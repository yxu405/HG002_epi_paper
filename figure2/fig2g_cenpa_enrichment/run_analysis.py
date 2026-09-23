#!/usr/bin/env python3
"""Run the CDR and non-CDR modified-base density workflow."""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import asdict
from pathlib import Path
from typing import Sequence

import pysam

from src.analysis import (
    coefficient_of_variation_between_chromosomes,
    fit_gini_threshold,
    fit_thresholds_by_chromosome,
    write_threshold_summary,
)
from src.config import AnalysisConfig, SampleConfig
from src.modifications import (
    ModificationSpec,
    adenine_spec,
    cpg_cytosine_spec,
    iter_region_densities,
    write_density_records,
)
from src.regions import (
    RegionMap,
    limit_regions,
    read_bed,
    split_regions,
    subtract_regions,
    validate_region_bounds,
)


def _sample_argument(value: str) -> SampleConfig:
    """Parse a sample argument written as NAME=/path/to/file.bam."""

    try:
        name, path = value.split("=", maxsplit=1)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "samples must use NAME=/path/to/file.bam"
        ) from error

    if not name.strip() or not path.strip():
        raise argparse.ArgumentTypeError(
            "samples must include both a name and BAM path"
        )

    return SampleConfig(
        name=name.strip(),
        bam_path=Path(path).expanduser(),
    )


def _slug(value: str) -> str:
    """Convert a sample name into a safe output filename component."""

    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return slug.strip("_") or "sample"


def _reference_lengths(reference_fasta: Path) -> dict[str, int]:
    """Read chromosome names and lengths from an indexed FASTA file."""

    try:
        with pysam.FastaFile(str(reference_fasta)) as reference:
            return dict(zip(reference.references, reference.lengths))
    except OSError as error:
        raise OSError(
            f"Could not open {reference_fasta}. Create its FASTA index with "
            f"'samtools faidx {reference_fasta}'."
        ) from error


def _prepare_regions(
    config: AnalysisConfig,
) -> tuple[RegionMap, RegionMap]:
    """Read, validate, subtract, and segment CDR and non-CDR regions."""

    cdr_regions = read_bed(config.cdr_bed)
    active_regions = read_bed(config.active_array_bed)

    reference_lengths = _reference_lengths(config.reference_fasta)
    validate_region_bounds(cdr_regions, reference_lengths)
    validate_region_bounds(active_regions, reference_lengths)

    non_cdr_regions = subtract_regions(
        containing_regions=active_regions,
        excluded_regions=cdr_regions,
    )

    segmented_cdr = split_regions(
        cdr_regions,
        segment_size=config.segment_size,
    )

    segmented_non_cdr = split_regions(
        non_cdr_regions,
        segment_size=config.segment_size,
    )

    segmented_non_cdr = limit_regions(
        segmented_non_cdr,
        limit=config.non_cdr_segments_per_chromosome,
    )

    return segmented_cdr, segmented_non_cdr


def _specs(args: argparse.Namespace) -> list[ModificationSpec]:
    """Construct the requested modified-base specifications."""

    requested = args.modification or ["mA", "mC"]
    specs: list[ModificationSpec] = []

    if "mA" in requested:
        specs.append(
            adenine_spec(
                minimum_score=args.ma_threshold,
            )
        )

    if "mC" in requested:
        specs.append(
            cpg_cytosine_spec(
                minimum_score=args.mc_threshold,
            )
        )

    return specs


def _summarize_group(records) -> dict[str, object]:
    """Create a JSON-compatible summary for one region group."""

    if not records:
        return {
            "record_count": 0,
        }

    variation = coefficient_of_variation_between_chromosomes(records)

    return {
        "record_count": len(records),
        "variation": asdict(variation),
    }


def run(
    config: AnalysisConfig,
    specs: Sequence[ModificationSpec],
) -> None:
    """Run every configured sample and modification type."""

    config.validate()
    output_dir = config.prepare_output_directory()

    cdr_regions, non_cdr_regions = _prepare_regions(config)

    analyses: dict[str, object] = {}

    summary: dict[str, object] = {
        "coordinate_convention": "BED zero-based half-open",
        "segment_size": config.segment_size,
        "non_cdr_segments_per_chromosome": (
            config.non_cdr_segments_per_chromosome
        ),
        "analyses": analyses,
    }

    for sample in config.samples:
        sample_slug = _slug(sample.name)

        for spec in specs:
            analysis_name = f"{sample.name}_{spec.name}"

            print(
                f"Processing {analysis_name}",
                flush=True,
            )

            cdr_records = list(
                iter_region_densities(
                    bam_path=sample.bam_path,
                    regions=cdr_regions,
                    spec=spec,
                    sample=sample.name,
                    minimum_mapping_quality=config.min_mapping_quality,
                )
            )

            non_cdr_records = list(
                iter_region_densities(
                    bam_path=sample.bam_path,
                    regions=non_cdr_regions,
                    spec=spec,
                    sample=sample.name,
                    minimum_mapping_quality=config.min_mapping_quality,
                )
            )

            prefix = f"{sample_slug}_{spec.name}"

            write_density_records(
                records=cdr_records,
                output_path=output_dir / f"{prefix}_cdr_densities.csv",
            )

            write_density_records(
                records=non_cdr_records,
                output_path=output_dir / f"{prefix}_non_cdr_densities.csv",
            )

            if not cdr_records:
                raise ValueError(
                    f"{analysis_name} produced no CDR density records"
                )

            if not non_cdr_records:
                raise ValueError(
                    f"{analysis_name} produced no non-CDR density records"
                )

            overall_threshold = fit_gini_threshold(
                non_cdr_densities=(
                    record.density for record in non_cdr_records
                ),
                cdr_densities=(
                    record.density for record in cdr_records
                ),
            )

            chromosome_thresholds = fit_thresholds_by_chromosome(
                non_cdr_records=non_cdr_records,
                cdr_records=cdr_records,
            )

            write_threshold_summary(
                thresholds=chromosome_thresholds,
                output_path=(
                    output_dir
                    / f"{prefix}_thresholds_by_chromosome.csv"
                ),
            )

            analyses[analysis_name] = {
                "specification": asdict(spec),
                "cdr": _summarize_group(cdr_records),
                "non_cdr": _summarize_group(non_cdr_records),
                "overall_threshold": asdict(overall_threshold),
                "chromosomes_with_thresholds": len(
                    chromosome_thresholds
                ),
            }

    summary_path = output_dir / "analysis_summary.json"

    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print(
        f"Finished. Results were written to {output_dir}",
        flush=True,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Calculate per-read mA and mC densities in segmented "
            "CDR and non-CDR regions, then fit Gini thresholds."
        )
    )

    parser.add_argument(
        "--reference",
        type=Path,
        required=True,
        help="Indexed reference FASTA file.",
    )

    parser.add_argument(
        "--cdr-bed",
        type=Path,
        required=True,
        help="BED file containing CDR intervals.",
    )

    parser.add_argument(
        "--active-bed",
        type=Path,
        required=True,
        help="BED file containing active-array intervals.",
    )

    parser.add_argument(
        "--sample",
        type=_sample_argument,
        action="append",
        required=True,
        metavar="NAME=BAM",
        help=(
            "Named BAM input written as NAME=/path/to/file.bam. "
            "Repeat this option for multiple samples."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for CSV and JSON results.",
    )

    parser.add_argument(
        "--segment-size",
        type=int,
        default=5_000,
        help="Size of each genomic segment in base pairs.",
    )

    parser.add_argument(
        "--non-cdr-limit",
        type=int,
        default=15,
        help=(
            "Maximum number of non-CDR segments retained per "
            "chromosome."
        ),
    )

    parser.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum read mapping quality.",
    )

    parser.add_argument(
        "--modification",
        action="append",
        choices=("mA", "mC"),
        help=(
            "Modification to process. Repeat for multiple types. "
            "Omit this option to process both mA and mC."
        ),
    )

    parser.add_argument(
        "--ma-threshold",
        type=int,
        default=230,
        help="Minimum modified-base score for mA calls.",
    )

    parser.add_argument(
        "--mc-threshold",
        type=int,
        default=210,
        help="Minimum modified-base score for mC calls.",
    )

    return parser.parse_args()


def main() -> None:
    """Build the configuration and run the analysis."""

    args = parse_args()

    config = AnalysisConfig(
        reference_fasta=args.reference.expanduser(),
        cdr_bed=args.cdr_bed.expanduser(),
        active_array_bed=args.active_bed.expanduser(),
        output_dir=args.output_dir.expanduser(),
        samples=tuple(args.sample),
        segment_size=args.segment_size,
        non_cdr_segments_per_chromosome=args.non_cdr_limit,
        min_mapping_quality=args.min_mapq,
        modification_score_thresholds={
            "mA": args.ma_threshold,
            "mC": args.mc_threshold,
        },
    )

    run(
        config=config,
        specs=_specs(args),
    )


if __name__ == "__main__":
    main()
