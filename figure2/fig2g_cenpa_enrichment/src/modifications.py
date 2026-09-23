"""Extract per-read modified-base densities from indexed BAM files.

Modification locations are mapped from forward-read query coordinates to genomic
coordinates before region filtering. This avoids manual gap-string slicing and
makes insertions, deletions, and reverse-strand reads explicit.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping

import pysam

from .regions import Region


@dataclass(frozen=True)
class ModificationSpec:
    """Describe how a modified base is represented in a BAM file."""

    name: str
    canonical_base: str
    modification_code: str
    opportunity_motif: str
    minimum_score: int
    strand: int = 0

    def __post_init__(self) -> None:
        if not self.canonical_base:
            raise ValueError("canonical_base cannot be empty")

        if not self.opportunity_motif:
            raise ValueError("opportunity_motif cannot be empty")

        if not 0 <= self.minimum_score <= 255:
            raise ValueError(
                "minimum_score must be between 0 and 255"
            )


@dataclass(frozen=True)
class DensityRecord:
    """Store modified-base counts for one read and genomic region."""

    sample: str
    read_name: str
    chromosome: str
    region_start: int
    region_end: int
    modification: str
    modified_count: int
    opportunity_count: int
    density: float


def adenine_spec(
    minimum_score: int = 230,
) -> ModificationSpec:
    """Return the BAM tag convention used for N6-methyladenine."""

    return ModificationSpec(
        name="mA",
        canonical_base="A",
        modification_code="a",
        opportunity_motif="A",
        minimum_score=minimum_score,
    )


def cpg_cytosine_spec(
    minimum_score: int = 210,
) -> ModificationSpec:
    """Return the BAM tag convention used for CpG 5-methylcytosine.

    The denominator is the number of aligned CG motifs. If the intended
    denominator is every cytosine, change opportunity_motif to C.
    """

    return ModificationSpec(
        name="mC",
        canonical_base="C",
        modification_code="m",
        opportunity_motif="CG",
        minimum_score=minimum_score,
    )


def _forward_query_to_reference(
    read: pysam.AlignedSegment,
) -> list[int | None]:
    """Map forward-read query positions to reference coordinates."""

    reference_positions = list(
        read.get_reference_positions(full_length=True)
    )

    if read.is_reverse:
        reference_positions.reverse()

    return reference_positions


def _modification_calls(
    read: pysam.AlignedSegment,
    spec: ModificationSpec,
) -> list[tuple[int, int]] | None:
    """Return positions and scores for the requested modification."""

    modifications = read.modified_bases_forward or {}

    exact_key = (
        spec.canonical_base,
        spec.strand,
        spec.modification_code,
    )

    if exact_key in modifications:
        return modifications[exact_key]

    for key, calls in modifications.items():
        canonical_base, strand, modification_code = key

        if isinstance(canonical_base, bytes):
            canonical_base = canonical_base.decode("ascii")

        if isinstance(modification_code, bytes):
            modification_code = modification_code.decode("ascii")

        if (
            canonical_base.upper()
            == spec.canonical_base.upper()
            and strand == spec.strand
            and modification_code == spec.modification_code
        ):
            return calls

    return None


def _count_modified_bases(
    read: pysam.AlignedSegment,
    region: Region,
    spec: ModificationSpec,
    query_to_reference: list[int | None],
) -> int | None:
    """Count high-confidence modified bases inside a region."""

    calls = _modification_calls(
        read=read,
        spec=spec,
    )

    if calls is None:
        return None

    count = 0

    for query_position, score in calls:
        if score < spec.minimum_score:
            continue

        if query_position < 0:
            continue

        if query_position >= len(query_to_reference):
            continue

        reference_position = query_to_reference[query_position]

        if reference_position is None:
            continue

        if region.start <= reference_position < region.end:
            count += 1

    return count


def _count_opportunities(
    read: pysam.AlignedSegment,
    region: Region,
    spec: ModificationSpec,
    query_to_reference: list[int | None],
) -> int:
    """Count aligned target motifs falling completely inside a region."""

    sequence = (
        read.get_forward_sequence() or ""
    ).upper()

    motif = spec.opportunity_motif.upper()
    motif_length = len(motif)
    count = 0

    final_start = len(sequence) - motif_length + 1

    for query_start in range(0, final_start):
        query_end = query_start + motif_length

        if sequence[query_start:query_end] != motif:
            continue

        reference_slice = query_to_reference[
            query_start:query_end
        ]

        if len(reference_slice) != motif_length:
            continue

        if any(
            position is None
            for position in reference_slice
        ):
            continue

        reference_coordinates = [
            position
            for position in reference_slice
            if position is not None
        ]

        if not all(
            region.start <= position < region.end
            for position in reference_coordinates
        ):
            continue

        if motif_length > 1:
            interrupted = any(
                abs(right - left) != 1
                for left, right in zip(
                    reference_coordinates,
                    reference_coordinates[1:],
                )
            )

            if interrupted:
                continue

        count += 1

    return count


def density_for_read(
    read: pysam.AlignedSegment,
    region: Region,
    spec: ModificationSpec,
    *,
    sample: str,
) -> DensityRecord | None:
    """Calculate one read's modified-base density inside one region."""

    query_to_reference = _forward_query_to_reference(read)

    opportunity_count = _count_opportunities(
        read=read,
        region=region,
        spec=spec,
        query_to_reference=query_to_reference,
    )

    if opportunity_count == 0:
        return None

    modified_count = _count_modified_bases(
        read=read,
        region=region,
        spec=spec,
        query_to_reference=query_to_reference,
    )

    if modified_count is None:
        return None

    return DensityRecord(
        sample=sample,
        read_name=read.query_name,
        chromosome=region.chromosome,
        region_start=region.start,
        region_end=region.end,
        modification=spec.name,
        modified_count=modified_count,
        opportunity_count=opportunity_count,
        density=modified_count / opportunity_count,
    )


def iter_region_densities(
    bam_path: Path | str,
    regions: Mapping[str, Iterable[Region]],
    spec: ModificationSpec,
    *,
    sample: str,
    minimum_mapping_quality: int = 0,
) -> Iterator[DensityRecord]:
    """Yield density records for primary mapped reads."""

    if not 0 <= minimum_mapping_quality <= 255:
        raise ValueError(
            "minimum_mapping_quality must be between 0 and 255"
        )

    with pysam.AlignmentFile(
        str(bam_path),
        "rb",
    ) as bam_file:
        for chromosome, chromosome_regions in regions.items():
            for region in chromosome_regions:
                if region.chromosome != chromosome:
                    raise ValueError(
                        f"region chromosome "
                        f"{region.chromosome!r} does not match "
                        f"mapping key {chromosome!r}"
                    )

                for read in bam_file.fetch(
                    chromosome,
                    region.start,
                    region.end,
                ):
                    if read.is_unmapped:
                        continue

                    if read.is_secondary:
                        continue

                    if read.is_supplementary:
                        continue

                    if read.is_duplicate:
                        continue

                    if (
                        read.mapping_quality
                        < minimum_mapping_quality
                    ):
                        continue

                    record = density_for_read(
                        read=read,
                        region=region,
                        spec=spec,
                        sample=sample,
                    )

                    if record is not None:
                        yield record


def write_density_records(
    records: Iterable[DensityRecord],
    output_path: Path | str,
) -> None:
    """Write density records as a tidy CSV table."""

    path = Path(output_path)
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fieldnames = [
        "sample",
        "read_name",
        "chromosome",
        "region_start",
        "region_end",
        "modification",
        "modified_count",
        "opportunity_count",
        "density",
    ]

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()

        for record in records:
            writer.writerow(
                {
                    "sample": record.sample,
                    "read_name": record.read_name,
                    "chromosome": record.chromosome,
                    "region_start": record.region_start,
                    "region_end": record.region_end,
                    "modification": record.modification,
                    "modified_count": (
                        record.modified_count
                    ),
                    "opportunity_count": (
                        record.opportunity_count
                    ),
                    "density": record.density,
                }
            )


def read_density_records(
    path: Path | str,
) -> list[DensityRecord]:
    """Load a density CSV written by write_density_records."""

    records: list[DensityRecord] = []

    with Path(path).open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle)

        for row in reader:
            records.append(
                DensityRecord(
                    sample=row["sample"],
                    read_name=row["read_name"],
                    chromosome=row["chromosome"],
                    region_start=int(
                        row["region_start"]
                    ),
                    region_end=int(
                        row["region_end"]
                    ),
                    modification=row["modification"],
                    modified_count=int(
                        row["modified_count"]
                    ),
                    opportunity_count=int(
                        row["opportunity_count"]
                    ),
                    density=float(row["density"]),
                )
            )

    return records
