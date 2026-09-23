"""Configuration objects for the modified-base density analysis.

Paths are intentionally relative to the project root. This keeps workstation and
cluster locations out of source control while still making the expected inputs
clear to other users.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping


@dataclass(frozen=True)
class SampleConfig:
    """A named BAM input."""

    name: str
    bam_path: Path


@dataclass(frozen=True)
class AnalysisConfig:
    """Inputs and analysis settings shared by the pipeline."""

    reference_fasta: Path
    cdr_bed: Path
    active_array_bed: Path
    output_dir: Path
    samples: tuple[SampleConfig, ...]
    segment_size: int = 5_000
    non_cdr_segments_per_chromosome: int = 15
    min_mapping_quality: int = 0
    modification_score_thresholds: Mapping[str, int] = field(
        default_factory=lambda: {"mA": 230, "mC": 210}
    )

    def validate(self) -> None:
        """Raise a helpful error when an input or setting is invalid."""

        if self.segment_size <= 0:
            raise ValueError("segment_size must be greater than zero")
        if self.non_cdr_segments_per_chromosome <= 0:
            raise ValueError(
                "non_cdr_segments_per_chromosome must be greater than zero"
            )
        if not 0 <= self.min_mapping_quality <= 255:
            raise ValueError("min_mapping_quality must be between 0 and 255")

        required_files = {
            "reference FASTA": self.reference_fasta,
            "CDR BED": self.cdr_bed,
            "active-array BED": self.active_array_bed,
            **{f"{sample.name} BAM": sample.bam_path for sample in self.samples},
        }
        missing = [
            f"{label} ({path})"
            for label, path in required_files.items()
            if not path.is_file()
        ]
        if missing:
            joined = "\n  ".join(missing)
            raise FileNotFoundError(f"Missing required input files\n  {joined}")

    def prepare_output_directory(self) -> Path:
        """Create and return the output directory."""

        self.output_dir.mkdir(parents=True, exist_ok=True)
        return self.output_dir


def default_config(project_root: Path | str = Path.cwd()) -> AnalysisConfig:
    """Return an editable, repository-friendly example configuration.

    Adjust the relative filenames to match the repository or construct an
    ``AnalysisConfig`` directly in a command-line entry point.
    """

    root = Path(project_root).expanduser().resolve()
    data_dir = root / "data"

    return AnalysisConfig(
        reference_fasta=data_dir / "reference" / "hg002v1.0.1.fasta",
        cdr_bed=data_dir / "regions" / "cdr_regions.bed",
        active_array_bed=data_dir / "regions" / "active_arrays.bed",
        output_dir=root / "results",
        samples=(
            SampleConfig("CENPA", data_dir / "bam" / "cenpa.bam"),
            SampleConfig("CENPC", data_dir / "bam" / "cenpc.bam"),
            SampleConfig("H3K9me3", data_dir / "bam" / "h3k9me3.bam"),
        ),
    )
