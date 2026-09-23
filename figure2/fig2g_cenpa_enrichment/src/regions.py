"""Read, validate, subtract, flank, and segment genomic regions.

All coordinates in this module follow the BED convention. Starts are zero-based,
ends are exclusive, and a region's length is therefore ``end - start``.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping


@dataclass(frozen=True, order=True)
class Region:
    """A validated half-open genomic interval."""

    chromosome: str
    start: int
    end: int

    def __post_init__(self) -> None:
        if not self.chromosome:
            raise ValueError("chromosome cannot be empty")
        if self.start < 0:
            raise ValueError(f"region start cannot be negative: {self}")
        if self.end <= self.start:
            raise ValueError(f"region end must be greater than start: {self}")

    @property
    def length(self) -> int:
        return self.end - self.start


RegionMap = dict[str, list[Region]]


def read_bed(path: Path | str) -> RegionMap:
    """Read the first three columns of a BED file into chromosome groups."""

    regions: dict[str, list[Region]] = defaultdict(list)
    bed_path = Path(path)

    with bed_path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue

            fields = line.split()
            if len(fields) < 3:
                raise ValueError(
                    f"{bed_path}:{line_number} has fewer than three columns"
                )

            chromosome = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{bed_path}:{line_number} contains non-integer coordinates"
                ) from error

            regions[chromosome].append(Region(chromosome, start, end))

    return {chromosome: sorted(items) for chromosome, items in regions.items()}


def merge_regions(regions: Mapping[str, Iterable[Region]]) -> RegionMap:
    """Merge overlapping or directly adjacent regions by chromosome."""

    merged: RegionMap = {}
    for chromosome, chromosome_regions in regions.items():
        ordered = sorted(chromosome_regions, key=lambda region: region.start)
        if not ordered:
            merged[chromosome] = []
            continue

        result = [ordered[0]]
        for current in ordered[1:]:
            previous = result[-1]
            if current.start <= previous.end:
                result[-1] = Region(
                    chromosome, previous.start, max(previous.end, current.end),
                )
            else:
                result.append(current)
        merged[chromosome] = result

    return merged


def subtract_regions(
    containing_regions: Mapping[str, Iterable[Region]],
    excluded_regions: Mapping[str, Iterable[Region]],
) -> RegionMap:
    """Return portions of containing regions not covered by excluded regions.

    This supports multiple active-array blocks per chromosome and clips exclusions
    at each containing interval's boundaries.
    """

    containers = merge_regions(containing_regions)
    exclusions = merge_regions(excluded_regions)
    result: RegionMap = {}

    for chromosome, chromosome_containers in containers.items():
        chromosome_result: list[Region] = []
        chromosome_exclusions = exclusions.get(chromosome, [])

        for container in chromosome_containers:
            cursor = container.start
            for exclusion in chromosome_exclusions:
                if exclusion.end <= cursor:
                    continue
                if exclusion.start >= container.end:
                    break

                clipped_start = max(exclusion.start, container.start)
                clipped_end = min(exclusion.end, container.end)
                if cursor < clipped_start:
                    chromosome_result.append(Region(chromosome, cursor, clipped_start))
                cursor = max(cursor, clipped_end)
                if cursor >= container.end:
                    break

            if cursor < container.end:
                chromosome_result.append(Region(chromosome, cursor, container.end))

        result[chromosome] = chromosome_result

    return result


def flank_regions(
    regions: Mapping[str, Iterable[Region]], flank_size: int = 1_000,
) -> RegionMap:
    """Create non-negative left and right flanks for each input region."""

    if flank_size <= 0:
        raise ValueError("flank_size must be greater than zero")

    flanks: RegionMap = {}
    for chromosome, chromosome_regions in regions.items():
        chromosome_flanks: list[Region] = []
        for region in chromosome_regions:
            left_start = max(0, region.start - flank_size)
            if left_start < region.start:
                chromosome_flanks.append(Region(chromosome, left_start, region.start))
            chromosome_flanks.append(
                Region(chromosome, region.end, region.end + flank_size)
            )
        flanks[chromosome] = chromosome_flanks
    return flanks


def split_regions(
    regions: Mapping[str, Iterable[Region]],
    segment_size: int = 5_000,
    *,
    include_partial: bool = False,
) -> RegionMap:
    """Split regions into fixed-size windows.

    Partial trailing windows are omitted by default so every returned segment has
    the same denominator length.
    """

    if segment_size <= 0:
        raise ValueError("segment_size must be greater than zero")

    segmented: RegionMap = {}
    for chromosome, chromosome_regions in regions.items():
        chromosome_segments: list[Region] = []
        for region in chromosome_regions:
            start = region.start
            while start + segment_size <= region.end:
                chromosome_segments.append(
                    Region(chromosome, start, start + segment_size)
                )
                start += segment_size

            if include_partial and start < region.end:
                chromosome_segments.append(Region(chromosome, start, region.end))

        segmented[chromosome] = chromosome_segments

    return segmented


def limit_regions(regions: Mapping[str, Iterable[Region]], limit: int) -> RegionMap:
    """Keep at most the first ``limit`` regions from each chromosome."""

    if limit < 0:
        raise ValueError("limit cannot be negative")
    return {
        chromosome: list(chromosome_regions)[:limit]
        for chromosome, chromosome_regions in regions.items()
    }


def count_regions(regions: Mapping[str, Iterable[Region]]) -> dict[str, int]:
    """Count regions by chromosome."""

    return {
        chromosome: sum(1 for _ in chromosome_regions)
        for chromosome, chromosome_regions in regions.items()
    }


def validate_region_bounds(
    regions: Mapping[str, Iterable[Region]], reference_lengths: Mapping[str, int],
) -> None:
    """Verify that every region exists within the reference assembly."""

    errors: list[str] = []
    for chromosome, chromosome_regions in regions.items():
        if chromosome not in reference_lengths:
            errors.append(f"{chromosome} is absent from the reference")
            continue

        chromosome_length = reference_lengths[chromosome]
        for region in chromosome_regions:
            if region.end > chromosome_length:
                errors.append(f"{region} exceeds chromosome length {chromosome_length}")
            if len(errors) >= 10:
                break
        if len(errors) >= 10:
            break

    if errors:
        details = "\n  ".join(errors)
        raise ValueError(f"Regions are incompatible with the reference\n  {details}")
