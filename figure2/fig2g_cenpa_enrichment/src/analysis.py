"""Statistical summaries for modified-base density records."""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import fmean, pstdev
from typing import Iterable, Mapping, Sequence

from .modifications import DensityRecord


@dataclass(frozen=True)
class VariationSummary:
    """Between-chromosome variation calculated from chromosome means."""

    chromosome_means: Mapping[str, float]
    mean_of_chromosome_means: float
    population_standard_deviation: float
    coefficient_of_variation: float


@dataclass(frozen=True)
class ThresholdResult:
    """A one-dimensional decision rule selected with weighted Gini impurity."""

    threshold: float
    weighted_gini: float
    left_label: str
    right_label: str
    left_purity: float
    right_purity: float
    sample_count: int

    def predict(self, density: float) -> str:
        return self.left_label if density <= self.threshold else self.right_label


def densities_by_chromosome(
    records: Iterable[DensityRecord],
) -> dict[str, list[float]]:
    """Collect record densities by chromosome."""

    grouped: dict[str, list[float]] = defaultdict(list)
    for record in records:
        grouped[record.chromosome].append(record.density)
    return dict(grouped)


def coefficient_of_variation_between_chromosomes(
    records: Iterable[DensityRecord],
) -> VariationSummary:
    """Calculate population CV across per-chromosome mean densities."""

    grouped = densities_by_chromosome(records)
    chromosome_means = {
        chromosome: fmean(densities)
        for chromosome, densities in grouped.items()
        if densities
    }
    if not chromosome_means:
        raise ValueError("at least one chromosome with density data is required")

    means = list(chromosome_means.values())
    overall_mean = fmean(means)
    standard_deviation = pstdev(means)
    if overall_mean == 0:
        raise ZeroDivisionError(
            "coefficient of variation is undefined when the mean is zero"
        )

    return VariationSummary(
        chromosome_means=chromosome_means,
        mean_of_chromosome_means=overall_mean,
        population_standard_deviation=standard_deviation,
        coefficient_of_variation=standard_deviation / overall_mean,
    )


def _gini(labels: Sequence[str]) -> float:
    if not labels:
        return 0.0
    counts: dict[str, int] = defaultdict(int)
    for label in labels:
        counts[label] += 1
    sample_count = len(labels)
    return 1.0 - sum((count / sample_count) ** 2 for count in counts.values())


def _majority_label(labels: Sequence[str]) -> tuple[str, float]:
    counts: dict[str, int] = defaultdict(int)
    for label in labels:
        counts[label] += 1
    # Sorting makes ties deterministic and therefore reproducible.
    label, count = min(counts.items(), key=lambda item: (-item[1], item[0]))
    return label, count / len(labels)


def fit_gini_threshold(
    non_cdr_densities: Iterable[float], cdr_densities: Iterable[float],
) -> ThresholdResult:
    """Fit a density threshold that separates non-CDR and CDR observations.

    The class assigned to each side is learned from the data instead of assuming
    in advance that CDR density must be higher or lower.
    """

    observations = [
        *((float(density), "non_CDR") for density in non_cdr_densities),
        *((float(density), "CDR") for density in cdr_densities),
    ]
    if not any(label == "non_CDR" for _, label in observations):
        raise ValueError("at least one non-CDR density is required")
    if not any(label == "CDR" for _, label in observations):
        raise ValueError("at least one CDR density is required")

    observations.sort(key=lambda item: item[0])
    unique_densities = sorted({density for density, _ in observations})
    if len(unique_densities) == 1:
        # No density threshold can split identical observations. Retaining one
        # deterministic rule makes the limitation visible in the result.
        threshold_candidates = unique_densities
    else:
        threshold_candidates = [
            (left + right) / 2
            for left, right in zip(unique_densities, unique_densities[1:])
        ]

    best_result: ThresholdResult | None = None
    total = len(observations)
    for threshold in threshold_candidates:
        left_labels = [label for density, label in observations if density <= threshold]
        right_labels = [label for density, label in observations if density > threshold]

        weighted_gini = len(left_labels) / total * _gini(left_labels) + len(
            right_labels
        ) / total * _gini(right_labels)
        left_label, left_purity = _majority_label(left_labels)
        if right_labels:
            right_label, right_purity = _majority_label(right_labels)
        else:
            right_label, right_purity = left_label, 0.0

        candidate = ThresholdResult(
            threshold=threshold,
            weighted_gini=weighted_gini,
            left_label=left_label,
            right_label=right_label,
            left_purity=left_purity,
            right_purity=right_purity,
            sample_count=total,
        )
        if best_result is None or (candidate.weighted_gini, candidate.threshold,) < (
            best_result.weighted_gini,
            best_result.threshold,
        ):
            best_result = candidate

    if best_result is None:
        raise RuntimeError("no threshold candidate could be evaluated")
    return best_result


def fit_thresholds_by_chromosome(
    non_cdr_records: Iterable[DensityRecord], cdr_records: Iterable[DensityRecord],
) -> dict[str, ThresholdResult]:
    """Fit thresholds for chromosomes represented in both groups."""

    non_cdr = densities_by_chromosome(non_cdr_records)
    cdr = densities_by_chromosome(cdr_records)
    shared_chromosomes = sorted(non_cdr.keys() & cdr.keys())
    return {
        chromosome: fit_gini_threshold(non_cdr[chromosome], cdr[chromosome])
        for chromosome in shared_chromosomes
    }


def write_threshold_summary(
    thresholds: Mapping[str, ThresholdResult], output_path: Path | str,
) -> None:
    """Write per-chromosome threshold results as CSV."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "chromosome",
                "threshold",
                "weighted_gini",
                "left_label",
                "right_label",
                "left_purity",
                "right_purity",
                "sample_count",
            ]
        )
        for chromosome, result in sorted(thresholds.items()):
            writer.writerow(
                [
                    chromosome,
                    result.threshold,
                    result.weighted_gini,
                    result.left_label,
                    result.right_label,
                    result.left_purity,
                    result.right_purity,
                    result.sample_count,
                ]
            )
