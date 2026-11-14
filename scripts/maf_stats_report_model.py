"""Pure cohort-classification model for :mod:`scripts.maf_stats_report`.

The report driver owns file schemas and presentation.  This module owns the
numeric parsing and robust cohort statistics that decide whether a sample is
normal, unusual, or needs review.  Keeping it free of TSV/HTML concerns makes
the classification contract independently testable and reusable.
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, field


NA_VALUES = frozenset({"", "NA", "NaN", "nan", "None"})
MAD_TO_SD = 1.4826
DEFAULT_MIN_SAMPLES = 5
DEFAULT_Z_THRESHOLD = 3.5
DEFAULT_RECURRENCE_WINDOW_BP = 500_000

# Direction in which a metric is concerning.  Insertion order controls report
# flag ordering and the corresponding strip plots. ``overlapping_reference_bp``
# and ``reference_contig_jumps`` are deliberately absent: AnchorWave's expected
# reference-chromosome-to-counterpart mapping makes both zero, so the report
# shows them without treating them as useful cohort flags.
FLAG_DIRECTIONS: dict[str, str] = {
    "aligned_reference_pct": "low",
    "aligned_query_pct": "low",
    "identity_pct": "low",
    "assembly_contig_n50_bp": "low",
    "overlapping_query_bp": "high",
    "breakpoints_per_gb_aligned": "high",
    "breakpoints_private": "high",
    "assembly_n_pct": "high",
    "unaligned_contig_bp": "high",
}

SCALE_ABSOLUTE_FLOORS: dict[str, float] = {
    "aligned_reference_pct": 1.0,
    "aligned_query_pct": 1.0,
    "identity_pct": 1.0,
    "assembly_n_pct": 0.5,
    "breakpoints_private": 2.0,
    "breakpoints_per_gb_aligned": 1.0,
    "overlapping_query_bp": 100_000.0,
    "assembly_contig_n50_bp": 1_000_000.0,
    "unaligned_contig_bp": 1_000_000.0,
}

SCALE_RELATIVE_FLOORS: dict[str, float] = {
    "breakpoints_per_gb_aligned": 0.05,
    "overlapping_query_bp": 0.05,
    "assembly_contig_n50_bp": 0.05,
    "unaligned_contig_bp": 0.05,
}

THRESHOLD_OPTIONS: dict[str, str] = {
    "flag_min_aligned_reference": "aligned_reference_pct",
    "flag_min_identity": "identity_pct",
    "flag_max_breakpoints_per_gb": "breakpoints_per_gb_aligned",
}

STATUS_OK = "ok"
STATUS_REVIEW = "review"
STATUS_UNUSUAL = "unusual"
STATUS_ICONS = {STATUS_OK: "\u2714", STATUS_REVIEW: "\u26a0", STATUS_UNUSUAL: "\u25c6"}
STATUS_LABELS = {STATUS_OK: "normal", STATUS_REVIEW: "needs review", STATUS_UNUSUAL: "unusual"}


def parse_number(raw: str | None) -> float | None:
    """Parse a TSV numeric cell; missing, non-finite, and invalid values are absent."""
    if raw is None:
        return None
    text = raw.strip()
    if text in NA_VALUES:
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


@dataclass
class MetricStats:
    n: int
    median: float
    mad: float
    scale: float
    relative_enabled: bool


@dataclass
class FlagResult:
    """Cohort flags and robust statistics, preserving sample input order."""

    flags: list[dict[str, list[str]]]
    stats: dict[str, MetricStats | None] = field(default_factory=dict)
    skipped: dict[str, int] = field(default_factory=dict)
    unusual: list[dict[str, list[str]]] = field(default_factory=list)
    z_scores: list[dict[str, float]] = field(default_factory=list)

    @property
    def status(self) -> list[str]:
        return compute_status(self)


def robust_scale(metric: str, median: float, mad: float) -> float:
    """Return ``1.4826*MAD`` after absolute and relative floor safeguards."""
    scale = MAD_TO_SD * mad
    scale = max(scale, SCALE_ABSOLUTE_FLOORS.get(metric, 0.0))
    scale = max(scale, SCALE_RELATIVE_FLOORS.get(metric, 0.0) * abs(median))
    return scale if scale > 0 else 1.0


def compute_flags(
    summaries: list[dict[str, str]],
    *,
    min_samples: int = DEFAULT_MIN_SAMPLES,
    z_threshold: float = DEFAULT_Z_THRESHOLD,
    thresholds: dict[str, float] | None = None,
) -> FlagResult:
    """Compute robust-z and absolute-threshold flags for each sample."""
    thresholds = thresholds or {}
    for metric in thresholds:
        if metric not in FLAG_DIRECTIONS:
            raise ValueError(f"No flag direction defined for threshold metric {metric!r}")
    result = FlagResult(
        flags=[{} for _ in summaries],
        unusual=[{} for _ in summaries],
        z_scores=[{} for _ in summaries],
    )

    for metric, direction in FLAG_DIRECTIONS.items():
        values = [parse_number(s.get(metric)) for s in summaries]
        present = [v for v in values if v is not None]

        stats: MetricStats | None = None
        if present:
            median = statistics.median(present)
            mad = statistics.median(abs(v - median) for v in present)
            stats = MetricStats(
                n=len(present),
                median=median,
                mad=mad,
                scale=robust_scale(metric, median, mad),
                relative_enabled=len(present) >= min_samples,
            )
        result.stats[metric] = stats
        if stats is None or not stats.relative_enabled:
            result.skipped[metric] = len(present)

        threshold = thresholds.get(metric)
        for idx, value in enumerate(values):
            if value is None:
                continue
            reasons: list[str] = []
            if stats is not None and stats.relative_enabled:
                z = (value - stats.median) / stats.scale
                result.z_scores[idx][metric] = z
                tag = "(mad0_floor)" if stats.mad == 0 else ""
                bad = z < -z_threshold if direction == "low" else z > z_threshold
                good = z > z_threshold if direction == "low" else z < -z_threshold
                if bad:
                    reasons.append(f"robust_z={z:.2f}{tag}")
                elif good:
                    result.unusual[idx][metric] = [f"robust_z={z:.2f}{tag}"]
            if threshold is not None:
                if direction == "low" and value < threshold:
                    reasons.append(f"below_threshold_{threshold:g}")
                elif direction == "high" and value > threshold:
                    reasons.append(f"above_threshold_{threshold:g}")
            if reasons:
                result.flags[idx][metric] = reasons
    return result


def sample_status(
    sample_flags: dict[str, list[str]], sample_unusual: dict[str, list[str]] | None = None
) -> str:
    """Return ``review`` before ``unusual`` before ``ok``."""
    if sample_flags:
        return STATUS_REVIEW
    if sample_unusual:
        return STATUS_UNUSUAL
    return STATUS_OK


def compute_status(flag_result: FlagResult) -> list[str]:
    """Return one status per sample in input order."""
    unusual = flag_result.unusual or [{} for _ in flag_result.flags]
    return [sample_status(flags, odd) for flags, odd in zip(flag_result.flags, unusual)]


def format_flags(sample_flags: dict[str, list[str]]) -> str:
    """Serialize flags in configured metric order."""
    entries: list[str] = []
    for metric in FLAG_DIRECTIONS:
        for reason in sample_flags.get(metric, []):
            entries.append(f"{metric}:{reason}")
    return ";".join(entries)
