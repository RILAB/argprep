#!/usr/bin/env python3
"""Cross-sample MAF QC report.

Merges the per-sample outputs of ``scripts/maf_stats.py`` into one TSV (one row
per sample plus a machine-readable ``flags`` column) and one self-contained
HTML report (embedded CSS, inline SVG, no JavaScript).
"""
from __future__ import annotations

import argparse
import html
import math
import statistics
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

try:
    from scripts.common import open_text
except ModuleNotFoundError:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from scripts.common import open_text


# ── input schemas ────────────────────────────────────────────────────────────

SUMMARY_COLUMNS = [
    "sample",
    "reference_length_bp",
    "covered_reference_bp",
    "reference_coverage_pct",
    "reference_bp_aligned_to_query_base",
    "query_length_bp",
    "query_length_source",
    "covered_query_bp",
    "query_coverage_pct",
    "identity_matches",
    "identity_compared_columns",
    "identity_pct",
    "gap_columns",
    "alignment_columns",
    "gap_fraction_pct",
    "blocks",
    "block_n50_bp",
    "overlapping_reference_bp",
    "strand_flips",
    "reference_contig_jumps",
    "out_of_order_adjacencies",
    "breakpoint_adjacencies",
    "breakpoints_per_gb_covered",
    "min_block_bp",
    "overlap_tolerance_bp",
]

BY_REFERENCE_COLUMNS = [
    "sample",
    "reference_contig",
    "reference_length_bp",
    "covered_reference_bp",
    "reference_coverage_pct",
    "reference_bp_aligned_to_query_base",
    "identity_matches",
    "identity_compared_columns",
    "identity_pct",
    "gap_columns",
    "alignment_columns",
    "gap_fraction_pct",
    "blocks",
    "block_n50_bp",
    "overlapping_reference_bp",
    "breakpoint_adjacencies",
    "dotplot",
]

BY_QUERY_COLUMNS = [
    "sample",
    "query_contig",
    "query_length_bp",
    "query_length_source",
    "covered_query_bp",
    "query_coverage_pct",
    "blocks",
    "blocks_considered",
    "strand_flips",
    "reference_contig_jumps",
    "out_of_order_adjacencies",
    "breakpoint_adjacencies",
    "breakpoints_per_gb_covered",
]

NA_VALUES = frozenset({"", "NA", "NaN", "nan", "None"})

# Columns rendered as text rather than numbers.
TEXT_COLUMNS = frozenset(
    {"sample", "reference_contig", "query_contig", "query_length_source", "dotplot"}
)
# Columns rendered with two decimals; all other numeric columns are integers.
FLOAT_COLUMNS = frozenset(
    {
        "reference_coverage_pct",
        "query_coverage_pct",
        "identity_pct",
        "gap_fraction_pct",
        "breakpoints_per_gb_covered",
    }
)

# ── flagging configuration (tune here) ───────────────────────────────────────

MAD_TO_SD = 1.4826
DEFAULT_MIN_SAMPLES = 5
DEFAULT_Z_THRESHOLD = 3.5

# Direction in which a metric is "bad". Order here is the order flags are listed.
FLAG_DIRECTIONS: dict[str, str] = {
    "reference_coverage_pct": "low",
    "query_coverage_pct": "low",
    "identity_pct": "low",
    "gap_fraction_pct": "high",
    "overlapping_reference_bp": "high",
    "blocks": "high",
    "breakpoints_per_gb_covered": "high",
}

# Absolute lower bound on the robust scale, in the metric's own units. Prevents
# division by zero and flags on trivially small deviations in tight cohorts.
SCALE_ABSOLUTE_FLOORS: dict[str, float] = {
    "reference_coverage_pct": 1.0,
    "query_coverage_pct": 1.0,
    "identity_pct": 1.0,
    "gap_fraction_pct": 1.0,
    "blocks": 10.0,
    "breakpoints_per_gb_covered": 1.0,
    "overlapping_reference_bp": 100_000.0,
}

# Additional scale floor as a fraction of |median| for count-like metrics.
SCALE_RELATIVE_FLOORS: dict[str, float] = {
    "blocks": 0.05,
    "breakpoints_per_gb_covered": 0.05,
    "overlapping_reference_bp": 0.05,
}

# CLI option dest -> metric for optional absolute thresholds. Direction comes
# from FLAG_DIRECTIONS.
THRESHOLD_OPTIONS: dict[str, str] = {
    "flag_min_reference_coverage": "reference_coverage_pct",
    "flag_min_identity": "identity_pct",
    "flag_max_gap_fraction": "gap_fraction_pct",
    "flag_max_breakpoints_per_gb": "breakpoints_per_gb_covered",
}

METRIC_LABELS: dict[str, str] = {
    "sample": "Sample",
    "reference_length_bp": "Reference length (bp)",
    "covered_reference_bp": "Covered ref bp",
    "reference_coverage_pct": "Ref coverage (%)",
    "reference_bp_aligned_to_query_base": "Ref bp aligned to query base",
    "query_length_bp": "Query length (bp)",
    "query_length_source": "Query length source",
    "covered_query_bp": "Covered query bp",
    "query_coverage_pct": "Query coverage (%)",
    "identity_matches": "Identity matches",
    "identity_compared_columns": "Identity compared columns",
    "identity_pct": "Identity (%)",
    "gap_columns": "Gap columns",
    "alignment_columns": "Alignment columns",
    "gap_fraction_pct": "Gap fraction (%)",
    "blocks": "Blocks",
    "blocks_considered": "Blocks considered",
    "block_n50_bp": "Block N50 (bp)",
    "overlapping_reference_bp": "Overlapping ref bp",
    "strand_flips": "Strand flips",
    "reference_contig_jumps": "Ref-contig jumps",
    "out_of_order_adjacencies": "Out-of-order adj.",
    "breakpoint_adjacencies": "Breakpoint adj.",
    "breakpoints_per_gb_covered": "Breakpoints / Gb covered",
    "min_block_bp": "Min block bp",
    "overlap_tolerance_bp": "Overlap tolerance bp",
    "reference_contig": "Reference contig",
    "query_contig": "Query contig",
    "dotplot": "Dotplot",
}

PRIMARY_TABLE_COLUMNS = [
    "reference_coverage_pct",
    "covered_reference_bp",
    "query_coverage_pct",
    "identity_pct",
    "gap_fraction_pct",
    "blocks",
    "block_n50_bp",
    "overlapping_reference_bp",
    "breakpoint_adjacencies",
    "breakpoints_per_gb_covered",
]

REFERENCE_TABLE_COLUMNS = [c for c in BY_REFERENCE_COLUMNS if c not in ("sample", "dotplot")]
QUERY_TABLE_COLUMNS = [c for c in BY_QUERY_COLUMNS if c != "sample"]

METRIC_DEFINITIONS: list[tuple[str, str]] = [
    ("reference_length_bp", "Total reference length from the reference <code>.fai</code>."),
    ("covered_reference_bp",
     "Union of reference intervals covered by any accepted block (merged block span); "
     "numerator of reference coverage."),
    ("reference_coverage_pct",
     "Covered reference bp divided by total reference length. Union-based: overlapping "
     "blocks are counted once."),
    ("reference_bp_aligned_to_query_base",
     "Reference bases opposite a non-gap query character; reported separately from "
     "covered reference bp."),
    ("query_length_bp", "Query genome size used as the query-coverage denominator."),
    ("query_length_source",
     "Source of the query-coverage denominator. <code>query_fai</code>: total length of "
     "all contigs in the query <code>.fai</code>. <code>maf_srcsize</code>: no query "
     "<code>.fai</code> was supplied, so only query contigs represented in the MAF (their "
     "<code>srcSize</code>) are counted. This <strong>overstates</strong> query coverage "
     "because unaligned query contigs are absent from the MAF; such values are marked "
     "&dagger; and are not comparable to <code>query_fai</code> values."),
    ("covered_query_bp", "Union of forward-coordinate query intervals covered by accepted blocks."),
    ("query_coverage_pct", "Covered query bp divided by query length (see denominator source)."),
    ("identity_matches", "Columns where both bases are A/C/G/T (case-insensitive) and equal."),
    ("identity_compared_columns",
     "Columns where both bases are A/C/G/T (case-insensitive); N and other ambiguity codes "
     "are excluded."),
    ("identity_pct", "Identity matches divided by identity compared columns. Raw block-column metric."),
    ("gap_columns", "Alignment columns with <code>-</code> in either row."),
    ("alignment_columns", "All alignment columns across accepted blocks."),
    ("gap_fraction_pct", "Gap columns divided by alignment columns. Raw block-column metric."),
    ("blocks", "Count of accepted pairwise MAF blocks. Raw block metric."),
    ("blocks_considered",
     "Blocks on this query contig at least <code>min_block_bp</code> long, i.e. those used "
     "for breakpoint classification."),
    ("block_n50_bp", "N50 of reference <code>size</code> across accepted blocks. Raw block metric."),
    ("overlapping_reference_bp",
     "Reference bp covered by more than one block (raw block span minus union span). "
     "Makes duplicated or secondary alignment content visible."),
    ("strand_flips",
     "Adjacent accepted blocks along a query contig whose query strands differ."),
    ("reference_contig_jumps",
     "Adjacent accepted blocks along a query contig aligned to different reference contigs."),
    ("out_of_order_adjacencies",
     "Same-reference, same-strand adjacent blocks whose reference coordinates contradict "
     "query order beyond <code>overlap_tolerance_bp</code>."),
    ("breakpoint_adjacencies",
     "Unique adjacent block pairs having any of the three breakpoint properties (strand "
     "flip, reference-contig jump, out-of-order). A pair with several properties counts once."),
    ("breakpoints_per_gb_covered",
     "Unique breakpoint adjacencies divided by covered reference bp, multiplied by 1e9."),
    ("min_block_bp",
     "Blocks shorter than this are excluded from breakpoint classification only, not from "
     "coverage or general alignment metrics."),
    ("overlap_tolerance_bp",
     "Reference overlap/backtrack allowed between adjacent same-strand blocks before the "
     "pair is classified out-of-order."),
    ("dotplot",
     "Structural dotplot of a reference contig against the query contigs aligned to it "
     "(forward blocks blue, reverse blocks red)."),
]


# ── reading ──────────────────────────────────────────────────────────────────


def parse_number(raw: str | None) -> float | None:
    """Parse a TSV numeric cell; ``NA``/empty/non-finite/garbage -> ``None``."""
    if raw is None:
        return None
    text = raw.strip()
    if text in NA_VALUES:
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    if not math.isfinite(value):
        return None
    return value


def read_tsv(path: Path, required: list[str]) -> list[dict[str, str]]:
    """Read a headered TSV into dicts, checking that ``required`` columns exist."""
    with open_text(path, "rt") as handle:
        header_line = handle.readline()
        if not header_line.strip():
            raise ValueError(f"{path}: empty file (missing header row)")
        header = header_line.rstrip("\r\n").split("\t")
        missing = [c for c in required if c not in header]
        if missing:
            raise ValueError(f"{path}: missing required column(s): {', '.join(missing)}")
        rows: list[dict[str, str]] = []
        for line_number, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) != len(header):
                raise ValueError(
                    f"{path}:{line_number}: expected {len(header)} columns, found {len(parts)}"
                )
            rows.append(dict(zip(header, parts)))
    return rows


def read_summary(path: Path) -> dict[str, str]:
    """Read a per-sample ``<sample>.maf_stats.tsv`` (header + exactly one row)."""
    rows = read_tsv(path, SUMMARY_COLUMNS)
    if len(rows) != 1:
        raise ValueError(f"{path}: expected exactly 1 data row, found {len(rows)}")
    row = rows[0]
    if not row["sample"].strip():
        raise ValueError(f"{path}: empty sample name")
    return {c: row[c] for c in SUMMARY_COLUMNS}


def read_contig_tables(
    paths: list[Path],
    required: list[str],
    samples: list[str],
    kind: str,
) -> dict[str, list[dict[str, str]]]:
    """Group rows of per-contig TSVs by their ``sample`` column.

    Raises if a row names a sample absent from the summaries, or if two files
    contribute rows for the same sample.
    """
    known = set(samples)
    grouped: dict[str, list[dict[str, str]]] = {}
    source: dict[str, Path] = {}
    for path in paths:
        rows = read_tsv(path, required)
        for row in rows:
            sample = row["sample"]
            if sample not in known:
                raise ValueError(
                    f"{path}: {kind} row for sample {sample!r}, which is not among the "
                    f"summary inputs"
                )
            if sample in source and source[sample] != path:
                raise ValueError(
                    f"{path}: duplicate {kind} input for sample {sample!r} "
                    f"(also in {source[sample]})"
                )
            source[sample] = path
            grouped.setdefault(sample, []).append({c: row[c] for c in required})
    return grouped


# ── flagging ─────────────────────────────────────────────────────────────────


@dataclass
class MetricStats:
    n: int
    median: float
    mad: float
    scale: float
    relative_enabled: bool


@dataclass
class FlagResult:
    # per sample (input order): metric -> list of reasons
    flags: list[dict[str, list[str]]]
    stats: dict[str, MetricStats | None] = field(default_factory=dict)
    # metrics whose relative flagging was skipped -> number of non-NA samples
    skipped: dict[str, int] = field(default_factory=dict)


def robust_scale(metric: str, median: float, mad: float) -> float:
    """``max(1.4826*MAD, absolute floor, relative floor * |median|)``."""
    scale = MAD_TO_SD * mad
    scale = max(scale, SCALE_ABSOLUTE_FLOORS.get(metric, 0.0))
    scale = max(scale, SCALE_RELATIVE_FLOORS.get(metric, 0.0) * abs(median))
    if not scale > 0:
        # Defensive: every flaggable metric has a positive absolute floor, but
        # never divide by zero if a new metric is added without one.
        scale = 1.0
    return scale


def compute_flags(
    summaries: list[dict[str, str]],
    *,
    min_samples: int = DEFAULT_MIN_SAMPLES,
    z_threshold: float = DEFAULT_Z_THRESHOLD,
    thresholds: dict[str, float] | None = None,
) -> FlagResult:
    """Compute relative (robust-z) and absolute flags for each sample.

    ``thresholds`` maps metric name -> absolute threshold; the direction is
    taken from ``FLAG_DIRECTIONS`` (low-is-bad -> flag when value < threshold).
    """
    thresholds = thresholds or {}
    for metric in thresholds:
        if metric not in FLAG_DIRECTIONS:
            raise ValueError(f"No flag direction defined for threshold metric {metric!r}")
    result = FlagResult(flags=[{} for _ in summaries])

    for metric, direction in FLAG_DIRECTIONS.items():
        values = [parse_number(s.get(metric)) for s in summaries]
        present = [v for v in values if v is not None]

        stats: MetricStats | None = None
        if present:
            med = statistics.median(present)
            mad = statistics.median(abs(v - med) for v in present)
            stats = MetricStats(
                n=len(present),
                median=med,
                mad=mad,
                scale=robust_scale(metric, med, mad),
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
                bad = z < -z_threshold if direction == "low" else z > z_threshold
                if bad:
                    tag = "(mad0_floor)" if stats.mad == 0 else ""
                    reasons.append(f"robust_z={z:.2f}{tag}")
            if threshold is not None:
                if direction == "low" and value < threshold:
                    reasons.append(f"below_threshold_{threshold:g}")
                elif direction == "high" and value > threshold:
                    reasons.append(f"above_threshold_{threshold:g}")
            if reasons:
                result.flags[idx][metric] = reasons
    return result


def format_flags(sample_flags: dict[str, list[str]]) -> str:
    """``metric:reason;metric:reason`` in FLAG_DIRECTIONS order; '' when none."""
    entries: list[str] = []
    for metric in FLAG_DIRECTIONS:
        for reason in sample_flags.get(metric, []):
            entries.append(f"{metric}:{reason}")
    return ";".join(entries)


# ── output TSV ───────────────────────────────────────────────────────────────


def write_summary_tsv(
    path: Path, summaries: list[dict[str, str]], flags: list[dict[str, list[str]]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(SUMMARY_COLUMNS + ["flags"]) + "\n")
        for row, sample_flags in zip(summaries, flags):
            fh.write(
                "\t".join([row[c] for c in SUMMARY_COLUMNS] + [format_flags(sample_flags)])
                + "\n"
            )


# ── HTML helpers ─────────────────────────────────────────────────────────────


def esc(text: object) -> str:
    return html.escape(str(text), quote=True)


def fmt_value(raw: str | None, column: str) -> str:
    """Format a cell for HTML (already escaped). NA -> muted 'NA'."""
    if column in TEXT_COLUMNS:
        return esc(raw or "")
    value = parse_number(raw)
    if value is None:
        if raw is not None and raw.strip() not in NA_VALUES:
            return esc(raw)
        return '<span class="na">NA</span>'
    if column in FLOAT_COLUMNS:
        return f"{value:,.2f}"
    if value.is_integer():
        return f"{int(value):,}"
    return f"{value:,.2f}"


def _fmt_axis(value: float, metric: str) -> str:
    if metric in FLOAT_COLUMNS:
        return f"{value:,.2f}"
    if abs(value - round(value)) < 1e-9:
        return f"{int(round(value)):,}"
    return f"{value:,.1f}"


def svg_strip_plot(
    metric: str,
    samples: list[str],
    values: list[float | None],
    flagged: list[bool],
    *,
    width: int = 900,
    height: int = 74,
) -> str:
    """One-dimensional strip plot: one dot per sample, median line, flagged dots red."""
    margin = {"left": 200, "right": 30, "top": 12, "bottom": 24}
    plot_w = width - margin["left"] - margin["right"]
    plot_h = height - margin["top"] - margin["bottom"]
    label = METRIC_LABELS.get(metric, metric)
    parts = [
        f'<svg class="strip" width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img" aria-label="{esc(label)}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{margin["left"] - 10}" y="{margin["top"] + plot_h / 2 + 4:.2f}" '
        f'text-anchor="end" font-size="12" font-family="sans-serif">{esc(label)}</text>',
    ]
    present = [v for v in values if v is not None]
    axis_y = margin["top"] + plot_h
    parts.append(
        f'<line x1="{margin["left"]}" y1="{axis_y}" x2="{margin["left"] + plot_w}" '
        f'y2="{axis_y}" stroke="#333" stroke-width="1"/>'
    )
    if not present:
        parts.append(
            f'<text x="{margin["left"] + plot_w / 2}" y="{margin["top"] + plot_h / 2 + 4:.2f}" '
            f'text-anchor="middle" font-size="11" font-family="sans-serif" fill="#777">'
            f"no values</text>"
        )
        parts.append("</svg>")
        return "\n".join(parts)

    lo, hi = min(present), max(present)
    if hi == lo:
        pad = abs(lo) * 0.05 or 1.0
        lo, hi = lo - pad, hi + pad

    def x_scale(v: float) -> float:
        return margin["left"] + (v - lo) / (hi - lo) * plot_w

    for frac in (0.0, 0.5, 1.0):
        val = lo + frac * (hi - lo)
        x = margin["left"] + frac * plot_w
        anchor = "start" if frac == 0 else ("end" if frac == 1 else "middle")
        parts.append(
            f'<line x1="{x:.2f}" y1="{axis_y}" x2="{x:.2f}" y2="{axis_y + 4}" '
            f'stroke="#333" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{x:.2f}" y="{axis_y + 16}" text-anchor="{anchor}" font-size="10" '
            f'font-family="sans-serif">{esc(_fmt_axis(val, metric))}</text>'
        )

    med = statistics.median(present)
    mx = x_scale(med)
    parts.append(
        f'<line x1="{mx:.2f}" y1="{margin["top"]}" x2="{mx:.2f}" y2="{axis_y}" '
        f'stroke="#888" stroke-width="1.5" stroke-dasharray="4 3">'
        f"<title>median {esc(_fmt_axis(med, metric))}</title></line>"
    )

    lanes = 5
    lane_h = plot_h / (lanes + 1)
    # Draw unflagged dots first so flagged ones sit on top.
    order = sorted(range(len(samples)), key=lambda i: flagged[i])
    for i in order:
        v = values[i]
        if v is None:
            continue
        cy = margin["top"] + lane_h * (1 + i % lanes)
        color = "#E45756" if flagged[i] else "#4C78A8"
        cls = "dot flagged" if flagged[i] else "dot"
        parts.append(
            f'<circle class="{cls}" cx="{x_scale(v):.2f}" cy="{cy:.2f}" r="4" fill="{color}" '
            f'fill-opacity="0.85" stroke="white" stroke-width="0.5">'
            f"<title>{esc(samples[i])}: {esc(_fmt_axis(v, metric))}</title></circle>"
        )
    parts.append("</svg>")
    return "\n".join(parts)


def _breakpoint_first(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Rows with breakpoint_adjacencies > 0 first (most first), others in input order."""
    def key(item: tuple[int, dict[str, str]]) -> tuple[int, float, int]:
        idx, row = item
        bp = parse_number(row.get("breakpoint_adjacencies")) or 0.0
        return (0 if bp > 0 else 1, -bp if bp > 0 else 0.0, idx)

    return [row for _, row in sorted(enumerate(rows), key=key)]


# Per-sample contig tables in the HTML stop here (breakpoint contigs come first);
# fragmented assemblies can have thousands of scaffolds. The TSVs keep every row.
MAX_HTML_CONTIG_ROWS = 50


def _contig_table(rows: list[dict[str, str]], columns: list[str], css_class: str) -> str:
    shown = rows[:MAX_HTML_CONTIG_ROWS]
    note = ""
    if len(rows) > len(shown):
        note = (
            f'<p class="na">Showing {len(shown):,} of {len(rows):,} contigs; '
            "the per-sample TSV has every row.</p>\n"
        )
    rows = shown
    out = [note, f'<table class="{css_class}">\n<tr>']
    out.extend(f"<th>{esc(METRIC_LABELS.get(c, c))}</th>" for c in columns)
    out.append("</tr>\n")
    for row in rows:
        bp = parse_number(row.get("breakpoint_adjacencies")) or 0.0
        cls = ' class="bp"' if bp > 0 else ""
        out.append(f"<tr{cls}>")
        for c in columns:
            out.append(f"<td>{fmt_value(row.get(c), c)}</td>")
        out.append("</tr>\n")
    out.append("</table>\n")
    return "".join(out)


CSS = """\
body{font-family:sans-serif;margin:24px;color:#111;max-width:1200px}
table{border-collapse:collapse;margin:12px 0}
th,td{border:1px solid #ccc;padding:4px 10px;text-align:right}
th{background:#f0f0f0;text-align:center}
td:first-child{text-align:left}
h1,h2,h3{margin-top:1.4em}
code{background:#f6f8fa;padding:0 3px;border-radius:3px}
details{margin:8px 0;border:1px solid #d0d7de;border-radius:6px;padding:4px 12px}
summary{cursor:pointer;font-size:1.05em;font-weight:bold;padding:6px 0;list-style:revert}
summary:hover{color:#0969da}
td.flag{background:#ffd7d7;font-weight:bold;cursor:help}
td.flags{text-align:left;font-family:monospace;font-size:0.85em;color:#a40000}
tr.bp td{background:#fff8dc}
.na{color:#999}
.note{background:#fff8dc;border:1px solid #e0c96b;border-radius:6px;padding:8px 12px}
.params td:first-child{font-family:monospace}
section.ref-contigs{border-left:4px solid #4C78A8;padding-left:12px;margin:12px 0}
section.query-contigs{border-left:4px solid #F58518;padding-left:12px;margin:12px 0}
section.dotplots{border-left:4px solid #54A24B;padding-left:12px;margin:12px 0}
.scroll{overflow-x:auto}
.thumbs{display:flex;flex-wrap:wrap;gap:12px}
.thumbs figure{margin:0;text-align:center;font-size:0.85em}
.thumbs img{border:1px solid #ccc}
dl.defs dt{font-family:monospace;font-weight:bold;margin-top:8px}
dl.defs dd{margin-left:20px}
.legend span{display:inline-block;width:10px;height:10px;border-radius:5px;margin:0 4px 0 12px}
"""


def render_html(
    summaries: list[dict[str, str]],
    flag_result: FlagResult,
    by_reference: dict[str, list[dict[str, str]]],
    by_query: dict[str, list[dict[str, str]]],
    *,
    min_samples: int,
    z_threshold: float,
    thresholds: dict[str, float],
    thumbnail_width: int = 280,
) -> str:
    samples = [s["sample"] for s in summaries]
    flags = flag_result.flags
    out: list[str] = []
    w = out.append

    w('<!doctype html>\n<html lang="en">\n<head>\n')
    w('<meta charset="utf-8" />\n<meta name="viewport" content="width=device-width, initial-scale=1" />\n')
    w("<title>ARGprep MAF QC</title>\n")
    w(f"<style>{CSS}</style>\n")
    w("</head>\n<body>\n")
    w("<h1>ARGprep MAF QC</h1>\n")
    w(
        f'<p>Generated: <code>{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</code>'
        f" &nbsp;|&nbsp; Samples: <code>{len(summaries):,}</code></p>\n"
    )

    # ── parameters ───────────────────────────────────────────────────────────
    w("<h2>Parameters</h2>\n<table class=\"params\">\n")
    w("<tr><th>Parameter</th><th>Value</th></tr>\n")
    w(f"<tr><td>min_samples</td><td>{min_samples:,}</td></tr>\n")
    w(f"<tr><td>z_threshold</td><td>{esc(f'{z_threshold:g}')}</td></tr>\n")
    for dest, metric in THRESHOLD_OPTIONS.items():
        option = "--" + dest.replace("_", "-")
        value = thresholds.get(metric)
        shown = esc(f"{value:g}") if value is not None else '<span class="na">off</span>'
        w(f"<tr><td>{esc(option)}</td><td>{shown}</td></tr>\n")
    for col in ("min_block_bp", "overlap_tolerance_bp"):
        distinct = list(dict.fromkeys(s[col] for s in summaries))
        shown = ", ".join(fmt_value(v, col) for v in distinct)
        if len(distinct) > 1:
            shown += ' &nbsp;<strong>(differs between samples)</strong>'
        w(f"<tr><td>{esc(col)}</td><td>{shown}</td></tr>\n")
    w("</table>\n")
    w(
        "<p>Relative flags: direction-aware robust z = (value &minus; median) / scale, "
        f"scale = max({MAD_TO_SD} &times; MAD, metric floor). A value is flagged when |z| "
        f"&gt; {esc(f'{z_threshold:g}')} in the bad direction. <code>(mad0_floor)</code> "
        "marks z-scores computed with MAD = 0, where the scale is the floor alone.</p>\n"
    )

    skipped = flag_result.skipped
    if len(summaries) < min_samples:
        w(
            f'<p class="note"><strong>Relative (cohort) flagging skipped:</strong> '
            f"only {len(summaries):,} sample(s), fewer than min_samples = "
            f"{min_samples:,}. Only absolute thresholds (if any) were applied.</p>\n"
        )
    elif skipped:
        items = ", ".join(f"<code>{esc(m)}</code> (n={n:,})" for m, n in skipped.items())
        w(
            f'<p class="note"><strong>Relative flagging skipped</strong> for metrics with '
            f"fewer than min_samples = {min_samples:,} non-NA values: {items}.</p>\n"
        )

    # ── primary table ────────────────────────────────────────────────────────
    w("<h2>Samples</h2>\n")
    w(
        "<p>Flagged cells are highlighted; hover for the reason. &dagger; marks query "
        "coverage computed against MAF <code>srcSize</code> (no query <code>.fai</code>), "
        "which overstates coverage.</p>\n"
    )
    w('<div class="scroll"><table class="samples">\n<tr><th>Sample</th>')
    for col in PRIMARY_TABLE_COLUMNS:
        w(f"<th>{esc(METRIC_LABELS.get(col, col))}</th>")
    w("<th>Flags</th></tr>\n")
    for row, sample_flags in zip(summaries, flags):
        sample = row["sample"]
        w(f'<tr><td><a href="#{esc(_sample_anchor(samples, sample))}">{esc(sample)}</a></td>')
        for col in PRIMARY_TABLE_COLUMNS:
            cell = fmt_value(row.get(col), col)
            if col == "query_coverage_pct" and row.get("query_length_source") == "maf_srcsize":
                cell += " &dagger;"
            reasons = sample_flags.get(col)
            if reasons:
                tip = "; ".join(f"{col}:{r}" for r in reasons)
                w(f'<td class="flag" title="{esc(tip)}">{cell}</td>')
            else:
                w(f"<td>{cell}</td>")
        w(f'<td class="flags">{esc(format_flags(sample_flags).replace(";", "; "))}</td></tr>\n')
    w("</table></div>\n")

    # ── strip plots ──────────────────────────────────────────────────────────
    w("<h2>Distributions across samples</h2>\n")
    w(
        '<p class="legend">One dot per sample (hover for name); dashed line = median.'
        '<span style="background:#4C78A8"></span>not flagged'
        '<span style="background:#E45756"></span>flagged</p>\n'
    )
    for metric in FLAG_DIRECTIONS:
        values = [parse_number(s.get(metric)) for s in summaries]
        flagged = [metric in f for f in flags]
        w(svg_strip_plot(metric, samples, values, flagged))
        w("\n")

    # ── per-sample details ───────────────────────────────────────────────────
    w("<h2>Per-sample detail</h2>\n")
    for row, sample_flags in zip(summaries, flags):
        sample = row["sample"]
        n_flags = sum(len(v) for v in sample_flags.values())
        w(f'<details id="{esc(_sample_anchor(samples, sample))}">\n')
        ref_cov = fmt_value(row.get("reference_coverage_pct"), "reference_coverage_pct")
        ident = fmt_value(row.get("identity_pct"), "identity_pct")
        pct = lambda cell: cell if "NA" in cell else cell + "%"  # noqa: E731
        w(
            f"<summary>{esc(sample)}"
            f" &mdash; ref coverage {pct(ref_cov)}"
            f" &nbsp;|&nbsp; identity {pct(ident)}"
            f" &nbsp;|&nbsp; {fmt_value(row.get('breakpoint_adjacencies'), 'breakpoint_adjacencies')}"
            f" breakpoint adj. &nbsp;|&nbsp; {n_flags} flag(s)</summary>\n"
        )

        w('<table class="sample-summary">\n<tr><th>Metric</th><th>Value</th></tr>\n')
        for col in SUMMARY_COLUMNS[1:]:
            reasons = sample_flags.get(col)
            cell = fmt_value(row.get(col), col)
            if reasons:
                tip = "; ".join(f"{col}:{r}" for r in reasons)
                w(f'<tr><td>{esc(col)}</td><td class="flag" title="{esc(tip)}">{cell}</td></tr>\n')
            else:
                w(f"<tr><td>{esc(col)}</td><td>{cell}</td></tr>\n")
        w("</table>\n")

        ref_rows = _breakpoint_first(by_reference.get(sample, []))
        w('<section class="ref-contigs">\n<h3>Reference contigs</h3>\n')
        if ref_rows:
            w("<p>Contigs with breakpoint adjacencies are listed first and shaded.</p>\n")
            w('<div class="scroll">')
            w(_contig_table(ref_rows, REFERENCE_TABLE_COLUMNS, "ref-table"))
            w("</div>\n")
        else:
            w('<p class="na">No reference-contig table supplied.</p>\n')
        w("</section>\n")

        query_rows = _breakpoint_first(by_query.get(sample, []))
        w('<section class="query-contigs">\n<h3>Query contigs</h3>\n')
        if query_rows:
            w("<p>Contigs with breakpoint adjacencies are listed first and shaded.</p>\n")
            w('<div class="scroll">')
            w(_contig_table(query_rows, QUERY_TABLE_COLUMNS, "query-table"))
            w("</div>\n")
        else:
            w('<p class="na">No query-contig table supplied.</p>\n')
        w("</section>\n")

        plots = [r for r in ref_rows if (r.get("dotplot") or "").strip()]
        w('<section class="dotplots">\n<h3>Dotplots</h3>\n')
        if plots:
            w('<div class="thumbs">\n')
            for r in plots:
                path = esc(r["dotplot"].strip())
                contig = esc(r["reference_contig"])
                bp = fmt_value(r.get("breakpoint_adjacencies"), "breakpoint_adjacencies")
                w(
                    f'<figure><a href="{path}"><img src="{path}" loading="lazy" '
                    f'width="{thumbnail_width}" alt="{contig} dotplot"></a>'
                    f"<figcaption>{contig} ({bp} breakpoint adj.)</figcaption></figure>\n"
                )
            w("</div>\n")
        else:
            w('<p class="na">No dotplots drawn for this sample.</p>\n')
        w("</section>\n")
        w("</details>\n")

    # ── definitions ──────────────────────────────────────────────────────────
    w("<h2>Metric definitions</h2>\n")
    w(
        '<p class="note"><strong>Raw versus union-based metrics.</strong> Identity, gap '
        "fraction, block count, and block N50 are raw block-column metrics: overlapping or "
        "secondary blocks are counted more than once. Reference and query coverage are "
        "union-based and count each base once. <code>overlapping_reference_bp</code> makes "
        "duplicated alignment content visible so it does not silently inflate raw "
        "metrics.</p>\n"
    )
    w('<dl class="defs">\n')
    for name, definition in METRIC_DEFINITIONS:
        direction = FLAG_DIRECTIONS.get(name)
        extra = ""
        if direction:
            extra = f" <em>Flagged when unusually {'low' if direction == 'low' else 'high'}.</em>"
        # Definitions are static trusted HTML; names are escaped.
        w(f"<dt>{esc(name)}</dt><dd>{definition}{extra}</dd>\n")
    w("</dl>\n")
    w("</body>\n</html>\n")
    return "".join(out)


def _sample_anchor(samples: list[str], sample: str) -> str:
    # Index-based anchors are unique and safe for any sample name.
    return f"s-{samples.index(sample)}"


# ── driver ───────────────────────────────────────────────────────────────────


def load_inputs(
    summary_paths: list[Path],
    by_reference_paths: list[Path],
    by_query_paths: list[Path],
) -> tuple[list[dict[str, str]], dict[str, list[dict[str, str]]], dict[str, list[dict[str, str]]]]:
    if not summary_paths:
        raise ValueError("No summary inputs given")
    summaries: list[dict[str, str]] = []
    seen: dict[str, Path] = {}
    for path in summary_paths:
        row = read_summary(path)
        sample = row["sample"]
        if sample in seen:
            raise ValueError(
                f"Duplicate sample name {sample!r} in {path} (also in {seen[sample]})"
            )
        seen[sample] = path
        summaries.append(row)
    samples = [s["sample"] for s in summaries]
    by_reference = read_contig_tables(
        by_reference_paths, BY_REFERENCE_COLUMNS, samples, "by-reference-contig"
    )
    by_query = read_contig_tables(by_query_paths, BY_QUERY_COLUMNS, samples, "by-query-contig")
    return summaries, by_reference, by_query


def build_report(
    summary_paths: list[Path],
    by_reference_paths: list[Path],
    by_query_paths: list[Path],
    out_tsv: Path,
    out_html: Path,
    *,
    min_samples: int = DEFAULT_MIN_SAMPLES,
    z_threshold: float = DEFAULT_Z_THRESHOLD,
    thresholds: dict[str, float] | None = None,
) -> FlagResult:
    thresholds = dict(thresholds or {})
    summaries, by_reference, by_query = load_inputs(
        summary_paths, by_reference_paths, by_query_paths
    )
    flag_result = compute_flags(
        summaries, min_samples=min_samples, z_threshold=z_threshold, thresholds=thresholds
    )
    write_summary_tsv(out_tsv, summaries, flag_result.flags)
    page = render_html(
        summaries,
        flag_result,
        by_reference,
        by_query,
        min_samples=min_samples,
        z_threshold=z_threshold,
        thresholds=thresholds,
    )
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(page, encoding="utf-8")
    return flag_result


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Build the cross-sample MAF QC TSV and HTML report.")
    ap.add_argument("--summaries", nargs="+", required=True,
                    help="Per-sample <sample>.maf_stats.tsv files.")
    ap.add_argument("--by-reference", nargs="*", default=[],
                    help="Per-sample <sample>.by_reference_contig.tsv files.")
    ap.add_argument("--by-query", nargs="*", default=[],
                    help="Per-sample <sample>.by_query_contig.tsv files.")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--min-samples", type=int, default=DEFAULT_MIN_SAMPLES,
                    help="Minimum samples with a value before relative flags are computed.")
    ap.add_argument("--z-threshold", type=float, default=DEFAULT_Z_THRESHOLD)
    ap.add_argument("--flag-min-reference-coverage", type=float, default=None, metavar="PCT")
    ap.add_argument("--flag-min-identity", type=float, default=None, metavar="PCT")
    ap.add_argument("--flag-max-gap-fraction", type=float, default=None, metavar="PCT")
    ap.add_argument("--flag-max-breakpoints-per-gb", type=float, default=None, metavar="X")
    args = ap.parse_args(argv)
    if args.min_samples < 1:
        ap.error("--min-samples must be >= 1")
    if not args.z_threshold > 0:
        ap.error("--z-threshold must be > 0")
    return args


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    thresholds = {
        metric: getattr(args, dest)
        for dest, metric in THRESHOLD_OPTIONS.items()
        if getattr(args, dest) is not None
    }
    try:
        build_report(
            [Path(p) for p in args.summaries],
            [Path(p) for p in args.by_reference],
            [Path(p) for p in args.by_query],
            Path(args.out_tsv),
            Path(args.out_html),
            min_samples=args.min_samples,
            z_threshold=args.z_threshold,
            thresholds=thresholds,
        )
    except ValueError as exc:
        raise SystemExit(f"maf_stats_report: error: {exc}") from exc


if __name__ == "__main__":
    main()
