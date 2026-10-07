import re
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.maf_stats_report import (  # noqa: E402
    BREAKPOINT_COLUMNS,
    BY_QUERY_COLUMNS,
    BY_REFERENCE_COLUMNS,
    FLAG_DIRECTIONS,
    NESTED_COLUMNS,
    NESTED_EXTRA_COLUMNS,
    OUTPUT_BREAKPOINT_COLUMNS,
    OUTPUT_NESTED_COLUMNS,
    OUTPUT_SUMMARY_COLUMNS,
    SUMMARY_COLUMNS,
    annotate_nested_recurrence,
    annotate_recurrence,
    build_report,
    compute_flags,
    format_flags,
    load_inputs,
    recurrence_counts,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "maf_stats_report.py"

ASSEMBLY_NA = {
    c: "NA"
    for c in SUMMARY_COLUMNS
    if c.startswith("assembly_") or c.startswith("unaligned_")
}


def summary_row(sample: str, *, assembly: bool = False, **overrides) -> dict[str, str]:
    row = {
        "sample": sample,
        "reference_length_bp": "1000000000",
        "aligned_reference_bp": "360000000",
        "aligned_reference_pct": "36.00",
        "block_span_reference_bp": "940000000",
        "block_span_reference_pct": "94.00",
        "query_length_bp": "1100000000",
        "query_length_source": "query_fai",
        "aligned_query_bp": "350000000",
        "aligned_query_pct": "31.82",
        "block_span_query_bp": "900000000",
        "block_span_query_pct": "81.82",
        "identity_matches": "340000000",
        "identity_compared_columns": "350000000",
        "identity_pct": "97.14",
        "alignment_columns": "1800000000",
        "insertion_columns": "700000000",
        "insertion_column_pct": "38.89",
        "deletion_columns": "740000000",
        "deletion_column_pct": "41.11",
        "query_n_bases_in_blocks": "1000",
        "blocks": "20000",
        "overlapping_reference_bp": "1000000",
        "overlapping_query_bp": "800000",
        "nested_blocks": "4",
        "nested_block_reference_bp": "200000",
        "aligned_query_contigs": "12",
        "aligned_query_contig_n50_bp": "90000000",
        "strand_flips": "10",
        "reference_contig_jumps": "5",
        "out_of_order_adjacencies": "3",
        "breakpoint_adjacencies": "15",
        "breakpoints_near_contig_end": "4",
        "breakpoints_near_n_gap": "NA",
        "breakpoints_interior": "11",
        "breakpoints_per_gb_aligned": "41.67",
        "min_block_bp": "0",
        "overlap_tolerance_bp": "0",
        "breakpoint_context_bp": "1000000",
    }
    row.update(ASSEMBLY_NA)
    if assembly:
        row.update(
            {
                "breakpoints_near_n_gap": "2",
                "assembly_length_bp": "1100000000",
                "assembly_sequences": "40",
                "assembly_scaffold_n50_bp": "120000000",
                "assembly_scaffold_l50": "5",
                "assembly_largest_sequence_bp": "150000000",
                "assembly_contig_pieces": "65",
                "assembly_contig_n50_bp": "80000000",
                "assembly_contig_l50": "6",
                "assembly_n_bp": "1100000",
                "assembly_n_pct": "0.10",
                "assembly_n_gaps": "25",
                "assembly_gc_pct": "46.50",
                "assembly_softmasked_pct": "70.00",
                "assembly_major_sequences": "10",
                "assembly_major_telomeric_ends": "14",
                "assembly_minor_sequences_with_telomere": "2",
                "unaligned_contigs": "28",
                "unaligned_contig_bp": "5000000",
                "unaligned_contig_n_pct": "1.00",
                "unaligned_contig_softmasked_pct": "90.00",
            }
        )
    row.update({k: str(v) for k, v in overrides.items()})
    assert set(row) == set(SUMMARY_COLUMNS)
    return row


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> Path:
    lines = ["\t".join(columns)]
    lines += ["\t".join(str(r.get(c, "")) for c in columns) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_summary(tmp_path: Path, sample_row: dict[str, str], name: str | None = None) -> Path:
    safe = re.sub(r"[^\w]", "_", sample_row["sample"])
    fname = name or f"{safe}.maf_stats.tsv"
    return write_tsv(tmp_path / fname, SUMMARY_COLUMNS, [sample_row])


def ref_row(sample: str, contig: str, breakpoints: int = 0, dotplot: str = "") -> dict[str, str]:
    return {
        "sample": sample,
        "reference_contig": contig,
        "reference_length_bp": "1000000",
        "aligned_reference_bp": "360000",
        "aligned_reference_pct": "36.00",
        "block_span_reference_bp": "940000",
        "block_span_reference_pct": "94.00",
        "identity_matches": "340000",
        "identity_compared_columns": "350000",
        "identity_pct": "97.14",
        "alignment_columns": "1800000",
        "insertion_columns": "700000",
        "insertion_column_pct": "38.89",
        "deletion_columns": "740000",
        "deletion_column_pct": "41.11",
        "blocks": "20",
        "query_contigs": "2",
        "overlapping_reference_bp": "0",
        "breakpoint_adjacencies": str(breakpoints),
        "dotplot": dotplot,
    }


def query_row(sample: str, contig: str, breakpoints: int = 0, aligned: bool = True) -> dict[str, str]:
    return {
        "sample": sample,
        "query_contig": contig,
        "query_length_bp": "1000000",
        "query_length_source": "query_fai",
        "aligned": "true" if aligned else "false",
        "aligned_query_bp": "300000" if aligned else "0",
        "block_span_query_bp": "800000" if aligned else "0",
        "block_span_query_pct": "80.00" if aligned else "0.00",
        "overlapping_query_bp": "0",
        "unaligned_start_bp": "1000",
        "unaligned_end_bp": "2000",
        "blocks": "20",
        "blocks_considered": "20",
        "nested_blocks": "0",
        "reference_contigs": "1",
        "strand_flips": str(breakpoints),
        "reference_contig_jumps": "0",
        "out_of_order_adjacencies": "0",
        "breakpoint_adjacencies": str(breakpoints),
        "n_bp": "NA",
        "n_gaps": "NA",
        "contig_pieces": "NA",
        "gc_pct": "NA",
        "softmasked_pct": "NA",
        "major": "NA",
        "telomere_start": "NA",
        "telomere_end": "NA",
    }


def bp_row(
    sample: str,
    left: tuple[str, int],
    right: tuple[str, int],
    *,
    query_contig: str = "q1",
    location: str = "interior",
    flip: bool = True,
    jump: bool | None = None,
    ooo: bool = False,
) -> dict[str, str]:
    if jump is None:
        jump = left[0] != right[0]
    tf = lambda b: "true" if b else "false"  # noqa: E731
    return {
        "sample": sample,
        "query_contig": query_contig,
        "query_contig_length_bp": "50000000",
        "query_junction_start": "1000000",
        "query_junction_end": "1000500",
        "distance_to_contig_end_bp": "1000000",
        "location": location,
        "near_contig_end": tf(location == "contig_end"),
        "near_n_gap": "NA",
        "strand_flip": tf(flip),
        "reference_contig_jump": tf(jump),
        "out_of_order": tf(ooo),
        "left_reference_contig": left[0],
        "left_reference_pos": str(left[1]),
        "left_strand": "+",
        "right_reference_contig": right[0],
        "right_reference_pos": str(right[1]),
        "right_strand": "-" if flip else "+",
    }


def nested_row(
    sample: str,
    ref: tuple[str, int, int],
    *,
    query_contig: str = "q1",
    container: tuple[str, int, int] = ("chr9", 50_000_000, 60_000_000),
    strand: str = "-",
) -> dict[str, str]:
    return {
        "sample": sample,
        "query_contig": query_contig,
        "query_start": "2000000",
        "query_end": "2050000",
        "strand": strand,
        "reference_contig": ref[0],
        "reference_start": str(ref[1]),
        "reference_end": str(ref[2]),
        "container_query_start": "1000000",
        "container_query_end": "9000000",
        "container_reference_contig": container[0],
        "container_reference_start": str(container[1]),
        "container_reference_end": str(container[2]),
        "container_strand": "+",
    }


def cohort(n: int, **per_sample_overrides) -> list[dict[str, str]]:
    """n samples with slightly varied metrics; overrides keyed by 's<index>'."""
    rows = []
    for i in range(n):
        rows.append(
            summary_row(
                f"S{i}",
                aligned_reference_pct=f"{36 + (i % 3) * 0.5:.2f}",
                identity_pct=f"{97 + (i % 3) * 0.2:.2f}",
                aligned_query_pct=f"{32 + (i % 3) * 0.5:.2f}",
                reference_contig_jumps=str(5 + i % 2),
            )
        )
    for idx, overrides in per_sample_overrides.items():
        rows[int(idx.lstrip("s"))].update({k: str(v) for k, v in overrides.items()})
    return rows


def read_out_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").rstrip("\n").split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def run_report(tmp_path, rows, ref=None, query=None, breakpoints=None, nested=None, **kwargs):
    paths = [write_summary(tmp_path, r) for r in rows]
    out_tsv = tmp_path / "out" / "maf_stats.tsv"
    out_html = tmp_path / "out" / "maf_stats.html"
    out_bp = tmp_path / "out" / "maf_breakpoints.tsv"
    out_nested = tmp_path / "out" / "maf_nested_blocks.tsv"
    result = build_report(
        paths, ref or [], query or [], out_tsv, out_html,
        out_breakpoints_tsv=out_bp, out_nested_tsv=out_nested,
        breakpoint_paths=breakpoints or [], nested_paths=nested or [], **kwargs,
    )
    return result, read_out_tsv(out_tsv), out_html.read_text(encoding="utf-8")


def write_breakpoints(tmp_path: Path, sample: str, rows: list[dict[str, str]]) -> Path:
    return write_tsv(tmp_path / f"{sample}.breakpoints.tsv", BREAKPOINT_COLUMNS, rows)


def write_nested(tmp_path: Path, sample: str, rows: list[dict[str, str]]) -> Path:
    return write_tsv(tmp_path / f"{sample}.nested_blocks.tsv", NESTED_COLUMNS, rows)


# ── aggregation ─────────────────────────────────────────────────────────────


def test_one_row_per_sample_in_input_order(tmp_path):
    rows = [summary_row(s) for s in ["zeta", "alpha", "mid"]]
    _, out, _ = run_report(tmp_path, rows)
    assert [r["sample"] for r in out] == ["zeta", "alpha", "mid"]
    header = (tmp_path / "out" / "maf_stats.tsv").read_text().split("\n")[0].split("\t")
    assert header == SUMMARY_COLUMNS + [
        "breakpoints_recurrent", "breakpoints_private", "nested_recurrent", "nested_private",
        "flags",
    ]
    assert header == OUTPUT_SUMMARY_COLUMNS + ["flags"]
    # Values pass through unchanged.
    assert out[1]["aligned_reference_bp"] == "360000000"
    assert out[1]["assembly_length_bp"] == "NA"
    assert all(r["flags"] == "" for r in out)
    # No breakpoint files: every sample has zero breakpoints.
    assert all(r["breakpoints_recurrent"] == "0" and r["breakpoints_private"] == "0" for r in out)
    bp_out = (tmp_path / "out" / "maf_breakpoints.tsv").read_text().split("\n")
    assert bp_out[0].split("\t") == OUTPUT_BREAKPOINT_COLUMNS
    assert bp_out[1:] == [""]
    nested_out = (tmp_path / "out" / "maf_nested_blocks.tsv").read_text().split("\n")
    assert nested_out[0].split("\t") == OUTPUT_NESTED_COLUMNS
    assert nested_out[1:] == [""]
    assert all(r["nested_recurrent"] == "0" and r["nested_private"] == "0" for r in out)


# ── relative flags ───────────────────────────────────────────────────────────


def test_low_is_bad_outlier_flagged_and_high_outlier_not(tmp_path):
    rows = cohort(8, s2={"aligned_reference_pct": "10.00"}, s5={"aligned_reference_pct": "90.00"})
    result, out, _ = run_report(tmp_path, rows)
    assert "aligned_reference_pct" in result.flags[2]
    assert out[2]["flags"].startswith("aligned_reference_pct:robust_z=-")
    # A high aligned-reference outlier is good, not flagged.
    assert "aligned_reference_pct" not in result.flags[5]
    assert "aligned_reference_pct" not in out[5]["flags"]


def test_high_is_bad_outlier_flagged_and_low_not():
    rows = cohort(
        8,
        s1={"overlapping_reference_bp": "50000000", "reference_contig_jumps": "40",
            "breakpoints_per_gb_aligned": "400"},
        s4={"overlapping_reference_bp": "0", "reference_contig_jumps": "0",
            "breakpoints_per_gb_aligned": "0"},
    )
    result = compute_flags(rows)
    for metric in ("overlapping_reference_bp", "reference_contig_jumps", "breakpoints_per_gb_aligned"):
        assert metric in result.flags[1], metric
        assert result.flags[1][metric][0].startswith("robust_z=")
        assert metric not in result.flags[4], metric
    z = float(result.flags[1]["reference_contig_jumps"][0].split("=")[1].split("(")[0])
    assert z > 3.5


def test_flag_directions_match_spec():
    assert FLAG_DIRECTIONS == {
        "aligned_reference_pct": "low",
        "aligned_query_pct": "low",
        "identity_pct": "low",
        "assembly_contig_n50_bp": "low",
        "overlapping_reference_bp": "high",
        "overlapping_query_bp": "high",
        "breakpoints_per_gb_aligned": "high",
        "reference_contig_jumps": "high",
        "breakpoints_private": "high",
        "assembly_n_pct": "high",
    }
    assert "nested_blocks" not in FLAG_DIRECTIONS
    assert "nested_private" not in FLAG_DIRECTIONS


def test_unflagged_metrics_are_not_flagged():
    # Huge outliers in metrics that are no longer flagged.
    rows = cohort(8, s0={"blocks": "9999999", "block_span_reference_pct": "1.00",
                         "insertion_column_pct": "99.00"})
    assert compute_flags(rows).flags[0] == {}


def test_min_samples_guard_skips_relative_flags(tmp_path):
    rows = cohort(4, s0={"aligned_reference_pct": "5.00"})
    result, out, page = run_report(tmp_path, rows)
    assert all(not f for f in result.flags)
    assert all(r["flags"] == "" for r in out)
    assert "Relative (cohort) flagging skipped" in page
    # Same data with a lower guard is flagged.
    assert "aligned_reference_pct" in compute_flags(rows, min_samples=3).flags[0]


def test_mad_zero_uses_floor():
    # Constant cohort except one trivially different value: not flagged.
    rows = [summary_row(f"S{i}", identity_pct="99.00") for i in range(6)]
    rows[3]["identity_pct"] = "98.50"
    result = compute_flags(rows)
    assert result.stats["identity_pct"].mad == 0
    assert result.stats["identity_pct"].scale == pytest.approx(1.0)
    assert "identity_pct" not in result.flags[3]

    # Hugely different value: flagged, with the floor made visible.
    rows[3]["identity_pct"] = "90.00"
    result = compute_flags(rows)
    reasons = result.flags[3]["identity_pct"]
    assert reasons == ["robust_z=-9.00(mad0_floor)"]
    assert format_flags(result.flags[3]) == "identity_pct:robust_z=-9.00(mad0_floor)"


def test_relative_floor_for_count_metrics():
    # breakpoints_per_gb MAD is 0; scale = max(1, 0.05*200=10) -> 10.
    rows = [summary_row(f"S{i}", breakpoints_per_gb_aligned="200") for i in range(6)]
    rows[0]["breakpoints_per_gb_aligned"] = "230"  # z = 3.0, below threshold
    rows[1]["breakpoints_per_gb_aligned"] = "240"  # z = 4.0, flagged
    result = compute_flags(rows)
    assert result.stats["breakpoints_per_gb_aligned"].scale == pytest.approx(10.0)
    assert "breakpoints_per_gb_aligned" not in result.flags[0]
    assert result.flags[1]["breakpoints_per_gb_aligned"] == ["robust_z=4.00(mad0_floor)"]


def test_assembly_n50_floor_and_direction():
    rows = [summary_row(f"S{i}", assembly=True) for i in range(6)]
    # median 80 Mb, MAD 0 -> scale = max(1e6, 0.05*80e6 = 4e6) = 4e6.
    rows[2]["assembly_contig_n50_bp"] = "60000000"  # z = -5
    rows[3]["assembly_contig_n50_bp"] = "200000000"  # high is good
    result = compute_flags(rows)
    assert result.stats["assembly_contig_n50_bp"].scale == pytest.approx(4_000_000)
    assert result.flags[2]["assembly_contig_n50_bp"] == ["robust_z=-5.00(mad0_floor)"]
    assert "assembly_contig_n50_bp" not in result.flags[3]


def test_breakpoints_private_flagged_high():
    rows = cohort(6)
    for i, r in enumerate(rows):
        r["breakpoints_private"] = str(3 + i % 2)
        r["breakpoints_recurrent"] = "5"
    rows[4]["breakpoints_private"] = "40"
    rows[5]["breakpoints_private"] = "0"  # low is fine
    result = compute_flags(rows)
    assert "breakpoints_private" in result.flags[4]
    assert "breakpoints_private" not in result.flags[5]


# ── absolute thresholds ─────────────────────────────────────────────────────


def test_absolute_thresholds(tmp_path):
    rows = [
        summary_row("A", identity_pct="94.00", aligned_reference_pct="20.00"),
        summary_row("B", identity_pct="95.00", aligned_reference_pct="30.00"),
        summary_row("C", identity_pct="99.00", breakpoints_per_gb_aligned="250"),
    ]
    thresholds = {
        "aligned_reference_pct": 30.0,
        "identity_pct": 95.0,
        "breakpoints_per_gb_aligned": 100.0,
    }
    result, out, page = run_report(tmp_path, rows, thresholds=thresholds)
    assert out[0]["flags"] == (
        "aligned_reference_pct:below_threshold_30;identity_pct:below_threshold_95"
    )
    assert out[1]["flags"] == ""  # equal to threshold is not flagged
    assert out[2]["flags"] == "breakpoints_per_gb_aligned:above_threshold_100"
    assert 'title="identity_pct:below_threshold_95"' in page
    # n=3 < 5: relative flagging skipped but absolute thresholds still applied.
    assert "Relative (cohort) flagging skipped" in page


def test_threshold_for_unknown_metric_errors():
    with pytest.raises(ValueError, match="No flag direction"):
        compute_flags([summary_row("A")], thresholds={"blocks": 1.0})


# ── NA handling ──────────────────────────────────────────────────────────────


def test_na_values_are_ignored(tmp_path):
    rows = cohort(6)
    rows[0]["identity_pct"] = "NA"
    rows[0]["breakpoints_per_gb_aligned"] = "NA"
    rows[1]["breakpoints_per_gb_aligned"] = "NA"
    result, out, page = run_report(
        tmp_path, rows, thresholds={"identity_pct": 99.9, "breakpoints_per_gb_aligned": 1.0}
    )
    assert "identity_pct" not in result.flags[0]
    assert result.stats["identity_pct"].n == 5
    # breakpoints has only 4 non-NA values -> relative flagging skipped for it.
    assert result.skipped["breakpoints_per_gb_aligned"] == 4
    assert out[0]["identity_pct"] == "NA"
    assert '<span class="na">NA</span>' in page
    assert "Relative flagging skipped</strong> for metrics" in page


def test_all_na_metric_does_not_crash(tmp_path):
    rows = cohort(5)
    for r in rows:
        r["aligned_query_pct"] = "NA"
    result, _, page = run_report(tmp_path, rows)
    assert result.stats["aligned_query_pct"] is None
    assert result.skipped["aligned_query_pct"] == 0
    assert "Not available for any sample" in page


def test_strip_plots_only_for_metrics_with_values(tmp_path):
    # No assembly data; 5 samples -> breakpoints_private is 0 (not NA).
    rows = cohort(5)
    _, _, page = run_report(tmp_path, rows)
    with_values = [m for m in FLAG_DIRECTIONS if not m.startswith("assembly_")]
    assert page.count('<svg class="strip"') == len(with_values)
    assert 'aria-label="Contig N50 (bp)"' not in page
    # One sample -> breakpoints_private NA -> no plot for it either.
    one = tmp_path / "one"
    one.mkdir()
    _, out, page = run_report(one, [summary_row("A")])
    assert out[0]["breakpoints_private"] == "NA"
    assert page.count('<svg class="strip"') == len(with_values) - 1


# ── assembly ─────────────────────────────────────────────────────────────────


def test_assembly_section_only_with_assembly_data(tmp_path):
    a = tmp_path / "a"
    a.mkdir()
    _, _, page = run_report(a, cohort(3))
    assert '<table class="assembly">' not in page
    assert "No sample has assembly statistics" in page
    # MAF-derived contiguity is always shown, with its caveat.
    assert '<table class="maf-contiguity">' in page
    assert "From the MAF: aligned contigs only, a lower bound on fragmentation" in page

    b = tmp_path / "b"
    b.mkdir()
    rows = [summary_row("A", assembly=True), summary_row("B")]
    _, out, page = run_report(b, rows)
    assert '<table class="assembly">' in page
    table = page[page.index('<table class="assembly">'):]
    table = table[: table.index("</table>")]
    assert "80,000,000" in table  # A's contig N50
    assert "120,000,000" in table  # A's scaffold N50
    assert "<td>14 / 20</td>" in table  # major telomeric ends / 2 x major sequences
    header = re.findall(r"<th>([^<]*)</th>", table[: table.index("</tr>")])
    assert header == [
        "Sample", "Assembly length (bp)", "Sequences", "Scaffold N50 (bp)", "Contig pieces",
        "Contig N50 (bp)", "N gaps", "N (%)", "GC (%)", "Soft-masked (%)",
        "Telomeric ends (major)", "Minor seqs with telomere", "Unaligned contigs",
        "Unaligned contig bp",
    ]
    b_row = table[table.index(">B</a>"):]
    assert '<span class="na">NA</span>' in b_row
    assert out[1]["assembly_contig_n50_bp"] == "NA"


def test_assembly_na_for_some_samples_does_not_crash_or_flag():
    rows = [summary_row(f"S{i}", assembly=i < 5) for i in range(8)]
    rows[0]["assembly_n_pct"] = "8.00"  # median 0.10, floor 0.5 -> z ~ 15.8
    result = compute_flags(rows)
    assert result.stats["assembly_n_pct"].n == 5
    assert "assembly_n_pct" in result.flags[0]
    for i in range(5, 8):
        assert not any(m.startswith("assembly_") for m in result.flags[i])
    # With only 3 assembly samples relative flagging is skipped for those metrics.
    result = compute_flags(rows[3:])
    assert result.skipped["assembly_n_pct"] == 2
    assert result.flags == [{} for _ in rows[3:]]


# ── recurrence ───────────────────────────────────────────────────────────────


def test_recurrence_same_orientation_within_window():
    bps = {
        "A": [bp_row("A", ("chr1", 1_000_000), ("chr1", 5_000_000))],
        "B": [bp_row("B", ("chr1", 1_400_000), ("chr1", 4_700_000))],
    }
    out = annotate_recurrence(["A", "B"], bps, 500_000)
    assert [r["recurrence_samples"] for r in out] == ["1", "1"]
    assert [r["recurrent"] for r in out] == ["true", "true"]


def test_recurrence_either_orientation():
    bps = {
        "A": [bp_row("A", ("chr1", 1_000_000), ("chr3", 7_000_000))],
        "B": [bp_row("B", ("chr3", 7_100_000), ("chr1", 900_000))],
        "C": [bp_row("C", ("chr1", 1_050_000), ("chr3", 6_950_000))],
    }
    out = annotate_recurrence(["A", "B", "C"], bps, 500_000)
    assert [r["recurrence_samples"] for r in out] == ["2", "2", "2"]


def test_recurrence_outside_window_is_private():
    bps = {
        "A": [bp_row("A", ("chr1", 1_000_000), ("chr1", 5_000_000))],
        # left end matches, right end 600 kb away
        "B": [bp_row("B", ("chr1", 1_000_000), ("chr1", 5_600_000))],
    }
    out = annotate_recurrence(["A", "B"], bps, 500_000)
    assert [r["recurrent"] for r in out] == ["false", "false"]
    assert [r["recurrence_samples"] for r in out] == ["0", "0"]
    # Exactly at the window boundary matches.
    out = annotate_recurrence(["A", "B"], bps, 600_000)
    assert [r["recurrent"] for r in out] == ["true", "true"]


def test_recurrence_different_contig_pair_no_match():
    bps = {
        "A": [bp_row("A", ("chr1", 1_000_000), ("chr2", 5_000_000))],
        "B": [bp_row("B", ("chr1", 1_000_000), ("chr3", 5_000_000))],
        # same positions but contigs swapped between ends: not a match either
        "C": [bp_row("C", ("chr2", 1_000_000), ("chr1", 5_000_000))],
    }
    out = annotate_recurrence(["A", "B", "C"], bps, 500_000)
    assert [r["recurrence_samples"] for r in out] == ["0", "0", "0"]


def test_recurrence_counts_distinct_other_samples_not_rows():
    bps = {
        "A": [bp_row("A", ("chr1", 100), ("chr1", 900_000))],
        "B": [
            bp_row("B", ("chr1", 200), ("chr1", 900_100)),
            bp_row("B", ("chr1", 300), ("chr1", 900_200)),
        ],
        # same-sample matches do not count
        "C": [
            bp_row("C", ("chr5", 100), ("chr6", 100)),
            bp_row("C", ("chr5", 150), ("chr6", 150)),
        ],
    }
    out = annotate_recurrence(["A", "B", "C"], bps, 1000)
    assert [r["recurrence_samples"] for r in out] == ["1", "1", "1", "0", "0"]
    counts = recurrence_counts(["A", "B", "C"], out)
    assert counts == {
        "A": {"breakpoints_recurrent": "1", "breakpoints_private": "0"},
        "B": {"breakpoints_recurrent": "2", "breakpoints_private": "0"},
        "C": {"breakpoints_recurrent": "0", "breakpoints_private": "2"},
    }


def test_recurrence_single_sample_is_na():
    bps = {"A": [bp_row("A", ("chr1", 100), ("chr1", 900_000))]}
    out = annotate_recurrence(["A"], bps, 500_000)
    assert out[0]["recurrence_samples"] == "NA" and out[0]["recurrent"] == "NA"
    assert recurrence_counts(["A"], out) == {
        "A": {"breakpoints_recurrent": "NA", "breakpoints_private": "NA"}
    }


def test_recurrence_na_position_never_matches():
    a = bp_row("A", ("chr1", 100), ("chr1", 900_000))
    b = bp_row("B", ("chr1", 100), ("chr1", 900_000))
    b["right_reference_pos"] = "NA"
    out = annotate_recurrence(["A", "B"], {"A": [a], "B": [b]}, 500_000)
    assert [r["recurrent"] for r in out] == ["false", "false"]


def test_combined_breakpoints_tsv_and_summary_counts(tmp_path):
    rows = cohort(3)
    files = [
        write_breakpoints(tmp_path, "S0", [
            bp_row("S0", ("chr1", 1_000_000), ("chr1", 5_000_000), query_contig="a"),
            bp_row("S0", ("chr2", 10), ("chr4", 10), query_contig="b", location="contig_end"),
        ]),
        write_breakpoints(tmp_path, "S1", [
            bp_row("S1", ("chr1", 5_100_000), ("chr1", 1_100_000), query_contig="c"),
        ]),
        write_breakpoints(tmp_path, "S2", []),  # header only
    ]
    # Files given out of sample order; output follows summary order.
    _, out, page = run_report(tmp_path, rows, breakpoints=list(reversed(files)))
    bp_out = read_out_tsv(tmp_path / "out" / "maf_breakpoints.tsv")
    header = (tmp_path / "out" / "maf_breakpoints.tsv").read_text().split("\n")[0].split("\t")
    assert header == BREAKPOINT_COLUMNS + ["recurrence_samples", "recurrent"]
    assert [(r["sample"], r["query_contig"]) for r in bp_out] == [("S0", "a"), ("S0", "b"), ("S1", "c")]
    assert [r["recurrence_samples"] for r in bp_out] == ["1", "0", "1"]
    assert [r["recurrent"] for r in bp_out] == ["true", "false", "true"]
    assert bp_out[1]["location"] == "contig_end"
    assert bp_out[1]["left_reference_pos"] == "10"

    assert [(r["breakpoints_recurrent"], r["breakpoints_private"]) for r in out] == [
        ("1", "1"), ("1", "0"), ("0", "0"),
    ]
    # Breakpoints section: explanations and the cross-sample table, recurrent first.
    assert "<h2>Breakpoints</h2>" in page
    assert "scaffold gap" in page and "cannot be distinguished without reads" in page
    table = page[page.index('<table class="bp-table">'):]
    table = table[: table.index("</table>")]
    assert table.index("chr1:1,000,000") < table.index("chr1:5,100,000") < table.index("chr2:10")
    assert "flip, jump" in table  # chr2->chr4 row
    assert '<tr class="recurrent">' in table
    # Summary breakpoint_adjacencies (15) do not match the rows: note is shown.
    assert "Breakpoint table does not match summary" in page


def test_single_sample_summary_recurrence_na(tmp_path):
    bp = write_breakpoints(tmp_path, "A", [bp_row("A", ("chr1", 1), ("chr1", 2))])
    _, out, page = run_report(tmp_path, [summary_row("A")], breakpoints=[bp])
    assert out[0]["breakpoints_recurrent"] == "NA"
    assert out[0]["breakpoints_private"] == "NA"
    bp_out = read_out_tsv(tmp_path / "out" / "maf_breakpoints.tsv")
    assert bp_out[0]["recurrence_samples"] == "NA" and bp_out[0]["recurrent"] == "NA"
    assert "recurrence not computed" in page


def test_breakpoint_tables_capped(tmp_path):
    rows = cohort(5)
    many = [bp_row("S0", ("chr1", i * 10_000_000), ("chr2", i * 10_000_000), query_contig=f"q{i}")
            for i in range(250)]
    bp = write_breakpoints(tmp_path, "S0", many)
    _, _, page = run_report(tmp_path, rows, breakpoints=[bp])
    assert "Showing 200 of 250 breakpoints; the combined breakpoints TSV has every row." in page
    assert "Showing 50 of 250 breakpoints" in page
    assert len(read_out_tsv(tmp_path / "out" / "maf_breakpoints.tsv")) == 250


# ── input errors ────────────────────────────────────────────────────────────


def test_duplicate_sample_error(tmp_path):
    a = write_summary(tmp_path, summary_row("dup"), name="a.maf_stats.tsv")
    b = write_summary(tmp_path, summary_row("dup"), name="b.maf_stats.tsv")
    with pytest.raises(ValueError, match="Duplicate sample name 'dup'"):
        load_inputs([a, b], [], [])


def test_multi_row_summary_error(tmp_path):
    path = write_tsv(tmp_path / "x.maf_stats.tsv", SUMMARY_COLUMNS, [summary_row("A"), summary_row("B")])
    with pytest.raises(ValueError, match="exactly 1 data row, found 2"):
        load_inputs([path], [], [])
    empty = write_tsv(tmp_path / "y.maf_stats.tsv", SUMMARY_COLUMNS, [])
    with pytest.raises(ValueError, match="found 0"):
        load_inputs([empty], [], [])


def test_missing_column_error(tmp_path):
    cols = [c for c in SUMMARY_COLUMNS if c != "identity_pct"]
    path = write_tsv(tmp_path / "x.maf_stats.tsv", cols, [summary_row("A")])
    with pytest.raises(ValueError, match="missing required column.*identity_pct"):
        load_inputs([path], [], [])


def test_old_schema_rejected(tmp_path):
    cols = [c for c in SUMMARY_COLUMNS if c != "aligned_reference_pct"] + ["reference_coverage_pct"]
    row = summary_row("A")
    row["reference_coverage_pct"] = "90"
    path = write_tsv(tmp_path / "x.maf_stats.tsv", cols, [row])
    with pytest.raises(ValueError, match="aligned_reference_pct"):
        load_inputs([path], [], [])


def test_by_file_with_unknown_sample_error(tmp_path):
    s = write_summary(tmp_path, summary_row("A"))
    ref = write_tsv(tmp_path / "B.by_reference_contig.tsv", BY_REFERENCE_COLUMNS, [ref_row("B", "chr1")])
    with pytest.raises(ValueError, match="'B'.*not among the summary inputs"):
        load_inputs([s], [ref], [])
    bp = write_breakpoints(tmp_path, "B", [bp_row("B", ("chr1", 1), ("chr1", 2))])
    with pytest.raises(ValueError, match="breakpoints row for sample 'B'"):
        load_inputs([s], [], [], [bp])


def test_duplicate_breakpoint_files_error(tmp_path):
    s = write_summary(tmp_path, summary_row("A"))
    one = write_tsv(tmp_path / "1.tsv", BREAKPOINT_COLUMNS, [bp_row("A", ("chr1", 1), ("chr1", 2))])
    two = write_tsv(tmp_path / "2.tsv", BREAKPOINT_COLUMNS, [bp_row("A", ("chr1", 1), ("chr1", 2))])
    with pytest.raises(ValueError, match="duplicate breakpoints input for sample 'A'"):
        load_inputs([s], [], [], [one, two])


def test_by_files_matched_by_sample_column_not_filename(tmp_path):
    s = write_summary(tmp_path, summary_row("A"))
    ref = write_tsv(tmp_path / "unrelated_name.tsv", BY_REFERENCE_COLUMNS, [ref_row("A", "chr1")])
    bp = write_tsv(tmp_path / "other.tsv", BREAKPOINT_COLUMNS, [bp_row("A", ("chr1", 1), ("chr1", 2))])
    nb = write_tsv(tmp_path / "n.tsv", NESTED_COLUMNS, [nested_row("A", ("chr1", 1, 2))])
    _, by_ref, _, bps, nested = load_inputs([s], [ref], [], [bp], [nb])
    assert [r["reference_contig"] for r in by_ref["A"]] == ["chr1"]
    assert len(bps["A"]) == 1
    assert len(nested["A"]) == 1


# ── HTML ─────────────────────────────────────────────────────────────────────


def test_html_escapes_sample_names(tmp_path):
    evil = "<b>&x"
    rows = cohort(5)
    rows[0]["sample"] = evil
    ref = write_tsv(
        tmp_path / "e.by_reference_contig.tsv",
        BY_REFERENCE_COLUMNS,
        [ref_row(evil, 'chr"<1>', breakpoints=2, dotplot='dotplots/x/chr"<1>.png')],
    )
    bp = write_breakpoints(
        tmp_path, "e", [bp_row(evil, ('chr"<1>', 5), ("chr2", 9), query_contig="<q&>")]
    )
    _, out, page = run_report(tmp_path, rows, ref=[ref], breakpoints=[bp])
    assert out[0]["sample"] == evil
    assert "<b>&x" not in page
    assert "&lt;b&gt;&amp;x" in page
    assert 'chr"<1>' not in page
    assert "<q&>" not in page and "&lt;q&amp;&gt;" in page
    assert "chr&quot;&lt;1&gt;:5" in page
    assert 'src="dotplots/x/chr&quot;&lt;1&gt;.png"' in page
    bp_out = read_out_tsv(tmp_path / "out" / "maf_breakpoints.tsv")
    assert bp_out[0]["sample"] == evil


def test_dotplot_links_and_separate_tables(tmp_path):
    rows = cohort(5)
    ref = write_tsv(
        tmp_path / "S0.by_reference_contig.tsv",
        BY_REFERENCE_COLUMNS,
        [
            ref_row("S0", "chr1", breakpoints=0, dotplot="dotplots/S0/chr1.png"),
            ref_row("S0", "chr2", breakpoints=0, dotplot=""),
            ref_row("S0", "chr3", breakpoints=7, dotplot="dotplots/S0/chr3.png"),
        ],
    )
    query = write_tsv(
        tmp_path / "S0.by_query_contig.tsv",
        BY_QUERY_COLUMNS,
        [query_row("S0", "scaf1"), query_row("S0", "scaf2", breakpoints=3),
         query_row("S0", "unal1", aligned=False)],
    )
    _, _, page = run_report(tmp_path, rows, ref=[ref], query=[query])

    assert '<a href="dotplots/S0/chr1.png"><img src="dotplots/S0/chr1.png" loading="lazy"' in page
    assert '<img src="dotplots/S0/chr3.png"' in page
    assert "dotplots/S0/chr2.png" not in page
    # Breakpoint contig first among thumbnails and in the reference table.
    assert page.index('src="dotplots/S0/chr3.png"') < page.index('src="dotplots/S0/chr1.png"')
    ref_table = page[page.index('<table class="ref-table">'):]
    ref_table = ref_table[: ref_table.index("</table>")]
    assert ref_table.index("chr3") < ref_table.index("chr1") < ref_table.index("chr2")
    assert "N50" not in ref_table
    assert "Aligned ref (%)" in ref_table

    # Separate reference, query and breakpoint sections with headings.
    assert '<section class="ref-contigs">' in page
    assert '<section class="query-contigs">' in page
    assert '<section class="breakpoints">' in page
    assert "<h3>Reference contigs</h3>" in page and "<h3>Query contigs</h3>" in page
    query_table = page[page.index('<table class="query-table">'):]
    query_table = query_table[: query_table.index("</table>")]
    assert "scaf1" in query_table and "scaf2" in query_table and "unal1" in query_table
    assert query_table.index("scaf2") < query_table.index("scaf1")
    assert "<td>false</td>" in query_table
    assert "chr1" not in query_table


def test_html_alignment_table_and_flag_highlighting(tmp_path):
    rows = cohort(8, s2={"aligned_reference_pct": "10.00"})
    _, _, page = run_report(tmp_path, rows)
    assert '<circle class="dot flagged"' in page
    assert "<title>S2: 10.00</title>" in page
    assert re.search(
        r'<td class="flag" title="aligned_reference_pct:robust_z=-[\d.]+">10\.00</td>', page
    )
    table = page[page.index('<table class="samples">'):]
    header = table[: table.index("</tr>")]
    labels = re.findall(r"<th>([^<]*)</th>", header)
    assert labels == [
        "Sample", "Aligned ref (%)", "Aligned query (%)", "Identity (%)", "Insertion cols (%)",
        "Deletion cols (%)", "Strand flips", "Out-of-order adj.", "Ref-contig jumps",
        "Breakpoint adj.", "Recurrent breakpoints", "Private breakpoints", "Nested blocks",
        "Recurrent nested", "Private nested", "Flags",
    ]
    assert "maf_srcsize" in page and "overstates" in page
    assert "block_span_reference_bp" in page and "AnchorWave" in page
    assert "TTTAGGG" in page


def test_cli_end_to_end(tmp_path):
    rows = cohort(6, s3={"identity_pct": "80.00"})
    rows[1]["query_length_source"] = "maf_srcsize"
    summaries = [write_summary(tmp_path, r) for r in rows]
    ref = write_tsv(
        tmp_path / "S3.by_reference_contig.tsv",
        BY_REFERENCE_COLUMNS,
        [ref_row("S3", "chr1", breakpoints=1, dotplot="dotplots/S3/chr1.png")],
    )
    query = write_tsv(tmp_path / "S3.by_query_contig.tsv", BY_QUERY_COLUMNS, [query_row("S3", "q1")])
    bps = [
        write_breakpoints(tmp_path, "S3", [bp_row("S3", ("chr1", 100), ("chr2", 100))]),
        write_breakpoints(tmp_path, "S4", [bp_row("S4", ("chr2", 900), ("chr1", 900))]),
    ]
    out_tsv = tmp_path / "results" / "maf_stats.tsv"
    out_html = tmp_path / "results" / "maf_stats.html"
    out_bp = tmp_path / "results" / "maf_breakpoints.tsv"
    out_nested = tmp_path / "results" / "maf_nested_blocks.tsv"
    nested = [
        write_nested(tmp_path, "S2", [nested_row("S2", ("chr3", 10_000, 20_000))]),
        write_nested(tmp_path, "S5", [nested_row("S5", ("chr3", 10_500, 20_500))]),
    ]
    proc = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--summaries", *map(str, summaries),
            "--by-reference", str(ref),
            "--by-query", str(query),
            "--breakpoints", *map(str, bps),
            "--nested-blocks", *map(str, nested),
            "--out-tsv", str(out_tsv),
            "--out-html", str(out_html),
            "--out-breakpoints-tsv", str(out_bp),
            "--out-nested-tsv", str(out_nested),
            "--recurrence-window-bp", "1000",
            "--min-samples", "5",
            "--z-threshold", "3.5",
            "--flag-min-aligned-reference", "10",
            "--flag-min-identity", "50",
            "--flag-max-breakpoints-per-gb", "1000",
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
    out = read_out_tsv(out_tsv)
    assert [r["sample"] for r in out] == [f"S{i}" for i in range(6)]
    assert out[3]["flags"].startswith("identity_pct:robust_z=-")
    assert out[3]["breakpoints_recurrent"] == "1"
    assert [r["recurrence_samples"] for r in read_out_tsv(out_bp)] == ["1", "1"]
    assert [r["recurrent"] for r in read_out_tsv(out_nested)] == ["true", "true"]
    assert out[2]["nested_recurrent"] == "1" and out[2]["nested_private"] == "0"
    assert out[0]["nested_recurrent"] == "0"
    page = out_html.read_text(encoding="utf-8")
    assert 'src="dotplots/S3/chr1.png"' in page
    assert "32.50 &dagger;" in page
    assert "--flag-min-aligned-reference" in page

    # Removed options are rejected.
    proc = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--summaries", str(summaries[0]),
            "--out-tsv", str(out_tsv), "--out-html", str(out_html),
            "--out-breakpoints-tsv", str(out_bp),
            "--out-nested-tsv", str(out_nested),
            "--flag-max-gap-fraction", "5",
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0

    # Error path: duplicate sample gives a clear non-zero exit.
    proc = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--summaries", str(summaries[0]), str(summaries[0]),
            "--out-tsv", str(out_tsv), "--out-html", str(out_html),
            "--out-breakpoints-tsv", str(out_bp), "--out-nested-tsv", str(out_nested),
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "Duplicate sample name" in proc.stderr


def test_html_contig_tables_are_capped_breakpoints_first(tmp_path):
    rows = cohort(5)
    contigs = [query_row("S0", f"scaffold_{i}") for i in range(120)]
    contigs.append(query_row("S0", "scaffold_jump", breakpoints=3))
    query = write_tsv(tmp_path / "S0.by_query_contig.tsv", BY_QUERY_COLUMNS, contigs)
    _, _, html = run_report(tmp_path, rows, query=[query])
    assert "Showing 50 of 121 contigs" in html
    assert "scaffold_jump" in html  # breakpoint contigs survive the cap
    assert "scaffold_119" not in html


# ── nested blocks ────────────────────────────────────────────────────────────


def test_nested_recurrence_match_within_window():
    nested = {
        "A": [nested_row("A", ("chr1", 1_000_000, 1_050_000))],
        "B": [nested_row("B", ("chr1", 1_400_000, 1_450_000))],
        "C": [nested_row("C", ("chr1", 1_000_000, 1_600_000))],  # end 550 kb off from A
    }
    out = annotate_nested_recurrence(["A", "B", "C"], nested, 500_000)
    assert [r["recurrence_samples"] for r in out] == ["1", "2", "1"]
    assert [r["recurrent"] for r in out] == ["true", "true", "true"]
    assert list(out[0]) == OUTPUT_NESTED_COLUMNS


def test_nested_recurrence_outside_window_is_private():
    nested = {
        "A": [nested_row("A", ("chr1", 1_000_000, 1_050_000))],
        "B": [nested_row("B", ("chr1", 1_600_000, 1_650_000))],
    }
    out = annotate_nested_recurrence(["A", "B"], nested, 500_000)
    assert [r["recurrent"] for r in out] == ["false", "false"]
    assert [r["recurrence_samples"] for r in out] == ["0", "0"]
    # Exactly at the window boundary matches.
    out = annotate_nested_recurrence(["A", "B"], nested, 600_000)
    assert [r["recurrent"] for r in out] == ["true", "true"]


def test_nested_recurrence_different_contig_or_same_sample_no_match():
    nested = {
        "A": [nested_row("A", ("chr1", 1_000_000, 1_050_000)),
              nested_row("A", ("chr1", 1_000_000, 1_050_000), query_contig="q2")],
        "B": [nested_row("B", ("chr2", 1_000_000, 1_050_000))],
    }
    out = annotate_nested_recurrence(["A", "B"], nested, 500_000)
    assert [r["recurrence_samples"] for r in out] == ["0", "0", "0"]
    counts = recurrence_counts(["A", "B"], out, NESTED_EXTRA_COLUMNS)
    assert counts == {
        "A": {"nested_recurrent": "0", "nested_private": "2"},
        "B": {"nested_recurrent": "0", "nested_private": "1"},
    }


def test_nested_recurrence_single_sample_na():
    nested = {"A": [nested_row("A", ("chr1", 1, 2))]}
    out = annotate_nested_recurrence(["A"], nested, 500_000)
    assert out[0]["recurrence_samples"] == "NA" and out[0]["recurrent"] == "NA"
    assert recurrence_counts(["A"], out, NESTED_EXTRA_COLUMNS) == {
        "A": {"nested_recurrent": "NA", "nested_private": "NA"}
    }


def test_nested_recurrence_na_position_never_matches():
    a = nested_row("A", ("chr1", 100, 200))
    b = nested_row("B", ("chr1", 100, 200))
    b["reference_end"] = "NA"
    out = annotate_nested_recurrence(["A", "B"], {"A": [a], "B": [b]}, 500_000)
    assert [r["recurrent"] for r in out] == ["false", "false"]


def test_combined_nested_tsv_summary_columns_and_html(tmp_path):
    rows = cohort(3)
    rows[0]["nested_blocks"] = "2"
    rows[1]["nested_blocks"] = "1"
    rows[2]["nested_blocks"] = "0"
    files = [
        write_nested(tmp_path, "S0", [
            nested_row("S0", ("chr1", 1_000_000, 1_050_000), query_contig="a"),
            nested_row("S0", ("chr2", 7_000_000, 7_010_000), query_contig="b",
                       container=("chr5", 1, 2)),
        ]),
        write_nested(tmp_path, "S1", [
            nested_row("S1", ("chr1", 1_100_000, 1_150_000), query_contig="c"),
        ]),
        write_nested(tmp_path, "S2", []),  # header only
    ]
    _, out, page = run_report(tmp_path, rows, nested=list(reversed(files)))
    path = tmp_path / "out" / "maf_nested_blocks.tsv"
    assert path.read_text().split("\n")[0].split("\t") == NESTED_COLUMNS + [
        "recurrence_samples", "recurrent"
    ]
    nested_out = read_out_tsv(path)
    assert [(r["sample"], r["query_contig"]) for r in nested_out] == [
        ("S0", "a"), ("S0", "b"), ("S1", "c")
    ]
    assert [r["recurrence_samples"] for r in nested_out] == ["1", "0", "1"]
    assert nested_out[1]["container_reference_contig"] == "chr5"
    assert [(r["nested_recurrent"], r["nested_private"]) for r in out] == [
        ("1", "1"), ("1", "0"), ("0", "0"),
    ]
    # Summary column order: breakpoint then nested recurrence columns before flags.
    header = (tmp_path / "out" / "maf_stats.tsv").read_text().split("\n")[0].split("\t")
    assert header[-5:] == [
        "breakpoints_recurrent", "breakpoints_private", "nested_recurrent", "nested_private",
        "flags",
    ]
    # HTML section: explanation, table with recurrent rows first, container shown.
    assert "<h2>Nested (secondary/transposed) alignments</h2>" in page
    assert "excluded from the breakpoint walk" in page
    assert "reference misassembly" in page
    table = page[page.index('<table class="nested-table">'):]
    table = table[: table.index("</table>")]
    assert table.index("chr1:1,000,000&ndash;1,050,000") < table.index("chr2:7,000,000")
    assert "chr5:1&ndash;2 (+)" in table
    assert '<tr class="recurrent">' in table
    assert "Nested-block table does not match summary" not in page
    assert '<section class="nested">' in page
    assert "No nested blocks for this sample." in page  # S2


def test_nested_mismatch_note_and_not_supplied(tmp_path):
    a = tmp_path / "a"
    a.mkdir()
    _, _, page = run_report(a, cohort(2))
    assert "No nested-block tables were supplied" in page
    b = tmp_path / "b"
    b.mkdir()
    nb = write_nested(b, "S0", [nested_row("S0", ("chr1", 1, 2))])
    _, _, page = run_report(b, cohort(2), nested=[nb])  # summaries say 4 nested blocks
    assert "Nested-block table does not match summary" in page


def test_nested_tables_capped(tmp_path):
    rows = cohort(5)
    many = [nested_row("S0", ("chr1", i * 10_000_000, i * 10_000_000 + 100), query_contig=f"q{i}")
            for i in range(250)]
    nb = write_nested(tmp_path, "S0", many)
    _, _, page = run_report(tmp_path, rows, nested=[nb])
    assert "Showing 200 of 250 nested blocks; the combined nested-blocks TSV has every row." in page
    assert "Showing 50 of 250 nested blocks" in page


def test_nested_html_escapes(tmp_path):
    evil = "<b>&x"
    rows = cohort(2)
    rows[0]["sample"] = evil
    nb = write_nested(tmp_path, "e", [nested_row(evil, ('chr"<1>', 5, 9), query_contig="<q&>",
                                                 container=("<c>", 1, 2))])
    _, _, page = run_report(tmp_path, rows, nested=[nb])
    assert "<q&>" not in page and "&lt;q&amp;&gt;" in page
    assert "chr&quot;&lt;1&gt;:5" in page
    assert "&lt;c&gt;:1" in page


def test_overlapping_query_bp_flagged_high_and_nested_never_flagged():
    rows = cohort(8)
    rows[1]["overlapping_query_bp"] = "90000000"
    rows[4]["overlapping_query_bp"] = "0"  # low is fine
    for i, r in enumerate(rows):
        r["nested_private"] = "1"
        r["nested_recurrent"] = "1"
    rows[3].update(nested_blocks="99999", nested_block_reference_bp="900000000",
                   nested_private="99999")
    result = compute_flags(rows)
    assert result.flags[1]["overlapping_query_bp"][0].startswith("robust_z=")
    assert "overlapping_query_bp" not in result.flags[4]
    assert result.flags[3] == {}
    # median 800 kb, MAD 0 -> scale = max(1e5, 0.05 * 8e5 = 4e4) = 1e5.
    assert result.stats["overlapping_query_bp"].scale == pytest.approx(100_000)


def test_query_table_new_columns(tmp_path):
    rows = cohort(5)
    q = query_row("S0", "chr1q")
    q.update(overlapping_query_bp="12345", nested_blocks="3", contig_pieces="7", major="true")
    query = write_tsv(tmp_path / "S0.by_query_contig.tsv", BY_QUERY_COLUMNS, [q])
    _, _, page = run_report(tmp_path, rows, query=[query])
    table = page[page.index('<table class="query-table">'):]
    table = table[: table.index("</table>")]
    labels = re.findall(r"<th>([^<]*)</th>", table)
    for label in ("Overlapping query bp", "Nested blocks", "Contig pieces", "Major"):
        assert label in labels
    assert "12,345" in table and "<td>true</td>" in table


def test_telomere_cell_na_safe(tmp_path):
    rows = [summary_row("A", assembly=True), summary_row("B", assembly=True)]
    rows[1]["assembly_major_telomeric_ends"] = "NA"
    _, _, page = run_report(tmp_path, rows)
    table = page[page.index('<table class="assembly">'):]
    table = table[: table.index("</table>")]
    assert "<td>14 / 20</td>" in table
    assert '<td><span class="na">NA</span> / 20</td>' in table


def test_definitions_cover_new_columns(tmp_path):
    _, _, page = run_report(tmp_path, [summary_row("A", assembly=True)])
    defs = page[page.index('<dl class="defs">'):]
    for name in ("overlapping_query_bp", "nested_blocks", "nested_block_reference_bp",
                 "nested_recurrent", "nested_private", "assembly_sequences",
                 "assembly_scaffold_n50_bp", "assembly_contig_pieces", "assembly_contig_n50_bp",
                 "assembly_major_sequences", "assembly_major_telomeric_ends",
                 "assembly_minor_sequences_with_telomere", "contig_pieces", "major"):
        assert f"<dt>{name}</dt>" in defs, name
    assert "1,000,000" in defs
    assert "CCCTAAA" in defs and "terminal 1 kb" in defs
    assert "no lowercase" in defs
    # nested_blocks definition has no "Flagged when" suffix.
    nb = defs[defs.index("<dt>nested_blocks</dt>"):]
    nb = nb[: nb.index("</dd>")]
    assert "Flagged when" not in nb
