import re
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.maf_stats_report import (  # noqa: E402
    BY_QUERY_COLUMNS,
    BY_REFERENCE_COLUMNS,
    SUMMARY_COLUMNS,
    build_report,
    compute_flags,
    format_flags,
    load_inputs,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "maf_stats_report.py"


def summary_row(sample: str, **overrides) -> dict[str, str]:
    row = {
        "sample": sample,
        "reference_length_bp": "1000000000",
        "covered_reference_bp": "900000000",
        "reference_coverage_pct": "90.00",
        "reference_bp_aligned_to_query_base": "850000000",
        "query_length_bp": "1100000000",
        "query_length_source": "query_fai",
        "covered_query_bp": "880000000",
        "query_coverage_pct": "80.00",
        "identity_matches": "800000000",
        "identity_compared_columns": "840000000",
        "identity_pct": "95.24",
        "gap_columns": "50000000",
        "alignment_columns": "950000000",
        "gap_fraction_pct": "5.26",
        "blocks": "20000",
        "block_n50_bp": "150000",
        "overlapping_reference_bp": "1000000",
        "strand_flips": "10",
        "reference_contig_jumps": "5",
        "out_of_order_adjacencies": "3",
        "breakpoint_adjacencies": "15",
        "breakpoints_per_gb_covered": "16.67",
        "min_block_bp": "0",
        "overlap_tolerance_bp": "0",
    }
    row.update({k: str(v) for k, v in overrides.items()})
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
        "covered_reference_bp": "900000",
        "reference_coverage_pct": "90.00",
        "reference_bp_aligned_to_query_base": "850000",
        "identity_matches": "800000",
        "identity_compared_columns": "840000",
        "identity_pct": "95.24",
        "gap_columns": "50000",
        "alignment_columns": "950000",
        "gap_fraction_pct": "5.26",
        "blocks": "20",
        "block_n50_bp": "150000",
        "overlapping_reference_bp": "0",
        "breakpoint_adjacencies": str(breakpoints),
        "dotplot": dotplot,
    }


def query_row(sample: str, contig: str, breakpoints: int = 0) -> dict[str, str]:
    return {
        "sample": sample,
        "query_contig": contig,
        "query_length_bp": "1000000",
        "query_length_source": "query_fai",
        "covered_query_bp": "800000",
        "query_coverage_pct": "80.00",
        "blocks": "20",
        "blocks_considered": "20",
        "strand_flips": str(breakpoints),
        "reference_contig_jumps": "0",
        "out_of_order_adjacencies": "0",
        "breakpoint_adjacencies": str(breakpoints),
        "breakpoints_per_gb_covered": "NA" if breakpoints == 0 else "1250.00",
    }


def cohort(n: int, **per_sample_overrides) -> list[dict[str, str]]:
    """n samples with slightly varied metrics; overrides keyed by sample index."""
    rows = []
    for i in range(n):
        rows.append(
            summary_row(
                f"S{i}",
                reference_coverage_pct=f"{90 + (i % 3) * 0.5:.2f}",
                identity_pct=f"{95 + (i % 3) * 0.2:.2f}",
                gap_fraction_pct=f"{5 + (i % 3) * 0.3:.2f}",
                query_coverage_pct=f"{80 + (i % 3) * 0.5:.2f}",
                blocks=str(20000 + (i % 3) * 100),
            )
        )
    for idx, overrides in per_sample_overrides.items():
        rows[int(idx.lstrip("s"))].update({k: str(v) for k, v in overrides.items()})
    return rows


def read_out_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").rstrip("\n").split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def run_report(tmp_path, rows, ref=None, query=None, **kwargs):
    paths = [write_summary(tmp_path, r) for r in rows]
    out_tsv = tmp_path / "out" / "maf_stats.tsv"
    out_html = tmp_path / "out" / "maf_stats.html"
    result = build_report(paths, ref or [], query or [], out_tsv, out_html, **kwargs)
    return result, read_out_tsv(out_tsv), out_html.read_text(encoding="utf-8")


# ── aggregation ─────────────────────────────────────────────────────────────


def test_one_row_per_sample_in_input_order(tmp_path):
    rows = [summary_row(s) for s in ["zeta", "alpha", "mid"]]
    _, out, _ = run_report(tmp_path, rows)
    assert [r["sample"] for r in out] == ["zeta", "alpha", "mid"]
    header = (tmp_path / "out" / "maf_stats.tsv").read_text().split("\n")[0].split("\t")
    assert header == SUMMARY_COLUMNS + ["flags"]
    # Values pass through unchanged.
    assert out[1]["covered_reference_bp"] == "900000000"
    assert all(r["flags"] == "" for r in out)


# ── relative flags ───────────────────────────────────────────────────────────


def test_low_is_bad_outlier_flagged_and_high_outlier_not(tmp_path):
    rows = cohort(8, s2={"reference_coverage_pct": "40.00"}, s5={"reference_coverage_pct": "99.90"})
    result, out, _ = run_report(tmp_path, rows)
    assert "reference_coverage_pct" in result.flags[2]
    assert out[2]["flags"].startswith("reference_coverage_pct:robust_z=-")
    # A high-coverage outlier is good, not flagged.
    assert "reference_coverage_pct" not in result.flags[5]
    assert "reference_coverage_pct" not in out[5]["flags"]


def test_high_is_bad_outlier_flagged_and_low_not():
    rows = cohort(
        8,
        s1={"gap_fraction_pct": "30.00", "blocks": "200000"},
        s4={"gap_fraction_pct": "0.10", "blocks": "100"},
    )
    result = compute_flags(rows)
    assert "gap_fraction_pct" in result.flags[1]
    assert "blocks" in result.flags[1]
    assert result.flags[1]["gap_fraction_pct"][0].startswith("robust_z=")
    assert float(result.flags[1]["gap_fraction_pct"][0].split("=")[1].split("(")[0]) > 3.5
    assert "gap_fraction_pct" not in result.flags[4]
    assert "blocks" not in result.flags[4]


def test_min_samples_guard_skips_relative_flags(tmp_path):
    rows = cohort(4, s0={"reference_coverage_pct": "5.00"})
    result, out, page = run_report(tmp_path, rows)
    assert all(not f for f in result.flags)
    assert all(r["flags"] == "" for r in out)
    assert "Relative (cohort) flagging skipped" in page
    # Same data with a lower guard is flagged.
    assert "reference_coverage_pct" in compute_flags(rows, min_samples=3).flags[0]


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
    # blocks MAD is 0; scale = max(10, 0.05*20000=1000) -> 1000.
    rows = [summary_row(f"S{i}", blocks="20000") for i in range(6)]
    rows[0]["blocks"] = "23000"  # z = 3.0, below threshold
    rows[1]["blocks"] = "24000"  # z = 4.0, flagged
    result = compute_flags(rows)
    assert result.stats["blocks"].scale == pytest.approx(1000.0)
    assert "blocks" not in result.flags[0]
    assert result.flags[1]["blocks"] == ["robust_z=4.00(mad0_floor)"]


# ── absolute thresholds ─────────────────────────────────────────────────────


def test_absolute_thresholds(tmp_path):
    rows = [
        summary_row("A", identity_pct="94.00", gap_fraction_pct="7.5"),
        summary_row("B", identity_pct="95.00", gap_fraction_pct="5.0"),
        summary_row("C", identity_pct="99.00", gap_fraction_pct="4.0", breakpoints_per_gb_covered="250"),
    ]
    thresholds = {
        "identity_pct": 95.0,
        "gap_fraction_pct": 5.0,
        "breakpoints_per_gb_covered": 100.0,
        "reference_coverage_pct": 50.0,
    }
    result, out, page = run_report(tmp_path, rows, thresholds=thresholds)
    assert out[0]["flags"] == (
        "identity_pct:below_threshold_95;gap_fraction_pct:above_threshold_5"
    )
    assert out[1]["flags"] == ""  # equal to threshold is not flagged
    assert out[2]["flags"] == "breakpoints_per_gb_covered:above_threshold_100"
    assert 'title="identity_pct:below_threshold_95"' in page
    # n=3 < 5: relative flagging skipped but absolute thresholds still applied.
    assert "Relative (cohort) flagging skipped" in page


# ── NA handling ──────────────────────────────────────────────────────────────


def test_na_values_are_ignored(tmp_path):
    rows = cohort(6)
    rows[0]["identity_pct"] = "NA"
    rows[0]["breakpoints_per_gb_covered"] = "NA"
    rows[1]["breakpoints_per_gb_covered"] = "NA"
    result, out, page = run_report(
        tmp_path, rows, thresholds={"identity_pct": 99.0, "breakpoints_per_gb_covered": 1.0}
    )
    assert "identity_pct" not in result.flags[0]
    assert result.stats["identity_pct"].n == 5
    # breakpoints has only 4 non-NA values -> relative flagging skipped for it.
    assert result.skipped["breakpoints_per_gb_covered"] == 4
    assert out[0]["identity_pct"] == "NA"
    assert '<span class="na">NA</span>' in page
    assert "Relative flagging skipped</strong> for metrics" in page


def test_all_na_metric_does_not_crash(tmp_path):
    rows = cohort(5)
    for r in rows:
        r["query_coverage_pct"] = "NA"
    result, _, page = run_report(tmp_path, rows)
    assert result.stats["query_coverage_pct"] is None
    assert "no values" in page


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


def test_by_file_with_unknown_sample_error(tmp_path):
    s = write_summary(tmp_path, summary_row("A"))
    ref = write_tsv(tmp_path / "B.by_reference_contig.tsv", BY_REFERENCE_COLUMNS, [ref_row("B", "chr1")])
    with pytest.raises(ValueError, match="'B'.*not among the summary inputs"):
        load_inputs([s], [ref], [])


def test_by_files_matched_by_sample_column_not_filename(tmp_path):
    s = write_summary(tmp_path, summary_row("A"))
    ref = write_tsv(tmp_path / "unrelated_name.tsv", BY_REFERENCE_COLUMNS, [ref_row("A", "chr1")])
    _, by_ref, _ = load_inputs([s], [ref], [])
    assert [r["reference_contig"] for r in by_ref["A"]] == ["chr1"]


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
    _, out, page = run_report(tmp_path, rows, ref=[ref])
    assert out[0]["sample"] == evil
    assert "<b>&x" not in page
    assert "&lt;b&gt;&amp;x" in page
    assert 'chr"<1>' not in page
    assert 'src="dotplots/x/chr&quot;&lt;1&gt;.png"' in page


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
        [query_row("S0", "scaf1"), query_row("S0", "scaf2", breakpoints=3)],
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

    # Separate reference and query sections with headings.
    assert '<section class="ref-contigs">' in page
    assert '<section class="query-contigs">' in page
    assert "<h3>Reference contigs</h3>" in page and "<h3>Query contigs</h3>" in page
    query_table = page[page.index('<table class="query-table">'):]
    query_table = query_table[: query_table.index("</table>")]
    assert "scaf1" in query_table and "scaf2" in query_table
    assert "chr1" not in query_table


def test_html_strip_plots_and_flag_highlighting(tmp_path):
    rows = cohort(8, s2={"reference_coverage_pct": "40.00"})
    _, _, page = run_report(tmp_path, rows)
    assert page.count('<svg class="strip"') == 7
    assert '<circle class="dot flagged"' in page
    assert "<title>S2: 40.00</title>" in page
    assert re.search(r'<td class="flag" title="reference_coverage_pct:robust_z=-[\d.]+">40\.00</td>', page)
    assert "maf_srcsize" in page and "overstates" in page
    assert "raw block-column metrics" in page


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
    out_tsv = tmp_path / "results" / "maf_stats.tsv"
    out_html = tmp_path / "results" / "maf_stats.html"
    proc = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--summaries", *map(str, summaries),
            "--by-reference", str(ref),
            "--by-query", str(query),
            "--out-tsv", str(out_tsv),
            "--out-html", str(out_html),
            "--min-samples", "5",
            "--z-threshold", "3.5",
            "--flag-min-reference-coverage", "50",
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
    out = read_out_tsv(out_tsv)
    assert [r["sample"] for r in out] == [f"S{i}" for i in range(6)]
    assert out[3]["flags"].startswith("identity_pct:robust_z=-")
    page = out_html.read_text(encoding="utf-8")
    assert 'src="dotplots/S3/chr1.png"' in page
    assert "80.50 &dagger;" in page

    # Error path: duplicate sample gives a clear non-zero exit.
    proc = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--summaries", str(summaries[0]), str(summaries[0]),
            "--out-tsv", str(out_tsv), "--out-html", str(out_html),
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
