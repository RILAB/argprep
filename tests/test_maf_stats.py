import gzip
import struct
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts import maf_stats
from scripts.maf_stats import (
    Block,
    MafValidationError,
    classify_breakpoints,
    column_counts,
    n50_l50,
    run,
    scan_maf,
    select_dotplot_contigs,
)

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "maf_stats.py"
DOTPLOT_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "maf_dotplot.py"


def _row(src, start, strand, src_size, text):
    size = len(text) - text.count("-")
    return f"s {src} {start} {size} {strand} {src_size} {text}\n"


def _write_maf(path: Path, blocks: list[list[tuple]]) -> Path:
    """Each block is a list of (src, start, strand, src_size, text) rows."""
    parts = ["##maf version=1\n\n"]
    for rows in blocks:
        parts.append("a score=0\n")
        parts.extend(_row(*row) for row in rows)
        parts.append("\n")
    text = "".join(parts)
    if path.name.endswith(".gz"):
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return path


def _write_fai(path: Path, lengths: dict[str, int]) -> Path:
    path.write_text("".join(f"{name}\t{length}\t0\t60\t61\n" for name, length in lengths.items()))
    return path


def _read_tsv(path: Path) -> list[dict[str, str]]:
    header, *rows = path.read_text().splitlines()
    columns = header.split("\t")
    return [dict(zip(columns, row.split("\t"))) for row in rows]


def _block(ref, rs, re_, qs, qe, strand="+", query="q1"):
    return Block(ref, rs, re_, query, qs, qe, strand)


# ---------------------------------------------------------------- column metrics


def test_column_counts_identity_indels_lowercase_and_n():
    counts = column_counts("ACGTAC-GT", "ACGAacgN-")
    assert counts.columns == 9
    assert counts.aligned == 7
    assert counts.insertions == 1  # query base opposite a reference gap
    assert counts.deletions == 1  # reference base opposite a query gap
    assert counts.compared == 6
    assert counts.matches == 5  # T/A mismatch; lowercase query bases still match
    assert counts.query_n == 1


def test_column_counts_excludes_ambiguity_codes_from_identity():
    counts = column_counts("ACGNTR", "ACGTTA")
    assert counts.compared == 4
    assert counts.matches == 4
    assert counts.aligned == 6


def test_column_counts_is_chunk_invariant(monkeypatch):
    ref, query = "ACGT--ACGTAC" * 50, "AC-TGGACnTA-" * 50
    whole = column_counts(ref, query)
    monkeypatch.setattr(maf_stats, "COLUMN_CHUNK", 7)
    assert column_counts(ref, query) == whole


def test_n50_l50():
    assert n50_l50([3, 10, 5]) == (10, 1)
    assert n50_l50([4, 4, 4, 4]) == (4, 2)
    assert n50_l50([]) == (None, None)


# ---------------------------------------------------------------- scan + summary


def test_summary_coverage_overlap_and_minus_strand_query(tmp_path):
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 20, "chr2": 10})
    maf = _write_maf(tmp_path / "s1.maf", [
        [("chr1", 0, "+", 20, "ACGTACGTAC"), ("q1", 0, "+", 30, "ACGTACGTAC")],
        [("chr1", 5, "+", 20, "CGTACGTACG"), ("q1", 2, "-", 30, "CGTACGTACG")],
    ])
    stats = run(maf, fai, "s1", tmp_path / "out", dotplots="false")
    s = stats.summary
    assert s["reference_length_bp"] == 30
    assert s["block_span_reference_bp"] == 15
    assert s["overlapping_reference_bp"] == 5
    assert s["block_span_reference_pct"] == pytest.approx(50.0)
    assert s["aligned_reference_bp"] == 20  # raw column count: overlapping blocks count twice
    # minus-strand query row start=2 size=10 srcSize=30 -> forward [18, 28)
    assert s["block_span_query_bp"] == 20
    assert s["query_length_source"] == "maf_srcsize"
    assert s["block_span_query_pct"] == pytest.approx(20 / 30 * 100)
    assert s["aligned_query_contigs"] == 1
    assert s["assembly_length_bp"] is None
    assert s["blocks"] == 2
    assert s["identity_pct"] == pytest.approx(100.0)

    by_ref = {row["reference_contig"]: row for row in _read_tsv(tmp_path / "out" / "s1.by_reference_contig.tsv")}
    assert set(by_ref) == {"chr1", "chr2"}  # uncovered contigs are reported too
    assert by_ref["chr2"]["block_span_reference_bp"] == "0"
    assert by_ref["chr2"]["identity_pct"] == "NA"
    summary_rows = _read_tsv(tmp_path / "out" / "s1.maf_stats.tsv")
    assert len(summary_rows) == 1 and summary_rows[0]["sample"] == "s1"


def test_query_fai_denominator_includes_unaligned_contigs(tmp_path):
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 10})
    qfai = _write_fai(tmp_path / "q.fai", {"q1": 10, "unaligned": 30})
    maf = _write_maf(tmp_path / "s1.maf", [
        [("chr1", 0, "+", 10, "ACGTACGTAC"), ("q1", 0, "+", 10, "ACGTACGTAC")],
    ])
    without = run(maf, fai, "s1", tmp_path / "a", dotplots="false").summary
    with_fai = run(maf, fai, "s1", tmp_path / "b", query_fai=qfai, dotplots="false").summary
    assert without["aligned_query_pct"] == pytest.approx(100.0)
    assert with_fai["query_length_bp"] == 40
    assert with_fai["query_length_source"] == "query_fai"
    assert with_fai["aligned_query_pct"] == pytest.approx(25.0)


def test_reference_names_are_normalized_only_when_unique(tmp_path):
    fai = _write_fai(tmp_path / "ref.fai", {"1": 10})
    maf = _write_maf(tmp_path / "s1.maf", [
        [("chr1", 0, "+", 10, "ACGTACGTAC"), ("q1", 0, "+", 10, "ACGTACGTAC")],
    ])
    scan = scan_maf(maf, {"1": 10})
    assert scan.reference["1"].block_sizes == [10]
    with pytest.raises(MafValidationError, match="ambiguous"):
        scan_maf(maf, {"1": 10, "chr01": 10})
    with pytest.raises(MafValidationError, match="not in the reference"):
        scan_maf(maf, {"chr2": 10})


def test_query_names_with_assembly_file_prefix_resolve_to_fasta_names(tmp_path):
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 10})
    maf = _write_maf(tmp_path / "s1.maf", [
        [("chr1", 0, "+", 10, "ACGTACGTAC"), ("Zm-X-1.0.fa.chr1", 0, "+", 10, "ACGTACGTAC")],
        [("chr1", 0, "+", 10, "ACGTACGTAC"), ("Zm-X-1.0.fa.scaf_7", 0, "+", 12, "ACGTACGTAC")],
    ])
    scan = scan_maf(maf, {"chr1": 10}, {"chr1": 10, "scaf_7": 12, "chr2": 5}, "query FASTA")
    assert set(scan.query) == {"chr1", "scaf_7"}
    with pytest.raises(MafValidationError, match="query contig 'Zm-X-1.0.fa.scaf_7' is not in the query FASTA"):
        scan_maf(maf, {"chr1": 10}, {"chr1": 10}, "query FASTA")


def test_gzipped_maf_is_supported(tmp_path):
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 10})
    maf = _write_maf(tmp_path / "s1.maf.gz", [
        [("chr1", 0, "+", 10, "ACGTACGTAC"), ("q1", 0, "+", 10, "ACGTACGTAC")],
    ])
    assert run(maf, fai, "s1", tmp_path / "out", dotplots="false").summary["blocks"] == 1


@pytest.mark.parametrize(
    "blocks, refs, message",
    [
        ([[("chr1", 0, "+", 10, "ACGT"), ("q1", 0, "+", 10, "ACGT"), ("q2", 0, "+", 10, "ACGT")]],
         {"chr1": 10}, "exactly 2 sequence rows"),
        ([[("chr1", 0, "+", 10, "ACGT")]], {"chr1": 10}, "exactly 2 sequence rows"),
        ([[("chr1", 8, "+", 10, "ACGT"), ("q1", 0, "+", 10, "ACGT")]], {"chr1": 10}, "exceeds srcSize"),
        ([[("chr1", 0, "*", 10, "ACGT"), ("q1", 0, "+", 10, "ACGT")]], {"chr1": 10}, "invalid strand"),
        ([[("chr1", 0, "-", 10, "ACGT"), ("q1", 0, "+", 10, "ACGT")]], {"chr1": 10}, "minus strand"),
        ([[("chr1", 0, "+", 12, "ACGT"), ("q1", 0, "+", 10, "ACGT")]], {"chr1": 10}, "does not match reference .fai"),
        ([[("chr1", 0, "+", 10, "ACGT"), ("q1", 0, "+", 10, "ACGT")],
          [("chr1", 4, "+", 10, "ACGT"), ("q1", 4, "+", 11, "ACGT")]], {"chr1": 10}, "inconsistent srcSize"),
    ],
)
def test_validation_errors_carry_block_context(tmp_path, blocks, refs, message):
    maf = _write_maf(tmp_path / "bad.maf", blocks)
    with pytest.raises(MafValidationError, match=message) as info:
        scan_maf(maf, refs)
    assert "block" in str(info.value)


def test_size_field_must_match_ungapped_length(tmp_path):
    maf = tmp_path / "bad.maf"
    maf.write_text("a\ns chr1 0 5 + 10 ACGT\ns q1 0 4 + 10 ACGT\n\n")
    with pytest.raises(MafValidationError, match="ungapped bases"):
        scan_maf(maf, {"chr1": 10})


def test_query_fai_mismatch_and_unknown_query_contig(tmp_path):
    maf = _write_maf(tmp_path / "s.maf", [[("chr1", 0, "+", 10, "ACGT"), ("q1", 0, "+", 10, "ACGT")]])
    with pytest.raises(MafValidationError, match="query .fai length"):
        scan_maf(maf, {"chr1": 10}, {"q1": 11})
    with pytest.raises(MafValidationError, match="is not in the query .fai"):
        scan_maf(maf, {"chr1": 10}, {"other": 10})


def test_empty_maf_is_an_error(tmp_path):
    maf = tmp_path / "empty.maf"
    maf.write_text("##maf version=1\n")
    with pytest.raises(MafValidationError, match="no alignment blocks"):
        scan_maf(maf, {"chr1": 10})


# ---------------------------------------------------------------- breakpoints

ORDER = {"chr1": 0, "chr2": 1}


def test_collinear_forward_and_reverse_runs_have_no_breakpoints():
    forward = [_block("chr1", 0, 10, 0, 10), _block("chr1", 20, 30, 10, 20)]
    reverse = [_block("chr1", 100, 110, 0, 10, "-"), _block("chr1", 90, 100, 10, 20, "-")]
    assert classify_breakpoints(forward, ORDER, 0, 0).breakpoint_adjacencies == 0
    assert classify_breakpoints(reverse, ORDER, 0, 0).breakpoint_adjacencies == 0


@pytest.mark.parametrize("strand, second", [("+", (50, 60)), ("-", (120, 130))])
def test_out_of_order_on_each_strand(strand, second):
    blocks = [_block("chr1", 100, 110, 0, 10, strand), _block("chr1", *second, 10, 20, strand)]
    counts = classify_breakpoints(blocks, ORDER, 0, 0)
    assert counts.out_of_order_adjacencies == 1
    assert counts.strand_flips == counts.reference_contig_jumps == 0
    assert counts.breakpoint_adjacencies == 1


def test_strand_flip_and_jump_are_recorded_independently():
    flip = [_block("chr1", 0, 10, 0, 10, "+"), _block("chr1", 20, 30, 10, 20, "-")]
    jump_and_flip = [_block("chr1", 0, 10, 0, 10, "+"), _block("chr2", 0, 10, 10, 20, "-")]
    counts = classify_breakpoints(flip, ORDER, 0, 0)
    assert (counts.strand_flips, counts.reference_contig_jumps, counts.breakpoint_adjacencies) == (1, 0, 1)
    counts = classify_breakpoints(jump_and_flip, ORDER, 0, 0)
    assert (counts.strand_flips, counts.reference_contig_jumps) == (1, 1)
    assert counts.out_of_order_adjacencies == 0
    assert counts.breakpoint_adjacencies == 1  # one adjacency, two properties


def test_overlap_tolerance():
    blocks = [_block("chr1", 0, 100, 0, 100), _block("chr1", 95, 200, 100, 205)]
    assert classify_breakpoints(blocks, ORDER, 0, 5).out_of_order_adjacencies == 0
    assert classify_breakpoints(blocks, ORDER, 0, 4).out_of_order_adjacencies == 1


def test_min_block_bp_skips_small_blocks_for_breakpoints_only():
    blocks = [
        _block("chr1", 0, 100, 0, 100),
        _block("chr2", 0, 5, 100, 105),  # tiny spurious hit
        _block("chr1", 100, 200, 105, 205),
    ]
    assert classify_breakpoints(blocks, ORDER, 0, 0).breakpoint_adjacencies == 2
    counts = classify_breakpoints(blocks, ORDER, 10, 0)
    assert counts.breakpoint_adjacencies == 0
    assert counts.blocks_considered == 2


def test_adjacency_order_is_deterministic_for_duplicated_query_regions():
    a = _block("chr1", 0, 10, 0, 10)
    b = _block("chr2", 0, 10, 0, 10)
    c = _block("chr1", 10, 20, 10, 20)
    first = classify_breakpoints([a, b, c], ORDER, 0, 0)
    second = classify_breakpoints([c, b, a], ORDER, 0, 0)
    assert first == second


def test_breakpoints_are_grouped_by_query_contig_and_credited_to_reference(tmp_path):
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 40, "chr2": 40})
    maf = _write_maf(tmp_path / "s1.maf", [
        [("chr1", 0, "+", 40, "ACGTACGTAC"), ("q1", 0, "+", 40, "ACGTACGTAC")],
        [("chr2", 0, "+", 40, "ACGTACGTAC"), ("q1", 10, "+", 40, "ACGTACGTAC")],
        # q2 is collinear on chr1: no breakpoint even though it interleaves with q1 in the file
        [("chr1", 20, "+", 40, "ACGTACGTAC"), ("q2", 0, "+", 40, "ACGTACGTAC")],
        [("chr1", 30, "+", 40, "ACGTACGTAC"), ("q2", 10, "+", 40, "ACGTACGTAC")],
    ])
    stats = run(maf, fai, "s1", tmp_path / "out", dotplots="false")
    assert stats.summary["reference_contig_jumps"] == 1
    assert stats.summary["breakpoint_adjacencies"] == 1
    by_query = {row["query_contig"]: row for row in stats.by_query}
    assert by_query["q1"]["breakpoint_adjacencies"] == 1
    assert by_query["q2"]["breakpoint_adjacencies"] == 0
    by_ref = {row["reference_contig"]: row for row in stats.by_reference}
    assert by_ref["chr1"]["breakpoint_adjacencies"] == 1
    assert by_ref["chr2"]["breakpoint_adjacencies"] == 1
    assert stats.summary["breakpoints_per_gb_aligned"] == pytest.approx(1 / 40 * 1e9)


# ---------------------------------------------------------------- dotplots


def _png_size(path: Path) -> tuple[int, int]:
    data = path.read_bytes()[:24]
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    return struct.unpack(">II", data[16:24])


def _structural_maf(tmp_path: Path) -> tuple[Path, Path]:
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 40, "chr2": 40, "chr3": 40})
    maf = _write_maf(tmp_path / "s1.maf", [
        [("chr1", 0, "+", 40, "ACGTACGTAC"), ("q1", 0, "+", 40, "ACGTACGTAC")],
        [("chr1", 10, "+", 40, "ACGTACGTAC"), ("q1", 10, "-", 40, "ACGTACGTAC")],
        [("chr2", 0, "+", 40, "ACGTACGTAC"), ("q2", 0, "+", 40, "ACGTACGTAC")],
        [("chr2", 10, "+", 40, "ACGTACGTAC"), ("q3", 0, "+", 20, "ACGTACGTAC")],
    ])
    return maf, fai


@pytest.mark.parametrize("mode, expected", [("all", {"chr1", "chr2"}), ("flagged", {"chr1"}), ("false", set())])
def test_dotplot_modes(tmp_path, mode, expected):
    maf, fai = _structural_maf(tmp_path)
    out = tmp_path / "out"
    stats = run(maf, fai, "s1", out, dotplots=mode)
    plotted = {row["reference_contig"] for row in stats.by_reference if row["dotplot"]}
    assert plotted == expected
    assert (out / "dotplots" / "s1").is_dir()
    for row in stats.by_reference:
        if row["dotplot"]:
            png = out / row["dotplot"]
            assert row["dotplot"] == f"dotplots/s1/{row['reference_contig']}.png"
            assert _png_size(png) == (800, 800)
    assert len(list((out / "dotplots" / "s1").glob("*.png"))) == len(expected)


def test_flagged_dotplots_are_capped_densest_first():
    rows = [
        {"reference_contig": "sparse", "blocks": 5, "breakpoint_adjacencies": 1, "block_span_reference_bp": 1000},
        {"reference_contig": "dense", "blocks": 5, "breakpoint_adjacencies": 5, "block_span_reference_bp": 1000},
        {"reference_contig": "clean", "blocks": 5, "breakpoint_adjacencies": 0, "block_span_reference_bp": 1000},
        {"reference_contig": "empty", "blocks": 0, "breakpoint_adjacencies": 0, "block_span_reference_bp": 0},
    ]
    assert select_dotplot_contigs(rows, "flagged", 1) == ["dense"]
    assert select_dotplot_contigs(rows, "flagged", 0) == ["dense", "sparse"]
    assert select_dotplot_contigs(rows, "all", 1) == ["sparse", "dense", "clean"]


def test_dotplot_segments_stack_query_contigs_in_bands():
    from scripts.maf_dotplot import segment

    forward = _block("chr1", 0, 1_000_000, 0, 1_000_000, "+", "q2")
    reverse = _block("chr1", 0, 1_000_000, 0, 1_000_000, "-", "q2")
    assert segment(forward, 5_000_000) == ((0.0, 5.0), (1.0, 6.0))
    assert segment(reverse, 5_000_000) == ((0.0, 6.0), (1.0, 5.0))


def test_dotplot_cli_writes_pngs_and_fails_loudly(tmp_path):
    maf, fai = _structural_maf(tmp_path)
    bad = tmp_path / "bad.maf"
    bad.write_text("a\ns chr1 0 4 + 40 ACGT\n\n")
    ok = subprocess.run(
        [sys.executable, str(DOTPLOT_SCRIPT), "--maf", str(maf), "--reference-fai", str(fai),
         "--out-dir", str(tmp_path / "plots")],
        capture_output=True, text=True,
    )
    assert ok.returncode == 0, ok.stderr
    assert {p.name for p in (tmp_path / "plots" / "dotplots" / "s1").glob("*.png")} == {"chr1.png", "chr2.png"}
    failed = subprocess.run(
        [sys.executable, str(DOTPLOT_SCRIPT), "--maf", str(maf), str(bad), "--reference-fai", str(fai),
         "--out-dir", str(tmp_path / "plots2")],
        capture_output=True, text=True,
    )
    assert failed.returncode == 1
    assert "FAIL" in failed.stderr


# ---------------------------------------------------------------- CLI


def test_cli_reports_validation_errors_without_traceback(tmp_path):
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 10})
    maf = tmp_path / "empty.maf"
    maf.write_text("")
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--maf", str(maf), "--reference-fai", str(fai),
         "--sample", "s1", "--out-dir", str(tmp_path / "out")],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "no alignment blocks" in result.stderr
    assert "Traceback" not in result.stderr


def test_cli_writes_three_tsvs(tmp_path):
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 10})
    maf = _write_maf(tmp_path / "s1.maf", [[("chr1", 0, "+", 10, "ACGTACGTAC"), ("q1", 0, "+", 10, "ACGTACGTAC")]])
    out = tmp_path / "out"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--maf", str(maf), "--reference-fai", str(fai),
         "--sample", "s1", "--out-dir", str(out), "--dotplots", "false"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    for name in ("s1.maf_stats.tsv", "s1.by_reference_contig.tsv", "s1.by_query_contig.tsv"):
        rows = _read_tsv(out / name)
        assert rows and all(len(row) == len(rows[0]) for row in rows)
    assert list(_read_tsv(out / "s1.maf_stats.tsv")[0]) == maf_stats.SUMMARY_COLUMNS


# ---------------------------------------------------------------- assembly FASTA

TELOMERE = "TTTAGGG" * 150  # 1,050 bp array: fills the terminal 1 kb window


def _write_fasta(path: Path, records: dict[str, str]) -> Path:
    text = "".join(f">{name} description\n{seq}\n" for name, seq in records.items())
    if path.name.endswith(".gz"):
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return path


def test_assembly_metrics_and_unaligned_contigs(tmp_path):
    # 5 Ns are an ambiguous stretch, not a gap; 12 Ns are a scaffold gap
    q1 = TELOMERE + "ACGTACGTAC" + "A" * 20 + "NNNNN" + "GGCC" * 5 + "C" * 20_000 + "N" * 12 + "acgtacgtac"
    unaligned = "N" * 10 + "acgt" * 5
    fasta = _write_fasta(tmp_path / "s1.fa.gz", {"q1": q1, "unaligned": unaligned})
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 20})
    maf = _write_maf(tmp_path / "s1.maf", [
        [("chr1", 0, "+", 20, "ACGTACGTAC"), ("q1", len(TELOMERE), "+", len(q1), "ACGTACGTAC")],
    ])
    stats = run(maf, fai, "s1", tmp_path / "out", fasta=fasta, dotplots="false")
    s = stats.summary
    assert s["query_length_source"] == "query_fasta"
    assert s["assembly_length_bp"] == len(q1) + len(unaligned)
    assert s["query_length_bp"] == s["assembly_length_bp"]
    assert s["assembly_sequences"] == 2
    assert s["assembly_scaffold_n50_bp"] == len(q1)
    assert s["assembly_scaffold_l50"] == 1
    assert s["assembly_n_bp"] == 27
    assert s["assembly_n_gaps"] == 2  # q1's 12-N run and the unaligned contig's 10-N run
    # q1 splits at its 12-N gap; the unaligned contig keeps one piece after its leading gap
    assert s["assembly_contig_pieces"] == 3
    assert s["assembly_contig_n50_bp"] == len(q1) - 12 - 10
    assert s["assembly_major_sequences"] == 1  # the unaligned contig is under 10% of q1's length
    assert s["assembly_major_telomeric_ends"] == 1
    assert s["assembly_minor_sequences_with_telomere"] == 0
    assert s["unaligned_contigs"] == 1
    assert s["unaligned_contig_bp"] == len(unaligned)
    assert s["unaligned_contig_n_pct"] == pytest.approx(10 / 30 * 100)
    assert s["assembly_softmasked_pct"] is not None
    assert s["unaligned_contig_softmasked_pct"] == pytest.approx(20 / 30 * 100)
    rows = {row["query_contig"]: row for row in _read_tsv(tmp_path / "out" / "s1.by_query_contig.tsv")}
    assert set(rows) == {"q1", "unaligned"}  # every assembly contig, aligned or not
    assert rows["q1"]["telomere_start"] == "true" and rows["q1"]["telomere_end"] == "false"
    assert rows["q1"]["major"] == "true" and rows["unaligned"]["major"] == "false"
    assert rows["q1"]["contig_pieces"] == "2"
    assert rows["unaligned"]["aligned"] == "false"
    assert rows["unaligned"]["unaligned_start_bp"] == "NA"
    assert rows["q1"]["unaligned_start_bp"] == str(len(TELOMERE))


def test_without_fasta_assembly_columns_are_na(tmp_path):
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 10})
    maf = _write_maf(tmp_path / "s1.maf", [[("chr1", 0, "+", 10, "ACGTACGTAC"), ("q1", 0, "+", 10, "ACGTACGTAC")]])
    run(maf, fai, "s1", tmp_path / "out", dotplots="false")
    summary = _read_tsv(tmp_path / "out" / "s1.maf_stats.tsv")[0]
    assert all(summary[c] == "NA" for c in summary if c.startswith(("assembly_", "unaligned_")))
    assert summary["breakpoints_near_n_gap"] == "NA"
    query_row = _read_tsv(tmp_path / "out" / "s1.by_query_contig.tsv")[0]
    assert query_row["n_bp"] == "NA" and query_row["telomere_start"] == "NA"


def test_fasta_length_must_match_maf_srcsize(tmp_path):
    fasta = _write_fasta(tmp_path / "s1.fa", {"q1": "ACGT" * 3})
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 10})
    maf = _write_maf(tmp_path / "s1.maf", [[("chr1", 0, "+", 10, "ACGT"), ("q1", 0, "+", 10, "ACGT")]])
    with pytest.raises(MafValidationError, match="query FASTA length 12"):
        run(maf, fai, "s1", tmp_path / "out", fasta=fasta, dotplots="false")


def test_fasta_and_query_fai_are_mutually_exclusive(tmp_path):
    with pytest.raises(ValueError, match="mutually exclusive"):
        run(tmp_path / "x.maf", tmp_path / "ref.fai", "s1", tmp_path / "out",
            query_fai=tmp_path / "q.fai", fasta=tmp_path / "q.fa")


# ---------------------------------------------------------------- breakpoint locations


def _breakpoint_maf(tmp_path: Path, q_length: int, junction: int) -> tuple[Path, Path]:
    """Two blocks on q1 joined at ``junction`` with a strand flip."""
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 1000})
    maf = _write_maf(tmp_path / "s1.maf", [
        [("chr1", 100, "+", 1000, "ACGTACGTAC"), ("q1", junction - 10, "+", q_length, "ACGTACGTAC")],
        # minus-strand query row: forward interval [junction, junction + 10)
        [("chr1", 500, "+", 1000, "ACGTACGTAC"), ("q1", q_length - junction - 10, "-", q_length, "ACGTACGTAC")],
    ])
    return maf, fai


def test_breakpoint_location_and_junction_reference_positions(tmp_path):
    maf, fai = _breakpoint_maf(tmp_path, q_length=1000, junction=500)
    stats = run(maf, fai, "s1", tmp_path / "out", dotplots="false", breakpoint_context_bp=50)
    (bp,) = stats.breakpoints
    assert bp["location"] == "interior"
    assert bp["distance_to_contig_end_bp"] == 500
    assert (bp["query_junction_start"], bp["query_junction_end"]) == (500, 500)
    assert bp["strand_flip"] is True and bp["near_n_gap"] is None
    # + block ends at ref 110; the - block's forward-query start sits opposite its ref end (510)
    assert (bp["left_reference_contig"], bp["left_reference_pos"], bp["left_strand"]) == ("chr1", 110, "+")
    assert (bp["right_reference_contig"], bp["right_reference_pos"], bp["right_strand"]) == ("chr1", 510, "-")
    assert stats.summary["breakpoints_interior"] == 1
    rows = _read_tsv(tmp_path / "out" / "s1.breakpoints.tsv")
    assert len(rows) == 1 and list(rows[0]) == maf_stats.BREAKPOINT_COLUMNS


def test_breakpoint_near_contig_end(tmp_path):
    maf, fai = _breakpoint_maf(tmp_path, q_length=1000, junction=40)
    stats = run(maf, fai, "s1", tmp_path / "out", dotplots="false", breakpoint_context_bp=50)
    assert stats.breakpoints[0]["location"] == "contig_end"
    assert stats.summary["breakpoints_near_contig_end"] == 1


def test_breakpoint_near_n_gap_needs_fasta(tmp_path):
    q = "A" * 480 + "N" * 15 + "A" * 505
    maf, fai = _breakpoint_maf(tmp_path, q_length=len(q), junction=500)
    fasta = _write_fasta(tmp_path / "s1.fa", {"q1": q})
    stats = run(maf, fai, "s1", tmp_path / "out", fasta=fasta, dotplots="false", breakpoint_context_bp=50)
    assert stats.breakpoints[0]["location"] == "n_gap"
    assert stats.breakpoints[0]["near_n_gap"] is True
    assert stats.summary["breakpoints_near_n_gap"] == 1


def test_empty_breakpoints_tsv_has_header(tmp_path):
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 10})
    maf = _write_maf(tmp_path / "s1.maf", [[("chr1", 0, "+", 10, "ACGTACGTAC"), ("q1", 0, "+", 10, "ACGTACGTAC")]])
    run(maf, fai, "s1", tmp_path / "out", dotplots="false")
    lines = (tmp_path / "out" / "s1.breakpoints.tsv").read_text().splitlines()
    assert lines == ["\t".join(maf_stats.BREAKPOINT_COLUMNS)]


def test_scattered_telomere_copies_do_not_make_a_telomeric_end():
    scattered = ("TTTAGGG" + "ACGTACGTACGTACGTACGT") * 40  # ~26% repeat in the last 1 kb
    contig = maf_stats.contig_stats((scattered + "ACGT" * 1000 + "TTTAGGG" * 150).encode())
    assert contig.telomere_start is False
    assert contig.telomere_end is True


def test_unmasked_fasta_reports_softmasking_as_na(tmp_path):
    fasta = _write_fasta(tmp_path / "s1.fa", {"q1": "ACGTACGTAC"})
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 10})
    maf = _write_maf(tmp_path / "s1.maf", [[("chr1", 0, "+", 10, "ACGTACGTAC"), ("q1", 0, "+", 10, "ACGTACGTAC")]])
    stats = run(maf, fai, "s1", tmp_path / "out", fasta=fasta, dotplots="false")
    assert stats.summary["assembly_softmasked_pct"] is None
    assert stats.summary["unaligned_contig_softmasked_pct"] is None
    assert stats.by_query[0]["softmasked_pct"] is None


# ---------------------------------------------------------------- nested blocks


def test_block_inside_another_blocks_query_span_is_nested_not_a_breakpoint(tmp_path):
    """The CML103 chr1 pattern: a small reversed block, aligned far away on the
    reference, whose query interval lies inside a large forward block's."""
    fai = _write_fai(tmp_path / "ref.fai", {"chr1": 2000})
    big_ref, big_query = "ACGT" * 250, "ACGT" * 100 + "-" * 0 + "ACGT" * 150
    maf = _write_maf(tmp_path / "s1.maf", [
        [("chr1", 0, "+", 2000, big_ref), ("q1", 0, "+", 1000, big_query)],
        # query forward interval [400, 410) on the minus strand, reference 1500-1510
        [("chr1", 1500, "+", 2000, "ACGTACGTAC"), ("q1", 590, "-", 1000, "ACGTACGTAC")],
    ])
    stats = run(maf, fai, "s1", tmp_path / "out", dotplots="false")
    assert stats.summary["breakpoint_adjacencies"] == 0
    assert stats.summary["strand_flips"] == 0
    assert stats.summary["nested_blocks"] == 1
    assert stats.summary["nested_block_reference_bp"] == 10
    assert stats.summary["overlapping_query_bp"] == 10
    (row,) = stats.nested
    assert (row["query_start"], row["query_end"], row["strand"]) == (400, 410, "-")
    assert (row["reference_contig"], row["reference_start"], row["reference_end"]) == ("chr1", 1500, 1510)
    assert (row["container_query_start"], row["container_query_end"]) == (0, 1000)
    rows = _read_tsv(tmp_path / "out" / "s1.nested_blocks.tsv")
    assert len(rows) == 1 and list(rows[0]) == maf_stats.NESTED_COLUMNS
    by_query = _read_tsv(tmp_path / "out" / "s1.by_query_contig.tsv")[0]
    assert by_query["nested_blocks"] == "1" and by_query["overlapping_query_bp"] == "10"


def test_partial_query_overlap_is_still_an_adjacency():
    blocks = [_block("chr1", 0, 100, 0, 100), _block("chr1", 500, 600, 90, 190, "-")]
    counts = classify_breakpoints(blocks, ORDER, 0, 0)
    assert counts.nested == []
    assert counts.strand_flips == 1


def test_identical_query_intervals_nest_the_later_block():
    a = _block("chr1", 0, 10, 0, 10)
    b = _block("chr2", 0, 10, 0, 10)
    counts = classify_breakpoints([b, a], ORDER, 0, 0)
    assert [(n.ref_contig, c.ref_contig) for n, c in counts.nested] == [("chr2", "chr1")]
    assert counts.breakpoint_adjacencies == 0


def test_major_sequences_ignore_scaffold_debris():
    contigs = {f"chr{i}": 300 - i * 10 for i in range(1, 11)}  # 290..200
    contigs.update({f"scaf_{i}": 5 for i in range(200)})  # lots of small unplaced sequence
    assembly = {name: maf_stats.contig_stats(b"A" * length) for name, length in contigs.items()}
    assert maf_stats.major_sequences(assembly) == {f"chr{i}" for i in range(1, 11)}

