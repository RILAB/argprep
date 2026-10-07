import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.maf_to_fasta import (
    build_columns,
    load_mask,
    load_sample_alignment,
    merge_insertion,
)
from scripts.maf_to_sites import read_contig_length, read_contig_region

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "maf_to_fasta.py"


@pytest.mark.parametrize("reverse", [False, True])
def test_insertion_conflicts_with_explicit_absence(tmp_path, reverse):
    blocks = [("AC--GT", "s1", "ACTTGT"), ("ACGT", "s1", "ACGT")]
    maf = tmp_path / "sample.maf"
    _write_maf(maf, "chr1", blocks[::-1] if reverse else blocks)
    anchors, insertions = load_sample_alignment(maf, "chr1", 0, 4)
    assert list(anchors) == [1, 2, 3, 4]
    assert insertions == {1: "??"}


@pytest.mark.parametrize("window_end", [2, 4])
def test_leading_block_insertion_kept_with_anchor_in_window(tmp_path, window_end):
    maf = tmp_path / "sample.maf"
    _write_maf_rows(maf, [("chr1", 2, "+", 4, "--GT"),
                          ("s1", 0, "+", 4, "TTGT")])
    _, insertions = load_sample_alignment(maf, "chr1", 0, window_end)
    assert insertions == {1: "TT"}
    _, insertions = load_sample_alignment(maf, "chr1", 2, 4)
    assert insertions == {}


@pytest.mark.parametrize("reverse", [False, True])
def test_block_boundary_does_not_establish_insertion_absence(tmp_path, reverse):
    maf = tmp_path / "sample.maf"
    blocks = [
        "a\ns chr1 0 2 + 4 AC\ns s1 0 2 + 6 AC\n\n",
        "a\ns chr1 2 2 + 4 GT\ns s1 4 2 + 6 GT\n\n",
        "a\ns chr1 0 4 + 4 AC--GT\ns s1 0 6 + 6 ACTTGT\n\n",
    ]
    maf.write_text("".join(blocks[::-1] if reverse else blocks))
    _, insertions = load_sample_alignment(maf, "chr1", 0, 4)
    assert insertions == {1: "TT"}


def test_reference_span_cap_precedes_sequence_and_alignment_loading(monkeypatch):
    import scripts.maf_to_fasta as module

    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "--maf-dir", "unused",
        "--samples", "s1", "--reference-fasta", "unused.fa", "--contig", "chr1",
        "--max-columns", "5", "--out", "unused.fa"])
    monkeypatch.setattr(module, "read_contig_length", lambda *args: 100)

    def unexpected(*args, **kwargs):
        pytest.fail("Oversized reference span must be rejected before loading")

    monkeypatch.setattr(module, "read_contig_region", unexpected)
    monkeypatch.setattr(module, "load_sample_alignment", unexpected)
    with pytest.raises(ValueError, match="Reference window alone.*over --max-columns"):
        module.main()


def _ungapped(text: str) -> str:
    return text.replace("-", "")


def _write_reference(tmp_path: Path, contig: str, seq: str, line_width: int = 0) -> Path:
    """Write a FASTA plus a matching .fai. ``line_width`` 0 writes one line."""
    width = line_width or len(seq)
    lines = [seq[i:i + width] for i in range(0, len(seq), width)]
    ref = tmp_path / "ref.fa"
    ref.write_text(f">{contig}\n" + "".join(f"{line}\n" for line in lines), encoding="utf-8")
    offset = len(contig) + 2
    (tmp_path / "ref.fa.fai").write_text(
        f"{contig}\t{len(seq)}\t{offset}\t{width}\t{width + 1}\n", encoding="utf-8"
    )
    return ref


def _write_maf(path: Path, contig: str, blocks: list[tuple[str, str, str]]) -> None:
    """Write a pairwise MAF. Each block is (ref_text, sample_name, sample_text),
    both rows starting at reference/sample offset 0."""
    parts = ["##maf version=1\n"]
    for ref_text, sample, sample_text in blocks:
        ref_len = len(_ungapped(ref_text))
        sample_len = len(_ungapped(sample_text))
        parts.append("a score=0\n")
        parts.append(f"s {contig} 0 {ref_len} + {ref_len} {ref_text}\n")
        parts.append(f"s {sample} 0 {sample_len} + {sample_len} {sample_text}\n")
        parts.append("\n")
    path.write_text("".join(parts), encoding="utf-8")


def _write_maf_rows(path: Path, rows: list[tuple[str, int, str, int, str]]) -> None:
    """Write a single-block pairwise MAF from explicit
    (src, start, strand, src_size, text) rows."""
    parts = ["##maf version=1\n", "a score=0\n"]
    for src, start, strand, src_size, text in rows:
        parts.append(f"s {src} {start} {len(_ungapped(text))} {strand} {src_size} {text}\n")
    path.write_text("".join(parts), encoding="utf-8")


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
            continue
        chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks)
    return records


def _run_script(tmp_path: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )


def test_insertion_adds_columns_that_other_samples_gap(tmp_path: Path) -> None:
    ref_seq = "ACGTACGT"
    _write_reference(tmp_path, "chr1", ref_seq)
    maf_dir = tmp_path / "mafs"
    maf_dir.mkdir()
    _write_maf(maf_dir / "s1.maf", "chr1", [("ACGT--ACGT", "s1", "ACGTGGACGT")])
    _write_maf(maf_dir / "s2.maf", "chr1", [("ACGTACGT", "s2", "ACGTACGT")])

    _run_script(
        tmp_path,
        "--maf-dir", "mafs",
        "--reference-fasta", "ref.fa",
        "--contig", "chr1",
        "--out", "win.fa",
    )
    records = _read_fasta(tmp_path / "win.fa")
    assert records == {
        "REF": "ACGT--ACGT",
        "s1": "ACGTGGACGT",
        "s2": "ACGT--ACGT",
    }


def test_insertions_of_unequal_length_pad_to_the_widest(tmp_path: Path) -> None:
    _write_reference(tmp_path, "chr1", "ACGT")
    maf_dir = tmp_path / "mafs"
    maf_dir.mkdir()
    # Both insert after reference base 2 (0-based index 1), but 3 bp vs 1 bp.
    _write_maf(maf_dir / "s1.maf", "chr1", [("AC---GT", "s1", "ACTTTGT")])
    _write_maf(maf_dir / "s2.maf", "chr1", [("AC-GT", "s2", "ACAGT")])

    _run_script(
        tmp_path,
        "--maf-dir", "mafs",
        "--reference-fasta", "ref.fa",
        "--contig", "chr1",
        "--out", "win.fa",
    )
    records = _read_fasta(tmp_path / "win.fa")
    # Shared columns are left-aligned; the narrower insertion pads on the right.
    assert records["REF"] == "AC---GT"
    assert records["s1"] == "ACTTTGT"
    assert records["s2"] == "ACA--GT"


def test_deletion_is_a_gap_and_unaligned_reference_is_n(tmp_path: Path) -> None:
    _write_reference(tmp_path, "chr1", "ACGTACGT")
    maf_dir = tmp_path / "mafs"
    maf_dir.mkdir()
    # Deletes reference base 5 and stops aligning after base 6.
    _write_maf(maf_dir / "s1.maf", "chr1", [("ACGTAC", "s1", "ACGT-C")])

    _run_script(
        tmp_path,
        "--maf-dir", "mafs",
        "--reference-fasta", "ref.fa",
        "--contig", "chr1",
        "--out", "win.fa",
    )
    records = _read_fasta(tmp_path / "win.fa")
    assert records["REF"] == "ACGTACGT"
    assert records["s1"] == "ACGT-CNN"


def test_insertion_before_the_first_window_base_is_dropped(tmp_path: Path) -> None:
    _write_reference(tmp_path, "chr1", "ACGTACGT")
    maf_dir = tmp_path / "mafs"
    maf_dir.mkdir()
    # Insertion sits between reference bases 4 and 5, i.e. before the window.
    _write_maf(maf_dir / "s1.maf", "chr1", [("ACGT--ACGT", "s1", "ACGTGGACGT")])

    _run_script(
        tmp_path,
        "--maf-dir", "mafs",
        "--reference-fasta", "ref.fa",
        "--contig", "chr1",
        "--start", "5", "--end", "8",
        "--out", "win.fa",
    )
    records = _read_fasta(tmp_path / "win.fa")
    assert records == {"REF": "ACGT", "s1": "ACGT"}


def test_insertion_anchored_to_the_last_window_base_is_kept(tmp_path: Path) -> None:
    _write_reference(tmp_path, "chr1", "ACGTACGT")
    maf_dir = tmp_path / "mafs"
    maf_dir.mkdir()
    _write_maf(maf_dir / "s1.maf", "chr1", [("ACGT--ACGT", "s1", "ACGTGGACGT")])

    _run_script(
        tmp_path,
        "--maf-dir", "mafs",
        "--reference-fasta", "ref.fa",
        "--contig", "chr1",
        "--start", "1", "--end", "4",
        "--out", "win.fa",
    )
    records = _read_fasta(tmp_path / "win.fa")
    assert records == {"REF": "ACGT--", "s1": "ACGTGG"}


def test_conflicting_blocks_degrade_anchors_and_insertions(tmp_path: Path) -> None:
    _write_reference(tmp_path, "chr1", "ACGT")
    maf_dir = tmp_path / "mafs"
    maf_dir.mkdir()
    # Two blocks cover the same reference span with different calls at base 4
    # and different insertions after base 2.
    _write_maf(
        maf_dir / "s1.maf",
        "chr1",
        [("AC--GT", "s1", "ACTTGT"), ("AC--GT", "s1", "ACAAGA")],
    )

    _run_script(
        tmp_path,
        "--maf-dir", "mafs",
        "--reference-fasta", "ref.fa",
        "--contig", "chr1",
        "--out", "win.fa",
    )
    records = _read_fasta(tmp_path / "win.fa")
    assert records["REF"] == "AC--GT"
    assert records["s1"] == "ACNNGN"


def test_quality_mask_renders_low_scoring_bases_as_n(tmp_path: Path) -> None:
    _write_reference(tmp_path, "chr1", "ACGT")
    maf_dir = tmp_path / "mafs"
    maf_dir.mkdir()
    _write_maf(maf_dir / "s1.maf", "chr1", [("AC-GT", "s1", "ACAGT")])
    qdir = tmp_path / "quality"
    qdir.mkdir()
    # Sample coordinates 2 and 3 cover the inserted base and reference base 3.
    (qdir / "s1.bed").write_text("s1\t2\t4\t0.1\n", encoding="utf-8")

    _run_script(
        tmp_path,
        "--maf-dir", "mafs",
        "--reference-fasta", "ref.fa",
        "--contig", "chr1",
        "--quality-bed-dir", "quality",
        "--quality-min", "0.5",
        "--out", "win.fa",
    )
    records = _read_fasta(tmp_path / "win.fa")
    assert records["REF"] == "AC-GT"
    assert records["s1"] == "ACNNT"


def test_quality_mask_uses_minus_strand_sample_coordinates(tmp_path: Path) -> None:
    _write_reference(tmp_path, "chr1", "ACGT")
    maf_dir = tmp_path / "mafs"
    maf_dir.mkdir()
    # Sample row on the minus strand: text position i maps to forward-strand
    # coordinate src_size - start - 1 - i = 3 - i.
    _write_maf_rows(
        maf_dir / "s1.maf",
        [("chr1", 0, "+", 4, "ACGT"), ("s1", 0, "-", 4, "ACGT")],
    )
    qdir = tmp_path / "quality"
    qdir.mkdir()
    (qdir / "s1.bed").write_text("s1\t3\t4\t0.1\n", encoding="utf-8")

    _run_script(
        tmp_path,
        "--maf-dir", "mafs",
        "--reference-fasta", "ref.fa",
        "--contig", "chr1",
        "--quality-bed-dir", "quality",
        "--quality-min", "0.5",
        "--out", "win.fa",
    )
    records = _read_fasta(tmp_path / "win.fa")
    assert records["s1"] == "NCGT"


def test_mask_bed_blanks_anchors_and_drops_their_insertions(tmp_path: Path) -> None:
    _write_reference(tmp_path, "chr1", "ACGT")
    maf_dir = tmp_path / "mafs"
    maf_dir.mkdir()
    _write_maf(maf_dir / "s1.maf", "chr1", [("AC-GT", "s1", "ACAGT")])
    # Mask reference base 2 (0-based [1, 2)), which anchors the insertion.
    (tmp_path / "mask.bed").write_text("chr1\t1\t2\n", encoding="utf-8")

    _run_script(
        tmp_path,
        "--maf-dir", "mafs",
        "--reference-fasta", "ref.fa",
        "--contig", "chr1",
        "--mask-bed", "mask.bed",
        "--out", "win.fa",
    )
    records = _read_fasta(tmp_path / "win.fa")
    assert records["REF"] == "ACGT"
    assert records["s1"] == "ANGT"


def test_chunk_root_discovers_samples_and_gzipped_chunks(tmp_path: Path) -> None:
    _write_reference(tmp_path, "chr1", "ACGT")
    root = tmp_path / "maf_by_contig"
    for sample, sample_text in (("s1", "ACGT"), ("s2", "ATGT")):
        sample_dir = root / sample
        sample_dir.mkdir(parents=True)
        _write_maf(sample_dir / "chr1.maf", "chr1", [("ACGT", sample, sample_text)])

    result = _run_script(
        tmp_path,
        "--maf-chunk-root", "maf_by_contig",
        "--reference-fasta", "ref.fa",
        "--contig", "chr1",
        "--out", "win.fa",
    )
    records = _read_fasta(tmp_path / "win.fa")
    assert records == {"REF": "ACGT", "s1": "ACGT", "s2": "ATGT"}
    assert "3 sequences" in result.stderr


def test_no_reference_and_line_width_control_output_shape(tmp_path: Path) -> None:
    _write_reference(tmp_path, "chr1", "ACGTACGTAC")
    maf_dir = tmp_path / "mafs"
    maf_dir.mkdir()
    _write_maf(maf_dir / "s1.maf", "chr1", [("ACGTACGTAC", "s1", "ACGTACGTAC")])

    _run_script(
        tmp_path,
        "--maf-dir", "mafs",
        "--reference-fasta", "ref.fa",
        "--contig", "chr1",
        "--no-reference",
        "--line-width", "4",
        "--out", "win.fa",
    )
    lines = (tmp_path / "win.fa").read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith(">s1 chr1:1-10")
    assert lines[1:] == ["ACGT", "ACGT", "AC"]


def test_max_columns_guard_refuses_wide_alignments(tmp_path: Path) -> None:
    _write_reference(tmp_path, "chr1", "ACGT")
    maf_dir = tmp_path / "mafs"
    maf_dir.mkdir()
    _write_maf(maf_dir / "s1.maf", "chr1", [("AC----GT", "s1", "ACTTTTGT")])

    with pytest.raises(subprocess.CalledProcessError) as excinfo:
        _run_script(
            tmp_path,
            "--maf-dir", "mafs",
            "--reference-fasta", "ref.fa",
            "--contig", "chr1",
            "--max-columns", "5",
            "--out", "win.fa",
        )
    assert "over --max-columns 5" in excinfo.value.stderr


@pytest.mark.parametrize(
    "start, end, message",
    [
        (0, 4, "--start must be >= 1"),
        (1, 99, "exceeds length 4"),
        (3, 2, "--end must be >= --start"),
    ],
)
def test_region_bounds_are_validated(tmp_path: Path, start: int, end: int, message: str) -> None:
    _write_reference(tmp_path, "chr1", "ACGT")
    maf_dir = tmp_path / "mafs"
    maf_dir.mkdir()
    _write_maf(maf_dir / "s1.maf", "chr1", [("ACGT", "s1", "ACGT")])

    with pytest.raises(subprocess.CalledProcessError) as excinfo:
        _run_script(
            tmp_path,
            "--maf-dir", "mafs",
            "--reference-fasta", "ref.fa",
            "--contig", "chr1",
            "--start", str(start), "--end", str(end),
            "--out", "win.fa",
        )
    assert message in excinfo.value.stderr


def test_quality_options_require_each_other(tmp_path: Path) -> None:
    _write_reference(tmp_path, "chr1", "ACGT")
    maf_dir = tmp_path / "mafs"
    maf_dir.mkdir()
    _write_maf(maf_dir / "s1.maf", "chr1", [("ACGT", "s1", "ACGT")])

    with pytest.raises(subprocess.CalledProcessError) as excinfo:
        _run_script(
            tmp_path,
            "--maf-dir", "mafs",
            "--reference-fasta", "ref.fa",
            "--contig", "chr1",
            "--quality-min", "0.5",
            "--out", "win.fa",
        )
    assert "--quality-min requires --quality-bed-dir" in excinfo.value.stderr


def test_load_sample_alignment_left_anchors_insertions(tmp_path: Path) -> None:
    maf = tmp_path / "s1.maf"
    _write_maf(maf, "chr1", [("AC--GT", "s1", "ACTTGT")])
    anchors, insertions = load_sample_alignment(maf, "chr1", 0, 4)
    assert list(anchors) == [1, 2, 3, 4]
    assert insertions == {1: "TT"}


def test_merge_insertion_degrades_conflicts_to_the_longer_width() -> None:
    insertions: dict[int, str] = {}
    merge_insertion(insertions, 3, "AC")
    merge_insertion(insertions, 3, "AC")
    assert insertions == {3: "AC"}
    merge_insertion(insertions, 3, "GGGG")
    assert insertions == {3: "????"}


def test_build_columns_uses_the_widest_insertion_per_anchor() -> None:
    widths, offsets, total = build_columns([{1: "AC"}, {1: "A", 2: "TTT"}], 4)
    assert widths == [0, 2, 3, 0]
    assert offsets == [0, 1, 4, 8]
    assert total == 9


def test_load_mask_clips_intervals_to_the_window(tmp_path: Path) -> None:
    bed = tmp_path / "mask.bed"
    bed.write_text("chr1\t2\t6\nchr2\t0\t10\n", encoding="utf-8")
    assert list(load_mask(bed, "chr1", 3, 8)) == [1, 1, 1, 0, 0]


def test_load_mask_rejects_malformed_rows(tmp_path: Path) -> None:
    bed = tmp_path / "mask.bed"
    bed.write_text("chr1\tnot_a_number\t6\n", encoding="utf-8")
    with pytest.raises(ValueError, match="must be integers"):
        load_mask(bed, "chr1", 0, 8)


def test_read_contig_region_matches_a_wrapped_fasta(tmp_path: Path) -> None:
    seq = "ACGTACGTACGTACGTACGT"
    ref = _write_reference(tmp_path, "chr1", seq, line_width=6)
    assert read_contig_length(ref, "chr1") == len(seq)
    for start, end in ((0, len(seq)), (0, 1), (5, 7), (6, 12), (13, 20)):
        assert read_contig_region(ref, "chr1", start, end) == seq[start:end]


def test_read_contig_region_falls_back_without_a_fai(tmp_path: Path) -> None:
    seq = "ACGTACGTAC"
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGTA\nCGTAC\n", encoding="utf-8")
    assert read_contig_length(ref, "chr1") == len(seq)
    assert read_contig_region(ref, "chr1", 3, 8) == seq[3:8]
    with pytest.raises(ValueError, match="exceeds length"):
        read_contig_region(ref, "chr1", 3, 99)


def test_minus_strand_reference_row_is_rejected(tmp_path: Path) -> None:
    maf = tmp_path / "s1.maf"
    _write_maf_rows(maf, [("chr1", 0, "-", 10, "ACGT"), ("s1", 0, "+", 10, "ACGT")])
    with pytest.raises(ValueError, match="reference rows are supported"):
        load_sample_alignment(maf, "chr1", 0, 10)
