import subprocess
import sys
import gzip
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.maf_to_sites import (
    choose_sample_record,
    discover_samples,
    iter_maf_blocks,
    load_quality_mask,
    load_sample_calls,
    maf_path_for_sample,
    maf_path_for_sample_with_map,
    missing_threshold,
    parse_maf_path_map,
    read_contig_sequence,
)
from scripts.summary_report import read_sample_missing_bp
from scripts.split_maf_by_contig import main as split_maf_main


def _run(cmd, cwd):
    subprocess.run(cmd, cwd=cwd, check=True)


def _write_pairwise_maf(path: Path, contig: str, ref_seq: str, sample: str, sample_seq: str) -> None:
    path.write_text(
        "##maf version=1\n"
        "a score=0\n"
        f"s {contig} 0 {len(ref_seq.replace('-', ''))} + {len(ref_seq.replace('-', ''))} {ref_seq}\n"
        f"s {sample} 0 {len(sample_seq.replace('-', ''))} + {len(sample_seq.replace('-', ''))} {sample_seq}\n",
        encoding="utf-8",
    )


def _write_pairwise_maf_strand(
    path: Path,
    contig: str,
    ref_seq: str,
    sample: str,
    sample_seq: str,
    strand: str,
    start: int,
    src_size: int,
) -> None:
    path.write_text(
        "##maf version=1\n"
        "a score=0\n"
        f"s {contig} 0 {len(ref_seq.replace('-', ''))} + {len(ref_seq.replace('-', ''))} {ref_seq}\n"
        f"s {sample} {start} {len(sample_seq.replace('-', ''))} {strand} {src_size} {sample_seq}\n",
        encoding="utf-8",
    )


def _write_same_src_pairwise_maf(path: Path, contig: str, ref_seq: str, sample_seq: str) -> None:
    path.write_text(
        "##maf version=1\n"
        "a score=0\n"
        f"s {contig} 0 {len(ref_seq.replace('-', ''))} + {len(ref_seq.replace('-', ''))} {ref_seq}\n"
        f"s {contig} 0 {len(sample_seq.replace('-', ''))} + {len(sample_seq.replace('-', ''))} {sample_seq}\n",
        encoding="utf-8",
    )


def _read_vcf_records(path: Path) -> list[list[str]]:
    records = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            records.append(line.rstrip("\n").split("\t"))
    return records


def _read_vcf_header_line(path: Path) -> str:
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                return line.rstrip("\n")
    raise AssertionError(f"No VCF header line found in {path}")


def _read_bed(path: Path) -> list[tuple]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if len(fields) >= 4:
                rows.append((chrom, start, end, fields[3]))
            else:
                rows.append((chrom, start, end))
    return rows


def _read_sites(path: Path):
    """Parse an ARGweaver .sites file into (names, region, [(pos, alleles)...])."""
    names: list[str] = []
    region: tuple[str, ...] = ()
    rows: list[tuple[str, str]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if fields[0] == "NAMES":
                names = fields[1:]
            elif fields[0] == "REGION":
                region = tuple(fields[1:])
            else:
                rows.append((fields[0], fields[1]))
    return names, region, rows


def test_maf_to_sites_emits_expected_records_and_mask(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGTACGT", "s1", "ACGTACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGTACGT", "s2", "AT-TCCGN")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--window-bp",
            "4",
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
            "--mask-indel-adjacent-snps",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".vcf")
    masked = Path(str(out_prefix) + ".mask.bed")

    all_records = _read_vcf_records(all_sites)
    variant_records = _read_vcf_records(variants)
    bed = _read_bed(masked)

    assert [record[1] for record in all_records] == ["1", "4", "5", "6", "7"]
    assert [record[1] for record in variant_records] == ["5"]
    assert bed == [("chr1", 1, 3), ("chr1", 7, 8)]
    assert Path(str(out_prefix) + ".report_stats.tsv").read_text(encoding="utf-8") == (
        "record_type\tcontig\tstart\tend\tsample\tinvariant\tvariant\tmasked\tcalled\tcarried_variant\n"
        "window\tchr1\t0\t4\t\t2\t0\t2\t\t\n"
        "window\tchr1\t4\t8\t\t2\t1\t1\t\t\n"
        "sample\tchr1\t\t\ts1\t\t\t\t5\t0\n"
        "sample\tchr1\t\t\ts2\t\t\t\t5\t1\n"
    )


def test_maf_to_sites_emits_per_sample_missing_masks(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    # s1: fully aligned, no missing data
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGTACGT", "s1", "ACGTACGT")
    # s2: deletion at ref pos 2 (idx 2), N at ref pos 7 (idx 7) treated as missing
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGTACGT", "s2", "AT-TCCGN")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(out_prefix),
            "--samples", "s1", "s2",
            "--max-missing-count", "0",
        ],
        cwd=Path.cwd(),
    )

    s1_mask = _read_bed(Path(str(out_prefix) + ".s1.missing.bed"))
    s2_mask = _read_bed(Path(str(out_prefix) + ".s2.missing.bed"))

    assert s1_mask == []
    assert s2_mask == [("chr1", 2, 3, "s2"), ("chr1", 7, 8, "s2")]


def test_maf_to_sites_preserves_explicit_sample_order(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nAC\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t2\t6\t2\t3\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "AC", "s1", "AT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "AC", "s2", "AC")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(out_prefix),
            "--samples", "s2", "s1",
            "--max-missing-count", "0",
            "--keep-indel-adjacent-snps",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    assert _read_vcf_header_line(all_sites).endswith("\ts2\ts1")
    records = _read_vcf_records(all_sites)
    assert [record[1] for record in records] == ["1", "2"]
    assert records[0][9:] == ["0", "0"]
    assert records[1][9:] == ["0", "1"]


def test_maf_to_sites_per_sample_missing_mask_includes_unaligned_positions(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    # s1: only covers the first base; positions 1-3 have no alignment (code 0)
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "A", "s1", "A")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "ACGT")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(out_prefix),
            "--samples", "s1", "s2",
            "--max-missing-count", "0",
        ],
        cwd=Path.cwd(),
    )

    s1_mask = _read_bed(Path(str(out_prefix) + ".s1.missing.bed"))
    s2_mask = _read_bed(Path(str(out_prefix) + ".s2.missing.bed"))

    assert s1_mask == [("chr1", 1, 4, "s1")]  # positions 1, 2, 3 unaligned → merged interval
    assert s2_mask == []


def test_maf_to_sites_keeps_multiallelic_sites_by_default(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nA\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t1\t6\t1\t2\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "A", "s1", "A")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "A", "s2", "C")
    _write_pairwise_maf(maf_dir / "s3.maf", "chr1", "A", "s3", "G")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "s3",
            "--max-missing-count",
            "0",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".vcf")
    masked = Path(str(out_prefix) + ".mask.bed")

    all_records = _read_vcf_records(all_sites)
    variant_records = _read_vcf_records(variants)

    assert len(all_records) == 1
    assert len(variant_records) == 1
    assert all_records[0][3] == "A"
    assert all_records[0][4] == "C,G"
    assert all_records[0][9:] == ["0", "1", "2"]
    assert variant_records == all_records
    assert _read_bed(masked) == []


def test_maf_to_sites_masks_multiallelic_when_flagged(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nA\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t1\t6\t1\t2\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "A", "s1", "A")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "A", "s2", "C")
    _write_pairwise_maf(maf_dir / "s3.maf", "chr1", "A", "s3", "G")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "s3",
            "--max-missing-count",
            "0",
            "--mask-multiallelic-snps",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".vcf")
    masked = Path(str(out_prefix) + ".mask.bed")
    summary = Path(str(out_prefix) + ".site_summary.tsv").read_text(encoding="utf-8")

    # The tri-allelic site (A/C/G) is masked out rather than emitted.
    assert _read_vcf_records(all_sites) == []
    assert _read_vcf_records(variants) == []
    assert _read_bed(masked) == [("chr1", 0, 1)]
    assert "masked_multiallelic\t1" in summary


def test_maf_to_sites_rejects_zero_length_contig(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    # chr1 has sequence; chr0 is an empty record (length 0).
    ref.write_text(">chr0\n\n>chr1\nA\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr0\t0\t6\t0\t1\nchr1\t1\t14\t1\t2\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "A", "s1", "A")

    out_prefix = tmp_path / "results" / "combined.chr0"
    proc = subprocess.run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr0",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "--max-missing-count",
            "0",
        ],
        cwd=Path.cwd(),
        text=True,
        capture_output=True,
    )
    assert proc.returncode != 0
    assert "Reference contig 'chr0' has length 0" in proc.stderr


def test_maf_to_sites_emits_argweaver_sites(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "AGGT")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(out_prefix),
            "--samples", "s1", "s2",
            "--max-missing-count", "0",
            "--emit-argweaver-sites",
        ],
        cwd=Path.cwd(),
    )

    names, region, rows = _read_sites(Path(str(out_prefix) + ".sites"))
    assert names == ["s1", "s2"]
    assert region == ("chr1", "1", "4")
    # Only the single variant site (pos 2: ref C, s2 G) is emitted, one real base per sample.
    assert rows == [("2", "CG")]
    # The .sites variant positions match the variant VCF exactly.
    variant_positions = [rec[1] for rec in _read_vcf_records(Path(str(out_prefix) + ".vcf"))]
    assert [pos for pos, _ in rows] == variant_positions


def test_maf_to_sites_omits_argweaver_sites_by_default(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "AGGT")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(out_prefix),
            "--samples", "s1", "s2",
            "--max-missing-count", "0",
        ],
        cwd=Path.cwd(),
    )

    assert not Path(str(out_prefix) + ".sites").exists()


def test_maf_to_sites_argweaver_sites_missing_becomes_n(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "AGGT")
    _write_pairwise_maf(maf_dir / "s3.maf", "chr1", "ACGT", "s3", "ANGT")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(out_prefix),
            "--samples", "s1", "s2", "s3",
            "--max-missing-count", "1",
            "--emit-argweaver-sites",
        ],
        cwd=Path.cwd(),
    )

    names, _region, rows = _read_sites(Path(str(out_prefix) + ".sites"))
    assert names == ["s1", "s2", "s3"]
    # s3's missing (N) call becomes 'N'; width stays one char per sample.
    assert rows == [("2", "CGN")]
    assert all(len(alleles) == len(names) for _pos, alleles in rows)


def test_maf_to_sites_argweaver_sites_add_ref_appends_ref_base(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "AGGT")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(out_prefix),
            "--samples", "s1", "s2",
            "--max-missing-count", "0",
            "--add-ref",
            "--emit-argweaver-sites",
        ],
        cwd=Path.cwd(),
    )

    names, _region, rows = _read_sites(Path(str(out_prefix) + ".sites"))
    assert names == ["s1", "s2", "REF"]
    # REF haplotype is appended last and carries the reference base (C at pos 2).
    assert rows == [("2", "CGC")]


def test_maf_to_sites_accepts_pairwise_blocks_with_same_src_names(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nA\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t1\t6\t1\t2\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_same_src_pairwise_maf(maf_dir / "s1.maf", "chr1", "A", "A")
    _write_same_src_pairwise_maf(maf_dir / "s2.maf", "chr1", "A", "C")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
        ],
        cwd=Path.cwd(),
    )

    all_records = _read_vcf_records(Path(str(out_prefix) + ".all_sites.vcf"))
    variant_records = _read_vcf_records(Path(str(out_prefix) + ".vcf"))
    summary = Path(str(out_prefix) + ".site_summary.tsv").read_text(encoding="utf-8")

    assert len(all_records) == 1
    assert len(variant_records) == 1
    assert all_records[0][3] == "A"
    assert all_records[0][4] == "C"
    assert "all_sites\t1" in summary
    assert "masked_no_alignment\t0" in summary


def test_maf_to_sites_masks_snp_adjacent_to_insertion(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "AC-GT", "s1", "ATCGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "AC-GT", "s2", "AC-GT")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
            "--mask-indel-adjacent-snps",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".vcf")
    masked = Path(str(out_prefix) + ".mask.bed")

    all_records = _read_vcf_records(all_sites)
    assert [record[1] for record in all_records] == ["1", "3", "4"]
    assert [record[1] for record in _read_vcf_records(variants)] == []
    assert _read_bed(masked) == [("chr1", 1, 2)]


def test_maf_to_sites_masks_snp_after_leading_insertion(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nAC\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t2\t6\t2\t3\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "-AC", "s1", "TGC")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "-AC", "s2", "-AC")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
            "--mask-indel-adjacent-snps",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".vcf")
    masked = Path(str(out_prefix) + ".mask.bed")

    assert [record[1] for record in _read_vcf_records(all_sites)] == ["2"]
    assert _read_vcf_records(variants) == []
    assert _read_bed(masked) == [("chr1", 0, 1)]


def test_maf_to_sites_can_keep_snp_adjacent_to_insertion(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "AC-GT", "s1", "ATCGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "AC-GT", "s2", "AC-GT")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
            "--keep-indel-adjacent-snps",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".vcf")
    masked = Path(str(out_prefix) + ".mask.bed")

    assert [record[1] for record in _read_vcf_records(all_sites)] == ["1", "2", "3", "4"]
    assert [record[1] for record in _read_vcf_records(variants)] == ["2"]
    assert _read_bed(masked) == []


def test_maf_to_sites_can_keep_snp_adjacent_to_indel(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t8\t6\t8\t9\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGTACGT", "s1", "ACGTACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGTACGT", "s2", "AT-TCCGN")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
            "--keep-indel-adjacent-snps",
        ],
        cwd=Path.cwd(),
    )

    all_sites = Path(str(out_prefix) + ".all_sites.vcf")
    variants = Path(str(out_prefix) + ".vcf")
    masked = Path(str(out_prefix) + ".mask.bed")

    assert [record[1] for record in _read_vcf_records(all_sites)] == ["1", "2", "4", "5", "6", "7"]
    assert [record[1] for record in _read_vcf_records(variants)] == ["2", "5"]
    assert _read_bed(masked) == [("chr1", 2, 3), ("chr1", 7, 8)]


def test_maf_to_sites_counts_all_deleted_sites_as_missingness(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nAAAA\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "AAAA", "s1", "A---")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "AAAA", "s2", "A---")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
        ],
        cwd=Path.cwd(),
    )

    summary = Path(str(out_prefix) + ".site_summary.tsv").read_text(encoding="utf-8")
    assert "masked_missingness\t3\n" in summary
    assert "masked_no_alignment\t0\n" in summary
    assert "masked_indel_adjacent\t0\n" in summary


def test_maf_to_sites_counts_missing_alignment_blocks_as_no_alignment(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "A", "s1", "A")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "ACGT")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
        ],
        cwd=Path.cwd(),
    )

    summary = Path(str(out_prefix) + ".site_summary.tsv").read_text(encoding="utf-8")
    assert "masked_no_alignment\t3\n" in summary
    assert "masked_missingness\t0\n" in summary


def test_maf_to_sites_matches_normalized_contig_names(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("1\t4\t3\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "ATGT")

    out_prefix = tmp_path / "results" / "combined.1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
        ],
        cwd=Path.cwd(),
    )

    all_records = _read_vcf_records(Path(str(out_prefix) + ".all_sites.vcf"))
    variant_records = _read_vcf_records(Path(str(out_prefix) + ".vcf"))
    summary = Path(str(out_prefix) + ".site_summary.tsv").read_text(encoding="utf-8")

    assert [record[1] for record in all_records] == ["1", "2", "3", "4"]
    assert [record[1] for record in variant_records] == ["2"]
    assert "masked_no_alignment\t0\n" in summary


def test_maf_to_sites_ignores_maf_columns_past_reference_contig_length(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nAC\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t2\t6\t2\t3\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "ATGT")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir",
            str(maf_dir),
            "--reference-fasta",
            str(ref),
            "--contig",
            "chr1",
            "--out-prefix",
            str(out_prefix),
            "--samples",
            "s1",
            "s2",
            "--max-missing-count",
            "0",
        ],
        cwd=Path.cwd(),
    )

    all_records = _read_vcf_records(Path(str(out_prefix) + ".all_sites.vcf"))
    variant_records = _read_vcf_records(Path(str(out_prefix) + ".vcf"))

    assert [record[1] for record in all_records] == ["1", "2"]
    assert [record[1] for record in variant_records] == ["2"]


def test_discover_samples_handles_maf_and_maf_gz(tmp_path: Path) -> None:
    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    (maf_dir / "sample_b.maf").write_text("", encoding="utf-8")
    (maf_dir / "sample_a.maf.gz").write_text("", encoding="utf-8")
    (maf_dir / "ignore.txt").write_text("", encoding="utf-8")

    assert discover_samples(maf_dir) == ["sample_a", "sample_b"]


def test_maf_path_for_sample_prefers_plain_then_gz(tmp_path: Path) -> None:
    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    plain = maf_dir / "sample1.maf"
    gz = maf_dir / "sample1.maf.gz"
    plain.write_text("", encoding="utf-8")
    gz.write_text("", encoding="utf-8")

    assert maf_path_for_sample(maf_dir, "sample1") == plain

    plain.unlink()
    assert maf_path_for_sample(maf_dir, "sample1") == gz


def test_maf_path_for_sample_raises_for_missing_sample(tmp_path: Path) -> None:
    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()

    with pytest.raises(FileNotFoundError, match="Missing MAF for sample 'sample1'"):
        maf_path_for_sample(maf_dir, "sample1")


def test_parse_maf_path_map_overrides_sample_paths(tmp_path: Path) -> None:
    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    fallback = maf_dir / "sample1.maf"
    fallback.write_text("", encoding="utf-8")
    sample2 = maf_dir / "sample2.maf"
    sample2.write_text("", encoding="utf-8")
    chunk = tmp_path / "chunks" / "sample1" / "chr1.maf"
    chunk.parent.mkdir(parents=True)
    chunk.write_text("", encoding="utf-8")

    paths = parse_maf_path_map([f"sample1={chunk}"])

    assert maf_path_for_sample_with_map(maf_dir, "sample1", paths) == chunk
    assert maf_path_for_sample_with_map(maf_dir, "sample2", paths) == sample2


@pytest.mark.parametrize(
    ("value", "message"),
    [
        ("sample1", "expected SAMPLE=PATH"),
        ("=/tmp/sample1.maf", "sample is empty"),
    ],
)
def test_parse_maf_path_map_rejects_invalid_entries(value: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        parse_maf_path_map([value])


def test_split_maf_by_contig_writes_only_requested_chunks(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    maf = tmp_path / "sample1.maf"
    maf.write_text(
        "##maf version=1\n"
        "a score=0\n"
        "s chr1 0 2 + 2 AC\n"
        "s sample1 0 2 + 2 AT\n"
        "\n"
        "a score=0\n"
        "s chr2 0 2 + 2 GG\n"
        "s sample1 0 2 + 2 GA\n"
        "\n"
        "a score=0\n"
        "s chr3 0 2 + 2 TT\n"
        "s sample1 0 2 + 2 TC\n",
        encoding="utf-8",
    )
    out_root = tmp_path / "chunks"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "split_maf_by_contig.py",
            "--maf", str(maf),
            "--sample", "sample1",
            "--out-root", str(out_root),
            "--contigs", "1", "chr2",
        ],
    )
    split_maf_main()

    chr1 = out_root / "sample1" / "1.maf.gz"
    chr2 = out_root / "sample1" / "chr2.maf.gz"
    with gzip.open(chr1, "rt", encoding="utf-8") as handle:
        chr1_text = handle.read()
    with gzip.open(chr2, "rt", encoding="utf-8") as handle:
        chr2_text = handle.read()
    assert "s chr1 0 2 + 2 AC" in chr1_text
    assert "s chr2 0 2 + 2 GG" in chr2_text
    assert "chr3" not in chr1_text
    assert "chr3" not in chr2_text


def test_split_maf_by_contig_reads_gzip_input(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    maf = tmp_path / "sample1.maf.gz"
    with gzip.open(maf, "wt", encoding="utf-8") as handle:
        handle.write(
            "##maf version=1\n"
            "a score=0\n"
            "s chr1 0 2 + 2 AC\n"
            "s sample1 0 2 + 2 AT\n",
        )
    out_root = tmp_path / "chunks"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "split_maf_by_contig.py",
            "--maf", str(maf),
            "--sample", "sample1",
            "--out-root", str(out_root),
            "--contigs", "chr1",
        ],
    )
    split_maf_main()

    chunk = out_root / "sample1" / "chr1.maf.gz"
    with gzip.open(chunk, "rt", encoding="utf-8") as handle:
        assert "s chr1 0 2 + 2 AC" in handle.read()


def test_split_maf_by_contig_rejects_ambiguous_normalized_contigs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    maf = tmp_path / "sample1.maf"
    maf.write_text(
        "##maf version=1\n"
        "a score=0\n"
        "s chr1 0 2 + 2 AC\n"
        "s sample1 0 2 + 2 AT\n",
        encoding="utf-8",
    )
    out_root = tmp_path / "chunks"

    # "chr1" and "1" normalize to the same key but map to different outputs.
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "split_maf_by_contig.py",
            "--maf", str(maf),
            "--sample", "sample1",
            "--out-root", str(out_root),
            "--contigs", "chr1", "1",
        ],
    )
    with pytest.raises(ValueError, match="Ambiguous normalized contig"):
        split_maf_main()


def test_read_sample_missing_bp_matches_longest_sample_suffix(tmp_path: Path) -> None:
    beds = []
    for name, body in [
        ("combined.chr1.s1.missing.bed", "chr1\t0\t2\ts1\n"),
        ("combined.chr1.xs1.missing.bed", "chr1\t2\t5\txs1\n"),
    ]:
        path = tmp_path / name
        path.write_text(body, encoding="utf-8")
        beds.append(path)

    data = read_sample_missing_bp(beds, ["s1", "xs1"], {"chr1": 10})

    assert data == {"s1": {"chr1": 2}, "xs1": {"chr1": 3}}


def test_read_contig_sequence_raises_when_contig_is_missing(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Contig 'chr2' not found"):
        read_contig_sequence(ref, "chr2")


@pytest.mark.parametrize("fraction", [-0.1, 1.1])
def test_missing_threshold_rejects_invalid_fraction(fraction: float) -> None:
    with pytest.raises(ValueError, match="--max-missing-fraction must be between 0 and 1"):
        missing_threshold(4, None, fraction)


def test_load_quality_mask_keeps_only_below_threshold(tmp_path: Path) -> None:
    bed = tmp_path / "s.bed"
    bed.write_text(
        "# comment\n"
        "track name=quality\n"
        "chrX\t10\t20\t0.95\n"  # passes threshold -> ignored
        "chrX\t30\t40\t0.5\n"  # below threshold -> low quality
        "chrY\t0\t5\t0.1\n",
        encoding="utf-8",
    )
    mask = load_quality_mask(bed, 0.9)

    assert mask.is_low("chrX", 35) is True
    assert mask.is_low("chrX", 15) is False  # passing interval not stored
    assert mask.is_low("chrX", 40) is False  # half-open end excluded
    assert mask.is_low("chrY", 3) is True
    assert mask.is_low("chrZ", 0) is False


@pytest.mark.parametrize(
    ("row", "message"),
    [
        ("chrX\tbad\trow\n", "expected at least 4 fields"),
        ("chrX\tbad\t20\t0.5\n", "start and end must be integers"),
        ("chrX\t10\t20\tbad\n", "score must be numeric"),
    ],
)
def test_load_quality_mask_rejects_malformed_data_rows(
    tmp_path: Path, row: str, message: str
) -> None:
    bed = tmp_path / "s.bed"
    bed.write_text("# comment\n" + row, encoding="utf-8")

    with pytest.raises(ValueError, match=message) as exc_info:
        load_quality_mask(bed, 0.9)

    assert f"{bed} at line 2" in str(exc_info.value)


def test_iter_maf_blocks_rejects_unequal_alignment_strings(tmp_path: Path) -> None:
    maf = tmp_path / "bad.maf"
    maf.write_text(
        "a score=0\n"
        "s chr1 0 2 + 2 AC\n"
        "s sample 0 1 + 1 A\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="alignment strings have unequal lengths"):
        list(iter_maf_blocks(maf))


def test_iter_maf_blocks_reports_malformed_numeric_field_context(tmp_path: Path) -> None:
    maf = tmp_path / "bad.maf"
    maf.write_text("a score=0\ns chr1 nope 2 + 2 AC\n", encoding="utf-8")

    with pytest.raises(ValueError, match=r"bad\.maf at line 2.*must be integers"):
        list(iter_maf_blocks(maf))


def test_iter_maf_blocks_rejects_short_sequence_row(tmp_path: Path) -> None:
    maf = tmp_path / "bad.maf"
    maf.write_text("a score=0\ns chr1 0 2 + 2\n", encoding="utf-8")

    with pytest.raises(ValueError, match=r"bad\.maf at line 2.*expected 7 fields"):
        list(iter_maf_blocks(maf))


def test_quality_mask_treats_plus_strand_bases_as_missing(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    # s2 carries a SNP (G) at ref idx 1; without masking this is a variant site.
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "AGGT")

    quality_dir = tmp_path / "quality"
    quality_dir.mkdir()
    # s2's own forward coordinate of the second aligned base is 1 (start 0, + strand).
    (quality_dir / "s2.bed").write_text("s2\t1\t2\t0.5\n", encoding="utf-8")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(out_prefix),
            "--samples", "s1", "s2",
            "--max-missing-count", "1",
            "--quality-bed-dir", str(quality_dir),
            "--quality-min", "0.9",
        ],
        cwd=Path.cwd(),
    )

    all_records = _read_vcf_records(Path(str(out_prefix) + ".all_sites.vcf"))
    variant_records = _read_vcf_records(Path(str(out_prefix) + ".vcf"))
    s2_mask = _read_bed(Path(str(out_prefix) + ".s2.missing.bed"))

    # The masked base makes ref idx 1 invariant, so no variant is emitted there.
    assert [record[1] for record in variant_records] == []
    pos2 = next(record for record in all_records if record[1] == "2")
    assert pos2[9:] == ["0", "."]  # s2 now missing at the masked position
    assert s2_mask == [("chr1", 1, 2, "s2")]


def test_quality_mask_handles_minus_strand_coordinates(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    # s2 aligns on the minus strand (start 0, src_size 4). The second aligned
    # base (ref idx 1) has forward coordinate 4 - 0 - 1 - 1 = 2.
    _write_pairwise_maf_strand(
        maf_dir / "s2.maf", "chr1", "ACGT", "s2", "AGGT", "-", 0, 4
    )

    quality_dir = tmp_path / "quality"
    quality_dir.mkdir()
    (quality_dir / "s2.bed").write_text("s2\t2\t3\t0.5\n", encoding="utf-8")

    out_prefix = tmp_path / "results" / "combined.chr1"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(out_prefix),
            "--samples", "s1", "s2",
            "--max-missing-count", "1",
            "--quality-bed-dir", str(quality_dir),
            "--quality-min", "0.9",
        ],
        cwd=Path.cwd(),
    )

    all_records = _read_vcf_records(Path(str(out_prefix) + ".all_sites.vcf"))
    variant_records = _read_vcf_records(Path(str(out_prefix) + ".vcf"))

    # Minus-strand coordinate 2 maps to ref idx 1, masking the SNP there.
    assert [record[1] for record in variant_records] == []
    pos2 = next(record for record in all_records if record[1] == "2")
    assert pos2[9:] == ["0", "."]


def test_quality_mask_requires_quality_min(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nAC\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t2\t6\t2\t3\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "AC", "s1", "AC")

    quality_dir = tmp_path / "quality"
    quality_dir.mkdir()

    proc = subprocess.run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(tmp_path / "results" / "combined.chr1"),
            "--samples", "s1",
            "--quality-bed-dir", str(quality_dir),
        ],
        cwd=Path.cwd(),
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "--quality-bed-dir requires --quality-min" in proc.stderr


def test_quality_min_requires_quality_bed_dir(tmp_path: Path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nAC\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr1\t2\t6\t2\t3\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "AC", "s1", "AC")

    proc = subprocess.run(
        [
            sys.executable,
            str(Path("scripts") / "maf_to_sites.py"),
            "--maf-dir", str(maf_dir),
            "--reference-fasta", str(ref),
            "--contig", "chr1",
            "--out-prefix", str(tmp_path / "results" / "combined.chr1"),
            "--samples", "s1",
            "--quality-min", "0.9",
        ],
        cwd=Path.cwd(),
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "--quality-min requires --quality-bed-dir" in proc.stderr


@pytest.mark.parametrize("strand", ["-", "?"])
def test_choose_sample_record_rejects_non_plus_reference_row(tmp_path: Path, strand: str) -> None:
    maf = tmp_path / "s1.maf"
    maf.write_text(f"a\ns chr1 0 4 {strand} 10 ACGT\ns s1 0 4 + 10 ACGT\n\n")
    block = next(iter_maf_blocks(maf))
    with pytest.raises(ValueError, match="only '\\+' reference rows are supported"):
        choose_sample_record(block, "chr1")
    with pytest.raises(ValueError, match="reference rows are supported"):
        load_sample_calls(maf, "chr1", 10)


def test_choose_sample_record_ignores_strand_of_other_contigs(tmp_path: Path) -> None:
    maf = tmp_path / "s1.maf"
    maf.write_text("a\ns chr2 0 4 - 10 ACGT\ns s1 0 4 + 10 ACGT\n\n")
    assert choose_sample_record(next(iter_maf_blocks(maf)), "chr1") is None
