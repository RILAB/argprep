import subprocess
import sys
from pathlib import Path


def _run(cmd, cwd):
    subprocess.run(cmd, cwd=cwd, check=True)


def _sum_bed(path: Path) -> int:
    total = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end = line.rstrip("\n").split("\t")[:3]
            total += int(end) - int(start)
    return total


def test_split_missing_with_fai(tmp_path: Path):
    ref_fai = tmp_path / "ref.fa.fai"
    ref_fai.write_text("1\t10\t0\t0\t0\n", encoding="utf-8")

    gvcf = tmp_path / "in.gvcf"
    gvcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=10>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "1\t3\t.\tA\t<NON_REF>\t.\t.\tEND=5\tGT:AD\t0/0:1,0\n"
        "1\t8\t.\tC\tT\t.\t.\tDP=5\tGT:AD\t0/1:1,1\n",
        encoding="utf-8",
    )

    prefix = tmp_path / "out"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "split.py"),
            "--out-prefix",
            str(prefix),
            "--fai",
            str(ref_fai),
            str(gvcf),
        ],
        cwd=Path.cwd(),
    )

    missing = tmp_path / "out.missing.bed"
    assert missing.exists()
    # Expected missing length: contig length (10) - covered (3-5 and 8) = 10 - 4 = 6
    assert _sum_bed(missing) == 6


def test_singer_dp_rewrite(tmp_path: Path):
    ref_fai = tmp_path / "ref.fa.fai"
    ref_fai.write_text("1\t5\t0\t0\t0\n", encoding="utf-8")

    gvcf = tmp_path / "in.gvcf"
    gvcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
        "1\t2\t.\tA\tG\t.\t.\tDP=30\tGT:AD\t0/1:0,1\t./.:.\t0/0:1,0\n",
        encoding="utf-8",
    )

    prefix = tmp_path / "out"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "split.py"),
            "--out-prefix",
            str(prefix),
            "--fai",
            str(ref_fai),
            str(gvcf),
        ],
        cwd=Path.cwd(),
    )

    clean = tmp_path / "out.clean"
    line = [l for l in clean.read_text(encoding="utf-8").splitlines() if not l.startswith("#")][0]
    info = line.split("\t")[7]
    # Two non-missing samples -> DP=2
    assert "DP=2" in info


def _count_vcf_records(path: Path) -> int:
    return len([l for l in path.read_text(encoding="utf-8").splitlines() if l and not l.startswith("#")])


def test_split_invariant_reference_block(tmp_path: Path):
    ref_fai = tmp_path / "ref.fa.fai"
    ref_fai.write_text("1\t5\t0\t0\t0\n", encoding="utf-8")

    gvcf = tmp_path / "in.gvcf"
    gvcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=5>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=5\tGT:AD\t0/0:5,0\n",
        encoding="utf-8",
    )

    prefix = tmp_path / "out"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "split.py"),
            "--out-prefix",
            str(prefix),
            "--fai",
            str(ref_fai),
            str(gvcf),
        ],
        cwd=Path.cwd(),
    )

    inv = tmp_path / "out.inv"
    missing = tmp_path / "out.missing.bed"

    assert inv.exists()
    assert missing.exists()
    assert _count_vcf_records(inv) == 5
    assert _sum_bed(missing) == 0


def test_split_alt_dot_requires_called_gt_for_inv(tmp_path: Path):
    ref_fai = tmp_path / "ref.fa.fai"
    ref_fai.write_text("1\t3\t0\t0\t0\n", encoding="utf-8")

    gvcf = tmp_path / "in.gvcf"
    gvcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "1\t1\t.\tA\t.\t.\t.\tDP=10\tGT\t.\n"
        "1\t2\t.\tC\t.\t.\t.\tDP=10\tGT\t0\n",
        encoding="utf-8",
    )

    prefix = tmp_path / "out"
    _run(
        [
            sys.executable,
            str(Path("scripts") / "split.py"),
            "--out-prefix",
            str(prefix),
            "--fai",
            str(ref_fai),
            str(gvcf),
        ],
        cwd=Path.cwd(),
    )

    inv = tmp_path / "out.inv"
    clean = tmp_path / "out.clean"

    inv_records = _count_vcf_records(inv)
    clean_records = _count_vcf_records(clean)

    # ALT="." records are invariant regardless of GT in current logic.
    assert inv_records == 2
    assert clean_records == 0
