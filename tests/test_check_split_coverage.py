import subprocess
import sys
from pathlib import Path


def _run(cmd, cwd):
    subprocess.run(cmd, cwd=cwd, check=True)


def test_check_split_coverage_pass(tmp_path: Path):
    prefix = tmp_path / "out"
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("1\t10\t0\t0\t0\n", encoding="utf-8")

    (tmp_path / "out.clean").write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1\t.\tA\tC\t.\t.\tDP=1\n"
        "1\t2\t.\tA\tG\t.\t.\tDP=1\n",
        encoding="utf-8",
    )
    (tmp_path / "out.inv").write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t3\t.\tA\t.\t.\t.\tDP=1\n"
        "1\t4\t.\tA\t.\t.\t.\tDP=1\n",
        encoding="utf-8",
    )
    (tmp_path / "out.filtered.bed").write_text(
        "1\t4\t10\n",
        encoding="utf-8",
    )

    _run(
        [
            sys.executable,
            str(Path("scripts") / "check_split_coverage.py"),
            str(prefix),
            "--fai",
            str(fai),
        ],
        cwd=Path.cwd(),
    )

    report = tmp_path / "out.coverage.txt"
    assert report.exists()


def test_check_split_coverage_overlap_fails(tmp_path: Path):
    prefix = tmp_path / "out"
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("1\t5\t0\t0\t0\n", encoding="utf-8")

    (tmp_path / "out.clean").write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1\t.\tA\tC\t.\t.\tDP=1\n",
        encoding="utf-8",
    )
    (tmp_path / "out.inv").write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t2\t.\tA\t.\t.\t.\tDP=1\n",
        encoding="utf-8",
    )
    # Overlaps clean at pos 1.
    (tmp_path / "out.filtered.bed").write_text(
        "1\t0\t2\n",
        encoding="utf-8",
    )

    try:
        _run(
            [
                sys.executable,
                str(Path("scripts") / "check_split_coverage.py"),
                str(prefix),
                "--fai",
                str(fai),
            ],
            cwd=Path.cwd(),
        )
    except subprocess.CalledProcessError:
        return
    raise AssertionError("check_split_coverage should have failed on overlap")
