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


def test_filt_to_bed_uses_ref_len(tmp_path: Path):
    prefix = tmp_path / "sample"

    filtered = tmp_path / "sample.filtered"
    filtered.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t10\t.\tAAA\tT\t.\t.\t.\n",
        encoding="utf-8",
    )

    missing = tmp_path / "sample.missing.bed"
    missing.write_text("1\t0\t1\n", encoding="utf-8")

    dropped = tmp_path / "dropped_indels.bed"
    dropped.write_text("", encoding="utf-8")

    (tmp_path / "sample.inv").write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", encoding="utf-8")
    (tmp_path / "sample.clean").write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", encoding="utf-8")

    _run(
        [
            sys.executable,
            str(Path("scripts") / "filt_to_bed.py"),
            str(prefix),
            "--dropped-bed",
            str(dropped),
        ],
        cwd=Path.cwd(),
    )

    out_bed = tmp_path / "sample.filtered.bed"
    assert out_bed.exists()

    # REF length 3 => interval length 3, plus missing interval length 1
    assert _sum_bed(out_bed) == 4
