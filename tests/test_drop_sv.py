import os
import subprocess
import shutil
from pathlib import Path

import pytest


def _conda_exe() -> str:
    return os.environ.get("CONDA_EXE") or "conda"


def _conda_cmd(*args: str) -> list[str]:
    return [_conda_exe(), "run", "-n", "argprep", *args]


def _require_conda_tools(*tools: str):
    if shutil.which(_conda_exe()) is None:
        pytest.skip("conda not available")
    missing = []
    for tool in tools:
        proc = subprocess.run(
            _conda_cmd(tool, "--version"),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if proc.returncode != 0:
            missing.append(tool)
    if missing:
        pytest.skip(f"Missing required tools in conda env argprep: {', '.join(missing)}")


def _run(cmd, cwd=None):
    subprocess.run(cmd, cwd=cwd, check=True)


def _write_vcf(path: Path):
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t2\t.\tA\tT\t.\t.\t.\n"
        "1\t5\t.\tAAAA\tA\t.\t.\t.\n",
        encoding="utf-8",
    )


def _count_records(path: Path) -> int:
    count = 0
    with subprocess.Popen(_conda_cmd("bcftools", "view", "-H", str(path)), stdout=subprocess.PIPE, text=True) as proc:
        assert proc.stdout is not None
        for _ in proc.stdout:
            line = _.rstrip("\n")
            if not line:
                continue
            # Only count true VCF records (tab-delimited with standard columns).
            if len(line.split("\t")) >= 8:
                count += 1
    return count


def test_drop_sv_filters_large_indel(tmp_path: Path):
    _require_conda_tools("bcftools", "bgzip", "tabix")
    vcf = tmp_path / "sample.gvcf"
    _write_vcf(vcf)

    _run(_conda_cmd("bgzip", "-f", str(vcf)))
    gvcf = tmp_path / "sample.gvcf.gz"
    _run(_conda_cmd("tabix", "-p", "vcf", str(gvcf)))

    _run(_conda_cmd(
        "python3",
        str(Path("scripts") / "dropSV.py"),
        "-d",
        str(tmp_path),
        "-c",
        "1",
    ))

    cleaned = tmp_path / "cleangVCF" / "sample.gvcf.gz"
    assert cleaned.exists()
    # Only SNP should remain.
    assert _count_records(cleaned) == 1

    dropped_bed = tmp_path / "cleangVCF" / "dropped_indels.bed"
    assert dropped_bed.exists()
    bed_lines = [l for l in dropped_bed.read_text(encoding="utf-8").splitlines() if l.strip()]
    assert len(bed_lines) == 1
