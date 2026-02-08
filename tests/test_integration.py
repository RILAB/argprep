import os
import subprocess
from pathlib import Path

import pytest


def _run(cmd, cwd):
    subprocess.run(cmd, cwd=cwd, check=True)


def _require_integration():
    if os.getenv("RUN_INTEGRATION") != "1":
        pytest.skip("Set RUN_INTEGRATION=1 to run integration tests.")


def _require_conda_tools(env: str, *tools: str):
    missing = []
    for tool in tools:
        proc = subprocess.run(
            ["conda", "run", "-n", env, tool, "--version"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if proc.returncode != 0:
            missing.append(tool)
    if missing:
        pytest.skip(f"Missing required tools in conda env {env}: {', '.join(missing)}")


def test_snakemake_summary():
    _require_integration()
    _require_conda_tools("argprep", "snakemake", "bcftools", "tabix", "bgzip", "gatk", "java")
    summary_target = str(Path.cwd() / "results" / "summary.html")
    _run(
        [
            "conda",
            "run",
            "-n",
            "argprep",
            "env",
            "PATH=/opt/anaconda3/envs/argprep/bin:/usr/bin:/bin:/usr/sbin:/sbin",
            "HOME=/tmp",
            "TMPDIR=/tmp",
            "XDG_CACHE_HOME=/tmp",
            "SNAKEMAKE_OUTPUT_CACHE=/tmp/snakemake",
            "SNAKEMAKE_SOURCE_CACHE=/tmp/snakemake",
            "snakemake",
            "-j",
            "2",
            summary_target,
        ],
        cwd=Path.cwd(),
    )


def test_maf_to_gvcf_single_sample(tmp_path: Path):
    _require_integration()
    _require_conda_tools("argprep", "java")
    # Minimal two-sequence MAF (ref + sample) for a small region.
    ref = "ACGTACGTACGTACGT"
    sample = "ACGTACGTACGTACGA"
    maf = tmp_path / "mini.maf"
    maf.write_text(
        "##maf version=1\n"
        "# synthetic alignment\n\n"
        "a score=0\n"
        f"s 1 0 {len(ref)} + {len(ref)} {ref}\n"
        f"s sample 0 {len(sample)} + {len(sample)} {sample}\n",
        encoding="utf-8",
    )
    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text(">1\n" + ref + "\n", encoding="utf-8")
    out = tmp_path / "out.gvcf"

    _run(
        [
            str(Path("tassel-5-standalone") / "run_pipeline.pl"),
            "-Xmx2G",
            "-debug",
            "-MAFToGVCFPlugin",
            "-referenceFasta",
            str(ref_fa),
            "-mafFile",
            str(maf),
            "-sampleName",
            "sample_anchorwave",
            "-gvcfOutput",
            str(out),
            "-fillGaps",
            "false",
        ],
        cwd=Path.cwd(),
    )
    assert out.exists() or (out.with_suffix(out.suffix + ".gz")).exists()
