import os
import subprocess
from pathlib import Path
import shutil

import pytest


def _conda_exe() -> str:
    return os.environ.get("CONDA_EXE") or "conda"


def _run(cmd, cwd):
    subprocess.run(cmd, cwd=cwd, check=True)


def _config_reference_fasta(config_path: Path) -> Path | None:
    if not config_path.exists():
        return None
    for raw in config_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("reference_fasta:"):
            value = line.split(":", 1)[1].strip().strip('"').strip("'")
            if value:
                return (config_path.parent / value).resolve()
    return None


def _require_integration():
    if os.getenv("RUN_INTEGRATION") != "1":
        pytest.skip("Set RUN_INTEGRATION=1 to run integration tests.")


def _require_conda_tools(env: str, *tools: str):
    missing = []
    for tool in tools:
        proc = subprocess.run(
            [_conda_exe(), "run", "-n", env, tool, "--version"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if proc.returncode != 0:
            missing.append(tool)
    if missing:
        pytest.skip(f"Missing required tools in conda env {env}: {', '.join(missing)}")


def _is_gzip(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def test_snakemake_summary(tmp_path: Path):
    _require_integration()
    _require_conda_tools("argprep", "snakemake", "bcftools", "tabix", "bgzip", "gatk", "java")
    repo = Path.cwd()
    ref = _config_reference_fasta(repo / "config.yaml")
    if ref is None or not ref.exists():
        pytest.skip(f"Missing reference FASTA from config.yaml: {ref}")
    maf_files = list((repo / "example_data").glob("*.maf*"))
    if not maf_files:
        pytest.skip("No example MAFs found under example_data/")
    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    for maf in maf_files:
        gzipped = _is_gzip(maf)
        if gzipped and maf.suffix != ".gz":
            dest_name = maf.name + ".gz"
        elif (not gzipped) and maf.name.endswith(".gz"):
            dest_name = maf.name[: -len(".gz")]
        else:
            dest_name = maf.name
        shutil.copyfile(maf, maf_dir / dest_name)
    summary_target = str(Path.cwd() / "results" / "summary.html")
    _run(
        [
            _conda_exe(),
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
            "--config",
            f"maf_dir={maf_dir}",
            "--",
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
