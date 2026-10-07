import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest


def _write_pairwise_maf(path: Path, contig: str, ref_seq: str, sample: str, sample_seq: str) -> None:
    path.write_text(
        "##maf version=1\n"
        "a score=0\n"
        f"s {contig} 0 {len(ref_seq.replace('-', ''))} + {len(ref_seq.replace('-', ''))} {ref_seq}\n"
        f"s {sample} 0 {len(sample_seq.replace('-', ''))} + {len(sample_seq.replace('-', ''))} {sample_seq}\n",
        encoding="utf-8",
    )


def _run_snakemake(tmp_path: Path, config: Path | None, *targets: str) -> subprocess.CompletedProcess[str]:
    repo_root = Path(__file__).resolve().parents[1]
    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(repo_root / "Snakefile"),
        "--directory",
        str(tmp_path),
        "--cores",
        "1",
    ]
    if config is not None:
        cmd.extend(["--configfile", str(config)])
    cmd.extend(["--", *targets])
    # Ensure tools (samtools, etc.) from the active conda env are on PATH
    env = os.environ.copy()
    env_bin = str(Path(sys.executable).parent)
    env["PATH"] = env_bin + os.pathsep + env.get("PATH", "")
    env["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    return subprocess.run(
        cmd,
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_remaps_requested_contigs_to_reference_names(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "chr1", "ACGT", "s2", "ATGT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'contigs: ["chr01"]',
                'samples: ["s1", "s2"]',
                "max_missing_count: 0",
                "mask_indel_adjacent_snps: false",
                "allow_multiallelic_snps: true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode == 0, result.stderr

    all_sites = tmp_path / "results" / "sites" / "combined.1.all_sites.vcf"
    assert all_sites.exists()
    assert (tmp_path / "results" / "maf_by_contig" / "s1" / "1.maf.gz").exists()
    assert (tmp_path / "results" / "maf_by_contig" / "s2" / "1.maf.gz").exists()
    assert not (tmp_path / "results" / "sites" / "combined.chr01.all_sites.vcf").exists()

    records = [
        line.rstrip("\n").split("\t")
        for line in all_sites.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert [record[0] for record in records] == ["1", "1", "1", "1"]
    assert [record[1] for record in records] == ["1", "2", "3", "4"]


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_default_contig_discovery_normalizes_maf_aliases(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "ACGT", "s2", "ATGT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["s1", "s2"]',
                "max_missing_count: 0",
                "mask_indel_adjacent_snps: false",
                "allow_multiallelic_snps: true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode == 0, result.stderr

    all_sites = tmp_path / "results" / "sites" / "combined.1.all_sites.vcf"
    assert all_sites.exists()
    records = [
        line.rstrip("\n").split("\t")
        for line in all_sites.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert [record[0] for record in records] == ["1", "1", "1", "1"]
    assert [record[1] for record in records] == ["1", "2", "3", "4"]


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_rejects_ambiguous_contig_remap(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n>chr1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "chr1", "ACGT", "s1", "ACGT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'contigs: ["chr01"]',
                'samples: ["s1"]',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode != 0
    assert "unmatched or ambiguous: chr01" in (
        result.stderr + result.stdout
    )


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
@pytest.mark.parametrize(
    ("fraction", "expect_retained"),
    # 4 samples, one missing at the site: int(4 * 0.25) == 1 allows it through,
    # int(4 * 0.1) == 0 does not. Guards the config -> CLI path for a threshold
    # that decides which sites are retained (cf. the allow_multiallelic_snps bug,
    # where a config value never reached the caller and the suite stayed green).
    [(0.25, True), (0.1, False)],
)
def test_workflow_propagates_max_missing_fraction(
    tmp_path: Path, fraction: float, expect_retained: bool
) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nA\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    for sample in ("s1", "s2", "s3"):
        _write_pairwise_maf(maf_dir / f"{sample}.maf", "1", "A", sample, "A")
    # s4 aligns but carries N, so it is missing at the only site.
    _write_pairwise_maf(maf_dir / "s4.maf", "1", "A", "s4", "N")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["s1", "s2", "s3", "s4"]',
                f"max_missing_fraction: {fraction}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    all_sites = tmp_path / "results" / "sites" / "combined.1.all_sites.vcf"
    result = _run_snakemake(tmp_path, config, str(all_sites))
    assert result.returncode == 0, result.stderr

    records = [
        line for line in all_sites.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert bool(records) is expect_retained


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_combines_per_sample_missing_masks_by_contig(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nAAAA\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "AAAA", "s1", "NAAA")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "AAAA", "s2", "AANA")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["s1", "s2"]',
                "max_missing_count: 2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    combined = tmp_path / "results" / "sites" / "combined.1.sample.mask.bed"
    result = _run_snakemake(tmp_path, config)
    assert result.returncode == 0, result.stderr
    assert combined.read_text(encoding="utf-8") == (
        "1\t0\t1\ts1\n"
        "1\t2\t3\ts2\n"
    )


def _write_config(path: Path, *lines: str) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _vcf_records(path: Path) -> list[str]:
    return [
        line for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]


def _vcf_sample_columns(path: Path) -> list[str]:
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#CHROM"):
            return line.split("\t")[9:]
    raise AssertionError(f"no #CHROM header in {path}")


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
@pytest.mark.parametrize(
    ("mask_adjacent", "expected_variants"),
    # Sample carries a deletion at ref position 2 and a SNP at position 3, so
    # the SNP is indel-adjacent: masked when the flag is on, retained when off.
    [(True, 0), (False, 1)],
)
def test_workflow_propagates_mask_indel_adjacent_snps(
    tmp_path: Path, mask_adjacent: bool, expected_variants: int
) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nAAAA\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "AAAA", "s1", "A-CA")

    config = tmp_path / "config.yaml"
    _write_config(
        config,
        "maf_dir: maf",
        "reference_fasta: ref.fa",
        "results_dir: results",
        'samples: ["s1"]',
        f"mask_indel_adjacent_snps: {str(mask_adjacent).lower()}",
    )

    variants = tmp_path / "results" / "sites" / "combined.1.vcf"
    result = _run_snakemake(tmp_path, config, str(variants))
    assert result.returncode == 0, result.stderr
    assert len(_vcf_records(variants)) == expected_variants


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
@pytest.mark.parametrize(
    ("add_ref", "expected_columns"),
    [(True, ["s1", "REF"]), (False, ["s1"])],
)
def test_workflow_propagates_add_ref(
    tmp_path: Path, add_ref: bool, expected_columns: list[str]
) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")

    config = tmp_path / "config.yaml"
    _write_config(
        config,
        "maf_dir: maf",
        "reference_fasta: ref.fa",
        "results_dir: results",
        'samples: ["s1"]',
        f"add_ref: {str(add_ref).lower()}",
    )

    all_sites = tmp_path / "results" / "sites" / "combined.1.all_sites.vcf"
    result = _run_snakemake(tmp_path, config, str(all_sites))
    assert result.returncode == 0, result.stderr
    assert _vcf_sample_columns(all_sites) == expected_columns


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
@pytest.mark.parametrize(
    ("use_quality", "expected_sites"),
    # The BED scores s1's first base below quality_min, making it missing; with
    # the default allowed_missing of 0 that masks ref position 1 of 4.
    [(True, 3), (False, 4)],
)
def test_workflow_propagates_quality_mask(
    tmp_path: Path, use_quality: bool, expected_sites: int
) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")

    qual_dir = tmp_path / "qual"
    qual_dir.mkdir()
    # Coordinates are in the sample's own genome, keyed by its MAF query name.
    (qual_dir / "s1.bed").write_text("s1\t0\t1\t0.1\n", encoding="utf-8")

    lines = [
        "maf_dir: maf",
        "reference_fasta: ref.fa",
        "results_dir: results",
        'samples: ["s1"]',
    ]
    if use_quality:
        # quality_bed_dir is resolved by the Snakefile, so pass it absolute.
        lines += [f"quality_bed_dir: {qual_dir}", "quality_min: 0.5"]
    config = tmp_path / "config.yaml"
    _write_config(config, *lines)

    all_sites = tmp_path / "results" / "sites" / "combined.1.all_sites.vcf"
    result = _run_snakemake(tmp_path, config, str(all_sites))
    assert result.returncode == 0, result.stderr
    assert len(_vcf_records(all_sites)) == expected_sites


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_rejects_any_unmatched_explicit_contig(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'contigs: ["1", "missing"]',
                'samples: ["s1"]',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode != 0
    assert "unmatched or ambiguous: missing" in (result.stderr + result.stdout)


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_default_discovery_warns_and_skips_unmatched_contigs(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    maf_text = (
        "##maf version=1\n"
        "a score=0\n"
        "s chr1 0 4 + 4 ACGT\n"
        "s s1 0 4 + 4 ACGT\n\n"
        "a score=0\n"
        "s chr2 0 4 + 4 ACGT\n"
        "s s1 0 4 + 4 ACGT\n"
    )
    (maf_dir / "s1.maf").write_text(maf_text, encoding="utf-8")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["s1"]',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    all_sites = tmp_path / "results" / "sites" / "combined.1.all_sites.vcf"
    result = _run_snakemake(tmp_path, config, str(all_sites))
    assert result.returncode == 0, result.stderr
    assert "Skipped contigs (no unambiguous match in reference .fai): 2" in (
        result.stderr + result.stdout
    )
    assert all_sites.exists()


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
@pytest.mark.parametrize(
    ("allow_multiallelic", "expected_records", "expected_mask"),
    [(True, 1, ""), (False, 0, "1\t0\t1")],
)
def test_workflow_propagates_multiallelic_policy(
    tmp_path: Path,
    allow_multiallelic: bool,
    expected_records: int,
    expected_mask: str,
) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nA\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "A", "s1", "A")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "A", "s2", "C")
    _write_pairwise_maf(maf_dir / "s3.maf", "1", "A", "s3", "G")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["s1", "s2", "s3"]',
                "max_missing_count: 0",
                f"allow_multiallelic_snps: {str(allow_multiallelic).lower()}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    all_sites = tmp_path / "results" / "sites" / "combined.1.all_sites.vcf"
    result = _run_snakemake(tmp_path, config, str(all_sites))
    assert result.returncode == 0, result.stderr

    records = [
        line for line in all_sites.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert len(records) == expected_records
    mask = tmp_path / "results" / "sites" / "combined.1.mask.bed"
    assert mask.read_text(encoding="utf-8").strip() == expected_mask


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_add_ref_appends_reference_sample_to_final_vcfs(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nAC\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "AC", "s1", "AC")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "AC", "s2", "AT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["s1", "s2"]',
                "max_missing_count: 0",
                "mask_indel_adjacent_snps: false",
                "allow_multiallelic_snps: true",
                "add_ref: true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode == 0, result.stderr

    all_sites = tmp_path / "results" / "sites" / "combined.1.all_sites.vcf"
    variants = tmp_path / "results" / "sites" / "combined.1.vcf"

    all_lines = all_sites.read_text(encoding="utf-8").splitlines()
    variant_lines = variants.read_text(encoding="utf-8").splitlines()

    assert next(line for line in all_lines if line.startswith("#CHROM\t")).endswith("\ts1\ts2\tREF")
    assert next(line for line in variant_lines if line.startswith("#CHROM\t")).endswith("\ts1\ts2\tREF")

    all_records = [line.split("\t") for line in all_lines if line and not line.startswith("#")]
    variant_records = [line.split("\t") for line in variant_lines if line and not line.startswith("#")]

    assert [record[1] for record in all_records] == ["1", "2"]
    assert [record[1] for record in variant_records] == ["2"]
    assert all_records[0][9:] == ["0", "0", "0"]
    assert all_records[1][9:] == ["0", "1", "0"]
    assert variant_records[0][9:] == ["0", "1", "0"]


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_emit_argweaver_sites(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "ACGT", "s2", "AGGT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["s1", "s2"]',
                "max_missing_count: 0",
                "emit_argweaver_sites: true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode == 0, result.stderr

    sites = tmp_path / "results" / "sites" / "combined.1.sites"
    assert sites.exists()
    lines = [ln for ln in sites.read_text(encoding="utf-8").splitlines() if ln]
    assert lines[0] == "NAMES\ts1\ts2"
    assert lines[1] == "REGION\t1\t1\t4"
    site_positions = [ln.split("\t")[0] for ln in lines[2:]]

    variants = tmp_path / "results" / "sites" / "combined.1.vcf"
    variant_positions = [
        ln.split("\t")[1]
        for ln in variants.read_text(encoding="utf-8").splitlines()
        if ln and not ln.startswith("#")
    ]
    assert site_positions == variant_positions
    assert lines[2:] == ["2\tCG"]


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_omits_argweaver_sites_by_default(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "ACGT", "s2", "AGGT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["s1", "s2"]',
                "max_missing_count: 0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "results" / "sites" / "combined.1.sites").exists()


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_default_target_emits_argweaver_sites_after_enabling(tmp_path: Path) -> None:
    """Enabling `.sites` in an existing results directory must produce the file.

    On a *fresh* run `.sites` appears regardless of `rule all`, because it is a
    declared output of `direct_maf_sites` and that rule has to run anyway. The
    case that depends on `_all_targets` listing `.sites` is enabling the flag
    when the other outputs are already up to date: without it, `direct_maf_sites`
    is not re-run and the file is never written. Since v1.9 removed `.sites` from
    the summary_report inputs, `rule all` is the only thing covering that path.
    """
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "ACGT", "s2", "AGGT")

    base = [
        "maf_dir: maf",
        "reference_fasta: ref.fa",
        "results_dir: results",
        'samples: ["s1", "s2"]',
        "max_missing_count: 0",
    ]
    config = tmp_path / "config.yaml"
    sites = tmp_path / "results" / "sites" / "combined.1.sites"

    # No explicit target throughout: exercises the default `all` rule.
    _write_config(config, *base, "emit_argweaver_sites: false")
    result = _run_snakemake(tmp_path, config)
    assert result.returncode == 0, result.stderr
    assert not sites.exists()

    _write_config(config, *base, "emit_argweaver_sites: true")
    result = _run_snakemake(tmp_path, config)
    assert result.returncode == 0, result.stderr
    assert sites.exists(), "rule all no longer pulls in the ARGweaver .sites output"

    lines = sites.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "NAMES\ts1\ts2"
    assert lines[1] == "REGION\t1\t1\t4"
    # One variant site: s2 carries G where the reference and s1 carry C.
    assert lines[2:] == ["2\tCG"]


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_disabling_argweaver_sites_leaves_stale_output(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "ACGT", "s2", "AGGT")

    config = tmp_path / "config.yaml"
    base_config = [
        "maf_dir: maf",
        "reference_fasta: ref.fa",
        "results_dir: results",
        'samples: ["s1", "s2"]',
        "max_missing_count: 0",
    ]
    config.write_text("\n".join(base_config + ["emit_argweaver_sites: true"]) + "\n", encoding="utf-8")

    target = str(tmp_path / "results" / "summary.html")
    result = _run_snakemake(tmp_path, config, target)
    assert result.returncode == 0, result.stderr
    sites = tmp_path / "results" / "sites" / "combined.1.sites"
    assert sites.exists()

    config.write_text("\n".join(base_config + ["emit_argweaver_sites: false"]) + "\n", encoding="utf-8")
    result = _run_snakemake(tmp_path, config, target)
    assert result.returncode == 0, result.stderr
    assert sites.exists()


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_enabling_argweaver_sites_creates_output_when_targeted(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "ACGT", "s2", "AGGT")

    config = tmp_path / "config.yaml"
    base_config = [
        "maf_dir: maf",
        "reference_fasta: ref.fa",
        "results_dir: results",
        'samples: ["s1", "s2"]',
        "max_missing_count: 0",
    ]
    config.write_text("\n".join(base_config + ["emit_argweaver_sites: false"]) + "\n", encoding="utf-8")

    target = str(tmp_path / "results" / "summary.html")
    result = _run_snakemake(tmp_path, config, target)
    assert result.returncode == 0, result.stderr
    sites = tmp_path / "results" / "sites" / "combined.1.sites"
    assert not sites.exists()

    config.write_text("\n".join(base_config + ["emit_argweaver_sites: true"]) + "\n", encoding="utf-8")
    result = _run_snakemake(tmp_path, config, str(sites))
    assert result.returncode == 0, result.stderr
    assert sites.exists()


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_handles_sample_names_with_spaces(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nAC\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "sample one.maf", "1", "AC", "sample_one", "AT")

    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(
            [
                "maf_dir: maf",
                "reference_fasta: ref.fa",
                "results_dir: results",
                'samples: ["sample one"]',
                "max_missing_count: 0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode == 0, result.stderr

    all_sites = tmp_path / "results" / "sites" / "combined.1.all_sites.vcf"
    assert next(
        line for line in all_sites.read_text(encoding="utf-8").splitlines()
        if line.startswith("#CHROM\t")
    ).endswith("\tsample one")


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_workflow_requires_explicit_configfile(tmp_path: Path) -> None:
    result = _run_snakemake(tmp_path, None, str(tmp_path / "results" / "summary.html"))
    assert result.returncode != 0
    assert "A config file is required" in (result.stderr + result.stdout)


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
@pytest.mark.parametrize("missing_key", ["maf_dir", "reference_fasta"])
def test_workflow_requires_expected_config_keys(tmp_path: Path, missing_key: str) -> None:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGT\n", encoding="utf-8")

    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGT", "s1", "ACGT")

    config_lines = [
        "maf_dir: maf",
        "reference_fasta: ref.fa",
        "results_dir: results",
        'samples: ["s1"]',
    ]
    config = tmp_path / "config.yaml"
    config.write_text(
        "\n".join(line for line in config_lines if not line.startswith(f"{missing_key}:"))
        + "\n",
        encoding="utf-8",
    )

    result = _run_snakemake(tmp_path, config, str(tmp_path / "results" / "summary.html"))
    assert result.returncode != 0
    assert f"Missing required config keys: {missing_key}" in (result.stderr + result.stdout)


def _maf_qc_fixture(tmp_path: Path, *extra_config: str) -> Path:
    ref = tmp_path / "ref.fa"
    ref.write_text(">1\nACGTACGTAC\n", encoding="utf-8")
    maf_dir = tmp_path / "maf"
    maf_dir.mkdir()
    _write_pairwise_maf(maf_dir / "s1.maf", "1", "ACGTACGTAC", "s1", "ACGTACGTAC")
    _write_pairwise_maf(maf_dir / "s2.maf", "1", "ACGTACGTAC", "s2", "ACGAACGTAC")
    config = tmp_path / "config.yaml"
    _write_config(
        config,
        "maf_dir: maf",
        "reference_fasta: ref.fa",
        "results_dir: results",
        'samples: ["s1", "s2"]',
        "maf_stats_dotplots: all",
        *extra_config,
    )
    return config


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_maf_qc_target_builds_stats_report_and_dotplots_only(tmp_path: Path) -> None:
    config = _maf_qc_fixture(tmp_path)
    result = _run_snakemake(tmp_path, config, "maf_qc")
    assert result.returncode == 0, result.stderr

    stats = tmp_path / "results" / "maf_stats"
    for sample in ("s1", "s2"):
        for suffix in ("maf_stats.tsv", "by_reference_contig.tsv", "by_query_contig.tsv"):
            assert (stats / f"{sample}.{suffix}").exists()
        assert (stats / "dotplots" / sample / "1.png").exists()
    rows = (stats / "maf_stats.tsv").read_text(encoding="utf-8").splitlines()
    assert len(rows) == 3 and rows[0].split("\t")[-1] == "flags"
    html = (stats / "maf_stats.html").read_text(encoding="utf-8")
    assert "dotplots/s1/1.png" in html
    # QC never starts the main analysis.
    assert not (tmp_path / "results" / "sites").exists()
    assert not (tmp_path / "results" / "maf_by_contig").exists()


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_default_target_does_not_run_maf_qc(tmp_path: Path) -> None:
    config = _maf_qc_fixture(tmp_path)
    result = _run_snakemake(tmp_path, config)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "results" / "summary.html").exists()
    assert not (tmp_path / "results" / "maf_stats").exists()


@pytest.mark.skipif(
    importlib.util.find_spec("snakemake") is None,
    reason="snakemake is not installed in the test environment",
)
def test_maf_qc_uses_query_fai_dir_and_requires_every_sample(tmp_path: Path) -> None:
    fai_dir = tmp_path / "query_fai"
    fai_dir.mkdir()
    (fai_dir / "s1.fa.fai").write_text("s1\t10\t4\t10\t11\nscaffold\t30\t20\t30\t31\n", encoding="utf-8")
    config = _maf_qc_fixture(tmp_path, "query_fai_dir: query_fai")

    missing = _run_snakemake(tmp_path, config, "maf_qc")
    assert missing.returncode != 0
    assert "No query .fai for sample 's2'" in (missing.stderr + missing.stdout)

    (fai_dir / "s2.fai").write_text("s2\t10\t4\t10\t11\n", encoding="utf-8")
    result = _run_snakemake(tmp_path, config, "maf_qc")
    assert result.returncode == 0, result.stderr
    header, row = (tmp_path / "results" / "maf_stats" / "s1.maf_stats.tsv").read_text(
        encoding="utf-8"
    ).splitlines()
    values = dict(zip(header.split("\t"), row.split("\t")))
    assert values["query_length_source"] == "query_fai"
    assert values["query_length_bp"] == "40"
