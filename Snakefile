import re
import sys
import shlex
from pathlib import Path

from snakemake.io import glob_wildcards

from scripts.common import normalize_contig, read_maf_contigs

if not workflow.configfiles:
    raise ValueError(
        "A config file is required. Run Snakemake with --configfile path/to/options.yaml."
    )

wildcard_constraints:
    contig="[^/]+"


def _config_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value != 0
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"1", "true", "t", "yes", "y", "on"}:
            return True
        if normalized in {"0", "false", "f", "no", "n", "off", ""}:
            return False
    raise ValueError(f"Invalid boolean config value: {value!r}")


REQUIRED_CONFIG_KEYS = ("maf_dir", "reference_fasta")
missing_config_keys = [key for key in REQUIRED_CONFIG_KEYS if key not in config]
if missing_config_keys:
    raise ValueError(
        "Missing required config keys: "
        + ", ".join(missing_config_keys)
        + ". Check the file passed via --configfile."
    )


MAF_DIR = Path(config["maf_dir"]).resolve()
ORIG_REF_FASTA = Path(config["reference_fasta"]).resolve()
RESULTS_DIR = Path(config.get("results_dir", "results")).resolve()

ALLOW_MULTIALLELIC = _config_bool(config.get("allow_multiallelic_snps", True))
MASK_INDEL_ADJACENT_SNPS = _config_bool(config.get("mask_indel_adjacent_snps", False))
ADD_REF = _config_bool(config.get("add_ref", False))
EMIT_ARGWEAVER_SITES = _config_bool(config.get("emit_argweaver_sites", False))
MAX_MISSING_COUNT = config.get("max_missing_count")
MAX_MISSING_FRACTION = config.get("max_missing_fraction")

QUALITY_BED_DIR = config.get("quality_bed_dir")
if QUALITY_BED_DIR in (None, ""):
    QUALITY_BED_DIR = None
else:
    QUALITY_BED_DIR = Path(QUALITY_BED_DIR).resolve()
QUALITY_MIN = config.get("quality_min")
if QUALITY_MIN not in (None, ""):
    QUALITY_MIN = float(QUALITY_MIN)
    if not 0 <= QUALITY_MIN <= 1:
        raise ValueError(f"quality_min must be between 0 and 1; got {QUALITY_MIN}")
    if QUALITY_BED_DIR is None:
        raise ValueError("quality_min is set but quality_bed_dir is missing; set both or neither.")
elif QUALITY_BED_DIR is not None:
    raise ValueError("quality_bed_dir is set but quality_min is missing; set both or neither.")
else:
    QUALITY_MIN = None

DEFAULT_MEM_MB = 48000
DEFAULT_TIME = "24:00:00"
SUMMARY_WINDOW_BP = int(config.get("summary_window_bp", 100000))
if SUMMARY_WINDOW_BP <= 0:
    raise ValueError(
        f"summary_window_bp must be a positive integer; got {SUMMARY_WINDOW_BP}"
    )

# First-stage MAF QC (`snakemake ... maf_qc`); never part of `rule all`.
MAF_STATS_DIR = RESULTS_DIR / "maf_stats"
QUERY_FAI_DIR = config.get("query_fai_dir")
QUERY_FAI_DIR = None if QUERY_FAI_DIR in (None, "") else Path(QUERY_FAI_DIR).resolve()
QUERY_FASTA_DIR = config.get("query_fasta_dir")
QUERY_FASTA_DIR = None if QUERY_FASTA_DIR in (None, "") else Path(QUERY_FASTA_DIR).resolve()
if QUERY_FASTA_DIR is not None and QUERY_FAI_DIR is not None:
    raise ValueError(
        "Set query_fasta_dir or query_fai_dir, not both; the FASTA already supplies contig lengths."
    )
MAF_STATS_BREAKPOINT_CONTEXT_BP = int(config.get("maf_stats_breakpoint_context_bp", 1000000))
MAF_STATS_RECURRENCE_WINDOW_BP = int(config.get("maf_stats_recurrence_window_bp", 500000))
MAF_STATS_MIN_BLOCK_BP = int(config.get("maf_stats_min_block_bp", 0))
MAF_STATS_OVERLAP_TOLERANCE_BP = int(config.get("maf_stats_overlap_tolerance_bp", 0))
MAF_STATS_DOTPLOTS = str(config.get("maf_stats_dotplots", "flagged"))
if MAF_STATS_DOTPLOTS in ("False", "false", "0"):  # YAML `false` arrives as a bool
    MAF_STATS_DOTPLOTS = "false"
if MAF_STATS_DOTPLOTS not in ("flagged", "all", "false"):
    raise ValueError(
        f"maf_stats_dotplots must be flagged, all, or false; got {MAF_STATS_DOTPLOTS!r}"
    )
MAF_STATS_DOTPLOT_MAX = int(config.get("maf_stats_dotplot_max", 20))
MAF_STATS_FLAG_OPTIONS = {
    "--flag-min-aligned-reference": config.get("maf_stats_flag_min_aligned_reference"),
    "--flag-min-identity": config.get("maf_stats_flag_min_identity"),
    "--flag-max-breakpoints-per-gb": config.get("maf_stats_flag_max_breakpoints_per_gb"),
}
# Placeholder resources until real MAFs are benchmarked; override in the config.
MAF_STATS_DEFAULT_MEM_MB = 16000
MAF_STATS_DEFAULT_TIME = "12:00:00"

DIRECT_REF_FASTA = RESULTS_DIR / "refs" / "reference_sites.fa"
REF_FAI = str(DIRECT_REF_FASTA) + ".fai"
MAF_CHUNK_ROOT = RESULTS_DIR / "maf_by_contig"

def _maf_path_for_sample(sample: str) -> Path:
    maf = MAF_DIR / f"{sample}.maf"
    maf_gz = MAF_DIR / f"{sample}.maf.gz"
    if maf.exists():
        return maf
    if maf_gz.exists():
        return maf_gz
    return maf


def _quality_bed_for_sample(sample: str) -> Path | None:
    if QUALITY_BED_DIR is None:
        return None
    plain = QUALITY_BED_DIR / f"{sample}.bed"
    gz = QUALITY_BED_DIR / f"{sample}.bed.gz"
    if plain.exists():
        return plain
    if gz.exists():
        return gz
    return None


QUERY_FAI_SUFFIXES = (".fai", ".fa.fai", ".fasta.fai", ".fna.fai", ".fa.gz.fai", ".fasta.gz.fai")


def _query_fai_for_sample(sample: str) -> list[str]:
    """`<sample>` + one of QUERY_FAI_SUFFIXES under query_fai_dir; [] when the
    option is unset. Missing or ambiguous files are errors, so query coverage
    never silently mixes denominator sources across samples."""
    if QUERY_FAI_DIR is None:
        return []
    found = [QUERY_FAI_DIR / f"{sample}{suffix}" for suffix in QUERY_FAI_SUFFIXES]
    found = [path for path in found if path.exists()]
    if len(found) != 1:
        problem = "No" if not found else "Multiple"
        raise ValueError(
            f"{problem} query .fai for sample '{sample}' in {QUERY_FAI_DIR} "
            f"(looked for {sample}{{{','.join(QUERY_FAI_SUFFIXES)}}})"
        )
    return [str(found[0])]


QUERY_FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")


def _query_fasta_for_sample(sample: str) -> list[str]:
    """`<sample>` + one of QUERY_FASTA_SUFFIXES under query_fasta_dir; [] when
    the option is unset. Missing or ambiguous files are errors."""
    if QUERY_FASTA_DIR is None:
        return []
    found = [QUERY_FASTA_DIR / f"{sample}{suffix}" for suffix in QUERY_FASTA_SUFFIXES]
    found = [path for path in found if path.exists()]
    if len(found) != 1:
        problem = "No" if not found else "Multiple"
        raise ValueError(
            f"{problem} query FASTA for sample '{sample}' in {QUERY_FASTA_DIR} "
            f"(looked for {sample}{{{','.join(QUERY_FASTA_SUFFIXES)}}})"
        )
    return [str(found[0])]


def _quality_bed_inputs():
    if QUALITY_BED_DIR is None:
        return []
    beds = []
    for sample in SAMPLES:
        bed = _quality_bed_for_sample(sample)
        if bed is not None:
            beds.append(str(bed))
    return beds


def _discover_samples():
    if "samples" in config:
        return list(config["samples"])
    # Constrain {sample} so it cannot span "/": glob_wildcards otherwise
    # descends into subdirectories of maf_dir (e.g. an example_data/ tree),
    # pulling in nested MAFs as bogus samples.
    maf_pattern = str(MAF_DIR / "{sample,[^/]+}.maf")
    maf_gz_pattern = str(MAF_DIR / "{sample,[^/]+}.maf.gz")
    samples = set(glob_wildcards(maf_pattern).sample)
    samples.update(glob_wildcards(maf_gz_pattern).sample)
    return sorted(samples)


def _read_maf_contig_sets(samples: list[str]) -> dict[str, set[str]]:
    contigs_by_sample: dict[str, set[str]] = {}
    for sample in samples:
        contigs_by_sample[sample] = read_maf_contigs(_maf_path_for_sample(sample))
    return contigs_by_sample


def _read_fai_contigs(fai: Path) -> list[str]:
    contigs = []
    with fai.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            contigs.append(line.split("\t", 1)[0])
    return contigs


def _resolve_requested_contigs(
    requested: list[str], available: list[str]
) -> tuple[list[str], list[str], list[tuple[str, str]]]:
    available_set = set(available)
    available_norm: dict[str, list[str]] = {}
    for name in available:
        available_norm.setdefault(normalize_contig(name), []).append(name)

    kept: list[str] = []
    dropped: list[str] = []
    remapped: list[tuple[str, str]] = []
    seen: set[str] = set()
    for raw in requested:
        req = str(raw)
        mapped = req
        if req in available_set:
            mapped = req
        else:
            candidates = available_norm.get(normalize_contig(req), [])
            if len(candidates) == 1:
                mapped = candidates[0]
                remapped.append((req, mapped))
            else:
                dropped.append(req)
                continue
        if mapped not in seen:
            kept.append(mapped)
            seen.add(mapped)
    return kept, dropped, remapped


SAMPLES = _discover_samples()
if not SAMPLES:
    raise ValueError(f"No MAF files found in {MAF_DIR}")

_MAF_CONTIG_INTERSECTION_CACHE: list[str] | None = None


def _maf_contig_intersection() -> list[str]:
    global _MAF_CONTIG_INTERSECTION_CACHE
    if _MAF_CONTIG_INTERSECTION_CACHE is None:
        contigs_by_sample = _read_maf_contig_sets(SAMPLES)
        if contigs_by_sample:
            normalized_sets = [
                {normalize_contig(contig) for contig in contigs}
                for contigs in contigs_by_sample.values()
            ]
            _MAF_CONTIG_INTERSECTION_CACHE = sorted(
                set.intersection(*normalized_sets)
            )
        else:
            _MAF_CONTIG_INTERSECTION_CACHE = []
    return _MAF_CONTIG_INTERSECTION_CACHE


def _active_contig_resolution() -> tuple[list[str], list[str], list[str], list[tuple[str, str]]]:
    ckpt = checkpoints.index_reference.get()
    fai = Path(str(ckpt.output.fai))
    available = _read_fai_contigs(fai)
    if "contigs" in config:
        requested = [str(c) for c in config["contigs"]]
        kept, dropped, remapped = _resolve_requested_contigs(requested, available)
        if dropped:
            raise ValueError(
                "Configured contigs must each have an exact or unambiguous normalized "
                "match in reference .fai; unmatched or ambiguous: "
                + ", ".join(dropped[:10])
            )
        return kept, dropped, requested, remapped

    requested = list(_maf_contig_intersection())
    if not requested:
        raise ValueError(
            "No contigs are shared across all MAF files. Set explicit 'contigs' in options.yaml to override."
        )
    kept, dropped, remapped = _resolve_requested_contigs(requested, available)
    if not kept:
        raise ValueError(
            "No shared MAF contigs are present in reference .fai. Set explicit 'contigs' in options.yaml to override."
        )
    return kept, dropped, requested, remapped


_CONTIG_RESOLUTION_LOGGED = False


def _active_contigs() -> list[str]:
    global _CONTIG_RESOLUTION_LOGGED
    kept, dropped, _requested, remapped = _active_contig_resolution()
    if not _CONTIG_RESOLUTION_LOGGED:
        _CONTIG_RESOLUTION_LOGGED = True
        if remapped:
            pairs = ", ".join(f"{r}->{m}" for r, m in remapped)
            print(f"[argprep] Remapped contigs to reference names: {pairs}", file=sys.stderr)
        if dropped:
            print(
                "[argprep] Skipped contigs (no unambiguous match in reference .fai): "
                + ", ".join(dropped),
                file=sys.stderr,
            )
    return kept


def _maf_input(sample):
    return str(_maf_path_for_sample(sample))


def _direct_prefix(contig):
    return RESULTS_DIR / "sites" / f"combined.{contig}"


def _direct_all_sites_out(contig):
    return Path(str(_direct_prefix(contig)) + ".all_sites.vcf")


def _direct_variants_out(contig):
    return Path(str(_direct_prefix(contig)) + ".vcf")


def _direct_mask_out(contig):
    return Path(str(_direct_prefix(contig)) + ".mask.bed")


def _direct_sites_out(contig):
    return Path(str(_direct_prefix(contig)) + ".sites")


def _direct_report_stats_out(contig):
    return Path(str(_direct_prefix(contig)) + ".report_stats.tsv")


def _direct_sample_missing_mask_out(contig, sample):
    return RESULTS_DIR / "sites" / f"combined.{contig}.{sample}.missing.bed"


def _direct_sample_mask_out(contig):
    return RESULTS_DIR / "sites" / f"combined.{contig}.sample.mask.bed"


def _split_sample_dir(sample):
    return MAF_CHUNK_ROOT / sample


def _split_sample_contig_maf(sample, contig):
    return _split_sample_dir(sample) / f"{contig}.maf.gz"


def _all_targets(_wc):
    contigs = _active_contigs()
    return (
        [str(_direct_all_sites_out(c)) for c in contigs]
        + [str(_direct_variants_out(c)) for c in contigs]
        + [str(_direct_mask_out(c)) for c in contigs]
        + [str(_direct_sample_missing_mask_out(c, s)) for c in contigs for s in SAMPLES]
        + [str(_direct_sample_mask_out(c)) for c in contigs]
        + ([str(_direct_sites_out(c)) for c in contigs] if EMIT_ARGWEAVER_SITES else [])
        + [str(RESULTS_DIR / "summary.html")]
    )


rule all:
    input: _all_targets


ruleorder: combine_sample_missing_masks > direct_maf_sites


rule prepare_reference:
    resources:
        mem_mb=1000,
        time="00:10:00",
    input:
        ref=str(ORIG_REF_FASTA),
    output:
        ref=str(DIRECT_REF_FASTA),
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.ref}")"
        ln -s "$(realpath "{input.ref}")" "{output.ref}"
        """


checkpoint index_reference:
    resources:
        mem_mb=4000,
        time="02:00:00",
    input:
        ref=str(DIRECT_REF_FASTA),
    output:
        fai=REF_FAI,
    shell:
        """
        set -euo pipefail
        samtools faidx "{input.ref}"
        """


rule split_sample_maf:
    input:
        maf=lambda wc: _maf_input(wc.sample),
        fai=REF_FAI,
    output:
        chunks=directory(str(MAF_CHUNK_ROOT / "{sample}")),
    params:
        out_root=str(MAF_CHUNK_ROOT),
        contigs=lambda wc: " ".join(shlex.quote(contig) for contig in _active_contigs()),
    shell:
        """
        set -euo pipefail
        python "{workflow.basedir}/scripts/split_maf_by_contig.py" \
          --maf "{input.maf}" \
          --sample "{wildcards.sample}" \
          --out-root "{params.out_root}" \
          --contigs {params.contigs}
        """


rule direct_maf_sites:
    resources:
        mem_mb=int(config.get("maf_mem_mb", DEFAULT_MEM_MB)),
        time=str(config.get("maf_time", DEFAULT_TIME))
    input:
        mafs=lambda wc: [str(_split_sample_dir(sample)) for sample in SAMPLES],
        ref=str(DIRECT_REF_FASTA),
        fai=REF_FAI,
        quality_beds=lambda wc: _quality_bed_inputs(),
    output:
        all_sites=str(RESULTS_DIR / "sites" / "combined.{contig}.all_sites.vcf"),
        variants=str(RESULTS_DIR / "sites" / "combined.{contig}.vcf"),
        mask=str(RESULTS_DIR / "sites" / "combined.{contig}.mask.bed"),
        summary=str(RESULTS_DIR / "sites" / "combined.{contig}.site_summary.tsv"),
        report_stats=str(RESULTS_DIR / "sites" / "combined.{contig}.report_stats.tsv"),
        sample_missing_masks=expand(
            str(RESULTS_DIR / "sites" / "combined.{{contig}}.{sample}.missing.bed"),
            sample=SAMPLES,
        ),
        **(
            {"sites": str(RESULTS_DIR / "sites" / "combined.{contig}.sites")}
            if EMIT_ARGWEAVER_SITES
            else {}
        ),
    params:
        maf_dir=str(MAF_DIR),
        maf_paths=lambda wc: " ".join(
            shlex.quote(f"{sample}={_split_sample_contig_maf(sample, wc.contig)}")
            for sample in SAMPLES
        ),
        samples=" ".join(shlex.quote(sample) for sample in SAMPLES),
        max_missing_count=(
            None if MAX_MISSING_COUNT in (None, "") else int(MAX_MISSING_COUNT)
        ),
        max_missing_fraction=(
            None
            if MAX_MISSING_FRACTION in (None, "")
            else float(MAX_MISSING_FRACTION)
        ),
        allow_multiallelic=ALLOW_MULTIALLELIC,
        mask_indel_adjacent_snps=MASK_INDEL_ADJACENT_SNPS,
        add_ref=ADD_REF,
        emit_argweaver_sites=EMIT_ARGWEAVER_SITES,
        quality_bed_dir=("" if QUALITY_BED_DIR is None else str(QUALITY_BED_DIR)),
        quality_min=("" if QUALITY_MIN is None else QUALITY_MIN),
        out_prefix=lambda wc: str(_direct_prefix(wc.contig)),
        window_bp=SUMMARY_WINDOW_BP,
    shell:
        """
        set -euo pipefail
        mkdir -p "{RESULTS_DIR}/sites"
        cmd=(python "{workflow.basedir}/scripts/maf_to_sites.py"
          --maf-dir "{params.maf_dir}"
          --reference-fasta "{input.ref}"
          --contig "{wildcards.contig}"
          --out-prefix "{params.out_prefix}"
          --window-bp "{params.window_bp}"
          --maf-paths {params.maf_paths}
          --samples {params.samples})
        if [ "{params.max_missing_count}" != "None" ]; then
          cmd+=(--max-missing-count "{params.max_missing_count}")
        fi
        if [ "{params.max_missing_fraction}" != "None" ]; then
          cmd+=(--max-missing-fraction "{params.max_missing_fraction}")
        fi
        if [ "{params.allow_multiallelic}" = "True" ]; then
          cmd+=(--allow-multiallelic-snps)
        else
          cmd+=(--mask-multiallelic-snps)
        fi
        if [ "{params.mask_indel_adjacent_snps}" = "True" ]; then
          cmd+=(--mask-indel-adjacent-snps)
        fi
        if [ "{params.add_ref}" = "True" ]; then
          cmd+=(--add-ref)
        fi
        if [ "{params.emit_argweaver_sites}" = "True" ]; then
          cmd+=(--emit-argweaver-sites)
        fi
        if [ -n "{params.quality_bed_dir}" ]; then
          cmd+=(--quality-bed-dir "{params.quality_bed_dir}" --quality-min "{params.quality_min}")
        fi
        "${{cmd[@]}}"
        """


rule combine_sample_missing_masks:
    input:
        lambda wc: [
            str(_direct_sample_missing_mask_out(wc.contig, sample))
            for sample in SAMPLES
        ],
    output:
        str(RESULTS_DIR / "sites" / "combined.{contig}.sample.mask.bed"),
    shell:
        "cat {input:q} > {output:q}"


rule summary_report:
    input:
        report_stats=lambda wc: [str(_direct_report_stats_out(c)) for c in _active_contigs()],
        summaries=lambda wc: [str(_direct_prefix(c)) + ".site_summary.tsv" for c in _active_contigs()],
        sample_missing_beds=lambda wc: [
            str(_direct_sample_missing_mask_out(c, s))
            for c in _active_contigs()
            for s in SAMPLES
        ],
        fai=REF_FAI,
        options_yaml=str(Path(workflow.configfiles[0]).resolve()),
    output:
        report=str(RESULTS_DIR / "summary.html"),
    params:
        window_bp=SUMMARY_WINDOW_BP,
    shell:
        """
        set -euo pipefail
        python "{workflow.basedir}/scripts/summary_report.py" \
          --fai "{input.fai}" \
          --window-bp "{params.window_bp}" \
          --report-out "{output.report}" \
          --report-stats {input.report_stats} \
          --site-summaries {input.summaries} \
          --sample-missing-beds {input.sample_missing_beds} \
          --options-yaml "{input.options_yaml}"
        """


rule maf_stats:
    threads: int(config.get("maf_stats_threads", 1))
    resources:
        mem_mb=int(config.get("maf_stats_mem_mb") or MAF_STATS_DEFAULT_MEM_MB),
        time=str(config.get("maf_stats_time") or MAF_STATS_DEFAULT_TIME),
    wildcard_constraints:
        sample="|".join(re.escape(sample) for sample in SAMPLES),
    input:
        maf=lambda wc: _maf_input(wc.sample),
        fai=REF_FAI,
        query_fai=lambda wc: _query_fai_for_sample(wc.sample),
        query_fasta=lambda wc: _query_fasta_for_sample(wc.sample),
    output:
        summary=str(MAF_STATS_DIR / "{sample}.maf_stats.tsv"),
        by_reference=str(MAF_STATS_DIR / "{sample}.by_reference_contig.tsv"),
        by_query=str(MAF_STATS_DIR / "{sample}.by_query_contig.tsv"),
        breakpoints=str(MAF_STATS_DIR / "{sample}.breakpoints.tsv"),
        nested=str(MAF_STATS_DIR / "{sample}.nested_blocks.tsv"),
        dotplots=directory(str(MAF_STATS_DIR / "dotplots" / "{sample}")),
    params:
        out_dir=str(MAF_STATS_DIR),
        min_block_bp=MAF_STATS_MIN_BLOCK_BP,
        overlap_tolerance_bp=MAF_STATS_OVERLAP_TOLERANCE_BP,
        dotplots=MAF_STATS_DOTPLOTS,
        dotplot_max=MAF_STATS_DOTPLOT_MAX,
        breakpoint_context_bp=MAF_STATS_BREAKPOINT_CONTEXT_BP,
    shell:
        """
        set -euo pipefail
        cmd=(python "{workflow.basedir}/scripts/maf_stats.py"
          --maf "{input.maf}"
          --reference-fai "{input.fai}"
          --sample "{wildcards.sample}"
          --out-dir "{params.out_dir}"
          --min-block-bp "{params.min_block_bp}"
          --overlap-tolerance-bp "{params.overlap_tolerance_bp}"
          --dotplots "{params.dotplots}"
          --dotplot-max "{params.dotplot_max}"
          --breakpoint-context-bp "{params.breakpoint_context_bp}")
        if [ -n "{input.query_fasta}" ]; then
          cmd+=(--fasta "{input.query_fasta}")
        elif [ -n "{input.query_fai}" ]; then
          cmd+=(--query-fai "{input.query_fai}")
        fi
        "${{cmd[@]}}"
        """


rule maf_stats_report:
    resources:
        mem_mb=int(config.get("maf_stats_report_mem_mb", 4000)),
        time=str(config.get("maf_stats_report_time", "00:30:00")),
    input:
        summaries=expand(str(MAF_STATS_DIR / "{sample}.maf_stats.tsv"), sample=SAMPLES),
        by_reference=expand(str(MAF_STATS_DIR / "{sample}.by_reference_contig.tsv"), sample=SAMPLES),
        by_query=expand(str(MAF_STATS_DIR / "{sample}.by_query_contig.tsv"), sample=SAMPLES),
        breakpoints=expand(str(MAF_STATS_DIR / "{sample}.breakpoints.tsv"), sample=SAMPLES),
        nested=expand(str(MAF_STATS_DIR / "{sample}.nested_blocks.tsv"), sample=SAMPLES),
        dotplots=expand(str(MAF_STATS_DIR / "dotplots" / "{sample}"), sample=SAMPLES),
    output:
        tsv=str(MAF_STATS_DIR / "maf_stats.tsv"),
        html=str(MAF_STATS_DIR / "maf_stats.html"),
        breakpoints=str(MAF_STATS_DIR / "maf_stats.breakpoints.tsv"),
        nested=str(MAF_STATS_DIR / "maf_stats.nested_blocks.tsv"),
    params:
        recurrence_window_bp=MAF_STATS_RECURRENCE_WINDOW_BP,
        flag_args=" ".join(
            f"{option} {shlex.quote(str(value))}"
            for option, value in MAF_STATS_FLAG_OPTIONS.items()
            if value not in (None, "")
        ),
    shell:
        """
        set -euo pipefail
        python "{workflow.basedir}/scripts/maf_stats_report.py" \
          --summaries {input.summaries:q} \
          --by-reference {input.by_reference:q} \
          --by-query {input.by_query:q} \
          --breakpoints {input.breakpoints:q} \
          --out-breakpoints-tsv "{output.breakpoints}" \
          --nested-blocks {input.nested:q} \
          --out-nested-tsv "{output.nested}" \
          --recurrence-window-bp "{params.recurrence_window_bp}" \
          --out-tsv "{output.tsv}" \
          --out-html "{output.html}" \
          {params.flag_args}
        """


# QC-only target: run and inspect before launching the main workflow.
rule maf_qc:
    input:
        str(MAF_STATS_DIR / "maf_stats.tsv"),
        str(MAF_STATS_DIR / "maf_stats.html"),
        str(MAF_STATS_DIR / "maf_stats.breakpoints.tsv"),
        str(MAF_STATS_DIR / "maf_stats.nested_blocks.tsv"),
        expand(str(MAF_STATS_DIR / "{sample}.breakpoints.tsv"), sample=SAMPLES),
        expand(str(MAF_STATS_DIR / "{sample}.nested_blocks.tsv"), sample=SAMPLES),
        expand(str(MAF_STATS_DIR / "{sample}.by_reference_contig.tsv"), sample=SAMPLES),
        expand(str(MAF_STATS_DIR / "{sample}.by_query_contig.tsv"), sample=SAMPLES),
        expand(str(MAF_STATS_DIR / "dotplots" / "{sample}"), sample=SAMPLES),
