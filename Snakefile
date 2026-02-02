import os
from pathlib import Path

from snakemake.io import glob_wildcards

configfile: "config.yaml"

wildcard_constraints:
    sample="[^/]+",
    gvcf_base="[^/]+",
    contig="[^/]+"

MAF_DIR = Path(config["maf_dir"]).resolve()
REF_FASTA = Path(config["reference_fasta"]).resolve()
GVCF_DIR = Path(config.get("gvcf_dir", "gvcf")).resolve()
RESULTS_DIR = Path(config.get("results_dir", "results")).resolve()
TASSEL_DIR = Path(config.get("tassel_dir", "tassel-5-standalone")).resolve()

SAMPLE_SUFFIX = config.get("sample_suffix", "_anchorwave")
FILL_GAPS = str(config.get("fill_gaps", "false")).lower()
DEPTH = int(config["depth"])
DROP_CUTOFF = config.get("drop_cutoff", "")
FILTER_MULTIALLELIC = bool(config.get("filter_multiallelic", False))
GZIP_OUTPUT = bool(config.get("gzip_output", False))
NO_MERGE = bool(config.get("no_merge", False))

REF_BASE = REF_FASTA.name
for ext in (".fa", ".fasta"):
    if REF_BASE.endswith(ext):
        REF_BASE = REF_BASE[: -len(ext)]

REF_FAI = str(REF_FASTA) + ".fai"
REF_DICT = str(REF_FASTA.with_suffix(".dict"))

def _discover_samples():
    if "samples" in config:
        return list(config["samples"])
    pattern = str(MAF_DIR / "{sample}.maf")
    return sorted(set(glob_wildcards(pattern).sample))


def _read_contigs():
    if "contigs" in config:
        return list(config["contigs"])
    fai = REF_FASTA.with_suffix(REF_FASTA.suffix + ".fai")
    if not fai.exists():
        raise ValueError(
            "Reference FASTA index (.fai) not found. "
            "Either run 'samtools faidx' on the reference or set 'contigs' in config.yaml."
        )
    contigs = []
    with fai.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            contigs.append(line.split("\t", 1)[0])
    return contigs


SAMPLES = _discover_samples()
CONTIGS = _read_contigs()

GVCF_BASES = [f"{sample}To{REF_BASE}" for sample in SAMPLES]


def _gvcf_out(base):
    return GVCF_DIR / f"{base}.gvcf.gz"


def _split_out(base, contig):
    return GVCF_DIR / "cleangVCF" / "split_gvcf" / f"{base}.{contig}.gvcf.gz"


def _combined_out(contig):
    return RESULTS_DIR / "combined" / f"combined.{contig}.gvcf.gz"


def _split_prefix(contig):
    return RESULTS_DIR / "split" / f"combined.{contig}"


SPLIT_SUFFIX = ".gz" if GZIP_OUTPUT else ""

rule all:
    # Final targets: merged gVCFs plus filtered bed masks per contig.
    input:
        [str(_combined_out(c)) for c in CONTIGS],
        [str(_split_prefix(c)) + ".filtered.bed" for c in CONTIGS],

rule index_reference:
    # Create reference FASTA index and sequence dictionary for GATK.
    input:
        ref=str(REF_FASTA),
    output:
        fai=REF_FAI,
        dict=REF_DICT,
    shell:
        """
        set -euo pipefail
        samtools faidx "{input.ref}"
        picard CreateSequenceDictionary R="{input.ref}" O="{output.dict}"
        """


rule maf_to_gvcf:
    # Convert each MAF to a gzipped gVCF using TASSEL.
    threads: 2
    resources:
        mem_mb=256000,
        time="24:00:00"
    input:
        maf=str(MAF_DIR / "{sample}.maf"),
        ref=str(REF_FASTA),
    output:
        gvcf=str(GVCF_DIR / (f"{{sample}}To{REF_BASE}.gvcf.gz")),
        tbi=str(GVCF_DIR / (f"{{sample}}To{REF_BASE}.gvcf.gz.tbi")),
    log:
        str(Path("logs") / "tassel" / "{sample}.log"),
    params:
        tassel_dir=str(TASSEL_DIR),
        sample_name=lambda wc: f"{wc.sample}{SAMPLE_SUFFIX}",
        fill_gaps=FILL_GAPS,
    shell:
        """
        set -euo pipefail
        mkdir -p "{GVCF_DIR}"
        mkdir -p "$(dirname "{log}")"
        out_base="{output.gvcf}"
        if [[ "$out_base" == *.gz ]]; then
          out_base="${{out_base%.gz}}"
        fi
        "{params.tassel_dir}/run_pipeline.pl" -Xmx256G -debug \
          -MAFToGVCFPlugin \
          -referenceFasta "{input.ref}" \
          -mafFile "{input.maf}" \
          -sampleName "{params.sample_name}" \
          -gvcfOutput "$out_base" \
          -fillGaps "{params.fill_gaps}" \
          > "{log}" 2>&1
        if [ ! -f "{output.gvcf}" ]; then
          echo "ERROR: TASSEL did not write {output.gvcf}" >&2
          tail -n 200 "{log}" >&2 || true
          exit 1
        fi
        if [ ! -f "{output.tbi}" ]; then
          tabix -p vcf "{output.gvcf}"
        fi
        """


rule drop_sv:
    # Remove large indels from all gVCFs in the directory.
    input:
        gvcfs=[str(_gvcf_out(b)) for b in GVCF_BASES],
    output:
        bed=str(GVCF_DIR / "cleangVCF" / "dropped_indels.bed"),
        gvcfs=[str(GVCF_DIR / "cleangVCF" / f"{b}.gvcf.gz") for b in GVCF_BASES],
        tbis=[str(GVCF_DIR / "cleangVCF" / f"{b}.gvcf.gz.tbi") for b in GVCF_BASES],
        clean_dir=directory(str(GVCF_DIR / "cleangVCF")),
    params:
        cutoff=DROP_CUTOFF,
    shell:
        """
        set -euo pipefail
        if [ -n "{params.cutoff}" ]; then
          python3 "{workflow.basedir}/dropSV.py" -d "{GVCF_DIR}" -c "{params.cutoff}"
        else
          python3 "{workflow.basedir}/dropSV.py" -d "{GVCF_DIR}"
        fi
        """


rule split_gvcf_by_contig:
    # Split each cleaned gVCF into per-contig gVCFs for merge.
    input:
        gvcf=lambda wc: str(GVCF_DIR / "cleangVCF" / f"{wc.gvcf_base}.gvcf.gz"),
        ref=str(REF_FASTA),
        fai=REF_FAI,
        dict=REF_DICT,
        bed=str(GVCF_DIR / "cleangVCF" / "dropped_indels.bed"),
    output:
        gvcf=str(GVCF_DIR / "cleangVCF" / "split_gvcf" / "{gvcf_base}.{contig}.gvcf.gz"),
        tbi=str(GVCF_DIR / "cleangVCF" / "split_gvcf" / "{gvcf_base}.{contig}.gvcf.gz.tbi"),
    shell:
        """
        set -euo pipefail
        mkdir -p "{GVCF_DIR}/cleangVCF/split_gvcf"
        gatk --java-options "-Xmx100g -Xms100g" SelectVariants \
          -R "{input.ref}" \
          -V "{input.gvcf}" \
          -L "{wildcards.contig}" \
          -O "{output.gvcf}"
        tabix -f -p vcf "{output.gvcf}"
        """


rule merge_contig:
    # Merge all samples for a contig with GenomicsDBImport + GenotypeGVCFs.
    input:
        gvcfs=lambda wc: [str(_split_out(b, wc.contig)) for b in GVCF_BASES],
        ref=str(REF_FASTA),
        fai=REF_FAI,
        dict=REF_DICT,
        bed=str(GVCF_DIR / "cleangVCF" / "dropped_indels.bed"),
    output:
        gvcf=str(RESULTS_DIR / "combined" / "combined.{contig}.gvcf.gz"),
        workspace=directory(str(RESULTS_DIR / "genomicsdb" / "{contig}")),
    params:
        gvcf_args=lambda wc: " ".join(
            f"-V {str(_split_out(b, wc.contig))}" for b in GVCF_BASES
        ),
    shell:
        """
        set -euo pipefail
        mkdir -p "{RESULTS_DIR}/combined"
        gatk --java-options "-Xmx100g -Xms100g" GenomicsDBImport \
          {params.gvcf_args} \
          --genomicsdb-workspace-path "{output.workspace}" \
          -L "{wildcards.contig}"
        gatk --java-options "-Xmx100g -Xms100g" GenotypeGVCFs \
          -R "{input.ref}" \
          -V "gendb://{output.workspace}" \
          -O "{output.gvcf}" \
          -L "{wildcards.contig}"
        """


rule split_gvcf:
    # Split merged contig gVCF into inv/filtered/clean plus missing bed.
    input:
        gvcf=lambda wc: str(_combined_out(wc.contig)),
        ref_fai=REF_FAI,
    output:
        inv=str(RESULTS_DIR / "split" / ("combined.{contig}.inv" + SPLIT_SUFFIX)),
        filt=str(RESULTS_DIR / "split" / ("combined.{contig}.filtered" + SPLIT_SUFFIX)),
        clean=str(RESULTS_DIR / "split" / ("combined.{contig}.clean" + SPLIT_SUFFIX)),
        missing=str(RESULTS_DIR / "split" / ("combined.{contig}.missing.bed" + SPLIT_SUFFIX)),
    params:
        depth=DEPTH,
        filter_multiallelic=FILTER_MULTIALLELIC,
        gzip_output=GZIP_OUTPUT,
        out_prefix=lambda wc: str(_split_prefix(wc.contig)),
    shell:
        """
        set -euo pipefail
        mkdir -p "{RESULTS_DIR}/split"
        cmd=(python3 "{workflow.basedir}/split.py" --depth="{params.depth}" --out-prefix "{params.out_prefix}" --fai "{input.ref_fai}")
        if [ "{params.filter_multiallelic}" = "True" ]; then
          cmd+=(--filter-multiallelic)
        fi
        if [ "{params.gzip_output}" = "True" ]; then
          cmd+=(--gzip-output)
        fi
        cmd+=("{input.gvcf}")
        "${{cmd[@]}}"
        """


rule mask_bed:
    # Build merged mask bed from filtered + missing + dropped indels.
    input:
        filt=lambda wc: str(_split_prefix(wc.contig)) + ".filtered" + SPLIT_SUFFIX,
        missing=lambda wc: str(_split_prefix(wc.contig)) + ".missing.bed" + SPLIT_SUFFIX,
        dropped=str(GVCF_DIR / "cleangVCF" / "dropped_indels.bed"),
    output:
        bed=str(RESULTS_DIR / "split" / "combined.{contig}.filtered.bed"),
    params:
        no_merge=NO_MERGE,
        prefix=lambda wc: str(_split_prefix(wc.contig)),
    shell:
        """
        set -euo pipefail
        cmd=(python3 "{workflow.basedir}/filt_to_bed.py" "{params.prefix}")
        cmd+=(--dropped-bed "{input.dropped}")
        if [ "{params.no_merge}" = "True" ]; then
          cmd+=(--no-merge)
        fi
        "${{cmd[@]}}"
        """
