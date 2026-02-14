import os
import re
import gzip
from pathlib import Path

from snakemake.io import glob_wildcards

configfile: "config.yaml"

wildcard_constraints:
    sample="[^/]+",
    gvcf_base="[^/]+",
    contig="[^/]+"

MAF_DIR = Path(config["maf_dir"]).resolve()
ORIG_REF_FASTA = Path(config["reference_fasta"]).resolve()
GVCF_DIR = Path(config.get("gvcf_dir", "gvcf")).resolve()
RESULTS_DIR = Path(config.get("results_dir", "results")).resolve()
TASSEL_DIR = Path(config.get("tassel_dir", "tassel-5-standalone")).resolve()

SAMPLE_SUFFIX = config.get("sample_suffix", "_anchorwave")
FILL_GAPS = str(config.get("fill_gaps", "false")).lower()
OUTPUT_JUST_GT = bool(config.get("outputJustGT", False))
DROP_CUTOFF = config.get("drop_cutoff", "")
FILTER_MULTIALLELIC = bool(config.get("filter_multiallelic", False))
BGZIP_OUTPUT = bool(config.get("bgzip_output", False))
GENOMICSDB_VCF_BUFFER_SIZE = int(config.get("genomicsdb_vcf_buffer_size", 1048576))
GENOMICSDB_SEGMENT_SIZE = int(config.get("genomicsdb_segment_size", 1048576))
MAF_TO_GVCF_THREADS = int(config.get("maf_to_gvcf_threads", 2))
MAF_TO_GVCF_MEM_MB = int(config.get("maf_to_gvcf_mem_mb", 256000))
MAF_TO_GVCF_TIME = str(config.get("maf_to_gvcf_time", "24:00:00"))
MAF_TO_GVCF_JAVA_MEM_MB = max(256, int(MAF_TO_GVCF_MEM_MB * 0.9))
MERGE_CONTIG_THREADS = int(config.get("merge_contig_threads", config.get("default_threads", 2)))
MERGE_CONTIG_MEM_MB = int(config.get("merge_contig_mem_mb", config.get("default_mem_mb", 48000)))
MERGE_CONTIG_JAVA_MEM_MB = max(256, int(MERGE_CONTIG_MEM_MB * 0.9))
MERGE_CONTIG_TIME = str(config.get("merge_contig_time", config.get("default_time", "48:00:00")))
DEFAULT_MEM_MB = int(config.get("default_mem_mb", 48000))
DEFAULT_JAVA_MEM_MB = max(256, int(DEFAULT_MEM_MB * 0.9))
PLOIDY = int(config.get("ploidy", 2))
VT_NORMALIZE = bool(config.get("vt_normalize", False))
VT_PATH = str(config.get("vt_path", "vt"))
MERGED_GENOTYPER = "selectvariants"

workflow.global_resources["merge_gvcf_jobs"] = int(config.get("merge_gvcf_max_jobs", 4))

REF_BASE = ORIG_REF_FASTA.name
if REF_BASE.endswith(".gz"):
    REF_BASE = REF_BASE[: -len(".gz")]
for ext in (".fa", ".fasta"):
    if REF_BASE.endswith(ext):
        REF_BASE = REF_BASE[: -len(ext)]

RENAMED_REF_FASTA = RESULTS_DIR / "refs" / "reference_gvcf.fa"
COMBINED_DIR = RESULTS_DIR / "combined"
COMBINED_RAW_DIR = RESULTS_DIR / "combined_raw"


def _normalize_contig(name: str) -> str:
    name = name.strip().lower()
    name = re.sub(r"^chr", "", name, flags=re.IGNORECASE)
    m = re.match(r"^(.*?)(\d+)$", name)
    if m:
        prefix, num = m.groups()
        num = num.lstrip("0") or "0"
        name = f"{prefix}{num}"
    else:
        name = name.lstrip("0")
    return name if name else "0"


def _read_maf_contigs() -> set[str]:
    contigs = set()
    maf_files = list(MAF_DIR.glob("*.maf")) + list(MAF_DIR.glob("*.maf.gz"))
    for maf in sorted(maf_files):
        try:
            if maf.name.endswith(".gz"):
                handle = gzip.open(maf, "rt", encoding="utf-8")
            else:
                handle = maf.open("r", encoding="utf-8")
            with handle:
                for line in handle:
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split()
                    if parts and parts[0] == "s" and len(parts) >= 2:
                        contigs.add(parts[1])
        except OSError:
            continue
    return contigs


def _open_fasta(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _read_fasta_contigs(path: Path) -> list[str]:
    contigs = []
    try:
        with _open_fasta(path) as handle:
            for line in handle:
                if line.startswith(">"):
                    contigs.append(line[1:].strip().split()[0])
    except OSError:
        pass
    return contigs


def _read_gvcf_contigs(path: Path) -> list[str]:
    contigs = []
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("##contig=<ID="):
                    entry = line.strip().split("ID=", 1)[1]
                    contig = entry.split(",", 1)[0].split(">", 1)[0]
                    contigs.append(contig)
                elif line.startswith("#CHROM"):
                    break
    except OSError:
        pass
    return contigs


def _contig_map_from_gvcf(gvcf_path: Path):
    gvcf_contigs = _read_gvcf_contigs(gvcf_path)
    fasta_contigs = _read_fasta_contigs(ORIG_REF_FASTA)
    if not gvcf_contigs or not fasta_contigs:
        raise ValueError(
            "Unable to read contigs from gVCF or reference; cannot rename reference."
        )

    gvcf_norm = {}
    for name in gvcf_contigs:
        gvcf_norm.setdefault(_normalize_contig(name), []).append(name)

    fasta_norm = {}
    for name in fasta_contigs:
        fasta_norm.setdefault(_normalize_contig(name), []).append(name)

    overlap = sorted(set(gvcf_norm.keys()) & set(fasta_norm.keys()))
    if not overlap:
        raise ValueError(
            "No overlapping contigs between gVCF and reference after normalization."
        )

    mapping = {}
    for key in overlap:
        if len(gvcf_norm[key]) != 1 or len(fasta_norm[key]) != 1:
            raise ValueError(
                "Ambiguous contig mapping between gVCF and reference; aborting."
            )
        mapping[fasta_norm[key][0]] = gvcf_norm[key][0]

    return mapping


REF_FASTA = ORIG_REF_FASTA
REF_FASTA_GATK = RENAMED_REF_FASTA
REF_FAI = str(REF_FASTA_GATK) + ".fai"
REF_DICT = str(REF_FASTA_GATK.with_suffix(".dict"))

def _discover_samples():
    if "samples" in config:
        return list(config["samples"])
    maf_pattern = str(MAF_DIR / "{sample}.maf")
    maf_gz_pattern = str(MAF_DIR / "{sample}.maf.gz")
    samples = set(glob_wildcards(maf_pattern).sample)
    samples.update(glob_wildcards(maf_gz_pattern).sample)
    return sorted(samples)


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


def _maf_input(sample):
    maf = MAF_DIR / f"{sample}.maf"
    maf_gz = MAF_DIR / f"{sample}.maf.gz"
    if maf.exists():
        return str(maf)
    if maf_gz.exists():
        return str(maf_gz)
    return str(maf)


def _split_out(base, contig):
    return GVCF_DIR / "cleangVCF" / "split_gvcf" / f"{base}.{contig}.gvcf.gz"


def _combined_out(contig):
    return COMBINED_DIR / f"combined.{contig}.gvcf.gz"


def _combined_raw_out(contig):
    return COMBINED_RAW_DIR / f"combined.{contig}.gvcf.gz"


def _split_prefix(contig):
    return RESULTS_DIR / "split" / f"combined.{contig}"


def _accessibility_out(contig):
    return RESULTS_DIR / "split" / f"combined.{contig}.accessible.npz"


SPLIT_SUFFIX = ".gz" if BGZIP_OUTPUT else ""

rule all:
    # Final targets: merged gVCFs plus filtered bed masks per contig.
    input:
        [str(_combined_out(c)) for c in CONTIGS],
        [str(_split_prefix(c)) + ".filtered.bed" for c in CONTIGS],
        [str(_split_prefix(c)) + ".coverage.txt" for c in CONTIGS],
        [str(_accessibility_out(c)) for c in CONTIGS],
        str(RESULTS_DIR / "summary.html"),

rule rename_reference:
    # Create a renamed reference FASTA to match gVCF contig names.
    input:
        ref=str(ORIG_REF_FASTA),
        gvcf=lambda wc: str(_gvcf_out(GVCF_BASES[0])),
    output:
        ref=str(RENAMED_REF_FASTA),
    run:
        from pathlib import Path

        out_path = Path(output.ref)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        mapping = _contig_map_from_gvcf(Path(input.gvcf))
        with _open_fasta(Path(input.ref)) as fin, open(
            output.ref, "w", encoding="utf-8"
        ) as fout:
            for line in fin:
                if line.startswith(">"):
                    name = line[1:].strip().split()[0]
                    new_name = mapping.get(name, name)
                    fout.write(f">{new_name}\n")
                else:
                    fout.write(line)

rule index_reference:
    # Create reference FASTA index and sequence dictionary for GATK.
    input:
        ref=str(REF_FASTA_GATK),
    output:
        fai=REF_FAI,
        dict=REF_DICT,
    shell:
        """
        set -euo pipefail
        samtools faidx "{input.ref}"
        picard CreateSequenceDictionary R="{input.ref}" O="{output.dict}"
        """


def _summary_jobs() -> list[tuple[str, list[str]]]:
    jobs = [
        ("index_reference", [REF_FAI, REF_DICT]),
        ("rename_reference", [str(RENAMED_REF_FASTA)]),
        ("maf_to_gvcf", [str(_gvcf_out(base)) for base in GVCF_BASES]),
        (
            "drop_sv",
            [str(GVCF_DIR / "cleangVCF" / f"{base}.gvcf.gz") for base in GVCF_BASES]
            + [str(GVCF_DIR / "cleangVCF" / "dropped_indels.bed")],
        ),
        (
            "split_gvcf_by_contig",
            [str(_split_out(base, contig)) for contig in CONTIGS for base in GVCF_BASES],
        ),
    ]
    if VT_NORMALIZE:
        jobs.append(("merge_contig_raw", [str(_combined_raw_out(c)) for c in CONTIGS]))
        jobs.append(("normalize_merged_gvcf", [str(_combined_out(c)) for c in CONTIGS]))
    else:
        jobs.append(("merge_contig", [str(_combined_out(c)) for c in CONTIGS]))
    jobs.extend(
        [
            (
                "split_gvcf",
                [
                    str(_split_prefix(c)) + suffix
                    for c in CONTIGS
                    for suffix in (".inv", ".filtered", ".clean", ".missing.bed")
                ],
            ),
            ("check_split_coverage", [str(_split_prefix(c)) + ".coverage.txt" for c in CONTIGS]),
            ("mask_bed", [str(_split_prefix(c)) + ".filtered.bed" for c in CONTIGS]),
            ("make_accessibility", [str(_accessibility_out(c)) for c in CONTIGS]),
        ]
    )
    return jobs


def _summary_temp_paths() -> set[str]:
    temp_paths = set()
    temp_paths.update(str(_gvcf_out(base)) for base in GVCF_BASES)
    temp_paths.update(str(GVCF_DIR / "cleangVCF" / f"{base}.gvcf.gz") for base in GVCF_BASES)
    temp_paths.update(str(_split_out(base, contig)) for contig in CONTIGS for base in GVCF_BASES)
    temp_paths.update(
        str(_split_out(base, contig)) + ".tbi"
        for contig in CONTIGS
        for base in GVCF_BASES
    )
    temp_paths.update(str(RESULTS_DIR / "genomicsdb" / f"{contig}") for contig in CONTIGS)
    if VT_NORMALIZE:
        temp_paths.update(str(_combined_raw_out(c)) for c in CONTIGS)
    return temp_paths


def _summary_arg_outputs() -> list[str]:
    return (
        [str(_split_prefix(c)) + ".clean" for c in CONTIGS]
        + [str(_split_prefix(c)) + ".filtered.bed" for c in CONTIGS]
        + [str(_accessibility_out(c)) for c in CONTIGS]
    )


rule summary_report:
    # Write an HTML summary of jobs, outputs, and warnings.
    input:
        combined=[str(_combined_out(c)) for c in CONTIGS],
        beds=[str(_split_prefix(c)) + ".filtered.bed" for c in CONTIGS],
        invs=[str(_split_prefix(c)) + ".inv" + SPLIT_SUFFIX for c in CONTIGS],
        filts=[str(_split_prefix(c)) + ".filtered" + SPLIT_SUFFIX for c in CONTIGS],
        cleans=[str(_split_prefix(c)) + ".clean" + SPLIT_SUFFIX for c in CONTIGS],
        dropped=str(GVCF_DIR / "cleangVCF" / "dropped_indels.bed"),
    output:
        report=str(RESULTS_DIR / "summary.html"),
    params:
        ref_fai=str(REF_FAI),
        maf_dir=str(MAF_DIR),
        orig_ref_fasta=str(ORIG_REF_FASTA),
        contigs=[str(c) for c in CONTIGS],
        jobs=_summary_jobs(),
        temp_paths=sorted(_summary_temp_paths()),
        arg_outputs=_summary_arg_outputs(),
        split_prefixes={str(c): str(_split_prefix(c)) for c in CONTIGS},
    script:
        "scripts/summary_report.py"

rule maf_to_gvcf:
    # Convert each MAF to a gzipped gVCF using TASSEL.
    threads: MAF_TO_GVCF_THREADS
    resources:
        mem_mb=MAF_TO_GVCF_MEM_MB,
        time=MAF_TO_GVCF_TIME
    input:
        maf=lambda wc: _maf_input(wc.sample),
        ref=str(REF_FASTA),
    output:
        gvcf=temp(str(GVCF_DIR / (f"{{sample}}To{REF_BASE}.gvcf.gz"))),
        tbi=temp(str(GVCF_DIR / (f"{{sample}}To{REF_BASE}.gvcf.gz.tbi"))),
    log:
        str(Path("logs") / "tassel" / "{sample}.log"),
    params:
        tassel_dir=str(TASSEL_DIR),
        sample_name=lambda wc: f"{wc.sample}{SAMPLE_SUFFIX}",
        fill_gaps=FILL_GAPS,
        output_just_gt=OUTPUT_JUST_GT,
    shell:
        """
        set -euo pipefail
        mkdir -p "{GVCF_DIR}"
        mkdir -p "$(dirname "{log}")"
        out_base="{output.gvcf}"
        if [[ "$out_base" == *.gz ]]; then
          out_base="${{out_base%.gz}}"
        fi
        "{params.tassel_dir}/run_pipeline.pl" -Xmx{MAF_TO_GVCF_JAVA_MEM_MB}m -debug \
          -MAFToGVCFPlugin \
          -referenceFasta "{input.ref}" \
          -mafFile "{input.maf}" \
          -sampleName "{params.sample_name}" \
          -gvcfOutput "$out_base" \
          -fillGaps "{params.fill_gaps}" \
          -outputJustGT "{params.output_just_gt}" \
          > "{log}" 2>&1
        if [ -f "$out_base" ]; then
          bgzip -f -c "$out_base" > "{output.gvcf}"
          rm -f "$out_base"
        elif [ -f "{output.gvcf}" ]; then
          if ! tabix -p vcf "{output.gvcf}" >/dev/null 2>&1; then
            gunzip -c "{output.gvcf}" | bgzip -c > "{output.gvcf}.tmp"
            mv "{output.gvcf}.tmp" "{output.gvcf}"
          fi
        else
          echo "ERROR: TASSEL did not write $out_base or {output.gvcf}" >&2
          tail -n 200 "{log}" >&2 || true
          exit 1
        fi
        tabix -f -p vcf "{output.gvcf}"
        """


rule drop_sv:
    # Remove large indels from all gVCFs in the directory.
    input:
        gvcfs=[str(_gvcf_out(b)) for b in GVCF_BASES],
    output:
        bed=str(GVCF_DIR / "cleangVCF" / "dropped_indels.bed"),
        gvcfs=[temp(str(GVCF_DIR / "cleangVCF" / f"{b}.gvcf.gz")) for b in GVCF_BASES],
        tbis=[str(GVCF_DIR / "cleangVCF" / f"{b}.gvcf.gz.tbi") for b in GVCF_BASES],
        clean_dir=directory(str(GVCF_DIR / "cleangVCF")),
    params:
        cutoff=DROP_CUTOFF,
    shell:
        """
        set -euo pipefail
        if [ -n "{params.cutoff}" ]; then
          python "{workflow.basedir}/scripts/dropSV.py" -d "{GVCF_DIR}" -c "{params.cutoff}"
        else
          python "{workflow.basedir}/scripts/dropSV.py" -d "{GVCF_DIR}"
        fi
        """


rule split_gvcf_by_contig:
    # Split each cleaned gVCF into per-contig gVCFs for merge.
    input:
        gvcf=lambda wc: str(GVCF_DIR / "cleangVCF" / f"{wc.gvcf_base}.gvcf.gz"),
        ref=str(REF_FASTA_GATK),
        fai=REF_FAI,
        dict=REF_DICT,
        bed=str(GVCF_DIR / "cleangVCF" / "dropped_indels.bed"),
    output:
        gvcf=temp(str(GVCF_DIR / "cleangVCF" / "split_gvcf" / "{gvcf_base}.{contig}.gvcf.gz")),
        tbi=temp(str(GVCF_DIR / "cleangVCF" / "split_gvcf" / "{gvcf_base}.{contig}.gvcf.gz.tbi")),
    shell:
        """
        set -euo pipefail
        mkdir -p "{GVCF_DIR}/cleangVCF/split_gvcf"
        tmp_vcf="{output.gvcf}.tmp.vcf"
        gatk --java-options "-Xmx{DEFAULT_JAVA_MEM_MB}m -Xms{DEFAULT_JAVA_MEM_MB}m" SelectVariants \
          -R "{input.ref}" \
          -V "{input.gvcf}" \
          -L "{wildcards.contig}" \
          -O "$tmp_vcf"
        bgzip -f -c "$tmp_vcf" > "{output.gvcf}"
        rm -f "$tmp_vcf"
        tabix -f -p vcf "{output.gvcf}"
        """


if VT_NORMALIZE:
    rule merge_contig_raw:
        # Merge all samples for a contig with GenomicsDBImport + GenotypeGVCFs.
        threads: MERGE_CONTIG_THREADS
        resources:
            mem_mb=MERGE_CONTIG_MEM_MB,
            time=MERGE_CONTIG_TIME,
            merge_gvcf_jobs=1,
        input:
            gvcfs=lambda wc: [str(_split_out(b, wc.contig)) for b in GVCF_BASES],
            tbis=lambda wc: [str(_split_out(b, wc.contig)) + ".tbi" for b in GVCF_BASES],
            ref=str(REF_FASTA_GATK),
            fai=REF_FAI,
            dict=REF_DICT,
            bed=str(GVCF_DIR / "cleangVCF" / "dropped_indels.bed"),
        output:
            gvcf=str(COMBINED_RAW_DIR / "combined.{contig}.gvcf.gz"),
            workspace=temp(directory(str(RESULTS_DIR / "genomicsdb" / "{contig}"))),
        params:
            gvcf_args=lambda wc: " ".join(
                f"-V {str(_split_out(b, wc.contig))}" for b in GVCF_BASES
            ),
            vcf_buffer_size=GENOMICSDB_VCF_BUFFER_SIZE,
            segment_size=GENOMICSDB_SEGMENT_SIZE,
        shell:
            """
            set -euo pipefail
            mkdir -p "{COMBINED_RAW_DIR}"
            gatk --java-options "-Xmx{MERGE_CONTIG_JAVA_MEM_MB}m -Xms{MERGE_CONTIG_JAVA_MEM_MB}m" GenomicsDBImport \
              {params.gvcf_args} \
              --genomicsdb-workspace-path "{output.workspace}" \
              -L "{wildcards.contig}" \
              --genomicsdb-vcf-buffer-size {params.vcf_buffer_size} \
              --genomicsdb-segment-size {params.segment_size}
            gatk --java-options "-Xmx{MERGE_CONTIG_JAVA_MEM_MB}m -Xms{MERGE_CONTIG_JAVA_MEM_MB}m" SelectVariants \
              -R "{input.ref}" \
              -V "gendb://{output.workspace}" \
              -O "{output.gvcf}" \
              -L "{wildcards.contig}" \
              --call-genotypes
            """
else:
    rule merge_contig:
        # Merge all samples for a contig with GenomicsDBImport + GenotypeGVCFs.
        threads: MERGE_CONTIG_THREADS
        resources:
            mem_mb=MERGE_CONTIG_MEM_MB,
            time=MERGE_CONTIG_TIME,
            merge_gvcf_jobs=1,
        input:
            gvcfs=lambda wc: [str(_split_out(b, wc.contig)) for b in GVCF_BASES],
            tbis=lambda wc: [str(_split_out(b, wc.contig)) + ".tbi" for b in GVCF_BASES],
            ref=str(REF_FASTA_GATK),
            fai=REF_FAI,
            dict=REF_DICT,
            bed=str(GVCF_DIR / "cleangVCF" / "dropped_indels.bed"),
        output:
            gvcf=str(COMBINED_DIR / "combined.{contig}.gvcf.gz"),
            workspace=temp(directory(str(RESULTS_DIR / "genomicsdb" / "{contig}"))),
        params:
            gvcf_args=lambda wc: " ".join(
                f"-V {str(_split_out(b, wc.contig))}" for b in GVCF_BASES
            ),
            vcf_buffer_size=GENOMICSDB_VCF_BUFFER_SIZE,
            segment_size=GENOMICSDB_SEGMENT_SIZE,
        shell:
            """
            set -euo pipefail
            mkdir -p "{COMBINED_DIR}"
            gatk --java-options "-Xmx{MERGE_CONTIG_JAVA_MEM_MB}m -Xms{MERGE_CONTIG_JAVA_MEM_MB}m" GenomicsDBImport \
              {params.gvcf_args} \
              --genomicsdb-workspace-path "{output.workspace}" \
              -L "{wildcards.contig}" \
              --genomicsdb-vcf-buffer-size {params.vcf_buffer_size} \
              --genomicsdb-segment-size {params.segment_size}
            gatk --java-options "-Xmx{MERGE_CONTIG_JAVA_MEM_MB}m -Xms{MERGE_CONTIG_JAVA_MEM_MB}m" SelectVariants \
              -R "{input.ref}" \
              -V "gendb://{output.workspace}" \
              -O "{output.gvcf}" \
              -L "{wildcards.contig}" \
              --call-genotypes
            """


if VT_NORMALIZE:
    rule normalize_merged_gvcf:
        # Normalize merged gVCFs with vt after GenotypeGVCFs.
        input:
            gvcf=str(COMBINED_RAW_DIR / "combined.{contig}.gvcf.gz"),
            ref=str(REF_FASTA_GATK),
        output:
            gvcf=str(COMBINED_DIR / "combined.{contig}.gvcf.gz"),
        shell:
            """
            set -euo pipefail
            mkdir -p "{COMBINED_DIR}"
            tmp_vcf="{output.gvcf}.tmp.vcf"
            "{VT_PATH}" normalize "{input.gvcf}" -r "{input.ref}" -o "$tmp_vcf"
            bgzip -f -c "$tmp_vcf" > "{output.gvcf}"
            rm -f "$tmp_vcf"
            tabix -f -p vcf "{output.gvcf}"
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
        filter_multiallelic=FILTER_MULTIALLELIC,
        bgzip_output=BGZIP_OUTPUT,
        out_prefix=lambda wc: str(_split_prefix(wc.contig)),
    shell:
        """
        set -euo pipefail
        mkdir -p "{RESULTS_DIR}/split"
        cmd=(python "{workflow.basedir}/scripts/split.py" --out-prefix "{params.out_prefix}" --fai "{input.ref_fai}")
        if [ "{params.filter_multiallelic}" = "True" ]; then
          cmd+=(--filter-multiallelic)
        fi
        if [ "{params.bgzip_output}" = "True" ]; then
          cmd+=(--bgzip-output)
        fi
        cmd+=("{input.gvcf}")
        "${{cmd[@]}}"
        """


rule check_split_coverage:
    # Validate that clean + inv + filtered bed sum to contig length.
    input:
        clean=lambda wc: str(_split_prefix(wc.contig)) + ".clean" + SPLIT_SUFFIX,
        inv=lambda wc: str(_split_prefix(wc.contig)) + ".inv" + SPLIT_SUFFIX,
        bed=lambda wc: str(_split_prefix(wc.contig)) + ".filtered.bed",
        fai=REF_FAI,
    output:
        report=str(RESULTS_DIR / "split" / "combined.{contig}.coverage.txt"),
    params:
        prefix=lambda wc: str(_split_prefix(wc.contig)),
    shell:
        """
        set -euo pipefail
        python "{workflow.basedir}/scripts/check_split_coverage.py" \
          "{params.prefix}" \
          --fai "{input.fai}"
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
        prefix=lambda wc: str(_split_prefix(wc.contig)),
    shell:
        """
        set -euo pipefail
        cmd=(python "{workflow.basedir}/scripts/filt_to_bed.py" "{params.prefix}")
        cmd+=(--dropped-bed "{input.dropped}")
        "${{cmd[@]}}"
        """


rule make_accessibility:
    # Build boolean accessibility array from clean + inv VCFs per contig.
    input:
        clean=lambda wc: str(_split_prefix(wc.contig)) + ".clean" + SPLIT_SUFFIX,
        inv=lambda wc: str(_split_prefix(wc.contig)) + ".inv" + SPLIT_SUFFIX,
        fai=REF_FAI,
    output:
        mask=str(RESULTS_DIR / "split" / "combined.{contig}.accessible.npz"),
    params:
        contig="{contig}",
    shell:
        """
        set -euo pipefail
        python "{workflow.basedir}/scripts/build_accessibility.py" \
          --clean "{input.clean}" \
          --inv "{input.inv}" \
          --fai "{input.fai}" \
          --contig "{params.contig}" \
          --output "{output.mask}"
        """
