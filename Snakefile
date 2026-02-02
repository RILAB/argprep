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
DROP_CUTOFF = config.get("drop_cutoff", "")
FILTER_MULTIALLELIC = bool(config.get("filter_multiallelic", False))
GZIP_OUTPUT = bool(config.get("gzip_output", False))
NO_MERGE = bool(config.get("no_merge", False))

REF_BASE = ORIG_REF_FASTA.name
if REF_BASE.endswith(".gz"):
    REF_BASE = REF_BASE[: -len(".gz")]
for ext in (".fa", ".fasta"):
    if REF_BASE.endswith(ext):
        REF_BASE = REF_BASE[: -len(ext)]

RENAMED_REF_FASTA = RESULTS_DIR / "refs" / "reference_renamed.fa"


def _normalize_contig(name: str) -> str:
    name = name.strip()
    name = re.sub(r"^chr", "", name, flags=re.IGNORECASE)
    name = name.lstrip("0")
    return name if name else "0"


def _read_maf_contigs() -> set[str]:
    contigs = set()
    for maf in sorted(MAF_DIR.glob("*.maf")):
        try:
            with maf.open("r", encoding="utf-8") as handle:
                for line in handle:
                    if not line or line.startswith("#"):
                        continue
                    if line.startswith("s "):
                        parts = line.split()
                        if len(parts) >= 2:
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


def _contig_map():
    maf_contigs = _read_maf_contigs()
    fasta_contigs = _read_fasta_contigs(ORIG_REF_FASTA)
    if not maf_contigs or not fasta_contigs:
        return False, {}

    maf_norm = {}
    for name in maf_contigs:
        maf_norm.setdefault(_normalize_contig(name), []).append(name)

    fasta_norm = {}
    for name in fasta_contigs:
        fasta_norm.setdefault(_normalize_contig(name), []).append(name)

    if set(maf_norm.keys()) != set(fasta_norm.keys()):
        print(
            "WARNING: MAF and FASTA contig sets differ. "
            "Only chr/Chr prefixes and leading zeros are auto-resolved.",
            file=os.sys.stderr,
        )
        return False, {}

    mapping = {}
    for key in maf_norm:
        if len(maf_norm[key]) != 1 or len(fasta_norm[key]) != 1:
            print(
                "WARNING: MAF/FASTA contig mapping is ambiguous; "
                "not renaming reference.",
                file=os.sys.stderr,
            )
            return False, {}
        mapping[fasta_norm[key][0]] = maf_norm[key][0]

    if any(k != v for k, v in mapping.items()):
        print(
            "WARNING: MAF/FASTA contigs differ only by chr prefix/zeros; "
            "a renamed reference will be generated to match MAF contigs.",
            file=os.sys.stderr,
        )
        return True, mapping

    return False, {}


NEED_RENAME, CONTIG_MAP = _contig_map()
REF_FASTA = RENAMED_REF_FASTA if NEED_RENAME else ORIG_REF_FASTA
REF_FAI = str(REF_FASTA) + ".fai"
REF_DICT = str(REF_FASTA.with_suffix(".dict"))

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


def _default_depth():
    depth_cfg = config.get("depth", None)
    if depth_cfg is None or str(depth_cfg).strip() == "":
        return max(1, len(SAMPLES))
    return int(depth_cfg)


DEPTH = _default_depth()


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
    return RESULTS_DIR / "combined" / f"combined.{contig}.gvcf.gz"


def _split_prefix(contig):
    return RESULTS_DIR / "split" / f"combined.{contig}"


SPLIT_SUFFIX = ".gz" if GZIP_OUTPUT else ""

rule all:
    # Final targets: merged gVCFs plus filtered bed masks per contig.
    input:
        [str(_combined_out(c)) for c in CONTIGS],
        [str(_split_prefix(c)) + ".filtered.bed" for c in CONTIGS],
        str(RESULTS_DIR / "summary.md"),

if NEED_RENAME:
    rule rename_reference:
        # Create a renamed reference FASTA to match MAF contig names.
        input:
            ref=str(ORIG_REF_FASTA),
        output:
            ref=str(RENAMED_REF_FASTA),
        run:
            from pathlib import Path

            out_path = Path(output.ref)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            with _open_fasta(Path(input.ref)) as fin, open(
                output.ref, "w", encoding="utf-8"
            ) as fout:
                for line in fin:
                    if line.startswith(">"):
                        name = line[1:].strip().split()[0]
                        new_name = CONTIG_MAP.get(name, name)
                        fout.write(f">{new_name}\n")
                    else:
                        fout.write(line)

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

rule summary_report:
    # Write a markdown summary of jobs, outputs, and warnings.
    input:
        combined=[str(_combined_out(c)) for c in CONTIGS],
        beds=[str(_split_prefix(c)) + ".filtered.bed" for c in CONTIGS],
        dropped=str(GVCF_DIR / "cleangVCF" / "dropped_indels.bed"),
    output:
        report=str(RESULTS_DIR / "summary.md"),
    run:
        from pathlib import Path

        report_path = Path(output.report)
        report_path.parent.mkdir(parents=True, exist_ok=True)

        jobs = [
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
            ("merge_contig", [str(_combined_out(c)) for c in CONTIGS]),
            (
                "split_gvcf",
                [
                    str(_split_prefix(c)) + suffix
                    for c in CONTIGS
                    for suffix in (".inv", ".filtered", ".clean", ".missing.bed")
                ],
            ),
            ("mask_bed", [str(_split_prefix(c)) + ".filtered.bed" for c in CONTIGS]),
        ]
        if NEED_RENAME:
            jobs.insert(0, ("rename_reference", [str(RENAMED_REF_FASTA)]))
        jobs.insert(0, ("index_reference", [REF_FAI, REF_DICT]))

        # Collect warnings from logs and snakemake logs.
        warnings = []
        log_paths = []
        log_paths.extend(sorted(Path("logs").rglob("*.log")))
        log_paths.extend(sorted(Path("logs").rglob("*.out")))
        log_paths.extend(sorted(Path("logs").rglob("*.err")))
        log_paths.extend(sorted(Path(".snakemake").rglob("*.log")))
        for log_path in log_paths:
            try:
                with log_path.open("r", encoding="utf-8", errors="ignore") as handle:
                    for line in handle:
                        if "WARNING" in line or "Warning" in line:
                            warnings.append(f"{log_path}: {line.rstrip()}")
            except OSError:
                continue

        with report_path.open("w", encoding="utf-8") as handle:
            handle.write("# Workflow summary\n\n")
            handle.write("## Jobs run\n")
            for job, outputs in jobs:
                handle.write(f"- {job}\n")
                for path in outputs:
                    handle.write(f"  - {path}\n")
            handle.write("\n## Files for ARG estimation\n")
            arg_outputs = (
                [str(_split_prefix(c)) + ".clean" for c in CONTIGS]
                + [str(_split_prefix(c)) + ".filtered.bed" for c in CONTIGS]
            )
            for path in arg_outputs:
                handle.write(f"- {path}\n")
            handle.write("\n## Warnings\n")
            if warnings:
                for line in warnings:
                    handle.write(f"- {line}\n")
            else:
                handle.write("- None found in logs\n")


rule maf_to_gvcf:
    # Convert each MAF to a gzipped gVCF using TASSEL.
    threads: 2
    resources:
        mem_mb=256000,
        time="24:00:00"
    input:
        maf=lambda wc: _maf_input(wc.sample),
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
          python3 "{workflow.basedir}/scripts/dropSV.py" -d "{GVCF_DIR}" -c "{params.cutoff}"
        else
          python3 "{workflow.basedir}/scripts/dropSV.py" -d "{GVCF_DIR}"
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
        cmd=(python3 "{workflow.basedir}/scripts/split.py" --depth="{params.depth}" --out-prefix "{params.out_prefix}" --fai "{input.ref_fai}")
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
        cmd=(python3 "{workflow.basedir}/scripts/filt_to_bed.py" "{params.prefix}")
        cmd+=(--dropped-bed "{input.dropped}")
        if [ "{params.no_merge}" = "True" ]; then
          cmd+=(--no-merge)
        fi
        "${{cmd[@]}}"
        """
