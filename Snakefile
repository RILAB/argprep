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


def _split_sanitized_out(base, contig):
    return GVCF_DIR / "cleangVCF" / "split_gvcf_sanitized" / f"{base}.{contig}.gvcf.gz"


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
    run:
        from pathlib import Path
        import gzip
        import html

        report_path = Path(output.report)
        report_path.parent.mkdir(parents=True, exist_ok=True)

        def _open_text(path: str):
            if str(path).endswith(".gz"):
                return gzip.open(path, "rt")
            return open(path, "r", encoding="utf-8", errors="ignore")

        def _read_fai_lengths(path: str) -> dict[str, int]:
            lengths: dict[str, int] = {}
            with open(path, "r", encoding="utf-8") as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2:
                        continue
                    try:
                        lengths[parts[0]] = int(parts[1])
                    except ValueError:
                        continue
            return lengths

        def _window_index(pos: int, window: int) -> int:
            return max((pos - 1) // window, 0)

        def _is_variant_alt(alt_field: str) -> bool:
            if not alt_field or alt_field == ".":
                return False
            alts = [a.strip() for a in alt_field.split(",") if a.strip()]
            alts = [a for a in alts if a != "<NON_REF>"]
            return len(alts) > 0

        def _add_bed_counts(
            path: str,
            counts: dict[str, list[int]],
            contig_lengths: dict[str, int],
            window: int,
        ) -> None:
            try:
                with _open_text(path) as f_in:
                    for line in f_in:
                        if not line or line.startswith("#"):
                            continue
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 3:
                            continue
                        contig = parts[0]
                        if contig not in counts:
                            continue
                        try:
                            start = int(parts[1])
                            end = int(parts[2])
                        except ValueError:
                            continue
                        if end <= start:
                            continue
                        # BED is 0-based half-open, convert to 1-based positions for windows.
                        pos = start + 1
                        while pos <= end:
                            idx = _window_index(pos, window)
                            if idx >= len(counts[contig]):
                                break
                            window_end = min((idx + 1) * window, contig_lengths[contig])
                            span_end = min(end, window_end)
                            counts[contig][idx] += max(span_end - pos + 1, 0)
                            pos = span_end + 1
            except OSError:
                pass

        def _svg_scatter_plot(
            xs: list[int],
            ys: list[int],
            width: int = 900,
            height: int = 240,
            x_label: str | None = None,
            y_label: str | None = None,
        ) -> str:
            if not xs or not ys or len(xs) != len(ys):
                return "<p>No data available.</p>"
            margin = {"left": 70, "right": 20, "top": 30, "bottom": 50}
            plot_w = width - margin["left"] - margin["right"]
            plot_h = height - margin["top"] - margin["bottom"]
            x_min = min(xs)
            x_max = max(xs)
            y_min = 0
            y_max = max(ys)
            if x_max == x_min:
                x_max = x_min + 1
            if y_max == y_min:
                y_max = y_min + 1

            def x_scale(x_val: int) -> float:
                return margin["left"] + (x_val - x_min) / (x_max - x_min) * plot_w

            def y_scale(y_val: int) -> float:
                return margin["top"] + plot_h - (y_val - y_min) / (y_max - y_min) * plot_h

            parts = [
                f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
                'xmlns="http://www.w3.org/2000/svg" role="img">',
                '<rect width="100%" height="100%" fill="white"/>',
            ]
            if y_label:
                parts.append(
                    f'<text x="16" y="{height/2}" text-anchor="middle" '
                    f'font-size="12" font-family="sans-serif" '
                    f'transform="rotate(-90 16 {height/2})">{html.escape(y_label)}</text>'
                )
            if x_label:
                parts.append(
                    f'<text x="{width/2}" y="{height-8}" text-anchor="middle" '
                    f'font-size="12" font-family="sans-serif">{html.escape(x_label)}</text>'
                )

            x0 = margin["left"]
            y0 = margin["top"] + plot_h
            parts.append(
                f'<line x1="{x0}" y1="{y0}" x2="{x0 + plot_w}" y2="{y0}" '
                'stroke="#333" stroke-width="1"/>'
            )
            parts.append(
                f'<line x1="{x0}" y1="{margin["top"]}" x2="{x0}" y2="{y0}" '
                'stroke="#333" stroke-width="1"/>'
            )
            # Y ticks
            for i in range(5):
                frac = i / 4
                y = y0 - frac * plot_h
                val = int(round(y_min + frac * (y_max - y_min)))
                parts.append(
                    f'<line x1="{x0 - 4}" y1="{y:.2f}" x2="{x0}" y2="{y:.2f}" '
                    'stroke="#333" stroke-width="1"/>'
                )
                parts.append(
                    f'<text x="{x0 - 8}" y="{y + 4:.2f}" text-anchor="end" '
                    f'font-size="10" font-family="sans-serif">{val:,}</text>'
                )
            # X ticks
            for i in range(5):
                frac = i / 4
                x = x0 + frac * plot_w
                val = int(round(x_min + frac * (x_max - x_min)))
                parts.append(
                    f'<line x1="{x:.2f}" y1="{y0}" x2="{x:.2f}" y2="{y0 + 4}" '
                    'stroke="#333" stroke-width="1"/>'
                )
                parts.append(
                    f'<text x="{x:.2f}" y="{y0 + 16}" text-anchor="middle" '
                    f'font-size="10" font-family="sans-serif">{val:,}</text>'
                )

            for x_val, y_val in zip(xs, ys):
                parts.append(
                    f'<circle cx="{x_scale(x_val):.2f}" cy="{y_scale(y_val):.2f}" '
                    'r="2.5" fill="#4C78A8" />'
                )
            parts.append("</svg>")
            return "\n".join(parts)

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

        temp_paths = set()
        temp_paths.update(str(_gvcf_out(base)) for base in GVCF_BASES)
        temp_paths.update(
            str(GVCF_DIR / "cleangVCF" / f"{base}.gvcf.gz") for base in GVCF_BASES
        )
        temp_paths.update(
            str(_split_out(base, contig))
            for contig in CONTIGS
            for base in GVCF_BASES
        )
        temp_paths.update(
            str(_split_out(base, contig)) + ".tbi"
            for contig in CONTIGS
            for base in GVCF_BASES
        )
        temp_paths.update(str(RESULTS_DIR / "genomicsdb" / f"{contig}") for contig in CONTIGS)
        if VT_NORMALIZE:
            temp_paths.update(str(_combined_raw_out(c)) for c in CONTIGS)

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
        try:
            maf_contigs = _read_maf_contigs()
            ref_contigs = set(_read_fasta_contigs(ORIG_REF_FASTA))
            missing_in_ref = sorted(set(maf_contigs) - ref_contigs)
            missing_in_maf = sorted(ref_contigs - set(maf_contigs))
            if missing_in_ref:
                warnings.append(
                    "WARNING: MAF contigs not present in reference (showing up to 5): "
                    + ", ".join(missing_in_ref[:5])
                )
            if missing_in_maf:
                warnings.append(
                    "WARNING: Reference contigs not present in MAFs (showing up to 5): "
                    + ", ".join(missing_in_maf[:5])
                )
        except Exception as exc:
            warnings.append(f"WARNING: Failed to compare MAF vs reference contigs: {exc}")

        with report_path.open("w", encoding="utf-8") as handle:
            handle.write("<!doctype html>\n")
            handle.write("<html lang=\"en\">\n")
            handle.write("<head>\n")
            handle.write("<meta charset=\"utf-8\" />\n")
            handle.write("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />\n")
            handle.write("<title>Workflow summary</title>\n")
            handle.write("<style>\n")
            handle.write("body { font-family: sans-serif; margin: 24px; color: #111; }\n")
            handle.write("h1, h2, h3 { margin-top: 1.4em; }\n")
            handle.write("code { background: #f6f6f6; padding: 0 4px; }\n")
            handle.write("table { border-collapse: collapse; margin: 12px 0; }\n")
            handle.write("th, td { border: 1px solid #ccc; padding: 4px 8px; }\n")
            handle.write(".temp { color: #555; }\n")
            handle.write("</style>\n")
            handle.write("</head>\n")
            handle.write("<body>\n")

            handle.write("<h1>Workflow summary</h1>\n")
            handle.write("<h2>Jobs run</h2>\n")
            handle.write("<ul>\n")
            for job, outputs in jobs:
                handle.write(f"<li><strong>{html.escape(job)}</strong>\n")
                handle.write("<ul>\n")
                for path in outputs:
                    mark = " *" if path in temp_paths else ""
                    cls = " class=\"temp\"" if path in temp_paths else ""
                    handle.write(
                        f"<li{cls}><code>{html.escape(path)}</code>{html.escape(mark)}</li>\n"
                    )
                handle.write("</ul>\n")
                handle.write("</li>\n")
            handle.write("</ul>\n")
            handle.write(
                "<p><em>Temporary outputs are marked with an asterisk and are removed "
                "after a successful run.</em></p>\n"
            )

            handle.write("<h2>Files for ARG estimation</h2>\n")
            arg_outputs = (
                [str(_split_prefix(c)) + ".clean" for c in CONTIGS]
                + [str(_split_prefix(c)) + ".filtered.bed" for c in CONTIGS]
                + [str(_accessibility_out(c)) for c in CONTIGS]
            )
            handle.write("<ul>\n")
            for path in arg_outputs:
                handle.write(f"<li><code>{html.escape(path)}</code></li>\n")
            handle.write("</ul>\n")
            handle.write(
                "<p><em>Accessibility arrays are provided to enable computing statistics "
                "with scikit-allel.</em></p>\n"
            )

            handle.write("<h2>Site density (100kb windows)</h2>\n")
            contig_lengths = _read_fai_lengths(REF_FAI)
            window = 100_000
            filtered_counts: dict[str, list[int]] = {}
            inv_counts: dict[str, list[int]] = {}
            variant_counts: dict[str, list[int]] = {}
            contig_names = [str(c) for c in CONTIGS]
            for contig in contig_names:
                length = contig_lengths.get(contig)
                if length is None:
                    continue
                n_windows = (length + window - 1) // window
                filtered_counts[contig] = [0 for _ in range(n_windows)]
                inv_counts[contig] = [0 for _ in range(n_windows)]
                variant_counts[contig] = [0 for _ in range(n_windows)]

            for path in input.beds:
                try:
                    _add_bed_counts(path, filtered_counts, contig_lengths, window)
                except OSError as exc:
                    handle.write(
                        f"<p>Failed to read filtered bed from {html.escape(path)}: "
                        f"{html.escape(str(exc))}</p>\n"
                    )

            # Prefer invariant BED if present; fallback to .inv VCF counts.
            for contig in contig_names:
                if contig not in inv_counts:
                    continue
                inv_bed = str(_split_prefix(contig)) + ".inv.bed"
                if Path(inv_bed).exists():
                    _add_bed_counts(inv_bed, inv_counts, contig_lengths, window)
                elif Path(inv_bed + ".gz").exists():
                    _add_bed_counts(inv_bed + ".gz", inv_counts, contig_lengths, window)

            for path in input.invs:
                try:
                    with _open_text(path) as f_in:
                        for line in f_in:
                            if not line or line.startswith("#"):
                                continue
                            parts = line.rstrip("\n").split("\t")
                            if len(parts) < 2:
                                continue
                            contig = parts[0]
                            if contig not in inv_counts:
                                continue
                            try:
                                pos = int(parts[1])
                            except ValueError:
                                continue
                            idx = _window_index(pos, window)
                            if idx < len(inv_counts[contig]):
                                inv_counts[contig][idx] += 1
                except OSError as exc:
                    handle.write(
                        f"<p>Failed to read invariant sites from {html.escape(path)}: "
                        f"{html.escape(str(exc))}</p>\n"
                    )

            for clean_path in input.cleans:
                try:
                    with _open_text(clean_path) as f_in:
                        for line in f_in:
                            if not line or line.startswith("#"):
                                continue
                            parts = line.rstrip("\n").split("\t")
                            if len(parts) < 5:
                                continue
                            contig = parts[0]
                            if contig not in variant_counts:
                                continue
                            if not _is_variant_alt(parts[4]):
                                continue
                            try:
                                pos = int(parts[1])
                            except ValueError:
                                continue
                            idx = _window_index(pos, window)
                            if idx < len(variant_counts[contig]):
                                variant_counts[contig][idx] += 1
                except OSError as exc:
                    handle.write(
                        f"<p>Failed to read variable sites from {html.escape(clean_path)}: "
                        f"{html.escape(str(exc))}</p>\n"
                    )

            handle.write("<h2>Summary Plots</h2>\n")
            for contig in contig_names:
                if contig not in filtered_counts or contig not in inv_counts or contig not in variant_counts:
                    continue
                length = contig_lengths.get(contig, 0)
                if length <= 0:
                    continue
                midpoints = []
                for idx in range(len(filtered_counts[contig])):
                    start = idx * window
                    end = min(start + window, length)
                    midpoints.append(start + (end - start) // 2)
                handle.write(f"<h3>{html.escape(contig)}</h3>\n")
                handle.write(f"<h4>Dropped + filtered bp for contig {html.escape(contig)}</h4>\n")
                handle.write(
                    _svg_scatter_plot(
                        midpoints,
                        filtered_counts[contig],
                        x_label="Window midpoint (bp)",
                        y_label="Base pairs",
                    )
                )
                handle.write(f"<h4>Invariant sites for contig {html.escape(contig)}</h4>\n")
                handle.write(
                    _svg_scatter_plot(
                        midpoints,
                        inv_counts[contig],
                        x_label="Window midpoint (bp)",
                        y_label="Site count",
                    )
                )
                handle.write(f"<h4>Variant sites for contig {html.escape(contig)}</h4>\n")
                handle.write(
                    _svg_scatter_plot(
                        midpoints,
                        variant_counts[contig],
                        x_label="Window midpoint (bp)",
                        y_label="Site count",
                    )
                )

            handle.write("<h2>Warnings</h2>\n")
            if warnings:
                handle.write("<ul>\n")
                for line in warnings:
                    handle.write(f"<li>{html.escape(line)}</li>\n")
                handle.write("</ul>\n")
            else:
                handle.write("<p>None found in logs</p>\n")

            handle.write("</body>\n</html>\n")


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


rule sanitize_split_gvcf:
    # Drop duplicate alleles/records before GenomicsDBImport.
    input:
        gvcf=lambda wc: str(_split_out(wc.gvcf_base, wc.contig)),
        tbi=lambda wc: str(_split_out(wc.gvcf_base, wc.contig)) + ".tbi",
    log:
        str(Path("logs") / "sanitize" / "{gvcf_base}.{contig}.log"),
    output:
        gvcf=temp(str(GVCF_DIR / "cleangVCF" / "split_gvcf_sanitized" / "{gvcf_base}.{contig}.gvcf.gz")),
        tbi=temp(str(GVCF_DIR / "cleangVCF" / "split_gvcf_sanitized" / "{gvcf_base}.{contig}.gvcf.gz.tbi")),
    shell:
        """
        set -euo pipefail
        mkdir -p "{GVCF_DIR}/cleangVCF/split_gvcf_sanitized"
        mkdir -p "$(dirname "{log}")"
        before=$(bcftools view -H "{input.gvcf}" | wc -l | tr -d ' ')
        bcftools norm -d both -Oz -o "{output.gvcf}" "{input.gvcf}"
        tabix -f -p vcf "{output.gvcf}"
        after=$(bcftools view -H "{output.gvcf}" | wc -l | tr -d ' ')
        removed=$((before - after))
        printf "input=%s\noutput=%s\nrecords_before=%s\nrecords_after=%s\nrecords_removed=%s\n" \
          "{input.gvcf}" "{output.gvcf}" "$before" "$after" "$removed" > "{log}"
        """


if VT_NORMALIZE:
    rule merge_contig_raw:
        # Merge all samples for a contig with GenomicsDBImport + GenotypeGVCFs.
        threads: MERGE_CONTIG_THREADS
        resources:
            mem_mb=MERGE_CONTIG_MEM_MB,
            time=MERGE_CONTIG_TIME,
        input:
            gvcfs=lambda wc: [str(_split_sanitized_out(b, wc.contig)) for b in GVCF_BASES],
            tbis=lambda wc: [str(_split_sanitized_out(b, wc.contig)) + ".tbi" for b in GVCF_BASES],
            ref=str(REF_FASTA_GATK),
            fai=REF_FAI,
            dict=REF_DICT,
            bed=str(GVCF_DIR / "cleangVCF" / "dropped_indels.bed"),
        output:
            gvcf=str(COMBINED_RAW_DIR / "combined.{contig}.gvcf.gz"),
            workspace=temp(directory(str(RESULTS_DIR / "genomicsdb" / "{contig}"))),
        params:
            gvcf_args=lambda wc: " ".join(
                f"-V {str(_split_sanitized_out(b, wc.contig))}" for b in GVCF_BASES
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
        input:
            gvcfs=lambda wc: [str(_split_sanitized_out(b, wc.contig)) for b in GVCF_BASES],
            tbis=lambda wc: [str(_split_sanitized_out(b, wc.contig)) + ".tbi" for b in GVCF_BASES],
            ref=str(REF_FASTA_GATK),
            fai=REF_FAI,
            dict=REF_DICT,
            bed=str(GVCF_DIR / "cleangVCF" / "dropped_indels.bed"),
        output:
            gvcf=str(COMBINED_DIR / "combined.{contig}.gvcf.gz"),
            workspace=temp(directory(str(RESULTS_DIR / "genomicsdb" / "{contig}"))),
        params:
            gvcf_args=lambda wc: " ".join(
                f"-V {str(_split_sanitized_out(b, wc.contig))}" for b in GVCF_BASES
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
