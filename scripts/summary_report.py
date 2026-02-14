#!/usr/bin/env python3
from pathlib import Path
import gzip
import html


def _open_text(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8", errors="ignore")


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


def _read_maf_contigs(maf_dir: Path) -> set[str]:
    contigs = set()
    maf_files = list(maf_dir.glob("*.maf")) + list(maf_dir.glob("*.maf.gz"))
    for maf in sorted(maf_files):
        try:
            if maf.name.endswith(".gz"):
                handle = gzip.open(maf, "rt", encoding="utf-8")
            else:
                handle = maf.open("r", encoding="utf-8")
            with handle:
                first_src = None
                for line in handle:
                    if not line or line.startswith("#"):
                        continue
                    stripped = line.strip()
                    if not stripped:
                        if first_src is not None:
                            contigs.add(first_src)
                            first_src = None
                        continue
                    parts = stripped.split()
                    if not parts:
                        continue
                    if parts[0] == "a":
                        if first_src is not None:
                            contigs.add(first_src)
                        first_src = None
                        continue
                    if parts[0] == "s" and len(parts) >= 2 and first_src is None:
                        # Use the first sequence source per alignment block (reference-side in our MAFs).
                        first_src = parts[1]
                if first_src is not None:
                    contigs.add(first_src)
        except OSError:
            continue
    return contigs


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


report_path = Path(snakemake.output.report)
report_path.parent.mkdir(parents=True, exist_ok=True)

contigs = [str(c) for c in snakemake.params.contigs]
jobs = [(str(job), [str(p) for p in outputs]) for job, outputs in snakemake.params.jobs]
temp_paths = set(str(p) for p in snakemake.params.temp_paths)
arg_outputs = [str(p) for p in snakemake.params.arg_outputs]
split_prefixes = {str(k): str(v) for k, v in dict(snakemake.params.split_prefixes).items()}

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
    maf_contigs = _read_maf_contigs(Path(snakemake.params.maf_dir))
    ref_contigs = set(_read_fasta_contigs(Path(snakemake.params.orig_ref_fasta)))
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
    handle.write('<html lang="en">\n')
    handle.write("<head>\n")
    handle.write('<meta charset="utf-8" />\n')
    handle.write('<meta name="viewport" content="width=device-width, initial-scale=1" />\n')
    handle.write("<title>Workflow summary</title>\n")
    handle.write("<style>\n")
    handle.write("body { font-family: sans-serif; margin: 24px; color: #111; }\n")
    handle.write("h1, h2, h3 { margin-top: 1.4em; }\n")
    handle.write("code { background: #f6f6f6; padding: 0 4px; }\n")
    handle.write("table { border-collapse: collapse; margin: 12px 0; }\n")
    handle.write("th, td { border: 1px solid #ccc; padding: 4px 8px; }\n")
    handle.write('.temp { color: #555; }\n')
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
            cls = ' class="temp"' if path in temp_paths else ""
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
    handle.write("<ul>\n")
    for path in arg_outputs:
        handle.write(f"<li><code>{html.escape(path)}</code></li>\n")
    handle.write("</ul>\n")
    handle.write(
        "<p><em>Accessibility arrays are provided to enable computing statistics "
        "with scikit-allel.</em></p>\n"
    )

    handle.write("<h2>Site density (100kb windows)</h2>\n")
    contig_lengths = _read_fai_lengths(str(snakemake.params.ref_fai))
    window = 100_000
    filtered_counts: dict[str, list[int]] = {}
    inv_counts: dict[str, list[int]] = {}
    variant_counts: dict[str, list[int]] = {}
    for contig in contigs:
        length = contig_lengths.get(contig)
        if length is None:
            continue
        n_windows = (length + window - 1) // window
        filtered_counts[contig] = [0 for _ in range(n_windows)]
        inv_counts[contig] = [0 for _ in range(n_windows)]
        variant_counts[contig] = [0 for _ in range(n_windows)]

    for path in snakemake.input.beds:
        try:
            _add_bed_counts(str(path), filtered_counts, contig_lengths, window)
        except OSError as exc:
            handle.write(
                f"<p>Failed to read filtered bed from {html.escape(str(path))}: "
                f"{html.escape(str(exc))}</p>\n"
            )

    for contig in contigs:
        if contig not in inv_counts:
            continue
        prefix = split_prefixes.get(contig)
        if prefix is None:
            continue
        inv_bed = f"{prefix}.inv.bed"
        if Path(inv_bed).exists():
            _add_bed_counts(inv_bed, inv_counts, contig_lengths, window)
        elif Path(inv_bed + ".gz").exists():
            _add_bed_counts(inv_bed + ".gz", inv_counts, contig_lengths, window)

    for path in snakemake.input.invs:
        try:
            with _open_text(str(path)) as f_in:
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
                f"<p>Failed to read invariant sites from {html.escape(str(path))}: "
                f"{html.escape(str(exc))}</p>\n"
            )

    for clean_path in snakemake.input.cleans:
        try:
            with _open_text(str(clean_path)) as f_in:
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
                f"<p>Failed to read variable sites from {html.escape(str(clean_path))}: "
                f"{html.escape(str(exc))}</p>\n"
            )

    handle.write("<h2>Summary Plots</h2>\n")
    for contig in contigs:
        if (
            contig not in filtered_counts
            or contig not in inv_counts
            or contig not in variant_counts
        ):
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
