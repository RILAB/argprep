#!/usr/bin/env python3
"""Make per-chromosome synteny dotplots from AnchorWave pairwise MAFs.

Each input MAF is one query sample aligned to the reference. Every alignment
block contributes one segment to the dotplot: reference interval on the x-axis,
query interval on the y-axis. Forward-strand blocks run along the diagonal;
reverse-strand blocks run anti-diagonal (inversions).

For each sample a single multi-page PDF is written, one chromosome per page.

Only the six leading coordinate fields of each `s` line are read (via awk), so
the large sequence columns are never pulled into Python -- this stays fast and
light even on multi-GB MAFs.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

MB = 1_000_000.0

# awk program: emit tab-separated coordinate fields for every `s` line,
# skipping the huge alignment sequence column entirely.
_AWK = r'/^s/ {print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}'


def _natural_key(name: str):
    """Sort chr1, chr2, ..., chr10 in numeric rather than lexical order."""
    m = re.search(r"(\d+)$", name)
    return (0, int(m.group(1))) if m else (1, name)


def parse_blocks(maf_path: Path):
    """Return {ref_contig: [block, ...]} where each block is a dict.

    Blocks are read in pairs of `s` lines: first line = reference, second =
    query. `size`/`start`/`srcSize` are ints; `strand` is '+' or '-'.
    """
    proc = subprocess.Popen(
        ["awk", _AWK, str(maf_path)],
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1 << 20,
    )
    by_chrom: dict[str, list] = {}
    pending = None  # holds the reference line until its query line arrives
    assert proc.stdout is not None
    for line in proc.stdout:
        name, start, size, strand, srcsize = line.split("\t")
        rec = (name, int(start), int(size), strand, int(srcsize))
        if pending is None:
            pending = rec
            continue
        r_name, r_start, r_size, _r_strand, r_src = pending
        q_name, q_start, q_size, q_strand, q_src = rec
        by_chrom.setdefault(r_name, []).append(
            {
                "r_start": r_start,
                "r_size": r_size,
                "r_src": r_src,
                "q_start": q_start,
                "q_size": q_size,
                "q_src": q_src,
                "q_strand": q_strand,
            }
        )
        pending = None
    ret = proc.wait()
    if ret != 0:
        raise RuntimeError(f"awk failed ({ret}) on {maf_path}")
    if pending is not None:
        sys.stderr.write(f"warning: odd number of s-lines in {maf_path}\n")
    return by_chrom


def _segment(block):
    """Endpoints (x0, y0, x1, y1) in bp of a block's dotplot segment.

    For reverse-strand blocks the MAF start is measured on the
    reverse-complemented source, so convert to forward coordinates and let the
    query coordinate decrease as the reference increases.
    """
    x0 = block["r_start"]
    x1 = block["r_start"] + block["r_size"]
    if block["q_strand"] == "+":
        y0 = block["q_start"]
        y1 = block["q_start"] + block["q_size"]
    else:
        y0 = block["q_src"] - block["q_start"]
        y1 = block["q_src"] - block["q_start"] - block["q_size"]
    return x0, y0, x1, y1


def plot_sample(maf_path: Path, out_dir: Path):
    sample = maf_path.name
    for suffix in (".maf.gz", ".maf"):
        if sample.endswith(suffix):
            sample = sample[: -len(suffix)]
            break
    by_chrom = parse_blocks(maf_path)
    if not by_chrom:
        return sample, 0, "no alignment blocks found"

    out_path = out_dir / f"{sample}.dotplot.pdf"
    chroms = sorted(by_chrom, key=_natural_key)
    with PdfPages(out_path) as pdf:
        for chrom in chroms:
            blocks = by_chrom[chrom]
            fig, ax = plt.subplots(figsize=(7.0, 7.0))
            r_src = blocks[0]["r_src"]
            q_src = max(b["q_src"] for b in blocks)
            for b in blocks:
                x0, y0, x1, y1 = _segment(b)
                color = "#1f77b4" if b["q_strand"] == "+" else "#d62728"
                ax.plot(
                    [x0 / MB, x1 / MB],
                    [y0 / MB, y1 / MB],
                    color=color,
                    linewidth=1.2,
                    solid_capstyle="round",
                )
            ax.set_xlim(0, r_src / MB)
            ax.set_ylim(0, q_src / MB)
            ax.set_xlabel(f"reference {chrom} (Mb)")
            ax.set_ylabel(f"{sample} {chrom} (Mb)")
            ax.set_title(f"{sample} — {chrom}  ({len(blocks)} blocks)")
            ax.grid(True, linewidth=0.3, alpha=0.4)
            ax.legend(
                handles=[
                    Line2D([0], [0], color="#1f77b4", label="forward (+)"),
                    Line2D([0], [0], color="#d62728", label="reverse (−)"),
                ],
                loc="lower right",
                fontsize=8,
                framealpha=0.9,
            )
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    return sample, len(chroms), str(out_path)


def _worker(args):
    maf_path, out_dir = args
    try:
        return plot_sample(Path(maf_path), Path(out_dir))
    except Exception as exc:  # keep one bad file from killing the batch
        return Path(maf_path).name, -1, f"ERROR: {exc}"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--maf-dir", help="directory of *.maf / *.maf.gz files")
    src.add_argument("--maf", nargs="+", help="explicit list of MAF files")
    ap.add_argument("--out-dir", required=True, help="output directory for PDFs")
    ap.add_argument("--procs", type=int, default=1, help="parallel workers")
    args = ap.parse_args()

    if args.maf_dir:
        d = Path(args.maf_dir)
        mafs = sorted(
            p for p in d.iterdir() if p.name.endswith((".maf", ".maf.gz"))
        )
    else:
        mafs = [Path(p) for p in args.maf]
    if not mafs:
        sys.exit("no MAF files found")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[maf_dotplot] {len(mafs)} samples -> {out_dir} ({args.procs} procs)")
    tasks = [(str(m), str(out_dir)) for m in mafs]
    with ProcessPoolExecutor(max_workers=args.procs) as ex:
        futures = {ex.submit(_worker, t): t[0] for t in tasks}
        for fut in as_completed(futures):
            sample, npages, info = fut.result()
            if npages < 0:
                print(f"  [FAIL] {sample}: {info}")
            else:
                print(f"  [ok]   {sample}: {npages} chromosomes -> {info}")
    print("[maf_dotplot] done")


if __name__ == "__main__":
    main()
