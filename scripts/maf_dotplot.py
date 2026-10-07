#!/usr/bin/env python3
"""Per-reference-chromosome synteny dotplots from pairwise reference-vs-query MAFs.

Every alignment block is one segment: reference interval on the x-axis, query
interval on the y-axis. Forward-strand blocks run along the diagonal (blue);
reverse-strand blocks run anti-diagonal (red, inversions).

One reference chromosome can align to several query contigs (translocations,
scaffolds). Each query contig gets its own labelled horizontal band on the
y-axis, so unrelated coordinate systems are never overlaid.

The MAF goes through the same validated parser as ``maf_stats.py`` (plain or
gzip-compressed input). The QC workflow calls :func:`render_dotplot` from
inside ``maf_stats.py``, so plotting shares the statistics pass; this CLI is for
ad-hoc use and writes one PNG per aligned reference contig.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

MB = 1_000_000.0
FORWARD_COLOR = "#1f77b4"
REVERSE_COLOR = "#d62728"
FIGSIZE_IN = (8.0, 8.0)
DPI = 100  # 800 x 800 px: predictable size for the HTML report
MAX_BAND_LABELS = 40


def _natural_key(name: str):
    """Sort chr1, chr2, ..., chr10 numerically rather than lexically."""
    m = re.search(r"(\d+)$", name)
    return (0, name[: m.start()], int(m.group(1))) if m else (1, name, 0)


def _band_order(
    blocks: Iterable, query_order: list[str] | None
) -> list[str]:
    present = {block.query_contig for block in blocks}
    if query_order is not None:
        ordered = [name for name in query_order if name in present]
        ordered += sorted(present.difference(ordered), key=_natural_key)
        return ordered
    return sorted(present, key=_natural_key)


def segment(block, band_offset: int) -> tuple[tuple[float, float], tuple[float, float]]:
    """Segment endpoints in Mb for one block. Query coordinates are already
    forward-strand; a reverse block runs from its query end down to its start."""
    x0, x1 = block.ref_start / MB, block.ref_end / MB
    if block.strand == "+":
        y0, y1 = block.query_start, block.query_end
    else:
        y0, y1 = block.query_end, block.query_start
    return (x0, (band_offset + y0) / MB), (x1, (band_offset + y1) / MB)


def render_dotplot(
    out_png: Path,
    *,
    sample: str,
    ref_contig: str,
    ref_length: int,
    blocks: list,
    query_lengths: dict[str, int],
    query_order: list[str] | None = None,
) -> None:
    """Write one PNG for one reference contig.

    ``blocks`` are ``maf_stats.Block`` objects (forward query coordinates) on
    ``ref_contig``. ``query_lengths`` sizes each query band; ``query_order``
    (e.g. query .fai order) orders the bands, falling back to natural sort.
    """
    bands = _band_order(blocks, query_order)
    offsets: dict[str, int] = {}
    total = 0
    for name in bands:
        offsets[name] = total
        total += query_lengths.get(name, 0)

    forward, reverse = [], []
    for block in blocks:
        seg = segment(block, offsets[block.query_contig])
        (forward if block.strand == "+" else reverse).append(seg)

    fig, ax = plt.subplots(figsize=FIGSIZE_IN, dpi=DPI)
    try:
        for segments, color in ((forward, FORWARD_COLOR), (reverse, REVERSE_COLOR)):
            if segments:
                collection = LineCollection(segments, colors=color, linewidths=1.2, capstyle="round")
                collection.set_rasterized(True)
                ax.add_collection(collection)

        ax.set_xlim(0, max(ref_length, 1) / MB)
        ax.set_ylim(0, max(total, 1) / MB)
        ax.set_xlabel(f"reference {ref_contig} (Mb)")
        ax.grid(True, linewidth=0.3, alpha=0.4)

        for name in bands[1:]:
            ax.axhline(offsets[name] / MB, color="0.6", linewidth=0.6)
        if len(bands) == 1:
            ax.set_ylabel(f"{sample} {bands[0]} (Mb)")
        else:
            ax.set_ylabel(f"{sample} query contigs (stacked, Mb)")
            labelled = sorted(bands, key=lambda n: -query_lengths.get(n, 0))[:MAX_BAND_LABELS]
            for name in labelled:
                mid = (offsets[name] + query_lengths.get(name, 0) / 2) / MB
                ax.text(
                    1.01, mid, name, transform=ax.get_yaxis_transform(),
                    fontsize=6, va="center", ha="left", clip_on=False,
                )
        ax.set_title(f"{sample} vs {ref_contig}  ({len(blocks)} blocks, {len(bands)} query contigs)")
        ax.legend(
            handles=[
                Line2D([0], [0], color=FORWARD_COLOR, label="forward (+)"),
                Line2D([0], [0], color=REVERSE_COLOR, label="reverse (−)"),
            ],
            loc="lower right",
            fontsize=8,
            framealpha=0.9,
        )
        fig.tight_layout()
        out_png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_png, dpi=DPI)
    finally:
        plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    try:
        from scripts import maf_stats
    except ModuleNotFoundError:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
        from scripts import maf_stats

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--maf", nargs="+", required=True, help="Pairwise MAF files (.maf / .maf.gz)")
    ap.add_argument("--reference-fai", required=True, help="Reference .fai index")
    ap.add_argument("--out-dir", required=True, help="Output directory; PNGs go to <out-dir>/dotplots/<sample>/")
    ap.add_argument("--query-fai", default=None, help="Query .fai (only valid with a single --maf)")
    args = ap.parse_args(argv)
    if args.query_fai and len(args.maf) > 1:
        ap.error("--query-fai applies to one query genome; pass a single --maf")

    reference_lengths = maf_stats.read_fai_lengths(Path(args.reference_fai))
    query_lengths = maf_stats.read_fai_lengths(Path(args.query_fai)) if args.query_fai else None
    failures = 0
    for maf in map(Path, args.maf):
        sample = re.sub(r"\.maf(\.gz)?$", "", maf.name)
        try:
            scan = maf_stats.scan_maf(maf, reference_lengths, query_lengths)
            contigs = [name for name, stats in scan.reference.items() if stats.block_sizes]
            paths = maf_stats.write_dotplots(Path(args.out_dir), sample, scan, contigs, query_lengths)
            print(f"[maf_dotplot] {sample}: {len(paths)} plots", file=sys.stderr)
        except Exception as exc:  # report every file, then fail overall
            failures += 1
            print(f"[maf_dotplot] FAIL {maf}: {exc}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
