#!/usr/bin/env python3
"""Draw the ARGprep graphical abstract (docs/graphical_abstract.{svg,png})."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

OUT = Path(__file__).resolve().parent

# Okabe-Ito based palette (colour-blind safe)
INV = "#0072B2"      # retained invariant
VAR = "#E69F00"      # retained variant
MASK = "#8C8C8C"     # masked
SNP_FILL = "#FBE3B0"  # base differing from reference
UNAL = "#E4E4E4"     # unaligned in this sample
INK = "#222222"
SOFT = "#555555"
PANEL = "#F6F7F9"

plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none"})

fig = plt.figure(figsize=(14, 7.2))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 140)
ax.set_ylim(0, 72)
ax.set_aspect("equal")
ax.axis("off")


def panel(x, y, w, h, title):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.5",
                                fc=PANEL, ec="#D5D8DD", lw=1))
    ax.text(x + w / 2, y + h - 2.6, title, ha="center", va="center",
            fontsize=12.5, weight="bold", color=INK)


def arrow(x0, y0, x1, y1, color=SOFT, lw=2.2):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=18,
                                 color=color, lw=lw, shrinkA=0, shrinkB=0))


def file_icon(x, y, w, h, label, sub=None, color=INK, dashed=False):
    fold = 1.4
    xs = [x, x + w - fold, x + w, x + w, x, x]
    ys = [y + h, y + h, y + h - fold, y, y, y + h]
    ax.fill(xs, ys, fc="white", ec=color, lw=1.4, ls="--" if dashed else "-")
    ax.plot([x + w - fold, x + w - fold, x + w], [y + h, y + h - fold, y + h - fold],
            color=color, lw=1.2)
    ax.text(x + w + 1.2, y + h / 2 + (0.9 if sub else 0), label, va="center",
            fontsize=10, weight="bold", color=INK, family="DejaVu Sans Mono")
    if sub:
        ax.text(x + w + 1.2, y + h / 2 - 1.5, sub, va="center", fontsize=8.5, color=SOFT)


# ---------------------------------------------------------------- title
ax.text(70, 68.6, "ARGprep", ha="center", va="center", fontsize=22, weight="bold", color=INK)
ax.text(70, 64.9,
        "Pairwise whole-genome alignments  →  reference-anchored all-sites data for ARG inference",
        ha="center", va="center", fontsize=12.5, color=SOFT)

# ---------------------------------------------------------------- inputs
panel(2, 13, 31, 48.5, "Input")

# pairwise alignment glyphs: reference bar over sample bar, with synteny blocks
ax.text(4.2, 54.3, "AnchorWave MAFs", fontsize=10.5, weight="bold", color=INK)
ax.text(4.2, 52.2, "one per sample, each vs. the same reference", fontsize=8.5, color=SOFT)
blocks = {
    "sample A": [(0.00, 0.30), (0.34, 0.72), (0.76, 1.00)],
    "sample B": [(0.00, 0.18), (0.25, 0.60), (0.66, 0.95)],
    "sample C": [(0.05, 0.42), (0.48, 1.00)],
}
gx, gw = 5.0, 23.0
for i, (name, segs) in enumerate(blocks.items()):
    top = 48.5 - i * 7.3
    bot = top - 3.6
    ax.add_patch(Rectangle((gx, top), gw, 0.9, fc=INK, ec="none"))
    ax.add_patch(Rectangle((gx, bot), gw, 0.9, fc=VAR, ec="none"))
    for a, b in segs:
        ax.fill([gx + a * gw, gx + b * gw, gx + b * gw, gx + a * gw],
                [top, top, bot + 0.9, bot + 0.9], fc="#C9D6E8", ec="none")
    ax.text(gx + gw + 0.8, top + 0.45, "ref", va="center", fontsize=8, color=SOFT)
    ax.text(gx + gw + 0.8, bot + 0.45, name.split()[1], va="center", fontsize=8,
            color=SOFT, weight="bold")
ax.text(gx + gw / 2, 26.4, "⋮", ha="center", va="center", fontsize=14, color=SOFT)

file_icon(5, 19.0, 3.2, 4.0, "reference.fa", "reference FASTA")
file_icon(5, 14.4, 3.2, 4.0, "<sample>.bed", "optional quality masks", color=SOFT, dashed=True)

arrow(33.6, 37, 37.4, 37)

# ---------------------------------------------------------------- core: per-site calling
panel(38, 13, 62, 48.5, "Project onto reference coordinates & call every site")

ref = "ACGTTAGCATNCGATC"
rows = {
    "A": "ACGTTAGCATNCGATC",
    "B": "ACTT....AATCAATC",
    "C": "ACGTTAGAATNCGA-C",
    "D": "ACGTT--CANNCGATT",
    "E": "ACGTTAGGATNCAATC",
}
# column status: (kind, tag)
status = [
    ("inv", ""), ("inv", ""), ("var", ""), ("inv", ""), ("inv", ""),
    ("mask", ""), ("mask", ""), ("mask", "multi-\nallelic*"), ("inv", ""),
    ("var", ""), ("mask", "ref N"), ("inv", ""), ("var", ""), ("inv", ""),
    ("inv", ""), ("mask", "indel-\nadjacent*"),
]
assert len(status) == len(ref) and all(len(r) == len(ref) for r in rows.values())

cw, ch = 3.2, 3.1
mx = 45.0
my_top = 53.0
ax.text(mx - 1.2, my_top + ch / 2, "ref", ha="right", va="center", fontsize=10,
        weight="bold", color=INK)
for j, b in enumerate(ref):
    ax.add_patch(Rectangle((mx + j * cw, my_top), cw, ch, fc="#DDE3EA", ec="white", lw=1.2))
    ax.text(mx + (j + 0.5) * cw, my_top + ch / 2, b, ha="center", va="center",
            fontsize=11, weight="bold", family="DejaVu Sans Mono", color=INK)
for i, (name, seq) in enumerate(rows.items()):
    y = my_top - (i + 1) * ch - 0.6
    ax.text(mx - 1.2, y + ch / 2, name, ha="right", va="center", fontsize=10, color=SOFT)
    for j, b in enumerate(seq):
        if b == ".":
            fc, txt = UNAL, ""
        elif b in "-N":
            fc, txt = UNAL, b
        elif b != ref[j] and ref[j] != "N":
            fc, txt = SNP_FILL, b
        else:
            fc, txt = "white", b
        ax.add_patch(Rectangle((mx + j * cw, y), cw, ch, fc=fc, ec="#D5D8DD", lw=0.8))
        if txt:
            ax.text(mx + (j + 0.5) * cw, y + ch / 2, txt, ha="center", va="center",
                    fontsize=10.5, family="DejaVu Sans Mono", color=INK)

# status strip
sy = my_top - 6 * ch - 1.8
kind_color = {"inv": INV, "var": VAR, "mask": MASK}
for j, (kind, tag) in enumerate(status):
    ax.add_patch(Rectangle((mx + j * cw + 0.15, sy), cw - 0.3, 1.6, fc=kind_color[kind], ec="none"))
    if tag:
        ax.text(mx + (j + 0.5) * cw, sy - 0.8, tag, ha="center", va="top", fontsize=7.2,
                color=SOFT, linespacing=0.95)
# columns 6-7 share one tag: B unaligned + D deleted = 2 missing
ax.text(mx + 6 * cw, sy - 0.8, "too many\nmissing", ha="center", va="top", fontsize=7.2,
        color=SOFT, linespacing=0.95)
ax.text(mx - 1.2, sy + 0.8, "call", ha="right", va="center", fontsize=9, weight="bold", color=INK)

# filter rules, below the call strip
fy = 27.6
ax.text(mx, fy, "Per-site filters", fontsize=10, weight="bold", color=INK, va="center")
rules = [
    ("reference base must be A/C/G/T", fy - 2.6),
    ("missing samples ≤ threshold (here 1)\n   unaligned, gap, N, low quality = missing", fy - 5.6),
]
opt = [
    ("*optional: mask SNPs flanking an indel", fy - 2.6),
    ("*optional: mask multiallelic sites", fy - 4.6),
]
for txt, yy in rules:
    ax.text(mx, yy, "• " + txt, fontsize=8.6, color=SOFT, va="center", linespacing=1.15)
for txt, yy in opt:
    ax.text(mx + 27, yy, txt, fontsize=8.6, color=SOFT, va="center")

# legend: cell shading, then call colours
legend = [
    [(SNP_FILL, "differs from ref", "#B8BCC2"), (UNAL, "unaligned / gap / N", "#B8BCC2")],
    [(INV, "retained invariant", "none"), (VAR, "retained variant", "none"), (MASK, "masked", "none")],
]
for r, items in enumerate(legend):
    yy = 17.6 - r * 2.9
    for k, (fc, lab, ec) in enumerate(items):
        lx = mx + k * 17.0
        ax.add_patch(Rectangle((lx, yy), 1.8, 1.8, fc=fc, ec=ec, lw=0.8))
        ax.text(lx + 2.4, yy + 0.9, lab, va="center", fontsize=8.3, color=SOFT)

arrow(100.6, 37, 104.4, 37)

# ---------------------------------------------------------------- outputs
panel(105, 13, 33, 48.5, "Per-contig outputs")
oy = 53.0
outs = [
    ("<contig>.all_sites.vcf", "invariant + variant sites", INV),
    ("<contig>.vcf", "variant sites only", VAR),
    ("<contig>.mask.bed", "every masked position", MASK),
    ("<contig>.sites", "ARGweaver format (optional)", VAR),
    ("summary.html", "genome-wide QC report", INK),
]
for k, (name, sub, col) in enumerate(outs):
    y = oy - k * 5.6
    file_icon(107.5, y - 2, 3.2, 4.0, name, sub, color=col, dashed=name.endswith(".sites"))

# tiling invariant
ax.text(121.5, 26.0, "retained sites + mask tile\neach contig exactly", ha="center", va="center",
        fontsize=8.5, color=SOFT, style="italic", linespacing=1.1)
tiles = [(INV, 4), (MASK, 1.5), (INV, 3), (VAR, 0.5), (INV, 2.5), (MASK, 3), (INV, 4),
         (VAR, 0.5), (INV, 3), (MASK, 1), (INV, 2.5)]
tx = 108.5
for col, w in tiles:
    ax.add_patch(Rectangle((tx, 21.6), w, 1.6, fc=col, ec="none"))
    tx += w
ax.text(121.5, 16.6, "→ Relate · ARGweaver · SINGER", ha="center", va="center",
        fontsize=10, weight="bold", color=INV)

# ---------------------------------------------------------------- workflow band
steps = [
    ("prepare_reference", "index FASTA, reconcile\ncontig names"),
    ("split_sample_maf", "partition each MAF\nby contig (gzip)"),
    ("direct_maf_sites", "one job per contig,\nsingle pass over MAFs"),
    ("summary_report", "aggregate stats,\nno VCF rescan"),
]
bw, gap, bx0, by = 27.5, 4.0, 4.0, 2.0
for k, (name, sub) in enumerate(steps):
    x = bx0 + k * (bw + gap)
    ax.add_patch(FancyBboxPatch((x, by), bw, 8.6, boxstyle="round,pad=0,rounding_size=1.2",
                                fc="white", ec="#B8BCC2", lw=1))
    ax.text(x + bw / 2, by + 6.3, name, ha="center", va="center", fontsize=10,
            weight="bold", family="DejaVu Sans Mono", color=INK)
    ax.text(x + bw / 2, by + 2.8, sub, ha="center", va="center", fontsize=8.3,
            color=SOFT, linespacing=1.1)
    if k < len(steps) - 1:
        arrow(x + bw + 0.4, by + 4.3, x + bw + gap - 0.4, by + 4.3, lw=1.6)
x_end = bx0 + 4 * (bw + gap) - gap
ax.text(x_end + 1.6, by + 6.6, "Snakemake", fontsize=9.5, weight="bold", color=INK)
ax.text(x_end + 1.6, by + 4.0, "local or SLURM\n(preemption-safe)", fontsize=8.3, color=SOFT,
        va="center", linespacing=1.1)


for ext in ("svg", "png"):
    fig.savefig(OUT / f"graphical_abstract.{ext}", dpi=200)
print("wrote", OUT / "graphical_abstract.svg", OUT / "graphical_abstract.png")
