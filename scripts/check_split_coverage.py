#!/usr/bin/env python3
"""
Check that split outputs cover the full contig length.

This script sums:
  - .clean records (1 bp each)
  - .inv records (1 bp each; END spans are expanded by split.py)
  - .filtered.bed intervals (0-based, half-open)

and compares the total to the contig length from a reference .fai.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import Dict, List, TextIO, Tuple


def open_maybe_gzip(path: Path, mode: str) -> TextIO:
    # Normalize IO for optional gzipped inputs.
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)  # type: ignore
    return path.open(mode, encoding="utf-8")


def extract_end(info: str) -> int | None:
    # Parse END= from INFO, if present.
    if info == ".":
        return None
    for field in info.split(";"):
        if field.startswith("END="):
            try:
                return int(field.split("=", 1)[1])
            except ValueError:
                return None
    return None


def load_fai_length(path: Path, chrom: str) -> int:
    # Resolve contig length from the .fai index.
    with path.open("r", encoding="utf-8") as fin:
        for raw in fin:
            if not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            if cols[0] == chrom:
                return int(cols[1])
    raise ValueError(f"Chromosome '{chrom}' not found in {path}")


def read_clean_inv_intervals(path: Path) -> tuple[str | None, List[Tuple[int, int]]]:
    # Convert 1-based VCF positions to 0-based half-open intervals.
    chrom = None
    intervals: List[Tuple[int, int]] = []
    with open_maybe_gzip(path, "rt") as fin:
        for raw in fin:
            if raw.startswith("#") or not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            if chrom is None:
                chrom = cols[0]
            try:
                pos = int(cols[1])
            except ValueError:
                continue
            intervals.append((pos - 1, pos))
    return chrom, intervals


def read_bed_intervals(path: Path) -> tuple[str | None, List[Tuple[int, int]]]:
    # Read BED intervals directly (already 0-based half-open).
    chrom = None
    intervals: List[Tuple[int, int]] = []
    with path.open("r", encoding="utf-8") as fin:
        for raw in fin:
            if not raw.strip() or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            if chrom is None:
                chrom = cols[0]
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                continue
            if end > start:
                intervals.append((start, end))
    return chrom, intervals


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    # Coalesce overlapping/adjacent intervals to simplify coverage math.
    if not intervals:
        return []
    intervals.sort(key=lambda x: (x[0], x[1]))
    merged: List[Tuple[int, int]] = []
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged


def overlap_bp(a: List[Tuple[int, int]], b: List[Tuple[int, int]]) -> int:
    # Compute total bp overlap between two interval sets.
    i = j = 0
    total = 0
    a = merge_intervals(a)
    b = merge_intervals(b)
    while i < len(a) and j < len(b):
        a_s, a_e = a[i]
        b_s, b_e = b[j]
        if a_e <= b_s:
            i += 1
            continue
        if b_e <= a_s:
            j += 1
            continue
        total += min(a_e, b_e) - max(a_s, b_s)
        if a_e <= b_e:
            i += 1
        else:
            j += 1
    return total


def overlap_intervals(
    a: List[Tuple[int, int]], b: List[Tuple[int, int]]
) -> List[Tuple[int, int]]:
    # Return explicit overlap intervals for reporting.
    i = j = 0
    overlaps: List[Tuple[int, int]] = []
    a = merge_intervals(a)
    b = merge_intervals(b)
    while i < len(a) and j < len(b):
        a_s, a_e = a[i]
        b_s, b_e = b[j]
        if a_e <= b_s:
            i += 1
            continue
        if b_e <= a_s:
            j += 1
            continue
        start = max(a_s, b_s)
        end = min(a_e, b_e)
        if end > start:
            overlaps.append((start, end))
        if a_e <= b_e:
            i += 1
        else:
            j += 1
    return overlaps


def format_overlap_report(
    chrom: str,
    label: str,
    overlaps: List[Tuple[int, int]],
    file_a: Path,
    file_b: Path,
    max_show: int = 20,
) -> str:
    # Pretty-print overlap intervals with file context for debugging.
    count = len(overlaps)
    if count == 0:
        return f"{label}: 0 overlaps"
    show = overlaps[:max_show]
    lines = [
        f"{label}: {count} overlap intervals (showing up to {max_show})",
        f"  files: {file_a} vs {file_b}",
    ]
    for start, end in show:
        lines.append(f"  {chrom}:{start}-{end} overlaps {label} ({file_a} vs {file_b})")
    if count > max_show:
        lines.append(f"  ... {count - max_show} more")
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description="Check split coverage against reference length.")
    ap.add_argument("prefix", help="Prefix for split outputs (e.g., results/split/combined.1)")
    ap.add_argument("--fai", required=True, help="Reference .fai path for contig lengths")
    args = ap.parse_args()

    prefix = Path(args.prefix)
    clean = prefix.with_suffix(prefix.suffix + ".clean")
    inv = prefix.with_suffix(prefix.suffix + ".inv")
    filtered_bed = prefix.with_suffix(prefix.suffix + ".filtered.bed")

    if clean.with_suffix(clean.suffix + ".gz").exists():
        clean = clean.with_suffix(clean.suffix + ".gz")
    if inv.with_suffix(inv.suffix + ".gz").exists():
        inv = inv.with_suffix(inv.suffix + ".gz")

    if not clean.exists():
        sys.exit(f"ERROR: clean file not found: {clean}")
    if not inv.exists():
        sys.exit(f"ERROR: inv file not found: {inv}")
    if not filtered_bed.exists():
        sys.exit(f"ERROR: filtered bed not found: {filtered_bed}")

    clean_chrom, clean_intervals = read_clean_inv_intervals(clean)
    inv_chrom, inv_intervals = read_clean_inv_intervals(inv)
    bed_chrom, bed_intervals = read_bed_intervals(filtered_bed)

    # Determine the contig name from whichever input provides one.
    chrom = clean_chrom or inv_chrom or bed_chrom
    if chrom is None:
        sys.exit("ERROR: unable to determine chromosome from inputs.")

    for name, c in (("clean", clean_chrom), ("inv", inv_chrom), ("bed", bed_chrom)):
        if c is not None and c != chrom:
            sys.exit(f"ERROR: mismatched chromosome in {name}: {c} != {chrom}")

    fai = Path(args.fai)
    chrom_len = load_fai_length(fai, chrom)

    # Overlaps indicate split outputs are not mutually exclusive.
    overlap_ci = overlap_bp(clean_intervals, inv_intervals)
    overlap_cb = overlap_bp(clean_intervals, bed_intervals)
    overlap_ib = overlap_bp(inv_intervals, bed_intervals)
    if overlap_ci or overlap_cb or overlap_ib:
        ci_intervals = overlap_intervals(clean_intervals, inv_intervals)
        cb_intervals = overlap_intervals(clean_intervals, bed_intervals)
        ib_intervals = overlap_intervals(inv_intervals, bed_intervals)
        report_lines = [
            f"ERROR: overlap detected for {chrom}: "
            f"clean vs inv={overlap_ci}, clean vs filtered_bed={overlap_cb}, "
            f"inv vs filtered_bed={overlap_ib}",
            format_overlap_report(chrom, "clean vs inv", ci_intervals, clean, inv),
            format_overlap_report(
                chrom, "clean vs filtered_bed", cb_intervals, clean, filtered_bed
            ),
            format_overlap_report(
                chrom, "inv vs filtered_bed", ib_intervals, inv, filtered_bed
            ),
        ]
        sys.exit("\n".join(report_lines))

    # Coverage should match the full contig length after unioning intervals.
    merged = merge_intervals(clean_intervals + inv_intervals + bed_intervals)
    total = sum(e - s for s, e in merged)
    if total != chrom_len:
        sys.exit(
            f"ERROR: coverage mismatch for {chrom}: "
            f"union={total}, chrom_len={chrom_len}"
        )

    report = prefix.with_suffix(prefix.suffix + ".coverage.txt")
    report.write_text(
        f"chrom={chrom}\n"
        f"clean_bp={sum(e - s for s, e in merge_intervals(clean_intervals))}\n"
        f"inv_bp={sum(e - s for s, e in merge_intervals(inv_intervals))}\n"
        f"filtered_bed_bp={sum(e - s for s, e in merge_intervals(bed_intervals))}\n"
        f"total_bp={total}\n"
        f"chrom_len={chrom_len}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
