#!/usr/bin/env python3
"""
filtered_vcf_to_bed.py

Read a *.filtered VCF (or VCF-like) file and write a BED file containing the
basepair positions of each record. Dropped-indels and missing BED files are
included based on the shared prefix.

- Skips VCF header lines beginning with '#'
- Uses CHROM (col 1) and POS (col 2)
- Outputs BED intervals that represent a single bp:
    start = POS-1
    end   = POS
  (BED is 0-based, half-open; a single base at POS becomes [POS-1, POS))

- Optionally sorts and merges overlapping/adjacent intervals per chromosome.

Example:
  python filtered_vcf_to_bed.py /path/to/sample.gvcf.gz
  python filtered_vcf_to_bed.py sample.gvcf --no-merge
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
from typing import TextIO, Dict, List, Tuple


def open_maybe_gzip(path: str, mode: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode)  # type: ignore
    return open(path, mode, encoding="utf-8")


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge overlapping or adjacent half-open intervals (start, end).
    Assumes intervals are sorted by (start, end).
    """
    if not intervals:
        return []
    merged: List[Tuple[int, int]] = []
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:  # overlap or adjacency (since half-open, adjacency means s == cur_e)
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged


def find_gvcf(path_or_prefix: str) -> str | None:
    """
    Resolve a gVCF path. Accepts:
      - an explicit file path, or
      - a prefix (tries .gvcf/.vcf with optional .gz).
    """
    if os.path.isfile(path_or_prefix):
        return path_or_prefix
    candidates = [
        path_or_prefix + ".gvcf.gz",
        path_or_prefix + ".gvcf",
        path_or_prefix + ".vcf.gz",
        path_or_prefix + ".vcf",
    ]
    for path in candidates:
        if os.path.isfile(path):
            return path
    return None


def prefix_from_gvcf(path_or_prefix: str) -> str:
    """
    Derive a shared prefix from a gVCF/VCF filename.
    If no known extension is present, return the string as-is.
    """
    path = path_or_prefix
    for suffix in (".gvcf.gz", ".gvcf", ".vcf.gz", ".vcf"):
        if path.endswith(suffix):
            return path[: -len(suffix)]
    return path


def read_bed_intervals(path: str, by_chrom: Dict[str, List[Tuple[int, int]]]) -> None:
    """
    Load a BED file into a per-chromosome interval map.
    Assumes 0-based half-open intervals [start, end).
    """
    with open(path, "rt", encoding="utf-8") as fin:
        for raw in fin:
            if not raw or raw.startswith("#"):
                continue
            line = raw.rstrip("\n")
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            chrom = cols[0]
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                continue
            if end <= start:
                continue
            by_chrom.setdefault(chrom, []).append((start, end))


def extract_end(info: str) -> int | None:
    """
    Parse END= from a VCF INFO field, if present.
    """
    if info == ".":
        return None
    for field in info.split(";"):
        if field.startswith("END="):
            try:
                return int(field.split("=", 1)[1])
            except ValueError:
                return None
    return None


def sum_merged_bp(by_chrom: Dict[str, List[Tuple[int, int]]]) -> int:
    """
    Merge intervals and return the total covered basepairs.
    """
    total = 0
    for chrom in by_chrom:
        intervals = by_chrom[chrom]
        if not intervals:
            continue
        intervals.sort(key=lambda x: (x[0], x[1]))
        merged = merge_intervals(intervals)
        total += sum(e - s for s, e in merged)
    return total


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert a .filtered VCF to a BED of bp positions.")
    ap.add_argument(
        "gvcf",
        help="Path to the gVCF/VCF (or its prefix) used to derive related files",
    )
    ap.add_argument(
        "--dropped-bed",
        default=None,
        help="Optional path to dropped_indels.bed (overrides the default cleangVCF/dropped_indels.bed).",
    )
    ap.add_argument(
        "--no-merge",
        action="store_true",
        help="Do not sort/merge overlapping or adjacent intervals.",
    )
    args = ap.parse_args()

    # Resolve the gVCF path and derive the shared prefix for all files.
    gvcf_arg = args.gvcf
    gvcf_path = find_gvcf(gvcf_arg)
    if gvcf_path is None:
        # Allow passing a direct prefix (e.g., results/split/combined.1).
        prefix = gvcf_arg
        filtered_path = prefix + ".filtered"
        filtered_gz = filtered_path + ".gz"
        if not (os.path.isfile(filtered_path) or os.path.isfile(filtered_gz)):
            sys.stderr.write(
                f"ERROR: gVCF not found: '{gvcf_arg}' (expected file or .gvcf/.vcf extension).\n"
            )
            sys.exit(1)
    else:
        prefix = prefix_from_gvcf(gvcf_path)
    filtered_path = prefix + ".filtered"
    filtered_gz = filtered_path + ".gz"
    if os.path.isfile(filtered_gz):
        filtered_path = filtered_gz
    out_path = prefix + ".filtered.bed"
    dropped_bed = args.dropped_bed or os.path.join(os.path.dirname(prefix), "cleangVCF", "dropped_indels.bed")
    missing_bed = prefix + ".missing.bed"

    # Collect intervals by chromosome from the filtered VCF and mask BEDs.
    by_chrom: Dict[str, List[Tuple[int, int]]] = {}

    if not os.path.isfile(filtered_path):
        sys.stderr.write(f"ERROR: filtered VCF not found: '{filtered_path}'.\n")
        sys.exit(1)
    if not os.path.isfile(missing_bed):
        sys.stderr.write(f"ERROR: missing BED not found: '{missing_bed}'.\n")
        sys.exit(1)
    if not os.path.isfile(dropped_bed):
        sys.stderr.write(
            f"ERROR: dropped indels BED not found: '{dropped_bed}'.\n"
        )
        sys.exit(1)

    with open_maybe_gzip(filtered_path, "rt") as fin:
        for raw in fin:
            if not raw or raw.startswith("#"):
                continue
            line = raw.rstrip("\n")
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            chrom = cols[0]
            try:
                pos = int(cols[1])
            except ValueError:
                continue
            start = pos - 1
            end = pos
            by_chrom.setdefault(chrom, []).append((start, end))

    # Add dropped-indel and missing-position masks.
    read_bed_intervals(dropped_bed, by_chrom)
    read_bed_intervals(missing_bed, by_chrom)

    # Default behavior is to sort + merge unless --no-merge is given.
    with open(out_path, "wt", encoding="utf-8") as fout:
        for chrom in sorted(by_chrom.keys()):
            intervals = by_chrom[chrom]
            intervals.sort(key=lambda x: (x[0], x[1]))
            output_intervals = intervals if args.no_merge else merge_intervals(intervals)
            for s, e in output_intervals:
                fout.write(f"{chrom}\t{s}\t{e}\n")


if __name__ == "__main__":
    main()
