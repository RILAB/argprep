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


def subtract_intervals(
    intervals: List[Tuple[int, int]],
    subtracts: List[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    """
    Subtract half-open intervals in 'subtracts' from 'intervals'.
    Inputs must be sorted and non-overlapping within each list.
    """
    if not intervals or not subtracts:
        return intervals
    result: List[Tuple[int, int]] = []
    i = 0
    j = 0
    while i < len(intervals):
        s, e = intervals[i]
        cur = s
        while j < len(subtracts) and subtracts[j][1] <= s:
            j += 1
        k = j
        while k < len(subtracts) and subtracts[k][0] < e:
            sub_s, sub_e = subtracts[k]
            if sub_s > cur:
                result.append((cur, min(sub_s, e)))
            cur = max(cur, sub_e)
            if cur >= e:
                break
            k += 1
        if cur < e:
            result.append((cur, e))
        i += 1
    return result


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


def read_vcf_intervals(path: str, by_chrom: Dict[str, List[Tuple[int, int]]]) -> None:
    """
    Load VCF records as half-open intervals, expanding END= spans or REF length.
    """
    with open_maybe_gzip(path, "rt") as fin:
        for raw in fin:
            if not raw or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            chrom = cols[0]
            try:
                pos = int(cols[1])
            except ValueError:
                continue
            ref = cols[3] if len(cols) >= 4 else "N"
            info = cols[7] if len(cols) >= 8 else "."
            end_val = extract_end(info)
            if end_val is None:
                span = max(len(ref), 1)
                end_val = pos + span - 1
            start = pos - 1
            end = end_val
            if end > start:
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
    inv_path = prefix + ".inv"
    inv_gz = inv_path + ".gz"
    if os.path.isfile(inv_gz):
        inv_path = inv_gz

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
            ref = cols[3] if len(cols) >= 4 else "N"
            info = cols[7] if len(cols) >= 8 else "."
            end_val = extract_end(info)
            if end_val is None:
                span = max(len(ref), 1)
                end_val = pos + span - 1
            start = pos - 1
            end = end_val
            by_chrom.setdefault(chrom, []).append((start, end))

    filtered_chroms: List[str] = []
    if by_chrom:
        filtered_chroms = sorted(by_chrom.keys())
    if len(filtered_chroms) > 1:
        sys.stderr.write(
            f"ERROR: filtered VCF has multiple chromosomes: {', '.join(filtered_chroms)}\n"
        )
        sys.exit(1)
    target_chrom = filtered_chroms[0] if filtered_chroms else None

    # If filtered is empty, try to infer target from other inputs.
    if target_chrom is None and os.path.isfile(inv_path):
        target_chrom = first_chrom_from_vcf(inv_path)
    clean_path = prefix + ".clean"
    clean_gz = clean_path + ".gz"
    if os.path.isfile(clean_gz):
        clean_path = clean_gz
    if target_chrom is None and os.path.isfile(clean_path):
        target_chrom = first_chrom_from_vcf(clean_path)
    if target_chrom is None and os.path.isfile(missing_bed):
        target_chrom = first_chrom_from_bed(missing_bed)

    if target_chrom is None:
        sys.stderr.write(
            "WARNING: unable to determine target chromosome; writing all intervals.\n"
        )

    # Add dropped-indel and missing-position masks.
    read_bed_intervals(dropped_bed, by_chrom)
    read_bed_intervals(missing_bed, by_chrom)

    # Subtract invariant and clean positions to prevent overlap.
    subtract_by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    if os.path.isfile(inv_path):
        read_vcf_intervals(inv_path, subtract_by_chrom)
    if os.path.isfile(clean_path):
        read_vcf_intervals(clean_path, subtract_by_chrom)

    subtract_by_chrom = filter_intervals(subtract_by_chrom, target_chrom)
    for chrom in subtract_by_chrom:
        intervals = subtract_by_chrom[chrom]
        intervals.sort(key=lambda x: (x[0], x[1]))
        subtract_by_chrom[chrom] = merge_intervals(intervals)

    by_chrom = filter_intervals(by_chrom, target_chrom)
    for chrom in by_chrom:
        intervals = by_chrom[chrom]
        if not intervals:
            continue
        intervals.sort(key=lambda x: (x[0], x[1]))
        sub = subtract_by_chrom.get(chrom, [])
        if sub:
            by_chrom[chrom] = subtract_intervals(intervals, sub)

    # Default behavior is to sort + merge unless --no-merge is given.
    with open(out_path, "wt", encoding="utf-8") as fout:
        for chrom in sorted(by_chrom.keys()):
            intervals = by_chrom[chrom]
            intervals.sort(key=lambda x: (x[0], x[1]))
            output_intervals = merge_intervals(intervals)
            for s, e in output_intervals:
                fout.write(f"{chrom}\t{s}\t{e}\n")


def first_chrom_from_vcf(path: str) -> str | None:
    with open_maybe_gzip(path, "rt") as fin:
        for raw in fin:
            if not raw or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) >= 1:
                return cols[0]
    return None


def first_chrom_from_bed(path: str) -> str | None:
    with open(path, "rt", encoding="utf-8") as fin:
        for raw in fin:
            if not raw or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) >= 1:
                return cols[0]
    return None


def filter_intervals(
    by_chrom: Dict[str, List[Tuple[int, int]]], target: str | None
) -> Dict[str, List[Tuple[int, int]]]:
    if target is None:
        return by_chrom
    filtered: Dict[str, List[Tuple[int, int]]] = {}
    if target in by_chrom:
        filtered[target] = by_chrom[target]
    return filtered


if __name__ == "__main__":
    main()
