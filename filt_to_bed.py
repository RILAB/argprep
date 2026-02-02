#!/usr/bin/env python3
"""
filtered_vcf_to_bed.py

Read a *.filtered VCF (or VCF-like) file and write a BED file containing the
basepair positions of each record.

- Skips VCF header lines beginning with '#'
- Uses CHROM (col 1) and POS (col 2)
- Outputs BED intervals that represent a single bp:
    start = POS-1
    end   = POS
  (BED is 0-based, half-open; a single base at POS becomes [POS-1, POS))

- Optionally sorts and merges overlapping/adjacent intervals per chromosome.

Example:
  python filtered_vcf_to_bed.py input.filtered --out input.filtered.bed
  python filtered_vcf_to_bed.py input.filtered --merge
"""

from __future__ import annotations

import argparse
import gzip
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


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert a .filtered VCF to a BED of bp positions.")
    ap.add_argument("filtered", help="Input .filtered file (optionally .gz)")
    ap.add_argument("--out", default=None, help="Output BED path (default: <input>.bed)")
    ap.add_argument(
        "--merge",
        action="store_true",
        help="Sort and merge overlapping/adjacent intervals per chromosome.",
    )
    args = ap.parse_args()

    in_path = args.filtered
    out_path = args.out if args.out else (in_path + ".bed")

    # If merging, collect intervals by chromosome
    by_chrom: Dict[str, List[Tuple[int, int]]] = {}

    with open_maybe_gzip(in_path, "rt") as fin:
        if args.merge:
            for raw in fin:
                if not raw or raw.startswith("#"):
                    continue
                line = raw.rstrip("\n")
                cols = line.split("\t")
                if len(cols) < 2:
                    continue
                chrom = cols[0]
                try:
                    pos = int(cols[1])
                except ValueError:
                    continue
                start = pos - 1
                end = pos
                by_chrom.setdefault(chrom, []).append((start, end))
        else:
            with open(out_path, "wt", encoding="utf-8") as fout:
                for raw in fin:
                    if not raw or raw.startswith("#"):
                        continue
                    line = raw.rstrip("\n")
                    cols = line.split("\t")
                    if len(cols) < 2:
                        continue
                    chrom = cols[0]
                    try:
                        pos = int(cols[1])
                    except ValueError:
                        continue
                    start = pos - 1
                    end = pos
                    fout.write(f"{chrom}\t{start}\t{end}\n")
            return

    # If we got here, --merge was enabled: sort + merge then write
    with open(out_path, "wt", encoding="utf-8") as fout:
        for chrom in sorted(by_chrom.keys()):
            intervals = by_chrom[chrom]
            intervals.sort(key=lambda x: (x[0], x[1]))
            merged = merge_intervals(intervals)
            for s, e in merged:
                fout.write(f"{chrom}\t{s}\t{e}\n")


if __name__ == "__main__":
    main()

