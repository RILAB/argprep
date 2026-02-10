#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np


def _open_text(path: Path):
    # Handle optional gzip input for split outputs.
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _iter_vcf_pos(path: Path):
    # Yield 1-based positions for each non-header VCF record.
    with _open_text(path) as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            try:
                yield int(fields[1])
            except ValueError:
                continue


def _contig_length_from_fai(fai_path: Path, contig: str) -> int:
    # Look up contig length from reference index.
    with fai_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == contig:
                return int(fields[1])
    raise ValueError(f"Contig {contig!r} not found in {fai_path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build a boolean accessibility array for a contig using clean + inv VCFs."
        )
    )
    parser.add_argument("--clean", required=True, type=Path, help="Clean VCF path")
    parser.add_argument("--inv", required=True, type=Path, help="Invariant VCF path")
    parser.add_argument("--fai", required=True, type=Path, help="Reference FASTA .fai")
    parser.add_argument("--contig", required=True, help="Contig name")
    parser.add_argument("--output", required=True, type=Path, help="Output .npz path")
    args = parser.parse_args()

    contig_len = _contig_length_from_fai(args.fai, args.contig)
    # Boolean mask indexed by 0-based position across the contig.
    mask = np.zeros(contig_len, dtype=bool)

    # Mark any position present in either clean or inv as accessible.
    for pos in _iter_vcf_pos(args.clean):
        if 1 <= pos <= contig_len:
            mask[pos - 1] = True
    for pos in _iter_vcf_pos(args.inv):
        if 1 <= pos <= contig_len:
            mask[pos - 1] = True

    # Persist as compressed numpy archive with a single "mask" array.
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, mask=mask)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
