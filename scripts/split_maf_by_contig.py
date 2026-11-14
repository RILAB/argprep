#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    from scripts.common import MafRecord, iter_maf_blocks, normalize_contig, open_text
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from scripts.common import MafRecord, iter_maf_blocks, normalize_contig, open_text


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Split one pairwise MAF into per-reference-contig MAF chunks."
    )
    ap.add_argument("--maf", required=True, help="Input .maf or .maf.gz")
    ap.add_argument("--sample", required=True, help="Sample name for output directory")
    ap.add_argument("--out-root", required=True, help="Output root directory")
    ap.add_argument("--contigs", nargs="+", required=True, help="Reference contigs to retain")
    return ap.parse_args()


def format_record(record: MafRecord) -> str:
    return (
        f"s {record.src} {record.start} {record.size} "
        f"{record.strand} {record.src_size} {record.text}"
    )


def main() -> None:
    args = parse_args()
    maf_path = Path(args.maf)
    sample_dir = Path(args.out_root) / args.sample
    sample_dir.mkdir(parents=True, exist_ok=True)

    contigs = [str(contig) for contig in args.contigs]
    normalized_to_output: dict[str, str] = {}
    for contig in contigs:
        normalized = normalize_contig(contig)
        if normalized in normalized_to_output and normalized_to_output[normalized] != contig:
            raise ValueError(
                f"Ambiguous normalized contig {normalized!r}: "
                f"{normalized_to_output[normalized]!r} and {contig!r}"
            )
        normalized_to_output[normalized] = contig

    handles = {}
    try:
        for contig in contigs:
            out_path = sample_dir / f"{contig}.maf.gz"
            handle = open_text(out_path, "wt")
            handle.write("##maf version=1\n")
            handles[contig] = handle

        for block in iter_maf_blocks(maf_path):
            if not block:
                continue
            out_contig = normalized_to_output.get(normalize_contig(block[0].src))
            if out_contig is None:
                continue
            handle = handles[out_contig]
            handle.write("\n")
            handle.write("a score=0\n")
            for record in block:
                handle.write(format_record(record))
                handle.write("\n")
    finally:
        for handle in handles.values():
            handle.close()


if __name__ == "__main__":
    main()
