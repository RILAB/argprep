#!/usr/bin/env python3
"""
dropSV.py

Remove large indels from gVCFs to avoid GATK merge failures.
Produces cleaned gVCFs in ./cleangVCF and writes per-sample
<prefix>.dropped_indels.bed for downstream masking.
"""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Tuple


def eprint(msg: str) -> None:
    # Centralized stderr printing for errors.
    sys.stderr.write(msg + "\n")


def require_tool(name: str) -> None:
    # Ensure external dependencies are available before doing work.
    if shutil.which(name) is None:
        eprint(f"Error: Please ensure {name} is installed and available.")
        sys.exit(1)


def parse_args() -> argparse.Namespace:
    # Mirror the original shell script CLI.
    ap = argparse.ArgumentParser(
        description="Remove large indels from gVCFs for downstream GATK parsing."
    )
    ap.add_argument(
        "-c",
        "--cutoff",
        type=int,
        default=9101264,
        help="INDEL size to remove (default: 9101264)",
    )
    ap.add_argument(
        "-d",
        "--directory",
        default=".",
        help="Source directory containing .gvcf.gz files (default: current directory)",
    )
    args = ap.parse_args()
    return args


def list_gvcfs(src_dir: Path) -> List[Path]:
    # Discover input gVCFs (gzipped only, consistent with the original script).
    return sorted(src_dir.glob("*.gvcf.gz"))


def run_cmd(cmd: List[str], cwd: Path | None = None) -> subprocess.CompletedProcess:
    # Run a subprocess and capture output for error reporting.
    return subprocess.run(cmd, cwd=cwd, check=False, text=True, capture_output=True)


def read_bcftools_query(path: Path) -> Iterable[Tuple[str, int, str, str]]:
    # Stream bcftools query output to avoid loading full files into memory.
    cmd = [
        "bcftools",
        "query",
        "-f",
        "%CHROM\t%POS\t%REF\t%ALT\n",
        str(path),
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) != 4:
            continue
        chrom, pos_s, ref, alt = parts
        try:
            pos = int(pos_s)
        except ValueError:
            continue
        yield chrom, pos, ref, alt
    proc.stdout.close()
    err = proc.stderr.read() if proc.stderr else ""
    code = proc.wait()
    if code != 0:
        eprint(err.strip() or f"bcftools query failed for {path}")
        sys.exit(1)


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    # Merge overlapping or adjacent half-open intervals.
    if not intervals:
        return []
    intervals.sort()
    merged: List[Tuple[int, int]] = []
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:
            if e > cur_e:
                cur_e = e
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged


def write_bed(path: Path, intervals_by_chrom: dict[str, List[Tuple[int, int]]]) -> None:
    # Write merged intervals in BED (0-based, half-open) format.
    with open(path, "wt", encoding="utf-8") as f:
        for chrom in sorted(intervals_by_chrom.keys()):
            merged = merge_intervals(intervals_by_chrom[chrom])
            for s, e in merged:
                f.write(f"{chrom}\t{s}\t{e}\n")


def count_data_lines_gz(path: Path) -> int:
    # Count non-header records to report how many lines were removed.
    count = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            count += 1
    return count


def main() -> None:
    args = parse_args()

    # External tools used for filtering/indexing.
    require_tool("bcftools")
    require_tool("tabix")

    src_dir = Path(args.directory).resolve()
    if not src_dir.is_dir():
        eprint(f"ERROR: Failed to change to {src_dir}. Please check directory exists")
        sys.exit(1)

    # Enumerate inputs and prepare output directory.
    gvcfs = list_gvcfs(src_dir)
    print(f"{len(gvcfs)} gVCF files found in {src_dir}")

    out_dir = src_dir / "cleangVCF"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Identify "super indels" across all gVCFs.
    super_indels: List[Tuple[str, str, int, int, int, int]] = []
    for gvcf in gvcfs:
        for chrom, pos, ref, alt in read_bcftools_query(gvcf):
            d = abs(len(ref) - len(alt))
            if d > args.cutoff:
                super_indels.append(
                    (str(gvcf), chrom, pos, d, len(ref), len(alt))
                )

    super_indels.sort(key=lambda x: x[0])
    super_indels_path = out_dir / "super_indels.txt"
    with open(super_indels_path, "wt", encoding="utf-8") as f:
        for row in super_indels:
            f.write(
                f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\t{row[5]}\n"
            )

    # Track which gVCFs actually need filtering.
    need_filter = sorted({row[0] for row in super_indels})
    need_filter_path = out_dir / "need_filter.txt"
    with open(need_filter_path, "wt", encoding="utf-8") as f:
        for path in need_filter:
            f.write(path + "\n")

    # Global merged BED across all dropped indels.
    global_intervals: dict[str, List[Tuple[int, int]]] = {}
    for _, chrom, pos, _, ref_len, alt_len in super_indels:
        start = pos - 1
        end = start + ref_len if ref_len > alt_len else pos
        global_intervals.setdefault(chrom, []).append((start, end))
    write_bed(out_dir / "dropped_indels.bed", global_intervals)

    # Process each gVCF: filter if needed, otherwise copy.
    for gvcf in gvcfs:
        base = gvcf.name
        print(f"Processing {base}")
        out_gvcf = out_dir / base
        out_tbi = out_dir / (base + ".tbi")

        if str(gvcf) in need_filter:
            # Filter out large indel positions using bcftools -T ^bad_sites.
            lines_before = count_data_lines_gz(gvcf)

            bad_sites = []
            for row in super_indels:
                if row[0] == str(gvcf):
                    bad_sites.append((row[1], row[2]))
            bad_sites_path = out_dir / "bad_sites.tmp"
            with open(bad_sites_path, "wt", encoding="utf-8") as f:
                for chrom, pos in bad_sites:
                    f.write(f"{chrom}\t{pos}\n")

            cmd = [
                "bcftools",
                "view",
                "-T",
                f"^{bad_sites_path}",
                "-Oz",
                "-o",
                str(out_gvcf),
                str(gvcf),
            ]
            proc = run_cmd(cmd)
            if proc.returncode != 0:
                eprint(proc.stderr.strip() or "bcftools view failed")
                sys.exit(1)

            # Index the filtered gVCF.
            proc = run_cmd(["tabix", "-p", "vcf", str(out_gvcf)])
            if proc.returncode != 0:
                eprint(proc.stderr.strip() or "tabix failed")
                sys.exit(1)

            lines_after = count_data_lines_gz(out_gvcf)
            print(
                f"  before: {lines_before}, after: {lines_after}, "
                f"removed: {lines_before - lines_after}"
            )

        else:
            # No large indels: copy input and emit an empty BED.
            shutil.copy2(gvcf, out_gvcf)
            if gvcf.with_suffix(gvcf.suffix + ".tbi").exists():
                shutil.copy2(gvcf.with_suffix(gvcf.suffix + ".tbi"), out_tbi)
            print("  no large indels found, file copied without changes.")

    print(
        f"Total super large indels identified across all gVCFs: {len(super_indels)}"
    )
    print(f"Total gVCFs with super large indels removed: {len(need_filter)}")
    print(f"Dropped indel positions BED written to {out_dir / 'dropped_indels.bed'}")

    # Cleanup temporary files created during filtering.
    for tmp in ("bad_sites.tmp", "need_filter.txt", "super_indels.txt"):
        tmp_path = out_dir / tmp
        if tmp_path.exists():
            tmp_path.unlink()

    print(f"All done. Cleaned gVCFs are in {out_dir}")


if __name__ == "__main__":
    # Keep the original behavior: print help if called with no args.
    if len(sys.argv) == 1:
        print("Usage: dropSV.py [OPTION]...")
        print("This script removes large indels for later parsing with GATK.")
        print("")
        print("Options:")
        print("  -h, --help       Display this help message and exit")
        print("  -c  --cutoff     INDEL size to remove (default: 9101264)")
        print(
            "  -d  --directory  Full path of the source directory containing gVCF files "
            "(default: current directory)"
        )
        sys.exit(1)
    main()
