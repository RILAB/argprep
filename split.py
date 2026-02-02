#!/usr/bin/env python3
"""
VCF splitter / filter with END-span expansion

This script reads a VCF (optionally gzipped) and writes each record
to EXACTLY ONE of three output files:

  <prefix>.inv
  <prefix>.filtered
  <prefix>.clean

Key feature:
- Any record with INFO containing END= is EXPANDED so that the record
  is replicated for every base position from POS to END (inclusive).

âš ï¸ WARNING:
Expanding END spans can massively increase file size if END-POS is large.
"""

from __future__ import annotations
import argparse
import gzip
from pathlib import Path
from typing import TextIO
import sys
import datetime
import getpass
import time
import platform
import os
from tqdm import tqdm

SCRIPT_NAME = "vcf-splitter"
SCRIPT_VERSION = "1.0.0"

# Valid single-base SNP alleles
VALID_BASES = {"A", "C", "G", "T"}

# ------------------------------------------------------------
# UTILITY: modify header to track provenance of split files!
# ------------------------------------------------------------
def build_provenance_headers(
    in_path: str,
    out_inv: str,
    out_filt: str,
    out_clean: str,
    depth: int,
    gzip_output: bool,
) -> list[str]:
    """
    Return a list of VCF meta-header lines (starting with '##') that describe
    how this file was produced.
    These lines must be written BEFORE the final '#CHROM\tPOS\t...' header line.
    """
    ts = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    user = getpass.getuser()
    host = platform.node()

    cmdline = " ".join([os.path.basename(sys.argv[0])] + sys.argv[1:])

    lines = [
        f"##source={SCRIPT_NAME} v{SCRIPT_VERSION}",
        f"##run.user={user}",
        f"##run.host={host}",
        f"##run.timestamp={ts}",
        f"##run.commandline={cmdline}",
        f"##input.file={in_path}",
        f"##output.inv={out_inv}",
        f"##output.filtered={out_filt}",
        f"##output.clean={out_clean}",
        f"##parameters.depth_threshold={depth}",
        f"##parameters.gzip_output={'true' if gzip_output else 'false'}",
    ]
    lines.append("##notes=Records with valid INFO/END are expanded across POS..END and routed to .inv; other filters apply as documented.")
    return [ln + "\n" for ln in lines]



# ------------------------------------------------------------
# Utility: open plain-text or gzipped files transparently
# ------------------------------------------------------------
def open_maybe_gzip(path: str, mode: str) -> TextIO:
    """
    Open a file that may or may not be gzipped.

    Parameters
    ----------
    path : str
        Path to file (.vcf or .vcf.gz)
    mode : str
        Text mode: "rt" for reading, "wt" for writing

    Returns
    -------
    TextIO
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode)  # type: ignore
    return open(path, mode, encoding="utf-8")


# ------------------------------------------------------------
# INFO field parsers
# ------------------------------------------------------------
def extract_dp(info: str) -> int | None:
    """
    Extract integer DP value from the INFO column.

    Returns
    -------
    int : DP value if present and parseable
    None : if INFO=".", DP missing, or malformed
    """
    if info == ".":
        return None

    for field in info.split(";"):
        if field.startswith("DP="):
            try:
                return int(field.split("=", 1)[1])
            except ValueError:
                return None

    return None


def extract_end(info: str) -> int | None:
    """
    Extract integer END value from the INFO column.

    Returns
    -------
    int : END value if present and parseable
    None : if INFO="." or END missing/malformed
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


# ------------------------------------------------------------
# Main driver
# ------------------------------------------------------------
def main() -> None:
    """
    Main entry point.
    Handles argument parsing, file routing, and record classification.
    """
    inv_bp = 0
    filtered_bp = 0
    clean_bp = 0


    # ---------------- Argument parsing ----------------
    ap = argparse.ArgumentParser(
        description="Split/filter a VCF into inv/filtered/clean, expanding END spans."
    )
    ap.add_argument(
        "--filter-multiallelic",
        action="store_true",
        help="Filter sites with multiple A/C/G/T alleles (e.g. ALT=T,G)"
    )
    ap.add_argument(
        "vcf",
        help="Input VCF (.vcf or .vcf.gz)"
    )
    ap.add_argument(
        "--depth",
        type=int,
        required=True,
        help="Depth threshold: DP < depth => filtered"
    )
    ap.add_argument(
        "--out-prefix",
        default=None,
        help="Output prefix (default: input filename without .vcf/.vcf.gz)"
    )
    ap.add_argument(
        "--gzip-output",
        action="store_true",
        help="Gzip all output files (.gz)"
    )
    args = ap.parse_args()

    in_path = args.vcf
    depth = args.depth

    # ---------------- Output prefix logic ----------------
    if args.out_prefix:
        prefix = args.out_prefix
    else:
        # Strip .vcf or .vcf.gz from input name
        p = Path(in_path)
        name = p.name
        if name.endswith(".vcf.gz"):
            name = name[:-7]
        elif name.endswith(".vcf"):
            name = name[:-4]
        prefix = str(p.with_name(name))

    out_inv = prefix + ".inv"
    out_filt = prefix + ".filtered"
    out_clean = prefix + ".clean"

    # ---------------- Open files ----------------
    def open_out(path, gzip_output):
        if gzip_output:
            return gzip.open(path + ".gz", "wt")
        return open(path, "wt", encoding="utf-8")

    
    with open_maybe_gzip(in_path, "rt") as fin, \
         open_out(out_inv, args.gzip_output) as f_inv, \
         open_out(out_filt, args.gzip_output) as f_filt, \
         open_out(out_clean, args.gzip_output) as f_clean:
    
        # --- NEW: buffer header lines until we see '#CHROM' ---
        header_buffer: list[str] = []
        headers_written = False
    
        # Prebuild provenance lines (to be inserted before '#CHROM')
        prov_lines = build_provenance_headers(
            in_path=in_path,
            out_inv=out_inv + (".gz" if args.gzip_output else ""),
            out_filt=out_filt + (".gz" if args.gzip_output else ""),
            out_clean=out_clean + (".gz" if args.gzip_output else ""),
            depth=depth,
            gzip_output=args.gzip_output,
        )
    
        record_count = 0
        last_report_time = time.time()
        for raw in fin:
            record_count += 1
            if record_count % 100_000 == 0:
                fields = raw.rstrip("\n").split("\t")
                if len(fields) >= 2:
                    chrom, pos = fields[0], fields[1]
                else:
                    chrom, pos = "?", "?"

                now = time.time()
                rate = record_count / max(now - last_report_time, 1e-6)

                sys.stderr.write(
                    f"\r[progress] {record_count:,} records  @ chrom {chrom} Mb {int(pos)/1E6}  "
                    f"({rate:,.0f} rec/s)"
                )
                sys.stderr.flush()

            # ---------------- Header handling ----------------
            if raw.startswith("#"):
                #print file format
                if raw.startswith("##fileformat"):
                    f_inv.write(raw)
                    f_filt.write(raw)
                    f_clean.write(raw)
                # Buffer meta/header lines until we see the '#CHROM ... INFO' line.
                elif raw.startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"):
                    # first: inject provenance headers (##...)
                    for ln in prov_lines:
                        f_inv.write(ln)
                        f_filt.write(ln)
                        f_clean.write(ln)
                    
                    # 2nd: write buffered input meta-headers (##...)
                    for h in header_buffer:
                        f_inv.write(h)
                        f_filt.write(h)
                        f_clean.write(h)

                    

                    # Third: write the #CHROM header line itself
                    f_inv.write(raw)
                    f_filt.write(raw)
                    f_clean.write(raw)
                    headers_written = True
                else:
                    header_buffer.append(raw)
                continue

            # If the input was missing '#CHROM' (malformed), we never wrote headers.
            # Flush buffered headers + provenance BEFORE the first data record.
            if not headers_written:
                for h in header_buffer:
                    f_inv.write(h)
                    f_filt.write(h)
                    f_clean.write(h)
                for ln in prov_lines:
                    f_inv.write(ln)
                    f_filt.write(ln)
                    f_clean.write(ln)
                headers_written = True

            # ---------------- Data line processing ----------------
            line = raw.rstrip("\n")
            cols = line.split("\t")

            # Malformed line (no INFO column) => filtered
            if len(cols) < 8:
                f_filt.write(line + "\n")
                filtered_bp += 1
                continue

            # VCF columns (0-based): 0 CHROM, 1 POS, 3 REF, 4 ALT, 7 INFO
            ref = cols[3]
            alt_field = cols[4]
            info = cols[7]

            # ============================================================
            # 1) .inv (highest priority): INFO == '.' OR contains END=
            # If END= exists, expand across POS..END inclusive.
            # ============================================================
            end_val = extract_end(info)
            if info == "." or end_val is not None:
                if end_val is None:
                    f_inv.write(line + "\n")
                    inv_bp += 1
                    continue

                try:
                    pos0 = int(cols[1])
                except ValueError:
                    f_inv.write(line + "\n")
                    inv_bp += 1
                    continue

                if end_val <= pos0:
                    f_inv.write(line + "\n")
                    inv_bp += 1
                    continue

                span = end_val - pos0 + 1
                inv_bp += span
                for pos in range(pos0, end_val + 1):
                    cols[1] = str(pos)
                    f_inv.write("\t".join(cols) + "\n")
                continue

            # ============================================================
            # 2) .filtered
            # ============================================================
            alts = [a.strip() for a in alt_field.split(",")] if alt_field != "." else []
            alts_no_nonref = [a for a in alts if a != "<NON_REF>"]

            dp = extract_dp(info)
            dp_low = (dp is not None and dp < depth)
            has_star = ("*" in alts_no_nonref)
            ref_long = (len(ref) > 1)

            valid_base_alts = [a for a in alts_no_nonref if a in VALID_BASES]
            multiple_valid_bases = (len(set(valid_base_alts)) > 1)

            has_non_acgt_nonstar = any(
                (a not in VALID_BASES) and (a != "*") for a in alts_no_nonref
            )

            if ( 
                dp_low or 
                has_star or 
                ref_long or 
                (args.filter_multiallelic and multiple_valid_bases) or
                has_non_acgt_nonstar):
                f_filt.write(line + "\n")
                filtered_bp += 1
                continue

            # ============================================================
            # 3) .clean
            # ============================================================
            f_clean.write(line + "\n")
            clean_bp += 1

        print("Output summary (non-header bp):",file=sys.stderr)
        print(f"  inv:      {inv_bp:,}",file=sys.stderr)
        print(f"  filtered: {filtered_bp:,}",file=sys.stderr)
        print(f"  clean:    {clean_bp:,}",file=sys.stderr)
    
if __name__ == "__main__":
    main()
