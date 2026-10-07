#!/usr/bin/env python3
"""
Directly project per-sample MAF alignments onto reference coordinates and emit:
  - an all-sites VCF
  - a variant-only VCF
  - a BED mask of excluded positions
  - a per-contig QC summary
"""

from __future__ import annotations

import argparse
import bisect
import contextlib
import mmap
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Iterable

import numpy as np

try:
    from scripts.common import (  # noqa: F401  (MafRecord/iter_maf_blocks re-exported)
        MafRecord,
        iter_maf_blocks,
        merge_intervals,
        normalize_contig,
        open_text,
    )
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from scripts.common import (  # noqa: F401  (MafRecord/iter_maf_blocks re-exported)
        MafRecord,
        iter_maf_blocks,
        merge_intervals,
        normalize_contig,
        open_text,
    )


NUC_TO_CODE = {
    "A": 1,
    "C": 2,
    "G": 3,
    "T": 4,
    "-": 5,
    "?": 7,
}
CODE_TO_BASE = {
    1: "A",
    2: "C",
    3: "G",
    4: "T",
}
VALID_BASES = {"A", "C", "G", "T"}
MISSING_CODES = {0, 5, 7}


class QualityMask:
    """Low-quality intervals per source sequence, in that sequence's own
    forward-strand coordinates.

    A per-sample BED file lists ``chrom start end score`` rows (score is a
    0-1 assembly-quality value). Any base whose score is *below* the configured
    threshold is considered low quality; aligned bases at those positions are
    treated as missing during parsing. Intervals with a passing score do not
    need to be listed.
    """

    def __init__(self, low_quality: dict[str, list[tuple[int, int]]]):
        # Merged, sorted low-quality intervals per source sequence, plus a
        # parallel list of interval starts for bisecting.
        self._intervals = low_quality
        self._starts = {
            src: [start for start, _ in intervals]
            for src, intervals in low_quality.items()
        }

    def is_low(self, src: str, pos: int) -> bool:
        intervals = self._intervals.get(src)
        if not intervals:
            return False
        starts = self._starts[src]
        i = bisect.bisect_right(starts, pos) - 1
        if i < 0:
            return False
        start, end = intervals[i]
        return start <= pos < end


def load_quality_mask(bed_path: Path, quality_min: float) -> QualityMask:
    """Build a :class:`QualityMask` from a per-sample quality BED file.

    Rows are ``chrom start end score`` (whitespace-separated); ``track``/
    ``browser``/comment lines are ignored. Malformed data rows raise an error.
    Only intervals whose score is below ``quality_min`` are retained.
    """
    raw: dict[str, list[tuple[int, int]]] = {}
    with open_text(bed_path, "rt", errors="ignore") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "track", "browser")):
                continue
            parts = stripped.split()
            if len(parts) < 4:
                raise ValueError(
                    f"Malformed quality BED row in {bed_path} at line {line_number}: "
                    "expected at least 4 fields"
                )
            try:
                start = int(parts[1])
                end = int(parts[2])
                score = float(parts[3])
            except ValueError as exc:
                raise ValueError(
                    f"Malformed quality BED row in {bed_path} at line {line_number}: "
                    "start and end must be integers and score must be numeric"
                ) from exc
            if score < quality_min:
                raw.setdefault(parts[0], []).append((start, end))
    merged = {src: merge_intervals(intervals) for src, intervals in raw.items()}
    return QualityMask(merged)


def quality_bed_for_sample(quality_bed_dir: Path, sample: str) -> Path | None:
    plain = quality_bed_dir / f"{sample}.bed"
    gz = quality_bed_dir / f"{sample}.bed.gz"
    if plain.exists():
        return plain
    if gz.exists():
        return gz
    return None


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--maf-dir", required=True, help="Directory containing per-sample .maf/.maf.gz files")
    ap.add_argument("--reference-fasta", required=True, help="Reference FASTA path")
    ap.add_argument("--contig", required=True, help="Reference contig to process")
    ap.add_argument("--out-prefix", required=True, help="Output prefix for this contig")
    ap.add_argument(
        "--window-bp",
        type=int,
        default=100000,
        help="Window size for compact report statistics (default: 100000)",
    )
    ap.add_argument(
        "--maf-paths",
        nargs="*",
        default=None,
        metavar="SAMPLE=PATH",
        help="Explicit per-sample MAF paths; overrides --maf-dir for listed samples.",
    )
    ap.add_argument(
        "--samples",
        nargs="*",
        default=None,
        help="Explicit sample names; defaults to inferring *.maf and *.maf.gz stems from --maf-dir",
    )
    ap.add_argument("--max-missing-count", type=int, default=None)
    ap.add_argument("--max-missing-fraction", type=float, default=None)
    ap.add_argument("--allow-multiallelic-snps", dest="allow_multiallelic_snps", action="store_true", default=True)
    ap.add_argument("--mask-multiallelic-snps", dest="allow_multiallelic_snps", action="store_false")
    ap.add_argument(
        "--mask-indel-adjacent-snps",
        dest="mask_indel_adjacent_snps",
        action="store_true",
        default=False,
    )
    ap.add_argument(
        "--keep-indel-adjacent-snps",
        dest="mask_indel_adjacent_snps",
        action="store_false",
    )
    ap.add_argument("--add-ref", action="store_true", default=False)
    ap.add_argument(
        "--emit-argweaver-sites",
        dest="emit_argweaver_sites",
        action="store_true",
        default=False,
        help="Also emit an ARGweaver .sites file (variant sites only) alongside the VCFs.",
    )
    ap.add_argument(
        "--quality-bed-dir",
        default=None,
        help=(
            "Directory of per-sample quality BED files (<sample>.bed / .bed.gz) in "
            "each sample's own genome coordinates. Aligned bases whose score is "
            "below --quality-min are treated as missing. Requires --quality-min."
        ),
    )
    ap.add_argument(
        "--quality-min",
        type=float,
        default=None,
        help="Quality threshold in [0, 1]; bases scoring below this are treated as missing.",
    )
    return ap.parse_args()


def discover_samples(maf_dir: Path) -> list[str]:
    samples = set()
    for path in maf_dir.iterdir():
        if path.name.endswith(".maf.gz"):
            samples.add(path.name[: -len(".maf.gz")])
        elif path.name.endswith(".maf"):
            samples.add(path.name[: -len(".maf")])
    return sorted(samples)


def maf_path_for_sample(maf_dir: Path, sample: str) -> Path:
    plain = maf_dir / f"{sample}.maf"
    gz = maf_dir / f"{sample}.maf.gz"
    if plain.exists():
        return plain
    if gz.exists():
        return gz
    raise FileNotFoundError(f"Missing MAF for sample '{sample}' under {maf_dir}")


def parse_maf_path_map(values: list[str] | None) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for value in values or []:
        if "=" not in value:
            raise ValueError(f"Invalid --maf-paths entry {value!r}; expected SAMPLE=PATH")
        sample, path = value.split("=", 1)
        if not sample:
            raise ValueError(f"Invalid --maf-paths entry {value!r}; sample is empty")
        paths[sample] = Path(path)
    return paths


def maf_path_for_sample_with_map(
    maf_dir: Path,
    sample: str,
    maf_paths: dict[str, Path],
) -> Path:
    if sample in maf_paths:
        return maf_paths[sample]
    return maf_path_for_sample(maf_dir, sample)


def _read_fai_entry(fai_path: Path, contig: str) -> tuple[int, int, int, int] | None:
    """Return (length, offset, linebases, linewidth) for a contig from a .fai file."""
    with fai_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5 and parts[0] == contig:
                return int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])
    return None


def read_contig_sequence(reference_fasta: Path, contig: str) -> str:
    fai_path = Path(str(reference_fasta) + ".fai")
    if fai_path.exists() and not str(reference_fasta).endswith(".gz"):
        entry = _read_fai_entry(fai_path, contig)
        if entry is not None:
            length, offset, linebases, linewidth = entry
            if length == 0:
                return ""
            full_lines, remainder = divmod(length, linebases)
            total_bytes = full_lines * linewidth + remainder
            with open(reference_fasta, "rb") as fh:
                fh.seek(offset)
                raw = fh.read(total_bytes)
            seq = raw.replace(b"\r", b"").replace(b"\n", b"")
            if len(seq) != length:
                raise ValueError(
                    f"FAI seek returned {len(seq)} bases for '{contig}', expected {length}"
                )
            return seq.decode("ascii").upper()
    # Fallback: linear scan (handles .gz or missing .fai)
    seq_parts: list[str] = []
    current: str | None = None
    with open_text(reference_fasta, "rt", errors="ignore") as handle:
        for raw in handle:
            if raw.startswith(">"):
                name = raw[1:].strip().split()[0]
                current = name
                continue
            if current == contig:
                seq_parts.append(raw.strip())
    if not seq_parts:
        raise ValueError(f"Contig '{contig}' not found in {reference_fasta}")
    return "".join(seq_parts).upper()


def read_contig_length(reference_fasta: Path, contig: str) -> int:
    """Length of a reference contig, from the .fai index when one exists."""
    fai_path = Path(str(reference_fasta) + ".fai")
    if fai_path.exists() and not str(reference_fasta).endswith(".gz"):
        entry = _read_fai_entry(fai_path, contig)
        if entry is not None:
            return entry[0]
    return len(read_contig_sequence(reference_fasta, contig))


def read_contig_region(reference_fasta: Path, contig: str, start: int, end: int) -> str:
    """Reference bases over the 0-based half-open interval ``[start, end)``.

    With a .fai index this seeks straight to the region, so pulling a kilobase
    window does not cost a whole-contig read; without one it falls back to
    :func:`read_contig_sequence` and slices.
    """
    if start < 0 or end < start:
        raise ValueError(f"Invalid region [{start}, {end}) for '{contig}'")
    if start == end:
        return ""
    fai_path = Path(str(reference_fasta) + ".fai")
    if fai_path.exists() and not str(reference_fasta).endswith(".gz"):
        entry = _read_fai_entry(fai_path, contig)
        if entry is not None:
            length, offset, linebases, linewidth = entry
            if end > length:
                raise ValueError(
                    f"Region [{start}, {end}) exceeds length {length} of '{contig}'"
                )
            first = offset + (start // linebases) * linewidth + (start % linebases)
            last = offset + ((end - 1) // linebases) * linewidth + ((end - 1) % linebases)
            with open(reference_fasta, "rb") as fh:
                fh.seek(first)
                raw = fh.read(last - first + 1)
            seq = raw.replace(b"\r", b"").replace(b"\n", b"")
            if len(seq) != end - start:
                raise ValueError(
                    f"FAI seek returned {len(seq)} bases for '{contig}':[{start}, {end}), "
                    f"expected {end - start}"
                )
            return seq.decode("ascii").upper()
    seq = read_contig_sequence(reference_fasta, contig)
    if end > len(seq):
        raise ValueError(
            f"Region [{start}, {end}) exceeds length {len(seq)} of '{contig}'"
        )
    return seq[start:end]


def choose_sample_record(block: list[MafRecord], contig: str) -> tuple[MafRecord, MafRecord] | None:
    if not block:
        return None
    ref = block[0]
    if normalize_contig(ref.src) != normalize_contig(contig):
        return None
    if ref.strand != "+":
        # Callers walk reference columns from ref.start as forward-strand
        # positions; a minus-strand reference row would be placed silently at
        # the wrong coordinates. AnchorWave never writes one, so reject it.
        raise ValueError(
            f"MAF reference row for '{ref.src}' at start {ref.start} is on strand "
            f"'{ref.strand}'; only '+' reference rows are supported"
        )
    if len(block) >= 2:
        # Some pairwise MAFs use the same contig name for both the reference
        # and query rows, so the second alignment row is still the sample.
        return ref, block[1]
    return None


def _assign_code(calls: bytearray, idx: int, code: int) -> None:
    existing = calls[idx]
    if existing == 0:
        calls[idx] = code
        return
    if existing == code:
        return
    calls[idx] = NUC_TO_CODE["?"]


def load_sample_calls(
    maf_path: Path,
    contig: str,
    contig_len: int,
    quality_mask: QualityMask | None = None,
    track_indel_adjacent: bool = False,
) -> tuple[bytearray, bytearray | None]:
    calls = bytearray(contig_len)
    # Only allocate/track the per-column indel-adjacency array when the
    # downstream masking feature is actually enabled (see main); when off it
    # is never read, so building it wastes O(contig_len) memory and time.
    adjacent_indel_flags = bytearray(contig_len) if track_indel_adjacent else None

    for block in iter_maf_blocks(maf_path):
        chosen = choose_sample_record(block, contig)
        if chosen is None:
            continue
        ref_record, sample_record = chosen
        ref_pos = ref_record.start
        # Track the sample's own-genome coordinate so per-sample quality masks
        # (keyed by the sample sequence name and forward-strand position) can be
        # applied. `sample_offset` counts aligning sample bases consumed so far.
        sample_src = sample_record.src
        sample_start = sample_record.start
        sample_src_size = sample_record.src_size
        sample_minus = sample_record.strand == "-"
        sample_offset = 0
        prev_ref_idx: int | None = None
        mark_next_ref_adjacent = False
        for ref_char, sample_char in zip(ref_record.text.upper(), sample_record.text.upper()):
            if ref_char == "-":
                if sample_char != "-":
                    if prev_ref_idx is not None:
                        if adjacent_indel_flags is not None:
                            adjacent_indel_flags[prev_ref_idx] = 1
                    mark_next_ref_adjacent = True
                    sample_offset += 1
                continue

            if ref_pos >= contig_len:
                break
            idx = ref_pos
            ref_pos += 1
            prev_ref_idx = idx

            if mark_next_ref_adjacent:
                if adjacent_indel_flags is not None:
                    adjacent_indel_flags[idx] = 1
                mark_next_ref_adjacent = False

            if sample_char == "-":
                _assign_code(calls, idx, NUC_TO_CODE["-"])
                if adjacent_indel_flags is not None:
                    adjacent_indel_flags[idx] = 1
                    if idx > 0:
                        adjacent_indel_flags[idx - 1] = 1
                mark_next_ref_adjacent = True
                continue

            # A sample base is present here, consuming one sample-genome
            # position. Resolve its forward-strand coordinate before advancing.
            if sample_minus:
                sample_coord = sample_src_size - sample_start - 1 - sample_offset
            else:
                sample_coord = sample_start + sample_offset
            sample_offset += 1

            if quality_mask is not None and quality_mask.is_low(sample_src, sample_coord):
                _assign_code(calls, idx, NUC_TO_CODE["?"])
                continue

            if sample_char in VALID_BASES:
                _assign_code(calls, idx, NUC_TO_CODE[sample_char])
                continue

            _assign_code(calls, idx, NUC_TO_CODE["?"])

    return calls, adjacent_indel_flags


def missing_threshold(sample_count: int, max_missing_count: int | None, max_missing_fraction: float | None) -> int:
    if max_missing_count is None and max_missing_fraction is None:
        return 0
    thresholds: list[int] = []
    if max_missing_count is not None:
        thresholds.append(max(0, max_missing_count))
    if max_missing_fraction is not None:
        if not 0 <= max_missing_fraction <= 1:
            raise ValueError("--max-missing-fraction must be between 0 and 1")
        thresholds.append(int(sample_count * max_missing_fraction))
    return min(thresholds)


def format_gt(call_code: int, alt_order: list[str]) -> str:
    if call_code in MISSING_CODES:
        return "."
    base = CODE_TO_BASE[call_code]
    if base not in alt_order:
        return "0"
    return str(alt_order.index(base) + 1)


def sites_allele_string(call_codes: list[int], ref_base: str, add_ref: bool) -> str:
    """One real base per (pseudo-haploid) sample for an ARGweaver .sites line.

    Non-missing codes are always 1-4 (ACGT); missing codes (0/5/7) become ``N``.
    When ``add_ref`` is set, the synthetic REF haplotype is appended last and
    always carries the reference base, mirroring the trailing ``0`` genotype the
    VCF writes for the REF sample.
    """
    chars = ["N" if code in MISSING_CODES else CODE_TO_BASE[code] for code in call_codes]
    if add_ref:
        chars.append(ref_base)
    return "".join(chars)


def vcf_header(contig: str, contig_len: int, samples: list[str]) -> list[str]:
    return [
        "##fileformat=VCFv4.2",
        "##source=argprep.maf_to_sites",
        f"##contig=<ID={contig},length={contig_len}>",
        '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of non-missing samples">',
        '##INFO=<ID=MS,Number=1,Type=Integer,Description="Number of missing samples">',
        '##INFO=<ID=SC,Number=1,Type=String,Description="Site class: invariant or variant">',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]


def write_lines(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        for line in lines:
            handle.write(line)
            handle.write("\n")


class IntervalBuilder:
    def __init__(self, contig: str):
        self.contig = contig
        self.intervals: list[tuple[int, int]] = []

    def add(self, pos: int) -> None:
        if self.intervals and self.intervals[-1][1] == pos:
            start, _end = self.intervals[-1]
            self.intervals[-1] = (start, pos + 1)
        else:
            self.intervals.append((pos, pos + 1))

    def lines(self, sample: str | None = None) -> list[str]:
        if sample is None:
            return [
                f"{self.contig}\t{start}\t{end}"
                for start, end in self.intervals
            ]
        return [
            f"{self.contig}\t{start}\t{end}\t{sample}"
            for start, end in self.intervals
        ]


def main() -> None:
    args = parse_args()
    maf_dir = Path(args.maf_dir)
    maf_paths = parse_maf_path_map(args.maf_paths)
    reference_fasta = Path(args.reference_fasta)
    contig = args.contig
    out_prefix = Path(args.out_prefix)
    if args.window_bp <= 0:
        raise ValueError("--window-bp must be a positive integer")
    samples = list(args.samples) if args.samples else discover_samples(maf_dir)
    if not samples:
        raise ValueError(f"No samples found under {maf_dir}")
    output_samples = [*samples, "REF"] if args.add_ref else samples

    quality_bed_dir: Path | None = None
    if args.quality_bed_dir is not None:
        if args.quality_min is None:
            raise ValueError("--quality-bed-dir requires --quality-min")
        if not 0 <= args.quality_min <= 1:
            raise ValueError("--quality-min must be between 0 and 1")
        quality_bed_dir = Path(args.quality_bed_dir)
    elif args.quality_min is not None:
        raise ValueError("--quality-min requires --quality-bed-dir")

    contig_seq = read_contig_sequence(reference_fasta, contig)
    contig_len = len(contig_seq)
    if contig_len == 0:
        raise ValueError(f"Reference contig '{contig}' has length 0")
    allowed_missing = missing_threshold(
        len(samples),
        args.max_missing_count,
        args.max_missing_fraction,
    )
    with tempfile.TemporaryDirectory(prefix="argprep-calls-") as call_tmp_dir:
        sample_arrays: list[mmap.mmap] = []
        sample_array_files = []
        track_indel_adjacent = args.mask_indel_adjacent_snps
        # Only allocated/merged when the indel-adjacent masking feature is on;
        # otherwise it is never read in the main loop.
        adjacent_indel_flags = bytearray(contig_len) if track_indel_adjacent else None
        merged_adjacent = (
            np.frombuffer(adjacent_indel_flags, dtype=np.uint8)
            if adjacent_indel_flags is not None
            else None
        )

        for sample in samples:
            quality_mask: QualityMask | None = None
            if quality_bed_dir is not None:
                bed_path = quality_bed_for_sample(quality_bed_dir, sample)
                if bed_path is not None:
                    quality_mask = load_quality_mask(bed_path, args.quality_min)
            calls, sample_adjacent_indels = load_sample_calls(
                maf_path_for_sample_with_map(maf_dir, sample, maf_paths),
                contig,
                contig_len,
                quality_mask,
                track_indel_adjacent=track_indel_adjacent,
            )
            call_file = tempfile.TemporaryFile(dir=call_tmp_dir)
            call_file.write(calls)
            call_file.flush()
            sample_arrays.append(
                mmap.mmap(call_file.fileno(), contig_len, access=mmap.ACCESS_READ)
            )
            sample_array_files.append(call_file)
            del calls
            if merged_adjacent is not None and sample_adjacent_indels is not None:
                # Vectorized OR-merge over a zero-copy uint8 view of the flags
                # buffer (replaces a pure-Python O(contig_len) loop per sample).
                merged_adjacent |= np.frombuffer(sample_adjacent_indels, dtype=np.uint8)

        all_sites_path = out_prefix.with_suffix(out_prefix.suffix + ".all_sites.vcf")
        variants_path = out_prefix.with_suffix(out_prefix.suffix + ".vcf")
        mask_path = out_prefix.with_suffix(out_prefix.suffix + ".mask.bed")
        summary_path = out_prefix.with_suffix(out_prefix.suffix + ".site_summary.tsv")
        report_stats_path = out_prefix.with_suffix(out_prefix.suffix + ".report_stats.tsv")
        sites_path = out_prefix.with_suffix(out_prefix.suffix + ".sites")

        all_sites_path.parent.mkdir(parents=True, exist_ok=True)
        with contextlib.ExitStack() as stack:
            all_sites = stack.enter_context(open_text(all_sites_path, "wt"))
            variants = stack.enter_context(open_text(variants_path, "wt"))
            sites = (
                stack.enter_context(open_text(sites_path, "wt"))
                if args.emit_argweaver_sites
                else None
            )
            for line in vcf_header(contig, contig_len, output_samples):
                all_sites.write(f"{line}\n")
                variants.write(f"{line}\n")
            if sites is not None:
                sites.write("NAMES\t" + "\t".join(output_samples) + "\n")
                sites.write(f"REGION\t{contig}\t1\t{contig_len}\n")

            mask_intervals = IntervalBuilder(contig)
            masked_total = 0
            retained_total = 0
            counts: Counter[str] = Counter()
            window_count = (contig_len + args.window_bp - 1) // args.window_bp
            invariant_by_window = np.zeros(window_count, dtype=np.int64)
            variant_by_window = np.zeros(window_count, dtype=np.int64)
            masked_by_window = np.zeros(window_count, dtype=np.int64)
            sample_called = np.zeros(len(samples), dtype=np.int64)
            sample_carried_variant = np.zeros(len(samples), dtype=np.int64)

            allow_multiallelic = args.allow_multiallelic_snps
            add_ref = args.add_ref
            # Zero-copy uint8 views over each sample's mmap-backed call buffer.
            sample_views = [np.frombuffer(calls, dtype=np.uint8) for calls in sample_arrays]
            num_samples = len(sample_views)
            # Codes 1-4 map to A, C, G, T (see CODE_TO_BASE); building the base
            # list in code order yields the same order as sorted(allele_set).
            code_bases = ("A", "C", "G", "T")
            # Process in position-chunks so the stacked per-sample array is
            # bounded to num_samples * CHUNK bytes rather than num_samples * L.
            CHUNK = 1_000_000

            for c0 in range(0, contig_len, CHUNK):
                c1 = min(c0 + CHUNK, contig_len)
                block = np.stack([view[c0:c1] for view in sample_views])
                missing_mask = (block == 0) | (block == 5) | (block == 7)
                missing_counts = missing_mask.sum(axis=0).tolist()
                has_unaligned = (block == 0).any(axis=0).tolist()
                present1 = (block == 1).any(axis=0).tolist()
                present2 = (block == 2).any(axis=0).tolist()
                present3 = (block == 3).any(axis=0).tolist()
                present4 = (block == 4).any(axis=0).tolist()

                for j in range(c1 - c0):
                    idx = c0 + j
                    ref_base = contig_seq[idx]
                    pos = idx + 1
                    if ref_base not in VALID_BASES:
                        mask_intervals.add(idx)
                        masked_total += 1
                        masked_by_window[idx // args.window_bp] += 1
                        counts["masked_ref_non_acgt"] += 1
                        continue

                    missing = missing_counts[j]
                    ns = num_samples - missing

                    if ns == 0:
                        mask_intervals.add(idx)
                        masked_total += 1
                        masked_by_window[idx // args.window_bp] += 1
                        if has_unaligned[j]:
                            counts["masked_no_alignment"] += 1
                        else:
                            counts["masked_missingness"] += 1
                        continue

                    present = (present1[j], present2[j], present3[j], present4[j])
                    allele_set = {code_bases[b] for b in range(4) if present[b]}
                    alt_order = [a for a in code_bases if a in allele_set and a != ref_base]

                    if track_indel_adjacent and alt_order and adjacent_indel_flags[idx]:
                        mask_intervals.add(idx)
                        masked_total += 1
                        masked_by_window[idx // args.window_bp] += 1
                        counts["masked_indel_adjacent"] += 1
                        continue

                    if missing > allowed_missing:
                        mask_intervals.add(idx)
                        masked_total += 1
                        masked_by_window[idx // args.window_bp] += 1
                        if has_unaligned[j]:
                            counts["masked_no_alignment"] += 1
                        else:
                            counts["masked_missingness"] += 1
                        continue

                    if len(alt_order) > 1 and not allow_multiallelic:
                        mask_intervals.add(idx)
                        masked_total += 1
                        masked_by_window[idx // args.window_bp] += 1
                        counts["masked_multiallelic"] += 1
                        continue

                    info = f"NS={ns};MS={missing}"
                    call_codes = block[:, j].tolist()
                    if not alt_order and allele_set == {ref_base}:
                        counts["all_sites"] += 1
                        counts["invariant"] += 1
                        retained_total += 1
                        invariant_by_window[idx // args.window_bp] += 1
                        sample_called += ~missing_mask[:, j]
                        genotypes = [format_gt(code, alt_order) for code in call_codes]
                        if add_ref:
                            genotypes.append("0")
                        record = (
                            f"{contig}\t{pos}\t.\t{ref_base}\t.\t.\tPASS\t{info};SC=invariant\tGT\t"
                            + "\t".join(genotypes)
                        )
                        all_sites.write(f"{record}\n")
                        continue

                    counts["all_sites"] += 1
                    counts["variants"] += 1
                    retained_total += 1
                    variant_by_window[idx // args.window_bp] += 1
                    called_at_site = ~missing_mask[:, j]
                    sample_called += called_at_site
                    sample_carried_variant += called_at_site & (
                        block[:, j] != NUC_TO_CODE[ref_base]
                    )
                    genotypes = [format_gt(code, alt_order) for code in call_codes]
                    if add_ref:
                        genotypes.append("0")
                    alt_field = ",".join(alt_order)
                    record = (
                        f"{contig}\t{pos}\t.\t{ref_base}\t{alt_field}\t.\tPASS\t{info};SC=variant\tGT\t"
                        + "\t".join(genotypes)
                    )
                    all_sites.write(f"{record}\n")
                    variants.write(f"{record}\n")
                    if sites is not None:
                        sites.write(
                            f"{pos}\t{sites_allele_string(call_codes, ref_base, add_ref)}\n"
                        )

        if retained_total + masked_total != contig_len:
            raise ValueError(
                f"coverage mismatch for {contig}: union={retained_total + masked_total}, chrom_len={contig_len}"
            )

        bed_lines = mask_intervals.lines()
        write_lines(mask_path, bed_lines)

        for sample, calls in zip(samples, sample_arrays):
            # Vectorized replacement for a pure-Python O(contig_len) scan: derive
            # missing-run intervals from a zero-copy uint8 view of the mmap.
            arr = np.frombuffer(calls, dtype=np.uint8)
            missing = (arr == 0) | (arr == 5) | (arr == 7)
            # Run boundaries: pad with False on both ends, then diff. +1 marks a
            # run start, -1 marks the position just past a run end. This yields
            # the same merged [start, end) intervals as IntervalBuilder.add.
            edges = np.diff(np.concatenate(([np.uint8(0)], missing.view(np.uint8), [np.uint8(0)])).astype(np.int8))
            starts = np.flatnonzero(edges == 1)
            ends = np.flatnonzero(edges == -1)
            sample_lines = [
                f"{contig}\t{int(start)}\t{int(end)}\t{sample}"
                for start, end in zip(starts, ends)
            ]
            write_lines(Path(str(out_prefix) + f".{sample}.missing.bed"), sample_lines)

        summary_lines = [
            "metric\tvalue",
            f"contig\t{contig}",
            f"contig_length\t{contig_len}",
            f"samples\t{len(samples)}",
            f"allowed_missing\t{allowed_missing}",
        ]
        for key in (
            "all_sites",
            "variants",
            "invariant",
            "masked_missingness",
            "masked_indel_adjacent",
            "masked_multiallelic",
            "masked_no_alignment",
            "masked_ref_non_acgt",
        ):
            summary_lines.append(f"{key}\t{counts.get(key, 0)}")
        summary_lines.append(f"masked_total\t{masked_total}")
        summary_lines.append(f"masked_intervals\t{len(bed_lines)}")
        write_lines(summary_path, summary_lines)

        report_stats_lines = [
            "record_type\tcontig\tstart\tend\tsample\tinvariant\tvariant\tmasked\tcalled\tcarried_variant"
        ]
        for window_idx in range(window_count):
            start = window_idx * args.window_bp
            end = min(start + args.window_bp, contig_len)
            report_stats_lines.append(
                f"window\t{contig}\t{start}\t{end}\t\t"
                f"{int(invariant_by_window[window_idx])}\t"
                f"{int(variant_by_window[window_idx])}\t"
                f"{int(masked_by_window[window_idx])}\t\t"
            )
        for sample_idx, sample in enumerate(samples):
            report_stats_lines.append(
                f"sample\t{contig}\t\t\t{sample}\t\t\t\t"
                f"{int(sample_called[sample_idx])}\t"
                f"{int(sample_carried_variant[sample_idx])}"
            )
        write_lines(report_stats_path, report_stats_lines)

        # Drop every numpy view derived from the mmaps: a live frombuffer view
        # keeps an exported pointer that would make mmap.close() raise
        # BufferError.
        sample_views.clear()
        arr = missing = edges = None
        for calls in sample_arrays:
            calls.close()
        for call_file in sample_array_files:
            call_file.close()


if __name__ == "__main__":
    main()
