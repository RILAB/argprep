#!/usr/bin/env python3
"""Build a gapped multiple-alignment FASTA for one reference window by merging
the per-sample pairwise MAFs transitively through the reference.

Unlike ``window_to_fasta.py`` -- which reconstructs a reference-anchored
substitution view from the all-sites VCF, so every row is exactly the window
length -- this reads the MAFs directly and keeps indels: a deletion in a sample
becomes a gap column, and an insertion adds columns that every other sample
pads with gaps.

Usage:
    python scripts/maf_to_fasta.py \\
        --maf-chunk-root results/maf_by_contig \\
        --reference-fasta /path/to/reference.fa \\
        --contig <contig> --start <start> --end <end> \\
        --out window.fa

See the "Auxiliary scripts" section of the README for the homology,
conflict, and masking caveats that come with a reference-anchored merge.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    from scripts.common import normalize_contig, open_text
    from scripts.maf_to_sites import (
        NUC_TO_CODE,
        QualityMask,
        VALID_BASES,
        _assign_code,
        choose_sample_record,
        discover_samples,
        iter_maf_blocks,
        load_quality_mask,
        maf_path_for_sample_with_map,
        parse_maf_path_map,
        quality_bed_for_sample,
        read_contig_length,
        read_contig_region,
    )
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from scripts.common import normalize_contig, open_text
    from scripts.maf_to_sites import (
        NUC_TO_CODE,
        QualityMask,
        VALID_BASES,
        _assign_code,
        choose_sample_record,
        discover_samples,
        iter_maf_blocks,
        load_quality_mask,
        maf_path_for_sample_with_map,
        parse_maf_path_map,
        quality_bed_for_sample,
        read_contig_length,
        read_contig_region,
    )


# Rendering of the per-anchor call codes that maf_to_sites already assigns.
# The one place this deliberately diverges from the VCF path is code 5: there a
# deletion is just another flavor of missing, here it is the gap it actually is.
CODE_TO_CHAR = {
    0: "N",  # no alignment block covered this reference position
    1: "A",
    2: "C",
    3: "G",
    4: "T",
    5: "-",  # sample carries a deletion over this reference base
    7: "N",  # ambiguous base, below --quality-min, or conflicting blocks
}
GAP = "-"
UNKNOWN = "?"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--maf-chunk-root",
        help=(
            "Root of the workflow's per-contig MAF chunks; sample MAFs are read "
            "from <root>/<sample>/<contig>.maf.gz"
        ),
    )
    source.add_argument(
        "--maf-dir",
        help="Directory of whole-genome per-sample MAFs (<sample>.maf / .maf.gz)",
    )
    ap.add_argument("--reference-fasta", required=True, help="Reference FASTA path")
    ap.add_argument("--contig", required=True, help="Reference contig to extract")
    ap.add_argument("--start", type=int, default=1, help="1-based inclusive window start (default: 1)")
    ap.add_argument(
        "--end",
        type=int,
        default=None,
        help="1-based inclusive window end (default: end of contig)",
    )
    ap.add_argument("--out", required=True, help="Output FASTA path (.gz is compressed)")
    ap.add_argument(
        "--samples",
        nargs="*",
        default=None,
        help="Explicit sample names; defaults to discovering them under the MAF source",
    )
    ap.add_argument(
        "--maf-paths",
        nargs="*",
        default=None,
        metavar="SAMPLE=PATH",
        help="Explicit per-sample MAF paths; overrides the MAF source for listed samples.",
    )
    ap.add_argument(
        "--quality-bed-dir",
        default=None,
        help=(
            "Directory of per-sample quality BED files (<sample>.bed / .bed.gz) in "
            "each sample's own genome coordinates, as in maf_to_sites.py. Requires "
            "--quality-min."
        ),
    )
    ap.add_argument(
        "--quality-min",
        type=float,
        default=None,
        help="Quality threshold in [0, 1]; bases scoring below this render as N.",
    )
    ap.add_argument(
        "--mask-bed",
        default=None,
        help=(
            "combined.<contig>.mask.bed from a completed run. Reference positions it "
            "masks render as N for every sample and their insertions are dropped, so "
            "the FASTA agrees with the VCF. Omitted, the FASTA can show bases at "
            "positions the VCF masked."
        ),
    )
    ap.add_argument(
        "--ref-name",
        default="REF",
        help="Sequence name for the reference row (default: REF)",
    )
    ap.add_argument(
        "--no-reference",
        dest="include_reference",
        action="store_false",
        default=True,
        help="Omit the reference row from the output",
    )
    ap.add_argument(
        "--line-width",
        type=int,
        default=80,
        help="FASTA wrap width; 0 writes each sequence on one line (default: 80)",
    )
    ap.add_argument(
        "--max-columns",
        type=int,
        default=5_000_000,
        help=(
            "Refuse to build an alignment wider than this. Insertion-rich regions "
            "across many samples can make a window far wider than its reference "
            "span, and memory is columns x samples (default: 5000000)."
        ),
    )
    return ap.parse_args()


def merge_insertion(insertions: dict[int, str], slot: int, seq: str) -> None:
    """Record ``seq`` as the sample's insertion at ``slot``.

    The analogue of :func:`~scripts.maf_to_sites._assign_code` for insertions:
    when two alignment blocks cover the same slot and disagree, there is no
    principled winner, so the slot degrades to unknown bases over the longer of
    the two lengths.
    """
    existing = insertions.get(slot)
    if existing is None or existing == seq:
        insertions[slot] = seq
        return
    insertions[slot] = UNKNOWN * max(len(existing), len(seq))


def load_sample_alignment(
    maf_path: Path,
    contig: str,
    win_start: int,
    win_end: int,
    quality_mask: QualityMask | None = None,
) -> tuple[bytearray, dict[int, str]]:
    """Project one sample's pairwise MAF onto the window ``[win_start, win_end)``.

    Returns the per-anchor call codes (window-local, one byte per reference
    base, using the same code space as ``maf_to_sites``) and the sample's
    insertions keyed by window-local anchor index.

    Insertions are left-anchored: a run of sample bases against reference gaps
    is attached to the reference base that *precedes* it, so a run preceding the
    window's first base anchors outside the window and is dropped.
    """
    length = win_end - win_start
    anchors = bytearray(length)
    insertions: dict[int, str] = {}
    absent = bytearray(length)

    def record_slot(slot: int, seq: str) -> None:
        if seq:
            merge_insertion(insertions, slot, seq)
            if absent[slot]:
                merge_insertion(insertions, slot, "")
        else:
            absent[slot] = 1
            if slot in insertions:
                merge_insertion(insertions, slot, "")

    for block in iter_maf_blocks(maf_path):
        chosen = choose_sample_record(block, contig)
        if chosen is None:
            continue
        ref_record, sample_record = chosen
        if ref_record.start > win_end or ref_record.start + ref_record.size <= win_start:
            continue

        ref_pos = ref_record.start
        # Sample-genome bookkeeping, as in maf_to_sites.load_sample_calls: the
        # quality mask is keyed by the sample's own forward-strand coordinate.
        sample_src = sample_record.src
        sample_start = sample_record.start
        sample_src_size = sample_record.src_size
        sample_minus = sample_record.strand == "-"
        sample_offset = 0
        leading_anchor = ref_record.start - win_start - 1
        prev_idx = leading_anchor if 0 <= leading_anchor < length else None
        previous_ref_seen = False
        pending: list[str] = []

        for ref_char, sample_char in zip(ref_record.text.upper(), sample_record.text.upper()):
            if ref_char == "-":
                if sample_char == GAP:
                    continue
                if sample_minus:
                    sample_coord = sample_src_size - sample_start - 1 - sample_offset
                else:
                    sample_coord = sample_start + sample_offset
                sample_offset += 1
                if prev_idx is None:
                    continue
                if quality_mask is not None and quality_mask.is_low(sample_src, sample_coord):
                    pending.append(UNKNOWN)
                elif sample_char in VALID_BASES:
                    pending.append(sample_char)
                else:
                    pending.append(UNKNOWN)
                continue

            idx = ref_pos - win_start
            ref_pos += 1
            if prev_idx is not None and (pending or previous_ref_seen):
                # Absence is evidence only between two reference bases in a
                # block; a block boundary alone says nothing about this slot.
                record_slot(prev_idx, "".join(pending))
                pending.clear()
            if idx >= length:
                break
            if idx < 0:
                if sample_char != GAP:
                    sample_offset += 1
                prev_idx = None
                continue

            prev_idx = idx
            previous_ref_seen = True
            if sample_char == GAP:
                _assign_code(anchors, idx, NUC_TO_CODE["-"])
                continue

            if sample_minus:
                sample_coord = sample_src_size - sample_start - 1 - sample_offset
            else:
                sample_coord = sample_start + sample_offset
            sample_offset += 1

            if quality_mask is not None and quality_mask.is_low(sample_src, sample_coord):
                _assign_code(anchors, idx, NUC_TO_CODE[UNKNOWN])
            elif sample_char in VALID_BASES:
                _assign_code(anchors, idx, NUC_TO_CODE[sample_char])
            else:
                _assign_code(anchors, idx, NUC_TO_CODE[UNKNOWN])

        if pending and prev_idx is not None:
            record_slot(prev_idx, "".join(pending))

    return anchors, insertions


def load_mask(mask_bed: Path, contig: str, win_start: int, win_end: int) -> bytearray:
    """Window-local flags for reference positions listed in a mask BED."""
    masked = bytearray(win_end - win_start)
    target = normalize_contig(contig)
    with open_text(mask_bed, "rt", errors="ignore") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "track", "browser")):
                continue
            parts = stripped.split()
            if len(parts) < 3:
                raise ValueError(
                    f"Malformed mask BED row in {mask_bed} at line {line_number}: "
                    "expected at least 3 fields"
                )
            if normalize_contig(parts[0]) != target:
                continue
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError as exc:
                raise ValueError(
                    f"Malformed mask BED row in {mask_bed} at line {line_number}: "
                    "start and end must be integers"
                ) from exc
            lo = max(start, win_start)
            hi = min(end, win_end)
            for pos in range(lo, hi):
                masked[pos - win_start] = 1
    return masked


def build_columns(
    insertions_by_sample: list[dict[int, str]],
    length: int,
) -> tuple[list[int], list[int], int]:
    """Lay the alignment out in columns.

    Each reference base gets one anchor column, followed by ``width[r]`` columns
    for insertions anchored to it -- the widest insertion any sample carries
    there. Returns the per-anchor widths, the column offset of each anchor, and
    the total column count.
    """
    widths = [0] * length
    for insertions in insertions_by_sample:
        for slot, seq in insertions.items():
            if len(seq) > widths[slot]:
                widths[slot] = len(seq)
    offsets = [0] * length
    column = 0
    for idx in range(length):
        offsets[idx] = column
        column += 1 + widths[idx]
    return widths, offsets, column


def render_row(
    anchors: bytearray,
    insertions: dict[int, str],
    offsets: list[int],
    total_columns: int,
) -> str:
    """One sample's alignment row: anchor calls in place, insertions
    left-aligned in their slot, every other column a gap."""
    row = bytearray(GAP.encode("ascii") * total_columns)
    for idx, code in enumerate(anchors):
        row[offsets[idx]] = ord(CODE_TO_CHAR[code])
    for slot, seq in insertions.items():
        base = offsets[slot] + 1
        for step, char in enumerate(seq):
            row[base + step] = ord("N" if char == UNKNOWN else char)
    return row.decode("ascii")


def render_reference_row(ref_seq: str, offsets: list[int], total_columns: int) -> str:
    row = bytearray(GAP.encode("ascii") * total_columns)
    for idx, char in enumerate(ref_seq):
        row[offsets[idx]] = ord(char)
    return row.decode("ascii")


def write_fasta(
    out_path: Path,
    rows: list[tuple[str, str]],
    contig: str,
    start: int,
    end: int,
    line_width: int,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(out_path, "wt") as handle:
        for name, seq in rows:
            handle.write(f">{name} {contig}:{start}-{end}\n")
            if line_width <= 0:
                handle.write(f"{seq}\n")
                continue
            for offset in range(0, len(seq), line_width):
                handle.write(f"{seq[offset:offset + line_width]}\n")


def maf_path_for_chunk_root(chunk_root: Path, sample: str, contig: str) -> Path:
    gz = chunk_root / sample / f"{contig}.maf.gz"
    plain = chunk_root / sample / f"{contig}.maf"
    if gz.exists():
        return gz
    if plain.exists():
        return plain
    raise FileNotFoundError(
        f"Missing MAF chunk for sample '{sample}', contig '{contig}' under {chunk_root}"
    )


def discover_chunk_samples(chunk_root: Path, contig: str) -> list[str]:
    samples = []
    for path in sorted(chunk_root.iterdir()):
        if not path.is_dir():
            continue
        if (path / f"{contig}.maf.gz").exists() or (path / f"{contig}.maf").exists():
            samples.append(path.name)
    return samples


def main() -> None:
    args = parse_args()
    reference_fasta = Path(args.reference_fasta)
    contig = args.contig
    out_path = Path(args.out)

    chunk_root = Path(args.maf_chunk_root) if args.maf_chunk_root else None
    maf_dir = Path(args.maf_dir) if args.maf_dir else None
    maf_paths = parse_maf_path_map(args.maf_paths)

    if args.samples:
        samples = list(args.samples)
    elif chunk_root is not None:
        samples = discover_chunk_samples(chunk_root, contig)
    else:
        samples = discover_samples(maf_dir)
    if not samples:
        source = chunk_root if chunk_root is not None else maf_dir
        raise ValueError(f"No samples found under {source}")

    quality_bed_dir: Path | None = None
    if args.quality_bed_dir is not None:
        if args.quality_min is None:
            raise ValueError("--quality-bed-dir requires --quality-min")
        if not 0 <= args.quality_min <= 1:
            raise ValueError("--quality-min must be between 0 and 1")
        quality_bed_dir = Path(args.quality_bed_dir)
    elif args.quality_min is not None:
        raise ValueError("--quality-min requires --quality-bed-dir")

    contig_len = read_contig_length(reference_fasta, contig)
    if contig_len == 0:
        raise ValueError(f"Reference contig '{contig}' has length 0")
    end = contig_len if args.end is None else args.end
    if args.start < 1:
        raise ValueError("--start must be >= 1 (coordinates are 1-based inclusive)")
    if end > contig_len:
        raise ValueError(
            f"--end {end} exceeds length {contig_len} of contig '{contig}'"
        )
    if end < args.start:
        raise ValueError("--end must be >= --start")
    win_start = args.start - 1
    win_end = end
    length = win_end - win_start

    if length > args.max_columns:
        raise ValueError(
            f"Reference window alone would be {length} columns, "
            f"over --max-columns {args.max_columns}. Narrow the window or raise the cap."
        )

    ref_seq = read_contig_region(reference_fasta, contig, win_start, win_end)

    anchors_by_sample: list[bytearray] = []
    insertions_by_sample: list[dict[int, str]] = []
    for sample in samples:
        quality_mask: QualityMask | None = None
        if quality_bed_dir is not None:
            bed_path = quality_bed_for_sample(quality_bed_dir, sample)
            if bed_path is not None:
                quality_mask = load_quality_mask(bed_path, args.quality_min)
        if sample in maf_paths:
            maf_path = maf_paths[sample]
        elif chunk_root is not None:
            maf_path = maf_path_for_chunk_root(chunk_root, sample, contig)
        else:
            maf_path = maf_path_for_sample_with_map(maf_dir, sample, maf_paths)
        anchors, insertions = load_sample_alignment(
            maf_path, contig, win_start, win_end, quality_mask
        )
        anchors_by_sample.append(anchors)
        insertions_by_sample.append(insertions)

    if args.mask_bed is not None:
        masked = load_mask(Path(args.mask_bed), contig, win_start, win_end)
        # A masked reference position has no trustworthy call, so neither does
        # anything inserted there: blank the anchor and drop the insertions,
        # which also removes those columns entirely.
        for anchors, insertions in zip(anchors_by_sample, insertions_by_sample):
            for idx in range(length):
                if masked[idx]:
                    anchors[idx] = NUC_TO_CODE[UNKNOWN]
                    insertions.pop(idx, None)

    _widths, offsets, total_columns = build_columns(insertions_by_sample, length)
    if total_columns > args.max_columns:
        raise ValueError(
            f"Alignment would be {total_columns} columns for a {length} bp window, "
            f"over --max-columns {args.max_columns}. Narrow the window, drop samples, "
            "or raise the cap."
        )

    rows: list[tuple[str, str]] = []
    if args.include_reference:
        rows.append((args.ref_name, render_reference_row(ref_seq, offsets, total_columns)))
    for sample, anchors, insertions in zip(samples, anchors_by_sample, insertions_by_sample):
        rows.append((sample, render_row(anchors, insertions, offsets, total_columns)))

    write_fasta(out_path, rows, contig, args.start, end, args.line_width)

    inserted_slots = sum(1 for width in _widths if width)
    print(
        f"wrote {out_path}: {len(rows)} sequences, {total_columns} columns "
        f"for a {length} bp window ({total_columns - length} insertion columns "
        f"at {inserted_slots} anchors)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
