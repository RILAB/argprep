#!/usr/bin/env python3
"""Per-sample alignment QC statistics for one pairwise reference-vs-query MAF.

One streaming pass over the MAF collects coverage, identity, gap, block and
structural-breakpoint statistics, plus the block coordinates needed to draw
dotplots, so a multi-GB MAF is never read twice.

Outputs, under ``--out-dir``:

- ``<sample>.maf_stats.tsv`` -- one genome-wide row.
- ``<sample>.by_reference_contig.tsv`` -- one row per reference contig in the
  reference ``.fai`` (uncovered contigs included, with zero coverage).
- ``<sample>.by_query_contig.tsv`` -- one row per query contig with alignments.
- ``dotplots/<sample>/<reference_contig>.png`` -- per ``--dotplots`` mode.

The parser is strict on purpose: anything that does not look like a pairwise
reference/query MAF (wrong row count, inconsistent coordinates or source
sizes, unknown reference contigs, a minus-strand reference row) is an error
with file and block context rather than being reinterpreted silently.

See reports/maf_stats_plan.md for metric definitions and design notes.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

try:
    from scripts.common import (
        MafRecord,
        iter_maf_blocks,
        merge_intervals,
        normalize_contig,
    )
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from scripts.common import (
        MafRecord,
        iter_maf_blocks,
        merge_intervals,
        normalize_contig,
    )


GAP = ord("-")
_ACGT = np.zeros(256, dtype=bool)
for _base in b"ACGT":
    _ACGT[_base] = True

DOTPLOT_MODES = ("flagged", "all", "false")

SUMMARY_COLUMNS = [
    "sample",
    "reference_length_bp",
    "covered_reference_bp",
    "reference_coverage_pct",
    "reference_bp_aligned_to_query_base",
    "query_length_bp",
    "query_length_source",
    "covered_query_bp",
    "query_coverage_pct",
    "identity_matches",
    "identity_compared_columns",
    "identity_pct",
    "gap_columns",
    "alignment_columns",
    "gap_fraction_pct",
    "blocks",
    "block_n50_bp",
    "overlapping_reference_bp",
    "strand_flips",
    "reference_contig_jumps",
    "out_of_order_adjacencies",
    "breakpoint_adjacencies",
    "breakpoints_per_gb_covered",
    "min_block_bp",
    "overlap_tolerance_bp",
]

BY_REFERENCE_COLUMNS = [
    "sample",
    "reference_contig",
    "reference_length_bp",
    "covered_reference_bp",
    "reference_coverage_pct",
    "reference_bp_aligned_to_query_base",
    "identity_matches",
    "identity_compared_columns",
    "identity_pct",
    "gap_columns",
    "alignment_columns",
    "gap_fraction_pct",
    "blocks",
    "block_n50_bp",
    "overlapping_reference_bp",
    "breakpoint_adjacencies",
    "dotplot",
]

BY_QUERY_COLUMNS = [
    "sample",
    "query_contig",
    "query_length_bp",
    "query_length_source",
    "covered_query_bp",
    "query_coverage_pct",
    "blocks",
    "blocks_considered",
    "strand_flips",
    "reference_contig_jumps",
    "out_of_order_adjacencies",
    "breakpoint_adjacencies",
    "breakpoints_per_gb_covered",
]


class MafValidationError(ValueError):
    """The MAF is not a valid pairwise reference/query alignment."""


@dataclass(frozen=True)
class Block:
    """Coordinates of one validated block. Query coordinates are forward-strand."""

    ref_contig: str
    ref_start: int
    ref_end: int
    query_contig: str
    query_start: int
    query_end: int
    strand: str

    @property
    def ref_size(self) -> int:
        return self.ref_end - self.ref_start


@dataclass
class ColumnCounts:
    matches: int = 0
    compared: int = 0
    gap_columns: int = 0
    columns: int = 0
    ref_aligned_to_query_base: int = 0

    def add(self, other: "ColumnCounts") -> None:
        self.matches += other.matches
        self.compared += other.compared
        self.gap_columns += other.gap_columns
        self.columns += other.columns
        self.ref_aligned_to_query_base += other.ref_aligned_to_query_base


@dataclass
class ReferenceContigStats:
    counts: ColumnCounts = field(default_factory=ColumnCounts)
    intervals: list[tuple[int, int]] = field(default_factory=list)
    block_sizes: list[int] = field(default_factory=list)
    breakpoint_adjacencies: int = 0


@dataclass
class QueryContigStats:
    src_size: int
    intervals: list[tuple[int, int]] = field(default_factory=list)
    blocks: list[Block] = field(default_factory=list)


@dataclass
class BreakpointCounts:
    blocks_considered: int = 0
    strand_flips: int = 0
    reference_contig_jumps: int = 0
    out_of_order_adjacencies: int = 0
    breakpoint_adjacencies: int = 0


@dataclass
class ScanResult:
    reference_lengths: dict[str, int]
    reference: dict[str, ReferenceContigStats]
    query: dict[str, QueryContigStats]
    block_count: int


def read_fai_lengths(path: Path) -> dict[str, int]:
    """Contig lengths from a ``.fai`` index, in file order."""
    lengths: dict[str, int] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                raise ValueError(f"Malformed .fai row in {path} at line {line_number}")
            try:
                lengths[parts[0]] = int(parts[1])
            except ValueError as exc:
                raise ValueError(
                    f"Malformed .fai row in {path} at line {line_number}: length must be an integer"
                ) from exc
    if not lengths:
        raise ValueError(f"No contigs in {path}")
    return lengths


class ReferenceResolver:
    """Map MAF reference names onto ``.fai`` names: exact match first, then the
    pipeline's contig normalization, but only where that mapping is unique."""

    def __init__(self, names: list[str]):
        self._exact = set(names)
        self._normalized: dict[str, list[str]] = {}
        for name in names:
            self._normalized.setdefault(normalize_contig(name), []).append(name)
        self._cache: dict[str, str] = {}

    def resolve(self, name: str) -> str:
        cached = self._cache.get(name)
        if cached is not None:
            return cached
        if name in self._exact:
            resolved = name
        else:
            candidates = self._normalized.get(normalize_contig(name), [])
            if len(candidates) != 1:
                reason = "is ambiguous" if candidates else "is not in the reference .fai"
                raise MafValidationError(f"reference contig '{name}' {reason}")
            resolved = candidates[0]
        self._cache[name] = resolved
        return resolved


def _ungapped_length(text: str) -> int:
    return len(text) - text.count("-")


def validate_pair(block: list[MafRecord]) -> tuple[MafRecord, MafRecord]:
    """Check that a block is one reference row plus one query row with
    self-consistent coordinates. Context (file, block number) is added by the
    caller."""
    if len(block) != 2:
        raise MafValidationError(
            f"expected exactly 2 sequence rows (reference, query), found {len(block)}"
        )
    for role, record in zip(("reference", "query"), block):
        if record.strand not in ("+", "-"):
            raise MafValidationError(f"{role} row has invalid strand '{record.strand}'")
        if record.start < 0 or record.size < 0 or record.src_size < 0:
            raise MafValidationError(f"{role} row has negative coordinates")
        if record.start + record.size > record.src_size:
            raise MafValidationError(
                f"{role} row interval {record.start}+{record.size} exceeds srcSize {record.src_size}"
            )
        ungapped = _ungapped_length(record.text)
        if ungapped != record.size:
            raise MafValidationError(
                f"{role} row has {ungapped} ungapped bases but size field {record.size}"
            )
    if block[0].strand != "+":
        raise MafValidationError("reference row is on the minus strand")
    return block[0], block[1]


def column_counts(ref_text: str, query_text: str) -> ColumnCounts:
    """Identity, gap, and aligned-base counts for one block's alignment columns."""
    ref = np.frombuffer(ref_text.upper().encode("ascii", "replace"), dtype=np.uint8)
    query = np.frombuffer(query_text.upper().encode("ascii", "replace"), dtype=np.uint8)
    ref_gap = ref == GAP
    query_gap = query == GAP
    both_acgt = _ACGT[ref] & _ACGT[query]
    return ColumnCounts(
        matches=int(np.count_nonzero(both_acgt & (ref == query))),
        compared=int(np.count_nonzero(both_acgt)),
        gap_columns=int(np.count_nonzero(ref_gap | query_gap)),
        columns=int(ref.size),
        ref_aligned_to_query_base=int(np.count_nonzero(~ref_gap & ~query_gap)),
    )


def forward_interval(record: MafRecord) -> tuple[int, int]:
    """A row's interval in forward-strand source coordinates."""
    if record.strand == "-":
        return record.src_size - record.start - record.size, record.src_size - record.start
    return record.start, record.start + record.size


def scan_maf(
    maf_path: Path,
    reference_lengths: dict[str, int],
    query_lengths: dict[str, int] | None = None,
) -> ScanResult:
    """Validate and accumulate every block of one pairwise MAF in a single pass."""
    resolver = ReferenceResolver(list(reference_lengths))
    reference = {name: ReferenceContigStats() for name in reference_lengths}
    query: dict[str, QueryContigStats] = {}
    block_count = 0

    for block_number, records in enumerate(iter_maf_blocks(maf_path), start=1):
        try:
            ref_record, query_record = validate_pair(records)
            ref_name = resolver.resolve(ref_record.src)
            if ref_record.src_size != reference_lengths[ref_name]:
                raise MafValidationError(
                    f"reference srcSize {ref_record.src_size} for '{ref_record.src}' does not "
                    f"match reference .fai length {reference_lengths[ref_name]}"
                )
            query_name = query_record.src
            if query_lengths is not None:
                if query_name not in query_lengths:
                    raise MafValidationError(f"query contig '{query_name}' is not in the query .fai")
                if query_record.src_size != query_lengths[query_name]:
                    raise MafValidationError(
                        f"query srcSize {query_record.src_size} for '{query_name}' does not "
                        f"match query .fai length {query_lengths[query_name]}"
                    )
            query_stats = query.get(query_name)
            if query_stats is None:
                query_stats = query[query_name] = QueryContigStats(src_size=query_record.src_size)
            elif query_stats.src_size != query_record.src_size:
                raise MafValidationError(
                    f"query contig '{query_name}' has inconsistent srcSize "
                    f"({query_stats.src_size} vs {query_record.src_size})"
                )
        except MafValidationError as exc:
            raise MafValidationError(f"{maf_path}, block {block_number}: {exc}") from None

        ref_start = ref_record.start
        ref_end = ref_record.start + ref_record.size
        query_start, query_end = forward_interval(query_record)

        ref_stats = reference[ref_name]
        ref_stats.counts.add(column_counts(ref_record.text, query_record.text))
        ref_stats.intervals.append((ref_start, ref_end))
        ref_stats.block_sizes.append(ref_record.size)
        query_stats.intervals.append((query_start, query_end))
        query_stats.blocks.append(
            Block(ref_name, ref_start, ref_end, query_name, query_start, query_end, query_record.strand)
        )
        block_count += 1

    if block_count == 0:
        raise MafValidationError(f"{maf_path}: no alignment blocks found")
    return ScanResult(reference_lengths, reference, query, block_count)


def union_length(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start for start, end in merge_intervals([iv for iv in intervals if iv[1] > iv[0]]))


def n50(sizes: list[int]) -> int | None:
    total = sum(sizes)
    if total <= 0:
        return None
    running = 0
    for size in sorted(sizes, reverse=True):
        running += size
        if running * 2 >= total:
            return size
    return None  # unreachable


def classify_breakpoints(
    blocks: list[Block],
    reference_order: dict[str, int],
    min_block_bp: int,
    overlap_tolerance_bp: int,
    reference_stats: dict[str, ReferenceContigStats] | None = None,
) -> BreakpointCounts:
    """Walk adjacent blocks along one query contig.

    Blocks shorter than ``min_block_bp`` (reference span) are skipped. Sorting is
    by forward query start with a full tie-break so duplicated query regions
    give the same adjacencies on every run. The three properties are recorded
    independently -- a jump that is also a strand flip counts towards both --
    and ``breakpoint_adjacencies`` counts pairs with any of them. When
    ``reference_stats`` is given, each breakpoint is also credited to the
    reference contig(s) on either side of it.
    """
    considered = sorted(
        (b for b in blocks if b.ref_size >= min_block_bp),
        key=lambda b: (
            b.query_start,
            b.query_end,
            reference_order[b.ref_contig],
            b.ref_start,
            b.ref_end,
            b.strand,
        ),
    )
    counts = BreakpointCounts(blocks_considered=len(considered))
    for previous, current in zip(considered, considered[1:]):
        jump = previous.ref_contig != current.ref_contig
        flip = previous.strand != current.strand
        out_of_order = False
        if not jump and not flip:
            if current.strand == "+":
                out_of_order = current.ref_start < previous.ref_end - overlap_tolerance_bp
            else:
                out_of_order = current.ref_end > previous.ref_start + overlap_tolerance_bp
        counts.reference_contig_jumps += jump
        counts.strand_flips += flip
        counts.out_of_order_adjacencies += out_of_order
        if jump or flip or out_of_order:
            counts.breakpoint_adjacencies += 1
            if reference_stats is not None:
                for name in {previous.ref_contig, current.ref_contig}:
                    reference_stats[name].breakpoint_adjacencies += 1
    return counts


def _pct(numerator: int, denominator: int) -> float | None:
    return None if denominator <= 0 else 100.0 * numerator / denominator


def _per_gb(count: int, denominator: int) -> float | None:
    return None if denominator <= 0 else count / denominator * 1e9


def _fmt(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def safe_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", name) or "_"


@dataclass
class SampleStats:
    summary: dict[str, object]
    by_reference: list[dict[str, object]]
    by_query: list[dict[str, object]]


def compute_stats(
    sample: str,
    scan: ScanResult,
    query_lengths: dict[str, int] | None,
    min_block_bp: int,
    overlap_tolerance_bp: int,
) -> SampleStats:
    reference_order = {name: i for i, name in enumerate(scan.reference_lengths)}
    for stats in scan.reference.values():
        stats.breakpoint_adjacencies = 0

    query_source = "query_fai" if query_lengths is not None else "maf_srcsize"
    total = BreakpointCounts()
    by_query = []
    covered_query_total = 0
    for name, stats in scan.query.items():
        counts = classify_breakpoints(
            stats.blocks, reference_order, min_block_bp, overlap_tolerance_bp, scan.reference
        )
        for attr in ("blocks_considered", "strand_flips", "reference_contig_jumps",
                     "out_of_order_adjacencies", "breakpoint_adjacencies"):
            setattr(total, attr, getattr(total, attr) + getattr(counts, attr))
        covered = union_length(stats.intervals)
        covered_query_total += covered
        by_query.append({
            "sample": sample,
            "query_contig": name,
            "query_length_bp": stats.src_size,
            "query_length_source": query_source,
            "covered_query_bp": covered,
            "query_coverage_pct": _pct(covered, stats.src_size),
            "blocks": len(stats.blocks),
            "blocks_considered": counts.blocks_considered,
            "strand_flips": counts.strand_flips,
            "reference_contig_jumps": counts.reference_contig_jumps,
            "out_of_order_adjacencies": counts.out_of_order_adjacencies,
            "breakpoint_adjacencies": counts.breakpoint_adjacencies,
            "breakpoints_per_gb_covered": _per_gb(counts.breakpoint_adjacencies, covered),
        })

    genome_counts = ColumnCounts()
    by_reference = []
    covered_reference_total = 0
    overlapping_total = 0
    all_sizes: list[int] = []
    for name, length in scan.reference_lengths.items():
        stats = scan.reference[name]
        genome_counts.add(stats.counts)
        covered = union_length(stats.intervals)
        overlapping = sum(stats.block_sizes) - covered
        covered_reference_total += covered
        overlapping_total += overlapping
        all_sizes.extend(stats.block_sizes)
        counts = stats.counts
        by_reference.append({
            "sample": sample,
            "reference_contig": name,
            "reference_length_bp": length,
            "covered_reference_bp": covered,
            "reference_coverage_pct": _pct(covered, length),
            "reference_bp_aligned_to_query_base": counts.ref_aligned_to_query_base,
            "identity_matches": counts.matches,
            "identity_compared_columns": counts.compared,
            "identity_pct": _pct(counts.matches, counts.compared),
            "gap_columns": counts.gap_columns,
            "alignment_columns": counts.columns,
            "gap_fraction_pct": _pct(counts.gap_columns, counts.columns),
            "blocks": len(stats.block_sizes),
            "block_n50_bp": n50(stats.block_sizes),
            "overlapping_reference_bp": overlapping,
            "breakpoint_adjacencies": stats.breakpoint_adjacencies,
            "dotplot": "",
        })

    if query_lengths is not None:
        query_length_total = sum(query_lengths.values())
    else:
        query_length_total = sum(stats.src_size for stats in scan.query.values())
    reference_length_total = sum(scan.reference_lengths.values())

    summary = {
        "sample": sample,
        "reference_length_bp": reference_length_total,
        "covered_reference_bp": covered_reference_total,
        "reference_coverage_pct": _pct(covered_reference_total, reference_length_total),
        "reference_bp_aligned_to_query_base": genome_counts.ref_aligned_to_query_base,
        "query_length_bp": query_length_total,
        "query_length_source": query_source,
        "covered_query_bp": covered_query_total,
        "query_coverage_pct": _pct(covered_query_total, query_length_total),
        "identity_matches": genome_counts.matches,
        "identity_compared_columns": genome_counts.compared,
        "identity_pct": _pct(genome_counts.matches, genome_counts.compared),
        "gap_columns": genome_counts.gap_columns,
        "alignment_columns": genome_counts.columns,
        "gap_fraction_pct": _pct(genome_counts.gap_columns, genome_counts.columns),
        "blocks": scan.block_count,
        "block_n50_bp": n50(all_sizes),
        "overlapping_reference_bp": overlapping_total,
        "strand_flips": total.strand_flips,
        "reference_contig_jumps": total.reference_contig_jumps,
        "out_of_order_adjacencies": total.out_of_order_adjacencies,
        "breakpoint_adjacencies": total.breakpoint_adjacencies,
        "breakpoints_per_gb_covered": _per_gb(total.breakpoint_adjacencies, covered_reference_total),
        "min_block_bp": min_block_bp,
        "overlap_tolerance_bp": overlap_tolerance_bp,
    }
    return SampleStats(summary, by_reference, by_query)


def select_dotplot_contigs(
    by_reference: list[dict[str, object]], mode: str, max_plots: int
) -> list[str]:
    """Reference contigs to plot.

    ``all``: every reference contig with at least one block. ``flagged``:
    contigs with breakpoint evidence, densest (breakpoints per covered Mb)
    first, capped at ``max_plots`` (0 = no cap) so a sample full of small-block
    noise cannot produce thousands of images. ``false``: none.
    """
    if mode == "false":
        return []
    aligned = [row for row in by_reference if row["blocks"]]
    if mode == "all":
        return [str(row["reference_contig"]) for row in aligned]
    flagged = [row for row in aligned if row["breakpoint_adjacencies"]]
    flagged.sort(
        key=lambda row: (
            -row["breakpoint_adjacencies"] / max(int(row["covered_reference_bp"]), 1),
            str(row["reference_contig"]),
        )
    )
    if max_plots > 0:
        flagged = flagged[:max_plots]
    return [str(row["reference_contig"]) for row in flagged]


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(_fmt(row[column]) for column in columns) + "\n")


def write_dotplots(
    out_dir: Path,
    sample: str,
    scan: ScanResult,
    contigs: list[str],
    query_lengths: dict[str, int] | None,
) -> dict[str, str]:
    """Render PNGs; return reference contig -> path relative to ``out_dir``."""
    plot_dir = out_dir / "dotplots" / sample  # matches the Snakefile directory() output
    plot_dir.mkdir(parents=True, exist_ok=True)
    if not contigs:
        return {}
    try:
        from scripts.maf_dotplot import render_dotplot
    except ModuleNotFoundError:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
        from scripts.maf_dotplot import render_dotplot

    wanted = set(contigs)
    blocks_by_reference: dict[str, list[Block]] = {name: [] for name in contigs}
    for stats in scan.query.values():
        for block in stats.blocks:
            if block.ref_contig in wanted:
                blocks_by_reference[block.ref_contig].append(block)
    query_sizes = query_lengths or {name: stats.src_size for name, stats in scan.query.items()}
    query_order = list(query_lengths) if query_lengths is not None else None

    paths: dict[str, str] = {}
    used: set[str] = set()
    for name in contigs:
        stem = safe_filename(name)
        candidate, suffix = stem, 1
        while candidate in used:
            suffix += 1
            candidate = f"{stem}_{suffix}"
        used.add(candidate)
        png = plot_dir / f"{candidate}.png"
        render_dotplot(
            png,
            sample=sample,
            ref_contig=name,
            ref_length=scan.reference_lengths[name],
            blocks=blocks_by_reference[name],
            query_lengths=query_sizes,
            query_order=query_order,
        )
        paths[name] = png.relative_to(out_dir).as_posix()
    return paths


def run(
    maf: Path,
    reference_fai: Path,
    sample: str,
    out_dir: Path,
    query_fai: Path | None = None,
    min_block_bp: int = 0,
    overlap_tolerance_bp: int = 0,
    dotplots: str = "flagged",
    dotplot_max: int = 20,
) -> SampleStats:
    if dotplots not in DOTPLOT_MODES:
        raise ValueError(f"--dotplots must be one of {', '.join(DOTPLOT_MODES)}")
    if min_block_bp < 0 or overlap_tolerance_bp < 0 or dotplot_max < 0:
        raise ValueError("--min-block-bp, --overlap-tolerance-bp and --dotplot-max must be >= 0")
    reference_lengths = read_fai_lengths(reference_fai)
    query_lengths = read_fai_lengths(query_fai) if query_fai is not None else None
    scan = scan_maf(maf, reference_lengths, query_lengths)
    stats = compute_stats(sample, scan, query_lengths, min_block_bp, overlap_tolerance_bp)

    out_dir.mkdir(parents=True, exist_ok=True)
    contigs = select_dotplot_contigs(stats.by_reference, dotplots, dotplot_max)
    plot_paths = write_dotplots(out_dir, sample, scan, contigs, query_lengths)
    for row in stats.by_reference:
        row["dotplot"] = plot_paths.get(str(row["reference_contig"]), "")

    write_tsv(out_dir / f"{sample}.maf_stats.tsv", SUMMARY_COLUMNS, [stats.summary])
    write_tsv(out_dir / f"{sample}.by_reference_contig.tsv", BY_REFERENCE_COLUMNS, stats.by_reference)
    write_tsv(out_dir / f"{sample}.by_query_contig.tsv", BY_QUERY_COLUMNS, stats.by_query)
    return stats


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--maf", required=True, help="Pairwise MAF (.maf or .maf.gz); reference row first")
    ap.add_argument("--reference-fai", required=True, help="Reference .fai index")
    ap.add_argument("--sample", required=True, help="Sample name used in output names and rows")
    ap.add_argument("--out-dir", required=True, help="Output directory")
    ap.add_argument(
        "--query-fai",
        default=None,
        help=(
            "Query genome .fai. Without it, query coverage is measured against only "
            "the query contigs that appear in the MAF, which overstates coverage."
        ),
    )
    ap.add_argument(
        "--min-block-bp",
        type=int,
        default=0,
        help="Ignore blocks with a shorter reference span when classifying breakpoints (default: 0)",
    )
    ap.add_argument(
        "--overlap-tolerance-bp",
        type=int,
        default=0,
        help="Reference overlap allowed between adjacent same-strand blocks before they count as out of order (default: 0)",
    )
    ap.add_argument(
        "--dotplots",
        choices=DOTPLOT_MODES,
        default="flagged",
        help="Which reference contigs to plot (default: flagged = contigs with breakpoints)",
    )
    ap.add_argument(
        "--dotplot-max",
        type=int,
        default=20,
        help="Maximum dotplots per sample in flagged mode, densest first; 0 = no cap (default: 20)",
    )
    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    try:
        stats = run(
            Path(args.maf),
            Path(args.reference_fai),
            args.sample,
            Path(args.out_dir),
            query_fai=Path(args.query_fai) if args.query_fai else None,
            min_block_bp=args.min_block_bp,
            overlap_tolerance_bp=args.overlap_tolerance_bp,
            dotplots=args.dotplots,
            dotplot_max=args.dotplot_max,
        )
    except (MafValidationError, ValueError, FileNotFoundError) as exc:
        sys.exit(f"[maf_stats] ERROR: {exc}")
    s = stats.summary
    print(
        f"[maf_stats] {args.sample}: {s['blocks']} blocks, reference coverage "
        f"{_fmt(s['reference_coverage_pct'])}%, identity {_fmt(s['identity_pct'])}%, "
        f"{s['breakpoint_adjacencies']} breakpoint adjacencies",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
