#!/usr/bin/env python3
"""Per-sample alignment QC statistics for one pairwise reference-vs-query MAF.

One streaming pass over the MAF collects alignment, indel, block and
structural-breakpoint statistics, plus the block coordinates needed to draw
dotplots, so a multi-GB MAF is never read twice.

Alignment quality and assembly quality are reported separately. Without
``--fasta`` only what the MAF can show is measured: contig lengths come from
``srcSize``, and query contigs with no alignment at all are invisible. With
``--fasta`` (the query assembly) the script also reports assembly contiguity,
N content and scaffold gaps, GC, soft-masking, telomeric contig ends and the
content of unaligned contigs, and classifies breakpoints next to N gaps.

Outputs, under ``--out-dir``:

- ``<sample>.maf_stats.tsv`` -- one genome-wide row.
- ``<sample>.by_reference_contig.tsv`` -- one row per reference contig in the
  reference ``.fai`` (uncovered contigs included).
- ``<sample>.by_query_contig.tsv`` -- one row per query contig with alignments;
  with ``--fasta``, every assembly contig.
- ``<sample>.breakpoints.tsv`` -- one row per breakpoint adjacency, with its
  location on the query contig and the reference coordinates at the junction.
- ``<sample>.nested_blocks.tsv`` -- blocks whose query interval lies inside a
  larger block's (secondary or transposed alignments), kept out of the
  breakpoint walk.
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
        open_text,
    )
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from scripts.common import (
        MafRecord,
        iter_maf_blocks,
        merge_intervals,
        normalize_contig,
        open_text,
    )


GAP = ord("-")
_ACGT = np.zeros(256, dtype=bool)
for _base in b"ACGT":
    _ACGT[_base] = True
_N = np.zeros(256, dtype=bool)
_N[ord("N")] = True
# Upper-case ASCII letters in place without a str copy: clear bit 0x20 on a-z.
_UPPER = np.arange(256, dtype=np.uint8)
_UPPER[ord("a"):ord("z") + 1] -= 32

# Columns per slice when comparing a block's two rows. AnchorWave blocks can be
# hundreds of millions of columns; slicing bounds the NumPy temporaries.
COLUMN_CHUNK = 8_000_000

# Plant telomere repeat (TTTAGGG)n. A sequence end counts as telomeric when at
# least TELOMERE_MIN_FRACTION of its terminal TELOMERE_WINDOW_BP is that repeat
# (either orientation): a real array sits at the very end and is kb-long, unlike
# scattered interstitial copies.
TELOMERE_WINDOW_BP = 1_000
TELOMERE_MIN_FRACTION = 0.5
# N runs at least this long are scaffold gaps; shorter runs are ambiguous bases.
MIN_GAP_N = 10
# "Major" sequences: at least this fraction of the longest sequence's length.
# Picks out the chromosomes of a chromosome-level assembly however much small
# unplaced-scaffold sequence it also carries (a share-of-assembly cut-off does
# not: 7% of debris pulled 74 scaffolds into a 95% cover).
MAJOR_MIN_FRACTION_OF_LONGEST = 0.10
_TELOMERE_FORWARD = b"TTTAGGG"
_TELOMERE_REVERSE = b"CCCTAAA"

DOTPLOT_MODES = ("flagged", "all", "false")

SUMMARY_COLUMNS = [
    "sample",
    "reference_length_bp",
    "aligned_reference_bp",
    "aligned_reference_pct",
    "block_span_reference_bp",
    "block_span_reference_pct",
    "query_length_bp",
    "query_length_source",
    "aligned_query_bp",
    "aligned_query_pct",
    "block_span_query_bp",
    "block_span_query_pct",
    "identity_matches",
    "identity_compared_columns",
    "identity_pct",
    "alignment_columns",
    "insertion_columns",
    "insertion_column_pct",
    "deletion_columns",
    "deletion_column_pct",
    "query_n_bases_in_blocks",
    "blocks",
    "overlapping_reference_bp",
    "overlapping_query_bp",
    "nested_blocks",
    "nested_block_reference_bp",
    "aligned_query_contigs",
    "aligned_query_contig_n50_bp",
    "strand_flips",
    "reference_contig_jumps",
    "out_of_order_adjacencies",
    "breakpoint_adjacencies",
    "breakpoints_near_contig_end",
    "breakpoints_near_n_gap",
    "breakpoints_interior",
    "breakpoints_per_gb_aligned",
    "assembly_length_bp",
    "assembly_sequences",
    "assembly_scaffold_n50_bp",
    "assembly_scaffold_l50",
    "assembly_largest_sequence_bp",
    "assembly_contig_pieces",
    "assembly_contig_n50_bp",
    "assembly_contig_l50",
    "assembly_n_bp",
    "assembly_n_pct",
    "assembly_n_gaps",
    "assembly_gc_pct",
    "assembly_softmasked_pct",
    "assembly_major_sequences",
    "assembly_major_telomeric_ends",
    "assembly_minor_sequences_with_telomere",
    "unaligned_contigs",
    "unaligned_contig_bp",
    "unaligned_contig_n_pct",
    "unaligned_contig_softmasked_pct",
    "min_block_bp",
    "overlap_tolerance_bp",
    "breakpoint_context_bp",
]

BY_REFERENCE_COLUMNS = [
    "sample",
    "reference_contig",
    "reference_length_bp",
    "aligned_reference_bp",
    "aligned_reference_pct",
    "block_span_reference_bp",
    "block_span_reference_pct",
    "identity_matches",
    "identity_compared_columns",
    "identity_pct",
    "alignment_columns",
    "insertion_columns",
    "insertion_column_pct",
    "deletion_columns",
    "deletion_column_pct",
    "blocks",
    "query_contigs",
    "overlapping_reference_bp",
    "breakpoint_adjacencies",
    "dotplot",
]

BY_QUERY_COLUMNS = [
    "sample",
    "query_contig",
    "query_length_bp",
    "query_length_source",
    "aligned",
    "aligned_query_bp",
    "block_span_query_bp",
    "block_span_query_pct",
    "overlapping_query_bp",
    "unaligned_start_bp",
    "unaligned_end_bp",
    "blocks",
    "blocks_considered",
    "nested_blocks",
    "reference_contigs",
    "strand_flips",
    "reference_contig_jumps",
    "out_of_order_adjacencies",
    "breakpoint_adjacencies",
    "n_bp",
    "n_gaps",
    "contig_pieces",
    "gc_pct",
    "softmasked_pct",
    "major",
    "telomere_start",
    "telomere_end",
]

NESTED_COLUMNS = [
    "sample",
    "query_contig",
    "query_start",
    "query_end",
    "strand",
    "reference_contig",
    "reference_start",
    "reference_end",
    "container_query_start",
    "container_query_end",
    "container_reference_contig",
    "container_reference_start",
    "container_reference_end",
    "container_strand",
]

BREAKPOINT_COLUMNS = [
    "sample",
    "query_contig",
    "query_contig_length_bp",
    "query_junction_start",
    "query_junction_end",
    "distance_to_contig_end_bp",
    "location",
    "near_contig_end",
    "near_n_gap",
    "strand_flip",
    "reference_contig_jump",
    "out_of_order",
    "left_reference_contig",
    "left_reference_pos",
    "left_strand",
    "right_reference_contig",
    "right_reference_pos",
    "right_strand",
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

    def ref_pos_at_query_end(self) -> int:
        """Reference coordinate opposite the block's forward-query end."""
        return self.ref_end if self.strand == "+" else self.ref_start

    def ref_pos_at_query_start(self) -> int:
        """Reference coordinate opposite the block's forward-query start."""
        return self.ref_start if self.strand == "+" else self.ref_end


@dataclass
class ColumnCounts:
    """Per-column tallies. ``aligned`` = columns with a base in both rows;
    ``insertions`` = query base opposite a reference gap; ``deletions`` =
    reference base opposite a query gap."""

    matches: int = 0
    compared: int = 0
    columns: int = 0
    aligned: int = 0
    insertions: int = 0
    deletions: int = 0
    query_n: int = 0

    def add(self, other: "ColumnCounts") -> None:
        self.matches += other.matches
        self.compared += other.compared
        self.columns += other.columns
        self.aligned += other.aligned
        self.insertions += other.insertions
        self.deletions += other.deletions
        self.query_n += other.query_n


@dataclass
class ReferenceContigStats:
    counts: ColumnCounts = field(default_factory=ColumnCounts)
    intervals: list[tuple[int, int]] = field(default_factory=list)
    block_sizes: list[int] = field(default_factory=list)
    query_contigs: set[str] = field(default_factory=set)
    breakpoint_adjacencies: int = 0


@dataclass
class QueryContigStats:
    src_size: int
    counts: ColumnCounts = field(default_factory=ColumnCounts)
    intervals: list[tuple[int, int]] = field(default_factory=list)
    blocks: list[Block] = field(default_factory=list)


@dataclass
class Breakpoint:
    previous: Block
    current: Block
    strand_flip: bool
    reference_contig_jump: bool
    out_of_order: bool


@dataclass
class BreakpointCounts:
    blocks_considered: int = 0
    strand_flips: int = 0
    reference_contig_jumps: int = 0
    out_of_order_adjacencies: int = 0
    breakpoint_adjacencies: int = 0
    breakpoints: list[Breakpoint] = field(default_factory=list)
    # (nested block, the block whose query interval contains it)
    nested: list[tuple[Block, Block]] = field(default_factory=list)


@dataclass
class ScanResult:
    reference_lengths: dict[str, int]
    reference: dict[str, ReferenceContigStats]
    query: dict[str, QueryContigStats]
    block_count: int


@dataclass
class AssemblyContig:
    length: int
    n_bp: int
    acgt_bp: int
    gc_bp: int
    softmasked_bp: int
    n_gaps: list[tuple[int, int]]
    pieces: list[int]
    telomere_start: bool
    telomere_end: bool


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


def iter_fasta(path: Path):
    """Yield ``(name, sequence_bytes)`` per record; plain or gzip-compressed."""
    name: str | None = None
    chunks: list[bytes] = []
    with open_text(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(chunks)
                header = line[1:].split()
                if not header:
                    raise ValueError(f"FASTA record with an empty name in {path}")
                name = header[0].decode("ascii", "replace")
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
            elif line.strip():
                raise ValueError(f"{path} does not start with a FASTA header")
    if name is not None:
        yield name, b"".join(chunks)


def _telomeric(window: bytes) -> bool:
    if not window:
        return False
    window = window.upper()
    repeat_bp = len(_TELOMERE_FORWARD) * max(
        window.count(_TELOMERE_FORWARD), window.count(_TELOMERE_REVERSE)
    )
    return repeat_bp >= TELOMERE_MIN_FRACTION * len(window)


_N_GAP = re.compile(rb"[Nn]{%d,}" % MIN_GAP_N)


def contig_stats(sequence: bytes) -> AssemblyContig:
    counts = np.bincount(np.frombuffer(sequence, dtype=np.uint8), minlength=256)
    n_bp = int(counts[ord("N")] + counts[ord("n")])
    gc_bp = int(sum(counts[ord(c)] for c in "GCgc"))
    acgt_bp = int(sum(counts[ord(c)] for c in "ACGTacgt"))
    softmasked_bp = int(counts[ord("a"):ord("z") + 1].sum())
    gaps = [(m.start(), m.end()) for m in _N_GAP.finditer(sequence)] if n_bp >= MIN_GAP_N else []
    pieces, previous_end = [], 0
    for start, end in gaps + [(len(sequence), len(sequence))]:
        if start > previous_end:
            pieces.append(start - previous_end)
        previous_end = end
    return AssemblyContig(
        length=len(sequence),
        n_bp=n_bp,
        acgt_bp=acgt_bp,
        gc_bp=gc_bp,
        softmasked_bp=softmasked_bp,
        n_gaps=gaps,
        pieces=pieces,
        telomere_start=_telomeric(sequence[:TELOMERE_WINDOW_BP]),
        telomere_end=_telomeric(sequence[-TELOMERE_WINDOW_BP:]),
    )


def read_assembly(path: Path) -> dict[str, AssemblyContig]:
    """Per-contig assembly statistics, one contig in memory at a time."""
    contigs: dict[str, AssemblyContig] = {}
    for name, sequence in iter_fasta(path):
        if name in contigs:
            raise ValueError(f"Duplicate contig name '{name}' in {path}")
        contigs[name] = contig_stats(sequence)
    if not contigs:
        raise ValueError(f"No sequences in {path}")
    return contigs


class ContigResolver:
    """Map MAF sequence names onto index (``.fai``/FASTA) names.

    Tried in order, each only when it identifies exactly one contig:

    1. exact match;
    2. an index name preceded by a dotted prefix -- aligners often write
       ``<assembly file>.<contig>``, e.g. ``Zm-CML103-REFERENCE-NAM-1.0.fa.chr6``
       for FASTA contig ``chr6`` (the longest such suffix wins);
    3. the pipeline's contig normalization (``chr05`` -> ``5``).
    """

    def __init__(self, names: list[str], role: str = "reference", source: str = "reference .fai"):
        self._exact = set(names)
        self._normalized: dict[str, list[str]] = {}
        for name in names:
            self._normalized.setdefault(normalize_contig(name), []).append(name)
        self._role = role
        self._source = source
        self._cache: dict[str, str] = {}

    def _by_suffix(self, name: str) -> str | None:
        parts = name.split(".")
        for i in range(1, len(parts)):
            candidate = ".".join(parts[i:])
            if candidate in self._exact:
                return candidate  # longest dotted suffix first
        return None

    def resolve(self, name: str) -> str:
        cached = self._cache.get(name)
        if cached is not None:
            return cached
        if name in self._exact:
            resolved = name
        else:
            resolved = self._by_suffix(name)
            if resolved is None:
                candidates = self._normalized.get(normalize_contig(name), [])
                if len(candidates) != 1:
                    reason = "is ambiguous in" if candidates else "is not in"
                    raise MafValidationError(f"{self._role} contig '{name}' {reason} the {self._source}")
                resolved = candidates[0]
        self._cache[name] = resolved
        return resolved


ReferenceResolver = ContigResolver  # backwards-compatible name


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
    """Identity, indel, and aligned-base counts for one block's columns.

    Works through the rows in slices of ``COLUMN_CHUNK`` columns so a
    chromosome-scale block never materialises full-length temporaries.
    """
    counts = ColumnCounts()
    for offset in range(0, len(ref_text), COLUMN_CHUNK):
        ref = _UPPER[np.frombuffer(
            ref_text[offset:offset + COLUMN_CHUNK].encode("ascii", "replace"), dtype=np.uint8
        )]
        query = _UPPER[np.frombuffer(
            query_text[offset:offset + COLUMN_CHUNK].encode("ascii", "replace"), dtype=np.uint8
        )]
        ref_gap = ref == GAP
        query_gap = query == GAP
        both_acgt = _ACGT[ref] & _ACGT[query]
        counts.matches += int(np.count_nonzero(both_acgt & (ref == query)))
        counts.compared += int(np.count_nonzero(both_acgt))
        counts.columns += int(ref.size)
        counts.aligned += int(np.count_nonzero(~ref_gap & ~query_gap))
        counts.insertions += int(np.count_nonzero(ref_gap & ~query_gap))
        counts.deletions += int(np.count_nonzero(~ref_gap & query_gap))
        counts.query_n += int(np.count_nonzero(_N[query]))
    return counts


def forward_interval(record: MafRecord) -> tuple[int, int]:
    """A row's interval in forward-strand source coordinates."""
    if record.strand == "-":
        return record.src_size - record.start - record.size, record.src_size - record.start
    return record.start, record.start + record.size


def scan_maf(
    maf_path: Path,
    reference_lengths: dict[str, int],
    query_lengths: dict[str, int] | None = None,
    query_length_label: str = "query .fai",
) -> ScanResult:
    """Validate and accumulate every block of one pairwise MAF in a single pass."""
    resolver = ContigResolver(list(reference_lengths))
    query_resolver = (
        ContigResolver(list(query_lengths), role="query", source=query_length_label)
        if query_lengths is not None
        else None
    )
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
            if query_resolver is not None:
                query_name = query_resolver.resolve(query_name)
                if query_record.src_size != query_lengths[query_name]:
                    raise MafValidationError(
                        f"query srcSize {query_record.src_size} for '{query_name}' does not "
                        f"match {query_length_label} length {query_lengths[query_name]}"
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
        counts = column_counts(ref_record.text, query_record.text)

        ref_stats = reference[ref_name]
        ref_stats.counts.add(counts)
        ref_stats.intervals.append((ref_start, ref_end))
        ref_stats.block_sizes.append(ref_record.size)
        ref_stats.query_contigs.add(query_name)
        query_stats.counts.add(counts)
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


def n50_l50(sizes: list[int]) -> tuple[int | None, int | None]:
    """N50 and L50 (count of pieces needed to reach half the total)."""
    total = sum(sizes)
    if total <= 0:
        return None, None
    running = 0
    for count, size in enumerate(sorted(sizes, reverse=True), start=1):
        running += size
        if running * 2 >= total:
            return size, count
    return None, None  # unreachable


def classify_breakpoints(
    blocks: list[Block],
    reference_order: dict[str, int],
    min_block_bp: int,
    overlap_tolerance_bp: int,
    reference_stats: dict[str, ReferenceContigStats] | None = None,
) -> BreakpointCounts:
    """Walk adjacent blocks along one query contig.

    Blocks shorter than ``min_block_bp`` (reference span) are skipped. A block
    whose query interval lies inside another block's is *nested* (the same
    query sequence aligned twice); nested blocks are returned separately and
    left out of the walk, which would otherwise pair them with their container
    as fake breakpoints spanning the whole container. Sorting is by forward
    query start with a full tie-break so the result is the same on every run. The three properties are recorded
    independently -- a jump that is also a strand flip counts towards both --
    and ``breakpoint_adjacencies`` counts pairs with any of them. When
    ``reference_stats`` is given, each breakpoint is also credited to the
    reference contig(s) on either side of it.
    """
    # Containing block first among equal starts, so anything whose query
    # interval lies inside an earlier block's is recognised as nested.
    ordered = sorted(
        (b for b in blocks if b.ref_size >= min_block_bp),
        key=lambda b: (
            b.query_start,
            -b.query_end,
            reference_order[b.ref_contig],
            b.ref_start,
            b.ref_end,
            b.strand,
        ),
    )
    counts = BreakpointCounts(blocks_considered=len(ordered))
    considered: list[Block] = []
    container: Block | None = None
    for block in ordered:
        if container is not None and block.query_end <= container.query_end:
            # The same query sequence is also inside a larger block: a secondary
            # or transposed alignment, not a step along the query contig.
            counts.nested.append((block, container))
            continue
        considered.append(block)
        if container is None or block.query_end > container.query_end:
            container = block
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
            counts.breakpoints.append(Breakpoint(previous, current, flip, jump, out_of_order))
            if reference_stats is not None:
                for name in {previous.ref_contig, current.ref_contig}:
                    reference_stats[name].breakpoint_adjacencies += 1
    return counts


def locate_breakpoint(
    breakpoint: Breakpoint,
    contig_length: int,
    n_gaps: list[tuple[int, int]] | None,
    context_bp: int,
) -> dict[str, object]:
    """Where a breakpoint sits on its query contig.

    The junction is the query interval between the two blocks. ``contig_end``
    when the junction is within ``context_bp`` of either contig end (assembly or
    scaffolding); otherwise ``n_gap`` when an N run lies within ``context_bp``
    (needs the FASTA; ``near_n_gap`` is None without it); otherwise
    ``interior``.
    """
    lo = min(breakpoint.previous.query_end, breakpoint.current.query_start)
    hi = max(breakpoint.previous.query_end, breakpoint.current.query_start)
    distance = max(0, min(lo, contig_length - hi))
    near_end = distance <= context_bp
    near_gap: bool | None = None
    if n_gaps is not None:
        near_gap = any(start < hi + context_bp and end > lo - context_bp for start, end in n_gaps)
    if near_end:
        location = "contig_end"
    elif near_gap:
        location = "n_gap"
    else:
        location = "interior"
    return {
        "query_junction_start": lo,
        "query_junction_end": hi,
        "distance_to_contig_end_bp": distance,
        "location": location,
        "near_contig_end": near_end,
        "near_n_gap": near_gap,
    }


def _pct(numerator: int, denominator: int) -> float | None:
    return None if denominator <= 0 else 100.0 * numerator / denominator


def _per_gb(count: int, denominator: int) -> float | None:
    return None if denominator <= 0 else count / denominator * 1e9


def _fmt(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
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
    breakpoints: list[dict[str, object]]
    nested: list[dict[str, object]]


def _assembly_summary(
    assembly: dict[str, AssemblyContig] | None, aligned_names: set[str]
) -> dict[str, object]:
    keys = [c for c in SUMMARY_COLUMNS if c.startswith(("assembly_", "unaligned_"))]
    if assembly is None:
        return {key: None for key in keys}
    contigs = list(assembly.values())
    scaffold_n50, scaffold_l50 = n50_l50([c.length for c in contigs])
    pieces = [piece for c in contigs for piece in c.pieces]
    contig_n50, contig_l50 = n50_l50(pieces)
    length = sum(c.length for c in contigs)
    acgt = sum(c.acgt_bp for c in contigs)
    masked = any(c.softmasked_bp for c in contigs)
    major = major_sequences(assembly)
    unaligned = [c for name, c in assembly.items() if name not in aligned_names]
    unaligned_bp = sum(c.length for c in unaligned)
    return {
        "assembly_length_bp": length,
        "assembly_sequences": len(contigs),
        "assembly_scaffold_n50_bp": scaffold_n50,
        "assembly_scaffold_l50": scaffold_l50,
        "assembly_largest_sequence_bp": max(c.length for c in contigs),
        "assembly_contig_pieces": len(pieces),
        "assembly_contig_n50_bp": contig_n50,
        "assembly_contig_l50": contig_l50,
        "assembly_n_bp": sum(c.n_bp for c in contigs),
        "assembly_n_pct": _pct(sum(c.n_bp for c in contigs), length),
        "assembly_n_gaps": sum(len(c.n_gaps) for c in contigs),
        "assembly_gc_pct": _pct(sum(c.gc_bp for c in contigs), acgt),
        # An all-uppercase FASTA was never soft-masked: 0% would mislead.
        "assembly_softmasked_pct": _pct(sum(c.softmasked_bp for c in contigs), length) if masked else None,
        "assembly_major_sequences": len(major),
        "assembly_major_telomeric_ends": sum(
            assembly[name].telomere_start + assembly[name].telomere_end for name in major
        ),
        "assembly_minor_sequences_with_telomere": sum(
            1 for name, c in assembly.items() if name not in major and (c.telomere_start or c.telomere_end)
        ),
        "unaligned_contigs": len(unaligned),
        "unaligned_contig_bp": unaligned_bp,
        "unaligned_contig_n_pct": _pct(sum(c.n_bp for c in unaligned), unaligned_bp),
        "unaligned_contig_softmasked_pct": (
            _pct(sum(c.softmasked_bp for c in unaligned), unaligned_bp) if masked else None
        ),
    }


def major_sequences(assembly: dict[str, AssemblyContig]) -> set[str]:
    """Sequences at least MAJOR_MIN_FRACTION_OF_LONGEST of the longest one's length."""
    longest = max(c.length for c in assembly.values())
    cutoff = MAJOR_MIN_FRACTION_OF_LONGEST * longest
    return {name for name, contig in assembly.items() if contig.length >= cutoff}


def compute_stats(
    sample: str,
    scan: ScanResult,
    query_lengths: dict[str, int] | None,
    min_block_bp: int,
    overlap_tolerance_bp: int,
    assembly: dict[str, AssemblyContig] | None = None,
    breakpoint_context_bp: int = 1_000_000,
) -> SampleStats:
    reference_order = {name: i for i, name in enumerate(scan.reference_lengths)}
    for stats in scan.reference.values():
        stats.breakpoint_adjacencies = 0

    if assembly is not None:
        query_source = "query_fasta"
        query_lengths = {name: contig.length for name, contig in assembly.items()}
    elif query_lengths is not None:
        query_source = "query_fai"
    else:
        query_source = "maf_srcsize"

    total = BreakpointCounts()
    breakpoint_rows: list[dict[str, object]] = []
    nested_rows: list[dict[str, object]] = []
    overlapping_query_total = 0
    query_rows: dict[str, dict[str, object]] = {}
    span_query_total = 0
    aligned_query_total = 0
    for name, stats in scan.query.items():
        counts = classify_breakpoints(
            stats.blocks, reference_order, min_block_bp, overlap_tolerance_bp, scan.reference
        )
        for attr in ("blocks_considered", "strand_flips", "reference_contig_jumps",
                     "out_of_order_adjacencies", "breakpoint_adjacencies"):
            setattr(total, attr, getattr(total, attr) + getattr(counts, attr))
        for block, holder in counts.nested:
            nested_rows.append({
                "sample": sample,
                "query_contig": name,
                "query_start": block.query_start,
                "query_end": block.query_end,
                "strand": block.strand,
                "reference_contig": block.ref_contig,
                "reference_start": block.ref_start,
                "reference_end": block.ref_end,
                "container_query_start": holder.query_start,
                "container_query_end": holder.query_end,
                "container_reference_contig": holder.ref_contig,
                "container_reference_start": holder.ref_start,
                "container_reference_end": holder.ref_end,
                "container_strand": holder.strand,
            })
        contig = assembly.get(name) if assembly is not None else None
        for bp in counts.breakpoints:
            located = locate_breakpoint(
                bp, stats.src_size, contig.n_gaps if contig is not None else None, breakpoint_context_bp
            )
            breakpoint_rows.append({
                "sample": sample,
                "query_contig": name,
                "query_contig_length_bp": stats.src_size,
                **located,
                "strand_flip": bp.strand_flip,
                "reference_contig_jump": bp.reference_contig_jump,
                "out_of_order": bp.out_of_order,
                "left_reference_contig": bp.previous.ref_contig,
                "left_reference_pos": bp.previous.ref_pos_at_query_end(),
                "left_strand": bp.previous.strand,
                "right_reference_contig": bp.current.ref_contig,
                "right_reference_pos": bp.current.ref_pos_at_query_start(),
                "right_strand": bp.current.strand,
            })
        span = union_length(stats.intervals)
        overlapping_query = sum(end - start for start, end in stats.intervals) - span
        span_query_total += span
        overlapping_query_total += overlapping_query
        aligned_query_total += stats.counts.aligned
        query_rows[name] = {
            "sample": sample,
            "query_contig": name,
            "query_length_bp": stats.src_size,
            "query_length_source": query_source,
            "aligned": True,
            "aligned_query_bp": stats.counts.aligned,
            "block_span_query_bp": span,
            "block_span_query_pct": _pct(span, stats.src_size),
            "overlapping_query_bp": overlapping_query,
            "unaligned_start_bp": min(start for start, _ in stats.intervals),
            "unaligned_end_bp": stats.src_size - max(end for _, end in stats.intervals),
            "blocks": len(stats.blocks),
            "blocks_considered": counts.blocks_considered,
            "nested_blocks": len(counts.nested),
            "reference_contigs": len({b.ref_contig for b in stats.blocks}),
            "strand_flips": counts.strand_flips,
            "reference_contig_jumps": counts.reference_contig_jumps,
            "out_of_order_adjacencies": counts.out_of_order_adjacencies,
            "breakpoint_adjacencies": counts.breakpoint_adjacencies,
        }

    empty_query = {
        "aligned": False, "aligned_query_bp": 0, "block_span_query_bp": 0, "block_span_query_pct": 0.0,
        "overlapping_query_bp": 0, "unaligned_start_bp": None, "unaligned_end_bp": None, "blocks": 0,
        "blocks_considered": 0, "nested_blocks": 0, "reference_contigs": 0, "strand_flips": 0, "reference_contig_jumps": 0,
        "out_of_order_adjacencies": 0, "breakpoint_adjacencies": 0,
    }
    by_query = []
    masked = assembly is not None and any(c.softmasked_bp for c in assembly.values())
    major = major_sequences(assembly) if assembly is not None else set()
    order = list(assembly) if assembly is not None else list(scan.query)
    for name in order:
        row = query_rows.get(name)
        if row is None:  # assembly contig with no alignment
            row = {"sample": sample, "query_contig": name, "query_length_bp": assembly[name].length,
                   "query_length_source": query_source, **empty_query}
        contig = assembly.get(name) if assembly is not None else None
        row.update({
            "n_bp": contig.n_bp if contig else None,
            "n_gaps": len(contig.n_gaps) if contig else None,
            "contig_pieces": len(contig.pieces) if contig else None,
            "gc_pct": _pct(contig.gc_bp, contig.acgt_bp) if contig else None,
            "softmasked_pct": _pct(contig.softmasked_bp, contig.length) if contig and masked else None,
            "major": (name in major) if contig else None,
            "telomere_start": contig.telomere_start if contig else None,
            "telomere_end": contig.telomere_end if contig else None,
        })
        by_query.append(row)

    genome_counts = ColumnCounts()
    by_reference = []
    span_reference_total = 0
    overlapping_total = 0
    for name, length in scan.reference_lengths.items():
        stats = scan.reference[name]
        genome_counts.add(stats.counts)
        span = union_length(stats.intervals)
        overlapping = sum(stats.block_sizes) - span
        span_reference_total += span
        overlapping_total += overlapping
        counts = stats.counts
        by_reference.append({
            "sample": sample,
            "reference_contig": name,
            "reference_length_bp": length,
            "aligned_reference_bp": counts.aligned,
            "aligned_reference_pct": _pct(counts.aligned, length),
            "block_span_reference_bp": span,
            "block_span_reference_pct": _pct(span, length),
            "identity_matches": counts.matches,
            "identity_compared_columns": counts.compared,
            "identity_pct": _pct(counts.matches, counts.compared),
            "alignment_columns": counts.columns,
            "insertion_columns": counts.insertions,
            "insertion_column_pct": _pct(counts.insertions, counts.columns),
            "deletion_columns": counts.deletions,
            "deletion_column_pct": _pct(counts.deletions, counts.columns),
            "blocks": len(stats.block_sizes),
            "query_contigs": len(stats.query_contigs),
            "overlapping_reference_bp": overlapping,
            "breakpoint_adjacencies": stats.breakpoint_adjacencies,
            "dotplot": "",
        })

    if query_lengths is not None:
        query_length_total = sum(query_lengths.values())
    else:
        query_length_total = sum(stats.src_size for stats in scan.query.values())
    reference_length_total = sum(scan.reference_lengths.values())
    aligned_contig_n50, _ = n50_l50([stats.src_size for stats in scan.query.values()])
    locations = [row["location"] for row in breakpoint_rows]

    summary = {
        "sample": sample,
        "reference_length_bp": reference_length_total,
        "aligned_reference_bp": genome_counts.aligned,
        "aligned_reference_pct": _pct(genome_counts.aligned, reference_length_total),
        "block_span_reference_bp": span_reference_total,
        "block_span_reference_pct": _pct(span_reference_total, reference_length_total),
        "query_length_bp": query_length_total,
        "query_length_source": query_source,
        "aligned_query_bp": aligned_query_total,
        "aligned_query_pct": _pct(aligned_query_total, query_length_total),
        "block_span_query_bp": span_query_total,
        "block_span_query_pct": _pct(span_query_total, query_length_total),
        "identity_matches": genome_counts.matches,
        "identity_compared_columns": genome_counts.compared,
        "identity_pct": _pct(genome_counts.matches, genome_counts.compared),
        "alignment_columns": genome_counts.columns,
        "insertion_columns": genome_counts.insertions,
        "insertion_column_pct": _pct(genome_counts.insertions, genome_counts.columns),
        "deletion_columns": genome_counts.deletions,
        "deletion_column_pct": _pct(genome_counts.deletions, genome_counts.columns),
        "query_n_bases_in_blocks": genome_counts.query_n,
        "blocks": scan.block_count,
        "overlapping_reference_bp": overlapping_total,
        "overlapping_query_bp": overlapping_query_total,
        "nested_blocks": len(nested_rows),
        "nested_block_reference_bp": sum(row["reference_end"] - row["reference_start"] for row in nested_rows),
        "aligned_query_contigs": len(scan.query),
        "aligned_query_contig_n50_bp": aligned_contig_n50,
        "strand_flips": total.strand_flips,
        "reference_contig_jumps": total.reference_contig_jumps,
        "out_of_order_adjacencies": total.out_of_order_adjacencies,
        "breakpoint_adjacencies": total.breakpoint_adjacencies,
        "breakpoints_near_contig_end": locations.count("contig_end"),
        "breakpoints_near_n_gap": locations.count("n_gap") if assembly is not None else None,
        "breakpoints_interior": locations.count("interior"),
        "breakpoints_per_gb_aligned": _per_gb(total.breakpoint_adjacencies, genome_counts.aligned),
        **_assembly_summary(assembly, set(scan.query)),
        "min_block_bp": min_block_bp,
        "overlap_tolerance_bp": overlap_tolerance_bp,
        "breakpoint_context_bp": breakpoint_context_bp,
    }
    return SampleStats(summary, by_reference, by_query, breakpoint_rows, nested_rows)


def select_dotplot_contigs(
    by_reference: list[dict[str, object]], mode: str, max_plots: int
) -> list[str]:
    """Reference contigs to plot.

    ``all``: every reference contig with at least one block. ``flagged``:
    contigs with breakpoint evidence, densest (breakpoints per spanned Mb)
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
            -row["breakpoint_adjacencies"] / max(int(row["block_span_reference_bp"]), 1),
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
    fasta: Path | None = None,
    breakpoint_context_bp: int = 1_000_000,
) -> SampleStats:
    if dotplots not in DOTPLOT_MODES:
        raise ValueError(f"--dotplots must be one of {', '.join(DOTPLOT_MODES)}")
    if min(min_block_bp, overlap_tolerance_bp, dotplot_max, breakpoint_context_bp) < 0:
        raise ValueError(
            "--min-block-bp, --overlap-tolerance-bp, --dotplot-max and "
            "--breakpoint-context-bp must be >= 0"
        )
    if fasta is not None and query_fai is not None:
        raise ValueError("--fasta and --query-fai are mutually exclusive (the FASTA supplies lengths)")
    reference_lengths = read_fai_lengths(reference_fai)
    assembly = read_assembly(fasta) if fasta is not None else None
    if assembly is not None:
        query_lengths = {name: contig.length for name, contig in assembly.items()}
        label = "query FASTA"
    else:
        query_lengths = read_fai_lengths(query_fai) if query_fai is not None else None
        label = "query .fai"
    scan = scan_maf(maf, reference_lengths, query_lengths, label)
    stats = compute_stats(
        sample, scan, query_lengths, min_block_bp, overlap_tolerance_bp, assembly, breakpoint_context_bp
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    contigs = select_dotplot_contigs(stats.by_reference, dotplots, dotplot_max)
    plot_paths = write_dotplots(out_dir, sample, scan, contigs, query_lengths)
    for row in stats.by_reference:
        row["dotplot"] = plot_paths.get(str(row["reference_contig"]), "")

    write_tsv(out_dir / f"{sample}.maf_stats.tsv", SUMMARY_COLUMNS, [stats.summary])
    write_tsv(out_dir / f"{sample}.by_reference_contig.tsv", BY_REFERENCE_COLUMNS, stats.by_reference)
    write_tsv(out_dir / f"{sample}.by_query_contig.tsv", BY_QUERY_COLUMNS, stats.by_query)
    write_tsv(out_dir / f"{sample}.breakpoints.tsv", BREAKPOINT_COLUMNS, stats.breakpoints)
    write_tsv(out_dir / f"{sample}.nested_blocks.tsv", NESTED_COLUMNS, stats.nested)
    return stats


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--maf", required=True, help="Pairwise MAF (.maf or .maf.gz); reference row first")
    ap.add_argument("--reference-fai", required=True, help="Reference .fai index")
    ap.add_argument("--sample", required=True, help="Sample name used in output names and rows")
    ap.add_argument("--out-dir", required=True, help="Output directory")
    query = ap.add_mutually_exclusive_group()
    query.add_argument(
        "--fasta",
        default=None,
        help=(
            "Query assembly FASTA (plain or .gz). Adds assembly metrics (contiguity, "
            "N gaps, GC, soft-masking, telomeric ends, unaligned contigs), uses its "
            "lengths for query coverage, and flags breakpoints next to N gaps."
        ),
    )
    query.add_argument(
        "--query-fai",
        default=None,
        help=(
            "Query genome .fai: lengths only. Without this or --fasta, query coverage "
            "is measured against only the query contigs in the MAF, which overstates it."
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
        "--breakpoint-context-bp",
        type=int,
        default=1_000_000,
        help="A breakpoint this close to a contig end (or, with --fasta, an N gap) is classed as contig_end / n_gap (default: 1000000)",
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
            fasta=Path(args.fasta) if args.fasta else None,
            breakpoint_context_bp=args.breakpoint_context_bp,
        )
    except (MafValidationError, ValueError, FileNotFoundError) as exc:
        sys.exit(f"[maf_stats] ERROR: {exc}")
    s = stats.summary
    print(
        f"[maf_stats] {args.sample}: {s['blocks']} blocks, aligned reference "
        f"{_fmt(s['aligned_reference_pct'])}%, identity {_fmt(s['identity_pct'])}%, "
        f"{s['breakpoint_adjacencies']} breakpoint adjacencies",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
