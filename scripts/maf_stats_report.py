#!/usr/bin/env python3
"""Cross-sample MAF QC report.

Merges the per-sample outputs of ``scripts/maf_stats.py`` into one TSV (one row
per sample plus machine-readable ``status``, ``unusual`` and ``flags`` columns),
one combined breakpoint TSV and one combined nested-block TSV, each annotated
with cross-sample recurrence, and one self-contained HTML report (embedded CSS,
inline SVG, and one tiny inline script that opens collapsed sections when an
in-page link points into them).

The HTML leads with a plain-language summary (status counts, cohort ranges,
one card per sample that needs review or is unusual, a compact sample table and
shared features) and puts everything else in collapsed drill-down sections.
"""
from __future__ import annotations

import argparse
import html
import math
import statistics
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

try:
    from scripts.common import open_text
except ModuleNotFoundError:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from scripts.common import open_text


# ── input schemas ────────────────────────────────────────────────────────────

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

# Columns this report adds.
RECURRENCE_COLUMNS = ["recurrence_samples", "recurrent"]
BREAKPOINT_EXTRA_COLUMNS = ["breakpoints_recurrent", "breakpoints_private"]
NESTED_EXTRA_COLUMNS = ["nested_recurrent", "nested_private"]
SUMMARY_EXTRA_COLUMNS = BREAKPOINT_EXTRA_COLUMNS + NESTED_EXTRA_COLUMNS
OUTPUT_SUMMARY_COLUMNS = SUMMARY_COLUMNS + SUMMARY_EXTRA_COLUMNS
OUTPUT_BREAKPOINT_COLUMNS = BREAKPOINT_COLUMNS + RECURRENCE_COLUMNS
OUTPUT_NESTED_COLUMNS = NESTED_COLUMNS + RECURRENCE_COLUMNS

NA_VALUES = frozenset({"", "NA", "NaN", "nan", "None"})

# Free-text names: always shown literally (escaped), never as a muted NA.
NAME_COLUMNS = frozenset(
    {
        "sample",
        "reference_contig",
        "query_contig",
        "left_reference_contig",
        "right_reference_contig",
        "container_reference_contig",
        "dotplot",
    }
)
# Categorical / boolean columns: shown as text; NA shown muted.
CATEGORY_COLUMNS = frozenset(
    {
        "query_length_source",
        "aligned",
        "major",
        "strand",
        "container_strand",
        "telomere_start",
        "telomere_end",
        "location",
        "near_contig_end",
        "near_n_gap",
        "strand_flip",
        "reference_contig_jump",
        "out_of_order",
        "left_strand",
        "right_strand",
        "recurrent",
    }
)
TEXT_COLUMNS = NAME_COLUMNS | CATEGORY_COLUMNS
# Columns rendered with two decimals; all other numeric columns are integers.
FLOAT_COLUMNS = frozenset(
    {
        "aligned_reference_pct",
        "block_span_reference_pct",
        "aligned_query_pct",
        "block_span_query_pct",
        "identity_pct",
        "insertion_column_pct",
        "deletion_column_pct",
        "breakpoints_per_gb_aligned",
        "assembly_n_pct",
        "assembly_gc_pct",
        "assembly_softmasked_pct",
        "unaligned_contig_n_pct",
        "unaligned_contig_softmasked_pct",
        "gc_pct",
        "softmasked_pct",
    }
)

# ── flagging configuration (tune here) ───────────────────────────────────────

MAD_TO_SD = 1.4826
DEFAULT_MIN_SAMPLES = 5
DEFAULT_Z_THRESHOLD = 3.5
DEFAULT_RECURRENCE_WINDOW_BP = 500_000

# Direction in which a metric is "bad". Order here is the order flags are listed,
# and each gets a strip plot. overlapping_reference_bp and reference_contig_jumps
# are deliberately absent: AnchorWave aligns each reference chromosome one-to-one
# with its counterpart, so both are 0 for every sample. They are reported but
# neither flagged nor plotted.
FLAG_DIRECTIONS: dict[str, str] = {
    "aligned_reference_pct": "low",
    "aligned_query_pct": "low",
    "identity_pct": "low",
    "assembly_contig_n50_bp": "low",
    "overlapping_query_bp": "high",
    "breakpoints_per_gb_aligned": "high",
    "breakpoints_private": "high",
    "assembly_n_pct": "high",
    "unaligned_contig_bp": "high",
}

# Absolute lower bound on the robust scale, in the metric's own units. Prevents
# division by zero and flags on trivially small deviations in tight cohorts.
SCALE_ABSOLUTE_FLOORS: dict[str, float] = {
    "aligned_reference_pct": 1.0,
    "aligned_query_pct": 1.0,
    "identity_pct": 1.0,
    "assembly_n_pct": 0.5,
    "breakpoints_private": 2.0,
    "breakpoints_per_gb_aligned": 1.0,
    "overlapping_query_bp": 100_000.0,
    "assembly_contig_n50_bp": 1_000_000.0,
    "unaligned_contig_bp": 1_000_000.0,
}

# Additional scale floor as a fraction of |median| for count-like metrics.
SCALE_RELATIVE_FLOORS: dict[str, float] = {
    "breakpoints_per_gb_aligned": 0.05,
    "overlapping_query_bp": 0.05,
    "assembly_contig_n50_bp": 0.05,
    "unaligned_contig_bp": 0.05,
}

# CLI option dest -> metric for optional absolute thresholds. Direction comes
# from FLAG_DIRECTIONS.
THRESHOLD_OPTIONS: dict[str, str] = {
    "flag_min_aligned_reference": "aligned_reference_pct",
    "flag_min_identity": "identity_pct",
    "flag_max_breakpoints_per_gb": "breakpoints_per_gb_aligned",
}

# Pseudo-column for the assembly table: "major telomeric ends / 2 x major sequences".
TELOMERE_DISPLAY_COLUMN = "assembly_major_telomeres"

METRIC_LABELS: dict[str, str] = {
    "sample": "Sample",
    "reference_length_bp": "Reference length (bp)",
    "aligned_reference_bp": "Aligned ref bp",
    "aligned_reference_pct": "Aligned ref (%)",
    "block_span_reference_bp": "Block-span ref bp",
    "block_span_reference_pct": "Block-span ref (%)",
    "query_length_bp": "Query length (bp)",
    "query_length_source": "Query length source",
    "aligned_query_bp": "Aligned query bp",
    "aligned_query_pct": "Aligned query (%)",
    "block_span_query_bp": "Block-span query bp",
    "block_span_query_pct": "Block-span query (%)",
    "identity_matches": "Identity matches",
    "identity_compared_columns": "Identity compared columns",
    "identity_pct": "Identity (%)",
    "alignment_columns": "Alignment columns",
    "insertion_columns": "Insertion columns",
    "insertion_column_pct": "Insertion cols (%)",
    "deletion_columns": "Deletion columns",
    "deletion_column_pct": "Deletion cols (%)",
    "query_n_bases_in_blocks": "Query N in blocks",
    "blocks": "Blocks",
    "blocks_considered": "Blocks considered",
    "overlapping_reference_bp": "Overlapping ref bp",
    "overlapping_query_bp": "Overlapping query bp",
    "nested_blocks": "Nested blocks",
    "nested_block_reference_bp": "Nested block ref bp",
    "nested_recurrent": "Recurrent nested",
    "nested_private": "Private nested",
    "aligned_query_contigs": "Aligned query contigs",
    "aligned_query_contig_n50_bp": "Aligned query contig N50 (bp)",
    "strand_flips": "Strand flips",
    "reference_contig_jumps": "Ref-contig jumps",
    "out_of_order_adjacencies": "Out-of-order adj.",
    "breakpoint_adjacencies": "Breakpoint adj.",
    "breakpoints_near_contig_end": "Breakpoints near contig end",
    "breakpoints_near_n_gap": "Breakpoints near N gap",
    "breakpoints_interior": "Interior breakpoints",
    "breakpoints_per_gb_aligned": "Breakpoints / Gb aligned",
    "breakpoints_recurrent": "Recurrent breakpoints",
    "breakpoints_private": "Private breakpoints",
    "assembly_length_bp": "Assembly length (bp)",
    "assembly_sequences": "Sequences",
    "assembly_scaffold_n50_bp": "Scaffold N50 (bp)",
    "assembly_scaffold_l50": "Scaffold L50",
    "assembly_largest_sequence_bp": "Largest sequence (bp)",
    "assembly_contig_pieces": "Contig pieces",
    "assembly_contig_n50_bp": "Contig N50 (bp)",
    "assembly_contig_l50": "Contig L50",
    "assembly_n_bp": "N bp",
    "assembly_n_pct": "N (%)",
    "assembly_n_gaps": "N gaps",
    "assembly_gc_pct": "GC (%)",
    "assembly_softmasked_pct": "Soft-masked (%)",
    "assembly_major_sequences": "Major sequences",
    "assembly_major_telomeric_ends": "Major telomeric ends",
    "assembly_minor_sequences_with_telomere": "Minor seqs with telomere",
    TELOMERE_DISPLAY_COLUMN: "Telomeric ends (major)",
    "unaligned_contigs": "Unaligned contigs",
    "unaligned_contig_bp": "Unaligned contig bp",
    "unaligned_contig_n_pct": "Unaligned contig N (%)",
    "unaligned_contig_softmasked_pct": "Unaligned contig soft-masked (%)",
    "min_block_bp": "Min block bp",
    "overlap_tolerance_bp": "Overlap tolerance bp",
    "breakpoint_context_bp": "Breakpoint context bp",
    "reference_contig": "Reference contig",
    "query_contig": "Query contig",
    "query_contigs": "Query contigs",
    "reference_contigs": "Ref contigs",
    "aligned": "Aligned",
    "unaligned_start_bp": "Unaligned start bp",
    "unaligned_end_bp": "Unaligned end bp",
    "n_bp": "N bp",
    "n_gaps": "N gaps",
    "contig_pieces": "Contig pieces",
    "major": "Major",
    "gc_pct": "GC (%)",
    "softmasked_pct": "Soft-masked (%)",
    "telomere_start": "Telomere start",
    "telomere_end": "Telomere end",
    "dotplot": "Dotplot",
    "query_contig_length_bp": "Query contig length (bp)",
    "query_junction_start": "Junction start",
    "query_junction_end": "Junction end",
    "distance_to_contig_end_bp": "Distance to contig end (bp)",
    "location": "Location",
    "recurrence_samples": "Other samples sharing",
    "recurrent": "Recurrent",
    "query_start": "Query start",
    "query_end": "Query end",
    "strand": "Strand",
    "reference_start": "Reference start",
    "reference_end": "Reference end",
    "container_query_start": "Container query start",
    "container_query_end": "Container query end",
    "container_reference_contig": "Container reference contig",
    "container_reference_start": "Container reference start",
    "container_reference_end": "Container reference end",
    "container_strand": "Container strand",
}

PRIMARY_TABLE_COLUMNS = [
    "aligned_reference_pct",
    "aligned_query_pct",
    "identity_pct",
    "insertion_column_pct",
    "deletion_column_pct",
    "strand_flips",
    "out_of_order_adjacencies",
    "reference_contig_jumps",
    "breakpoint_adjacencies",
    "breakpoints_recurrent",
    "breakpoints_private",
    "nested_blocks",
    "nested_recurrent",
    "nested_private",
]

ASSEMBLY_TABLE_COLUMNS = [
    "assembly_length_bp",
    "assembly_sequences",
    "assembly_scaffold_n50_bp",
    "assembly_contig_pieces",
    "assembly_contig_n50_bp",
    "assembly_n_gaps",
    "assembly_n_pct",
    "assembly_gc_pct",
    "assembly_softmasked_pct",
    TELOMERE_DISPLAY_COLUMN,
    "assembly_minor_sequences_with_telomere",
    "unaligned_contigs",
    "unaligned_contig_bp",
]

MAF_CONTIGUITY_COLUMNS = ["aligned_query_contigs", "aligned_query_contig_n50_bp"]

REFERENCE_TABLE_COLUMNS = [
    "reference_contig",
    "reference_length_bp",
    "aligned_reference_bp",
    "aligned_reference_pct",
    "block_span_reference_pct",
    "identity_pct",
    "insertion_column_pct",
    "deletion_column_pct",
    "blocks",
    "query_contigs",
    "overlapping_reference_bp",
    "breakpoint_adjacencies",
]

QUERY_TABLE_COLUMNS = [
    "query_contig",
    "query_length_bp",
    "major",
    "aligned",
    "aligned_query_bp",
    "block_span_query_pct",
    "overlapping_query_bp",
    "unaligned_start_bp",
    "unaligned_end_bp",
    "blocks",
    "nested_blocks",
    "reference_contigs",
    "strand_flips",
    "reference_contig_jumps",
    "out_of_order_adjacencies",
    "breakpoint_adjacencies",
    "n_bp",
    "contig_pieces",
    "gc_pct",
    "softmasked_pct",
    "telomere_start",
    "telomere_end",
]

METRIC_DEFINITIONS: list[tuple[str, str]] = [
    ("reference_length_bp", "Total reference length from the reference <code>.fai</code>."),
    ("aligned_reference_bp",
     "Reference bases in alignment columns where <em>both</em> rows have a base (no "
     "<code>-</code>), counted once per reference position. This is the honest measure of "
     "how much of the reference is opposite query sequence."),
    ("aligned_reference_pct", "Aligned reference bp divided by reference length."),
    ("block_span_reference_bp",
     "Union of reference intervals spanned by accepted blocks. Includes long indels inside "
     "blocks: AnchorWave blocks can span whole chromosome arms between collinear anchors, "
     "so block-span coverage <strong>overstates</strong> alignment. Shown for comparison only."),
    ("block_span_reference_pct", "Block-span reference bp divided by reference length."),
    ("query_length_bp", "Query genome size used as the query-percentage denominator."),
    ("query_length_source",
     "Source of the query denominator. <code>query_fasta</code>: the query assembly given "
     "with <code>--fasta</code>. <code>query_fai</code>: all contigs in the query "
     "<code>.fai</code>. <code>maf_srcsize</code>: only query contigs present in the MAF "
     "(their <code>srcSize</code>) are counted, which <strong>overstates</strong> query "
     "percentages because unaligned contigs are absent; such values are marked &dagger; and "
     "are not comparable to the other sources."),
    ("aligned_query_bp", "Query bases in columns where both rows have a base."),
    ("aligned_query_pct", "Aligned query bp divided by query length (see denominator source)."),
    ("block_span_query_bp", "Union of forward-coordinate query intervals spanned by accepted blocks."),
    ("block_span_query_pct", "Block-span query bp divided by query length."),
    ("identity_pct",
     "Matching columns divided by columns where both bases are A/C/G/T (case-insensitive); "
     "N and other ambiguity codes are excluded. Raw block-column metric."),
    ("alignment_columns", "All alignment columns across accepted blocks."),
    ("insertion_columns", "Columns with a query base opposite a reference gap."),
    ("insertion_column_pct", "Insertion columns divided by alignment columns."),
    ("deletion_columns", "Columns with a reference base opposite a query gap."),
    ("deletion_column_pct", "Deletion columns divided by alignment columns."),
    ("query_n_bases_in_blocks", "Query N bases inside accepted blocks."),
    ("blocks", "Count of accepted pairwise MAF blocks. Raw block metric."),
    ("overlapping_reference_bp",
     "Reference bp covered by more than one block. Makes duplicated or secondary alignment "
     "content visible."),
    ("overlapping_query_bp",
     "Query bp (forward coordinates) covered by more than one block: the same query "
     "sequence aligned to more than one reference location."),
    ("nested_blocks",
     "Blocks whose forward query interval lies inside a larger block's query interval: the "
     "same query sequence aligned in two places (a secondary or transposed alignment). "
     "Nested blocks are <strong>excluded from the breakpoint walk</strong> (walking them as "
     "adjacencies produced spurious strand flips and very long junctions) and listed "
     "separately in <code>&lt;sample&gt;.nested_blocks.tsv</code>. Reported, not flagged: "
     "they may be biology (duplications, transpositions) or artefacts."),
    ("nested_block_reference_bp", "Total reference span of nested blocks."),
    ("nested_recurrent",
     "This sample's nested blocks matched by a nested block in at least one other sample "
     "(same reference contig, start and end each within the recurrence window). NA with "
     "fewer than two samples."),
    ("nested_private",
     "This sample's nested blocks not seen in any other sample. NA with fewer than two "
     "samples."),
    ("aligned_query_contigs",
     "Query contigs with at least one accepted block. From the MAF only, so unaligned "
     "contigs are invisible."),
    ("aligned_query_contig_n50_bp",
     "N50 of the lengths of aligned query sequences. A <strong>lower bound on "
     "fragmentation</strong>: unaligned (usually small) sequences are missing. For "
     "chromosome-scale assemblies it is essentially a chromosome length and says little "
     "about contiguity; prefer <code>assembly_contig_n50_bp</code>."),
    ("strand_flips", "Adjacent blocks along a query contig whose query strands differ."),
    ("reference_contig_jumps",
     "Adjacent blocks along a query contig aligned to different reference contigs."),
    ("out_of_order_adjacencies",
     "Same-reference, same-strand adjacent blocks whose reference coordinates contradict "
     "query order beyond <code>overlap_tolerance_bp</code>."),
    ("breakpoint_adjacencies",
     "Unique adjacent block pairs having any breakpoint property (strand flip, "
     "reference-contig jump, out-of-order). A pair with several properties counts once."),
    ("breakpoints_near_contig_end",
     "Breakpoints whose junction lies within <code>breakpoint_context_bp</code> of a query "
     "contig end."),
    ("breakpoints_near_n_gap",
     "Breakpoints whose junction lies within <code>breakpoint_context_bp</code> of an N gap "
     "in the query. Needs <code>--fasta</code>."),
    ("breakpoints_interior", "Breakpoints neither near a contig end nor near an N gap."),
    ("breakpoints_per_gb_aligned",
     "Breakpoint adjacencies divided by <code>aligned_reference_bp</code>, times 1e9."),
    ("breakpoints_recurrent",
     "This sample's breakpoints matched by a breakpoint in at least one other sample (see "
     "Breakpoints). NA with fewer than two samples."),
    ("breakpoints_private",
     "This sample's breakpoints not seen in any other sample. NA with fewer than two samples."),
    ("assembly_length_bp",
     "Total length of the query assembly FASTA. All <code>assembly_*</code> and "
     "<code>unaligned_*</code> metrics need <code>--fasta</code> and are NA without it."),
    ("assembly_sequences", "Number of FASTA sequences (scaffolds/chromosomes) in the query assembly."),
    ("assembly_scaffold_n50_bp",
     "N50 of the FASTA sequence lengths (scaffold N50). For a chromosome-scale assembly "
     "this is a chromosome length."),
    ("assembly_scaffold_l50", "Number of sequences needed to reach the scaffold N50."),
    ("assembly_largest_sequence_bp", "Length of the longest FASTA sequence."),
    ("assembly_contig_pieces",
     "Contigs after splitting every sequence at N runs of at least 10 bp: maximal "
     "stretches between such gaps."),
    ("assembly_contig_n50_bp",
     "True <strong>contig</strong> N50: N50 of the contig pieces (sequences split at N runs "
     "&ge; 10 bp). Unlike scaffold N50 this reflects contiguity of the sequence itself."),
    ("assembly_contig_l50", "Number of contig pieces needed to reach the contig N50."),
    ("assembly_n_pct", "N bases divided by assembly length."),
    ("assembly_n_gaps",
     "Runs of at least 10 N (scaffold gaps). Isolated or short N runs are not counted."),
    ("assembly_gc_pct", "G+C over non-N bases."),
    ("assembly_softmasked_pct",
     "Lowercase (soft-masked) bases divided by assembly length. NA when the FASTA has no "
     "lowercase at all (not soft-masked)."),
    ("assembly_major_sequences",
     "Sequences at least 10% as long as the longest one, i.e. the chromosome-scale "
     "sequences of a chromosome-level assembly (unplaced scaffolds are excluded)."),
    ("assembly_major_telomeric_ends",
     "Telomeric ends among the major sequences, out of 2 &times; "
     "<code>assembly_major_sequences</code> (shown as &ldquo;ends / possible&rdquo;). An end "
     "is telomeric when at least 50% of its terminal 1 kb is plant telomere repeat "
     "(<code>TTTAGGG</code>/<code>CCCTAAA</code>)."),
    ("assembly_minor_sequences_with_telomere",
     "Sequences outside the major set with at least one telomeric end (e.g. unplaced "
     "telomeric fragments)."),
    ("unaligned_contigs", "Assembly contigs with no accepted block."),
    ("unaligned_contig_bp", "Total length of unaligned contigs."),
    ("unaligned_contig_n_pct", "N content of unaligned contigs."),
    ("unaligned_contig_softmasked_pct",
     "Soft-masked fraction of unaligned contigs; high values suggest repeats. NA when the "
     "FASTA has no lowercase."),
    ("min_block_bp",
     "Blocks shorter than this are excluded from breakpoint classification only."),
    ("overlap_tolerance_bp",
     "Reference overlap/backtrack allowed between adjacent same-strand blocks before the "
     "pair is classified out-of-order."),
    ("breakpoint_context_bp",
     "Distance used to decide whether a breakpoint junction is near a contig end or an N gap "
     "(<code>maf_stats.py</code> default 1,000,000)."),
    ("recurrence_samples",
     "Per breakpoint or nested block: number of <em>other</em> samples with a matching one."),
    ("contig_pieces",
     "Per query contig: pieces after splitting at N runs &ge; 10 bp (1 = no gap)."),
    ("major",
     "Per query contig: whether it is one of the assembly's major (chromosome-scale) "
     "sequences."),
    ("dotplot",
     "Structural dotplot of a reference contig against the query contigs aligned to it "
     "(forward blocks blue, reverse blocks red)."),
]


# ── reading ──────────────────────────────────────────────────────────────────


def parse_number(raw: str | None) -> float | None:
    """Parse a TSV numeric cell; ``NA``/empty/non-finite/garbage -> ``None``."""
    if raw is None:
        return None
    text = raw.strip()
    if text in NA_VALUES:
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    if not math.isfinite(value):
        return None
    return value


def read_tsv(path: Path, required: list[str]) -> list[dict[str, str]]:
    """Read a headered TSV into dicts, checking that ``required`` columns exist."""
    with open_text(path, "rt") as handle:
        header_line = handle.readline()
        if not header_line.strip():
            raise ValueError(f"{path}: empty file (missing header row)")
        header = header_line.rstrip("\r\n").split("\t")
        missing = [c for c in required if c not in header]
        if missing:
            raise ValueError(f"{path}: missing required column(s): {', '.join(missing)}")
        rows: list[dict[str, str]] = []
        for line_number, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) != len(header):
                raise ValueError(
                    f"{path}:{line_number}: expected {len(header)} columns, found {len(parts)}"
                )
            rows.append(dict(zip(header, parts)))
    return rows


def read_summary(path: Path) -> dict[str, str]:
    """Read a per-sample ``<sample>.maf_stats.tsv`` (header + exactly one row)."""
    rows = read_tsv(path, SUMMARY_COLUMNS)
    if len(rows) != 1:
        raise ValueError(f"{path}: expected exactly 1 data row, found {len(rows)}")
    row = rows[0]
    if not row["sample"].strip():
        raise ValueError(f"{path}: empty sample name")
    return {c: row[c] for c in SUMMARY_COLUMNS}


def read_contig_tables(
    paths: list[Path],
    required: list[str],
    samples: list[str],
    kind: str,
) -> dict[str, list[dict[str, str]]]:
    """Group rows of per-sample TSVs by their ``sample`` column.

    Raises if a row names a sample absent from the summaries, or if two files
    contribute rows for the same sample.
    """
    known = set(samples)
    grouped: dict[str, list[dict[str, str]]] = {}
    source: dict[str, Path] = {}
    for path in paths:
        rows = read_tsv(path, required)
        for row in rows:
            sample = row["sample"]
            if sample not in known:
                raise ValueError(
                    f"{path}: {kind} row for sample {sample!r}, which is not among the "
                    f"summary inputs"
                )
            if sample in source and source[sample] != path:
                raise ValueError(
                    f"{path}: duplicate {kind} input for sample {sample!r} "
                    f"(also in {source[sample]})"
                )
            source[sample] = path
            grouped.setdefault(sample, []).append({c: row[c] for c in required})
    return grouped


# ── breakpoint recurrence ────────────────────────────────────────────────────


def _breakpoint_ends(row: dict[str, str]) -> tuple[tuple[str, float], tuple[str, float]] | None:
    """``((left_contig, left_pos), (right_contig, right_pos))`` or None if incomplete."""
    lc = (row.get("left_reference_contig") or "").strip()
    rc = (row.get("right_reference_contig") or "").strip()
    lp = parse_number(row.get("left_reference_pos"))
    rp = parse_number(row.get("right_reference_pos"))
    if lc in NA_VALUES or rc in NA_VALUES or lp is None or rp is None:
        return None
    return (lc, lp), (rc, rp)


def _ends_match(
    a: tuple[tuple[str, float], tuple[str, float]],
    b: tuple[tuple[str, float], tuple[str, float]],
    window_bp: float,
) -> bool:
    def near(x: tuple[str, float], y: tuple[str, float]) -> bool:
        return x[0] == y[0] and abs(x[1] - y[1]) <= window_bp

    return (near(a[0], b[0]) and near(a[1], b[1])) or (near(a[0], b[1]) and near(a[1], b[0]))


def annotate_recurrence(
    samples: list[str],
    breakpoints: dict[str, list[dict[str, str]]],
    window_bp: int = DEFAULT_RECURRENCE_WINDOW_BP,
) -> list[dict[str, str]]:
    """Flatten breakpoints (sample input order) and add recurrence columns.

    A breakpoint in sample A matches one in sample B (B != A) when their
    unordered junction-end pairs match in either orientation, each paired end on
    the same reference contig and within ``window_bp``. ``recurrence_samples``
    counts the distinct other samples with at least one match. With fewer than
    two samples both added columns are ``NA``.

    Breakpoints are bucketed by unordered reference-contig pair and, inside a
    bucket, indexed by ``position // window`` of each end, so only candidates in
    neighbouring bins are compared.
    """
    flat: list[dict[str, str]] = []
    for sample in samples:
        for row in breakpoints.get(sample, []):
            out = {c: row.get(c, "") for c in BREAKPOINT_COLUMNS}
            flat.append(out)

    if len(samples) < 2:
        for row in flat:
            row["recurrence_samples"] = "NA"
            row["recurrent"] = "NA"
        return flat

    bin_width = max(int(window_bp), 1)
    ends = [_breakpoint_ends(r) for r in flat]
    # bucket key (sorted contig pair) -> (contig, bin) -> indices of breakpoints
    buckets: dict[tuple[str, str], dict[tuple[str, int], list[int]]] = {}
    for idx, e in enumerate(ends):
        if e is None:
            continue
        key = tuple(sorted((e[0][0], e[1][0])))
        index = buckets.setdefault(key, {})  # type: ignore[arg-type]
        for contig, pos in e:
            slot = (contig, int(pos // bin_width))
            lst = index.setdefault(slot, [])
            if not lst or lst[-1] != idx:
                lst.append(idx)

    for idx, e in enumerate(ends):
        matched: set[str] = set()
        if e is not None:
            key = tuple(sorted((e[0][0], e[1][0])))
            index = buckets[key]  # type: ignore[index]
            contig, pos = e[0]
            b = int(pos // bin_width)
            sample = flat[idx]["sample"]
            seen: set[int] = set()
            for nb in (b - 1, b, b + 1):
                for j in index.get((contig, nb), ()):
                    if j in seen or j == idx:
                        continue
                    seen.add(j)
                    other = flat[j]["sample"]
                    if other == sample or other in matched:
                        continue
                    if _ends_match(e, ends[j], window_bp):  # type: ignore[arg-type]
                        matched.add(other)
        flat[idx]["recurrence_samples"] = str(len(matched))
        flat[idx]["recurrent"] = "true" if matched else "false"
    return flat


def recurrence_counts(
    samples: list[str],
    annotated: list[dict[str, str]],
    columns: tuple[str, str] | list[str] = tuple(BREAKPOINT_EXTRA_COLUMNS),
) -> dict[str, dict[str, str]]:
    """Per sample ``{recurrent_col: n, private_col: n}`` as strings.

    ``columns`` defaults to ``("breakpoints_recurrent", "breakpoints_private")``;
    pass ``NESTED_EXTRA_COLUMNS`` for nested blocks. NA with fewer than two samples.
    """
    rec_col, priv_col = columns
    if len(samples) < 2:
        return {s: {rec_col: "NA", priv_col: "NA"} for s in samples}
    counts = {s: [0, 0] for s in samples}
    for row in annotated:
        counts[row["sample"]][0 if row["recurrent"] == "true" else 1] += 1
    return {s: {rec_col: str(r), priv_col: str(p)} for s, (r, p) in counts.items()}


def _nested_key(row: dict[str, str]) -> tuple[str, float, float] | None:
    """``(reference_contig, reference_start, reference_end)`` or None if incomplete."""
    contig = (row.get("reference_contig") or "").strip()
    start = parse_number(row.get("reference_start"))
    end = parse_number(row.get("reference_end"))
    if contig in NA_VALUES or start is None or end is None:
        return None
    return contig, start, end


def annotate_nested_recurrence(
    samples: list[str],
    nested: dict[str, list[dict[str, str]]],
    window_bp: int = DEFAULT_RECURRENCE_WINDOW_BP,
) -> list[dict[str, str]]:
    """Flatten nested blocks (sample input order) and add recurrence columns.

    A nested block in sample A matches one in sample B (B != A) when both lie
    on the same reference contig and their ``reference_start`` and
    ``reference_end`` each differ by at most ``window_bp``.
    ``recurrence_samples`` counts distinct other samples with at least one
    match; both added columns are ``NA`` with fewer than two samples.
    Candidates are indexed by ``(contig, reference_start // window)``.
    """
    flat: list[dict[str, str]] = []
    for sample in samples:
        for row in nested.get(sample, []):
            flat.append({c: row.get(c, "") for c in NESTED_COLUMNS})

    if len(samples) < 2:
        for row in flat:
            row["recurrence_samples"] = "NA"
            row["recurrent"] = "NA"
        return flat

    bin_width = max(int(window_bp), 1)
    keys = [_nested_key(r) for r in flat]
    index: dict[tuple[str, int], list[int]] = {}
    for idx, k in enumerate(keys):
        if k is not None:
            index.setdefault((k[0], int(k[1] // bin_width)), []).append(idx)

    for idx, k in enumerate(keys):
        matched: set[str] = set()
        if k is not None:
            contig, start, end = k
            sample = flat[idx]["sample"]
            b = int(start // bin_width)
            for nb in (b - 1, b, b + 1):
                for j in index.get((contig, nb), ()):
                    other = flat[j]["sample"]
                    if j == idx or other == sample or other in matched:
                        continue
                    _, o_start, o_end = keys[j]  # type: ignore[misc]
                    if abs(o_start - start) <= window_bp and abs(o_end - end) <= window_bp:
                        matched.add(other)
        flat[idx]["recurrence_samples"] = str(len(matched))
        flat[idx]["recurrent"] = "true" if matched else "false"
    return flat


# ── flagging ─────────────────────────────────────────────────────────────────


@dataclass
class MetricStats:
    n: int
    median: float
    mad: float
    scale: float
    relative_enabled: bool


STATUS_OK = "ok"
STATUS_REVIEW = "review"
STATUS_UNUSUAL = "unusual"
STATUS_ICONS = {STATUS_OK: "\u2714", STATUS_REVIEW: "\u26a0", STATUS_UNUSUAL: "\u25c6"}
STATUS_LABELS = {STATUS_OK: "normal", STATUS_REVIEW: "needs review", STATUS_UNUSUAL: "unusual"}


@dataclass
class FlagResult:
    # per sample (input order): metric -> list of reasons (bad direction or
    # absolute-threshold failures; these make a sample "review")
    flags: list[dict[str, list[str]]]
    stats: dict[str, MetricStats | None] = field(default_factory=dict)
    # metrics whose relative flagging was skipped -> number of non-NA samples
    skipped: dict[str, int] = field(default_factory=dict)
    # per sample: metric -> reasons for |z| > threshold in the GOOD direction
    unusual: list[dict[str, list[str]]] = field(default_factory=list)
    # per sample: metric -> robust z (only metrics with relative flagging enabled)
    z_scores: list[dict[str, float]] = field(default_factory=list)

    @property
    def status(self) -> list[str]:
        return compute_status(self)


def robust_scale(metric: str, median: float, mad: float) -> float:
    """``max(1.4826*MAD, absolute floor, relative floor * |median|)``."""
    scale = MAD_TO_SD * mad
    scale = max(scale, SCALE_ABSOLUTE_FLOORS.get(metric, 0.0))
    scale = max(scale, SCALE_RELATIVE_FLOORS.get(metric, 0.0) * abs(median))
    if not scale > 0:
        # Defensive: every flaggable metric has a positive absolute floor, but
        # never divide by zero if a new metric is added without one.
        scale = 1.0
    return scale


def compute_flags(
    summaries: list[dict[str, str]],
    *,
    min_samples: int = DEFAULT_MIN_SAMPLES,
    z_threshold: float = DEFAULT_Z_THRESHOLD,
    thresholds: dict[str, float] | None = None,
) -> FlagResult:
    """Compute relative (robust-z) and absolute flags for each sample.

    ``thresholds`` maps metric name -> absolute threshold; the direction is
    taken from ``FLAG_DIRECTIONS`` (low-is-bad -> flag when value < threshold).
    Metrics that are NA (or absent) for a sample are skipped for that sample.

    |z| above ``z_threshold`` in the bad direction goes in ``flags``; in the
    good direction it goes in ``unusual`` (never in ``flags``). Absolute
    thresholds only ever produce flags.
    """
    thresholds = thresholds or {}
    for metric in thresholds:
        if metric not in FLAG_DIRECTIONS:
            raise ValueError(f"No flag direction defined for threshold metric {metric!r}")
    result = FlagResult(
        flags=[{} for _ in summaries],
        unusual=[{} for _ in summaries],
        z_scores=[{} for _ in summaries],
    )

    for metric, direction in FLAG_DIRECTIONS.items():
        values = [parse_number(s.get(metric)) for s in summaries]
        present = [v for v in values if v is not None]

        stats: MetricStats | None = None
        if present:
            med = statistics.median(present)
            mad = statistics.median(abs(v - med) for v in present)
            stats = MetricStats(
                n=len(present),
                median=med,
                mad=mad,
                scale=robust_scale(metric, med, mad),
                relative_enabled=len(present) >= min_samples,
            )
        result.stats[metric] = stats
        if stats is None or not stats.relative_enabled:
            result.skipped[metric] = len(present)

        threshold = thresholds.get(metric)
        for idx, value in enumerate(values):
            if value is None:
                continue
            reasons: list[str] = []
            if stats is not None and stats.relative_enabled:
                z = (value - stats.median) / stats.scale
                result.z_scores[idx][metric] = z
                tag = "(mad0_floor)" if stats.mad == 0 else ""
                bad = z < -z_threshold if direction == "low" else z > z_threshold
                good = z > z_threshold if direction == "low" else z < -z_threshold
                if bad:
                    reasons.append(f"robust_z={z:.2f}{tag}")
                elif good:
                    result.unusual[idx][metric] = [f"robust_z={z:.2f}{tag}"]
            if threshold is not None:
                if direction == "low" and value < threshold:
                    reasons.append(f"below_threshold_{threshold:g}")
                elif direction == "high" and value > threshold:
                    reasons.append(f"above_threshold_{threshold:g}")
            if reasons:
                result.flags[idx][metric] = reasons
    return result


def sample_status(
    sample_flags: dict[str, list[str]], sample_unusual: dict[str, list[str]] | None = None
) -> str:
    """``review`` if any flag, else ``unusual`` if any unusual metric, else ``ok``."""
    if sample_flags:
        return STATUS_REVIEW
    if sample_unusual:
        return STATUS_UNUSUAL
    return STATUS_OK


def compute_status(flag_result: FlagResult) -> list[str]:
    """Per-sample status (input order): review > unusual > ok."""
    unusual = flag_result.unusual or [{} for _ in flag_result.flags]
    return [sample_status(f, u) for f, u in zip(flag_result.flags, unusual)]


def format_flags(sample_flags: dict[str, list[str]]) -> str:
    """``metric:reason;metric:reason`` in FLAG_DIRECTIONS order; '' when none.

    Also used for the ``unusual`` column.
    """
    entries: list[str] = []
    for metric in FLAG_DIRECTIONS:
        for reason in sample_flags.get(metric, []):
            entries.append(f"{metric}:{reason}")
    return ";".join(entries)


# ── output TSVs ──────────────────────────────────────────────────────────────


STATUS_COLUMNS = ["status", "unusual", "flags"]


def write_summary_tsv(
    path: Path,
    summaries: list[dict[str, str]],
    flags: list[dict[str, list[str]]],
    unusual: list[dict[str, list[str]]] | None = None,
) -> None:
    """Write summary rows (input + recurrence columns) plus status/unusual/flags."""
    unusual = unusual or [{} for _ in summaries]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(OUTPUT_SUMMARY_COLUMNS + STATUS_COLUMNS) + "\n")
        for row, sample_flags, sample_unusual in zip(summaries, flags, unusual):
            fh.write(
                "\t".join([row.get(c, "NA") for c in OUTPUT_SUMMARY_COLUMNS]
                          + [sample_status(sample_flags, sample_unusual),
                             format_flags(sample_unusual),
                             format_flags(sample_flags)])
                + "\n"
            )


def _write_rows(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(c, "NA") for c in columns) + "\n")


def write_breakpoints_tsv(path: Path, annotated: list[dict[str, str]]) -> None:
    _write_rows(path, OUTPUT_BREAKPOINT_COLUMNS, annotated)


def write_nested_tsv(path: Path, annotated: list[dict[str, str]]) -> None:
    _write_rows(path, OUTPUT_NESTED_COLUMNS, annotated)


# ── plain-language reasons ───────────────────────────────────────────────────


def fmt_bp(value: float) -> str:
    """Human-readable length: ``812 bp``, ``45.3 kb``, ``48.1 Mb``, ``2.13 Gb``."""
    a = abs(value)
    if a >= 1e9:
        return f"{value / 1e9:.2f} Gb"
    if a >= 1e6:
        return f"{value / 1e6:.1f} Mb"
    if a >= 1e3:
        return f"{value / 1e3:.1f} kb"
    return f"{value:,.0f} bp"


def _scaled(value: float, unit: float, decimals: int) -> str:
    return f"{value / unit:,.{decimals}f}"


def fmt_span(contig: str, start: float, end: float) -> str:
    """``chr2:122.0–148.0 Mb`` (enough decimals that start and end differ)."""
    top = max(abs(start), abs(end))
    if top >= 1e6:
        unit, suffix = 1e6, " Mb"
    elif top >= 1e3:
        unit, suffix = 1e3, " kb"
    else:
        return f"{contig}:{start:,.0f}–{end:,.0f}"
    for decimals in (1, 2, 3):
        a, b = _scaled(start, unit, decimals), _scaled(end, unit, decimals)
        if a != b:
            break
    return f"{contig}:{a}–{b}{suffix}"


def _fmt_pos(value: float) -> str:
    if abs(value) >= 1e6:
        return f"{value / 1e6:,.2f} Mb"
    if abs(value) >= 1e3:
        return f"{value / 1e3:,.1f} kb"
    return f"{value:,.0f} bp"


def _fmt_metric_value(metric: str, value: float) -> str:
    if metric in ("aligned_reference_pct", "aligned_query_pct"):
        return f"{value:.1f}%"
    if metric in ("identity_pct", "assembly_n_pct"):
        return f"{value:.2f}%"
    if metric.endswith("_bp"):
        return fmt_bp(value)
    if metric == "breakpoints_per_gb_aligned":
        return f"{value:.1f}"
    return f"{value:g}"


# metric -> (subject with {v}, comparison with {m})
SENTENCE_TEMPLATES: dict[str, tuple[str, str]] = {
    "aligned_reference_pct": ("Aligned reference {v}", "vs cohort median {m}"),
    "aligned_query_pct": ("Aligned query {v}", "vs cohort median {m}"),
    "identity_pct": ("Identity {v}", "vs cohort median {m}"),
    "assembly_contig_n50_bp": ("Contig N50 {v}", "vs cohort median {m}"),
    "overlapping_query_bp": ("{v} of query sequence aligned more than once", "vs median {m}"),
    "breakpoints_per_gb_aligned": ("{v} breakpoints per Gb aligned", "vs median {m}"),
    "breakpoints_private": ("{v} private breakpoints (not seen in any other sample)",
                            "vs median {m}"),
    "assembly_n_pct": ("N content {v}", "vs cohort median {m}"),
    "unaligned_contig_bp": ("{v} in assembly contigs with no alignment", "vs median {m}"),
}

THRESHOLD_OPTION_FOR_METRIC = {
    metric: "--" + dest.replace("_", "-") for dest, metric in THRESHOLD_OPTIONS.items()
}


@dataclass
class Reason:
    """One plain-language reason a sample needs review or is unusual."""

    metric: str
    kind: str  # STATUS_REVIEW or STATUS_UNUSUAL
    sentence: str  # plain text (escape when rendering)
    where: str = ""  # plain text pointer to the region(s), or ""
    dotplots: list[tuple[str, str]] = field(default_factory=list)  # (ref contig, path)


def metric_sentence(
    metric: str,
    value: float,
    stats: MetricStats | None,
    z: float | None,
    reasons: list[str],
    z_threshold: float,
) -> str:
    """``Aligned reference 52.6% vs cohort median 38.8% (much higher than the rest; z=+4.1)``."""
    subject, comparison = SENTENCE_TEMPLATES.get(metric, (f"{metric} {{v}}", "vs median {m}"))
    text = subject.format(v=_fmt_metric_value(metric, value))
    if stats is not None:
        text += " " + comparison.format(m=_fmt_metric_value(metric, stats.median))
    notes: list[str] = []
    if z is not None and abs(z) > z_threshold:
        notes.append("much higher than the rest" if z > 0 else "much lower than the rest")
        notes.append(f"z={z:+.1f}")
        if stats is not None and stats.mad == 0:
            notes.append("other samples nearly identical")
    for reason in reasons:
        if reason.startswith(("below_threshold_", "above_threshold_")):
            word, limit = reason.split("_threshold_")
            option = THRESHOLD_OPTION_FOR_METRIC.get(metric, "threshold")
            notes.append(f"{word} the {option} threshold of {limit}")
    if notes:
        text += f" ({'; '.join(notes)})"
    return text


def _num(row: dict[str, str], col: str) -> float:
    return parse_number(row.get(col)) or 0.0


def _dotplot_for(ref_rows: list[dict[str, str]], contigs: list[str]) -> list[tuple[str, str]]:
    paths = {r.get("reference_contig", ""): (r.get("dotplot") or "").strip() for r in ref_rows}
    out = []
    for c in contigs:
        p = paths.get(c, "")
        if p and (c, p) not in out:
            out.append((c, p))
    return out


def _strand_word(strand: str) -> str:
    return {"-": "inverted", "+": "forward"}.get(strand.strip(), strand.strip() or "?")


def _pointer_overlap(
    query_rows: list[dict[str, str]],
    nested_rows: list[dict[str, str]],
    ref_rows: list[dict[str, str]],
    bp_rows: list[dict[str, str]] | None = None,
) -> tuple[str, list[tuple[str, str]]]:
    parts: list[str] = []
    contigs: list[str] = []
    candidates = [r for r in query_rows if _num(r, "overlapping_query_bp") > 0]
    top = max(candidates, key=lambda r: _num(r, "overlapping_query_bp"), default=None)
    pool = nested_rows
    if top is not None:
        name = top.get("query_contig", "")
        overlap = _num(top, "overlapping_query_bp")
        length = parse_number(top.get("query_length_bp"))
        of = f" of {fmt_bp(length)}" if length else ""
        parts.append(f"most on query contig {name} ({fmt_bp(overlap)}{of})")
        on_top = [r for r in nested_rows if r.get("query_contig") == name]
        pool = on_top or nested_rows
        junctions = sorted(
            v for v in (parse_number(r.get("query_junction_start"))
                        for r in (bp_rows or []) if r.get("query_contig") == name)
            if v is not None
        )
        if junctions:
            shown = ", ".join(_fmt_pos(v) for v in junctions[:4])
            more = f" and {len(junctions) - 4} more" if len(junctions) > 4 else ""
            parts.append(f"breakpoints on {name} at query {shown}{more}")
            for r in bp_rows or []:
                if r.get("query_contig") == name:
                    for side in ("left", "right"):
                        c = (r.get(f"{side}_reference_contig") or "").strip()
                        if c and c not in NA_VALUES:
                            contigs.append(c)

    def span(r: dict[str, str]) -> float:
        a, b = parse_number(r.get("query_start")), parse_number(r.get("query_end"))
        return (b - a) if a is not None and b is not None else -1.0

    nb = max(pool, key=span, default=None)
    if nb is not None and span(nb) >= 0:
        qs, qe = parse_number(nb.get("query_start")), parse_number(nb.get("query_end"))
        text = (f"largest nested block ({fmt_bp(qe - qs)}): query "
                f"{fmt_span(nb.get('query_contig', ''), qs, qe)}")
        rs, re_ = parse_number(nb.get("reference_start")), parse_number(nb.get("reference_end"))
        rc = (nb.get("reference_contig") or "").strip()
        if rc not in NA_VALUES and rs is not None and re_ is not None:
            text += (f" also aligns {_strand_word(nb.get('strand') or '')} to reference "
                     f"{fmt_span(rc, rs, re_)}")
            contigs.append(rc)
        cc = (nb.get("container_reference_contig") or "").strip()
        cs = parse_number(nb.get("container_reference_start"))
        ce = parse_number(nb.get("container_reference_end"))
        if cc not in NA_VALUES and cs is not None and ce is not None:
            text += (f", inside a {_strand_word(nb.get('container_strand') or '')} block "
                     f"on {fmt_span(cc, cs, ce)}")
            contigs.append(cc)
        parts.append(text)
    return "; ".join(parts), _dotplot_for(ref_rows, contigs)


def _pointer_breakpoints(
    bp_rows: list[dict[str, str]],
    query_rows: list[dict[str, str]],
    ref_rows: list[dict[str, str]],
    *,
    private_only: bool,
) -> tuple[str, list[tuple[str, str]]]:
    rows = [r for r in bp_rows if not private_only or r.get("recurrent") != "true"]
    counts: dict[str, int] = {}
    refs: dict[str, list[str]] = {}
    for r in rows:
        q = r.get("query_contig", "")
        counts[q] = counts.get(q, 0) + 1
        for side in ("left", "right"):
            c = (r.get(f"{side}_reference_contig") or "").strip()
            if c and c not in NA_VALUES and c not in refs.setdefault(q, []):
                refs[q].append(c)
    if not counts and not private_only:
        for r in query_rows:
            n = int(_num(r, "breakpoint_adjacencies"))
            if n > 0:
                counts[r.get("query_contig", "")] = n
    if not counts:
        return "", []
    top = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:3]
    text = "on query contig" + ("s " if len(top) > 1 else " ")
    text += ", ".join(f"{q} ({n})" for q, n in top)
    if len(counts) > 3:
        text += f" and {len(counts) - 3} more"
    involved: list[str] = []
    for q, _ in top:
        for c in refs.get(q, []):
            if c not in involved:
                involved.append(c)
    if involved:
        text += "; reference contig" + ("s " if len(involved) > 1 else " ") + ", ".join(involved)
    return text, _dotplot_for(ref_rows, involved[:3])


def contig_medians(
    by_reference: dict[str, list[dict[str, str]]], column: str, min_samples: int = 3
) -> dict[str, float]:
    """Per reference contig: cohort median of ``column`` (contigs seen in >= min_samples)."""
    values: dict[str, list[float]] = {}
    for rows in by_reference.values():
        for r in rows:
            v = parse_number(r.get(column))
            if v is not None:
                values.setdefault(r.get("reference_contig", ""), []).append(v)
    return {c: statistics.median(v) for c, v in values.items() if len(v) >= min_samples}


def _pointer_contigs(
    column: str,
    sign: float,
    ref_rows: list[dict[str, str]],
    medians: dict[str, float],
) -> tuple[str, list[tuple[str, str]]]:
    """Reference contigs that deviate most (in direction ``sign``) from their cohort median.

    Contigs shorter than 1% of the total reference length (unplaced scaffolds)
    are ignored.
    """
    total = sum(_num(r, "reference_length_bp") for r in ref_rows)
    devs: list[tuple[float, str, float, float]] = []
    for r in ref_rows:
        c = r.get("reference_contig", "")
        length = parse_number(r.get("reference_length_bp"))
        if total > 0 and length is not None and length < 0.01 * total:
            continue
        v = parse_number(r.get(column))
        if v is None or c not in medians:
            continue
        devs.append((sign * (v - medians[c]), c, v, medians[c]))
    if not devs:
        return "", []
    devs.sort(key=lambda d: (-d[0], d[1]))
    floor = 0.5 if column == "identity_pct" else 1.0
    if devs[0][0] <= floor:
        return "", []
    decimals = 2 if column == "identity_pct" else 1
    word = "higher" if sign > 0 else "lower"

    def show(d: tuple[float, str, float, float]) -> str:
        return f"{d[1]} {d[2]:.{decimals}f}% vs {d[3]:.{decimals}f}%"

    same = sum(1 for d in devs if d[0] > floor)
    typical = statistics.median(d[0] for d in devs)
    if same >= 0.75 * len(devs) and len(devs) >= 2:
        text = (f"{word} on {same} of {len(devs)} reference contigs (genome-wide); "
                f"most: {show(devs[0])} (contig median)")
        top = [devs[0]]
    elif devs[0][0] >= 2 * max(typical, 0.0):
        top = [d for d in devs[:3] if d[0] > floor and d[0] >= devs[0][0] / 2]
        text = "mostly " + ", ".join(show(d) for d in top) + " (contig medians)"
    else:
        return "", []
    return text, _dotplot_for(ref_rows, [d[1] for d in top])


def explain_sample(
    row: dict[str, str],
    sample_flags: dict[str, list[str]],
    sample_unusual: dict[str, list[str]],
    z_scores: dict[str, float],
    stats: dict[str, MetricStats | None],
    *,
    ref_rows: list[dict[str, str]],
    query_rows: list[dict[str, str]],
    bp_rows: list[dict[str, str]],
    nested_rows: list[dict[str, str]],
    medians: dict[str, dict[str, float]],
    z_threshold: float,
) -> list[Reason]:
    """Plain-language reasons (review first, then unusual) in FLAG_DIRECTIONS order.

    ``medians`` maps a by-reference column (``aligned_reference_pct``,
    ``identity_pct``) to per-contig cohort medians (see ``contig_medians``).
    """
    out: list[Reason] = []
    for kind, source in ((STATUS_REVIEW, sample_flags), (STATUS_UNUSUAL, sample_unusual)):
        for metric in FLAG_DIRECTIONS:
            if metric not in source:
                continue
            value = parse_number(row.get(metric))
            if value is None:
                continue
            sentence = metric_sentence(
                metric, value, stats.get(metric), z_scores.get(metric), source[metric],
                z_threshold,
            )
            where, plots = "", []
            if metric == "overlapping_query_bp":
                where, plots = _pointer_overlap(query_rows, nested_rows, ref_rows, bp_rows)
            elif metric in ("breakpoints_per_gb_aligned", "breakpoints_private"):
                where, plots = _pointer_breakpoints(
                    bp_rows, query_rows, ref_rows,
                    private_only=metric == "breakpoints_private",
                )
            elif metric in ("aligned_reference_pct", "aligned_query_pct", "identity_pct"):
                column = "identity_pct" if metric == "identity_pct" else "aligned_reference_pct"
                st = stats.get(metric)
                sign = 1.0 if (st is not None and value > st.median) else -1.0
                where, plots = _pointer_contigs(column, sign, ref_rows, medians.get(column, {}))
            out.append(Reason(metric, kind, sentence, where, plots))
    return out


# ── shared features (cross-sample clusters) ──────────────────────────────────

MAX_HTML_SHARED_FEATURES = 30


@dataclass
class SharedFeature:
    kind: str  # "breakpoint" or "nested block"
    region: str  # plain text
    samples: list[str]  # distinct samples, input order
    rows: int
    detail: str = ""  # e.g. breakpoint types
    sort_key: tuple = ()


class _UnionFind:
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))

    def find(self, i: int) -> int:
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def _contig_sort_key(contig: str) -> tuple:
    digits = "".join(ch for ch in contig if ch.isdigit())
    return (int(digits) if digits else float("inf"), contig)


def _range_text(contig: str, lo: float, hi: float) -> str:
    if hi - lo < 1:
        return f"{contig} {_fmt_pos(lo)}"
    return fmt_span(contig, lo, hi).replace(":", " ", 1)


def cluster_shared_features(
    samples: list[str],
    breakpoints: list[dict[str, str]],
    nested: list[dict[str, str]],
    window_bp: int = DEFAULT_RECURRENCE_WINDOW_BP,
) -> list[SharedFeature]:
    """Group breakpoints / nested blocks that match across samples into clusters.

    Matching uses the same rules as the recurrence annotation (breakpoints:
    both ends on the same reference contigs within ``window_bp``, either
    orientation; nested blocks: same reference contig, start and end each
    within ``window_bp``), linked transitively (single linkage). Only clusters
    with at least two distinct samples are returned, largest first.
    """
    if len(samples) < 2:
        return []
    order = {s: i for i, s in enumerate(samples)}
    bin_width = max(int(window_bp), 1)
    features: list[SharedFeature] = []

    # breakpoints
    ends = [_breakpoint_ends(r) for r in breakpoints]
    uf = _UnionFind(len(breakpoints))
    index: dict[tuple[str, str, str, int], list[int]] = {}
    for i, e in enumerate(ends):
        if e is None:
            continue
        key = tuple(sorted((e[0][0], e[1][0])))
        for contig, pos in e:
            index.setdefault((*key, contig, int(pos // bin_width)), []).append(i)  # type: ignore[arg-type]
    for i, e in enumerate(ends):
        if e is None:
            continue
        key = tuple(sorted((e[0][0], e[1][0])))
        contig, pos = e[0]
        b = int(pos // bin_width)
        for nb in (b - 1, b, b + 1):
            for j in index.get((*key, contig, nb), ()):  # type: ignore[arg-type]
                if j != i and _ends_match(e, ends[j], window_bp):  # type: ignore[arg-type]
                    uf.union(i, j)
    groups: dict[int, list[int]] = {}
    for i, e in enumerate(ends):
        if e is not None:
            groups.setdefault(uf.find(i), []).append(i)
    for members in groups.values():
        names = sorted({breakpoints[i]["sample"] for i in members}, key=lambda s: order.get(s, 0))
        if len(names) < 2:
            continue
        normed = [sorted(ends[i], key=lambda x: (_contig_sort_key(x[0]), x[1])) for i in members]  # type: ignore[arg-type]
        a_contig, b_contig = normed[0][0][0], normed[0][1][0]
        a_pos = [n[0][1] for n in normed if n[0][0] == a_contig]
        b_pos = [n[1][1] for n in normed if n[1][0] == b_contig]
        region = (f"{_range_text(a_contig, min(a_pos), max(a_pos))} ↔ "
                  f"{_range_text(b_contig, min(b_pos), max(b_pos))}")
        types: list[str] = []
        for i in members:
            for t in _breakpoint_types(breakpoints[i]).split(", "):
                if t and t not in types:
                    types.append(t)
        features.append(SharedFeature(
            "breakpoint", region, names, len(members), ", ".join(types),
            (-len(names), 0, _contig_sort_key(a_contig), min(a_pos)),
        ))

    # nested blocks
    keys = [_nested_key(r) for r in nested]
    uf = _UnionFind(len(nested))
    nindex: dict[tuple[str, int], list[int]] = {}
    for i, k in enumerate(keys):
        if k is not None:
            nindex.setdefault((k[0], int(k[1] // bin_width)), []).append(i)
    for i, k in enumerate(keys):
        if k is None:
            continue
        b = int(k[1] // bin_width)
        for nb in (b - 1, b, b + 1):
            for j in nindex.get((k[0], nb), ()):
                o = keys[j]
                if j != i and abs(o[1] - k[1]) <= window_bp and abs(o[2] - k[2]) <= window_bp:  # type: ignore[index]
                    uf.union(i, j)
    ngroups: dict[int, list[int]] = {}
    for i, k in enumerate(keys):
        if k is not None:
            ngroups.setdefault(uf.find(i), []).append(i)
    for members in ngroups.values():
        names = sorted({nested[i]["sample"] for i in members}, key=lambda s: order.get(s, 0))
        if len(names) < 2:
            continue
        contig = keys[members[0]][0]  # type: ignore[index]
        lo = min(keys[i][1] for i in members)  # type: ignore[index]
        hi = max(keys[i][2] for i in members)  # type: ignore[index]
        strands = sorted({(nested[i].get("strand") or "").strip() for i in members} - {""})
        features.append(SharedFeature(
            "nested block", fmt_span(contig, lo, hi), names, len(members),
            "strand " + "/".join(strands) if strands else "",
            (-len(names), 1, _contig_sort_key(contig), lo),
        ))

    features.sort(key=lambda f: f.sort_key)
    return features


# ── HTML helpers ─────────────────────────────────────────────────────────────


NA_HTML = '<span class="na">NA</span>'
COLOR_REVIEW = "#D6423F"
COLOR_UNUSUAL = "#7B5EA7"
COLOR_NORMAL = "#7A93AC"


def esc(text: object) -> str:
    return html.escape(str(text), quote=True)


def fmt_value(raw: str | None, column: str) -> str:
    """Format a cell for HTML (already escaped). NA -> muted 'NA'."""
    if column in NAME_COLUMNS:
        return esc(raw or "")
    if column in CATEGORY_COLUMNS:
        if raw is None or raw.strip() in NA_VALUES:
            return NA_HTML
        return esc(raw)
    value = parse_number(raw)
    if value is None:
        if raw is not None and raw.strip() not in NA_VALUES:
            return esc(raw)
        return NA_HTML
    if column in FLOAT_COLUMNS:
        return f"{value:,.2f}"
    if value.is_integer():
        return f"{int(value):,}"
    return f"{value:,.2f}"


def _fmt_axis(value: float, metric: str) -> str:
    if metric in FLOAT_COLUMNS:
        return f"{value:,.2f}"
    if abs(value - round(value)) < 1e-9:
        return f"{int(round(value)):,}"
    return f"{value:,.1f}"


def _is_na(raw: str | None) -> bool:
    return raw is None or raw.strip() in NA_VALUES


def svg_strip_plot(
    metric: str,
    samples: list[str],
    values: list[float | None],
    flagged: list[bool],
    unusual: list[bool] | None = None,
    *,
    width: int = 900,
    height: int = 74,
) -> str:
    """One-dimensional strip plot: one dot per sample, median line.

    Flagged (review) dots are red, unusual dots purple, the rest grey-blue.
    """
    unusual = unusual or [False] * len(samples)
    margin = {"left": 200, "right": 30, "top": 12, "bottom": 24}
    plot_w = width - margin["left"] - margin["right"]
    plot_h = height - margin["top"] - margin["bottom"]
    label = METRIC_LABELS.get(metric, metric)
    parts = [
        f'<svg class="strip" width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img" aria-label="{esc(label)}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{margin["left"] - 10}" y="{margin["top"] + plot_h / 2 + 4:.2f}" '
        f'text-anchor="end" font-size="12" font-family="sans-serif">{esc(label)}</text>',
    ]
    present = [v for v in values if v is not None]
    axis_y = margin["top"] + plot_h
    parts.append(
        f'<line x1="{margin["left"]}" y1="{axis_y}" x2="{margin["left"] + plot_w}" '
        f'y2="{axis_y}" stroke="#333" stroke-width="1"/>'
    )
    if not present:
        parts.append(
            f'<text x="{margin["left"] + plot_w / 2}" y="{margin["top"] + plot_h / 2 + 4:.2f}" '
            f'text-anchor="middle" font-size="11" font-family="sans-serif" fill="#777">'
            f"no values</text>"
        )
        parts.append("</svg>")
        return "\n".join(parts)

    lo, hi = min(present), max(present)
    if hi == lo:
        pad = abs(lo) * 0.05 or 1.0
        lo, hi = lo - pad, hi + pad

    def x_scale(v: float) -> float:
        return margin["left"] + (v - lo) / (hi - lo) * plot_w

    for frac in (0.0, 0.5, 1.0):
        val = lo + frac * (hi - lo)
        x = margin["left"] + frac * plot_w
        anchor = "start" if frac == 0 else ("end" if frac == 1 else "middle")
        parts.append(
            f'<line x1="{x:.2f}" y1="{axis_y}" x2="{x:.2f}" y2="{axis_y + 4}" '
            f'stroke="#333" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{x:.2f}" y="{axis_y + 16}" text-anchor="{anchor}" font-size="10" '
            f'font-family="sans-serif">{esc(_fmt_axis(val, metric))}</text>'
        )

    med = statistics.median(present)
    mx = x_scale(med)
    parts.append(
        f'<line x1="{mx:.2f}" y1="{margin["top"]}" x2="{mx:.2f}" y2="{axis_y}" '
        f'stroke="#888" stroke-width="1.5" stroke-dasharray="4 3">'
        f"<title>median {esc(_fmt_axis(med, metric))}</title></line>"
    )

    lanes = 5
    lane_h = plot_h / (lanes + 1)
    # Draw normal dots first so flagged/unusual ones sit on top.
    order = sorted(range(len(samples)), key=lambda i: (flagged[i], unusual[i]))
    for i in order:
        v = values[i]
        if v is None:
            continue
        cy = margin["top"] + lane_h * (1 + i % lanes)
        if flagged[i]:
            color, cls = COLOR_REVIEW, "dot flagged"
        elif unusual[i]:
            color, cls = COLOR_UNUSUAL, "dot unusual"
        else:
            color, cls = COLOR_NORMAL, "dot"
        parts.append(
            f'<circle class="{cls}" cx="{x_scale(v):.2f}" cy="{cy:.2f}" r="4" fill="{color}" '
            f'fill-opacity="0.85" stroke="white" stroke-width="0.5">'
            f"<title>{esc(samples[i])}: {esc(_fmt_axis(v, metric))}</title></circle>"
        )
    parts.append("</svg>")
    return "\n".join(parts)


def _breakpoint_first(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Rows with breakpoint_adjacencies > 0 first (most first), others in input order."""
    def key(item: tuple[int, dict[str, str]]) -> tuple[int, float, int]:
        idx, row = item
        bp = parse_number(row.get("breakpoint_adjacencies")) or 0.0
        return (0 if bp > 0 else 1, -bp if bp > 0 else 0.0, idx)

    return [row for _, row in sorted(enumerate(rows), key=key)]


# Per-sample contig/breakpoint tables in the HTML stop here (breakpoint contigs
# come first); fragmented assemblies can have thousands of scaffolds. The TSVs
# keep every row.
MAX_HTML_CONTIG_ROWS = 50
# Cross-sample breakpoint table cap (recurrent breakpoints first).
MAX_HTML_BREAKPOINT_ROWS = 200
# Cross-sample nested-block table cap (recurrent first).
MAX_HTML_NESTED_ROWS = 200


def _cap_note(shown: int, total: int, noun: str, where: str) -> str:
    if total <= shown:
        return ""
    return (
        f'<p class="na">Showing {shown:,} of {total:,} {noun}; '
        f"{where} has every row.</p>\n"
    )


def _contig_table(rows: list[dict[str, str]], columns: list[str], css_class: str) -> str:
    shown = rows[:MAX_HTML_CONTIG_ROWS]
    out = [
        _cap_note(len(shown), len(rows), "contigs", "the per-sample TSV"),
        f'<table class="{css_class}">\n<tr>',
    ]
    out.extend(f"<th>{esc(METRIC_LABELS.get(c, c))}</th>" for c in columns)
    out.append("</tr>\n")
    for row in shown:
        bp = parse_number(row.get("breakpoint_adjacencies")) or 0.0
        cls = ' class="bp"' if bp > 0 else ""
        out.append(f"<tr{cls}>")
        for c in columns:
            out.append(f"<td>{fmt_value(row.get(c), c)}</td>")
        out.append("</tr>\n")
    out.append("</table>\n")
    return "".join(out)


def _breakpoint_types(row: dict[str, str]) -> str:
    types = [
        name
        for col, name in (
            ("strand_flip", "flip"),
            ("reference_contig_jump", "jump"),
            ("out_of_order", "out_of_order"),
        )
        if (row.get(col) or "").strip().lower() == "true"
    ]
    return ", ".join(types)


def _ref_position(row: dict[str, str], side: str) -> str:
    contig = row.get(f"{side}_reference_contig") or ""
    pos = parse_number(row.get(f"{side}_reference_pos"))
    if contig.strip() in NA_VALUES or pos is None:
        return NA_HTML
    shown = f"{int(pos):,}" if pos.is_integer() else f"{pos:,.1f}"
    strand = (row.get(f"{side}_strand") or "").strip()
    suffix = f" ({esc(strand)})" if strand and strand not in NA_VALUES else ""
    return f"{esc(contig)}:{shown}{suffix}"


def _recurrence_sort(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Recurrent first (most shared first); otherwise input order."""
    def key(item: tuple[int, dict[str, str]]) -> tuple[float, int]:
        idx, row = item
        n = parse_number(row.get("recurrence_samples")) or 0.0
        return (-n, idx)

    return [row for _, row in sorted(enumerate(rows), key=key)]


def _breakpoint_table(
    rows: list[dict[str, str]], *, cap: int, with_sample: bool, samples: list[str], where: str
) -> str:
    shown = rows[:cap]
    header = (["Sample"] if with_sample else []) + [
        "Query contig", "Junction (query)", "Location", "Types",
        "Left ref position", "Right ref position", METRIC_LABELS["recurrence_samples"],
    ]
    out = [
        _cap_note(len(shown), len(rows), "breakpoints", where),
        '<table class="bp-table">\n<tr>',
    ]
    out.extend(f"<th>{esc(h)}</th>" for h in header)
    out.append("</tr>\n")
    for row in shown:
        recurrent = row.get("recurrent") == "true"
        out.append('<tr class="recurrent">' if recurrent else "<tr>")
        if with_sample:
            sample = row["sample"]
            anchor = _sample_anchor(samples, sample) if sample in samples else ""
            out.append(f'<td><a href="#{esc(anchor)}">{esc(sample)}</a></td>')
        start = fmt_value(row.get("query_junction_start"), "query_junction_start")
        end = fmt_value(row.get("query_junction_end"), "query_junction_end")
        loc = fmt_value(row.get("location"), "location")
        loc_raw = (row.get("location") or "").strip()
        loc_cls = f' class="loc-{esc(loc_raw)}"' if loc_raw in ("contig_end", "n_gap", "interior") else ""
        out.append(
            f"<td>{fmt_value(row.get('query_contig'), 'query_contig')}</td>"
            f"<td>{start}&ndash;{end}</td>"
            f"<td{loc_cls}>{loc}</td>"
            f"<td>{esc(_breakpoint_types(row))}</td>"
            f"<td>{_ref_position(row, 'left')}</td>"
            f"<td>{_ref_position(row, 'right')}</td>"
            f"<td>{fmt_value(row.get('recurrence_samples'), 'recurrence_samples')}</td>"
        )
        out.append("</tr>\n")
    out.append("</table>\n")
    return "".join(out)


def _telomere_cell(row: dict[str, str]) -> str:
    """``major_telomeric_ends / 2 x major_sequences`` (e.g. ``14 / 20``); NA-safe."""
    ends = parse_number(row.get("assembly_major_telomeric_ends"))
    major = parse_number(row.get("assembly_major_sequences"))
    if ends is None and major is None:
        return NA_HTML
    left = f"{int(ends):,}" if ends is not None else NA_HTML
    right = f"{int(2 * major):,}" if major is not None else NA_HTML
    return f"{left} / {right}"


def _interval(start: str | None, end: str | None) -> str:
    a = fmt_value(start, "query_start")
    b = fmt_value(end, "query_end")
    return f"{a}&ndash;{b}"


def _ref_interval(row: dict[str, str], prefix: str, *, strand_col: str | None = None) -> str:
    contig = row.get(f"{prefix}reference_contig") or ""
    if contig.strip() in NA_VALUES:
        return NA_HTML
    text = f"{esc(contig)}:{_interval(row.get(f'{prefix}reference_start'), row.get(f'{prefix}reference_end'))}"
    if strand_col:
        strand = (row.get(strand_col) or "").strip()
        if strand and strand not in NA_VALUES:
            text += f" ({esc(strand)})"
    return text


def _nested_table(
    rows: list[dict[str, str]], *, cap: int, with_sample: bool, samples: list[str], where: str
) -> str:
    shown = rows[:cap]
    header = (["Sample"] if with_sample else []) + [
        "Query contig", "Query interval", "Strand", "Reference",
        "Container reference (strand)", METRIC_LABELS["recurrence_samples"],
    ]
    out = [
        _cap_note(len(shown), len(rows), "nested blocks", where),
        '<table class="nested-table">\n<tr>',
    ]
    out.extend(f"<th>{esc(h)}</th>" for h in header)
    out.append("</tr>\n")
    for row in shown:
        recurrent = row.get("recurrent") == "true"
        out.append('<tr class="recurrent">' if recurrent else "<tr>")
        if with_sample:
            sample = row["sample"]
            anchor = _sample_anchor(samples, sample) if sample in samples else ""
            out.append(f'<td><a href="#{esc(anchor)}">{esc(sample)}</a></td>')
        out.append(
            f"<td>{fmt_value(row.get('query_contig'), 'query_contig')}</td>"
            f"<td>{_interval(row.get('query_start'), row.get('query_end'))}</td>"
            f"<td>{fmt_value(row.get('strand'), 'strand')}</td>"
            f"<td>{_ref_interval(row, '')}</td>"
            f"<td>{_ref_interval(row, 'container_', strand_col='container_strand')}</td>"
            f"<td>{fmt_value(row.get('recurrence_samples'), 'recurrence_samples')}</td>"
        )
        out.append("</tr>\n")
    out.append("</table>\n")
    return "".join(out)


def _count_mismatches(
    summaries: list[dict[str, str]], rows: list[dict[str, str]], column: str
) -> list[tuple[str, int, int]]:
    """``(sample, summary count, row count)`` where ``column`` disagrees with the rows."""
    per_sample: dict[str, int] = {}
    for r in rows:
        per_sample[r["sample"]] = per_sample.get(r["sample"], 0) + 1
    mismatched = []
    for s in summaries:
        expected = parse_number(s.get(column))
        got = per_sample.get(s["sample"], 0)
        if expected is not None and int(expected) != got:
            mismatched.append((s["sample"], int(expected), got))
    return mismatched


def _mismatch_note(mismatched: list[tuple[str, int, int]], what: str, files: str) -> str:
    if not mismatched:
        return ""
    items = ", ".join(
        f"{esc(name)} ({exp:,} in summary, {got:,} rows)" for name, exp, got in mismatched
    )
    return (
        f'<p class="note"><strong>{what} table does not match summary</strong> for: '
        f"{items}. Missing {files} count as zero.</p>\n"
    )


def _sample_table(
    summaries: list[dict[str, str]],
    flags: list[dict[str, list[str]]],
    samples: list[str],
    columns: list[str],
    css_class: str,
    *,
    flags_column: bool,
    unusual: list[dict[str, list[str]]] | None = None,
) -> str:
    """Full table; flagged cells (review) and unusual cells highlighted, tooltip = raw reason."""
    unusual = unusual or [{} for _ in summaries]
    out = [f'<div class="scroll"><table class="{css_class}">\n<tr><th>Sample</th>']
    out.extend(f"<th>{esc(METRIC_LABELS.get(c, c))}</th>" for c in columns)
    if flags_column:
        out.append("<th>Unusual</th><th>Flags</th>")
    out.append("</tr>\n")
    for row, sample_flags, sample_unusual in zip(summaries, flags, unusual):
        sample = row["sample"]
        out.append(f'<tr><td><a href="#{esc(_sample_anchor(samples, sample))}">{esc(sample)}</a></td>')
        for col in columns:
            if col == TELOMERE_DISPLAY_COLUMN:
                out.append(f"<td>{_telomere_cell(row)}</td>")
                continue
            cell = fmt_value(row.get(col), col)
            if col == "aligned_query_pct" and row.get("query_length_source") == "maf_srcsize":
                cell += " &dagger;"
            reasons = sample_flags.get(col)
            odd = sample_unusual.get(col)
            if reasons:
                tip = "; ".join(f"{col}:{r}" for r in reasons)
                out.append(f'<td class="flag" title="{esc(tip)}">{cell}</td>')
            elif odd:
                tip = "; ".join(f"{col}:{r}" for r in odd)
                out.append(f'<td class="unusual" title="{esc(tip)}">{cell}</td>')
            else:
                out.append(f"<td>{cell}</td>")
        if flags_column:
            out.append(
                f'<td class="flags unusual-list">'
                f'{esc(format_flags(sample_unusual).replace(";", "; "))}</td>'
                f'<td class="flags">{esc(format_flags(sample_flags).replace(";", "; "))}</td>'
            )
        out.append("</tr>\n")
    out.append("</table></div>\n")
    return "".join(out)


CSS = """\
:root{--review:#D6423F;--review-bg:#fde2e1;--unusual:#7B5EA7;--unusual-bg:#ece6f6;
--ok:#2e7d4f;--muted:#6b7280;--line:#d0d7de}
body{font-family:system-ui,-apple-system,"Segoe UI",sans-serif;margin:24px auto;padding:0 16px;
color:#111;max-width:1200px;line-height:1.45;background:#fff}
table{border-collapse:collapse;margin:12px 0}
th,td{border:1px solid #ccc;padding:4px 10px;text-align:right}
th{background:#f0f0f0;text-align:center}
td:first-child{text-align:left}
h1{margin:0.2em 0 0.1em;font-size:1.6em}
h2{margin-top:1.6em;font-size:1.25em;border-bottom:1px solid var(--line);padding-bottom:4px}
h3{margin-top:1.2em}
p.meta{color:var(--muted);margin-top:0}
code{background:#f6f8fa;padding:0 3px;border-radius:3px}
details{margin:8px 0;border:1px solid var(--line);border-radius:6px;padding:4px 12px}
details.drill>summary{font-size:1.1em}
summary{cursor:pointer;font-weight:bold;padding:6px 0;list-style:revert}
summary:hover{color:#0969da}
.summary-box{border:1px solid var(--line);border-radius:8px;padding:12px 16px;background:#f8fafc}
.counts{font-size:1.25em;font-weight:600;margin:0 0 6px}
.counts .ok{color:var(--ok)} .counts .review{color:var(--review)} .counts .unusual{color:var(--unusual)}
.cohort{margin:4px 0}
.card{border:1px solid var(--line);border-left:6px solid var(--review);border-radius:6px;
padding:8px 14px;margin:10px 0;background:#fff}
.card.unusual{border-left-color:var(--unusual)}
.card h3{margin:4px 0 6px;font-size:1.05em}
.card ul{margin:4px 0 6px 0;padding-left:20px}
.card li{margin:3px 0}
.where{color:#374151;font-size:0.92em}
.card .thumbs img{border:1px solid #ccc}
.card-link{font-size:0.9em}
td.status{text-align:center;cursor:help}
td.status.review,.icon.review{color:var(--review)}
td.status.unusual,.icon.unusual{color:var(--unusual)}
td.status.ok,.icon.ok{color:var(--ok)}
td.flag{background:var(--review-bg);font-weight:bold;cursor:help}
td.unusual{background:var(--unusual-bg);cursor:help}
td.flags{text-align:left;font-family:monospace;font-size:0.85em;color:#a40000}
td.flags.unusual-list{color:var(--unusual)}
td.left{text-align:left}
tr.bp td{background:#fff8dc}
tr.recurrent td{background:#e8f4ea}
td.loc-contig_end{color:#8a5a00}
td.loc-n_gap{color:#6a3d9a}
td.loc-interior{color:#a40000}
.na{color:#999}
.note{background:#fff8dc;border:1px solid #e0c96b;border-radius:6px;padding:8px 12px}
.params td:first-child{font-family:monospace}
section.ref-contigs{border-left:4px solid #4C78A8;padding-left:12px;margin:12px 0}
section.query-contigs{border-left:4px solid #F58518;padding-left:12px;margin:12px 0}
section.breakpoints{border-left:4px solid #B279A2;padding-left:12px;margin:12px 0}
section.nested{border-left:4px solid #9D755D;padding-left:12px;margin:12px 0}
section.reasons{border-left:4px solid var(--review);padding-left:12px;margin:12px 0}
table.maf-contiguity{font-size:0.85em;color:#555}
section.dotplots{border-left:4px solid #54A24B;padding-left:12px;margin:12px 0}
.scroll{overflow-x:auto}
.thumbs{display:flex;flex-wrap:wrap;gap:12px}
.thumbs figure{margin:0;text-align:center;font-size:0.85em}
.thumbs img{border:1px solid #ccc;max-width:100%;height:auto}
dl.defs dt{font-family:monospace;font-weight:bold;margin-top:8px}
dl.defs dd{margin-left:20px}
.legend span{display:inline-block;width:10px;height:10px;border-radius:5px;margin:0 4px 0 12px}
svg.strip{max-width:100%;height:auto}
"""

# Opens every collapsed <details> that contains the target of an in-page link,
# so links from the summary into the drill-down sections land on visible content.
OPEN_TARGET_JS = """\
function argprepOpenTarget(){var id=decodeURIComponent(location.hash.slice(1));
if(!id)return;var t=document.getElementById(id);if(!t)return;
for(var el=t;el;el=el.parentElement){if(el.tagName==='DETAILS')el.setAttribute('open','');}
t.scrollIntoView();}
window.addEventListener('hashchange',argprepOpenTarget);
window.addEventListener('DOMContentLoaded',argprepOpenTarget);
"""

# Compact sample table: (column key, label). Column keys are summary columns or
# TELOMERE_DISPLAY_COLUMN.
COMPACT_COLUMNS_ASSEMBLY: list[tuple[str, str]] = [
    ("aligned_reference_pct", "Aligned ref %"),
    ("identity_pct", "Identity %"),
    ("breakpoint_adjacencies", "Breakpoints"),
    ("nested_blocks", "Nested blocks"),
    ("assembly_contig_n50_bp", "Contig N50"),
    (TELOMERE_DISPLAY_COLUMN, "Telomeres"),
]
COMPACT_COLUMNS_MAF_ONLY: list[tuple[str, str]] = [
    ("aligned_reference_pct", "Aligned ref %"),
    ("aligned_query_pct", "Aligned query %"),
    ("identity_pct", "Identity %"),
    ("breakpoint_adjacencies", "Breakpoints"),
    ("breakpoints_private", "Private breakpoints"),
    ("nested_blocks", "Nested blocks"),
]
# Compact column -> flagged metrics whose reasons highlight that cell.
COMPACT_CELL_METRICS: dict[str, tuple[str, ...]] = {
    "aligned_reference_pct": ("aligned_reference_pct",),
    "aligned_query_pct": ("aligned_query_pct",),
    "identity_pct": ("identity_pct",),
    "breakpoint_adjacencies": ("breakpoints_per_gb_aligned",),
    "breakpoints_private": ("breakpoints_private",),
    "nested_blocks": ("overlapping_query_bp",),
    "assembly_contig_n50_bp": ("assembly_contig_n50_bp",),
}

STATUS_DEFINITIONS: list[tuple[str, str]] = [
    ("status",
     "Per-sample verdict: <code>review</code> if any metric is flagged, else "
     "<code>unusual</code> if any metric is unusual, else <code>ok</code> (shown as normal)."),
    ("review",
     "A metric is flagged when its robust z-score is beyond the threshold in the "
     "<em>bad</em> direction (e.g. far less of the reference aligned than the rest of the "
     "cohort), or when it fails an absolute threshold given on the command line. Listed in "
     "the <code>flags</code> column."),
    ("unusual",
     "A metric is unusual when its robust z-score is beyond the threshold in the "
     "<em>good</em> direction (e.g. far more similar to the reference than the rest). Not "
     "a failure, but worth a look: it can mean admixture, contamination, a mislabelled "
     "sample, or a sample much closer to the reference than the others. Listed in the "
     "<code>unusual</code> column; never in <code>flags</code>."),
]


def _compact_cell(row: dict[str, str], col: str) -> str:
    if col == TELOMERE_DISPLAY_COLUMN:
        return _telomere_cell(row)
    if col == "assembly_contig_n50_bp":
        v = parse_number(row.get(col))
        return esc(fmt_bp(v)) if v is not None else NA_HTML
    cell = fmt_value(row.get(col), col)
    if col == "aligned_query_pct" and row.get("query_length_source") == "maf_srcsize":
        cell += " &dagger;"
    return cell


def _range_phrase(
    summaries: list[dict[str, str]], metric: str, fmt: str, unit: str = ""
) -> str | None:
    values = [parse_number(s.get(metric)) for s in summaries]
    present = [v for v in values if v is not None]
    if not present:
        return None
    lo, hi = min(present), max(present)
    a, b = format(lo, fmt), format(hi, fmt)
    return f"{a}{unit}" if a == b else f"{a}–{b}{unit}"


def _thumb(path: str, contig: str, width: int, caption: str = "") -> str:
    p, c = esc(path), esc(contig)
    cap = caption or c
    return (
        f'<figure><a href="{p}"><img src="{p}" loading="lazy" width="{width}" '
        f'alt="{c} dotplot"></a><figcaption>{cap}</figcaption></figure>\n'
    )


def _reason_list(reasons: list[Reason]) -> str:
    out = ["<ul>\n"]
    for r in reasons:
        icon = STATUS_ICONS[r.kind]
        out.append(f'<li><span class="icon {r.kind}">{icon}</span> {esc(r.sentence)}')
        if r.where:
            out.append(f'<br><span class="where">Where: {esc(r.where)}</span>')
        out.append("</li>\n")
    out.append("</ul>\n")
    return "".join(out)


def _card(sample: str, anchor: str, status: str, reasons: list[Reason]) -> str:
    out = [
        f'<div class="card {status}">\n<h3><span class="icon {status}">'
        f'{STATUS_ICONS[status]}</span> <a href="#{esc(anchor)}">{esc(sample)}</a></h3>\n'
    ]
    out.append(_reason_list(reasons))
    plots: list[tuple[str, str]] = []
    for r in reasons:
        for p in r.dotplots:
            if p not in plots:
                plots.append(p)
    if plots:
        out.append('<div class="thumbs">\n')
        for contig, path in plots[:4]:
            out.append(_thumb(path, contig, 200, f"{esc(contig)} dotplot"))
        out.append("</div>\n")
    out.append(
        f'<p class="card-link"><a href="#{esc(anchor)}">Full detail for {esc(sample)} '
        f"&darr;</a></p>\n</div>\n"
    )
    return "".join(out)


def render_html(
    summaries: list[dict[str, str]],
    flag_result: FlagResult,
    by_reference: dict[str, list[dict[str, str]]],
    by_query: dict[str, list[dict[str, str]]],
    breakpoints: list[dict[str, str]] | None = None,
    nested: list[dict[str, str]] | None = None,
    *,
    min_samples: int,
    z_threshold: float,
    thresholds: dict[str, float],
    recurrence_window_bp: int = DEFAULT_RECURRENCE_WINDOW_BP,
    breakpoints_supplied: bool = True,
    nested_supplied: bool = True,
    thumbnail_width: int = 280,
) -> str:
    """Render the report: plain-language summary first, drill-down sections collapsed.

    ``breakpoints`` are annotated rows from ``annotate_recurrence`` and
    ``nested`` annotated rows from ``annotate_nested_recurrence``.
    """
    breakpoints = breakpoints or []
    nested = nested or []
    samples = [s["sample"] for s in summaries]
    n = len(samples)
    flags = flag_result.flags
    unusual = flag_result.unusual or [{} for _ in summaries]
    z_scores = flag_result.z_scores or [{} for _ in summaries]
    statuses = [sample_status(f, u) for f, u in zip(flags, unusual)]

    bp_by_sample: dict[str, list[dict[str, str]]] = {}
    for r in breakpoints:
        bp_by_sample.setdefault(r["sample"], []).append(r)
    nested_by_sample: dict[str, list[dict[str, str]]] = {}
    for r in nested:
        nested_by_sample.setdefault(r["sample"], []).append(r)
    medians = {
        col: contig_medians(by_reference, col) for col in ("aligned_reference_pct", "identity_pct")
    }
    reasons: list[list[Reason]] = [
        explain_sample(
            row, flags[i], unusual[i], z_scores[i], flag_result.stats,
            ref_rows=by_reference.get(row["sample"], []),
            query_rows=by_query.get(row["sample"], []),
            bp_rows=bp_by_sample.get(row["sample"], []),
            nested_rows=nested_by_sample.get(row["sample"], []),
            medians=medians,
            z_threshold=z_threshold,
        )
        for i, row in enumerate(summaries)
    ]

    n_assembly = sum(1 for s in summaries if not _is_na(s.get("assembly_length_bp")))
    has_assembly = n_assembly > 0
    shared = cluster_shared_features(samples, breakpoints, nested, recurrence_window_bp)

    out: list[str] = []
    w = out.append

    title = f"ARGprep MAF QC — {n:,} sample{'s' if n != 1 else ''}"
    w('<!doctype html>\n<html lang="en">\n<head>\n')
    w('<meta charset="utf-8" />\n<meta name="viewport" content="width=device-width, initial-scale=1" />\n')
    w(f"<title>{esc(title)}</title>\n")
    w(f"<style>{CSS}</style>\n")
    w(f"<script>{OPEN_TARGET_JS}</script>\n")
    w("</head>\n<body>\n")

    # ── 1. title line ────────────────────────────────────────────────────────
    w(f"<h1>{esc(title)}</h1>\n")
    if n_assembly == n:
        fasta = "assembly FASTAs used"
    elif n_assembly == 0:
        fasta = "no assembly FASTAs (alignment metrics only)"
    else:
        fasta = f"assembly FASTAs used for {n_assembly:,} of {n:,} samples"
    w(
        f'<p class="meta">Generated {esc(datetime.now().strftime("%Y-%m-%d %H:%M"))} '
        f"&middot; {esc(fasta)}</p>\n"
    )

    # ── 2. summary box ───────────────────────────────────────────────────────
    n_ok = statuses.count(STATUS_OK)
    n_review = statuses.count(STATUS_REVIEW)
    n_unusual = statuses.count(STATUS_UNUSUAL)
    w('<div class="summary-box">\n')
    w(
        f'<p class="counts"><span class="ok">{STATUS_ICONS[STATUS_OK]} {n_ok:,} normal</span>'
        f' &middot; <span class="review">{STATUS_ICONS[STATUS_REVIEW]} {n_review:,} need review</span>'
        f' &middot; <span class="unusual">{STATUS_ICONS[STATUS_UNUSUAL]} {n_unusual:,} unusual</span></p>\n'
    )
    cohort_parts: list[str] = []
    phrase = _range_phrase(summaries, "aligned_reference_pct", ".1f", "%")
    if phrase:
        cohort_parts.append(f"aligned reference {phrase}")
    phrase = _range_phrase(summaries, "identity_pct", ".2f", "%")
    if phrase:
        cohort_parts.append(f"identity {phrase}")
    phrase = _range_phrase(summaries, "breakpoint_adjacencies", ",.0f")
    if phrase:
        cohort_parts.append(f"{phrase} breakpoints per sample")
    nested_values = [parse_number(s.get("nested_blocks")) for s in summaries]
    if any(v for v in nested_values if v is not None):
        phrase = _range_phrase(summaries, "nested_blocks", ",.0f")
        cohort_parts.append(f"{phrase} nested blocks per sample")
    if cohort_parts:
        w(f'<p class="cohort">Cohort: {esc(" · ".join(cohort_parts))}.</p>\n')
    if shared:
        n_bp_shared = sum(1 for f in shared if f.kind == "breakpoint")
        n_nb_shared = len(shared) - n_bp_shared
        bits = []
        if n_bp_shared:
            bits.append(f"{n_bp_shared:,} breakpoint{'s' if n_bp_shared != 1 else ''}")
        if n_nb_shared:
            bits.append(f"{n_nb_shared:,} nested block{'s' if n_nb_shared != 1 else ''}")
        w(
            f'<p class="cohort">Shared by two or more samples: {esc(" and ".join(bits))} '
            '(see <a href="#shared">Shared features</a>).</p>\n'
        )
    if n < min_samples:
        w(
            f'<p class="note"><strong>Relative (cohort) flagging skipped:</strong> only '
            f"{n:,} sample(s), fewer than min_samples = {min_samples:,}, so samples were not "
            "compared with each other. Only absolute thresholds (if any) were applied.</p>\n"
        )
    w("</div>\n")

    # ── 3. needs review / unusual cards ──────────────────────────────────────
    if n_review:
        w("<h2>Needs review</h2>\n")
        for i, row in enumerate(summaries):
            if statuses[i] == STATUS_REVIEW:
                w(_card(row["sample"], _sample_anchor(samples, row["sample"]),
                        STATUS_REVIEW, reasons[i]))
    if n_unusual:
        w("<h2>Unusual</h2>\n")
        w(
            '<p class="where">Unusual means a large deviation in the &ldquo;good&rdquo; '
            "direction, e.g. far more of the reference aligned or far more similar to the "
            "reference than the rest of the cohort. Not a failure, but check for admixture, "
            "contamination or a mislabelled sample.</p>\n"
        )
        for i, row in enumerate(summaries):
            if statuses[i] == STATUS_UNUSUAL:
                w(_card(row["sample"], _sample_anchor(samples, row["sample"]),
                        STATUS_UNUSUAL, reasons[i]))

    # ── 4. compact sample table ──────────────────────────────────────────────
    columns = COMPACT_COLUMNS_ASSEMBLY if has_assembly else COMPACT_COLUMNS_MAF_ONLY
    w("<h2>Samples</h2>\n")
    w('<div class="scroll"><table class="samples">\n<tr><th>Status</th><th>Sample</th>')
    w("".join(f"<th>{esc(label)}</th>" for _, label in columns))
    w("</tr>\n")
    for i, row in enumerate(summaries):
        sample = row["sample"]
        st = statuses[i]
        tip = " | ".join(r.sentence for r in reasons[i]) or STATUS_LABELS[st]
        w(
            f'<tr><td class="status {st}" title="{esc(tip)}">{STATUS_ICONS[st]}</td>'
            f'<td><a href="#{esc(_sample_anchor(samples, sample))}">{esc(sample)}</a></td>'
        )
        by_metric = {r.metric: r for r in reasons[i]}
        for col, _ in columns:
            cell = _compact_cell(row, col)
            hits = [by_metric[m] for m in COMPACT_CELL_METRICS.get(col, ()) if m in by_metric]
            if col == "breakpoint_adjacencies" and has_assembly and "breakpoints_private" in by_metric:
                hits.append(by_metric["breakpoints_private"])
            if hits:
                kind = STATUS_REVIEW if any(h.kind == STATUS_REVIEW for h in hits) else STATUS_UNUSUAL
                cls = "flag" if kind == STATUS_REVIEW else "unusual"
                ttl = esc(" | ".join(h.sentence for h in hits))
                w(f'<td class="{cls}" title="{ttl}">{cell}</td>')
            else:
                w(f"<td>{cell}</td>")
        w("</tr>\n")
    w("</table></div>\n")
    legend = (
        f'<p class="where"><span class="icon ok">{STATUS_ICONS[STATUS_OK]}</span> normal '
        f'&nbsp; <span class="icon review">{STATUS_ICONS[STATUS_REVIEW]}</span> needs review '
        f'&nbsp; <span class="icon unusual">{STATUS_ICONS[STATUS_UNUSUAL]}</span> unusual. '
        "Shaded cells: red = flagged, purple = unusual; hover for the reason. Aligned "
        "percentages count only columns where both genomes have a base."
    )
    if has_assembly:
        legend += " Telomeres = telomeric ends of the major (chromosome-scale) sequences / possible ends."
    elif any(s.get("query_length_source") == "maf_srcsize" for s in summaries):
        legend += (
            " &dagger; query % computed against MAF <code>srcSize</code> (no query FASTA or "
            "<code>.fai</code>), which overstates it."
        )
    w(legend + "</p>\n")

    # ── 5. shared features ───────────────────────────────────────────────────
    w('<h2 id="shared">Shared features</h2>\n')
    w(
        "<p>Breakpoints and nested blocks found in two or more samples at the same reference "
        f"location (ends within {recurrence_window_bp:,} bp, grouped transitively). Shared "
        "features usually reflect real shared structure (e.g. a common inversion) or "
        "structure or misassembly specific to the reference, rather than the quality of "
        "individual samples.</p>\n"
    )
    if n < 2:
        w('<p class="na">Needs at least two samples.</p>\n')
    elif not shared:
        w('<p class="na">No breakpoint or nested block is shared between samples.</p>\n')
    else:
        shown = shared[:MAX_HTML_SHARED_FEATURES]
        if len(shared) > len(shown):
            w(
                f'<p class="na">Showing the {len(shown):,} largest of {len(shared):,} shared '
                "features; the All breakpoints / All nested blocks sections and the combined "
                "TSVs list every row.</p>\n"
            )
        w('<div class="scroll"><table class="shared">\n<tr><th>Type</th><th>Reference region</th>'
          "<th>Samples</th><th>Which</th></tr>\n")
        for f in shown:
            kind = f.kind + (f" ({f.detail})" if f.detail else "")
            if len(f.samples) == n:
                which = "all"
            else:
                which = ", ".join(
                    f'<a href="#{esc(_sample_anchor(samples, s))}">{esc(s)}</a>' for s in f.samples
                )
            w(
                f'<tr><td class="left">{esc(kind)}</td><td class="left">{esc(f.region)}</td>'
                f"<td>{len(f.samples):,}/{n:,}</td><td class=\"left\">{which}</td></tr>\n"
            )
        w("</table></div>\n")

    # ── 6. drill-down sections (collapsed) ───────────────────────────────────
    w("<h2>Details</h2>\n")

    # Distributions
    w('<details class="drill" id="distributions">\n<summary>Distributions</summary>\n')
    w(
        '<p class="legend">One dot per sample (hover for name); dashed line = median.'
        f'<span style="background:{COLOR_NORMAL}"></span>normal'
        f'<span style="background:{COLOR_REVIEW}"></span>flagged (review)'
        f'<span style="background:{COLOR_UNUSUAL}"></span>unusual</p>\n'
    )
    skipped = flag_result.skipped
    if n >= min_samples:
        partial = {m: k for m, k in skipped.items() if k > 0}
        empty = [m for m, k in skipped.items() if k == 0]
        if partial:
            items = ", ".join(f"<code>{esc(m)}</code> (n={k:,})" for m, k in partial.items())
            w(
                f'<p class="note"><strong>Relative flagging skipped</strong> for metrics with '
                f"fewer than min_samples = {min_samples:,} non-NA values: {items}.</p>\n"
            )
        if empty:
            items = ", ".join(f"<code>{esc(m)}</code>" for m in empty)
            w(f'<p class="na">Not available for any sample (not flagged): {items}.</p>\n')
    for metric in FLAG_DIRECTIONS:
        values = [parse_number(s.get(metric)) for s in summaries]
        if all(v is None for v in values):
            continue
        w(svg_strip_plot(metric, samples, values, [metric in f for f in flags],
                         [metric in u for u in unusual]))
        w("\n")
    w("</details>\n")

    # All breakpoints
    w('<details class="drill" id="all-breakpoints">\n<summary>All breakpoints</summary>\n')
    w(
        "<p>A breakpoint is an adjacent pair of blocks along a query contig with a strand "
        "flip, a reference-contig jump, or out-of-order reference coordinates. "
        "Location on the query contig hints at the cause:</p>\n<ul>\n"
        "<li><code>contig_end</code>: near a query contig end &rArr; points to assembly "
        "or scaffolding (a join or contig boundary).</li>\n"
        "<li><code>n_gap</code>: next to an N run &rArr; a scaffold gap, i.e. a "
        "scaffolding join.</li>\n"
        "<li><code>interior</code>: inside contiguous sequence &rArr; real structural "
        "variation or misalignment; the two cannot be distinguished without reads.</li>\n"
        "</ul>\n"
        "<p><strong>Recurrence.</strong> A breakpoint is <em>recurrent</em> when another "
        "sample has one joining the same pair of reference positions (either orientation, "
        f"each end within {recurrence_window_bp:,} bp on the same reference contig). Shared "
        "breakpoints are more likely biological (e.g. a known inversion); private ones are "
        "more likely assembly or alignment artifacts. Recurrence needs at least two "
        "samples.</p>\n"
    )
    if not breakpoints_supplied:
        w('<p class="na">No breakpoint tables were supplied (<code>--breakpoints</code>).</p>\n')
    else:
        w(_mismatch_note(
            _count_mismatches(summaries, breakpoints, "breakpoint_adjacencies"),
            "Breakpoint", "breakpoint files",
        ))
    if breakpoints:
        n_recurrent = sum(1 for r in breakpoints if r.get("recurrent") == "true")
        if n >= 2:
            w(
                f"<p>{len(breakpoints):,} breakpoint(s) across samples; "
                f"{n_recurrent:,} recurrent. Recurrent rows are listed first and shaded.</p>\n"
            )
        else:
            w(f"<p>{len(breakpoints):,} breakpoint(s); recurrence not computed (one sample).</p>\n")
        w('<div class="scroll">')
        w(_breakpoint_table(
            _recurrence_sort(breakpoints), cap=MAX_HTML_BREAKPOINT_ROWS, with_sample=True,
            samples=samples, where="the combined breakpoints TSV",
        ))
        w("</div>\n")
    elif breakpoints_supplied:
        w('<p class="na">No breakpoints in any sample.</p>\n')
    w("</details>\n")

    # All nested blocks
    w('<details class="drill" id="all-nested">\n<summary>All nested blocks</summary>\n')
    w(
        "<p>A <em>nested block</em> is a block whose forward query interval lies inside a "
        "larger block's query interval: the same query sequence aligned in two places "
        "(a secondary or transposed alignment). "
        "Nested blocks are <strong>excluded from the breakpoint walk</strong> (otherwise they "
        "appear as spurious strand flips with very long junctions) and listed here, each "
        "with the reference location of its containing block.</p>\n"
        "<p><strong>Recurrence.</strong> A nested block is <em>recurrent</em> when another "
        "sample has a nested block on the same reference contig with start and end each "
        f"within {recurrence_window_bp:,} bp. Nested blocks shared across samples at the "
        "same reference coordinates may reflect reference-specific structure (e.g. a "
        "duplication or transposition in the reference) or a reference misassembly; private "
        "ones are more likely repeats or alignment artefacts. Nested blocks are reported, "
        "not flagged.</p>\n"
    )
    if not nested_supplied:
        w('<p class="na">No nested-block tables were supplied (<code>--nested-blocks</code>).</p>\n')
    else:
        w(_mismatch_note(
            _count_mismatches(summaries, nested, "nested_blocks"),
            "Nested-block", "nested-block files",
        ))
    if nested:
        n_recurrent = sum(1 for r in nested if r.get("recurrent") == "true")
        if n >= 2:
            w(
                f"<p>{len(nested):,} nested block(s) across samples; {n_recurrent:,} "
                "recurrent. Recurrent rows are listed first and shaded.</p>\n"
            )
        else:
            w(f"<p>{len(nested):,} nested block(s); recurrence not computed (one sample).</p>\n")
        w('<div class="scroll">')
        w(_nested_table(
            _recurrence_sort(nested), cap=MAX_HTML_NESTED_ROWS, with_sample=True,
            samples=samples, where="the combined nested-blocks TSV",
        ))
        w("</div>\n")
    elif nested_supplied:
        w('<p class="na">No nested blocks in any sample.</p>\n')
    w("</details>\n")

    # Assembly metrics (or MAF-only contiguity)
    heading = "Assembly metrics" if has_assembly else "Contiguity (from the MAF only)"
    w(f'<details class="drill" id="assembly">\n<summary>{esc(heading)}</summary>\n')
    if has_assembly:
        w(
            "<p>Query assembly statistics computed from the FASTA given to "
            "<code>maf_stats.py --fasta</code>; independent of alignment. Samples run "
            "without a FASTA show NA.</p>\n"
        )
        w(_sample_table(summaries, flags, samples, ASSEMBLY_TABLE_COLUMNS, "assembly",
                        flags_column=False, unusual=unusual))
    else:
        w(
            '<p class="na">No sample has assembly statistics (run <code>maf_stats.py '
            "--fasta</code> with the query assembly to get contiguity, N content, GC, "
            "soft-masking, telomeres and unaligned contigs).</p>\n"
        )
    w(
        '<h3>Contiguity from the MAF</h3>\n<p class="na">From the MAF: aligned contigs only, a '
        "lower bound on fragmentation. Unaligned contigs are invisible here. For "
        "chromosome-scale assemblies the aligned-contig N50 is essentially a chromosome "
        "length and says little about contiguity; use the contig N50 above.</p>\n"
    )
    w(_sample_table(summaries, flags, samples, MAF_CONTIGUITY_COLUMNS, "maf-contiguity",
                    flags_column=False))
    w("</details>\n")

    # Per-sample detail
    w('<details class="drill" id="per-sample">\n<summary>Per-sample detail</summary>\n')
    w(
        "<p>All alignment metrics, one row per sample. Flagged cells red, unusual cells "
        "purple; hover for the raw reason. &dagger; marks query percentages computed "
        "against MAF <code>srcSize</code> (no query FASTA or <code>.fai</code>), which "
        "overstates them.</p>\n"
    )
    w(_sample_table(summaries, flags, samples, PRIMARY_TABLE_COLUMNS, "alignment",
                    flags_column=True, unusual=unusual))
    for i, row in enumerate(summaries):
        sample = row["sample"]
        st = statuses[i]
        w(f'<details id="{esc(_sample_anchor(samples, sample))}">\n')
        aligned = fmt_value(row.get("aligned_reference_pct"), "aligned_reference_pct")
        ident = fmt_value(row.get("identity_pct"), "identity_pct")
        pct = lambda cell: cell if "NA" in cell else cell + "%"  # noqa: E731
        w(
            f'<summary><span class="icon {st}">{STATUS_ICONS[st]}</span> {esc(sample)}'
            f" &mdash; {esc(STATUS_LABELS[st])}"
            f" &nbsp;|&nbsp; aligned ref {pct(aligned)}"
            f" &nbsp;|&nbsp; identity {pct(ident)}"
            f" &nbsp;|&nbsp; {fmt_value(row.get('breakpoint_adjacencies'), 'breakpoint_adjacencies')}"
            f" breakpoint adj.</summary>\n"
        )
        if reasons[i]:
            w('<section class="reasons">\n<h3>Why</h3>\n')
            w(_reason_list(reasons[i]))
            w("</section>\n")

        w('<table class="sample-summary">\n<tr><th>Metric</th><th>Value</th></tr>\n')
        w(f'<tr><td>status</td><td class="status {st}">{esc(st)}</td></tr>\n')
        for col in OUTPUT_SUMMARY_COLUMNS[1:]:
            flagged = flags[i].get(col)
            odd = unusual[i].get(col)
            cell = fmt_value(row.get(col), col)
            if flagged:
                tip = "; ".join(f"{col}:{r}" for r in flagged)
                w(f'<tr><td>{esc(col)}</td><td class="flag" title="{esc(tip)}">{cell}</td></tr>\n')
            elif odd:
                tip = "; ".join(f"{col}:{r}" for r in odd)
                w(f'<tr><td>{esc(col)}</td><td class="unusual" title="{esc(tip)}">{cell}</td></tr>\n')
            else:
                w(f"<tr><td>{esc(col)}</td><td>{cell}</td></tr>\n")
        w("</table>\n")

        sample_bps = bp_by_sample.get(sample, [])
        w('<section class="breakpoints">\n<h3>Breakpoints</h3>\n')
        if sample_bps:
            w('<div class="scroll">')
            w(_breakpoint_table(
                _recurrence_sort(sample_bps), cap=MAX_HTML_CONTIG_ROWS, with_sample=False,
                samples=samples, where="the combined breakpoints TSV",
            ))
            w("</div>\n")
        else:
            w('<p class="na">No breakpoints for this sample.</p>\n')
        w("</section>\n")

        sample_nested = nested_by_sample.get(sample, [])
        w('<section class="nested">\n<h3>Nested alignments</h3>\n')
        if sample_nested:
            w('<div class="scroll">')
            w(_nested_table(
                _recurrence_sort(sample_nested), cap=MAX_HTML_CONTIG_ROWS, with_sample=False,
                samples=samples, where="the combined nested-blocks TSV",
            ))
            w("</div>\n")
        else:
            w('<p class="na">No nested blocks for this sample.</p>\n')
        w("</section>\n")

        ref_rows = _breakpoint_first(by_reference.get(sample, []))
        w('<section class="ref-contigs">\n<h3>Reference contigs</h3>\n')
        if ref_rows:
            w("<p>Contigs with breakpoint adjacencies are listed first and shaded.</p>\n")
            w('<div class="scroll">')
            w(_contig_table(ref_rows, REFERENCE_TABLE_COLUMNS, "ref-table"))
            w("</div>\n")
        else:
            w('<p class="na">No reference-contig table supplied.</p>\n')
        w("</section>\n")

        query_rows = _breakpoint_first(by_query.get(sample, []))
        w('<section class="query-contigs">\n<h3>Query contigs</h3>\n')
        if query_rows:
            w("<p>Contigs with breakpoint adjacencies are listed first and shaded.</p>\n")
            w('<div class="scroll">')
            w(_contig_table(query_rows, QUERY_TABLE_COLUMNS, "query-table"))
            w("</div>\n")
        else:
            w('<p class="na">No query-contig table supplied.</p>\n')
        w("</section>\n")

        plots = [r for r in ref_rows if (r.get("dotplot") or "").strip()]
        w('<section class="dotplots">\n<h3>Dotplots</h3>\n')
        if plots:
            w('<div class="thumbs">\n')
            for r in plots:
                bp = fmt_value(r.get("breakpoint_adjacencies"), "breakpoint_adjacencies")
                contig = esc(r["reference_contig"])
                w(_thumb(r["dotplot"].strip(), r["reference_contig"], thumbnail_width,
                         f"{contig} ({bp} breakpoint adj.)"))
            w("</div>\n")
        else:
            w('<p class="na">No dotplots drawn for this sample.</p>\n')
        w("</section>\n")
        w("</details>\n")
    w("</details>\n")

    # Metric definitions
    w('<details class="drill" id="definitions">\n<summary>Metric definitions</summary>\n')
    w('<dl class="defs">\n')
    for name, definition in STATUS_DEFINITIONS:
        w(f"<dt>{esc(name)}</dt><dd>{definition}</dd>\n")
    w("</dl>\n")
    w(
        '<p class="note"><strong>Aligned versus block-span metrics.</strong> '
        "<code>aligned_*</code> counts bases in columns where both rows have a base. "
        "<code>block_span_*</code> is the union of block intervals and includes long indels "
        "inside blocks; AnchorWave blocks span chromosomes between collinear anchors, so "
        "block-span coverage overstates alignment and is not flagged. Identity and indel "
        "column percentages are raw block-column metrics: overlapping blocks count more "
        "than once (<code>overlapping_reference_bp</code> makes that visible). "
        "<strong>Alignment versus assembly.</strong> Assembly metrics describe the query "
        "FASTA itself and need <code>--fasta</code>; the MAF alone cannot see unaligned "
        "contigs.</p>\n"
    )
    w('<dl class="defs">\n')
    for name, definition in METRIC_DEFINITIONS:
        direction = FLAG_DIRECTIONS.get(name)
        extra = ""
        if direction:
            other = "high" if direction == "low" else "low"
            extra = (f" <em>Flagged when unusually {direction}; marked unusual when "
                     f"unusually {other}.</em>")
        # Definitions are static trusted HTML; names are escaped.
        w(f"<dt>{esc(name)}</dt><dd>{definition}{extra}</dd>\n")
    w("</dl>\n")
    w("</details>\n")

    # Parameters
    w('<details class="drill" id="parameters">\n<summary>Parameters</summary>\n')
    w('<table class="params">\n<tr><th>Parameter</th><th>Value</th></tr>\n')
    w(f"<tr><td>min_samples</td><td>{min_samples:,}</td></tr>\n")
    w(f"<tr><td>z_threshold</td><td>{esc(f'{z_threshold:g}')}</td></tr>\n")
    w(f"<tr><td>recurrence_window_bp</td><td>{recurrence_window_bp:,}</td></tr>\n")
    for dest, metric in THRESHOLD_OPTIONS.items():
        option = "--" + dest.replace("_", "-")
        value = thresholds.get(metric)
        shown = esc(f"{value:g}") if value is not None else '<span class="na">off</span>'
        w(f"<tr><td>{esc(option)}</td><td>{shown}</td></tr>\n")
    for col in ("min_block_bp", "overlap_tolerance_bp", "breakpoint_context_bp"):
        distinct = list(dict.fromkeys(s[col] for s in summaries))
        shown = ", ".join(fmt_value(v, col) for v in distinct)
        if len(distinct) > 1:
            shown += ' &nbsp;<strong>(differs between samples)</strong>'
        w(f"<tr><td>{esc(col)}</td><td>{shown}</td></tr>\n")
    w("</table>\n")
    w(
        "<p>Relative flags: direction-aware robust z = (value &minus; median) / scale, "
        f"scale = max({MAD_TO_SD} &times; MAD, metric floor). A value is flagged (review) "
        f"when |z| &gt; {esc(f'{z_threshold:g}')} in the bad direction and marked unusual "
        "when |z| exceeds it in the good direction. <code>(mad0_floor)</code> "
        "marks z-scores computed with MAD = 0, where the scale is the floor alone. Metrics "
        "that are NA for a sample (e.g. assembly metrics without <code>--fasta</code>) are "
        "not flagged for that sample.</p>\n"
    )
    w("</details>\n")
    w("</body>\n</html>\n")
    return "".join(out)


def _sample_anchor(samples: list[str], sample: str) -> str:
    # Index-based anchors are unique and safe for any sample name.
    return f"s-{samples.index(sample)}"


# ── driver ───────────────────────────────────────────────────────────────────


def load_inputs(
    summary_paths: list[Path],
    by_reference_paths: list[Path],
    by_query_paths: list[Path],
    breakpoint_paths: list[Path] | None = None,
    nested_paths: list[Path] | None = None,
) -> tuple[
    list[dict[str, str]],
    dict[str, list[dict[str, str]]],
    dict[str, list[dict[str, str]]],
    dict[str, list[dict[str, str]]],
    dict[str, list[dict[str, str]]],
]:
    """Return ``(summaries, by_reference, by_query, breakpoints, nested)``.

    Per-sample tables are grouped by their ``sample`` column; a sample with no
    breakpoint (nested-block) file has no breakpoints (nested blocks).
    """
    if not summary_paths:
        raise ValueError("No summary inputs given")
    summaries: list[dict[str, str]] = []
    seen: dict[str, Path] = {}
    for path in summary_paths:
        row = read_summary(path)
        sample = row["sample"]
        if sample in seen:
            raise ValueError(
                f"Duplicate sample name {sample!r} in {path} (also in {seen[sample]})"
            )
        seen[sample] = path
        summaries.append(row)
    samples = [s["sample"] for s in summaries]
    by_reference = read_contig_tables(
        by_reference_paths, BY_REFERENCE_COLUMNS, samples, "by-reference-contig"
    )
    by_query = read_contig_tables(by_query_paths, BY_QUERY_COLUMNS, samples, "by-query-contig")
    breakpoints = read_contig_tables(
        list(breakpoint_paths or []), BREAKPOINT_COLUMNS, samples, "breakpoints"
    )
    nested = read_contig_tables(
        list(nested_paths or []), NESTED_COLUMNS, samples, "nested-blocks"
    )
    return summaries, by_reference, by_query, breakpoints, nested


def build_report(
    summary_paths: list[Path],
    by_reference_paths: list[Path],
    by_query_paths: list[Path],
    out_tsv: Path,
    out_html: Path,
    *,
    out_breakpoints_tsv: Path,
    out_nested_tsv: Path,
    breakpoint_paths: list[Path] | None = None,
    nested_paths: list[Path] | None = None,
    recurrence_window_bp: int = DEFAULT_RECURRENCE_WINDOW_BP,
    min_samples: int = DEFAULT_MIN_SAMPLES,
    z_threshold: float = DEFAULT_Z_THRESHOLD,
    thresholds: dict[str, float] | None = None,
) -> FlagResult:
    thresholds = dict(thresholds or {})
    breakpoint_paths = list(breakpoint_paths or [])
    nested_paths = list(nested_paths or [])
    summaries, by_reference, by_query, breakpoints, nested = load_inputs(
        summary_paths, by_reference_paths, by_query_paths, breakpoint_paths, nested_paths
    )
    samples = [s["sample"] for s in summaries]
    annotated = annotate_recurrence(samples, breakpoints, recurrence_window_bp)
    nested_annotated = annotate_nested_recurrence(samples, nested, recurrence_window_bp)
    counts = recurrence_counts(samples, annotated)
    nested_counts = recurrence_counts(samples, nested_annotated, NESTED_EXTRA_COLUMNS)
    for row in summaries:
        row.update(counts[row["sample"]])
        row.update(nested_counts[row["sample"]])

    flag_result = compute_flags(
        summaries, min_samples=min_samples, z_threshold=z_threshold, thresholds=thresholds
    )
    write_summary_tsv(out_tsv, summaries, flag_result.flags, flag_result.unusual)
    write_breakpoints_tsv(out_breakpoints_tsv, annotated)
    write_nested_tsv(out_nested_tsv, nested_annotated)
    page = render_html(
        summaries,
        flag_result,
        by_reference,
        by_query,
        annotated,
        nested_annotated,
        min_samples=min_samples,
        z_threshold=z_threshold,
        thresholds=thresholds,
        recurrence_window_bp=recurrence_window_bp,
        breakpoints_supplied=bool(breakpoint_paths),
        nested_supplied=bool(nested_paths),
    )
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(page, encoding="utf-8")
    return flag_result


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Build the cross-sample MAF QC TSVs and HTML report.")
    ap.add_argument("--summaries", nargs="+", required=True,
                    help="Per-sample <sample>.maf_stats.tsv files.")
    ap.add_argument("--by-reference", nargs="*", default=[],
                    help="Per-sample <sample>.by_reference_contig.tsv files.")
    ap.add_argument("--by-query", nargs="*", default=[],
                    help="Per-sample <sample>.by_query_contig.tsv files.")
    ap.add_argument("--breakpoints", nargs="*", default=[],
                    help="Per-sample <sample>.breakpoints.tsv files (matched by sample column).")
    ap.add_argument("--nested-blocks", nargs="*", default=[],
                    help="Per-sample <sample>.nested_blocks.tsv files (matched by sample column).")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--out-breakpoints-tsv", required=True,
                    help="Combined breakpoint table with recurrence columns.")
    ap.add_argument("--out-nested-tsv", required=True,
                    help="Combined nested-block table with recurrence columns.")
    ap.add_argument("--recurrence-window-bp", type=int, default=DEFAULT_RECURRENCE_WINDOW_BP,
                    help="Max distance between paired breakpoint ends (or nested-block "
                         "reference starts/ends) on the same reference contig for two "
                         "samples' breakpoints (nested blocks) to match.")
    ap.add_argument("--min-samples", type=int, default=DEFAULT_MIN_SAMPLES,
                    help="Minimum samples with a value before relative flags are computed.")
    ap.add_argument("--z-threshold", type=float, default=DEFAULT_Z_THRESHOLD)
    ap.add_argument("--flag-min-aligned-reference", type=float, default=None, metavar="PCT")
    ap.add_argument("--flag-min-identity", type=float, default=None, metavar="PCT")
    ap.add_argument("--flag-max-breakpoints-per-gb", type=float, default=None, metavar="X")
    args = ap.parse_args(argv)
    if args.min_samples < 1:
        ap.error("--min-samples must be >= 1")
    if not args.z_threshold > 0:
        ap.error("--z-threshold must be > 0")
    if args.recurrence_window_bp < 0:
        ap.error("--recurrence-window-bp must be >= 0")
    return args


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    thresholds = {
        metric: getattr(args, dest)
        for dest, metric in THRESHOLD_OPTIONS.items()
        if getattr(args, dest) is not None
    }
    try:
        build_report(
            [Path(p) for p in args.summaries],
            [Path(p) for p in args.by_reference],
            [Path(p) for p in args.by_query],
            Path(args.out_tsv),
            Path(args.out_html),
            out_breakpoints_tsv=Path(args.out_breakpoints_tsv),
            out_nested_tsv=Path(args.out_nested_tsv),
            breakpoint_paths=[Path(p) for p in args.breakpoints],
            nested_paths=[Path(p) for p in args.nested_blocks],
            recurrence_window_bp=args.recurrence_window_bp,
            min_samples=args.min_samples,
            z_threshold=args.z_threshold,
            thresholds=thresholds,
        )
    except ValueError as exc:
        raise SystemExit(f"maf_stats_report: error: {exc}") from exc


if __name__ == "__main__":
    main()
