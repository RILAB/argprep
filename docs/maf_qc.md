# MAF alignment QC (`maf_qc`)

Run this QC-only target before the more expensive per-contig workflow:

```bash
sbatch profiles/slurm/run-controller.sbatch options.yaml maf_qc
# Inspect results/maf_stats/maf_stats.html, then start the main run separately:
sbatch profiles/slurm/run-controller.sbatch options.yaml
```

`maf_qc` runs one job per sample with `scripts/maf_stats.py`, followed by one
small cross-sample report job using `scripts/maf_stats_report.py`. It is not
part of `rule all`: completing QC does not start the main analysis, and the
main analysis neither runs nor waits for QC. Constructing the QC DAG does not
read the MAFs. The main workflow's automatic contig discovery does (see
[issue #36](https://github.com/RILAB/argprep/issues/36)).

Each sample's MAF is read once. QC keeps **alignment quality** and **assembly
quality** separate. Alignment metrics come from the MAF alone. Assembly metrics
require the query assembly FASTA through `query_fasta_dir`; without it they are
`NA` and the MAF-only view is used.

## Alignment metrics

- **Aligned bases:** reference and query bases opposite a base in the other
  row, in bp and as a percentage of the genome. This is the headline coverage
  figure. **Block span** (the union of block intervals) is also reported, but
  AnchorWave blocks can span long indels and substantially overstate alignment.
  In a test of two teosinte genomes against B73 v5, blocks spanned about 94% of
  the reference while only about 36% of reference bases aligned base-to-base.
- **Identity:** matches divided by columns where both bases are A/C/G/T,
  case-insensitively.
- **Insertion and deletion columns:** query bases opposite reference gaps and
  reference bases opposite query gaps, each as a percentage of alignment
  columns.
- **Blocks and overlap:** block count plus `overlapping_reference_bp` and
  `overlapping_query_bp`, the bases covered by more than one block.
- **Nested blocks:** blocks whose query interval lies inside a larger block's
  interval. These are listed in `<sample>.nested_blocks.tsv` and excluded from
  the breakpoint walk to avoid creating false strand flips. The report also
  checks whether nested blocks recur across samples at the same reference
  coordinates.
- **Structural breakpoints:** strand flips, reference-contig jumps, and
  out-of-order blocks are counted separately and together per Gb aligned.
  Breakpoints are classified as `contig_end`, `n_gap` (requires FASTA), or
  `interior`, and checked for recurrence at matching reference coordinates in
  other samples. Shared breakpoints are more likely biological; private ones
  are more likely alignment or assembly artifacts, although read evidence is
  needed to distinguish them confidently.
- **MAF-only contiguity:** number and N50 of query contigs with at least one
  alignment. Because unaligned contigs are absent from the MAF, this is a lower
  bound on assembly fragmentation.
- **Dotplots:** optional per-reference-contig alignment plots.

Identity, aligned-base, and insertion/deletion counts are summed over columns,
so overlapping or secondary blocks are counted more than once. Block span is
union-based; `overlapping_reference_bp` reports the duplicated coverage.

## Assembly metrics

These metrics require `query_fasta_dir`:

- **Size and contiguity:** total length and record count. Scaffold N50/L50 is
  calculated over FASTA records; contig N50/L50 is calculated after splitting
  records at N gaps.
- **Ns:** total N bases and the number of N gaps (runs of at least 10 Ns).
- **GC and soft masking:** soft-masked percentage is `NA` when the FASTA
  contains no lowercase sequence.
- **Telomeres:** a sequence end is telomeric when at least half of its terminal
  1 kb is plant telomere repeat (`TTTAGGG`/`CCCTAAA`). Ends are summarized on
  major sequences (at least 10% as long as the longest sequence); telomeric
  minor sequences are counted separately.
- **Unaligned sequences:** size and N/soft-masked content of sequences without
  alignments. FASTA lengths also become the query-coverage denominator.

## Strict input validation

The parser stops with file and block context for:

- blocks without exactly one reference and one query row;
- size fields that disagree with the ungapped sequence;
- coordinates outside `srcSize`;
- inconsistent `srcSize` values;
- unknown or ambiguous reference contigs after exact or normalized matching;
- minus-strand reference rows; or
- an empty MAF.

## Outputs

Outputs are written under `results/maf_stats/`:

- `maf_stats.html`: flagged alignment and assembly tables, recurrent
  breakpoints, metric plots, per-sample drill-downs, dotplots, and definitions.
- `maf_stats.tsv`: one row per sample with genome-wide metrics, recurrent and
  private breakpoint counts, and flags.
- `maf_stats.breakpoints.tsv`: every breakpoint with location and recurrence.
- `maf_stats.nested_blocks.tsv`: every nested block with recurrence.
- `<sample>.maf_stats.tsv`, `<sample>.by_reference_contig.tsv`,
  `<sample>.by_query_contig.tsv`, `<sample>.breakpoints.tsv`, and
  `<sample>.nested_blocks.tsv`: per-sample details.
- `dotplots/<sample>/<reference_contig>.png`: 800×800 dotplots, with forward
  blocks in blue, reverse blocks in red, and query contigs in labelled bands.

## Flagging

**Relative flags:** a sample is flagged when its robust z-score (median/MAD)
passes 3.5 in the adverse direction. Low values are adverse for aligned
reference percentage, aligned query percentage, identity, and assembly contig
N50. High values are adverse for overlapping query bp, breakpoints per Gb,
private breakpoints, assembly N percentage, and unaligned-contig bp. Scale
floors prevent near-constant cohorts from being flagged for trivial
differences. Relative flags require at least five samples with a value for that
metric.

**Absolute thresholds:** optional and disabled by default because expected
identity and fragmentation depend on divergence.

## Configuration

All keys are optional:

| Key | Default | Meaning |
|---|---|---|
| `query_fasta_dir` | unset | Query FASTAs named `<sample>.fa`, `.fasta`, or `.fna`, optionally gzipped. Enables assembly metrics and `n_gap` breakpoint classification. Every sample must have exactly one matching FASTA. |
| `query_fai_dir` | unset | Length-only alternative to FASTAs. Accepts `<sample>.fai`, `<sample>.fa.fai`, `<sample>.fasta.fai`, `<sample>.fna.fai`, `<sample>.fa.gz.fai`, or `<sample>.fasta.gz.fai`, and is mutually exclusive with `query_fasta_dir`. Without either option, query coverage uses only contigs present in the MAF and therefore overstates coverage; the report marks this with †. |
| `maf_stats_min_block_bp` | `0` | Ignore shorter reference spans when counting breakpoints. |
| `maf_stats_overlap_tolerance_bp` | `0` | Allowed reference overlap between adjacent same-strand blocks before they count as out of order. |
| `maf_stats_breakpoint_context_bp` | `1000000` | Distance from a query contig end or N gap used for breakpoint classification. |
| `maf_stats_recurrence_window_bp` | `500000` | Coordinate tolerance at both junction ends when matching breakpoints across samples. |
| `maf_stats_dotplots` | `flagged` | `flagged`, `all`, or `false`. |
| `maf_stats_dotplot_max` | `20` | Maximum flagged plots per sample; `0` means unlimited. |
| `maf_stats_flag_min_aligned_reference`, `maf_stats_flag_min_identity`, `maf_stats_flag_max_breakpoints_per_gb` | unset | Optional absolute thresholds. |
| `maf_stats_threads`, `maf_stats_mem_mb`, `maf_stats_time` | `1`, `16000`, `12:00:00` | Per-sample resources. These remain provisional: before block comparisons were chunked, one 7.7 GB AnchorWave MAF took about one minute and peaked near 10 GB; the current implementation has not yet been remeasured, with or without FASTA input. |
| `maf_stats_report_mem_mb`, `maf_stats_report_time` | `4000`, `00:30:00` | Cross-sample report resources. |

## Standalone commands

Run statistics for one sample (`--fasta` is optional):

```bash
python scripts/maf_stats.py --maf S.maf.gz --reference-fai ref.fa.fai \
  --sample S --out-dir qc [--fasta S.fa.gz]
```

Generate dotplots directly:

```bash
python scripts/maf_dotplot.py --maf S.maf.gz --reference-fai ref.fa.fai \
  --out-dir plots
```

The dotplot command exits nonzero if any input fails.
