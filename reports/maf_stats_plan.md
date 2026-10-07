# Plan: first-stage per-sample MAF quality control

Written 2026-10-06 and revised after design review. File and line references must be
rechecked before implementation.

## Goal and user workflow

Add a standalone first-stage QC workflow that scans each input MAF, computes alignment
statistics, creates structural dotplots, and produces one cross-sample TSV and HTML report.
Users run and inspect this stage before launching the main per-contig pipeline.

```bash
mkdir -p logs/slurm
sbatch profiles/slurm/run-controller.sbatch options.yaml maf_qc
# Inspect results/maf_stats/maf_stats.html, then run the main workflow separately:
sbatch profiles/slurm/run-controller.sbatch options.yaml
```

`maf_qc` is intentionally not part of `rule all`, and the main workflow does not depend on
it. The two commands are separate user decisions: a completed QC run never automatically
starts the expensive analysis. The recommended controller wrapper runs Snakemake itself in a
non-preemptable Slurm allocation and forwards `maf_qc` as a target. Within that controller,
only genuinely trivial metadata operations may be `localrules`; every MAF/BED scan,
concatenation, index, report, and plotting operation is a separate rule with explicit Slurm
resources. This policy also keeps direct controller invocation from moving substantive work
onto a login node.

## User-facing outputs

The primary outputs are deliberately simple:

- `maf_stats/maf_stats.tsv` — one row per sample with genome-wide metrics and flags.
- `maf_stats/maf_stats.html` — the same summary as a readable report, with expandable
  per-sample details, metric definitions, flagged values, and dotplot thumbnails that link
  to full-resolution PNGs.

Diagnostic files support drill-down without mixing reference- and query-contig row types:

- `maf_stats/<sample>.maf_stats.tsv` — one genome-wide row for one sample.
- `maf_stats/<sample>.by_reference_contig.tsv` — coverage, identity, gaps, and block metrics
  grouped by reference contig.
- `maf_stats/<sample>.by_query_contig.tsv` — alignment and breakpoint metrics grouped by
  query contig.
- `maf_stats/dotplots/<sample>/<reference_contig>.png` — structural dotplots selected by
  the configured dotplot mode.

The HTML report is the normal entry point. The two contig-level TSVs exist to answer such
questions as “which reference chromosome has low coverage?” and “which query scaffold
contains the chromosome jump?”

## Metric definitions

### Genome-wide summary (`maf_stats.tsv`)

| Metric | Definition |
|---|---|
| Reference coverage (%) | Union of reference intervals covered by any accepted block, divided by total reference length from the reference `.fai`. |
| Covered reference bp | Numerator of reference coverage; preferred name for the merged block span. |
| Reference bp aligned to a query base | Reference bases opposite a non-gap query character; reported separately from covered reference bp. |
| Query coverage (%) | Union of forward-coordinate query intervals divided by query genome size from a query `.fai`. If no query `.fai` is supplied, use only distinct query contigs represented in the MAF and label the denominator basis explicitly. |
| Sequence identity (%) | Matches divided by columns where both bases are A/C/G/T, case-insensitive; N and other ambiguity codes are excluded. |
| Gap fraction (%) | Columns with `-` in either row divided by all alignment columns. |
| Number of blocks | Count of accepted pairwise MAF blocks. |
| Block N50 | N50 of reference `size` across accepted blocks. |
| Strand flips | Adjacent accepted blocks along a query contig whose query strands differ. |
| Reference-contig jumps | Adjacent accepted blocks along a query contig aligned to different reference contigs. |
| Out-of-order adjacencies | Same-reference, same-strand adjacent blocks whose reference coordinates contradict query order beyond the configured tolerance. |
| Unique breakpoint adjacencies | Number of adjacent block pairs having any of the three breakpoint properties. |
| Breakpoints per Gb covered | Unique breakpoint adjacencies divided by covered reference bp, multiplied by 1e9. |

Genome-wide identity, gap fraction, block count, and N50 are raw block-column metrics and
can count overlapping or secondary blocks more than once. Coverage is union-based and does
not. The TSV and HTML definitions must state this distinction; the implementation should
also report `overlapping_reference_bp` so unexpectedly duplicated alignment content is
visible rather than silently inflating raw metrics.

### Reference-contig detail

One row per reference contig: reference length, covered reference bp, coverage percentage,
reference bp aligned to a query base, identity numerator/denominator and percentage, gap
columns/total columns and percentage, block count, block N50, and overlapping reference bp.

### Query-contig detail and breakpoint classification

Convert every query interval to forward coordinates. For a `-` row:

```text
forward_start = srcSize - start - size
forward_end   = srcSize - start
```

Group accepted blocks by query contig and sort them deterministically by forward query start,
forward query end, reference contig, reference start, reference end, query strand, and finally
the original block ordinal. Walk adjacent pairs in that order. Record the three breakpoint
properties independently so a chromosome jump combined with a strand flip is not hidden.
Also record a single `has_any_breakpoint` boolean for the unique total.

For same-reference, same-strand blocks, classify out-of-order/overlap beyond tolerance as:

- `+`: `next_ref_start < previous_ref_end - overlap_tolerance_bp`
- `-`: `next_ref_end > previous_ref_start + overlap_tolerance_bp`

The query-contig TSV contains query length and denominator source, covered query bp, block
count, each breakpoint-property count, unique breakpoint adjacencies, and breakpoints per Gb.
Blocks shorter than `min_block_bp` are excluded from breakpoint classification only, not
from coverage or general alignment metrics. Defaults remain zero so small-block noise is
visible until real-data validation supports a different cutoff.

## Shared parsing versus QC validation

Move/refactor the reusable MAF record and block iterator from `scripts/maf_to_sites.py` into
`scripts/common.py`; do not maintain two parsers. The shared iterator must stream plain and
gzip-compressed MAFs, preserve every `s` row, and retain its existing basic syntax and
equal-alignment-length checks. Moving it must not make the main pipeline reject inputs it
currently accepts.

QC applies a separate strict pairwise-validation layer after parsing. That layer validates
rather than silently reinterpreting unexpected data:

- Each accepted block must contain exactly the expected pairwise reference and query `s`
  rows. Extra or missing sequence rows are an error with file and block context.
- Aligned strings must have equal lengths, and each row's ungapped sequence length must
  equal its MAF `size` field.
- Coordinates must be non-negative and contained within `srcSize`; strand must be `+` or
  `-`.
- Reference names are checked against the reference `.fai`, using the pipeline's existing
  contig normalization rules only where the mapping is unique. A MAF reference contig absent
  from the `.fai` is an error, not a warning.
- Repeated `srcSize` values for the same source must agree. Query `.fai` lengths, when
  provided, must agree with the corresponding MAF sources.
- Empty MAFs and MAFs with zero accepted blocks produce clear QC failures, not empty reports
  that appear successful.

## 1. Shared contig-manifest preflight

Resolve the production-scale planning problem tracked separately in issue #36 before building
QC-specific discovery. Add a small per-sample Slurm rule that scans each raw MAF for reference
contig names and writes a compact contig-list file, followed by a manifest rule/checkpoint that
intersects and resolves those lists against the reference `.fai`.

Both the default workflow and `maf_qc` consume the resolved manifest. The manifest stage is
independent infrastructure, not a QC output: running the main workflow must not cause the QC
statistics or report to run. This removes whole-MAF discovery scans from Snakemake DAG
construction and avoids implementing two different contig-discovery paths. The lightweight
manifest scan and the later full QC scan both read each MAF; that deliberate duplicate pass is
preferable to reading all MAFs on the controller host.

## 2. `scripts/maf_stats.py` — one sample, one streaming pass

```bash
python scripts/maf_stats.py \
    --maf S.maf[.gz] \
    --reference-fai ref.fa.fai \
    --sample S \
    --out-dir results/maf_stats \
    [--query-fai S.fa.fai] \
    [--min-block-bp 0] \
    [--overlap-tolerance-bp 0] \
    [--dotplots flagged|all|false]
```

Use the shared iterator plus the QC-only strict validator described above. Sequence identity
and gap calculations may compare upper-cased NumPy byte arrays, but benchmark against direct
byte operations before assuming NumPy is faster for typical blocks.

Collect metric counters, mergeable reference/query intervals, breakpoint coordinates, and
dotplot segments during the same scan. Dotplot generation must not rescan a multi-GB MAF.
Expected memory is the largest alignment block plus stored interval/segment coordinates;
benchmark this on representative MAFs before choosing resource defaults. If coordinate
storage is too large, spill compact per-contig coordinates to node-local `$TMPDIR` and render
one contig at a time.

This deliberately scans the raw MAF during the separate QC stage. The later main workflow
will read it again to split by contig; that duplicate read is the cost of keeping QC a clear,
independent approval point. Within QC, statistics and plots still share one scan.

## 3. Dotplots — refactor `scripts/maf_dotplot.py`

Use the existing untracked `scripts/maf_dotplot.py` as starter code, not as-is. Preserve its
forward/reverse coordinate logic and blue/red visual convention, while addressing these
requirements:

- Use the shared MAF parser followed by the QC-only validator, and support `.maf.gz`; do not
  invoke `awk` directly on compressed files.
- Do not assume every two `s` lines form a valid pair without block validation.
- Produce browser-friendly PNGs with a predictable pixel size. Optionally retain a PDF output
  later, but it is not required for the first version.
- Rasterize large block collections and use batched line collections rather than one plotting
  call per block.
- A reference contig may align to multiple query contigs. Stack query contigs into labeled
  y-axis bands (using query `.fai` order when available) rather than overlaying unrelated
  coordinate systems.
- `flagged` mode ranks reference contigs by breakpoint rate, applies both a configurable
  minimum absolute breakpoint count and minimum rate, and renders at most
  `maf_stats_max_dotplots_per_sample` contigs in deterministic rank order. Reference contigs
  on either side of a jump are eligible. `all` renders every aligned reference contig and
  `false` disables plotting. Start with bounded `flagged` mode as the provisional default;
  real-MAF validation must set meaningful count/rate floors so tiny-block noise does not make
  `flagged` equivalent to `all`.

Matplotlib must be added to `argprep.yml`. Plotting occurs in each sample's Slurm QC job so
samples are parallelized by the scheduler. Do not create an internal process pool by default;
one Snakemake job should use its declared `threads` and memory without oversubscribing a node.

Fully interactive block-level plots are out of scope for the first version. Static PNG
thumbnails with links to full-resolution images provide useful drill-down while keeping the
HTML small and dependency-free.

## 4. `scripts/maf_stats_report.py` — cross-sample report

Merge per-sample summary TSVs into `maf_stats/maf_stats.tsv`, one row per sample, plus a
machine-readable `flags` column. Generate `maf_stats/maf_stats.html` with embedded CSS and
inline summary SVGs using the style of `scripts/summary_report.py`, with no JavaScript library;
dotplot PNGs remain linked assets beside the report.

The HTML contains:

- A primary one-row-per-sample table with highlighted flagged cells and explanations in
  tooltips.
- A compact strip plot for each genome-wide metric across samples.
- Native HTML `<details>` drill-down sections per sample.
- Reference- and query-contig diagnostic tables kept visually separate.
- Dotplot thumbnails, with contigs having breakpoint evidence shown first, linked to the
  full-resolution PNGs.
- Metric definitions, denominator sources, filtering parameters, and warnings about raw
  versus union-based metrics.

Relative flags use direction-aware robust z-scores based on median and MAD × 1.4826 and are
enabled only with at least five samples. Every metric also has a minimum absolute-difference
floor so statistically unusual but scientifically trivial deviations are not flagged. When
MAD is zero, do not emit a relative flag unless that metric has an explicitly configured
zero-MAD difference floor and the observation exceeds it; never divide by zero. Low
coverage/identity are bad-direction outliers; high gap fraction, overlapping bp, block count,
and breakpoints/Gb are high-direction outliers. Absolute thresholds are optional and off by
default because expected identity and fragmentation depend on divergence and data type.

## 5. Snakefile and configuration

- `rule list_maf_contigs`: one small Slurm job per sample producing the shared contig-list
  input used by the resolved manifest for both QC and the default workflow.
- `rule maf_stats`: one job per sample, input raw MAF plus `REF_FAI` and optional query `.fai`;
  outputs the sample TSVs and configured plots. Resources come from `maf_stats_threads`,
  `maf_stats_mem_mb`, and `maf_stats_time`.
- `rule maf_stats_report`: one small Slurm job aggregating all samples and creating the main
  TSV and HTML. It receives explicit time and memory rather than sharing the controller
  allocation.
- `rule maf_qc`: the QC-only user target depending on the report and all declared diagnostic
  outputs.
- Do not add QC outputs to `rule all`, `summary.html`, `split_sample_maf`, or the main calling
  DAG. The user explicitly launches the main workflow only after inspecting QC.
- Every per-sample output path is unique, and Slurm output/error names retain `%j`, satisfying
  Farm/Quobyte's one-writer-per-file requirement.
- Restrict `localrules` to trivial metadata operations such as creating a symlink. Any rule
  that scans, indexes, concatenates, reports on, or plots biological data gets explicit small
  resources and runs as its own Slurm job.

Proposed config keys:

```yaml
query_fai_dir: null
maf_stats_min_block_bp: 0
maf_stats_overlap_tolerance_bp: 0
maf_stats_dotplots: flagged       # flagged | all | false
maf_stats_dotplot_min_breakpoints: 5           # provisional; tune after validation
maf_stats_dotplot_min_breakpoints_per_gb: 10.0 # provisional; tune after validation
maf_stats_max_dotplots_per_sample: 10
maf_stats_flag_min_reference_coverage: null
maf_stats_flag_min_identity: null
maf_stats_flag_max_gap_fraction: null
maf_stats_flag_max_breakpoints_per_gb: null
maf_stats_flag_min_absolute_difference: {}    # optional per-metric floors
maf_stats_flag_zero_mad_difference: {}        # optional per-metric floors
maf_stats_threads: 1
maf_stats_mem_mb: 16000           # provisional; replace after benchmarking
maf_stats_time: "08:00:00"        # provisional; replace after benchmarking
maf_stats_report_mem_mb: 4000
maf_stats_report_time: "00:30:00"
maf_contigs_mem_mb: 2000           # provisional shared-manifest scan
maf_contigs_time: "04:00:00"       # provisional; compressed MAFs may be I/O-bound
```

Do not add a general `maf_stats: true` switch: `maf_qc` itself opts into this stage, while the
default main target remains unchanged.

## 6. Tests

Run non-trivial workflow tests through `hpc_run` in the `argprep` Conda environment.

- Parser validation: malformed row counts, unequal alignment lengths, invalid coordinates or
  strands, inconsistent `srcSize`, unknown/ambiguous reference names, empty MAF, and `.maf.gz`.
- Metric tests: identity/gap counting with N and lowercase, overlapping interval unions and
  `overlapping_reference_bp`, query `-` coordinate conversion, raw block N50, and query coverage
  with and without a query `.fai`.
- Breakpoint tests: each property on both strands, combinations where more than one property
  is true, unique adjacency counting, overlap tolerance, minimum block size, duplicated query
  intervals, and deterministic tie-breaking independent of input order where keys differ.
- Dotplot tests: forward/reverse segment endpoints, multiple query-contig y-axis bands,
  `flagged|all|false`, count/rate floors, per-sample plot cap, deterministic ranking and
  filenames, and bounded image dimensions.
- Report tests: one-row-per-sample aggregation, direction-aware outliers, minimum-sample guard,
  MAD = 0, minimum absolute-difference floors, absolute thresholds, separate reference/query
  drill-downs, escaped HTML, missing optional query FAI, and dotplot links.
- Workflow tests: `maf_qc` produces all declared outputs; the default target neither runs nor
  requires QC; both targets reuse the contig manifest; substantive jobs have explicit resources
  under the Slurm profile; and the controller wrapper forwards the `maf_qc` target.

## 7. Documentation and real-data validation

Document the two-command workflow, metrics, output files, configuration, query-coverage
denominator caveat, dotplot interpretation, and the fact that QC does not automatically gate
or start the main workflow. Update `README.md`, `example_data/options.yaml`, and `changelog.md`.

Before selecting defaults or treating flags as scientifically meaningful, run `maf_qc` on a
small representative set of real project MAFs through Farm Slurm and check:

- Coverage against independently merged reference intervals for selected chromosomes.
- Identity and gap counts in manually inspected regions.
- Whether overlapping/secondary blocks inflate raw metrics as expected.
- High-breakpoint samples and contigs against their dotplots.
- Query coverage with and without query `.fai` files.
- Whether cohort outlier flags identify plausible problems rather than normal biological
  divergence.
- Runtime, peak RSS, output size, and PNG count, then set conservative resource defaults.

Real-data validation reads inputs and writes new QC outputs only; it does not modify source
MAFs.

## Remaining decisions before implementation

1. Confirm whether per-sample query `.fai` files are normally available and how filenames map
   to sample IDs.
2. Benchmark bounded `flagged` dotplots and choose meaningful minimum count/rate thresholds;
   retain the provisional cap of 10 unless output volume supports another value.
3. Decide whether a future version should optionally emit multipage PDFs in addition to PNGs.
4. Replace provisional `maf_stats_mem_mb` and `maf_stats_time` using real MaxRSS and elapsed
   measurements.
