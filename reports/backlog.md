# ARGprep backlog

The living list of known-open work. Previously this lived at the bottom of
[code_review/2026-07.md](code_review/2026-07.md); that file and
[code_review/2026-08-09.md](code_review/2026-08-09.md) are historical records of
what those reviews found and should not be edited further. Track open work here.

Last reconciled against the code at v1.9. Items marked *(verified open)* were
re-checked in the source at that point, not carried over on faith.

## Features

- *(Nothing open.)* The "real alignment FASTA from the MAFs" item shipped as
  [scripts/maf_to_fasta.py](../scripts/maf_to_fasta.py): a second consumer of the
  core calling logic that keeps indels, with insertions left-anchored and
  left-aligned, conflicting blocks degraded to `N`, VCF agreement behind
  `--mask-bed`, and a `--max-columns` guard against width blowup. The homology
  caveat is inherent to a reference-anchored merge and is documented rather than
  solved. `scripts/window_to_fasta.py` keeps the cheap reference-anchored view.

## Bugs and hardening

- **Low — `report_stats.tsv` has no record of the `summary_window_bp` it was
  written with** *(verified open, new at v1.9)*. `read_report_stats` computes
  `idx = int(start) // window_bp` and *assigns* into the window list. Running
  `summary_report.py` by hand with a different `--window-bp` than the stats file
  used silently collapses multiple windows onto one index instead of failing.
  The workflow always passes `summary_window_bp` to both sides, so this only
  bites manual invocation. Fix: record the window size in the file and validate
  it on read.
- **Feature question — should soft-masked repeats be excludable?** The behavior
  (upper-case everything, so lowercase repeat bases are genotyped normally;
  `maf_to_sites.py:398` and `read_contig_sequence`) is documented in the README
  as of v1.9, which closes the "confirm intent" half. What remains is a product
  decision: whether to add a config flag that treats soft-masked reference bases
  as missing. Today the workarounds are a hard-masked reference or repeat
  intervals via `quality_bed_dir`.

  Narrowed 2026-08-09: only the **reference FASTA** path is live. Reference bases
  come from `read_contig_sequence` (`maf_to_sites.py:555,662`), not from the MAF,
  so a soft-masked `reference_fasta` has its repeat content genotyped regardless
  of the aligner. The query path is moot in practice — scanning a real AnchorWave
  MAF found no lowercase in the sequence rows, so case is already lost upstream
  and the `.upper()` on sample rows is a no-op. (Sampled, not exhaustive, and a
  query genome that was never soft-masked would look the same.)

  Two follow-ups, both optional:
  - The README note says "reference and query sequences are upper-cased", which
    is true but implies a query-side concern that does not arise for AnchorWave
    input. Consider narrowing it to the reference, or keep it as documentation of
    what happens with MAFs from other aligners. Deliberately left as-is for v1.9.
  - Note that `.upper()` is load-bearing, not cosmetic: without it a lowercase
    base fails the `VALID_BASES` check and is assigned `?`, i.e. **missing**. Any
    "exclude soft-masked repeats" flag should mask them explicitly rather than
    reach that behavior by dropping the upper-casing, so the two causes stay
    distinguishable in the summary counters.
- **Low — split normalization collision aborts with no override** *(verified
  open)*. `split_maf_by_contig.py:47` raises on an ambiguous normalized contig
  with no way to force a mapping.

## Tests (coverage gaps)

All *(verified open)* — no test currently exercises these:

- Contig-end overrun (a MAF ref row extending past the `.fai` length).
- Empty / header-only MAF (expect all `masked_no_alignment`).
- The coverage-invariant `ValueError` path in `main`.
- `read_contig_sequence` FAI-mismatch and gzip-fallback branches.

**Config propagation is now fully covered** (closed in v1.9). Every config key
that changes which sites are retained has a both-states end-to-end workflow test
proving the YAML value reaches the caller — the bug class that hid the
`allow_multiallelic_snps` defect behind a green suite: `allow_multiallelic_snps`,
`max_missing_fraction`, `mask_indel_adjacent_snps`, `add_ref`, and the
`quality_bed_dir` / `quality_min` pair. All were verified correct. The tests were
mutation-checked by breaking each `cmd+=(...)` line in the Snakefile and
confirming the enabled-state case fails.

## Docs

- **Low — changelog v1.0 lists filenames that never shipped**
  (`combined.<contig>.variants.vcf`, `.masked.bed`; actual names are `.vcf` and
  `.mask.bed`). Deferred by maintainer decision; recorded so it is not
  rediscovered as a bug.

## Deferred by decision

- **`temp()` on the per-contig MAF chunks** *(verified open)*. The gzip half of
  2026-07 efficiency #6 shipped in v1.8; marking
  `results/maf_by_contig/<sample>` as `temp()` is still undone and is DAG-safe,
  but `test_workflow_remaps_requested_contigs_to_reference_names` asserts the
  chunks persist, so it needs a test update. Maintainer decision pending.
- **2026-07 efficiency #5** — str→bytes reference conversion. Intentionally
  skipped; low value once the classification loop was vectorized.
- **`summary_report.py` remains a single ~810-line module** doing parsing,
  aggregation, HTML, CSS, and SVG. v1.9 removed its VCF parser and mask
  re-binning, which was the bulk of the 2026-08-09 recommendation; splitting the
  renderer further has not been judged worth it.

## Closed since the 2026-07 review

Recorded so these are not re-opened from the old list:

- Zero-length contig handling — v1.7 added an empty-output path; v1.9 replaced it
  with an explicit rejection.
- Leading-insertion indel-adjacent flagging — fixed in v1.8.
- `intervals_from_positions`, `summarize_site_and_mask_coverage` — deleted.
- "Unaligned" genome-overview segment always 0 — the segment is gone.
- gzip per-contig MAF chunks (half of efficiency #6) — shipped in v1.8.
- Test gaps for `parse_maf_path_map` error branches, split `.maf.gz` input, and
  `read_sample_missing_bp` longest-suffix matching — all now covered.
- Hardcoded personal paths in the two ad-hoc analysis scripts — generalized to
  argparse; `chr1_variable_plot.py` was deleted outright in v1.9 once the report's
  y-axis scaling made it unnecessary.
