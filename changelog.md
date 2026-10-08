# Changelog

Versions are git tags. For production runs, check out the most recent **stable**
release (`git checkout v1.10`). Tags with a `-beta.N` suffix are pre-releases
for testing new features: `v1.11-beta.1` adds the `maf_qc` alignment-QC stage,
and feedback on it is welcome.
v1.0 was a full rewrite from the legacy TASSEL/gVCF/GATK pipeline — the
pre-v1.0 (`v0.x`) entries at the bottom describe that older lineage and do not
carry forward to the v1.x series.

## v1.11-beta.1 (pre-release)

The new `maf_qc` stage is in **beta**: metric definitions, flag thresholds, the
report layout and the output columns may change after user feedback. The rest of
this release is stable.

### First-stage MAF QC: `maf_qc` target ([#33](https://github.com/RILAB/argprep/issues/33)) — beta

- New QC-only target, `snakemake ... maf_qc`. It runs [scripts/maf_stats.py](scripts/maf_stats.py) as one SLURM job per sample, then [scripts/maf_stats_report.py](scripts/maf_stats_report.py) as one small report job. Outputs go to `results/maf_stats/`.
- It is not part of `rule all`. The main workflow neither runs it nor depends on it, so inspecting QC and starting the main run are separate decisions.
- Strict validation of pairwise MAF structure, coordinates, source sizes and reference names.
- Cross-sample flags: direction-aware robust z-scores with per-metric scale floors (needs at least 5 samples), plus optional absolute thresholds.
- [scripts/maf_dotplot.py](scripts/maf_dotplot.py) is refactored and now used by QC:
  - reads `.maf.gz`;
  - pairs rows by block through the shared validated parser;
  - stacks query contigs in labelled bands;
  - writes fixed-size PNGs with batched line drawing;
  - exits non-zero when an input fails.
- `MafRecord` and `iter_maf_blocks` moved from `maf_to_sites.py` to `common.py` so there is one MAF parser; `maf_to_sites` re-exports both.
- Metrics, with alignment and assembly quality reported separately:
  - Coverage headlines are base-aligned (`aligned_reference_pct`, `aligned_query_pct`). Block-span coverage is kept as a secondary metric, because AnchorWave blocks span chromosomes between anchors: in a test on two teosinte genomes, blocks spanned about 94% of B73 v5 but about 36% of reference bases were aligned.
  - Identity over A/C/G/T pairs; insertion and deletion columns reported separately (a single gap fraction was dominated by AnchorWave's long within-block indels).
  - Strand flips, reference-contig jumps and out-of-order adjacencies counted independently, plus unique breakpoint adjacencies per Gb aligned.
  - Optional query assembly FASTA (`--fasta`, `query_fasta_dir`) adds scaffold and contig N50 (contigs split at runs of at least 10 N), N gaps, GC, soft-masking (`NA` for unmasked FASTAs), telomeric ends on chromosome-scale sequences, and the content of unaligned sequences. Without it, MAF-only contiguity (aligned contigs only) is reported. Query names such as `<assembly file>.chr6` resolve to FASTA name `chr6`.
  - Nested blocks (a query interval inside a larger block's) are reported separately and kept out of the breakpoint walk, where each would have produced two fake strand flips. Also `overlapping_query_bp`.
  - Each breakpoint is located on its query contig (`contig_end`, `n_gap` with the FASTA, or `interior`) and checked for recurrence across samples (`maf_stats_recurrence_window_bp`). New outputs: `<sample>.breakpoints.tsv` and `maf_stats.breakpoints.tsv`.
  - Block comparisons are chunked, to bound memory on chromosome-scale blocks.
- `matplotlib` added to `argprep.yml` (and so to the containers).
- Report leads with a short summary:
  - Per-sample status: ✔ normal, ⚠ review (bad-direction outlier or absolute threshold), ◆ unusual (large deviation in the good direction, e.g. far closer to the reference than the rest; check for admixture, contamination or mislabelling).
  - Plain-language reasons that point to contigs, coordinates and dotplots.
  - A compact sample table and a shared-features table (breakpoints and nested blocks recurring across samples).
  - Everything else in collapsed sections.
  - `maf_stats.tsv` gains `status` and `unusual` columns.
  - `overlapping_reference_bp` and `reference_contig_jumps` are reported but not flagged (0 for every sample under AnchorWave's one-to-one alignment); unaligned sequence is flagged when high.
- The QC rules pass a hash of their scripts as a param, so updating ARGprep reruns QC instead of keeping stale outputs.
- Tested on 8 NAM maize lines (with assembly FASTAs) and 8 teosintes (MAF only) against B73 v5.

### Other changes

- New output `sites/combined.<contig>.sample.mask.bed`: the per-sample missing BEDs for each contig concatenated into one file, labelled by sample ([#32](https://github.com/RILAB/argprep/issues/32)). It is grouped by sample, not coordinate-sorted, and is not a merged mask. Follow-ups are tracked in [#34](https://github.com/RILAB/argprep/issues/34).
- MAF blocks whose reference row is on the minus strand are now rejected with an error. `maf_to_sites.py` and `maf_to_fasta.py` would previously have placed them at the wrong reference coordinates without warning; AnchorWave never writes them ([#35](https://github.com/RILAB/argprep/issues/35)).
- `prepare_reference` and `index_reference` now request explicit small SLURM resources instead of the profile's 20 GB / 48 h default.
- Open work moved from `reports/backlog.md` to [GitHub issues](https://github.com/RILAB/argprep/issues); see `reports/code_review/2026-10.md` for the latest review.

## v1.10

### New auxiliary script: gapped alignment FASTA from the MAFs

- Added [scripts/maf_to_fasta.py](scripts/maf_to_fasta.py), which builds a true
  gapped multiple alignment for a reference window by merging the per-sample
  pairwise MAFs through the reference. Deletions render as gap columns and
  insertions add columns the other samples pad with gaps, so it is not limited to
  the substitution-only view `window_to_fasta.py` reconstructs from the
  `all_sites.vcf`. It reads the per-contig MAF chunks under
  `results/maf_by_contig/` (or whole-genome MAFs via `--maf-dir`) and reuses the
  existing MAF parser, block selection, conflict resolution, and quality-mask
  machinery, so its calls agree with the VCF path base for base.
- Like `window_to_fasta.py` it is standalone: no workflow rule invokes it, and the
  pipeline's deliverables are unchanged. The window is user-chosen, so there is
  nothing for `rule all` to target.
- Conventions, all documented in the README: insertions are left-anchored to the
  preceding reference base and left-aligned within its columns (which asserts a
  homology the pairwise data cannot establish — inherent to any reference-anchored
  merge); conflicting alignment blocks degrade an anchor to `N` and an insertion to
  a run of `N`; `--mask-bed` reconciles the output with `mask.bed` rather than
  shipping two contradictory views by default; `--max-columns` fails loudly instead
  of exhausting memory in insertion-rich regions.
- `maf_to_sites.py` gained `read_contig_length` and `read_contig_region`, which use
  the `.fai` index to seek straight to a window so extracting a kilobase region does
  not cost a whole-contig read. Existing callers are unaffected.
- Closes the only open item under "Features" in [reports/backlog.md](reports/backlog.md).

## v1.9

Acts on the [2026-08-09 code review](reports/code_review/2026-08-09.md). The
primary scientific outputs (`all_sites.vcf`, `combined.<contig>.vcf`,
`mask.bed`, `site_summary.tsv`, `*.missing.bed`) keep the same formats and
semantics; the changes are one correctness fix, a stricter configuration and
input-validation contract, and a reporting-path simplification.

### Behavior changes (review your config when upgrading)

- Fixed `allow_multiallelic_snps: false`, which previously never passed the CLI's masking flag and therefore still retained multiallelic sites. **Runs that set this to `false` did not get what they asked for and will now produce fewer retained sites.** `true` (the default) is unaffected.
- Explicit `contigs` entries now all have to resolve to the reference, exactly or through one unambiguous normalized alias. A partially invalid explicit list now fails and names the unresolved entries, instead of silently running a subset. Automatic contig discovery retains its skip-and-warn behavior.
- Contig-name normalization now strips an assembly/genome prefix when it precedes a `chr` token, so `chr2`, `Zm-B73v5.chr2`, and `Zx-TIL25.chr2` all resolve to contig `2`. The rule is gated on the `chr` token, so accession-style names such as `NC_050096.1` (where `.1` is a version suffix) are left intact. This can make previously unmatched contigs resolve, and MAFs whose reference contigs differ only by assembly prefix now intersect.
- Sample auto-discovery is now non-recursive: `*.maf` / `*.maf.gz` are matched only directly in `maf_dir`, not in nested subdirectories. Previously the `{sample}` wildcard could span `/`, so pointing `maf_dir` at a directory that also contained an `example_data/example.maf/` tree (or any nested MAFs) pulled those in as bogus samples, whose contigs failed to intersect and triggered "No contigs are shared across all MAF files."
- Removed the unused `maf_threads` setting; site calling is single-threaded and each per-contig job requests one core. The key is now ignored if left in an existing `options.yaml`; drop it and keep tuning `maf_mem_mb` / `maf_time`.

### Fail fast on malformed or degenerate input

Previously these cases were silently absorbed, turning bad input into missing
data, skipped contigs, or partial summaries. They now raise with the file and
line number.

- Zero-length reference contigs are rejected rather than taking a duplicated empty-output path (`write_empty_contig_outputs` and its tests are gone). A 1 bp contig still works through the normal path.
- MAF `s` rows with fewer than 7 fields, or non-integer `start`/`size`/`srcSize`, now raise. Alignment blocks whose sequence strings have unequal lengths now raise instead of being silently truncated by `zip`.
- Quality-BED rows with fewer than 4 fields, or non-numeric coordinates/score, now raise. Comment, `track`, and `browser` lines are still skipped.
- Documented that **soft masking is ignored**: reference and query sequences are upper-cased before calling, so lowercase repeat-masked bases are genotyped like any other. A soft-masked reference — the default distribution format for most released genomes — has its repeat content called, not excluded. Hard masking still works, since `N` is non-ACGT and counts as missing. No behavior change; this was previously undocumented, and it is a scientific choice rather than an implementation detail.
- `read_maf_contigs` no longer swallows `OSError`; an unreadable or missing MAF fails during contig discovery rather than contributing an empty contig set.

### Reporting path

- `all_sites.vcf` remains a required scientific output — it is *not* being made optional or removed. Instead, site calling now accumulates report counters during the pass it was already making and writes `sites/combined.<contig>.report_stats.tsv`: one `window` row per `summary_window_bp` window (invariant / variant / masked counts) and one `sample` row per sample (called and variant-carrying retained sites).
- `summary_report.py` reads those counters instead of reparsing every genotype cell of every `all_sites.vcf` and re-binning every `mask.bed` interval. On chromosome-scale data this removes a second full pass over the largest artifact. Reported values are unchanged.
- `maf_to_sites.py` gained `--window-bp` (default 100000); the workflow passes `summary_window_bp` so the counters and the report always agree on window boundaries.
- `summary_report.py`'s `--all-sites` and `--masked-beds` arguments were replaced by `--report-stats`. Anyone invoking the script outside the workflow must update their command.
- Report-stats parsing is strict: an unexpected header, column count, or record type raises instead of being skipped.

### Optional ARGweaver `.sites` output

- Disabling `emit_argweaver_sites` no longer deletes a stale `.sites` file as a side effect of rebuilding `summary.html`. Cleanup was unrelated to report generation; use a clean results directory when changing output modes if stale optional files matter.
- Consequently `.sites` files are no longer inputs to the `summary_report` rule. Ordinary runs are unaffected: `.sites` is a declared output of `direct_maf_sites`, which has to run anyway, so any fresh run produces it — whether you use the default `all` target or ask only for `results/summary.html`. The one case that changed is enabling `emit_argweaver_sites` in a results directory whose other outputs are already up to date. The default target still produces `.sites` there, because `rule all` lists it; asking only for `results/summary.html` does not, because nothing in that sub-DAG forces `direct_maf_sites` to re-run. Use the default target when toggling the flag in place.

### Summary report plot scaling

- Split the per-contig plot in two. Invariant and missing are complementary shares of a window and keep the fixed 0–100% axis; variable sites, typically a fraction of a percent, were an unreadable flat line against that axis and now get their own plot.
- The variable-site axis is auto-scaled but **shared across all contigs**, rounded up to a readable 1/2/2.5/5 × 10ᵏ bound. Scaling each contig independently would have made the y-axis mean something different in each section and silently broken cross-contig comparison — the reason the axis was pinned to 0–100% originally.
- Axis ticks now carry enough decimals to be exact. A quartered axis produces steps like 1.25, which the previous integer formatting rendered as `0, 1, 2, 4, 5`.
- This removes the need for the deleted `chr1_variable_plot.py` (below), which existed only to work around the flattening.

### Removed dead code

None of these had callers in the workflow, but they were importable, so remove them from any external scripts:

- `scripts/common.py`: `read_fasta_contigs`, `extract_info_int`, `extract_end`, `extract_dp`.
- `scripts/summary_report.py`: `read_mask_percentages`, and the never-consumed `per_sample_missing_retained_by_contig` aggregation.
- `scripts/simulate_msprime_indels.py`: `apply_indel_events` (exercised only by its own unit test, not by the simulation CLI).
- Deleted the `scripts/chr1_variable_plot.py` helper. It re-parsed SVG polyline coordinates back out of a generated `summary.html` to redraw the chromosome-1 variable-site line on an auto-scaled y-axis, working around the combined plot's fixed 0–100% scale. Reading a rendered chart back in as a data source was brittle and assumed chromosome-specific element IDs and chart geometry; the fix belongs in the report's y-axis scaling, not in a downstream re-render.
- Kept `scripts/window_to_fasta.py`; retaining `all_sites.vcf` as a required output keeps its design valid, since reconstructing sequence from the reference is only sound when invariant positions were positively called. The README section is now titled "Auxiliary scripts" and states explicitly that no workflow rule invokes them. Its entry now documents that output is a reference-anchored substitution view, not a true alignment: insertions consume no reference coordinate and are dropped, deletions are recorded as missing and render as `N`, and an `N` therefore conflates unaligned, deleted, ambiguous, and below-`quality_min` bases.

### Tests

- 109 tests pass (was 77). Every config key that changes which sites are retained now has a both-states end-to-end workflow test proving the YAML value actually reaches the caller — `allow_multiallelic_snps`, `max_missing_fraction`, `mask_indel_adjacent_snps`, `add_ref`, and `quality_bed_dir`/`quality_min`. Only the first was broken; the rest were verified correct. This closes the gap that let the multiallelic defect sit behind a green suite, since the CLI flags themselves had unit tests all along. Also new: explicit-`contigs` hard failure vs. discovery skip-and-warn, non-recursive sample discovery, assembly-prefix normalization, strict MAF/BED/report-stats parsing, the report built from `report_stats.tsv`, and axis-bound rounding, tick precision, and out-of-range clamping for the plots. Tests for the deleted zero-length-contig writer and `apply_indel_events` were removed with the code.

## v1.8

- Added optional ARGweaver `.sites` output. Set `emit_argweaver_sites: true` to also emit `sites/combined.<contig>.sites` (variant sites only; one real base per pseudo-haploid sample, `N` for missing) alongside the VCFs. Disabled by default, so existing runs are unchanged. When `add_ref` is set, the synthetic REF haplotype is appended as the final column. The variant-site set matches `sites/combined.<contig>.vcf`.
- `split_sample_maf` now writes per-contig MAF chunks gzip-compressed (`results/maf_by_contig/<sample>/<contig>.maf.gz`) instead of uncompressed, so `.maf.gz` inputs no longer get re-materialized as a full uncompressed copy of the alignment corpus under `results/`.
- Fixed indel-adjacent SNP masking (`mask_indel_adjacent_snps`) for insertions at the start of an alignment block: a SNP immediately following such a leading insertion was not being flagged as indel-adjacent, so it slipped through unmasked. It is now flagged like any other post-insertion site.
- `summary_report` now treats the config file as a tracked input (so the summary rebuilds when it changes) and removes stale `sites/combined.<contig>.sites` files left over from a previous run when `emit_argweaver_sites` is disabled.
- Documented preparing inputs for downstream ARG inference (Relate via `RelateFileFormats`, native ARGweaver `.sites`, and SINGER's direct VCF input) in the README.

## v1.7

- Added a `split_sample_maf` rule (`scripts/split_maf_by_contig.py`) that partitions each per-sample pairwise MAF into per-reference-contig chunks under `results/maf_by_contig/<sample>/<contig>.maf`. These are regenerable intermediates.
- `direct_maf_sites` now consumes the per-contig chunks via a new `--maf-paths SAMPLE=PATH` argument to `maf_to_sites.py`, so each contig reads only its own slice of every MAF instead of rescanning each full MAF once per contig.
- Reduced peak memory in `maf_to_sites.py`: per-sample call arrays are now mmap-backed temp files and mask/missing intervals are streamed rather than materialized in memory. Added an explicit coverage invariant (`retained_total + masked_total == contig_len`) to the internal contig-span check.
- Faster site calling: the main all-sites classification, the per-sample missing-BED pass, and the indel-adjacent flag merge are now vectorized with NumPy over the mmap-backed call arrays, and indel-adjacent tracking is skipped entirely unless `--mask-indel-adjacent-snps` is set.
- Fixed a crash on zero-length reference contigs (empty FASTA records): `maf_to_sites.py` now emits header-only VCFs and empty BED/summary outputs instead of failing to `mmap` an empty buffer.

## v1.6

- Added optional per-sample assembly-quality masking. Set `quality_bed_dir` (a directory of `<sample>.bed`/`.bed.gz` files in each sample's own genome coordinates, rows `chrom start end score`) and `quality_min` to treat aligned bases scoring below the threshold as missing. Disabled by default; masked bases flow through the existing `max_missing_*` thresholds and per-sample `*.missing.bed` outputs. Both `+` and `-` strand query alignments are supported.

## v1.5

No changes to pipeline outputs since v1.4. SLURM preemption safety only:

- Added a cluster-generic status command (`profiles/slurm/status-sacct.sh`) that maps `PREEMPTED` to *running* and genuine failures (`TIMEOUT`, `OUT_OF_MEMORY`, `NODE_FAIL`, `scancel`) to *failed*. Rule jobs can now run on a preemptable queue (set `slurm_partition: low`): a preempted job is auto-requeued by SLURM and the controller waits for the rerun, while real failures surface cleanly instead of hanging in "wait for output" mode.
- Added `profiles/slurm/run-controller.sbatch` to run the long-lived Snakemake controller on a non-preemptable partition so it outlives the rule jobs.
- Documented the recommended `sbatch` run and how to override the two cluster-specific values (`--partition`/`--account`) for other clusters.

## v1.4

Behavior changes (may change variant counts on the same `options.yaml` — review your config when upgrading):

- Changed the default for `mask_indel_adjacent_snps` from `true` to `false`. Set it explicitly to `true` in `options.yaml` (or pass `--mask-indel-adjacent-snps`) to opt back in.
- Removed the `mask_indels` config option. Deletion gaps (`-`) already count toward per-sample missingness, so setting `max_missing_count: 0` reproduces the old `mask_indels: true` behavior.
- Renamed the `masked_indel` site-summary counter to `masked_indel_adjacent` (it now reflects only the adjacent-SNP masking it still controls).

## v1.3

- Redesigned the HTML summary report for clarity and browser performance.
- Removed the `treat_n_as_missing` flag; `N` is now always treated as missing.
- Per-individual missing BED files now include the sample name in the filename.
- Improved contig-resolution UX with a lazy MAF scan; default contig discovery normalizes contig names across per-sample MAFs (e.g. `chr1` and `1` intersect), and `maf_to_sites` preserves explicit `--samples` order instead of sorting.
- Clarified `summary_window_bp` scope and documented skipped-contig behavior.
- Split README/NOTES, added a citation note, and fixed citation author formatting.

## v1.2

- Added per-sample missing-data BED masks.
- Documented new config keys, TSV outputs, and defaults; fixed the `treat_n_as_missing` default and updated the README to match current pipeline behavior.

## v1.1

- Added the `add_ref` option so final `all_sites` and `variants` VCFs can include a synthetic `REF` sample with genotype `0` at every retained site.
- Skip empty scaffolds/chromosomes in `summary.html` plots.
- Embedded the config file and pipeline version in the summary report.
- Fixed bugs and inefficiencies identified in a pipeline review.

## v1.0 — new direct-MAF pipeline

- Replaced the legacy TASSEL/gVCF/GATK workflow with a direct MAF-to-sites pipeline.
- The workflow now reads per-sample MAFs directly and emits three per-contig outputs:
  - `combined.<contig>.all_sites.vcf`
  - `combined.<contig>.variants.vcf`
  - `combined.<contig>.masked.bed`
- Site classification is now done in reference coordinates from the MAF alignment itself rather than from merged gVCF records.
- Missingness filtering is configurable with either `max_missing_count` or `max_missing_fraction`; if both are set, the stricter threshold is used.
- Indel handling is configurable: `mask_indel_adjacent_snps` optionally masks SNPs adjacent to insertions or deletions.
- Multiallelic SNPs are retained by default and can be disabled with `allow_multiallelic_snps`.
- Coverage validation was simplified: `maf_to_sites.py` performs the retained-sites-plus-mask contig-span check internally without a separate `coverage.txt` artifact.
- Added a direct-pipeline `summary.html` report with 100 kb window plots showing percent invariant, percent variable, and percent missing along each contig.
- Added an msprime-based simulation helper that generates reference FASTA, sample FASTA sequences, pairwise MAFs against the reference, realized indel tables, and summary truth counts for SNPs and indel-affected reference positions.

---

# Legacy (pre-v1.0 TASSEL/gVCF pipeline)

These entries describe the older pipeline that v1.0 replaced. They are retained
for historical reference and do not apply to the v1.x series above.

## Changes since v0.4

- Clean VCF handling was tightened: `<NON_REF>` is now always stripped from `ALT`, and reference-donor `AD/PL/DP` values are copied where needed for cleaner downstream outputs.
- Summary warning collection was fixed to use only current-run logs passed by Snakemake, avoiding leakage of historical warnings from old `.snakemake` and SLURM logs.
- Pipeline defaults were updated in `config.yaml`/`Snakefile`, including example-data-oriented path defaults and associated resource/default cleanup.
- Integration tests were adjusted to be more portable (not relying on hardcoded `PATH` behavior), and test coverage was updated across split/integration/summary-related tests.
- Random-position integration test paths were aligned with example-output conventions and related test defaults were refreshed.
- README/changelog documentation was refreshed with expanded SLURM usage notes, formatting fixes, and changelog updates.

## Changes since v0.3

- Added/expanded contig handling across the workflow: reference renaming from gVCF/MAF contigs, numeric-suffix normalization, configured-contig remapping (e.g. `chr1` -> `1`), filtering contigs to the reference `.fai`, and defaulting to the shared MAF-contig intersection when `contigs:` is not set.
- Improved missing-contig behavior and reporting: split/merge steps now skip missing per-sample contigs more safely, preserve valid placeholders, and surface clearer warnings in the HTML summary.
- Expanded HTML summary reporting with stronger warning parsing, MAF-vs-reference contig mismatch warnings, ploidy reporting, and a missing-genotype exclusion histogram.
- Fixed summary warning de-duplication so contigs reported in "Configured contigs were remapped..." are not also repeated in generic MAF/reference mismatch warnings.
- Improved split/filter/mask robustness: refined split classification and filtered BED span handling, gzip/bgzip split-output support, gzipped missing-BED input support, and better empty-coverage handling in `check_split_coverage.py`.
- `.clean.vcf` output now always strips `<NON_REF>` from `ALT`, including GT-only (`outputJustGT=true`) runs.
- Added optional `add_reference` split output behavior: when enabled, `.clean.vcf` appends a synthetic `REF` sample with reference-only genotypes (`0`, `0/0`, etc.) and fallback ploidy matching the pipeline-resolved ploidy.
- Fixed boolean config parsing in `Snakefile` so CLI overrides like `--config bgzip_output=false add_reference=true` are interpreted correctly instead of treating non-empty strings as truthy.
- Added accessibility mask outputs (`combined.<contig>.accessible.npz`) for scikit-allel workflows and documented downstream use.
- Added ploidy inference from MAF blocks (with summary reporting) while retaining config override support.
- Added and refined resource/concurrency controls: configurable `maf_to_gvcf_*` and `merge_contig_*` resources, `merge_gvcf_max_jobs`, `maf_to_gvcf_array_max_jobs`, GenomicsDB buffer tuning, and `default_mem_mb` fallback behavior.
- Improved SLURM profile/runtime behavior: numeric default resource values in the profile, SLURM stdout/stderr written to `logs/slurm/`, and clearer documentation of profile `jobs` (Snakemake parallelism) vs per-job CPUs.
- Added substantial test coverage: pytest scaffolding, split/filt/dropSV unit tests, coverage checks, integration tests (including contig mismatch/default-contig behavior), random-position validation, summary-report regressions, Snakefile contig-resolution unit tests, and a SLURM-profile dry-run smoke test.
- Documentation updates: `logic.md`, expanded README notes/options/outputs/testing instructions, and a standalone `changelog.md`.

## Changes since v0.2

- Moved HTML summary generation into `scripts/summary_report.py` and simplified `Snakefile`.
- Corrected example MAF inputs so `example_data/*.maf.gz` are valid gzip files.
- Updated split classification so `ALT=<NON_REF>`-only records are treated as invariant.
- Updated SINGER clean-output formatting to strip `<NON_REF>` while preserving genotype/sample fields.
- Added `logic.md` with detailed site-routing/filtering logic and concrete examples.
- Added `results/split/combined.<contig>.coverage.txt` to documented workflow outputs.
- Added split-test coverage for invariant/nonref and genotype-preserving clean formatting.
- Updated SLURM profile default resources to numeric values to avoid resource conversion/submission errors.
- Added `merge_gvcf_max_jobs` pipeline concurrency control.

## Changes since v0.1

- Added HTML summary report with embedded SVG histograms and expanded output details.
- Split logic tightened: clean sites now require all samples called; missing GTs are routed to filtered.
- Invariant/filtered/clean outputs are enforced as mutually exclusive per position; filtered BED spans now respect END/REF lengths and subtract inv/clean.
- Merged gVCFs are produced via GATK SelectVariants with genotype calling; TASSEL `outputJustGT` default set to `false` to retain likelihoods for calling.
- Added accessibility mask generation (`combined.<contig>.accessible.npz`) for scikit-allel workflows.
- New/expanded validation and tests: split coverage checks, filtered-bed tests, integration tests gated by `RUN_INTEGRATION=1`.
- Example data regenerated via msprime with indels and missing data and AnchorWave-style MAF formatting.
- `check_split_coverage.py` now reports overlap intervals with file names to aid debugging.
- `filt_to_bed.py` filters masks to the target contig, preventing cross-contig lines in `combined.<contig>.filtered.bed`.
- SLURM default resources now read `default_*` from `config.yaml` instead of hardcoded profile values.
