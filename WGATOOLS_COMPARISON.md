# MAF→VCF: ARGprep vs. wgatools

A summary of the main algorithmic differences between how ARGprep
([scripts/maf_to_sites.py](scripts/maf_to_sites.py)) and
[wjwei-handsome/wgatools](https://github.com/wjwei-handsome/wgatools)
(`wgatools call`, [src/tools/caller.rs](https://github.com/wjwei-handsome/wgatools/blob/main/src/tools/caller.rs))
generate VCFs from MAF alignments.

The two tools solve different problems with the same input format, so the
"differences" below are mostly differences in goal that drive different
algorithms, rather than competing implementations of the same algorithm.

## TL;DR

- **ARGprep** projects *N* per-sample pairwise MAFs onto a shared reference
  coordinate grid and emits an **all-sites multi-sample VCF** (plus a
  variant-only VCF and a BED mask) suitable for ARG estimation. It tracks
  missingness explicitly and masks indel-affected positions rather than
  calling indels.
- **wgatools `call`** walks a **single pairwise MAF** (one target vs. one
  query) column-by-column and emits a **variant-only single-sample VCF**
  containing SNPs, large insertions, large deletions, and inversions. There
  is no all-sites output and no missingness model.

## Side-by-side

| Aspect | ARGprep `maf_to_sites.py` | wgatools `call` (caller.rs) |
|---|---|---|
| Input cardinality | N per-sample MAFs against a common reference | One pairwise MAF (or PAF) |
| Output samples | Multi-sample VCF (one column per sample) | Single-sample VCF (one query column) |
| Coordinate frame | Reference-anchored, per-contig call arrays of length `contig_len` | Walk alignment columns, maintaining target/query offsets |
| Sites emitted | Every reference base is either retained (in VCF) or masked (in BED); the union spans the contig exactly | Only positions inside aligned blocks; unaligned regions are silently absent |
| Variant types | SNPs only (multi-allelic optional). Indels are a masking signal, not a variant. | SNPs + INS + DEL (length > `--svlen`, default 50) + INV (`<INV>` symbolic) |
| Missingness | Explicit per-site count; site rejected if `missing > threshold`. Missing = unaligned, `-`, `N`, or non-ACGT. | No concept — only aligned bases are considered |
| Indel encoding | Indel-overlapped (and optionally adjacent) reference positions go to BED mask | VCF records with anchor-base REF/ALT encoding plus `SVTYPE/SVLEN/END` |
| Strand / inversion | Not handled at the calling layer — bases are projected as AnchorWave aligned them | Negative-strand blocks emit a leading `<INV>` record and add `INV_NEST=TRUE` to nested variants |
| Chunking | Snakemake first splits each MAF by reference contig; site classification then processes 1 Mb reference-position chunks over mmap-backed per-sample arrays | Default 1 Mb alignment-column chunks with an "SV-safe boundary" search |
| Genotype | Haploid integers (`0`, `1`, ..., `.` for missing) | Diploid phased (`1\|1`) plus a custom `QI` FORMAT field encoding query coords + strand |
| VCF builder | Hand-written `##` header and tab-joined record lines | `noodles::vcf` typed builder |
| QC outputs | Per-contig `site_summary.tsv` and `report_stats.tsv`, per-sample missing BEDs (plus a per-contig `sample.mask.bed`), BED mask, genome-wide `summary.html`; optional `maf_qc` per-sample MAF QC report | None beyond the VCF |

## How the core loop differs

### ARGprep — reference-anchored, per-position scan

For each sample, [load_sample_calls](scripts/maf_to_sites.py) walks the sample's
per-contig MAF and writes integer base codes into a `bytearray` indexed by
reference position. Reference-aligned `-` columns advance neither reference
index and flag neighboring reference positions for optional indel-adjacent
masking; gap characters in the sample are recorded as missing.

The main loop then scans the contig in 1 Mb chunks by stacking zero-copy NumPy
views over the mmap-backed sample call arrays. At every reference base it:

1. Reject if `ref_base` is non-ACGT.
2. Counts missing samples from the chunked call matrix.
3. If indel-overlapped (or adjacent and `mask_indel_adjacent_snps`), mask.
4. If `missing > allowed_missing`, mask (distinguishing unaligned vs. missing).
5. Otherwise emit either an invariant record (`ALT=.`, `SC=invariant`) into
   the all-sites VCF, or a variant record into both VCFs.

Because the loop visits every reference position exactly once, retained
sites + BED mask provably tile the contig. Production code checks this directly
with the inline invariant `retained_total + masked_total == contig_len`.

### wgatools — alignment-column walk with run-length categories

[`call_within_var`](https://github.com/wjwei-handsome/wgatools/blob/main/src/tools/caller.rs#L388)
zips the target and query alignment strings and groups consecutive columns
by [`cigar_cat_ext_caller`](https://github.com/wjwei-handsome/wgatools/blob/main/src/parser/cigar.rs#L314):

| Category | target, query | Action |
|---|---|---|
| `=` | base, same base | advance both offsets, set `after_m=true` |
| `X` | base, different base | per column emit one SNP record (if `--snp`) |
| `I` | `-`, base | emit one INS record iff `len > svlen_cutoff` and `after_m` |
| `D` | base, `-` | emit one DEL record iff `len > svlen_cutoff` and `after_m` |
| `W` | `-`, `-` | ignore |

The `after_m` guard prevents emitting an indel that lacks an anchor base
(i.e. when an INS/DEL appears at the very start of a block or directly after
another gap run); the corresponding offsets are still advanced. Each block
is processed in chunks of `--chunk-size` (default 1 Mb), and
[`find_safe_chunk_boundary`](https://github.com/wjwei-handsome/wgatools/blob/main/src/tools/caller.rs#L159)
extends the chunk past any gap run ≥ `svlen_cutoff` so SVs are never split.

Negative-strand blocks emit a synthetic `<INV>` record at the block start
and tag every nested variant with `INV_NEST=TRUE`.

## Example-data comparison rerun

Rerun on 2026-11-13 with the bundled `example_data/example.maf/*.maf` files
(8 samples, one 1 Mb contig).

| Software | Version | How it was run |
|---|---|---|
| ARGprep | v1.12 (`scripts/maf_to_sites.py` at commit `0204aef`), Python 3.11 | all 8 samples together, settings below |
| wgatools | 1.1.0 (bioconda) | `wgatools call --snp --svlen 50`, one run per sample MAF |

wgatools 1.1.1 (released 2026-11-11) changes only `dotplot`; its `call` code
is identical to 1.1.0.

To reproduce, run on a compute node:

```bash
python scripts/compare_wgatools.py --wgatools /path/to/wgatools --out-dir OUT
```

It writes `comparison.md` and `comparison.json` (including both tools'
versions) to `OUT` and exits non-zero if any shared call disagrees on REF/ALT.

Two ARGprep settings are useful for interpreting the result:

- **Example/default missingness** (`max_missing_fraction: 0.0`,
  `mask_indel_adjacent_snps: false`): ARGprep emitted 9,578 multi-sample
  variant records. Across the eight samples, every ARGprep non-reference SNP
  call was present in the corresponding wgatools single-sample VCF
  (25,680 / 25,680 calls, zero REF/ALT mismatches). wgatools emitted 5,114
  additional SNP calls. These are expected: wgatools calls from one pairwise
  alignment at a time, while ARGprep masks a site if any other sample is
  missing at that reference base.
- **Permissive missingness** (`max_missing_count: 8`,
  `mask_indel_adjacent_snps: false`): ARGprep emitted 11,867 multi-sample
  variant records. ARGprep and wgatools made exactly the same SNP calls for
  every sample: 30,794 / 30,794 shared, with zero REF/ALT mismatches and no
  calls unique to either tool.

Per-sample SNP counts from the strict example setting:

| Sample | ARGprep non-ref SNPs | wgatools SNPs | Shared | wgatools-only |
|---|---:|---:|---:|---:|
| sample1 | 3,254 | 3,923 | 3,254 | 669 |
| sample2 | 3,194 | 3,872 | 3,194 | 678 |
| sample3 | 3,226 | 3,919 | 3,226 | 693 |
| sample4 | 3,160 | 3,768 | 3,160 | 608 |
| sample5 | 3,186 | 3,760 | 3,186 | 574 |
| sample6 | 3,238 | 3,826 | 3,238 | 588 |
| sample7 | 3,232 | 3,839 | 3,232 | 607 |
| sample8 | 3,190 | 3,887 | 3,190 | 697 |
| **Total** | **25,680** | **30,794** | **25,680** | **5,114** |

The SNP counts are identical to those from the earlier rerun (2026-07-05,
wgatools 1.1.0).

wgatools also emitted 522 INS and 665 DEL records (all longer than 50 bp)
across the eight example VCFs. The 2026-07-05 write-up reported 550 INS and
684 DEL with the same wgatools version; how those were counted was not
recorded, and lowering `--svlen` to 49 or 45 does not reproduce them. ARGprep
intentionally represents indel-affected reference positions in the BED masks
rather than as VCF indel records.

## Would they make the same SNP calls?

Restricted to the same single pairwise MAF (one sample vs. reference), the
two tools will not necessarily produce identical SNP call sets. The disagreements are
concentrated at:

1. **Indel-adjacent SNPs.** With `mask_indel_adjacent_snps: true` (off by
   default), ARGprep drops any SNP whose reference position immediately
   flanks an insertion or deletion in any sample. wgatools emits every `X`
   column regardless of neighbors, so those SNPs survive there but are
   filtered out by ARGprep when this option is enabled.
2. **Non-ACGT bases.** ARGprep masks any site where the reference is
   non-ACGT and treats `N` or ambiguity in a sample as *missing*, so it
   cannot drive a SNP call. wgatools' classifier only checks `c1 == c2`,
   so `N` vs. `A` is an `X` column and gets emitted as a SNP record
   (with `N` in REF or ALT).
3. **Overlapping or duplicate alignment blocks.** ARGprep's
   [`_assign_code`](scripts/maf_to_sites.py) collapses conflicting
   calls at the same reference position into `?` (missing), so the site
   cannot be called as a SNP. wgatools walks each block independently;
   overlapping blocks can yield duplicate or contradictory records.
4. **Missingness threshold (multi-sample runs only).** ARGprep can mask
   a site because *another* sample is missing there, suppressing a SNP
   that wgatools would still call from the single pair you fed it.
5. **Multi-allelic vs. pairwise.** With multiple samples, ARGprep can
   emit one record with multiple ALT alleles; running wgatools on each
   pair independently produces one ALT per record. The *positions* can
   still agree, but per-record content will not.

Where they agree: at any reference position covered by exactly one
alignment block, with ACGT bases in both reference and query, and not
adjacent to an indel, both tools call the same SNP with the same REF/ALT
bases. The disagreements cluster near indels, ambiguous bases, and
overlapping blocks.

## Why the algorithms diverge

The differences flow from a single design choice: **what is a "site"?**

- ARGprep treats the *reference contig* as the site space. Every reference
  base must end up either in the VCF (with a confident multi-sample call)
  or in the BED mask (with a reason). This is what an ARG estimator needs:
  invariant bases anchor the tree, missing/uncertain bases must be told
  apart from confidently-invariant bases, and indels are nuisance signal to
  be excluded rather than variants to be called.
- wgatools treats *alignment columns* as the site space. Only columns that
  exist in the pairwise alignment can produce records, and each variant
  type (SNP, INS, DEL, INV) is encoded directly in standard VCF/SV
  conventions. This is what a pairwise variant-discovery pipeline needs:
  identify what is different between target and query, including SVs.

As a result:

- ARGprep can answer "is this reference base invariant across all samples?"
  while wgatools cannot.
- wgatools can answer "what is the structural variation between this pair?"
  while ARGprep deliberately discards that information.
- ARGprep's missingness/threshold filters have no analogue in wgatools, and
  wgatools' SV-length cutoff and inversion handling have no analogue in
  ARGprep.
