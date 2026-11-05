# Helper scripts

These are standalone utilities that are **not** part of the Snakemake workflow. No rule
invokes them, they are not run by `snakemake`, and their outputs are not part of the
pipeline's deliverables. Run them by hand against a completed run's outputs.

- [popstats.py](#popstatspy) — pi, Watterson's theta, Tajima's D and a downsampled SFS from a FASTA with missing data
- [window_to_fasta.py](#window_to_fastapy) — reference-anchored window FASTA from the all-sites VCF
- [maf_to_fasta.py](#maf_to_fastapy) — gapped multiple alignment for a window, from the MAFs
- [simulate_msprime_indels.py](#simulate_msprime_indelspy) — msprime test-data simulator

The workflow's own scripts (`split_maf_by_contig.py`, `maf_to_sites.py`, `summary_report.py`, the MAF-QC scripts and `common.py`) are described in the main [README](../README.md).

## `popstats.py`

Computes nucleotide diversity (pi), Watterson's theta, Tajima's D and a folded,
downsampled site frequency spectrum (SFS) from an aligned FASTA that contains
missing data. It handles variable per-site sample size, so it suits alignments
where indels or masking leave different numbers of sequences at different sites.

```
python scripts/popstats.py aln.fa --min-seqs 20 [--sfs-n 16] \
    [--window 10000 --overlap 5000] [--out mystats]
```

### Input

An aligned FASTA (optionally `.gz`) in which every sequence has the same length.
`A/C/G/T` (either case) are called bases. Anything else (`N`, `-`, `?`, IUPAC
codes) is treated as missing, so alignments with gap characters, such as the output of
[maf_to_fasta.py](#maf_to_fastapy), can be used directly (a gap counts as missing, not as a
fifth allele).

### Flags

| Flag | Default | Meaning |
|---|---|---|
| `fasta` (positional) | required | Aligned FASTA, plain or `.gz`. |
| `--min-seqs X` | required | Keep only sites with at least X called sequences. Applies to every statistic. Must be >= 2. |
| `--sfs-n N` | `--min-seqs` | Number of chromosomes to downsample the SFS to. Must satisfy 2 <= N <= `--min-seqs`, so every retained site can be projected. |
| `--window W` | off | Window size in alignment columns. Without it the whole alignment is one window. |
| `--overlap Y` | `0` | Overlap between consecutive windows, so windows advance by W - Y. Requires 0 <= Y < W. |
| `--out PREFIX` | `popstats` | Output prefix. |

### How missing data are handled

Let n_i be the number of called sequences at site i. Only sites with
n_i >= `--min-seqs` are retained, and the number of retained sites (including
monomorphic ones) is the denominator for per-site values.

- **pi:** at each site, the mean over the n_i called sequences of pairwise
  differences, `1 - sum_a c_a(c_a-1) / (n_i(n_i-1))`, where c_a is the count of
  allele a. So the number of pairwise comparisons follows the data at each site.
  `pi_sum` is the sum over retained sites. `pi` is `pi_sum` divided by the number
  of retained sites.
- **Watterson's theta:** `S / sum_i a1(n_i)` over retained sites, where
  `a1(n) = sum_{k=1}^{n-1} 1/k` and S is the number of segregating sites (any
  retained site with two or more alleles). The result is per site.
- **Tajima's D:** `(pi_sum - S/A) / sqrt(e1*S + e2*S*(S-1))`, with A and the
  variance coefficients computed from the per-site n_i, averaged across retained
  sites (using the mean n_i in place of n). This is exact when n_i is constant
  across retained sites and an approximation when it varies. It is `NaN` when
  S = 0.
- **SFS:** folded (no outgroup is used) and restricted to biallelic sites. Each
  site with minor-allele count k out of n_i is projected down to `--sfs-n`
  chromosomes by hypergeometric sampling (X drawn without replacement from the
  n_i called chromosomes), so each retained biallelic site contributes a total
  weight of 1 in expected counts. Sites with three or more alleles are excluded
  from the SFS but are included in pi and S.

### Windows

Windows are in **alignment columns**, not genomic coordinates, and are reported
1-based inclusive. With `--window W --overlap Y` they start every W - Y columns.
The last window may be shorter than W. A trailing window that lies entirely
inside the previous one is dropped. A window with no retained sites is reported
with `n_sites` 0 and `NaN` statistics.

### Output

`PREFIX.stats.tsv`, one row per window:

| Column | Meaning |
|---|---|
| `start`, `end` | Window bounds (1-based, inclusive). |
| `n_sites` | Retained sites (n_i >= `--min-seqs`). |
| `n_biallelic` | Retained sites with exactly two alleles. |
| `mean_n` | Mean n_i over retained sites. |
| `S` | Segregating sites among retained sites. |
| `pi_sum` | Summed per-site pi. |
| `pi` | `pi_sum / n_sites`. |
| `theta_w` | Watterson's theta per site. |
| `tajima_d` | Tajima's D. |

`PREFIX.sfs.tsv`, one row per window with the same `start`/`end`:

| Column | Meaning |
|---|---|
| `proj_invariant` | Expected number of biallelic sites that look monomorphic after projecting to `--sfs-n` chromosomes (minor count 0). Fully monomorphic sites are not included. |
| `k1` ... `kM` | Expected number of sites with minor-allele count 1..M in the downsampled sample, where M = floor(`--sfs-n` / 2). |

### Examples

Whole-alignment statistics, keeping sites with at least 10 called sequences:

```
python scripts/popstats.py aln.fa --min-seqs 10 --out all
```

10 kb windows overlapping by 5 kb, SFS downsampled to 12 chromosomes:

```
python scripts/popstats.py aln.fa --min-seqs 16 --sfs-n 12 \
    --window 10000 --overlap 5000 --out win10k
```

### Tests

```
pytest tests/test_popstats.py            # small hand-checked cases
python tests/test_popstats_msprime.py    # compares to tskit on msprime loci with varying n
                                         # (also runs under pytest)
```

[tests/test_popstats_msprime.py](../tests/test_popstats_msprime.py) simulates
ten 1 kb loci with different sample sizes, concatenates them into one alignment
with `N` filling the missing rows, and checks `popstats.py` window output
against tskit statistics computed on each locus alone. Run both through
`hpc_run` on the cluster.

Good inputs are [window_to_fasta.py](#window_to_fastapy) output (reference-anchored, deletions
as `N`) and [maf_to_fasta.py](#maf_to_fastapy) output (gapped alignment).

## `window_to_fasta.py`

[window_to_fasta.py](window_to_fasta.py) builds a reference-anchored
multi-FASTA for a single window, so you can eyeball a region as sequence rather than as
a variant table:

```bash
python scripts/window_to_fasta.py \
  --vcf results/sites/combined.<contig>.all_sites.vcf \
  --bed-dir results/sites \
  --reference /path/to/reference.fa \
  --contig <contig> --start <start> --end <end> \
  --out window.fa
```

`--start`/`--end` are 1-based inclusive. Output is the reference sequence followed by one
sequence per sample. It builds each sample by starting from the reference, overwriting
positions where that sample's `GT` names a different allele, then overwriting every
interval in `combined.<contig>.<sample>.missing.bed` with `N`.

It requires `all_sites.vcf` rather than the variant-only VCF: starting from the reference
and patching in differences is only sound if every invariant position was positively
called as matching reference. The variant-only VCF cannot distinguish "called, matches
reference" from "never assessed".

**Indels are not represented, in either direction.** Every output sequence is exactly the
window length and contains no gap characters, because the whole pipeline works in
reference coordinates:

- An **insertion** in a sample consumes no reference position, so the inserted bases are
  dropped and leave no trace in the FASTA.
- A **deletion** in a sample is recorded as *missing*, not as a gap, so it appears as `N`.

An `N` is therefore ambiguous: it may mean the sample never aligned there, that the
sample carries a deletion, that the base was `N`/ambiguous in the MAF, or that it fell
below `quality_min`. The per-sample `.missing.bed` merges those same cases and cannot
separate them either; recovering that distinction means going back to the MAF. Treat the
output as a reference-anchored substitution view, not as a true alignment.

The window is held in memory once per sample, so this is intended for kilobase-scale
windows, not whole chromosomes.

## `maf_to_fasta.py`

[maf_to_fasta.py](maf_to_fasta.py) builds a **gapped multiple alignment**
for a window by going back to the MAFs, so indels survive: a deletion becomes a gap
column and an insertion adds columns that every other sample pads with gaps. This is the
thing `window_to_fasta.py` deliberately is not. It reads the per-contig MAF chunks the
workflow already writes under `results/maf_by_contig/`, so it does not need the
`all_sites.vcf`:

```bash
python scripts/maf_to_fasta.py \
  --maf-chunk-root results/maf_by_contig \
  --reference-fasta /path/to/reference.fa \
  --contig <contig> --start <start> --end <end> \
  --out window.fa
```

`--start`/`--end` are 1-based inclusive and default to the whole contig. Pass
`--maf-dir` instead of `--maf-chunk-root` to read whole-genome per-sample MAFs directly
(`<sample>.maf` / `.maf.gz`), and `--quality-bed-dir` / `--quality-min` to apply the same
per-sample quality masking the workflow uses. `--max-columns` caps the alignment width.

The merge is **reference-anchored**: each reference base gets one anchor column, and the
pairwise alignments are related to each other only through it. Four consequences are
worth stating outright.

- **Insertion homology is asserted, not established.** Insertions are left-anchored to
  the preceding reference base and left-aligned within that anchor's columns, padded on
  the right to the widest insertion any sample carries there. Two samples inserting at
  the same anchor therefore share columns, but the pairwise MAFs cannot show that they
  are the same event. An insertion preceding the window's first base anchors outside the
  window and is dropped. Leading insertions at a block boundary are retained when
  their preceding reference coordinate lies inside the window.
- **Conflicting blocks degrade to `N`.** Where two alignment blocks cover one reference
  position with different calls, the anchor becomes `N`, matching the `?` that
  `maf_to_sites.py` assigns. The analogous rule for insertions is that conflicting
  sequences at one anchor become a run of `N` as long as the longer of the two.
  This includes an insertion conflicting with explicit absence in another block
  spanning that slot. A block ending or starting at the slot does not establish absence.
- **By default it does not agree with the VCF.** It applies no missingness threshold and
  does not know which sites the workflow masked, so it can show a base where
  `mask.bed` masked the site. Pass `--mask-bed results/sites/combined.<contig>.mask.bed`
  to reconcile them: masked reference positions render as `N` for every sample, and
  insertions anchored there are dropped.
- **Width can blow up.** An insertion-rich region across many samples can make a 10 kb
  window far wider than 10 kb, and memory is columns × samples. `--max-columns` (default
  5,000,000) rejects reference spans above the cap before loading sample alignments,
  and checks the insertion-expanded width before rendering rows. It is a width limit,
  not a memory guarantee: input parsing, stored projections, and sample count still
  consume memory. Use an explicit small window for large genomes.

Characters mean: `A`/`C`/`G`/`T` a called base; `-` a gap, meaning either a deletion in
that sample or padding beside another sample's insertion; `N` no information — the sample
never aligned there, the MAF base was ambiguous, it fell below `--quality-min`, blocks
conflicted, or `--mask-bed` masked the position. Unlike `window_to_fasta.py`, `N` no
longer conflates deletions with the rest.

## `simulate_msprime_indels.py`

The repository includes [simulate_msprime_indels.py](simulate_msprime_indels.py) for generating haploid test datasets with msprime SNP variation plus branch-based indels on the tree sequence. Note that these simulations are not intended to be evolutionarily accurate, but simply to give a reasonable example data.

Example:

```bash
python scripts/simulate_msprime_indels.py \
  --sequence-length 1000000 \
  --num-samples 8 \
  --theta 0.01 \
  --rho 0.01 \
  --ne 10000 \
  --indel-rate 1e-8 \
  --indel-lambda 0.001 \
  --seed 8675309 \
  --out-prefix example_data/example
```

Outputs:

- `<prefix>.reference.fa`
- `<prefix>.samples.fa`
- `<prefix>.indels.tsv`
- `<prefix>.summary.tsv`
- `<prefix>.maf/`

Summary fields include:

- `seed`
- `sequence_length`
- `reference_bp_with_indel_in_ge1_sample`
- `total_snps`
- `snps_without_overlapping_indel`
