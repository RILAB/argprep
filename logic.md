# Pipeline Logic

This file documents how sites are routed through the workflow and where they end up.

## High-Level Flow

1. `maf_to_gvcf` converts each sample MAF to sample gVCF.
2. `drop_sv.py` removes super-large indel records from sample gVCFs.
3. `split_gvcf_by_contig` extracts per-contig gVCFs.
4. `merge_contig` merges sample gVCFs into `results/combined/combined.<contig>.gvcf.gz`.
5. `split.py` splits merged gVCF into:
   - `combined.<contig>.inv`
   - `combined.<contig>.clean`
   - `combined.<contig>.filtered`
   - `combined.<contig>.missing.bed`
6. `filt_to_bed.py` builds `combined.<contig>.filtered.bed` from filtered/missing/dropped-indel masks (minus any clean/inv coverage).
7. `check_split_coverage.py` validates that clean + inv + filtered.bed exactly covers the contig.

## What Happens to Sites Not in the MAF

Positions with no merged gVCF record are written to `combined.<contig>.missing.bed` by gap tracking in `split.py`.

Example:
- If contig length is 10 and records only cover positions 3-5 and 8,
- missing positions are 1-2, 6-7, 9-10,
- written in BED (0-based, half-open), e.g.:

```bed
1\t0\t2
1\t5\t7
1\t8\t10
```

These missing intervals are included in `combined.<contig>.filtered.bed` (mask output).

## What Happens to Large Indels

`drop_sv.py` removes records where `abs(len(REF) - len(ALT)) > cutoff` (default cutoff: `9101264`).

Outputs:
- cleaned per-sample gVCFs in `gvcf/cleangVCF/`
- removed positions merged into `gvcf/cleangVCF/dropped_indels.bed`

`dropped_indels.bed` is later added to `combined.<contig>.filtered.bed` so downstream masking excludes those positions.

## Site Classification in `split.py`

Priority is important:

1. Invariant (`.inv`) has highest priority.
2. Then filtered (`.filtered`).
3. Remaining called variant sites go to clean (`.clean`).

For each position, records are grouped and made mutually exclusive.

### Invariant (`combined.<contig>.inv`)

A record is invariant if either is true:
- `INFO` contains `END=` (expanded across `POS..END` into one record per bp in `.inv`)
- `ALT == "."`
- `ALT` contains only `<NON_REF>` (for example `ALT=<NON_REF>`)

Example (goes to `.inv`):

```vcf
1\t180\t.\tT\t<NON_REF>\t.\t.\t.\tGT:AD\t0/0:30,0\t0/0:30,0\t0/0:30,0
```

### Filtered (`combined.<contig>.filtered`)

A record is filtered if any of these are true:
- `ALT` includes `*`
- `REF` length > 1
- `ALT` has non-ACGT allele (excluding `*` and `<NON_REF>`)
- `--filter-multiallelic` is enabled and there are multiple distinct A/C/G/T ALTs
- any sample GT is missing (`.`)

Examples (go to `.filtered`):

```vcf
# missing GT
1\t200\t.\tA\tG\t.\t.\tDP=20\tGT:AD\t0/1:10,10\t.:.\t0/0:20,0

# star allele
1\t300\t.\tC\t*,<NON_REF>\t.\t.\tDP=8\tGT\t0/1\t0/0\t0/0

# symbolic/non-ACGT ALT
1\t400\t.\tG\t<DEL>,<NON_REF>\t.\t.\tDP=12\tGT\t0/1\t0/0\t0/0
```

### Clean (`combined.<contig>.clean`)

A record goes to clean when it is not invariant and not filtered, and all samples are called.

Example (goes to `.clean`):

```vcf
1\t500\t.\tA\tG,<NON_REF>\t.\t.\tDP=30\tGT:AD\t0/1:15,15\t0/0:30,0\t1/1:0,30
```

In clean output, `<NON_REF>` is removed from ALT, but genotype/sample fields are preserved.
So this becomes:

```vcf
1\t500\t.\tA\tG\t.\t.\tDP=30\tGT:AD\t0/1:15,15\t0/0:30,0\t1/1:0,30
```

## SNP-Type Examples and Destination

- Biallelic called SNP (`A -> G`, no missing GT): `clean`
- SNP with missing sample GT: `filtered`
- SNP with `*` ALT: `filtered`
- SNP with non-ACGT symbolic ALT (`<DEL>`, etc.): `filtered`
- Multiallelic A/C/G/T SNP:
  - `clean` if `--filter-multiallelic` is off
  - `filtered` if `--filter-multiallelic` is on
- gVCF non-variant SNP-position record (`ALT=<NON_REF>` only or `ALT=.`): `inv`

## Final Mask (`combined.<contig>.filtered.bed`)

`filtered.bed` is built from:
- split filtered records
- `combined.<contig>.missing.bed` (sites absent from merged gVCF)
- `gvcf/cleangVCF/dropped_indels.bed` (large indel masks)

Then any intervals overlapping clean/inv are subtracted, so outputs remain mutually exclusive.
