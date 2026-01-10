
## 1 Align genomes to reference

Align each assembly to the reference using [anchorwave](https://github.com/baoxingsong/AnchorWave).

## 2 Assemble a gvcf

Individual `.maf` files need to be converted to `.gvcf` and then combined to a single `.gvcf`. 
We recommend doing this separately by chromosome. 
Instructions for these steps are [here](https://github.com/baoxingsong/AnchorWave/blob/master/doc/GATK.md).

## 3 Clean gvcf

Run `split.py` using `python3 split.py --depth=<depth> <filename.vcf>`. 

This script writes three files, `.inv`, `.filtered`, and `.clean`. Each includes the regular header.

##### `.inv` 
Contains lines from vcf where `INFO` is `.` OR `INFO` contains `END=`, and expands `END=` lines for each bp. 
This means the number of non-header lines (e.g. from `wc -l`) should represent the number of invariant bp in the file.

##### `.filtered`
Contains lines from vcf where any of the following occur:

- `DP` is less than the `depth` parameter given
- the line contains `*` as an allele
- multi-bp REF allele
- multiallelic SNPs
- symbolic / non-ACGT alleles

##### `.clean`
Should contain only biallelic SNPs in vcf passing all checks.
