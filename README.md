## 1 Align genomes to reference

Align each assembly to the reference using [anchorwave](https://github.com/baoxingsong/AnchorWave).

## 2 Assemble a gvcf

Individual `.maf` files need to be converted to `.gvcf` and then combined to a single `.gvcf`. 
We recommend doing this separately by chromosome. 
Instructions for these steps are [here](https://github.com/baoxingsong/AnchorWave/blob/master/doc/GATK.md).

## 3 Clean gvcf

Assuming your gvcf is formatted like the [example file](https://github.com/RILAB/arg-ne/blob/main/test.vcf.gz).
Run `split.py` using `python3 split.py --depth=<depth> <filename.vcf>`. 

This script writes three files, `.inv`, `.filtered`, and `.clean`. Each includes the regular header.
File outputs may be large when unzipped, so run with `--gzip-output` to automatically zip output files.

##### `.inv` 
Contains lines from vcf where:
- `INFO` is `.`
-  `INFO` contains `END=`

For multi-bp lines (with `END=`, the script copies the line once for each bp. 
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

## 4 Prep for SINGER

`.clean` will be the SNP data you give to SINGER. You will also need a `.bed` format file of bp that are masked. 
Usually these are everything in your `.filtered` file.
`filt_to_bed.py` will take a vcf and make a bedfile. 
Run using: `python3 filt_to_bed.py <filtered vcf file> --merge`. 
Dropping the `--merge` will result in a bigger bedfile with many small, contiguous regions and is not recommended.
