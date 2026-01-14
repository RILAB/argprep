## 1 Align genomes to reference

Align each assembly to the reference using [anchorwave](https://github.com/baoxingsong/AnchorWave).

## 2 Assemble a gvcf

Individual `.maf` files need to be converted to `.gvcf` and then combined to a single `.gvcf`. 
We recommend doing this separately by chromosome. 
Instructions for these steps are [here](https://github.com/baoxingsong/AnchorWave/blob/master/doc/GATK.md).

## 3 Clean gvcf & prep for SINGER

### 3A Clean gvcf 
Assuming your gvcf is formatted like the [example file](https://github.com/RILAB/arg-ne/blob/main/test.vcf.gz) and is for a single chromosome. 
Please split any multi-chromosome gvcfs into individual chromosomes before continuing.
The script can read both gzipped and unzipped vcfs.

Run `split.py` using `python3 split.py --depth=<depth> <filename.vcf>`. 
Normally you will want to set depth equal to your sample size. 
In some files, for example, depth is recorded as 30 for each individual, so you should set depth to 30 x sample size.

This script writes three files, `.inv`, `.filtered`, and `.clean`. Each includes the regular header.
File outputs will be large when unzipped, it is recommended to run with `--gzip-output` to automatically zip output files.
Writes to stderr log of how many bp (expanding `END` segments) were written to each file.
Your gvcf **must** have invariant sites. If there are no invariant sites, go back to [step 2](https://github.com/RILAB/arg-ne/blob/main/README.md#2-assemble-a-gvcf).

##### `.inv` 
Contains lines from vcf where:
- `INFO` is `.`
-  `INFO` contains `END=`

For multi-bp lines (with `END=`, the script copies the line once for each bp. 
This means the number of non-header lines (e.g. from `wc -l`) should represent the number of invariant bp in the file.

##### `.filtered`
Contains lines from vcf where any of the following occur:

- the line contains `*` as an allele
- symbolic / non-ACGT alleles. Among other things, this removes bp with `N`. This is OK if you have assemblies and you assume genotyping error is effectively zero. If you cannot assume that, then these alleles would need to go into `.clean` but the script cannot currently do that.
- `DP` is less than the `depth` parameter given. Since we are using whole genome alignment we are assuming sites with DP < depth but no explicit indels "*" have missing data still due to structural variation of some sort and should be removed.

##### `.clean`
Should contain only biallelic SNPs in vcf passing all checks as well as mutliallelic SNPs with no indels that will be filtered out by SINGER snakemake pipleine later (and used to adjust the mutation rate). 
If you run the script with `--filter-multiallelic` this will send multi-allelic site to the `.filtered` file instead. 

### 3B Prep for SINGER

`.clean` will be the SNP data you give to SINGER. 
You will also need a `.bed` format file of bp that are masked. 
Usually these are everything in your `.filtered` file.
`filt_to_bed.py` will take a vcf and make a bedfile. 

Run using: `python3 filt_to_bed.py <vcf file of filtered snps> --merge`. 
Dropping the `--merge` will result in a bigger bedfile with many small, contiguous regions and is not recommended.

### validate

As alignment software and GATK versions may produce gvcfs of different formats, you should validate your output makes sense. 
Some suggestions include:

- `grep -v "#" test.inv | cut -f 5 | sort -n | uniq` -- check that invariant sites file only has "NON_REF" as an ALT allele
- `grep -v "#"  test.clean | cut -f 8 | sort -n | uniq` -- check that the INFO field has all "DP=<depth>" values
- `grep -v "#" test.filtered | grep -v "*" | less -S` -- scroll through some filtered records (excluding indels) and check they should all be removed and none are invariant or good SNPs

## 4 Run SINGER

Use Nate Pope's [snakemake pipeline](https://github.com/nspope/singer-snakemake/tree/main).
Using the steps above, there is no need to have a filter file. 
Use the bedfile made in step 4 as the mask bedfile. 

## 5 Edit and validate ARG, estimate and validate Ne models

This step requires interactive use of the jupyter notebook [argcheck.ipynb](argcheck.ipynb). You will need to have `msprime` `tskit` `tszip` `demes` `demesdraw` and `yaml` among other packages. There are many places in the notebook where you will need to **modify** the code or variables to suite your files/system. Example outputs for maize are shown in the notebook preview.
