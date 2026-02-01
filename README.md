## Workflow Overview

This workflow documents taking `.maf` files of assemblies aligned to a reference, creating a `.gvcf`, filtering, and then estimating an ARG. It assumes this is done one chromosome at a time.

## 1 Assemble a gvcf

### 1A Align genomes to reference

Align each assembly to the reference using [anchorwave](https://github.com/baoxingsong/AnchorWave). This step produces per-sample MAF files and should be done separately by chromosome.

### 1B Convert to gvcf
Use `maf_to_gvcf.sh` to run tassel on a directory of `.maf` files.
It expects a reference FASTA located in the same directory as the `.maf` files
(either `reference.fa` or the only `.fa/.fasta` file in that directory).

Submit with: `sbatch --array=0-<N-1>%4 maf_to_gvcf.sh -m /path/to/maf_dir -d /path/to/out_dir` where `N` is the number of `.maf` files in the directory. The `%4` cap limits to 4 concurrent array tasks.
You can also run `maf_to_gvcf.sh -m /path/to/maf_dir -d /path/to/out_dir` directly; it will auto-submit an array sized to the number of `.maf` files with the same 4-task cap. These scripts should be launched from the repo root so `logs/` exists for SLURM output.
Logs are written to `logs/` alongside the SLURM script.

### 1C Make joint gvcf
To merge per-sample gVCFs after filtering large indels, use `gatk_merge_gvcf.sh`.
It runs `dropSV.py`, bgzip/tabix, and then GATK GenomicsDBImport + GenotypeGVCFs.
Submit with: `sbatch gatk_merge_gvcf.sh -g /path/to/gvcf_dir -r /path/to/reference.fa [-l interval] [-c max_indel_len]` (default cutoff: 9101264). Run from the repo root so SLURM can write to `logs/`.

      OPTIONAL: `dropSV.py` can be run separately to remove large indels. Run `./dropSV.py -h` for options. This writes       `cleangVCF/dropped_indels.bed` (full-span intervals). 

## 2 GVCF parsing
### 2A Clean gvcf 
 
Use `split_and_mask.sh` to run `split.py` followed by `filt_to_bed.py` in one SLURM job:
`sbatch split_and_mask.sh -p /path/to/prefix -d <depth> [--filter-multiallelic] [--no-gzip] [--no-merge]` (run from the repo root so `logs/` is present).
Normally you will want to set depth equal to your sample size. 
In some files, for example, depth is recorded as 30 for each individual, so you should set depth to 30 x sample size.
If your samples are not inbred, you may need to change this by a factor of two. 
In addition, if you run the script with `--filter-multiallelic`, this will send multi-allelic sites to the `.filtered` file described below. 

VCFs from anchorwave often have genotypes using depth like:

`13      216881  .       G       A,<NON_REF>     .       .       DP=120  GT:AD:PL:DP     .:.:.:. .:30,0,0:0,90,90:30     .:.:.:. .:30,0,0:0,90,90:30     .:30,0,0:0,90,90:30     .:0,30,0:90,90,0:30     .:.:.:.`

`split.py` will automatically reformat to genotypes for SINGER when needed based on this pattern.

This script writes three files, `.inv`, `.filtered`, and `.clean`. Each includes the regular header.
It also writes a `.missing.bed` file (no header) listing bp positions absent from the input VCF (i.e., not covered by invariant END ranges, indel spans, or variable sites), in BED 0-based half-open format.
File outputs will be large when unzipped, it is recommended to run with `--gzip-output` to automatically zip output files.
Writes to stderr log of how many bp (expanding `END` segments) were written to each file.
Your gvcf **must** have invariant sites. If there are no invariant sites, go back to [step 2](https://github.com/rossibarra/jri-arg/blob/main/README.md#2-gvcf-parsing).

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
**Note:** If you run the script with `--filter-multiallelic` this will send multi-allelic site to the `.filtered` file instead. 

The bash script then takes your `.filtered` file plus any large indels removed by `dropSV.py` and any missing positions identified by `split.py`. `filt_to_bed.py` takes a gVCF/VCF filename (or its prefix) and builds a merged mask from `<prefix>.filtered`, `<prefix>.missing.bed`, and a dropped-indels bedfile. By default it looks for `cleangVCF/dropped_indels.bed` relative to the prefix directory; you can override this with `--dropped-bed /path/to/dropped_indels.bed`. It exits with an error if required files are missing. It also checks that filtered bed bp + `.inv` bp + `.clean` bp equals the chromosome length inferred from the gVCF header (or last bp if no header length).

Using `--no-merge` will result in a bigger bedfile with many small, contiguous regions and is not recommended.

Alternatively, you can run these individually.

            Run `split.py` using `python3 split.py --depth=<depth> <filename.vcf>`. 

            Run usfilt_to_bed.py using: `python3 filt_to_bed.py /path/to/<prefix>.gvcf[.gz] [--no-merge]`. Output is `<prefix>.filtered.bed`.

### 2B validate

As alignment software and GATK versions may produce gvcfs of different formats, you should validate your output makes sense. 
Some suggestions include:

- `grep -v "#" test.inv | cut -f 5 | sort -n | uniq` -- check that invariant sites file only has "NON_REF" as an ALT allele
- `grep -v "#"  test.clean | cut -f 8 | sort -n | uniq` -- check that the INFO field has all "DP=`depth`" values
- `grep -v "#" test.filtered | grep -v "*" | less -S` -- scroll through some filtered records (excluding indels) and check they should all be removed and none are invariant or good SNPs

## 3 ARG estimation

Use Nate Pope's [snakemake pipeline](https://github.com/nspope/singer-snakemake/tree/main).
Using the steps above, there is no need to have a filter file. 
Use the `.clean` vcf and `.filtered.bed` made above as the mask bedfile. 
Please make sure your recombination 'hapmap' file extends to the end of the chromosome. 

## 4 ARG processing and 5 Ne modeling (under construction)

Separate scripts are being developed for these steps. In the meantime, one approach to the functionality here can be found in this interactive jupyter notebook [argcheck.ipynb](argcheck.ipynb). You will need to have `msprime` `tskit` `tszip` `demes` `demesdraw` and `yaml` among other packages. There are many places in the notebook where you will need to **modify** the code or variables to suite your files/system. Example outputs for maize are shown in the notebook preview.
