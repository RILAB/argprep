# ARG Pipeline (Snakemake)

This repository provides a Snakemake workflow that converts AnchorWave MAFs into per-contig merged gVCFs, splits them into clean/filtered/invariant sets, and produces mask bedfiles for downstream ARG inference.

## Requirements

- Conda (you may need to `module load conda` on your cluster)
- TASSEL (provided via the `tassel-5-standalone` submodule)
- GATK, Picard, htslib (installed in the conda env defined below)

## Setup

Clone with submodules so TASSEL is available:

```bash
git clone --recurse-submodules <repo>
```

Create and activate the environment (do this before running Snakemake):

```bash
module load conda
conda env create -f argprep.yml
conda activate argprep
```

## Configure

Edit `config.yaml` to point at your MAF directory and reference FASTA. At minimum you must set:

- `maf_dir`: directory containing `*.maf` or `*.maf.gz`
- `reference_fasta`: reference FASTA path (plain `.fa/.fasta` or bgzipped `.fa.gz`)
- `depth`: depth cutoff for `split.py` (defaults to the number of MAF files)

If your reference FASTA does **not** have an index (`.fai`), either create one (`samtools faidx`) or set `contigs:` explicitly in `config.yaml`.

## Run 

The pipeline can be run one of two ways, both from the repo root:

### On Slurm 

A default SLURM profile is provided under `profiles/slurm/`. Edit `profiles/slurm/config.yaml` to customize sbatch options if needed.
Defaults for account/partition and baseline resources are set in `config.yaml` (`slurm_account`, `slurm_partition`, `default_*`).

```bash
snakemake --profile profiles/slurm
```

### Locally

```bash
snakemake -j 8 
```

Common options:

- `-j <N>`: number of parallel jobs
- `--rerun-incomplete`: clean up partial outputs
- `--printshellcmds`: show executed commands

## Workflow Outputs

By default the workflow uses these locations (override in `config.yaml`):

- `gvcf/` : TASSEL gVCFs (`*.gvcf.gz`) from MAFs
- `gvcf/cleangVCF/` : cleaned gVCFs from `scripts/dropSV.py`
- `gvcf/cleangVCF/dropped_indels.bed` : bedfile of large indels
- `gvcf/cleangVCF/split_gvcf/` : per-contig gVCFs for merging
- `results/combined/combined.<contig>.gvcf.gz` : merged gVCF per contig
- `results/split/combined.<contig>.inv` : invariant sites
- `results/split/combined.<contig>.filtered` : filtered sites
- `results/split/combined.<contig>.clean` : clean sites
- `results/split/combined.<contig>.missing.bed` : missing positions
- `results/split/combined.<contig>.filtered.bed` : merged mask bed
- `results/summary.md` : markdown summary of jobs run, outputs created, and warnings
  - If `gzip_output: true`, the `.inv`, `.filtered`, `.clean`, and `.missing.bed` files will have a `.gz` suffix.

## Notes

- `scripts/dropSV.py` removes indels larger than `drop_cutoff` (if set in `config.yaml`).
- `scripts/split.py` supports `--filter-multiallelic` and `--gzip-output` (toggle via `config.yaml`).
- `scripts/filt_to_bed.py` merges `<prefix>.filtered`, `<prefix>.missing.bed`, and `dropped_indels.bed` into a final mask bed.

## Downstream ARG estimation

Use Nate Pope's Snakemake pipeline:

- `singer-snakemake` (included as a submodule)

For SINGER, use `combined.<contig>.clean` and `combined.<contig>.filtered.bed` as inputs.
