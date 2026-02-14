# ARG Pipeline (Snakemake)

This repository provides a Snakemake workflow that converts AnchorWave MAFs into per-contig merged gVCFs, splits them into clean/filtered/invariant sets, and produces mask bedfiles for downstream ARG inference (see `logic.md` for details).

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

If your reference FASTA does **not** have an index (`.fai`), either create one (`samtools faidx`) or set `contigs:` explicitly in `config.yaml`.

## Run 

The pipeline can be run one of two ways, both from the repo root. **It is recommended you first run on the example data provided to ensure the pipeline works on your system.**

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
- `--notemp`: keep temporary intermediates (see Notes)

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
- `results/split/combined.<contig>.coverage.txt` : split coverage validation summary
- `results/split/combined.<contig>.accessible.npz` : boolean accessibility array (union of clean + invariant sites), for scikit-allel statistics
- `results/summary.html` : HTML summary of jobs run, outputs created, and warnings

## Notes

- If `bgzip_output: true`, the `.inv`, `.filtered`, `.clean`, and `.missing.bed` files will have a `.gz` suffix.
- All gzipped outputs in this pipeline use bgzip (required for `tabix`).
- `scripts/dropSV.py` removes indels larger than `drop_cutoff` (if set in `config.yaml`).
- `scripts/split.py` supports `--filter-multiallelic` and `--bgzip-output` (toggle via `config.yaml`).
- `scripts/filt_to_bed.py` merges `<prefix>.filtered`, `<prefix>.missing.bed`, and `dropped_indels.bed` into a final mask bed.
- `make_accessibility` builds a per-contig accessibility array from the union of `combined.<contig>.clean` and `combined.<contig>.inv` using the reference `.fai` to size the array. The output is a compressed NumPy archive containing a boolean array named `mask`, intended for scikit-allel statistics.
- Warning: The workflow defaults to haploid genotyping (`ploidy: 1`) and has not been validated on diploid genome assemblies.
- Optional: enable `vt_normalize: true` in `config.yaml` to normalize merged gVCFs with `vt normalize` after `SelectVariants`.
- If GenomicsDBImport fails with a buffer-size error, increase `genomicsdb_vcf_buffer_size` and `genomicsdb_segment_size` in `config.yaml` (set them above your longest gVCF line length).
- Large intermediate files are marked as temporary and removed after a successful run (per-sample gVCFs, cleaned gVCFs, per-contig split gVCFs, and the GenomicsDB workspace). Use `snakemake --notemp` if you want to preserve them for debugging or reruns.
- Resource knobs (memory/threads/time) and GenomicsDB buffer sizes are configurable in `config.yaml` (e.g., `merge_contig_mem_mb`, `maf_to_gvcf_*`, `genomicsdb_*`).
- To cap concurrent contig-merge jobs on SLURM, set `merge_gvcf_max_jobs` in `config.yaml` (used by the profile as a global `merge_gvcf_jobs` resource limit).
- To cap the SLURM array concurrency for `scripts/maf_to_gvcf.sh`, set `maf_to_gvcf_array_max_jobs` in `config.yaml` (default 4).

## Changes since v0.1

- Added HTML summary report with embedded SVG histograms and expanded output details.
- Split logic tightened: clean sites now require all samples called; missing GTs are routed to filtered.
- Invariant/filtered/clean outputs are enforced as mutually exclusive per position; filtered BED spans now respect END/REF lengths and subtract inv/clean.
- Merged gVCFs are produced via GATK SelectVariants with genotype calling; TASSEL `outputJustGT` default set to `false` to retain likelihoods for calling.
- Added accessibility mask generation (`combined.<contig>.accessible.npz`) for scikit‑allel workflows.
- New/expanded validation and tests: split coverage checks, filtered‑bed tests, integration tests gated by `RUN_INTEGRATION=1`.
- Example data regenerated via msprime with indels and missing data and AnchorWave‑style MAF formatting.
- `check_split_coverage.py` now reports overlap intervals with file names to aid debugging.
- `filt_to_bed.py` filters masks to the target contig, preventing cross‑contig lines in `combined.<contig>.filtered.bed`.
- SLURM default resources now read `default_*` from `config.yaml` instead of hardcoded profile values.

## Changes since v0.2

- Moved HTML summary generation into `scripts/summary_report.py` and simplified `Snakefile`.
- Corrected example MAF inputs so `example_data/*.maf.gz` are valid gzip files.
- Updated split classification so `ALT=<NON_REF>`-only records are treated as invariant.
- Updated SINGER clean-output formatting to strip `<NON_REF>` while preserving genotype/sample fields.
- Added `logic.md` with detailed site-routing/filtering logic and concrete examples.
- Added `results/split/combined.<contig>.coverage.txt` to documented workflow outputs.
- Added split-test coverage for invariant/nonref and genotype-preserving clean formatting.
- Updated SLURM profile default resources to numeric values to avoid resource conversion/submission errors.
- Added merge_gvcf_max_jobs pipeline concurrency control

## Downstream Uses

### ARG estimation

Use Nate Pope's SINGER Snakemake [pipeline](https://github.com/nspope/singer-snakemake) with `combined.<contig>.clean` and `combined.<contig>.filtered.bed` as inputs.

### Population genetic statistics 

If you use scikit-allel, you can use the `combined.<contig>.clean` VCF and load the accessibility mask like this:

```python
import numpy as np
mask = np.load("results/split/combined.1.accessible.npz")["mask"]
```
