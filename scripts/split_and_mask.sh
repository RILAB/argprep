#!/bin/bash
#SBATCH --job-name=split_and_mask
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err

set -euo pipefail

# One job: split gVCF into .inv/.filtered/.clean and build mask bed.
usage() {
  echo "Usage: $0 -p <prefix> -d <depth> [--filter-multiallelic] [--no-gzip] [--no-merge]"
  echo "  -p  Prefix for .inv/.filtered/.clean outputs"
  echo "  -d  Depth cutoff passed to split.py"
  echo "  --filter-multiallelic  Filter multi-allelic SNPs in split.py"
  echo "  --no-gzip              Do not gzip split.py outputs (default)"
  echo "  --no-merge             Do not merge intervals in filt_to_bed.py"
  exit 1
}

PREFIX=""
DEPTH=""
FILTER_MULTIALLELIC="false"
NO_GZIP="true"
NO_MERGE="false"

while [ "$#" -gt 0 ]; do
  case "$1" in
    -p) PREFIX="$2"; shift 2 ;;
    -d) DEPTH="$2"; shift 2 ;;
    --filter-multiallelic) FILTER_MULTIALLELIC="true"; shift ;;
    --no-gzip) NO_GZIP="true"; shift ;;
    --no-merge) NO_MERGE="true"; shift ;;
    *) usage ;;
  esac
done

if [ -z "$PREFIX" ] || [ -z "$DEPTH" ]; then
  usage
fi

# Logs go to submission directory logs/.
LOG_DIR="${SLURM_SUBMIT_DIR:-.}/logs"
mkdir -p "$LOG_DIR"

# Resolve script directory so we can find helper scripts when SLURM runs from a spool dir.
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  REPO_ROOT="$SLURM_SUBMIT_DIR"
else
  REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
fi
SCRIPT_DIR="$REPO_ROOT/scripts"

# Find input VCF/gVCF from prefix.
INPUT_VCF=""
for ext in ".vcf.gz" ".vcf" ".gvcf.gz" ".gvcf"; do
  if [ -f "${PREFIX}${ext}" ]; then
    INPUT_VCF="${PREFIX}${ext}"
    break
  fi
done

if [ -z "$INPUT_VCF" ]; then
  echo "ERROR: could not find input VCF/GVCF for prefix '${PREFIX}'"
  exit 1
fi

# Build split.py args based on requested options.
SPLIT_ARGS=(python3 "${SCRIPT_DIR}/split.py" --depth="$DEPTH")
if [ "$FILTER_MULTIALLELIC" = "true" ]; then
  SPLIT_ARGS+=(--filter-multiallelic)
fi
if [ "$NO_GZIP" != "true" ]; then
  SPLIT_ARGS+=(--gzip-output)
fi
SPLIT_ARGS+=("$INPUT_VCF")

# Run split.py to produce .inv/.filtered/.clean/.missing.bed.
"${SPLIT_ARGS[@]}"

# Build mask bed from filtered + missing + dropped indels.
FILT_ARGS=(python3 "${SCRIPT_DIR}/filt_to_bed.py" "$PREFIX")
if [ "$NO_MERGE" = "true" ]; then
  FILT_ARGS+=(--no-merge)
fi

"${FILT_ARGS[@]}"
