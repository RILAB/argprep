#!/bin/bash
#SBATCH --job-name=maf_to_gvcf
#SBATCH --array=0-0
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

MAF_DIR=""
OUT_DIR=""

while getopts "m:d:" opt; do
  case "$opt" in
    m) MAF_DIR="$OPTARG" ;;
    d) OUT_DIR="$OPTARG" ;;
    *) ;;
  esac
done

if [ -z "$MAF_DIR" ] || [ -z "$OUT_DIR" ]; then
  echo "Usage: $0 -m <maf_dir> -d <out_dir>"
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="$SCRIPT_DIR/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

REPO_ROOT="$SCRIPT_DIR"
TASSEL_DIR="$REPO_ROOT/tassel-5-standalone"
FILL_GAPS="false"
SAMPLE_SUFFIX="_anchorwave"

# Reference FASTA must live in the MAF directory.
if [ -f "$MAF_DIR/reference.fa" ]; then
  REF_FASTA="$MAF_DIR/reference.fa"
else
  REF_CANDIDATES=()
  while IFS= read -r line; do
    REF_CANDIDATES+=("$line")
  done < <(find "$MAF_DIR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \) | sort)
  if [ "${#REF_CANDIDATES[@]}" -ne 1 ]; then
    echo "ERROR: expected a single reference FASTA in $MAF_DIR (or reference.fa)."
    exit 1
  fi
  REF_FASTA="${REF_CANDIDATES[0]}"
fi

MAF_FILES=()
while IFS= read -r line; do
  MAF_FILES+=("$line")
done < <(find "$MAF_DIR" -maxdepth 1 -type f -name "*.maf" | sort)
if [ "${#MAF_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no .maf files found in $MAF_DIR"
  exit 1
fi

# If not running as a SLURM array task, submit ourselves with an array sized to the MAF count.
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
  N="${#MAF_FILES[@]}"
  ARRAY="0-$((N-1))%4"
  echo "Submitting SLURM array with $N tasks (max 4 running at a time): $ARRAY"
  sbatch --array="${ARRAY}" "$0" -m "$MAF_DIR" -d "$OUT_DIR"
  exit 0
fi

MAF_FILE="${MAF_FILES[$SLURM_ARRAY_TASK_ID]}"
BASE="$(basename "$MAF_FILE" .maf)"
SAMPLE_NAME="${BASE}${SAMPLE_SUFFIX}"
REF_BASE="$(basename "$REF_FASTA")"
REF_BASE="${REF_BASE%.fa}"
REF_BASE="${REF_BASE%.fasta}"

GVCF_OUT="${OUT_DIR}/${BASE}To${REF_BASE}.gvcf"
LOG_OUT="${LOG_DIR}/${BASE}_outputMafToGVCF.txt"

"$TASSEL_DIR/run_pipeline.pl" -Xmx100G -debug \
  -MAFToGVCFPlugin \
  -referenceFasta "$REF_FASTA" \
  -mafFile "$MAF_FILE" \
  -sampleName "$SAMPLE_NAME" \
  -gvcfOutput "$GVCF_OUT" \
  -fillGaps "$FILL_GAPS" \
  > "$LOG_OUT"
