#!/bin/bash
#SBATCH --job-name=maf_to_gvcf
#SBATCH --array=0-0
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

MAF_DIR=""
OUT_DIR=""
ACCOUNT=""
PARTITION=""
MAX_CONCURRENT=""

while getopts "m:o:A:P:J:" opt; do
  case "$opt" in
    m) MAF_DIR="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    A) ACCOUNT="$OPTARG" ;;
    P) PARTITION="$OPTARG" ;;
    J) MAX_CONCURRENT="$OPTARG" ;;
    *) ;;
  esac
done

if [ -z "$MAF_DIR" ] || [ -z "$OUT_DIR" ]; then
  echo "Usage: $0 -m <maf_dir> -o <out_dir> [-A <account>] [-P <partition>] [-J <max_concurrent>]"
  exit 1
fi

# Determine repository root: prefer the SLURM submit directory when available.
# When SLURM runs a job, the job script may be executed from a spool directory
# where the user cannot create new directories. Use `SLURM_SUBMIT_DIR` so
# logs are created next to the repository (the submission location).
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  REPO_ROOT="$SLURM_SUBMIT_DIR"
else
  REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
fi

TASSEL_DIR="$REPO_ROOT/tassel-5-standalone"
LOG_DIR="$REPO_ROOT/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"
FILL_GAPS="false"
SAMPLE_SUFFIX="_anchorwave"
CONFIG_PATH="$REPO_ROOT/config.yaml"

if [ -z "$MAX_CONCURRENT" ] && [ -f "$CONFIG_PATH" ]; then
  MAX_CONCURRENT="$(awk -F: '
    /^[[:space:]]*maf_to_gvcf_array_max_jobs:/ {
      val=$2
      sub(/#.*/, "", val)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", val)
      print val
      exit
    }' "$CONFIG_PATH")"
fi
if [ -z "$MAX_CONCURRENT" ]; then
  MAX_CONCURRENT=4
fi

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

# Build list of MAFs for the array.
MAF_FILES=()
while IFS= read -r line; do
  MAF_FILES+=("$line")
done < <(find "$MAF_DIR" -maxdepth 1 -type f \( -name "*.maf" -o -name "*.maf.gz" \) | sort)
if [ "${#MAF_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no .maf or .maf.gz files found in $MAF_DIR"
  exit 1
fi

# If not running as a SLURM array task, submit ourselves with an array sized to the MAF count.
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
  N="${#MAF_FILES[@]}"
  ARRAY="0-$((N-1))%${MAX_CONCURRENT}"
  echo "Submitting SLURM array with $N tasks (max ${MAX_CONCURRENT} running at a time): $ARRAY"
  SBATCH_ARGS=(--array="${ARRAY}" --output "${LOG_DIR}/%x_%A_%a.out" --error "${LOG_DIR}/%x_%A_%a.err")
  if [ -n "$ACCOUNT" ]; then
    SBATCH_ARGS+=(--account "$ACCOUNT")
  fi
  if [ -n "$PARTITION" ]; then
    SBATCH_ARGS+=(--partition "$PARTITION")
  fi
  SCRIPT_ARGS=(-m "$MAF_DIR" -o "$OUT_DIR")
  if [ -n "$ACCOUNT" ]; then
    SCRIPT_ARGS+=(-A "$ACCOUNT")
  fi
  if [ -n "$PARTITION" ]; then
    SCRIPT_ARGS+=(-P "$PARTITION")
  fi
  sbatch "${SBATCH_ARGS[@]}" "$0" "${SCRIPT_ARGS[@]}"
  exit 0
fi

MAF_FILE="${MAF_FILES[$SLURM_ARRAY_TASK_ID]}"
BASE="$(basename "$MAF_FILE")"
BASE="${BASE%.maf.gz}"
BASE="${BASE%.maf}"
SAMPLE_NAME="${BASE}${SAMPLE_SUFFIX}"
REF_BASE="$(basename "$REF_FASTA")"
REF_BASE="${REF_BASE%.fa}"
REF_BASE="${REF_BASE%.fasta}"

GVCF_OUT="${OUT_DIR}/${BASE}To${REF_BASE}.gvcf"
LOG_OUT="${LOG_DIR}/${BASE}_outputMafToGVCF.txt"

TMP_MAF=""
MAF_INPUT="$MAF_FILE"
if [[ "$MAF_FILE" == *.gz ]]; then
  TMP_MAF="$(mktemp "${TMPDIR:-/tmp}/maf_to_gvcf.XXXXXX.maf")"
  gunzip -c "$MAF_FILE" > "$TMP_MAF"
  MAF_INPUT="$TMP_MAF"
fi

# Run tassel MAF→GVCF for this array shard.
"$TASSEL_DIR/run_pipeline.pl" -Xmx230G -debug \
  -MAFToGVCFPlugin \
  -referenceFasta "$REF_FASTA" \
  -mafFile "$MAF_INPUT" \
  -sampleName "$SAMPLE_NAME" \
  -gvcfOutput "$GVCF_OUT" \
  -fillGaps "$FILL_GAPS" \
  > "$LOG_OUT"

if [ -n "$TMP_MAF" ] && [ -f "$TMP_MAF" ]; then
  rm -f "$TMP_MAF"
fi
