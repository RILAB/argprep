#!/bin/bash
#SBATCH --job-name=gatk_merge_gvcf
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err

set -euo pipefail

# Validate environment modules and load required tools.
usage() {
  echo "Usage: $0 -g <gvcf_dir> -r <reference_fasta> [-l interval] [-c cutoff]"
  echo "  -g  Directory containing per-sample .gvcf.gz files"
  echo "  -r  Reference FASTA (used for GATK and indexing)"
  echo "  -l  Optional interval (e.g., chr1) passed to GenomicsDBImport/GenotypeGVCFs"
  echo "  -c  Optional indel length cutoff passed to dropSV.py (default: dropSV.py default)"
  exit 1
}

LOG_DIR="${SLURM_SUBMIT_DIR:-.}/logs"
mkdir -p "$LOG_DIR"

module load conda || { echo "ERROR: failed to load conda module"; exit 1; }
conda activate argprep || { echo "ERROR: failed to activate conda env: argprep"; exit 1; }


echo "Tool versions:"
picard MarkDuplicates --version 2>&1 | head -n 1 || true
tabix 2>&1 | head -n 3 | tail -1 || true
gatk --version 2>&1 | head -n 1 || true

# Parse inputs for gVCF dir, reference, optional interval, optional cutoff.
GVCF_DIR=""
REF_FASTA=""
INTERVAL=""
DROP_CUTOFF=""

while getopts "g:r:l:c:" opt; do
  case "$opt" in
    g) GVCF_DIR="$OPTARG" ;;
    r) REF_FASTA="$OPTARG" ;;
    l) INTERVAL="$OPTARG" ;;
    c) DROP_CUTOFF="$OPTARG" ;;
    *) usage ;;
  esac
done

if [ -z "$GVCF_DIR" ] || [ -z "$REF_FASTA" ]; then
  usage
fi

# Prefer the submit directory so we can find repo files even when SLURM runs from a spool dir.
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  REPO_ROOT="$SLURM_SUBMIT_DIR"
else
  REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
fi
SCRIPT_DIR="$REPO_ROOT/scripts"
DROP_SV="${SCRIPT_DIR}/dropSV.py"

if [ ! -x "$DROP_SV" ]; then
  echo "ERROR: dropSV.py not found or not executable at $DROP_SV"
  exit 1
fi

if [ ! -f "$REF_FASTA" ]; then
  echo "ERROR: reference FASTA not found: $REF_FASTA"
  exit 1
fi

# Required external tools (after conda activation).
# command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH"; exit 1; }
# command -v gatk >/dev/null 2>&1 || { echo "ERROR: gatk not found in PATH"; exit 1; }
# command -v bgzip >/dev/null 2>&1 || { echo "ERROR: bgzip not found in PATH"; exit 1; }
# command -v tabix >/dev/null 2>&1 || { echo "ERROR: tabix not found in PATH"; exit 1; }
# command -v picard >/dev/null 2>&1 || { echo "ERROR: picard not found in PATH"; exit 1; }


# Ensure reference index and dictionary exist.
if [ ! -f "${REF_FASTA}.fai" ]; then
  samtools faidx "$REF_FASTA"
fi

DICT="${REF_FASTA%.*}.dict"
if [ ! -f "$DICT" ]; then
  picard CreateSequenceDictionary R="$REF_FASTA" O="$DICT"
fi


# Run dropSV.py on the input directory (optional cutoff override).
if [ -n "$DROP_CUTOFF" ]; then
  python3 "$DROP_SV" -d "$GVCF_DIR" -c "$DROP_CUTOFF"
else
  python3 "$DROP_SV" -d "$GVCF_DIR"
fi

CLEAN_DIR="${GVCF_DIR%/}/cleangVCF"
if [ ! -d "$CLEAN_DIR" ]; then
  echo "ERROR: expected cleaned gVCFs in $CLEAN_DIR"
  exit 1
fi

# Ensure all cleaned gVCFs are bgzipped + indexed.
for f in "$CLEAN_DIR"/*.gvcf; do
  [ -e "$f" ] || continue
  bgzip -c "$f" > "${f}.gz"
  tabix -p vcf "${f}.gz"
done

for f in "$CLEAN_DIR"/*.gvcf.gz; do
  [ -e "$f" ] || continue
  if [ ! -f "${f}.tbi" ]; then
    tabix -p vcf "$f"
  fi
done

GVCF_FILES=()
while IFS= read -r line; do
  GVCF_FILES+=("$line")
done < <(find "$CLEAN_DIR" -maxdepth 1 -type f -name "*.gvcf.gz" | sort)
if [ "${#GVCF_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no .gvcf.gz files found in $CLEAN_DIR"
  exit 1
fi

# If an interval is provided (e.g., chr1 or chr1:1-1000), only process that contig.
INTERVAL_CONTIG=""
if [ -n "$INTERVAL" ]; then
  INTERVAL_CONTIG="${INTERVAL%%:*}"
fi

SPLIT_DIR="${CLEAN_DIR%/}/split_gvcf"
mkdir -p "$SPLIT_DIR"

declare -A CONTIG_SET=()

for f in "${GVCF_FILES[@]}"; do
  CONTIGS=()
  while IFS= read -r line; do
    CONTIGS+=("$line")
  done < <(tabix -l "$f")
  if [ "${#CONTIGS[@]}" -eq 0 ]; then
    echo "ERROR: no contigs found in $f (missing or corrupt index?)"
    exit 1
  fi

  if [ -n "$INTERVAL_CONTIG" ]; then
    WANT_CONTIGS=()
    for c in "${CONTIGS[@]}"; do
      if [ "$c" = "$INTERVAL_CONTIG" ]; then
        WANT_CONTIGS+=("$c")
        break
      fi
    done
  else
    WANT_CONTIGS=("${CONTIGS[@]}")
  fi

  if [ "${#WANT_CONTIGS[@]}" -eq 0 ]; then
    continue
  fi

  base="$(basename "$f" .gvcf.gz)"
  if [ "${#CONTIGS[@]}" -eq 1 ]; then
    contig="${CONTIGS[0]}"
    if [ -n "$INTERVAL_CONTIG" ] && [ "$contig" != "$INTERVAL_CONTIG" ]; then
      continue
    fi
    ln -sf "$f" "$SPLIT_DIR/${base}.${contig}.gvcf.gz"
    ln -sf "${f}.tbi" "$SPLIT_DIR/${base}.${contig}.gvcf.gz.tbi"
    CONTIG_SET["$contig"]=1
    continue
  fi

  for contig in "${WANT_CONTIGS[@]}"; do
    out="${SPLIT_DIR}/${base}.${contig}.gvcf.gz"
    if [ ! -f "$out" ]; then
      gatk --java-options "-Xmx100g -Xms100g" SelectVariants \
        -R "$REF_FASTA" -V "$f" -L "$contig" -O "$out"
      tabix -p vcf "$out"
    fi
    CONTIG_SET["$contig"]=1
  done
done

if [ "${#CONTIG_SET[@]}" -eq 0 ]; then
  echo "ERROR: no contigs to process after splitting."
  exit 1
fi

for contig in "${!CONTIG_SET[@]}"; do
  CONTIG_FILES=()
  while IFS= read -r line; do
    CONTIG_FILES+=("$line")
  done < <(find "$SPLIT_DIR" -maxdepth 1 -type f -name "*.${contig}.gvcf.gz" | sort)
  if [ "${#CONTIG_FILES[@]}" -eq 0 ]; then
    continue
  fi

  GVCF_ARGS=()
  for f in "${CONTIG_FILES[@]}"; do
    GVCF_ARGS+=( -V "$f" )
  done

  WORKSPACE="genomicsdb_workspace_${contig}"
  OUT_GVCF="${PWD}/combined.${contig}.gvcf.gz"

  if [ -d "$WORKSPACE" ]; then
    echo "ERROR: $WORKSPACE already exists; remove it before re-running."
    exit 1
  fi

  IMPORT_CMD=(gatk --java-options "-Xmx100g -Xms100g" GenomicsDBImport "${GVCF_ARGS[@]}" \
    --genomicsdb-workspace-path "$WORKSPACE")
  GENO_CMD=(gatk --java-options "-Xmx100g -Xms100g" GenotypeGVCFs \
    -R "$REF_FASTA" -V "gendb://${WORKSPACE}" -O "$OUT_GVCF")

  if [ -n "$INTERVAL" ]; then
    IMPORT_CMD+=( -L "$INTERVAL" )
    GENO_CMD+=( -L "$INTERVAL" )
  else
    IMPORT_CMD+=( -L "$contig" )
    GENO_CMD+=( -L "$contig" )
  fi

  # Build GenomicsDB and emit merged gVCF for this contig.
  "${IMPORT_CMD[@]}"
  "${GENO_CMD[@]}"

  echo "Combined gVCF written to $OUT_GVCF"
done
