#!/usr/bin/env bash
#SBATCH --job-name=genomeark_annot_gff
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-1

set -euo pipefail

# -------- user inputs --------
ARK="${ARK:-ARK}"  # allow override: sbatch --export=ARK=/path/to/ARK this_script.sh
BASE="/data/Wilson_Lab/data/VGP_genomes_phase1"
WORKDIR="/data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations"
# -----------------------------

mkdir -p "${WORKDIR}/logs"
cd "${WORKDIR}"
module load agat/1.4.0

# Determine number of tasks if you didn't set it right:
#   N=$(($(wc -l < "$ARK") - 1)); sbatch --array=1-$N this_script.sh
# We assume --array is already set correctly.

# Pull the line corresponding to this array index (skip header)
LINE="$(tail -n +2 "${ARK}" | sed -n "${SLURM_ARRAY_TASK_ID}p" || true)"
if [[ -z "${LINE}" ]]; then
  echo "ERROR: No line for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} in ${ARK}" >&2
  exit 1
fi

# Parse CSV: GenomeArkID,Species
SP_ID="$(echo "${LINE}" | awk -F',' '{print $1}')"
SPECIES="$(echo "${LINE}" | awk -F',' '{print $2}')"

if [[ -z "${SP_ID}" || -z "${SPECIES}" ]]; then
  echo "ERROR: Failed to parse SP_ID/SPECIES from line: ${LINE}" >&2
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: SP_ID=${SP_ID} SPECIES=${SPECIES}"

mkdir -p "${SP_ID}"
cd "${SP_ID}"

# Download annotation files

BASE_URL="https://genomeark.s3.amazonaws.com/downstream_analyses/EGAPx_annotations/${SP_ID}_output"

gff="${SP_ID}.gff"
gtf="${SP_ID}.gtf"

have_10_lines() {
  local f="$1"
  [[ -s "$f" ]] && [[ "$(wc -l < "$f")" -ge 10 ]]
}

# 1) Try to get a .gff directly
# 1) Try to get a .gff directly (URL A, then fallback URL B)
SRC_URL_A="${BASE_URL}/complete_genomic_gff.gff"
SRC_URL_B="${BASE_URL}/complete.genomic.gff"   # fallback (example)

if ! have_10_lines "$gff"; then
  wget -q -O "${SP_ID}.tmp.gff" "$SRC_URL_A" || true

  if have_10_lines "${SP_ID}.tmp.gff"; then
    mv -f "${SP_ID}.tmp.gff" "$gff"
  else
    rm -f "${SP_ID}.tmp.gff"

    wget -q -O "${SP_ID}.tmp.gff" "$SRC_URL_B" || true
    if have_10_lines "${SP_ID}.tmp.gff"; then
      mv -f "${SP_ID}.tmp.gff" "$gff"
    else
      rm -f "${SP_ID}.tmp.gff"
    fi
  fi
fi

# 2) If .gff still missing/too small, try to get a .gtf
SRC_URL_A="${BASE_URL}/complete_genomic_gtf.gtf"
SRC_URL_B="${BASE_URL}/complete.genomic.gtf"   # fallback (example)

if ! have_10_lines "$gff"; then
  wget -q -O "${SP_ID}.tmp.gtf" "$SRC_URL_A" || true

  if have_10_lines "${SP_ID}.tmp.gtf"; then
    mv -f "${SP_ID}.tmp.gtf" "$gtf"
  else
    rm -f "${SP_ID}.tmp.gtf"

    wget -q -O "${SP_ID}.tmp.gtf" "$SRC_URL_B" || true
    if have_10_lines "${SP_ID}.tmp.gtf"; then
      mv -f "${SP_ID}.tmp.gtf" "$gtf"
    else
      rm -f "${SP_ID}.tmp.gtf"
    fi
  fi
fi