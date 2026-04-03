#!/usr/bin/env bash
#SBATCH --job-name=genomeark_annot_faa
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


out="${SP_ID}.cds"
url_fasta="${BASE_URL}/complete_proteins.fasta"
url_faa="${BASE_URL}/complete.proteins.faa"

have_10_lines() {
  local f="$1"
  [[ -s "$f" ]] && [[ "$(wc -l < "$f")" -ge 10 ]]
}

# Try .fasta first, then .faa
if ! have_10_lines "$out"; then
  wget -q -O "${SP_ID}.tmp.proteins" "$url_fasta" || true
  if have_10_lines "${SP_ID}.tmp.proteins"; then
    mv -f "${SP_ID}.tmp.proteins" "$out"
  else
    rm -f "${SP_ID}.tmp.proteins"
    wget -q -O "${SP_ID}.tmp.proteins" "$url_faa" || true
    if have_10_lines "${SP_ID}.tmp.proteins"; then
      mv -f "${SP_ID}.tmp.proteins" "$out"
    else
      rm -f "${SP_ID}.tmp.proteins"
    fi
  fi
fi

if have_10_lines "$out"; then
  echo "OK: using $out"
else
  echo "ERROR: could not obtain usable proteins from:"
  echo "  $url_fasta"
  echo "  $url_faa"
  exit 1
fi
