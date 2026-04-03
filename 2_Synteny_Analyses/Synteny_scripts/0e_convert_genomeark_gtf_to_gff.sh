#!/usr/bin/env bash
#SBATCH --job-name=genomeark_annot_gff
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
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

# If no usable .gff but a usable .gtf exists, convert gtf -> gff
if ! have_10_lines "$gff" && have_10_lines "$gtf"; then
  agat_convert_sp_gxf2gxf.pl -g "$gtf" -o "$gff.tmp"
  mv "$gff.tmp" "$gff"
fi

# 4) Final sanity check
if have_10_lines "$gff"; then
  echo "OK: using $gff"
elif have_10_lines "$gtf"; then
  echo "WARN: only $gtf is usable; conversion may have failed"
  exit 1
else
  echo "ERROR: could not obtain usable $gff or $gtf from: $SRC_URL"
  exit 1
fi