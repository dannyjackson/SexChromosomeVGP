#!/usr/bin/env bash
#SBATCH --job-name=vgp_dl
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_logs/vgp_dl_%A_%a.out
#SBATCH --error=slurm_logs/vgp_dl_%A_%a.err
#SBATCH --mail-type=ALL

set -euo pipefail

source myconda
conda activate lifton

while getopts "s:a:" opt; do
  case "$opt" in
    s) SPECIES="$OPTARG" ;;
    a) ACCESSION="$OPTARG" ;;
    *) echo "Usage: $0 -s SPECIES -a ACCESSION" >&2; exit 1 ;;
  esac
done

if [[ -z "${SPECIES:-}" || -z "${ACCESSION:-}" ]]; then
  echo "Error: both -s and -a are required" >&2
  exit 1
fi

SOURCE_DIR="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Anolis_sagrei/ncbi_dataset/data/GCF_037176765.1"
FNA_REF="${SOURCE_DIR}/GCF_037176765.1_rAnoSag1.mat_genomic.fna"
FAA_REF="${SOURCE_DIR}/protein.faa"
GFF_REF_DB="${SOURCE_DIR}/genomic.gff_db"

GFF_QRY="/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_nonVGP_lifted_gffs/${SPECIES}/lifted.${ACCESSION}.0_5.gff"
FNA_QRY_GLOB="/data/Wilson_Lab/data/VGP_genomes_phase1/squamate_nonVGP_genomes/${SPECIES}/ncbi_dataset/data/${ACCESSION}/${ACCESSION}"*genomic.fna

mkdir -p "$(dirname "$GFF_QRY")" slurm_logs

qry_files=( $FNA_QRY_GLOB )
if [[ ${#qry_files[@]} -ne 1 ]]; then
  echo "Error: expected exactly one query FASTA, found ${#qry_files[@]}" >&2
  printf '%s\n' "${qry_files[@]}" >&2
  exit 1
fi
FNA_QRY="${qry_files[0]}"

SC=0.95

# Use node-local scratch if available
JOB_TMP="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/lifton_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$JOB_TMP"

# Copy the prebuilt DB to a job-local path
LOCAL_GFF_DB="${JOB_TMP}/genomic.gff_db"
cp "$GFF_REF_DB" "$LOCAL_GFF_DB"

lifton \
    -g "$LOCAL_GFF_DB" \
    -o "${GFF_QRY}" \
    -P "${FAA_REF}" \
    -t "${SLURM_CPUS_PER_TASK}" \
    -copies \
    -sc "${SC}" \
    "${FNA_QRY}" \
    "${FNA_REF}"