#!/usr/bin/env bash
#SBATCH --job-name=vgp_dl
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH --output=slurm_logs/vgp_dl_%j.out
#SBATCH --error=slurm_logs/vgp_dl_%j.err
#SBATCH --mail-type=ALL

set -euo pipefail
shopt -s nullglob

source myconda
conda activate lifton

while getopts "s:r:" opt; do
  case "$opt" in
    s) SPECIES="$OPTARG" ;;
    r) REFGENOME="$OPTARG" ;;
    *) echo "Usage: $0 -s SPECIES -r REFGENOME" >&2; exit 1 ;;
  esac
done

if [[ -z "${SPECIES:-}" || -z "${REFGENOME:-}" ]]; then
  echo "Error: both -s and -r are required" >&2
  exit 1
fi

GENOME_DIR="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes"
OUT_DIR="/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/${SPECIES}"

mkdir -p "$OUT_DIR" slurm_logs

# Find the single GC* directory for the reference
ref_data_dirs=( "${GENOME_DIR}/${REFGENOME}/ncbi_dataset/data/"GC* )
if [[ ${#ref_data_dirs[@]} -ne 1 ]]; then
  echo "Error: expected exactly one reference data dir, found ${#ref_data_dirs[@]}" >&2
  printf '%s\n' "${ref_data_dirs[@]}" >&2
  exit 1
fi
REFGENOME_DIR="${ref_data_dirs[0]}"

# Find the single GC* directory for the query
qry_data_dirs=( "${GENOME_DIR}/${SPECIES}/ncbi_dataset/data/"GC* )
if [[ ${#qry_data_dirs[@]} -ne 1 ]]; then
  echo "Error: expected exactly one query data dir, found ${#qry_data_dirs[@]}" >&2
  printf '%s\n' "${qry_data_dirs[@]}" >&2
  exit 1
fi
SPECIES_DIR="${qry_data_dirs[0]}"

if [[ -f "${REFGENOME_DIR}/genomic.gff_db" ]]; then
  GFF_REF="${REFGENOME_DIR}/genomic.gff_db"
else
  GFF_REF="${REFGENOME_DIR}/genomic.gff"
fi

FAA_REF="${REFGENOME_DIR}/protein.faa"
GFF_QRY="${OUT_DIR}/lifted.${SPECIES}.0_5.gff"

ref_files=( "${REFGENOME_DIR}/"GC*genomic.fna )
if [[ ${#ref_files[@]} -ne 1 ]]; then
  echo "Error: expected exactly one reference FASTA, found ${#ref_files[@]}" >&2
  printf '%s\n' "${ref_files[@]}" >&2
  exit 1
fi
FNA_REF="${ref_files[0]}"

qry_files=( "${SPECIES_DIR}/"GC*genomic.fna )
if [[ ${#qry_files[@]} -ne 1 ]]; then
  echo "Error: expected exactly one query FASTA, found ${#qry_files[@]}" >&2
  printf '%s\n' "${qry_files[@]}" >&2
  exit 1
fi
FNA_QRY="${qry_files[0]}"

SC=0.95

echo "Reference dir:   $REFGENOME_DIR"
echo "Reference FASTA: $FNA_REF"
echo "Query dir:       $SPECIES_DIR"
echo "Query FASTA:     $FNA_QRY"


# Use node-local scratch if available
JOB_TMP="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/lifton_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$JOB_TMP"

# Copy the prebuilt DB to a job-local path



if [[ -f "${REFGENOME_DIR}/genomic.gff_db" ]]; then
  LOCAL_GFF_DB="${JOB_TMP}/genomic.gff_db"

  cp "$GFF_REF" "$LOCAL_GFF_DB"
  
  lifton \
    -g "$LOCAL_GFF_DB" \
    -o "${GFF_QRY}" \
    -P "${FAA_REF}" \
    -t "${SLURM_CPUS_PER_TASK:-12}" \
    -copies \
    -sc "${SC}" \
    "${FNA_QRY}" \
    "${FNA_REF}"
else
  lifton \
    -g "$GFF_REF" \
    -o "${GFF_QRY}" \
    -P "${FAA_REF}" \
    -t "${SLURM_CPUS_PER_TASK:-12}" \
    -copies \
    -sc "${SC}" \
    "${FNA_QRY}" \
    "${FNA_REF}"
fi

