#!/usr/bin/env bash
#SBATCH --job-name=vgp_download
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=slurm_logs/vgp_download_%A_%a.out
#SBATCH --error=slurm_logs/vgp_download_%A_%a.err

set -euo pipefail

# ---- required array guard ----
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: This script must be submitted as a SLURM array."
  echo "Example:"
  echo '  N=$(( $(wc -l < /path/to/species_latest_accessions.csv) - 1 ))'
  echo '  sbatch --array=1-"$N"%25 download_vgp_array.sh'
  exit 1
fi

# ---- env setup ----
source myconda
conda activate ncbi_datasets

GENOME_DIR="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes"
ACCESSION_LIST="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species_latest_accessions.csv"
FILES_TO_DOWNLOAD="gff3,rna,cds,protein,genome,seq-report"

mkdir -p "$GENOME_DIR"
mkdir -p "$GENOME_DIR/slurm_logs"

# ---- fetch the (headerless) line corresponding to this array index ----
# array task 1 => CSV line 2 (skip header)
LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))
LINE="$(sed -n "${LINE_NUM}p" "$ACCESSION_LIST" || true)"

if [[ -z "${LINE}" ]]; then
  echo "ERROR: No line found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} (expected line ${LINE_NUM})."
  exit 1
fi

# ---- parse CSV line (expects: Species,Accession_GCF,Accession_GCA) ----
IFS=, read -r SPECIES_NAME ACC_GCF ACC_GCA _REST <<< "$LINE"

# trim possible CR (Windows newlines)
SPECIES_NAME="${SPECIES_NAME//$'\r'/}"
ACC_GCF="${ACC_GCF//$'\r'/}"
ACC_GCA="${ACC_GCA//$'\r'/}"

if [[ -z "${SPECIES_NAME}" ]]; then
  echo "ERROR: Empty species name on line ${LINE_NUM}"
  exit 1
fi

# ---- pick accession: prefer GCF; else GCA; else error ----
ACCESSION=""
if [[ -n "${ACC_GCF}" ]]; then
  ACCESSION="${ACC_GCF}"
elif [[ -n "${ACC_GCA}" ]]; then
  ACCESSION="${ACC_GCA}"
else
  echo "WARN: ${SPECIES_NAME} has no GCF or GCA; skipping."
  exit 0
fi

OUTDIR="${GENOME_DIR}/${SPECIES_NAME}"
mkdir -p "$OUTDIR"

ZIPFILE="${OUTDIR}/${ACCESSION}.zip"

echo "=== Task ${SLURM_ARRAY_TASK_ID}: ${SPECIES_NAME} -> ${ACCESSION}"
echo "Downloading: ${FILES_TO_DOWNLOAD}"
datasets download genome accession "${ACCESSION}" \
  --include "${FILES_TO_DOWNLOAD}" \
  --filename "${ZIPFILE}"

echo "Unzipping to: ${OUTDIR}"
unzip -o -q "${ZIPFILE}" -d "${OUTDIR}"
rm -f "${ZIPFILE}"

echo "Done: ${SPECIES_NAME}"
