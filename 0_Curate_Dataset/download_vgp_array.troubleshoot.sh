#!/usr/bin/env bash
set -euo pipefail

# --- env setup ---
source myconda
conda activate ncbi_datasets

GENOME_DIR="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes"
ACCESSION_LIST="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/did_not_run.csv"

FILES_TO_DOWNLOAD="gff3,rna,cds,protein,genome,seq-report"

mkdir -p "$GENOME_DIR"
mkdir -p "$GENOME_DIR/slurm_logs"

# --- main loop ---
# Expects CSV columns: Species,Accession_GCF,Accession_GCA
# Rule: if GCF present use it (and ignore GCA); else use GCA.
tail -n +2 "$ACCESSION_LIST" | while IFS=, read -r SPECIES_NAME ACC_GCF ACC_GCA REST; do
  # trim possible CR (Windows newlines)
  SPECIES_NAME="${SPECIES_NAME//$'\r'/}"
  ACC_GCF="${ACC_GCF//$'\r'/}"
  ACC_GCA="${ACC_GCA//$'\r'/}"

  # skip empty lines
  [[ -z "${SPECIES_NAME}" ]] && continue

  # pick accession per rule
  ACCESSION=""
  if [[ -n "${ACC_GCF}" ]]; then
    ACCESSION="${ACC_GCF}"
  elif [[ -n "${ACC_GCA}" ]]; then
    ACCESSION="${ACC_GCA}"
  else
    echo "WARN: ${SPECIES_NAME} has no GCF or GCA; skipping." >&2
    continue
  fi

  OUTDIR="${GENOME_DIR}/${SPECIES_NAME}"
  mkdir -p "$OUTDIR"

  ZIPFILE="${OUTDIR}/${ACCESSION}.zip"

  echo "=== ${SPECIES_NAME} -> ${ACCESSION}"
  echo "Downloading: ${FILES_TO_DOWNLOAD}"
  datasets download genome accession "${ACCESSION}" \
    --include "${FILES_TO_DOWNLOAD}" \
    --filename "${ZIPFILE}"

  echo "Unzipping to: ${OUTDIR}"
  unzip -o -q "${ZIPFILE}" -d "${OUTDIR}"

  # optional: remove the zip after successful unzip
  rm -f "${ZIPFILE}"

  echo "Done: ${SPECIES_NAME}"
done
