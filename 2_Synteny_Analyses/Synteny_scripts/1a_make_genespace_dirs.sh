#!/bin/bash

REF="${1:?Usage: $0 <species>}"

OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${REF}

mkdir -p ${OUTDIR}/sexshared

SEX_CHR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv"

SS_BED_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed"
SS_PEP_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide"

REF_BED="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/${REF}.bed"
REF_PEP="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/${REF}.fa"

tail -n +2 "$SEX_CHR" | while IFS=',' read -r SPECIES CHR_TYPE ACC; do
  ACC="${ACC%$'\r'}"

  # Make simple per-species structure
  mkdir -p "${OUTDIR}/sexshared/${SPECIES}/bed" "${OUTDIR}/sexshared/${SPECIES}/peptide"

  # Symlink species inputs
  ln -sf "${SS_BED_DIR}/${SPECIES}.bed" "${OUTDIR}/sexshared/${SPECIES}/bed/${SPECIES}.bed"
  ln -sf "${SS_PEP_DIR}/${SPECIES}.fa"  "${OUTDIR}/sexshared/${SPECIES}/peptide/${SPECIES}.fa"

  # Symlink references (handy to keep alongside each species)
  ln -sf "${SS_BED_DIR}/${REF}.bed" "${OUTDIR}/sexshared/${SPECIES}/bed/${REF}.bed"
  ln -sf "${SS_PEP_DIR}/${REF}.fa" "${OUTDIR}/sexshared/${SPECIES}/peptide/${REF}.fa"
done
