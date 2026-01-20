#!/usr/bin/env bash

SPECIES_LIST="cetacean.species.txt"

# TODO: set this to the correct genome FASTA per species (see notes below).
# Option A (recommended): a TSV mapping species -> genome.fa path.
GENOME_MAP_TSV="cetacean.genomes.tsv"

BASE_ANALYSIS_DIR="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/cetaceans"
SEXCHR_DIR="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/sex_chrs"

if [[ ! -f "$SPECIES_LIST" ]]; then
  echo "ERROR: species list not found: $SPECIES_LIST" >&2
  exit 1
fi

if [[ ! -f "$GENOME_MAP_TSV" ]]; then
  echo "ERROR: genome map TSV not found: $GENOME_MAP_TSV" >&2
  echo "Create it like:  SPECIES<TAB>/path/to/genome.fa" >&2
  exit 1
fi

# load genome paths into an associative array
declare -A GENOME
while IFS=$'\t' read -r sp fa; do
  [[ -z "${sp:-}" ]] && continue
  [[ "${sp:0:1}" == "#" ]] && continue
  GENOME["$sp"]="$fa"
done < "$GENOME_MAP_TSV"

while IFS= read -r SPECIES; do
  [[ -z "${SPECIES:-}" ]] && continue
  [[ "${SPECIES:0:1}" == "#" ]] && continue

  echo "==> Processing: $SPECIES"

  GENOME_FA="${GENOME[$SPECIES]:-}"
  if [[ -z "${GENOME_FA}" ]]; then
    echo "WARNING: No genome.fa path found for $SPECIES in $GENOME_MAP_TSV; skipping." >&2
    continue
  fi
  if [[ ! -f "$GENOME_FA" ]]; then
    echo "WARNING: genome.fa not found for $SPECIES at $GENOME_FA; skipping." >&2
    continue
  fi

  SPECIES_DIR="${BASE_ANALYSIS_DIR}/${SPECIES}"
  mkdir -p "$SPECIES_DIR"
  cd "$SPECIES_DIR"

  # Ensure fasta index exists for this genome
  samtools faidx "$GENOME_FA"

  # Get X chromosome
  samtools faidx "$GENOME_FA" chrX > "${SEXCHR_DIR}/${SPECIES}.X.fa"
  bgzip -f "${SEXCHR_DIR}/${SPECIES}.X.fa"

  # Get Y chromosome
  samtools faidx "$GENOME_FA" chrY > "${SEXCHR_DIR}/${SPECIES}.Y.fa"
  bgzip -f "${SEXCHR_DIR}/${SPECIES}.Y.fa"

  # Symlinks in species analysis dir
  ln -sf "${SEXCHR_DIR}/${SPECIES}.X.fa.gz" "${SPECIES}.hap1.fasta.gz"
  ln -sf "${SEXCHR_DIR}/${SPECIES}.Y.fa.gz" "${SPECIES}.hap2.fasta.gz"

  # Index gz fasta
  samtools faidx "${SPECIES}.hap1.fasta.gz"
  samtools faidx "${SPECIES}.hap2.fasta.gz"

  # Write config file (expands species correctly)
  cat > "${SPECIES_DIR}/config.${SPECIES}.json" <<EOF
{
  "prefix": "${SPECIES}"
}
EOF

  # Submit job
  sbatch scinkd_submit.cetaceans.sh -s "$SPECIES"

done < "$SPECIES_LIST"

echo "Done."
