#!/usr/bin/env bash

SPECIES_LIST="cetacean.species.txt"

BASE_ANALYSIS_DIR="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/cetaceans"
SEXCHR_DIR="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/sex_chrs"
GENOME_ROOT="/data/Wilson_Lab/data/VGP_genomes"

SEXCHROM_CSV="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/sexchrom_accessions.csv"

if [[ ! -f "$SPECIES_LIST" ]]; then
  echo "ERROR: species list not found: $SPECIES_LIST" >&2
  exit 1
fi
if [[ ! -f "$SEXCHROM_CSV" ]]; then
  echo "ERROR: mapping CSV not found: $SEXCHROM_CSV" >&2
  exit 1
fi

mkdir -p "$SEXCHR_DIR"

get_acc() {
  # usage: get_acc <species> <X|Y>
  local sp="$1"
  local chrom="$2"
  awk -F',' -v sp="$sp" -v c="$chrom" '
    $1==sp && $2==c { gsub(/\r/,"",$3); print $3; exit }
  ' "$SEXCHROM_CSV"
}

while IFS= read -r SPECIES; do
  [[ -z "${SPECIES:-}" ]] && continue
  [[ "${SPECIES:0:1}" == "#" ]] && continue

  GENOME_FA="${GENOME_ROOT}/${SPECIES}/${SPECIES}.fna"
  if [[ ! -f "$GENOME_FA" ]]; then
    echo "WARNING: genome not found for $SPECIES at $GENOME_FA; skipping." >&2
    continue
  fi

  X_ACC="$(get_acc "$SPECIES" "X")"
  Y_ACC="$(get_acc "$SPECIES" "Y")"

  if [[ -z "${X_ACC:-}" || -z "${Y_ACC:-}" ]]; then
    echo "WARNING: missing X/Y accession in $SEXCHROM_CSV for $SPECIES (X='${X_ACC:-}', Y='${Y_ACC:-}'); skipping." >&2
    continue
  fi

  echo "==> Processing: $SPECIES  (X=$X_ACC, Y=$Y_ACC)"

  SPECIES_DIR="${BASE_ANALYSIS_DIR}/${SPECIES}"
  mkdir -p "$SPECIES_DIR"
  cd "$SPECIES_DIR"

  # Ensure genome is indexed
  samtools faidx "$GENOME_FA"

  # Extract X and Y by accession/contig name from CSV
  samtools faidx "$GENOME_FA" "$X_ACC" > "${SEXCHR_DIR}/${SPECIES}.X.fa"
  bgzip -f "${SEXCHR_DIR}/${SPECIES}.X.fa"

  samtools faidx "$GENOME_FA" "$Y_ACC" > "${SEXCHR_DIR}/${SPECIES}.Y.fa"
  bgzip -f "${SEXCHR_DIR}/${SPECIES}.Y.fa"

  # Symlinks to hap FASTAs
  ln -sf "${SEXCHR_DIR}/${SPECIES}.X.fa.gz" "${SPECIES}.hap1.fasta.gz"
  ln -sf "${SEXCHR_DIR}/${SPECIES}.Y.fa.gz" "${SPECIES}.hap2.fasta.gz"

  # Index gz FASTAs
  samtools faidx "${SPECIES}.hap1.fasta.gz"
  samtools faidx "${SPECIES}.hap2.fasta.gz"

  # Write config
  cat > "${SPECIES_DIR}/config.${SPECIES}.json" <<EOF
{
  "prefix": "${SPECIES}"
}
EOF

done < "$SPECIES_LIST"

echo "Done."
