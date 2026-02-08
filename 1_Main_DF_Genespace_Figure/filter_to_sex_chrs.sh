#!/usr/bin/env bash

# ---- inputs ----
SEX_CHR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv

# Root where the original per-species folders live (edit if different)
IN_ROOT=/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks

# Outputs
ONLY_YW_ROOT="${IN_ROOT}/Only_YW"
ONLY_XZ_ROOT="${IN_ROOT}/Only_XZ"

mkdir -p "$ONLY_YW_ROOT" "$ONLY_XZ_ROOT"

# If you only want to process a subset, fill this list with species names (e.g. Homo_sapiens).
# Leave empty to auto-detect species from $SEX_CHR.
SPECIES_FILTER=(
  "Homo_sapiens"
  "Gallus_gallus"
  "Anolis_sagrei"
  "Podarcis_raffonei"
  "Hoplias_malabaricus"
  "Narcine_bancroftii"
)

# ---- helpers ----
get_acc_re () {
  local species="$1"
  local chr_re="$2"   # e.g. '^(Y|Y[0-9]+|W|W[0-9]+)$'

  awk -F',' -v sp="$species" -v chr_re="$chr_re" '
    BEGIN{IGNORECASE=1}
    $1==sp && $2 ~ chr_re {print $3}
  ' "$SEX_CHR" | paste -sd'|' -
}

# Keep FASTA records that match accession OR obvious name patterns (chrY, Y1, etc.)
keep_fasta () {
  local infile="$1"
  local outfile="$2"
  local acc_re="$3"
  local name_re="$4"

  awk -v acc_re="$acc_re" -v name_re="$name_re" '
    BEGIN{IGNORECASE=1; keep=0}
    /^>/{
      h=$0
      keep=0
      if (acc_re != "" && h ~ acc_re) keep=1
      if (h ~ name_re) keep=1
    }
    { if(keep) print }
  ' "$infile" > "$outfile"
}

# Keep GFF features whose seqid matches accession OR name patterns
keep_gff () {
  local infile="$1"
  local outfile="$2"
  local acc_re="$3"
  local name_re="$4"

  awk -v acc_re="$acc_re" -v name_re="$name_re" '
    BEGIN{FS=OFS="\t"; IGNORECASE=1}
    /^#/ {print; next}
    {
      seqid=$1
      keep=0
      if (acc_re != "" && seqid ~ ("^(" acc_re ")$")) keep=1
      if (seqid ~ name_re) keep=1
      if (keep) print
    }
  ' "$infile" > "$outfile"
}

# ---- regexes for name-based matching (covers Y, Y1.., W, W1.. and X, X1.., Z, Z1..) ----
# Matches tokens like "Y", "chrY", "Y1", "chrY2", "W", "chrW3", and also "chromosome Y1" etc.
NAME_RE_YW='(^|[^A-Za-z0-9])(chr)?(Y([0-9]+)?|W([0-9]+)?)([^A-Za-z0-9]|$)|chromosome[[:space:]]*(Y([0-9]+)?|W([0-9]+)?)|chrom[[:space:]]*(Y([0-9]+)?|W([0-9]+)?)'
NAME_RE_XZ='(^|[^A-Za-z0-9])(chr)?(X([0-9]+)?|Z([0-9]+)?)([^A-Za-z0-9]|$)|chromosome[[:space:]]*(X([0-9]+)?|Z([0-9]+)?)|chrom[[:space:]]*(X([0-9]+)?|Z([0-9]+)?)'

CHR_RE_YW='^(Y|Y[0-9]+|W|W[0-9]+)$'
CHR_RE_XZ='^(X|X[0-9]+|Z|Z[0-9]+)$'

# ---- determine species list ----
if ((${#SPECIES_FILTER[@]} > 0)); then
  SPECIES_LIST=("${SPECIES_FILTER[@]}")
else
  mapfile -t SPECIES_LIST < <(awk -F',' 'NR>1{print $1}' "$SEX_CHR" | sort -u)
fi

echo "Processing ${#SPECIES_LIST[@]} species..."

for SPECIES in "${SPECIES_LIST[@]}"; do
  INDIR="${IN_ROOT}/${SPECIES}"
  [ -d "$INDIR" ] || { echo "Skipping $SPECIES (no dir: $INDIR)"; continue; }

  OUT_YW="${ONLY_YW_ROOT}/${SPECIES}"
  OUT_XZ="${ONLY_XZ_ROOT}/${SPECIES}"
  mkdir -p "$OUT_YW" "$OUT_XZ"

  YW_ACC_RE=$(get_acc_re "$SPECIES" "$CHR_RE_YW")
  XZ_ACC_RE=$(get_acc_re "$SPECIES" "$CHR_RE_XZ")

  echo "==> $SPECIES"
  echo "    IN:  $INDIR"
  echo "    Y/W accessions: ${YW_ACC_RE:-<none>}"
  echo "    X/Z accessions: ${XZ_ACC_RE:-<none>}"

  # FASTA-like
  for f in "$INDIR"/*.fna "$INDIR"/*.faa "$INDIR"/*.cds "$INDIR"/*.translated.cds; do
    [ -e "$f" ] || continue
    base=$(basename "$f")

    # Only Y/W
    keep_fasta "$f" "$OUT_YW/$base" "$YW_ACC_RE" "$NAME_RE_YW"

    # Only X/Z
    keep_fasta "$f" "$OUT_XZ/$base" "$XZ_ACC_RE" "$NAME_RE_XZ"
  done

  # GFF
  gff_in="$INDIR/${SPECIES}.gff"
  if [ -e "$gff_in" ]; then
    keep_gff "$gff_in" "$OUT_YW/${SPECIES}.gff" "$YW_ACC_RE" "$NAME_RE_YW"
    keep_gff "$gff_in" "$OUT_XZ/${SPECIES}.gff" "$XZ_ACC_RE" "$NAME_RE_XZ"
  fi

  # FAI rebuild for fna outputs (if present)
  if command -v samtools >/dev/null 2>&1; then
    [ -e "$OUT_YW/${SPECIES}.fna" ] && samtools faidx "$OUT_YW/${SPECIES}.fna"
    [ -e "$OUT_XZ/${SPECIES}.fna" ] && samtools faidx "$OUT_XZ/${SPECIES}.fna"
  fi
done

echo "Done."
echo "Only Y/W outputs: $ONLY_YW_ROOT"
echo "Only X/Z outputs: $ONLY_XZ_ROOT"
