#!/bin/bash

SEX_CHR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv

noYW_DIR=/data/Wilson_Lab/data/VGP_genomes_phase1/No_YW/

# species list: "SpeciesName  input_dir"
SPEC_LIST=(
  "Homo_sapiens        /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Homo_sapiens"
  "Gallus_gallus       /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Gallus_gallus"
  "Anolis_sagrei       /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Anolis_sagrei"
  "Podarcis_raffonei   /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Podarcis_raffonei"
  "Hoplias_malabaricus /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Hoplias_malabaricus"
  "Narcine_bancroftii  /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Narcine_bancroftii"
  "Pseudacris_triseriata  /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Pseudacris_triseriata"
  "Hyla_sarda/        /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Hyla_sarda/"
)

mkdir -p "$noYW_DIR"

filter_fasta_noYW () {
  local infile="$1"
  local outfile="$2"
  local acc_re="$3"

  awk -v acc_re="$acc_re" '
    BEGIN{IGNORECASE=1}
    /^>/{
      header=$0
      drop=0
      # accession-based drop
      if (acc_re != "" && header ~ acc_re) drop=1
      # name-based fallback drop (common conventions)
      if (header ~ /(^|[[:space:];|,])(chr)?[YW]([[:space:];|,]|$)/) drop=1
      if (header ~ /(chromosome|chrom)[[:space:]]*[YW]\b/) drop=1
    }
    { if(!drop) print }
  ' "$infile" > "$outfile"
}

filter_gff_noYW () {
  local infile="$1"
  local outfile="$2"
  local acc_re="$3"

  awk -v acc_re="$acc_re" '
    BEGIN{FS=OFS="\t"; IGNORECASE=1}
    /^#/ {print; next}
    {
      seqid=$1
      drop=0
      # seqid often exactly equals an accession
      if (acc_re != "" && seqid ~ ("^(" acc_re ")$")) drop=1
      # common seqid conventions
      if (seqid ~ /^(chr)?[YW]$/) drop=1
      if (seqid ~ /(chromosome|chrom)[[:space:]]*[YW]\b/) drop=1
      if (!drop) print
    }
  ' "$infile" > "$outfile"
}

for item in "${SPEC_LIST[@]}"; do
  # split "Species  /path"
  SPECIES=$(awk '{print $1}' <<<"$item")
  INDIR=$(awk '{print $2}' <<<"$item")
  OUTDIR="$noYW_DIR/$SPECIES"

  mkdir -p "$OUTDIR"

  # build regex like: NC_060948.1|CM123.4 (may be empty)
  YW_ACC_RE=$(awk -F',' -v sp="$SPECIES" '$1==sp && ($2=="Y" || $2=="W"){print $3}' "$SEX_CHR" | paste -sd'|' -)

  echo "==> $SPECIES"
  echo "    in : $INDIR"
  echo "    out: $OUTDIR"
  echo "    Y/W accessions: ${YW_ACC_RE:-<none>}"

  # FASTA-like files
  for f in "$INDIR"/*.fna "$INDIR"/*.faa "$INDIR"/*.cds "$INDIR"/*.translated.cds; do
    [ -e "$f" ] || continue
    base=$(basename "$f")
    filter_fasta_noYW "$f" "$OUTDIR/$base" "$YW_ACC_RE"
  done

  # GFF
  shopt -s nullglob
  gff_in=( "$INDIR/${SPECIES}"*.gff )
  shopt -u nullglob 

  if [ -e "$gff_in" ]; then
    filter_gff_noYW "$gff_in" "$OUTDIR/${SPECIES}.gff" "$YW_ACC_RE"
  else
    echo "    (no GFF found at $gff_in)"
  fi

  # Rebuild FAI if possible
  if command -v samtools >/dev/null 2>&1; then
    if [ -e "$OUTDIR/${SPECIES}.fna" ]; then
      samtools faidx "$OUTDIR/${SPECIES}.fna"
    fi
  fi

done
