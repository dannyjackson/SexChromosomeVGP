#!/usr/bin/env bash
#SBATCH -J sexchr_split
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-330

set -euo pipefail

mkdir -p logs

module load samtools

# ---- inputs ----
SEX_CHR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv

# Root where the original per-species folders live
IN_ROOT=/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks

# Outputs
ONLY_YW_ROOT="${IN_ROOT}/Only_YW"
ONLY_XZ_ROOT="${IN_ROOT}/Only_XZ"

mkdir -p "$ONLY_YW_ROOT" "$ONLY_XZ_ROOT"

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

# ---- regexes for name-based matching ----
NAME_RE_YW='(^|[^A-Za-z0-9])(chr)?(Y([0-9]+)?|W([0-9]+)?)([^A-Za-z0-9]|$)|chromosome[[:space:]]*(Y([0-9]+)?|W([0-9]+)?)|chrom[[:space:]]*(Y([0-9]+)?|W([0-9]+)?)'
NAME_RE_XZ='(^|[^A-Za-z0-9])(chr)?(X([0-9]+)?|Z([0-9]+)?)([^A-Za-z0-9]|$)|chromosome[[:space:]]*(X([0-9]+)?|Z([0-9]+)?)|chrom[[:space:]]*(X([0-9]+)?|Z([0-9]+)?)'

CHR_RE_YW='^(Y|Y[0-9]+|W|W[0-9]+)$'
CHR_RE_XZ='^(X|X[0-9]+|Z|Z[0-9]+)$'

# ---- build species list, pick one by array index ----
# NOTE: SLURM_ARRAY_TASK_ID is 1-based by default in the --array range below.
mapfile -t SPECIES_LIST < <(awk -F',' 'NR>1{print $1}' "$SEX_CHR" | sort -u)

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
IDX=$((TASK_ID - 1))

if (( IDX < 0 || IDX >= ${#SPECIES_LIST[@]} )); then
  echo "ERROR: SLURM_ARRAY_TASK_ID=$TASK_ID out of range (1..${#SPECIES_LIST[@]})"
  exit 1
fi

SPECIES="${SPECIES_LIST[$IDX]}"
echo "Array task ${TASK_ID}/${#SPECIES_LIST[@]} -> SPECIES=$SPECIES"

# ---- per-species processing ----
INDIR="${IN_ROOT}/${SPECIES}"
if [ ! -d "$INDIR" ]; then
  echo "Skipping $SPECIES (no dir: $INDIR)"
  exit 0
fi

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

echo "Done."
echo "Only Y/W outputs: $ONLY_YW_ROOT"
echo "Only X/Z outputs: $ONLY_XZ_ROOT"
