# Panthera onca synteny
https://genomeark.s3.amazonaws.com/index.html?prefix=downstream_analyses/EGAPx_annotations/mPanOnc1.hap2_output/

# Remove old genome (hap1, lacks sex chrs)
```
cd /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Panthera_onca
rm -r *
```
# Download hap2
```
sinteractive
source myconda
conda activate ncbi_datasets

FILES_TO_DOWNLOAD="gff3,rna,cds,protein,genome,seq-report"
ACCESSION="GCA_046562875.2"
datasets download genome accession "${ACCESSION}" \
  --include "${FILES_TO_DOWNLOAD}" 

unzip ncbi_dataset.zip
rm ncbi_dataset.zip
```
# Download annotation file
```
cd /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations
mkdir mPanOnc1
cd mPanOnc1
wget https://genomeark.s3.amazonaws.com/downstream_analyses/EGAPx_annotations/mPanOnc1.hap2_output/complete.proteins.faa
wget https://genomeark.s3.amazonaws.com/downstream_analyses/EGAPx_annotations/mPanOnc1.hap2_output/complete.genomic.gff

mv complete.proteins.faa mPanOnc1.cds
mv complete.genomic.gff mPanOnc1.gff
```
# Make symlink for FASTA
```
ln -sf /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Panthera_onca/ncbi_dataset/data/GCA_046562875.2/GCA_046562875.2_mPanOnc1_haplotype_2_genomic.fna mPanOnc1.fa 
```
# Make bed and peptide files
```
library(GENESPACE)

SPECIES <- "mPanOnc1"
outdir <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark/"
repo <- "/data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/"

parsedPaths2 <- parse_annotations(
  rawGenomeRepo = repo,
  genomeDirs    = SPECIES,
  genomeIDs     = SPECIES,
  gffString     = "gff",
  faString      = "cds",
  presets       = "none",
  genespaceWd   = outdir,
  gffIdColumn        = "protein_id",
  headerSep          = " ",
  headerEntryIndex   = 1,
  headerStripText    = ".*\\|",   # strip everything up to last '|'
  gffStripText       = ".*\\|"    # same normalization on protein_id values
)
```
## Rename all bed and peptide files to species names, not IDs
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark/bed
mv mPanOnc1.bed Panthera_onca.bed
cd ../peptide
mv mPanOnc1.fa Panthera_onca.fa
```
# Reformat bed files to have chromosome numbers, not accession numbers
## Create a reformating index
```
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark
FNA="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Panthera_onca/ncbi_dataset/data/GCA_046562875.2/GCA_046562875.2_mPanOnc1_haplotype_2_genomic.fna"
MAPDIR="${OUTDIR}/reformatting_chr"
mkdir -p "$MAPDIR"

SPECIES="Panthera_onca"
MAP="${MAPDIR}/${SPECIES}.chr_reformat.txt"

[[ ! -s "$FNA" ]] && continue

awk '
  BEGIN{ OFS="," }
  function trim(s){ gsub(/^[[:space:]]+|[[:space:]]+$/, "", s); return s }
  function up(s){ return toupper(s) }

  /^>/{
    line=$0
    gsub(/\r/,"",line)

    acc=$1; sub(/^>/,"",acc); acc=trim(acc)
    u=up(line)

    # same exclusion as above (optional but recommended)
    if (u ~ /(SCAFFOLD|CONTIG|UNPLACED|UNLOCALIZED|UNLOC|PATCH|ALT_REF_LOCI|ALTERNATE|HAPLOTYP)/) next

    chr=""
    if (match(u, /CHROMOSOME[:[:space:]]+([A-Z][0-9]*|[0-9]+)/, m)) chr=m[1]
    else if (match(u, /(^|[^A-Z0-9])CHR[[:space:]]*([A-Z][0-9]*|[0-9]+)/, m)) chr=m[2]
    else if (match(u, /LINKAGE[[:space:]]+GROUP[:[:space:]]+([0-9]+)/, m)) chr=m[1]
    else if (match(u, /(^|[^A-Z0-9])LG[[:space:]]*([0-9]+)/, m)) chr=m[2]

    if (chr!="") print acc, chr
  }
' "$FNA" | sort -u > "$MAP"

[[ -s "$MAP" ]] && echo "Wrote $MAP" >&2

grep 'CM' /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark/reformatting_chr/Panthera_onca.chr_reformat.txt > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark/reformatting_chr/Panthera_onca.chr_reformat.txt.tmp
mv /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark/reformatting_chr/Panthera_onca.chr_reformat.txt.tmp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark/reformatting_chr/Panthera_onca.chr_reformat.txt
```
# Output list of instances in which chrs in reformating index do not match any chrs in bed
This suggests a mismatch between versions of genome between annotation and assembly
```

SPECIES="Panthera_onca"

OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark
MAP="${MAPDIR}/${SPECIES}.chr_reformat.txt"
INBED_DIR="${OUTDIR}/bed"
OUTBED_DIR="${INBED_DIR}_chrfixed"
BED="${INBED_DIR}/${SPECIES}.bed"
OUT="${OUTBED_DIR}/${SPECIES}.bed"
MISSING_REPORT=missingreport.txt

# ---- 1) Report: chromosomes present in MAP (keys) that do NOT occur in BED col1 ----
# Output columns: species<TAB>map_chr
awk -v FS="\t" -v OFS="\t" -v sp="${SPECIES}" '
  # First pass: BED -> collect chromosomes present
  FNR==NR { bed[$1]=1; next }

  # Second pass: MAP -> parse "old,new" and report old not in BED
  {
    split($0,a,",");
    old=a[1];
    if (!(old in bed)) print sp, old;
  }
' "${BED}" "${MAP}" >> "${MISSING_REPORT}"

# (Optional) also echo a per-species summary to stderr
n_missing=$(awk -v FS="\t" -v sp="${SPECIES}" '$1==sp{c++} END{print c+0}' "${MISSING_REPORT}")
if (( n_missing > 0 )); then
  echo "[${SPECIES}] MAP keys missing from BED: ${n_missing} (see ${MISSING_REPORT})" >&2
fi

# ---- 2) Do the chr reformat on BED using MAP ----
awk -v FS="\t" -v OFS="\t" '
  FNR==NR { split($0,a,","); map[a[1]]=a[2]; next }
  ($1 in map) { $1=map[$1] }
  { print }
' "${MAP}" "${BED}" > "${OUT}"

[[ -s "${OUT}" ]] && echo "Wrote ${OUT}" >&2

echo "Missing-chrom report written to: ${MISSING_REPORT}" >&2
```
## Remove W and Y chromosomes from bed files, and all scaffolds
```
INDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/


INBED="${INDIR}/bed_chrfixed"
INPEP="${INDIR}/peptide"

OUTBED="${OUTDIR}/sexshared/usable_bed"
OUTPEP="${OUTDIR}/sexshared/usable_peptide"

bed_in="${INBED}/Panthera_onca.bed"
bed_out="${OUTBED}/Panthera_onca.bed"

pep_in="${INPEP}/Panthera_onca.fa"
pep_out="${OUTPEP}/Panthera_onca.fa"

# 1) Filter BED
if [[ ! -s "$bed_in" ]]; then
  echo "ERROR: missing/empty BED: $bed_in" >&2
else
  awk -F'\t' 'BEGIN{OFS="\t"}
    # drop Y, Y1.. and W, W1..
    $1 ~ /^(Y[0-9]*|W[0-9]*)$/ { next }

    # drop if chrom name is longer than 3 chars total
    length($1) > 3 { next }

    { print }
  ' "$bed_in" > "$bed_out"

  echo "Wrote $bed_out"
fi

# 2) Filter peptide FASTA to match filtered BED
if [[ ! -s "$pep_in" ]]; then
  echo "WARN: missing/empty peptide FASTA: $pep_in" >&2
elif [[ ! -s "$bed_out" ]]; then
  echo "WARN: filtered BED missing/empty: $bed_out" >&2
else
  keep_ids="$(mktemp)"
  awk -F'\t' '$4!=""{print $4}' "$bed_out" | sort -u > "$keep_ids"

  seqkit grep -n -f "$keep_ids" "$pep_in" > "$pep_out"

  rm -f "$keep_ids"
  echo "Wrote $pep_out"
fi
```
# Set up genespace
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/Panthera_onca
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/Panthera_onca
mkdir -p bed peptide
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace//sexshared/usable_bed/Gallus_gallus.bed bed/
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace//sexshared/usable_bed/Panthera_onca.bed bed/

cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace//sexshared/usable_peptide/Gallus_gallus.fa peptide/
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace//sexshared/usable_peptide/Panthera_onca.fa peptide/