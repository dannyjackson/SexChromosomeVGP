# Prepare genespace for liftover files
ls /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/scripts

# Create a list of all species with a lifted gff
```
find /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/ -type d -maxdepth 1


find /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs \
  -mindepth 1 -maxdepth 1 -type d \
  -regextype posix-extended \
  -regex '.*/[A-Z][a-z]+_[a-z]+' \
  -printf '%f\n' > /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/species_with_lifted_genomes.txt
```

# Test species
SPECIES=Rousettus_aegyptiacus


## Remove scaffolds from FNA files
```
cd /data/Wilson_Lab/data/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts

find /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/${SPECIES}/ncbi_dataset/data/ -mindepth 2 -maxdepth 2 -name "*.fna" -print0 |
  while IFS= read -r -d '' f; do
    out="${f%.fna}.contigs_only.fna"
    awk 'BEGIN{IGNORECASE=1}
         /^>/{
           drop = ($0 ~ /(scaffold|unplaced|unlocalized|unloc|mitochon)/)
         }
         !drop{print}
        ' "$f" > "$out"
    echo "$f -> $out"
  done
```

## Filter /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/${SPECIES}/lifted.${SPECIES}.0_5.gff to just regions kept in file ${SPECIES}/${SPECIES}.contigs_only.fna
```
for fasta in /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/${SPECIES}/ncbi_dataset/data/*/*.contigs_only.fna; do
  gff="/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/${SPECIES}/lifted.${SPECIES}.0_5.gff"
  out="/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/${SPECIES}/lifted.${SPECIES}.0_5.contigs_only.gff"

  if [[ ! -s "$gff" ]]; then
    echo "WARN: missing/empty GFF for ${SPECIES}: $gff (skipping)" >&2
    continue
  fi

  awk '
    BEGIN { FS=OFS="\t" }

    # File 1: FASTA -> allowed seqids
    FNR==NR {
      if ($0 ~ /^>/) {
        h=$0
        sub(/^>/, "", h)
        sub(/[ \t].*$/, "", h)   # first token only
        keep[h]=1
      }
      next
    }

    # File 2: GFF -> keep header/comments, and features whose seqid is allowed
    /^#/ { print; next }
    ($1 in keep) { print }
  ' "$fasta" "$gff" > "$out"

  echo "OK: ${SPECIES} -> $out"
done
```
## Translate nucleic acid sequences for all species
```
gff="/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/${SPECIES}/lifted.${SPECIES}.0_5.contigs_only.gff"
genome="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/${SPECIES}/ncbi_dataset/data/*/*.contigs_only.fna"
cds="/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/${SPECIES}/lifted.${SPECIES}.0_5.contigs_only.cds"
trans_cds="/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/${SPECIES}/lifted.${SPECIES}.0_5.contigs_only.translated.cds"

gffread -x "$cds" -g $genome "$gff"

transeq -sequence "$cds" -outseq "$trans_cds"
```
## Prepare all bed and peptide files
```
SYMDIR=/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/symlinks/${SPECIES}
mkdir -p ${SYMDIR}

ln -sf "$trans_cds" "${SYMDIR}/${SPECIES}.translated.cds"
ln -sf "$gff" "${SYMDIR}/${SPECIES}.gff"
ln -sf $genome "${SYMDIR}/${SPECIES}.fa"

OUTDIR=/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/Genespace_input
Rscript 0b_prepare_genespace.lifted.R "$OUTDIR" /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/symlinks/
```
# rename chromosome numbers in bed file
```
Make a mapping file to replace Pseudacris triseriata chromosome names:
```
grep '^>' ${genome} |
sed 's/^>//' |
awk '
{
  acc=$1
  chr=""
  if (match($0, /chromosome[[:space:]]+([0-9XYWZ]+)/, m)) chr=m[1]
  if (chr!="") print acc "\t" chr
}' >  /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/${SPECIES}.chromosome_mapping.tsv

# Using mapping file, replace chromosome names in bed with numbers/X
```
Rscript remap_bed.liftover.R ${SPECIES}

mv /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/Genespace_input/bed/${SPECIES}.remapped.bed /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/Genespace_input/bed/${SPECIES}.bed
```

## Remove W and Y chromosomes from bed files, and all scaffolds
```
INDIR=/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/Genespace_input/
OUTDIR=/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/Genespace_input_filtered/


INBED="${INDIR}/bed"
INPEP="${INDIR}/peptide"

OUTBED="${OUTDIR}/usable_bed"
OUTPEP="${OUTDIR}/usable_peptide"

mkdir -p "$OUTBED" "$OUTPEP"

shopt -s nullglob

# 1) Filter BEDs: drop rows where chrom (col 1) is Y*, W*, and drop rows where chrom (col 1) has > 3 total characters
for bed in "$INBED"/*.bed; do
  outbed="${OUTBED}/$(basename "$bed")"

  awk -F'\t' 'BEGIN{OFS="\t"}
    # drop Y, Y1.. and W, W1..
    $1 ~ /^(Y[0-9]*|W[0-9]*)$/ { next }

    # drop if chrom name is longer than 3 chars total
    length($1) > 3 { next }

    { print }
  ' "$bed" > "$outbed"

  echo "Wrote $outbed"
done

# 2) Filter peptide FASTAs to match filtered BEDs (keep only genes still present)
#    - assumes peptide file name matches bed base name: SPECIES.bed <-> SPECIES.fa
for bed in "$OUTBED"/*.bed; do
  base="$(basename "$bed" .bed)"
  pep_in="${INPEP}/${base}.fa"
  pep_out="${OUTPEP}/${base}.fa"

  if [[ ! -s "$pep_in" ]]; then
    echo "WARN: missing/empty peptide FASTA: $pep_in (skipping)" >&2
    continue
  fi

  # IDs to keep = col 4 from filtered BED
  keep_ids="$(mktemp)"
  awk -F'\t' '$4!=""{print $4}' "$bed" | sort -u > "$keep_ids"

  # Subset FASTA by header name (requires seqkit)
  seqkit grep -n -f "$keep_ids" "$pep_in" > "$pep_out"

  rm -f "$keep_ids"
  echo "Wrote $pep_out"
done
```
# Create requisite genespace directory structure
```
BED_DIR="/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/Genespace_input_filtered/usable_bed"
PEP_DIR="/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/Genespace_input_filtered/usable_peptide"

REF_BED="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/Option2/bed/Gallus_gallus.bed"
REF_PEP="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/Option2/peptide/Gallus_gallus.fa"

GENESPACEDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared_liftover

mkdir -p ${GENESPACEDIR}

cd ${GENESPACEDIR}

# Make simple per-species structure
mkdir -p "${SPECIES}/bed" "${SPECIES}/peptide"

# Symlink species inputs
ln -sf "${BED_DIR}/${SPECIES}.bed" "${SPECIES}/bed/${SPECIES}.bed"
ln -sf "${PEP_DIR}/${SPECIES}.fa"  "${SPECIES}/peptide/${SPECIES}.fa"

# Symlink references (handy to keep alongside each species)
ln -sf "$REF_BED" "${SPECIES}/bed/Gallus_gallus_REF.bed"
ln -sf "$REF_PEP" "${SPECIES}/peptide/Gallus_gallus_REF.fa"
```
# Submit genespace
# Remake sexchrom_accessions.csv
```
export TSV=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/VGP_freeze_hap1_combined_sexchroms_seq_reports.tsv
export GENOMEDIR=/data/Wilson_Lab/data/VGP_genomes

GITREP="~/dannys_githubs/SexChromosomeVGP/"

bash $GITREPO/make_sexchrom_accessions.sh $OUTDIR/sexchrom_accessions.csv
```
N=$(wc -l < /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/species_with_lifted_genomes.txt)
N=1
sbatch --array=1-${N} 1b_genespace_array.lifted.sh Gallus_gallus
```
