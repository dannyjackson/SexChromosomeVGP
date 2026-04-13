
# Create a list of all species with a lifted gff
```
find /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/ \
  -mindepth 1 -maxdepth 1 -type d \
  -regextype posix-extended \
  -regex '.*/[A-Z][a-z]+_[a-z]+' \
  -printf '%f\n' > /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/species_with_lifted_genomes.squamates_VGP.txt

```
rm -r /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Anniella_stebbinsi/
mv /data/Wilson_Lab/data/VGP_genomes_phase1/squamate_nonVGP_genomes/Anniella_stebbinsi /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/

## Eliminate overlapping cds; Translate nucleic acid sequences for all species
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts

while IFS= read -r SPECIES; do
    echo "${SPECIES}"
    gff_matches=(/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/${SPECIES}/lifted.GC*.0_5.gff)
    cds="/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/${SPECIES}/lifted.${SPECIES}.0_5.cds"
    trans_cds="/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/${SPECIES}/lifted.${SPECIES}.0_5.translated.cds"
    gff="${gff_matches[0]}"

    mapfile -t genomes < <(
    find "/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/${SPECIES}/ncbi_dataset/data" \
        -mindepth 2 -maxdepth 2 -name "G*genomic.fna"
    )

    if [[ ${#genomes[@]} -eq 0 ]]; then
    echo "ERROR: no genome FASTA found for ${SPECIES}" >&2
    exit 1
    elif [[ ${#genomes[@]} -gt 1 ]]; then
    echo "ERROR: multiple genome FASTAs found for ${SPECIES}:" >&2
    printf '  %s\n' "${genomes[@]}" >&2
    exit 1
    fi

    genome="${genomes[0]}"

    printf 'genome=%s\n' "$genome"

    awk -F'\t' '
    BEGIN{OFS="\t"}

    # pass 1: identify parents with overlapping CDS
    FNR==NR {
    if($0 ~ /^#/ || $3!="CDS") next
    parent=""
    n=split($9,a,";")
    for(i=1;i<=n;i++){
        if(a[i] ~ /^Parent=/){
        parent=a[i]
        sub(/^Parent=/,"",parent)
        break
        }
    }
    if(parent=="") next

    key=parent SUBSEP $1
    if((key in prev_end) && $4 <= prev_end[key]) bad[parent]=1
    if(!(key in prev_end) || $5 > prev_end[key]) prev_end[key]=$5
    next
    }

    # pass 2: print everything except bad transcripts and their children
    {
    if($0 ~ /^#/) { print; next }

    id=""
    parent=""
    n=split($9,a,";")
    for(i=1;i<=n;i++){
        if(a[i] ~ /^ID=/){ id=a[i]; sub(/^ID=/,"",id) }
        if(a[i] ~ /^Parent=/){ parent=a[i]; sub(/^Parent=/,"",parent) }
    }

    if(id in bad) next
    if(parent in bad) next
    print
    }
    ' "$gff" "$gff" > "${gff%.gff}.no_overlapping_cds.gff"

    gffread -x "$cds" -g "$genome" "${gff%.gff}.no_overlapping_cds.gff"

    transeq -sequence "$cds" -outseq "$trans_cds"

done < /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/species_with_lifted_genomes.squamates_VGP.txt

```
## Prepare all bed and peptide files
```
while IFS= read -r SPECIES; do

  SYMDIR=/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/symlinks/${SPECIES}
  gff_matches=(/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/${SPECIES}/lifted.GC*.0_5.no_overlapping_cds.gff)
  gff="${gff_matches[0]}"
  trans_cds="/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/${SPECIES}/lifted.${SPECIES}.0_5.translated.cds"

  mapfile -t genomes < <(
    find "/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/${SPECIES}/ncbi_dataset/data" \
      -mindepth 2 -maxdepth 2 \
      -name "G*_genomic.fna"
  )

  if [[ ${#genomes[@]} -eq 0 ]]; then
    echo "ERROR: no *_genomic.fna found for ${SPECIES}" >&2
    exit 1
  elif [[ ${#genomes[@]} -gt 1 ]]; then
    echo "ERROR: multiple *_genomic.fna files found for ${SPECIES}:" >&2
    printf '  %s\n' "${genomes[@]}" >&2
    exit 1
  fi

  genome="${genomes[0]}"
  printf 'genome=%s\n' "$genome"

  mkdir -p ${SYMDIR}

  ln -sf "$trans_cds" "${SYMDIR}/${SPECIES}.translated.cds"
  ln -sf "$gff" "${SYMDIR}/${SPECIES}.gff"
  ln -sf $genome "${SYMDIR}/${SPECIES}.fa"

done < /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/species_with_lifted_genomes.squamates_VGP.txt

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts

OUTDIR=/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/Genespace_input_unfiltered
mkdir -p $OUTDIR
Rscript 0b_prepare_genespace.lifted.R "$OUTDIR" /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/symlinks/
```
# rename chromosome numbers in bed file
Make a mapping file to replace chromosome names:
```
while IFS= read -r SPECIES; do

  mapfile -t genomes < <(
    find "/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/${SPECIES}/ncbi_dataset/data" \
      -mindepth 2 -maxdepth 2 \
      -name "G*_genomic.fna"
  )

  if [[ ${#genomes[@]} -eq 0 ]]; then
    echo "ERROR: no *_genomic.fna found for ${SPECIES}" >&2
    exit 1
  elif [[ ${#genomes[@]} -gt 1 ]]; then
    echo "ERROR: multiple *_genomic.fna files found for ${SPECIES}:" >&2
    printf '  %s\n' "${genomes[@]}" >&2
    exit 1
  fi

  genome="${genomes[0]}"
  printf 'genome=%s\n' "$genome"

  grep '^>' "$genome" |
  sed 's/^>//' |
  awk '
  /^>/ {
    line=$0
    sub(/^>/, "", line)

    split(line, a, /[[:space:]]+/)
    acc=a[1]

    if (line ~ /chromosome:?[[:space:]]*[0-9XYWZ]+/) {
      chr=line
      sub(/^.*chromosome:?[[:space:]]*/, "", chr)
      sub(/[^0-9XYWZ].*$/, "", chr)
      print acc "\t" chr
    }
  }' "$genome" > "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/${SPECIES}.chromosome_mapping.tsv"

  echo -e 'File length:'
  wc -l /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/${SPECIES}.chromosome_mapping.tsv

  echo -e '\n'

done < /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/species_with_lifted_genomes.squamates_VGP.txt
```
# Using mapping file, replace chromosome names in bed with numbers/X
```
mkdir -p /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/Genespace_input_unfiltered/bed_reformatted/

while IFS= read -r SPECIES; do
    echo -e 'Starting' ${SPECIES}
    Rscript remap_bed.liftover.squamates_VGP.R ${SPECIES}
    echo -e 'Finished' ${SPECIES} '\n\n'
done < /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/species_with_lifted_genomes.squamates_VGP.txt

while IFS= read -r SPECIES; do
  echo -e '\n' $SPECIES
  awk '{print $1}' /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/Genespace_input_unfiltered/bed_reformatted/${SPECIES}.remapped.bed | sort -u | head
  echo -e '\n'
  awk '{print $1}' /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/Genespace_input_unfiltered/bed/${SPECIES}.bed | sort -u | head
  echo -e '\n'
done < /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/species_with_lifted_genomes.squamates_VGP.txt

while IFS= read -r SPECIES; do
echo -e '\n' $SPECIES
awk '{print $1}' /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/Genespace_input_unfiltered/bed_reformatted/${SPECIES}.remapped.bed | sort -u | wc -l
done < /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/species_with_lifted_genomes.squamates_VGP.txt



while IFS= read -r SPECIES; do
mv /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/Genespace_input_unfiltered/bed_reformatted/${SPECIES}.remapped.bed \
   /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/Genespace_input_unfiltered/bed_reformatted/${SPECIES}.bed
done < /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/species_with_lifted_genomes.squamates_VGP.txt

```

## Remove W and Y chromosomes from bed files, and all scaffolds
```
INDIR=/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/Genespace_input_unfiltered/
OUTDIR=/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/Genespace_input_filtered/


INBED="${INDIR}/bed_reformatted"
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
BED_DIR="/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/Genespace_input_filtered/usable_bed"
PEP_DIR="/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_VGP_lifted_gffs/Genespace_input_filtered/usable_peptide"

REF_BED="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/Option2/bed/Gallus_gallus.bed"
REF_PEP="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/Option2/peptide/Gallus_gallus.fa"

GENESPACEDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared_squamates_VGP

mkdir -p ${GENESPACEDIR}

cd ${GENESPACEDIR}

while IFS= read -r SPECIES; do

    # Make simple per-species structure
    mkdir -p "${SPECIES}/bed" "${SPECIES}/peptide"

    # Symlink species inputs
    ln -sf "${BED_DIR}/${SPECIES}.bed" "${SPECIES}/bed/${SPECIES}.bed"
    ln -sf "${PEP_DIR}/${SPECIES}.fa"  "${SPECIES}/peptide/${SPECIES}.fa"

    # Symlink references (handy to keep alongside each species)
    ln -sf "$REF_BED" "${SPECIES}/bed/Gallus_gallus.bed"
    ln -sf "$REF_PEP" "${SPECIES}/peptide/Gallus_gallus.fa"
    
done < /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/species_with_lifted_genomes.squamates_VGP.txt

```
# Submit genespace
# replace chr 3 with Z in Shinisaurus
awk 'BEGIN{OFS="\t"} {$1=($1=="3"?"Z":$1); print}' \
../Gallus_gallus/sexshared_squamates_VGP/Shinisaurus_crocodilurus/bed/Shinisaurus_crocodilurus.bed \
> tmp && mv tmp ../Gallus_gallus/sexshared_squamates_VGP/Shinisaurus_crocodilurus/bed/Shinisaurus_crocodilurus.bed

# Remake sexchrom_accessions.csv
```
find /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared_squamates_VGP -path "*/bed/*.bed" | while read -r bed; do
    species=$(basename "$bed" .bed)
    awk -v sp="$species" '
        !seen[$1]++ && $1 !~ /^[0-9]+$/ {
            print sp "\t" $1
        }
    ' "$bed"
done > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.lifted.squamates_VGP.csv

sort -u  /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.lifted.squamates_VGP.csv | grep -v Gallus_gallus > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.lifted.squamates_VGP.csv.tmp
mv /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.lifted.squamates_VGP.csv.tmp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.lifted.squamates_VGP.csv

echo -e 'Shinisaurus_crocodilurus\tZ' > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.lifted.squamates_VGP.csv
```
# Submit genespace
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts
N=$(wc -l < /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/species_with_lifted_genomes.squamates_VGP.txt)
sbatch --array=1-${N} 1b_genespace_array.squamates_VGP.sh Gallus_gallus
```

# Species with only bed and peptide subdirs
# All are missing from /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.lifted.csv
Antrozous_pallidus -- no annotated sex chr in genome; prob chr 2 but chr 4 is also synt with Gg1&4
Doryrhina_cyclops -- no annotated sex chr in genome; prob chr 7
Eptesicus_fuscus -- no annotated sex chr in genome; prob chr 1 but chr 2 is also synt with Gg1&4
Rhinolophus_yonghoiseni -- no annotated sex chr in genome; Chr 8 is syntenic with Chicken 1&4, likely X

/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.lifted.No_Sex_Chr_Annot.csv

sbatch --array=1-4 1b_genespace_array.lifted.No_Sex_Chr_Annot.sh Gallus_gallus

# Copy finished directories into full genespace directory
cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared_liftover/* /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/
