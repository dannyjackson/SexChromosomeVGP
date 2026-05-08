# Prepare Pseudacris triseriata genome
# Liftover annotations from Hyla sarda
```
/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/analyses/lacks_gff/Pseudacris_triseriata
```
source myconda
mamba activate genespace
FROG=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/analyses/lacks_gff/Pseudacris_triseriata
cd $FROG
```
## Remove scaffolds from FNA files
```
out="Pseudacris_triseriata.contigs_only.fna"
awk 'BEGIN{IGNORECASE=1}
        /^>/{
        drop = ($0 ~ /(scaffold|unplaced|unlocalized|unloc|mitochon)/)
        }
        !drop{print}
    ' Pseudacris_triseriata.fna > "$out"
echo "Pseudacris_triseriata.fna -> $out"
```

## Filter ${SPECIES}/${SPECIES}.gff to just regions kept in file ${SPECIES}/${SPECIES}.contigs_only.fna
```

OUT_SUFFIX=".contigs_only.gff"   # output name: ${SPECIES}/${SPECIES}${OUT_SUFFIX}
fasta="Pseudacris_triseriata.contigs_only.fna"
species="Pseudacris_triseriata"
gff="Pseudacris_triseriata.lifton.REF_Hyla_sarda.SC_0_5.gff"
out="${species}${OUT_SUFFIX}"

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

echo "OK: ${species} -> $out"

```
# Translate cds files from nucleotides to amino acids
gffread -x Pseudacris_triseriata.cds -g Pseudacris_triseriata.contigs_only.fna Pseudacris_triseriata.contigs_only.gff

transeq -sequence Pseudacris_triseriata.cds -outseq Pseudacris_triseriata.translated.cds
```
# Parse annotation files
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure

R

library(GENESPACE)

SPECIES=c("Pseudacris_triseriata")

parsedPaths <- parse_annotations(
  rawGenomeRepo = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/analyses/lacks_gff/", 
  genomeDirs = SPECIES,
  genomeIDs = SPECIES,
  gffString = "contigs_only.gff",
  faString = "translated.cds",
  headerEntryIndex = 1,
  headerStripText = "rna-|_[0-9]+$",
  genespaceWd = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure")

```
Make a mapping file to replace Pseudacris triseriata chromosome names:
```
grep '^>' /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/analyses/lacks_gff/Pseudacris_triseriata/Pseudacris_triseriata.contigs_only.fna |
sed 's/^>//' |
awk '
{
  acc=$1
  chr=""
  if (match($0, /chromosome[[:space:]]+([0-9XYWZ]+)/, m)) chr=m[1]
  if (chr!="") print acc "\t" chr
}' >  /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/Pseudacris_triseriata.chromosome_mapping.tsv

# Replace chromosome 1 with X in mapping file
awk 'BEGIN{OFS="\t"} $2==1{$2="X"} {print}' \
  /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/Pseudacris_triseriata.chromosome_mapping.tsv \
  > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/Pseudacris_triseriata.chromosome_mapping.tsv.tmp \
&& mv /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/Pseudacris_triseriata.chromosome_mapping.tsv.tmp \
      /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/Pseudacris_triseriata.chromosome_mapping.tsv

# Using mapping file, replace chromosome names in bed with numbers/X

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure

R

library(readr)

bed <- read_tsv("bed/Pseudacris_triseriata.bed", col_names = FALSE, show_col_types = FALSE)

map_tbl <- read_tsv(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/Pseudacris_triseriata.chromosome_mapping.tsv",
  col_names = FALSE,
  show_col_types = FALSE
)

# make a named vector: names = accessions, values = chrom labels
map_vec <- setNames(map_tbl$X2, map_tbl$X1)

# remap, leaving anything unmapped unchanged
bed$X1 <- ifelse(bed$X1 %in% names(map_vec), unname(map_vec[bed$X1]), bed$X1)

write_tsv(bed, "bed/Pseudacris_triseriata.remapped.bed", col_names = FALSE)

q()

mv bed/Pseudacris_triseriata.remapped.bed bed/Pseudacris_triseriata.bed
