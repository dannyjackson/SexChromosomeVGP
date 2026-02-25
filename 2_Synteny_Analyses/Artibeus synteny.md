# Artibeus synteny

## Set up directory
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus
```
## Identify dataset genomes
```
GENOME_DIR=/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks
OUTGROUP=${GENOME_DIR}/Homo_sapiens/
ART_LIT=${GENOME_DIR}/Artibeus_lituratus/
ART_INT=${GENOME_DIR}/Artibeus_intermedius/
```
## Lift over annotatiosn for both species

```
source myconda
mamba activate lifton
cd $ART_LIT
minimap2 -d Artibeus_lituratus.fna.mmi Artibeus_lituratus.fna
cd $ART_INT
minimap2 -d Artibeus_intermedius.fna.mmi Artibeus_intermedius.fna

```

### Prepare reference genome for liftover
Remove mitochondrial genes (these are often in a format that breaks lifton, i.e. "gene-12")
```
GFF_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Desmodus_rotundus/ncbi_dataset/data/GCF_022682495.2/genomic.gff
awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $1!="NC_022423.1"  {print}' "$GFF_REF" > /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Desmodus_rotundus/ncbi_dataset/data/GCF_022682495.2/genomic.noMT.gff
```

### Lift over for Artibeus literatus
```
#!/bin/bash

cd /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Artibeus_lituratus

source myconda
mamba activate lifton

# Use Desmodus rotundus; the closest relative with not a CLR based assembly
FNA_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Desmodus_rotundus/ncbi_dataset/data/GCF_022682495.2/GCF_022682495.2_HLdesRot8A.1_genomic.fna

GFF_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Desmodus_rotundus/ncbi_dataset/data/GCF_022682495.2/genomic.noMT.gff

FAA_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Desmodus_rotundus/ncbi_dataset/data/GCF_022682495.2/protein.faa

GFF_QRY=/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Artibeus_lituratus/Artibeus_lituratus.lifton.REF_Desmodus_rotundus.SC_0_5.gff

FNA_QRY=/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Artibeus_lituratus/Artibeus_lituratus.fna


SC=0.5

lifton \
    -g ${GFF_REF} \
    -o ${GFF_QRY} \
    -P ${FAA_REF} \
    -t 16 \
    -sc ${SC} \
    ${FNA_QRY} \
    ${FNA_REF}
```

```
sbatch \
  -c 16 \
  -t 6:00:00 \
  --mem-per-cpu=24G \
  --mail-user=jacksondan@nih.gov \
  --mail-type=ALL \
  --output=slurm_output/lifton_Artibeus_literatus.%j \
  submit_lifton.Artibeus_literatus.sh 
```

### Lift over for Artibeus intermedius
```
#!/bin/bash

cd /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Artibeus_intermedius

source myconda
mamba activate lifton

# Use Desmodus rotundus; the closest relative with not a CLR based assembly
FNA_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Desmodus_rotundus/ncbi_dataset/data/GCF_022682495.2/GCF_022682495.2_HLdesRot8A.1_genomic.fna

GFF_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Desmodus_rotundus/ncbi_dataset/data/GCF_022682495.2/genomic.noMT.gff

FAA_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Desmodus_rotundus/ncbi_dataset/data/GCF_022682495.2/protein.faa

GFF_QRY=/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Artibeus_intermedius/Artibeus_intermedius.lifton.REF_Desmodus_rotundus.SC_0_5.gff

FNA_QRY=/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Artibeus_intermedius/Artibeus_intermedius.fna

SC=0.5

lifton \
    -g ${GFF_REF} \
    -o ${GFF_QRY} \
    -P ${FAA_REF} \
    -t 16 \
    -sc ${SC} \
    ${FNA_QRY} \
    ${FNA_REF}
```

```
sbatch \
  -c 16 \
  -t 6:00:00 \
  --mem-per-cpu=24G \
  --mail-user=jacksondan@nih.gov \
  --mail-type=ALL \
  --output=slurm_output/lifton_Artibeus_intermedius.%j \
  submit_lifton.Artibeus_intermedius.sh 
```

# Translate cds files from nucleotides to amino acids
```
source myconda
mamba activate genespace

cd /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Artibeus_intermedius

gffread -x Artibeus_intermedius.lifton.REF_Desmodus_rotundus.SC_0_5.cds -g Artibeus_intermedius.fna Artibeus_intermedius.lifton.REF_Desmodus_rotundus.SC_0_5.gff

transeq -sequence Artibeus_intermedius.lifton.REF_Desmodus_rotundus.SC_0_5.cds -outseq Artibeus_intermedius.translated.cds

cd /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Artibeus_lituratus

gffread -x Artibeus_lituratus.lifton.REF_Desmodus_rotundus.SC_0_5.cds -g Artibeus_lituratus.fna Artibeus_lituratus.lifton.REF_Desmodus_rotundus.SC_0_5.gff

transeq -sequence Artibeus_lituratus.lifton.REF_Desmodus_rotundus.SC_0_5.cds -outseq Artibeus_lituratus.translated.cds


```

# Parse annotation files
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/figure

R

library(GENESPACE)

SPECIES=c("Homo_sapiens")

parsedPaths <- parse_annotations(
  rawGenomeRepo = "/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/", 
  genomeDirs = SPECIES,
  genomeIDs = SPECIES,
  gffString = "gff",
  faString = "translated.cds",
  presets = "ncbi", 
  genespaceWd = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/figure")

SPECIES=c("Artibeus_lituratus", "Artibeus_intermedius")

parsedPaths <- parse_annotations(
  rawGenomeRepo = "/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/", 
  genomeDirs = SPECIES,
  genomeIDs = SPECIES,
  gffString = "\\.gff$",
  faString = "translated.cds",
  headerEntryIndex = 1,
  headerStripText = "rna-|_[0-9]+$",
  genespaceWd = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/figure")

```
## Make a mapping file to replace chromosome names in each Artibeus:
First, Artibeus lituratus
```
grep '^>' /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Artibeus_lituratus/Artibeus_lituratus.fna |
sed 's/^>//' |
awk '
{
  acc=$1
  chr=""
  if (match($0, /chromosome[[:space:]]+([0-9XYWZ]+)/, m)) chr=m[1]
  if (chr!="") print acc "\t" chr
}' >  /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/Artibeus_lituratus.chromosome_mapping.tsv

# Using mapping file, replace chromosome names in bed with numbers/X

R

library(readr)

bed <- read_tsv("bed/Artibeus_lituratus.bed", col_names = FALSE, show_col_types = FALSE)

map_tbl <- read_tsv(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/Artibeus_lituratus.chromosome_mapping.tsv",
  col_names = FALSE,
  show_col_types = FALSE
)

# make a named vector: names = accessions, values = chrom labels
map_vec <- setNames(map_tbl$X2, map_tbl$X1)

# remap, leaving anything unmapped unchanged
bed$X1 <- ifelse(bed$X1 %in% names(map_vec), unname(map_vec[bed$X1]), bed$X1)

write_tsv(bed, "bed/Artibeus_lituratus.remapped.bed", col_names = FALSE)

q()

mv bed/Artibeus_lituratus.remapped.bed bed/Artibeus_lituratus.bed

```
Second, Artibeus intermedius
```
grep '^>' /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Artibeus_intermedius/Artibeus_intermedius.fna |
sed 's/^>//' |
awk '
{
  acc=$1
  chr=""
  if (match($0, /chromosome[[:space:]]+([0-9XYWZ]+)/, m)) chr=m[1]
  if (chr!="") print acc "\t" chr
}' >  /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/Artibeus_intermedius.chromosome_mapping.tsv

# Using mapping file, replace chromosome names in bed with numbers/X

R

library(readr)

bed <- read_tsv("bed/Artibeus_intermedius.bed", col_names = FALSE, show_col_types = FALSE)

map_tbl <- read_tsv(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/Artibeus_intermedius.chromosome_mapping.tsv",
  col_names = FALSE,
  show_col_types = FALSE
)

# make a named vector: names = accessions, values = chrom labels
map_vec <- setNames(map_tbl$X2, map_tbl$X1)

# remap, leaving anything unmapped unchanged
bed$X1 <- ifelse(bed$X1 %in% names(map_vec), unname(map_vec[bed$X1]), bed$X1)

write_tsv(bed, "bed/Artibeus_intermedius.remapped.bed", col_names = FALSE)

q()

mv bed/Artibeus_intermedius.remapped.bed bed/Artibeus_intermedius.bed
```
# Run Genespace
Save this as submit_genespace.Artibeus.sh:
```
#!/bin/bash

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/

source myconda

mamba activate genespace

Rscript /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/scripts/genespace.Artibeus.r
```
Then run it with:
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/scripts

sbatch \
  -c 16 \
  -t 1:00:00 \
  --mem-per-cpu=20G \
  --mail-user=jacksondan@nih.gov \
  --mail-type=ALL \
  --gres=lscratch:500 \
  --output=slurm_output/genespace_Artibeus.%j \
  submit_genespace.Artibeus.sh 
```
# Redo with version 2
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2/genomes
mamba activate ncbi_datasets
FILES_TO_DOWNLOAD="gff3,rna,cds,protein,genome,seq-report"

datasets download genome accession GCA_038363095.2 \
--include "${FILES_TO_DOWNLOAD}" 
unzip ncbi_dataset.zip

mamba activate lifton
minimap2 -d ncbi_dataset/data/GCA_038363095.2/GCA_038363095.2_mArtLit.hap1_genomic.fna.mmi ncbi_dataset/data/GCA_038363095.2/GCA_038363095.2_mArtLit.hap1_genomic.fna

```
## Liftover
```

#!/bin/bash

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2/genomes

source myconda
mamba activate lifton

# Use Desmodus rotundus; the closest relative with not a CLR based assembly
FNA_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Desmodus_rotundus/ncbi_dataset/data/GCF_022682495.2/GCF_022682495.2_HLdesRot8A.1_genomic.fna

GFF_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Desmodus_rotundus/ncbi_dataset/data/GCF_022682495.2/genomic.noMT.gff

FAA_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Desmodus_rotundus/ncbi_dataset/data/GCF_022682495.2/protein.faa

GFF_QRY=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2/genomes/Artibeus_lituratus.v2.lifton.REF_Desmodus_rotundus.SC_0_5.gff

FNA_QRY=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2/genomes/ncbi_dataset/data/GCA_038363095.2/GCA_038363095.2_mArtLit.hap1_genomic.fna


SC=0.5

lifton \
    -g ${GFF_REF} \
    -o ${GFF_QRY} \
    -P ${FAA_REF} \
    -t 16 \
    -sc ${SC} \
    ${FNA_QRY} \
    ${FNA_REF}
    ```

sbatch \
  -c 16 \
  -t 3:00:00 \
  --mem-per-cpu=24G \
  --mail-user=jacksondan@nih.gov \
  --mail-type=ALL \
  --output=slurm_output/lifton_Artibeus_lituratus.v2.%j \
  submit_lifton.Artibeus_lituratus.v2.sh 
```

## Translate cds files from nucleotides to amino acids
```
source myconda
mamba activate genespace


cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2/symlinks/Artibeus_lituratus

gffread -x Artibeus_lituratus.lifton.REF_Desmodus_rotundus.SC_0_5.cds -g /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2/genomes/ncbi_dataset/data/GCA_038363095.2/GCA_038363095.2_mArtLit.hap1_genomic.fna /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2/genomes/Artibeus_lituratus.v2.lifton.REF_Desmodus_rotundus.SC_0_5.gff

transeq -sequence Artibeus_lituratus.lifton.REF_Desmodus_rotundus.SC_0_5.cds -outseq Artibeus_lituratus.translated.cds

ln /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2/genomes/ncbi_dataset/data/GCA_038363095.2/GCA_038363095.2_mArtLit.hap1_genomic.fna Artibeus_lituratus.fna
ln /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2/genomes/Artibeus_lituratus.v2.lifton.REF_Desmodus_rotundus.SC_0_5.gff Artibeus_lituratus.gff
cd ..
```

## Parse annotation files
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2/symlinks/

R

library(GENESPACE)

SPECIES=c("Artibeus_lituratus")

parsedPaths <- parse_annotations(
  rawGenomeRepo = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2/symlinks/", 
  genomeDirs = SPECIES,
  genomeIDs = SPECIES,
  gffString = "\\.gff$",
  faString = "translated.cds",
  headerEntryIndex = 1,
  headerStripText = "rna-|_[0-9]+$",
  genespaceWd = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2/figure")

```
## Make a mapping file to replace chromosome names in Artibeus lituratus v.2:

```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2/figure

# Using mapping file, replace chromosome names in bed with numbers/X

grep '^>' /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2/symlinks/Artibeus_lituratus/Artibeus_lituratus.fna |
sed 's/^>//' |
awk '
{
  acc=$1
  chr=""
  if (match($0, /chromosome[[:space:]]+([0-9XYWZ]+)/, m)) chr=m[1]
  if (chr!="") print acc "\t" chr
}' >  /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/Artibeus_lituratus.chromosome_mapping.V2.tsv

# Using mapping file, replace chromosome names in bed with numbers/X

R

library(readr)

bed <- read_tsv("bed/Artibeus_lituratus.bed", col_names = FALSE, show_col_types = FALSE)

map_tbl <- read_tsv(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/Artibeus_lituratus.chromosome_mapping.V2.tsv",
  col_names = FALSE,
  show_col_types = FALSE
)

# make a named vector: names = accessions, values = chrom labels
map_vec <- setNames(map_tbl$X2, map_tbl$X1)

# remap, leaving anything unmapped unchanged
bed$X1 <- ifelse(bed$X1 %in% names(map_vec), unname(map_vec[bed$X1]), bed$X1)

write_tsv(bed, "bed/Artibeus_lituratus.remapped.bed", col_names = FALSE)

q()

mv bed/Artibeus_lituratus.remapped.bed bed/Artibeus_lituratus.bed


```
## Copy Artibeus intermedius into wd
```
cp ../../figure/bed/Artibeus_intermedius.bed bed/
cp ../../figure/peptide/Artibeus_intermedius.fa peptide/


cp ../../figure/bed/Homo_sapiens.bed bed/
cp ../../figure/peptide/Homo_sapiens.fa peptide/
```
# Run Genespace
Save this as submit_genespace.Artibeus_lituratus.v2.sh:
```
#!/bin/bash

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/version2

source myconda

mamba activate genespace

Rscript /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/scripts/genespace.Artibeus.v2.r
```
#Then run it with:
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/scripts

sbatch \
  -c 16 \
  -t 1:00:00 \
  --mem-per-cpu=20G \
  --mail-user=jacksondan@nih.gov \
  --mail-type=ALL \
  --gres=lscratch:500 \
  --output=slurm_output/genespace_Artibeus.v2.%j \
  submit_genespace.Artibeus_lituratus.v2.sh