# VGP main data freeze paper GeneSpace plot
https://github.com/jtlovell/GENESPACE
## Set up environment
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/yamls

source myconda

mamba env create -f /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/yamls/genespace.yml

mamba env remove -n genespace_clean
mamba list envs
mamba activate genespace

mamba install -c bioconda mcscanx
mamba install bioconda::emboss
mamba install -c bioconda bioconductor-biostrings
mamba install -c conda-forge r-igraph
mamba install -c bioconda gffread
mamba install -c bioconda samtools


mamba env export > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/yamls/genespace.exported.yml


R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jtlovell/GENESPACE")


```
## Set up directory
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure
```
## Identify dataset genomes
```
GENOME_DIR=/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks
MAMMAL=${GENOME_DIR}/Homo_sapiens/
BIRD=${GENOME_DIR}/Gallus_gallus/
SQUA_ONE=${GENOME_DIR}/Anolis_sagrei/
SQUA_TWO=${GENOME_DIR}/Podarcis_raffonei/
CARTILAGE=${GENOME_DIR}/Hoplias_malabaricus/
RAY=${GENOME_DIR}/Narcine_bancroftii/
FROG=${GENOME_DIR}/Pseudacris_triseriata/
ls $MAMMAL
ls $BIRD
ls $SQUA_ONE
ls $SQUA_TWO
ls $CARTILAGE
ls $RAY
ls $FROG
```
# Translate cds files from nucleotides to amino acids
```
cd $MAMMAL
transeq -sequence Homo_sapiens.cds -outseq Homo_sapiens.translated.cds

cd $BIRD
transeq -sequence Gallus_gallus.cds -outseq Gallus_gallus.translated.cds

cd $SQUA_ONE 
transeq -sequence Anolis_sagrei.cds -outseq Anolis_sagrei.translated.cds

cd $SQUA_TWO 
transeq -sequence Podarcis_raffonei.cds -outseq Podarcis_raffonei.translated.cds

cd $CARTILAGE
transeq -sequence Hoplias_malabaricus.cds -outseq Hoplias_malabaricus.translated.cds

cd $RAY
transeq -sequence Narcine_bancroftii.cds -outseq Narcine_bancroftii.translated.cds

cd $FROG
gffread -x Pseudacris_triseriata.lifton.REF_Hyla_sarda.SC_0_5.cds -g Pseudacris_triseriata.fna Pseudacris_triseriata.lifton.REF_Hyla_sarda.SC_0_5.gff

transeq -sequence Pseudacris_triseriata.lifton.REF_Hyla_sarda.SC_0_5.cds -outseq Pseudacris_triseriata.translated.cds

cd /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Hyla_sarda
transeq -sequence Hyla_sarda.cds -outseq Hyla_sarda.translated.cds

```
# Remove Y and W chromosomes
cd /data/Wilson_Lab/data/VGP_genomes_phase1

echo 'Pseudacris_triseriata,X,CM130676.1' >> /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv

chmod +x filter_YW.sh
./filter_YW.sh

# Remove everything but sex chromosomes
chmod +x filter_to_sex_chrs.sh
./filter_to_sex_chrs.sh

# Parse annotation files
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure

R

library(GENESPACE)

SPECIES=c("Homo_sapiens", "Gallus_gallus", "Anolis_sagrei", "Podarcis_raffonei", "Hoplias_malabaricus", "Narcine_bancroftii")

parsedPaths <- parse_annotations(
  rawGenomeRepo = "/data/Wilson_Lab/data/VGP_genomes_phase1/No_YW", 
  genomeDirs = SPECIES,
  genomeIDs = SPECIES,
  gffString = "gff",
  faString = "translated.cds",
  presets = "ncbi", 
  genespaceWd = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure")

SPECIES=c("Pseudacris_triseriata")

parsedPaths <- parse_annotations(
  rawGenomeRepo = "/data/Wilson_Lab/data/VGP_genomes_phase1/No_YW", 
  genomeDirs = SPECIES,
  genomeIDs = SPECIES,
  gffString = "gff",
  faString = "translated.cds",
  headerEntryIndex = 1,
  headerStripText = "rna-|_[0-9]+$",
  genespaceWd = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure")

```
Make a mapping file to replace Pseudacris triseriata chromosome names:
```
grep '^>' /data/Wilson_Lab/data/VGP_genomes_phase1/No_YW/Pseudacris_triseriata/Pseudacris_triseriata.fna |
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

# Replace column 12 with X in Narcine bancroftii

awk 'BEGIN{OFS="\t"} $1=="12"{$1="X"} {print}' bed/Narcine_bancroftii.bed > bed/Narcine_bancroftii.bed.tmp \
&& mv bed/Narcine_bancroftii.bed.tmp bed/Narcine_bancroftii.bed

```
# Run Genespace
Save this as submit_genespace.VGP_MainDF_Figure.sh:
```
#!/bin/bash

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace

source myconda

mamba activate genespace

Rscript /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts/genespace.VGP_MainDF_Figure.r
```
Then run it with:
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts

sbatch \
  -c 16 \
  -t 1:00:00 \
  --mem-per-cpu=20G \
  --mail-user=jacksondan@nih.gov \
  --mail-type=ALL \
  --gres=lscratch:500 \
  --output=slurm_output/genespace_VGP_MainDF_Figure.%j \
  submit_genespace.VGP_MainDF_Figure.sh 
```

       Flagging over-dispered OGs
        ...Anolis_sagrei      :  368 genes in 13 OGs hit > 8 unique places 
        ...Gallus_gallus      :  294 genes in  7 OGs hit > 8 unique places 
        ...Homo_sapiens       :  298 genes in  7 OGs hit > 8 unique places 
        ...Hoplias_malabaricus:  776 genes in 22 OGs hit > 8 unique places 
        ...Narcine_bancroftii : 1010 genes in 27 OGs hit > 8 unique places ***
        ...Podarcis_raffonei  :  638 genes in 12 OGs hit > 8 unique places 
        NOTE! Genomes flagged *** have > 5% of genes in over-dispersed
                orthogroups. These are likely not great annotations, or
                the synteny run contains un-specified WGDs. Regardless,
                these should be examined carefully

# Prune and plot phylogeny
```
module load R/4.5.0

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure

Rscript phylogeny.R
```

## Rerun synteny with just anurans
```
cp ../VGP_DF_Figure/bed/Pseudacris_triseriata.bed ../anuran_plot/bed/
cp ../VGP_DF_Figure/peptide/Pseudacris_triseriata.fa ../anuran_plot/peptide/

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/anuran_plot

R

library(GENESPACE)

SPECIES=c("Hyla_sarda")

parsedPaths <- parse_annotations(
  rawGenomeRepo = "/data/Wilson_Lab/data/VGP_genomes_phase1/No_YW", 
  genomeDirs = SPECIES,
  genomeIDs = SPECIES,
  gffString = "gff",
  faString = "translated.cds",
  presets = "ncbi", 
  genespaceWd = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure")

```
## Same the following as submit_anurans.sh
```
#!/bin/bash

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace

source myconda

mamba activate genespace

Rscript /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts/genespace.anurans.r
```

```
sbatch \
  -c 16 \
  -t 0:30:00 \
  --mem-per-cpu=20G \
  --mail-user=jacksondan@nih.gov \
  --mail-type=ALL \
  --gres=lscratch:500 \
  --output=slurm_output/genespace_VGP_MainDF_Figure.%j \
  submit_anurans.sh
```