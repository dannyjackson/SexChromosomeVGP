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
ls $MAMMAL
ls $BIRD
ls $SQUA_ONE
ls $SQUA_TWO
ls $CARTILAGE
ls $RAY
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


```

# Parse annotation files
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/

R

library(GENESPACE)

SPECIES=c("Homo_sapiens", "Gallus_gallus") 
# SPECIES=c("Homo_sapiens", "Gallus_gallus", "Anolis_sagrei", "Podarcis_raffonei", "Hoplias_malabaricus", "Narcine_bancroftii")

parsedPaths <- parse_annotations(
  rawGenomeRepo = "/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks", 
  genomeDirs = SPECIES,
  genomeIDs = SPECIES,
  gffString = "gff",
  faString = "translated.cds",
  presets = "ncbi", 
  genespaceWd = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/formatted_annotations")
```
# Run Genespace
Save this as submit_genespace.sh:
```
#!/bin/bash

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace

source myconda

mamba activate genespace

Rscript genespace.r
```
Then run it with:
```
sbatch \
  -c 16 \
  -t 1:00:00 \
  --mem-per-cpu=20G \
  --mail-user=jacksondan@nih.gov \
  --mail-type=ALL \
  --gres=lscratch:500 \
  --output=slurm_output/genespace.%j \
  submit_genespace.sh 
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