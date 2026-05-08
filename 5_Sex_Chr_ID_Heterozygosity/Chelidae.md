# Can we use heterozygosity to ID sex chr in a turtle?
1. Download WGS data
2. Align to Chelidae genome
3. Compare patterns of heterozygosity across the genome


source myconda
mamba activate ncbi_datasets
1. SRA datasets

```
#!/bin/bash
#SBATCH --job-name=sra_fasta
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-35
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Mobula/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Mobula/logs/%x_%A_%a.err

module load sratoolkit/3.3.0

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae || exit 1

mkdir -p logs

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

mkdir -p "${run}_out"

fasterq-dump "$run" -O "${run}_out"

```
# submit
```
sbatch download_fasta_array.sh
```
# accessions.txt
```
SRR25606201
SRR25606202
SRR25606203
SRR25606204
SRR25606205
SRR32517510
SRR32517511
SRR32517512
SRR32517513
SRR32517514
SRR32517515
SRR32517516
SRR32517517
SRR32517518
SRR32517519
SRR32517520
SRR32517521
SRR32517522
SRR32517523
SRR32517524
SRR32517525
SRR32517526
SRR32517527
SRR32517528
SRR32517529
SRR32517530
SRR32517531
SRR32517532
SRR32517533
SRR32517534
SRR32517535
SRR32517536
SRR32517537
SRR32517538
SRR32517539
SRR32517540
SRR32517541
SRR32517542
SRR32517543
SRR32517544
```
# Index Myuchelys georgesi genome
```
#!/bin/bash
#SBATCH --job-name=index_genome
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --gres=lscratch:500
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae/logs/index.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae/logs/index.err

module load bwa
REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Myuchelys_georgesi/ncbi_dataset/data/GCA_040894355.2/GCA_040894355.2_rMyuGeo1.pri_genomic.fna"

bwa index $REF_GENOME
```
# Fast QC
```
#!/bin/bash
#SBATCH --job-name=raw_QC
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-35
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae/logs/%x_%A_%a.err

module load fastqc

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae || exit 1

mkdir -p trimmed_fastas raw_QC_files

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

while read -r run; do
    fastqc \
    ${run}_out/${run}_1.fastq ${run}_out/${run}_2.fastq \
    --outdir raw_QC_files/
    echo $run" raw QC done" >> trim_and_QC_log.txt
done < accessions.txt
```
# align
```
#!/bin/bash
#SBATCH --job-name=align
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-35
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae/logs/%x_%A_%a.err

module load bwa
module load samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae || exit 1

mkdir -p logs bam_files

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Myuchelys_georgesi/ncbi_dataset/data/GCA_040894355.2/GCA_040894355.2_rMyuGeo1.pri_genomic.fna"
R1="${run}_out/${run}_1.fastq"
R2="${run}_out/${run}_2.fastq"

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
  echo "Missing FASTQ files for $run"
  exit 1
fi

bwa mem -t "${SLURM_CPUS_PER_TASK}" "$REF_GENOME" "$R1" "$R2" | \
  samtools view -@ "${SLURM_CPUS_PER_TASK}" -o "bam_files/${run}.bam" -S
```
#
```
#!/bin/bash
#SBATCH --job-name=sortbams
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-95
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae/logs/%x_%A_%a.err

module load samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae || exit 1

mkdir -p logs bam_files

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Myuchelys_georgesi/ncbi_dataset/data/GCA_040894355.2/GCA_040894355.2_rMyuGeo1.pri_genomic.fna"

mkdir -p sorted_bam_files

samtools sort "bam_files/${run}.bam" -o "sorted_bam_files/${run}.sorted.bam"

samtools index "sorted_bam_files/${run}.sorted.bam"
```
# Call variants
## First, generate list of chromosomes
```
module load samtools
REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Myuchelys_georgesi/ncbi_dataset/data/GCA_040894355.2/GCA_040894355.2_rMyuGeo1.pri_genomic.fna"
samtools faidx "$REF_GENOME"

cut -f1 "${REF_GENOME}.fai" > contigs.txt
```
## Then, call variants in parallel using all chrs
```
#!/bin/bash
#SBATCH --job-name=callvariants
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=48:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-109%12
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae/logs/callvariants_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae/logs/callvariants_%A_%a.err

set -euo pipefail

module load bcftools
module load samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae

bamdir="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae/sorted_bam_files"
REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Myuchelys_georgesi/ncbi_dataset/data/GCA_040894355.2/GCA_040894355.2_rMyuGeo1.pri_genomic.fna"

contig=$(sed -n "${SLURM_ARRAY_TASK_ID}p" contigs.txt)

mkdir -p parts

bcftools mpileup \
  --threads "$SLURM_CPUS_PER_TASK" \
  -Ou \
  -f "$REF_GENOME" \
  -r "$contig" \
  -a FORMAT/AD,DP,INFO/AD,SP \
  "${bamdir}"/*.bam \
| bcftools call \
  --threads "$SLURM_CPUS_PER_TASK" \
  -mv -V indels \
  -Ob \
  -o "parts/${contig}.bcf"

bcftools index --threads "$SLURM_CPUS_PER_TASK" "parts/${contig}.bcf"
```
## Concatenate after array finshes
```
bcftools concat \
  -Ob \
  -o Myuchelys_georgesi_snps_multiallelic.bcf \
  $(awk '{print "parts/"$1".bcf"}' contigs.txt)

bcftools index --threads 12 Myuchelys_georgesi_snps_multiallelic.bcf
```
# Rename samples and zip vcf
```
sed -E '/^#CHROM/ s#(/[^[:space:]]*/)(SRR[0-9]+)\.sorted\.bam#\2#g' Myuchelys_georgesi_snps_multiallelic.vcf > Myuchelys_georgesi_snps_multiallelic.renamed.vcf
bgzip Myuchelys_georgesi_snps_multiallelic.renamed.vcf
```
# Remove scaffolds
```
cat Myuchelys_georgesi_snps_multiallelic.vcf | grep -v 'JAV' > Myuchelys_georgesi_snps_multiallelic.contigs.vcf
bgzip Myuchelys_georgesi_snps_multiallelic.contigs.vcf
```
# Compute heterozygosity
zcat Myuchelys_georgesi_snps_multiallelic.contigs.vcf.gz | grep -v "#" | awk '{print $1}' | sort -u > CONTIGS.txt
mkdir -p heterozygosity
module load vcftools

while read -r CHROM; do
    VCF=Myuchelys_georgesi_snps_multiallelic.contigs.vcf.gz
    vcftools --gzvcf "$VCF" --chr "$CHROM" --het --out heterozygosity/"${CHROM}"
done < CONTIGS.txt

# Compute average heterozygosity from these output files using this equation: Observed heterozygosity H_O = (N_Sites - O(HOM)) / N_Sites (from the .het file)

while read -r CHROM; do
  infile=heterozygosity/"${CHROM}.het"
  outfile=heterozygosity/"${CHROM}.with_HO.tsv"

  echo "Processing ${infile}..."
  awk 'BEGIN {OFS="\t"} 
       NR==1 {print $0, "H_O"; next} 
       {
         H_O = ($4 - $2) / $4;  # (N_SITES - O(HOM)) / N_SITES
         print $0, H_O
       }' "$infile" > "$outfile"
done < CONTIGS.txt

# Create a general matrix of all chromosomes and individuals for plotting.

awk 'FNR==1 && NR!=1 { next }
     {
       split(FILENAME,a,"[./]");
       CHROM=a[1]"."a[2];
       if($1=="INDV") next;
       print CHROM, IND, $0
     }' heterozygosity/*.with_HO.tsv \
| awk 'BEGIN{OFS="\t"; print "CHROM","IND","O(HOM)","E(HOM)","N_SITES","F", "H_O"}1' \
> heterozygosity/all_het_matrix.tsv

awk '{$1=$1; gsub(/ +/, "\t"); print}' heterozygosity/all_het_matrix.tsv > heterozygosity/all_het_matrix.tsv.tmp
mv heterozygosity/all_het_matrix.tsv.tmp heterozygosity/all_het_matrix.tsv

# Generate plots of heterozygosity by sex by chromosome (run interactively).
## Make chrom conversion file
```
CM105055,1
CM105056,2
CM105057,3
CM105058,4
CM105059,5
CM105060,6
CM105061,7
CM105062,8
CM105063,9
CM105064,10
CM105065,11
CM105066,12
CM105067,13
CM105068,14
CM105069,15
CM105070,16
CM105071,17
CM105072,18
CM105073,19
CM105074,20
CM105075,21
CM105076,22
CM105077,23
CM105078,24
CM105079,25
```
# Plot heterozygosity using a gradient division of sex
```
library(tidyverse)

het_file   <- "heterozygosity/all_het_matrix.tsv"
chrom_file <- "chrom_conversion.txt"

het <- read_tsv(het_file, na = c("", "NA", ".", "NaN"), show_col_types = FALSE) %>%
  select(CHROM, IND, H_O) %>%
  mutate(
    CHROM = str_remove(CHROM, "^heterozygosity[._]"),
    CHROM = str_replace(CHROM, "^(NC_[0-9]+)$", "\\1.1")
  )

chrom_conversion <- read_csv(
  chrom_file,
  col_names = c("chrom_accession", "chrom_label"),
  show_col_types = FALSE
)

het <- het %>%
  left_join(chrom_conversion, by = c("CHROM" = "chrom_accession")) %>%
  mutate(chrom_plot = chrom_label)

lev_base  <- str_sort(unique(het$chrom_plot), numeric = TRUE)
specials  <- c("Z", "W", "X", "Y", "MT", "M", "Un", "Unplaced")
lev_order <- c(setdiff(lev_base, specials), specials[specials %in% lev_base])

het <- het %>%
  mutate(chrom_plot = factor(chrom_plot, levels = lev_order))


het <- het %>%
  filter(!is.na(chrom_plot), chrom_plot != "13")

# get chromosome 1 heterozygosity per individual
chr1_ref <- het %>%
  filter(chrom_plot == "1") %>%
  group_by(IND) %>%
  summarize(H_chr1 = mean(H_O, na.rm = TRUE), .groups = "drop")

# standardize all chromosomes by chr1 for each individual
het <- het %>%
  left_join(chr1_ref, by = "IND") %>%
  mutate(std_H = H_O / H_chr1)

# use standardized X heterozygosity to define groups
x_std <- het %>%
  filter(chrom_plot == "X") %>%
  group_by(IND) %>%
  summarize(std_H_X = mean(std_H, na.rm = TRUE), .groups = "drop") %>%
  arrange(std_H_X)

n <- nrow(x_std)

x_std <- x_std %>%
  mutate(group = if_else(row_number() <= floor(n / 2), "M", "F"))

het <- het %>%
  left_join(x_std, by = "IND") %>%
  mutate(group = factor(group, levels = c("M", "F")))

pd <- position_dodge(width = 0.8)

# one heterozygosity value per individual from chromosome X
x_colors <- het %>%
  filter(chrom_plot == "X") %>%
  group_by(IND) %>%
  summarize(H_X = mean(H_O, na.rm = TRUE), .groups = "drop")

# join back so each individual's X value is used for all chromosomes
het <- het %>%
  left_join(x_colors, by = "IND")

MIDPOINT <- het %>%
  filter(chrom_plot == "X") %>%
  pull(H_X) %>%
  median(na.rm = TRUE)

  
gg <- ggplot(het, aes(x = chrom_plot, y = std_H, group = IND, color = H_X)) +
  geom_line(alpha = 0.4, linewidth = 0.25) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0),
    alpha = 0.8,
    size = 0.9
  ) +
  scale_color_gradient2(low = "blue", high = "black", midpoint = MIDPOINT, na.value = "grey70") +
  labs(
    x = "Chromosome",
    y = NULL,
    color = "H on X"
  ) +
  theme_bw(base_size = 1) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 10),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )

ggsave("heterozygosity_by_chromosome.Xgradient.pdf", gg, width = 12, height = 3, dpi = 300)

```
# Make a PCA of chromosomes based on heterozygosity score
```
library(tidyverse)

# make chromosome x individual matrix
het_wide <- het %>%
  select(chrom_plot, IND, H_O) %>%
  filter(!is.na(chrom_plot), chrom_plot != "13") %>%
  pivot_wider(names_from = IND, values_from = H_O)

# save chromosome names
chrom_names <- het_wide$chrom_plot

# numeric matrix for PCA
mat <- het_wide %>%
  select(-chrom_plot) %>%
  as.data.frame()

# remove chromosomes with any missing values
keep <- complete.cases(mat)
mat <- mat[keep, , drop = FALSE]
chrom_names <- chrom_names[keep]

# remove chromosomes with no variation across individuals
var_keep <- apply(mat, 1, function(x) sd(x, na.rm = TRUE) > 0)
mat <- mat[var_keep, , drop = FALSE]
chrom_names <- chrom_names[var_keep]

# PCA: chromosomes are observations, individuals are variables
pca <- prcomp(mat, center = TRUE, scale. = TRUE)

# variance explained
var_explained <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

# scores
pca_df <- data.frame(
  CHROM = chrom_names,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

write_tsv(pca_df, "chromosome_pca_scores.tsv")

# plot
gg <- ggplot(pca_df, aes(x = PC1, y = PC2, label = CHROM)) +
  geom_point(color = "black", size = 2) +
  geom_text(vjust = -0.6, size = 3) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_bw(base_size = 12)

ggsave("chromosome_PCA.pdf", gg, width = 6, height = 5, dpi = 300)
```
# Use a binary sex assignment based on top and bottom 50% heterozygosity score on X chromosome
```
library(tidyverse)

het_file   <- "heterozygosity/all_het_matrix.tsv"
chrom_file <- "chrom_conversion.txt"

het <- read_tsv(het_file, na = c("", "NA", ".", "NaN"), show_col_types = FALSE) %>%
  select(CHROM, IND, H_O) %>%
  mutate(
    CHROM = str_remove(CHROM, "^heterozygosity[._]"),
    CHROM = str_replace(CHROM, "^(NC_[0-9]+)$", "\\1.1")
  )

chrom_conversion <- read_csv(
  chrom_file,
  col_names = c("chrom_accession", "chrom_label"),
  show_col_types = FALSE
)

het <- het %>%
  left_join(chrom_conversion, by = c("CHROM" = "chrom_accession")) %>%
  mutate(chrom_plot = chrom_label)

lev_base  <- str_sort(unique(het$chrom_plot), numeric = TRUE)
specials  <- c("Z", "W", "X", "Y", "MT", "M", "Un", "Unplaced")
lev_order <- c(setdiff(lev_base, specials), specials[specials %in% lev_base])

het <- het %>%
  mutate(chrom_plot = factor(chrom_plot, levels = lev_order))

keep_inds <- c(
  "ERR12115347","ERR12115348","ERR12115349","ERR12115350","ERR12115351",
  "ERR12115352","ERR12115353","ERR12115354","ERR12115355","ERR12115356",
  "ERR12115357","ERR12115358","ERR12115359","ERR12115360","ERR12115361",
  "ERR12115362","ERR12115363","ERR12115364","ERR12115365","ERR12115366",
  "ERR12115367","ERR12115368","ERR12115369","ERR12115370","ERR12115371",
  "ERR12115372","ERR12115373","ERR12115374","ERR12115375","ERR12115376",
  "ERR12115377","ERR12115378","ERR12115379","ERR12115380","ERR12115381",
  "ERR12115382","ERR12115383","ERR12115384","ERR12115385","ERR12115386",
  "ERR12115387","ERR12115388","ERR12115389","ERR12115390","ERR12115391",
  "ERR12115392","ERR12115393","ERR12115394","ERR12115395","ERR12115396",
  "ERR12115397","ERR12115398","ERR12115399","ERR12115400","ERR12115401",
  "ERR12115402","ERR12115403","ERR12115404","ERR12115405","ERR12115406",
  "ERR12115407","ERR12115408","ERR12115409","ERR12115410","ERR12115411",
  "ERR12115412","ERR12115413","ERR12115414","ERR12115415","ERR12115416",
  "ERR12115417","ERR12115418","ERR12115419","ERR12115420","ERR12115421",
  "ERR12115422","ERR12115423","ERR12115424","ERR12115425","ERR12115426",
  "ERR12115427","ERR12115428"
)

het <- het %>%
  filter(IND %in% keep_inds)

het <- het %>%
  filter(chrom_plot != "13")

# X-chromosome heterozygosity per individual
x_het <- het %>%
  filter(chrom_plot == "X") %>%
  group_by(IND) %>%
  summarize(H_X = mean(H_O, na.rm = TRUE), .groups = "drop") %>%
  arrange(H_X)

# split into lower half = M, upper half = F
n <- nrow(x_het)

x_het <- x_het %>%
  mutate(group = if_else(row_number() <= floor(n / 2), "M", "F"))

# if you want exact median split with ties handled by <= median:
# med_x <- median(x_het$H_X, na.rm = TRUE)
# x_het <- x_het %>% mutate(group = if_else(H_X <= med_x, "M", "F"))

# join groups back onto full data
het <- het %>%
  left_join(x_het %>% select(IND, H_X, group), by = "IND") %>%
  mutate(group = factor(group, levels = c("M", "F")))

pd <- position_dodge(width = 0.8)

gg <- ggplot(het, aes(x = chrom_plot, y = H_O, fill = group)) +
  geom_violin(
    position = pd,
    width = 0.75,
    trim = FALSE,
    scale = "width",
    color = NA
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    position = pd,
    shape = 95,
    size = 5,
    color = "black"
  ) +
  geom_point(
    aes(color = group),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    alpha = 0.25,
    size = 0.6,
    stroke = 0,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(M = "#4F7CAC", F = "#E26D5A")) +
  scale_color_manual(values = c(M = "#4F7CAC", F = "#E26D5A")) +
  labs(
    x = "Chromosome",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    legend.position = "right"
  )

ggsave("heterozygosity_by_chromosome_MF_violin.pdf", gg, width = 12, height = 3, dpi = 300)
```
# Standardized heterozygosity violin plot
```
library(tidyverse)

het_file   <- "heterozygosity/all_het_matrix.tsv"
chrom_file <- "chrom_conversion.txt"

het <- read_tsv(het_file, na = c("", "NA", ".", "NaN"), show_col_types = FALSE) %>%
  select(CHROM, IND, H_O) %>%
  mutate(
    CHROM = str_remove(CHROM, "^heterozygosity[._]"),
    CHROM = str_replace(CHROM, "^(NC_[0-9]+)$", "\\1.1")
  )

chrom_conversion <- read_csv(
  chrom_file,
  col_names = c("chrom_accession", "chrom_label"),
  show_col_types = FALSE
)

het <- het %>%
  left_join(chrom_conversion, by = c("CHROM" = "chrom_accession")) %>%
  mutate(chrom_plot = chrom_label)

lev_base  <- str_sort(unique(het$chrom_plot), numeric = TRUE)
specials  <- c("Z", "W", "X", "Y", "MT", "M", "Un", "Unplaced")
lev_order <- c(setdiff(lev_base, specials), specials[specials %in% lev_base])

het <- het %>%
  mutate(chrom_plot = factor(chrom_plot, levels = lev_order))

keep_inds <- c(
  "ERR12115347","ERR12115348","ERR12115349","ERR12115350","ERR12115351",
  "ERR12115352","ERR12115353","ERR12115354","ERR12115355","ERR12115356",
  "ERR12115357","ERR12115358","ERR12115359","ERR12115360","ERR12115361",
  "ERR12115362","ERR12115363","ERR12115364","ERR12115365","ERR12115366",
  "ERR12115367","ERR12115368","ERR12115369","ERR12115370","ERR12115371",
  "ERR12115372","ERR12115373","ERR12115374","ERR12115375","ERR12115376",
  "ERR12115377","ERR12115378","ERR12115379","ERR12115380","ERR12115381",
  "ERR12115382","ERR12115383","ERR12115384","ERR12115385","ERR12115386",
  "ERR12115387","ERR12115388","ERR12115389","ERR12115390","ERR12115391",
  "ERR12115392","ERR12115393","ERR12115394","ERR12115395","ERR12115396",
  "ERR12115397","ERR12115398","ERR12115399","ERR12115400","ERR12115401",
  "ERR12115402","ERR12115403","ERR12115404","ERR12115405","ERR12115406",
  "ERR12115407","ERR12115408","ERR12115409","ERR12115410","ERR12115411",
  "ERR12115412","ERR12115413","ERR12115414","ERR12115415","ERR12115416",
  "ERR12115417","ERR12115418","ERR12115419","ERR12115420","ERR12115421",
  "ERR12115422","ERR12115423","ERR12115424","ERR12115425","ERR12115426",
  "ERR12115427","ERR12115428"
)

het <- het %>%
  filter(IND %in% keep_inds)


het <- het %>%
  filter(chrom_plot != "13")

# get chromosome 1 heterozygosity per individual
chr1_ref <- het %>%
  filter(chrom_plot == "1") %>%
  group_by(IND) %>%
  summarize(H_chr1 = mean(H_O, na.rm = TRUE), .groups = "drop")

# standardize all chromosomes by chr1 for each individual
het <- het %>%
  left_join(chr1_ref, by = "IND") %>%
  mutate(std_H = H_O / H_chr1)

# use standardized X heterozygosity to define groups
x_std <- het %>%
  filter(chrom_plot == "X") %>%
  group_by(IND) %>%
  summarize(std_H_X = mean(std_H, na.rm = TRUE), .groups = "drop") %>%
  arrange(std_H_X)

n <- nrow(x_std)

x_std <- x_std %>%
  mutate(group = if_else(row_number() <= floor(n / 2), "M", "F"))

het <- het %>%
  left_join(x_std, by = "IND") %>%
  mutate(group = factor(group, levels = c("M", "F")))

pd <- position_dodge(width = 0.8)

gg <- ggplot(het, aes(x = chrom_plot, y = std_H, fill = group)) +
  geom_violin(
    position = pd,
    width = 0.75,
    trim = FALSE,
    scale = "width",
    color = NA
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    position = pd,
    shape = 95,
    size = 5,
    color = "black"
  ) +
  geom_point(
    aes(fill = group),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.8,
    size = 0.6,
    show.legend = FALSE
    ) +
  scale_fill_manual(values = c(M = "#4F7CAC", F = "#E26D5A")) +
  scale_color_manual(values = c(M = "#4F7CAC", F = "#E26D5A")) +
  labs(
    x = "Chromosome",
    y = "Standardized heterozygosity (H / H_chr1)",
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    legend.position = "right"
  )

ggsave("heterozygosity_by_chromosome_MF_violin_stdH.pdf", gg, width = 12, height = 3, dpi = 300)
```






# Trim?
```
#!/bin/bash
#SBATCH --job-name=trim
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-95
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae/logs/%x_%A_%a.err

module load trimmmomatic

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Chelidae || exit 1

while read -r ID; do
  java -jar $TRIMMOJAR PE -phred33 -threads 12 \
    ${ID}_out/${ID}_1.fastq \
    ${ID}_out/${ID}_2.fastq \
    -baseout trimmed_fastas/"$ID"_trimmed.fq.gz \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90 >> trim_log.txt
done < accessions.txt

```