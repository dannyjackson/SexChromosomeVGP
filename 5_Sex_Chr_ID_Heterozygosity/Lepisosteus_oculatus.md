# What is the sex chromosome in the gar?
1. Download ddRAD data
2. Align to gar genome
3. Assign sex to each sample
5. Infer if any chromosome exhibits sex linked patterns

mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus
```
#!/bin/bash
#SBATCH --job-name=sra_fasta
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-95
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/logs/%x_%A_%a.err

module load sratoolkit/3.3.0

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus || exit 1

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
SRR9702659
SRR9702660
SRR9702661
SRR9702662
SRR9702663
SRR9702664
SRR9702665
SRR9702666
SRR9702667
SRR9702668
SRR9702669
SRR9702670
SRR9702671
SRR9702672
SRR9702673
SRR9702674
SRR9702675
SRR9702676
SRR9702677
SRR9702678
SRR9702679
SRR9702680
SRR9702681
SRR9702682
SRR9702683
SRR9702684
SRR9702685
SRR9702686
SRR9702687
SRR9702688
SRR9702689
SRR9702690
SRR9702691
SRR9702692
SRR9702693
SRR9702694
SRR9702695
SRR9702696
SRR9702697
SRR9702698
SRR9702699
SRR9702700
SRR9702701
SRR9702702
SRR9702703
SRR9702704
SRR9702705
SRR9702706
SRR9702707
SRR9702708
SRR9702709
SRR9702710
SRR9702711
SRR9702712
SRR9702713
SRR9702714
SRR9702715
SRR9702716
SRR9702717
SRR9702718
SRR9702719
SRR9702720
SRR9702721
SRR9702722
SRR9702723
SRR9702724
SRR9702725
SRR9702726
SRR9702727
SRR9702728
SRR9702729
SRR9702730
SRR9702731
SRR9702732
SRR9702733
SRR9702734
SRR9702735
SRR9702736
SRR9702737
SRR9702738
SRR9702739
```
# Index Gar genome
```
#!/bin/bash
#SBATCH --job-name=index_genome
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --gres=lscratch:500
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/logs/index.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/logs/index.err

module load bwa

# Hap with Y chrs
REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Lepisosteus_oculatus/ncbi_dataset/data/GCF_040954835.1/GCF_040954835.1_fLepOcu1.hap2_genomic.fna"

bwa index $REF_GENOME

# Hap with X chrs
source myconda
mamba activate ncbi_datasets


FILES_TO_DOWNLOAD="gff3,rna,cds,protein,genome,seq-report"
ACCESSION=GCA_040954845.1
datasets download genome accession "${ACCESSION}" \
  --include "${FILES_TO_DOWNLOAD}" 

REF_GENOME="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/ncbi_dataset/data/GCA_040954845.1/GCA_040954845.1_fLepOcu1.hap1_genomic.fna"

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
#SBATCH --array=1-95
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/logs/%x_%A_%a.err

module load bwa
module load samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus || exit 1

mkdir -p trimmed_fastas raw_QC_files

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

while read -r run; do
    fastqc \
    ${run}_out/${run}.fastq \
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
#SBATCH --time=4:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-95
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/logs/%x_%A_%a.err

module load bwa
module load samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus || exit 1

mkdir -p logs bam_files

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

# REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Lepisosteus_oculatus/ncbi_dataset/data/GCF_040954835.1/GCF_040954835.1_fLepOcu1.hap2_genomic.fna"

REF_GENOME="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/ncbi_dataset/data/GCA_040954845.1/GCA_040954845.1_fLepOcu1.hap1_genomic.fna"

R1="${run}_out/${run}.fastq"

if [ ! -f "$R1" ] ; then
  echo "Missing FASTQ files for $run"
  exit 1
fi

bwa mem -t "${SLURM_CPUS_PER_TASK}" "$REF_GENOME" "$R1" | \
  samtools view -@ "${SLURM_CPUS_PER_TASK}" -o "bam_files/${run}.bam" -S
```
# Sort and index bams
```
#!/bin/bash
#SBATCH --job-name=sortbams
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-95
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/logs/%x_%A_%a.err

module load samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus || exit 1

mkdir -p logs bam_files

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

# REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Lepisosteus_oculatus/ncbi_dataset/data/GCF_040954835.1/GCF_040954835.1_fLepOcu1.hap2_genomic.fna"

REF_GENOME="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/ncbi_dataset/data/GCA_040954845.1/GCA_040954845.1_fLepOcu1.hap1_genomic.fna"

mkdir -p sorted_bam_files

samtools sort "bam_files/${run}.bam" -o "sorted_bam_files/${run}.sorted.bam"
samtools index "sorted_bam_files/${run}.sorted.bam"

```
# Call variants
```
#!/bin/bash
#SBATCH --job-name=callvariants
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --gres=lscratch:500
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/logs/callvariants.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/logs/callvariants.err

module load bcftools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus || exit 1

bamdir="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/sorted_bam_files/"

# REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Lepisosteus_oculatus/ncbi_dataset/data/GCF_040954835.1/GCF_040954835.1_fLepOcu1.hap2_genomic.fna"

REF_GENOME="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/Sex_Chr_ID/Lepisosteus_oculatus/ncbi_dataset/data/GCA_040954845.1/GCA_040954845.1_fLepOcu1.hap1_genomic.fna"

bcftools mpileup -Ou -f "$REF_GENOME" -a FORMAT/AD,DP,INFO/AD,SP "$bamdir"*.bam | bcftools call -mv -V indels > Lepisosteus_oculatus_snps_multiallelic.vcf
```
# Rename samples and zip vcf
```
sed -E '/^#CHROM/ s#(/[^[:space:]]*/)(SRR[0-9]+)\.sorted\.bam#\2#g' Lepisosteus_oculatus_snps_multiallelic.vcf > Lepisosteus_oculatus_snps_multiallelic.renamed.vcf
bgzip Lepisosteus_oculatus_snps_multiallelic.renamed.vcf
```
# Remove scaffolds
```
zcat Lepisosteus_oculatus_snps_multiallelic.renamed.vcf.gz | grep -v 'NW_' > Lepisosteus_oculatus_snps_multiallelic.renamed.contigs.vcf
bgzip Lepisosteus_oculatus_snps_multiallelic.renamed.contigs.vcf
```
# Compute heterozygosity
zcat Lepisosteus_oculatus_snps_multiallelic.renamed.contigs.vcf.gz | grep -v "#" | awk '{print $1}' | sort -u > CONTIGS.txt
mkdir -p heterozygosity
module load vcftools

while read -r CHROM; do
    VCF=Lepisosteus_oculatus_snps_multiallelic.renamed.contigs.vcf.gz
    vcftools --gzvcf "$VCF" --chr "$CHROM" --het --out heterozygosity/"${CHROM}"
done < CONTIGS.txt

# Didn't properly remove scaffolds before, do so now
rm heterozygosity/JBDI*
grep -v 'JBDI' CONTIGS.txt > CONTIGS.tmp
mv CONTIGS.tmp CONTIGS.txt

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
## Make chrom conversion file (chrom_conversion.txt)
### For Y haplotype
```
NC_090696.1,1
NC_090697.1,2
NC_090698.1,3
NC_090699.1,4
NC_090700.1,5
NC_090701.1,6
NC_090702.1,7
NC_090703.1,8
NC_090704.1,9
NC_090705.1,10
NC_090706.1,11
NC_090707.1,12
NC_090708.1,13
NC_090709.1,14
NC_090710.1,15
NC_090711.1,16
NC_090712.1,17
NC_090713.1,18
NC_090714.1,19
NC_090715.1,20
NC_090716.1,21
NC_090717.1,22
NC_090718.1,23
NC_090719.1,24
NC_090720.1,25
NC_090721.1,26
NC_090722.1,27
NC_090723.1,28
NC_090724.1,29
```
### For X haplotype
```
CM082521,1
CM082522,2
CM082523,3
CM082524,4
CM082525,5
CM082526,6
CM082527,7
CM082528,8
CM082529,9
CM082530,10
CM082531,11
CM082532,12
CM082533,13
CM082534,14
CM082535,15
CM082536,16
CM082537,17
CM082538,18
CM082539,19
CM082540,20
CM082541,21
CM082542,22
CM082543,23
CM082544,24
CM082545,25
CM082546,26
CM082547,27
CM082548,28
CM082549,29
```
# Pull sex ID information
```
outfile="accession_sex.csv"
printf 'run_accession,sex\n' > "$outfile"

while read -r run sample study species; do
  sex=$(
    curl -s -H 'Accept: application/json' \
      "https://www.ebi.ac.uk/biosamples/samples/${sample}" | \
    jq -r '
      [
        .characteristics.sex[]?.text,
        .characteristics.gender[]?.text,
        .characteristics["sample sex"][]?.text,
        .characteristics["donor sex"][]?.text,
        .characteristics.Sex[]?.text,
        .characteristics.Gender[]?.text
      ]
      | map(select(. != null and . != ""))
      | unique
      | if length==0 then "NA" else join(";") end
    '
  )
  printf '%s,%s\n' "$run" "$sex"
done < <(
  curl -sG 'https://www.ebi.ac.uk/ena/portal/api/search' \
    --data-urlencode 'result=read_run' \
    --data-urlencode 'fields=run_accession,sample_accession,study_accession,scientific_name' \
    --data-urlencode "query=$QUERY" \
    --data-urlencode 'format=tsv' |
  tail -n +2
) >> "$outfile"
```
# Plot heterozygosity using assigned sex file (accession_sex.csv)
```
library(tidyverse)

het_file   <- "heterozygosity/all_het_matrix.tsv"
chrom_file <- "chrom_conversion.txt"
sex_file <- "accession_sex.csv"

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

sex <- read_csv(
  sex_file,
  col_names = c("accession", "sex"),
  show_col_types = FALSE
)

het <- het %>%
  left_join(sex, by = c("IND" = "accession")) 
  
lev_base  <- str_sort(unique(het$chrom_plot), numeric = TRUE)
specials  <- c("Z", "W", "X", "Y", "MT", "M", "Un", "Unplaced")
lev_order <- c(setdiff(lev_base, specials), specials[specials %in% lev_base])

het <- het %>%
  mutate(chrom_plot = factor(chrom_plot, levels = lev_order))

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


pd <- position_dodge(width = 0.8)

  
gg <- ggplot(het, aes(x = chrom_plot, y = H_O, group = IND, color = sex)) +
  geom_line(alpha = 0.4, linewidth = 0.25) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0),
    alpha = 0.8,
    size = 0.9
  ) +
  scale_color_manual(
    values = c(
      male = "blue",
      female = "black"
    ),
    na.value = "grey70"
  ) +
  labs(
    x = "Chromosome",
    y = NULL,
    color = "sex"
  ) +
  theme_bw(base_size = 1) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 10),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )
  
ggsave("heterozygosity_by_chromosome.sex.pdf", gg, width = 12, height = 3, dpi = 300)
  
gg <- ggplot(het, aes(x = chrom_plot, y = std_H, group = IND, color = sex)) +
  geom_line(alpha = 0.4, linewidth = 0.25) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0),
    alpha = 0.8,
    size = 0.9
  ) +
  scale_color_manual(
    values = c(
      male = "blue",
      female = "black"
    ),
    na.value = "grey70"
  ) +
  labs(
    x = "Chromosome",
    y = NULL,
    color = "sex"
  ) +
  theme_bw(base_size = 1) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 10),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )
  
ggsave("std_heterozygosity_by_chromosome.sex.pdf", gg, width = 12, height = 3, dpi = 300)

```
# Make a PCA of chromosomes based on heterozygosity score
```
library(tidyverse)

# make chromosome x individual matrix
het_wide <- het %>%
  select(chrom_plot, IND, H_O) %>%
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
# Use a binary sex assignment
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

sex_file <- "accession_sex.csv"
sex <- read_csv(
  sex_file,
  col_names = c("accession", "sex"),
  show_col_types = FALSE
)

het <- het %>%
  left_join(sex, by = c("IND" = "accession")) 

pd <- position_dodge(width = 0.8)

gg <- ggplot(het, aes(x = chrom_plot, y = H_O, fill = sex)) +
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
    aes(color = sex),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    alpha = 0.25,
    size = 0.6,
    stroke = 0,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(male = "#4F7CAC", female = "#E26D5A")) +
  scale_color_manual(values = c(male = "#4F7CAC", female = "#E26D5A")) +
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

# plot standardized heterozygosity

# get chromosome 1 heterozygosity per individual
chr1_ref <- het %>%
  filter(chrom_plot == "1") %>%
  group_by(IND) %>%
  summarize(H_chr1 = mean(H_O, na.rm = TRUE), .groups = "drop")

# standardize all chromosomes by chr1 for each individual
het <- het %>%
  left_join(chr1_ref, by = "IND") %>%
  mutate(std_H = H_O / H_chr1)

write.csv(het, "std_het_sex.csv")

gg <- ggplot(het, aes(x = chrom_plot, y = std_H, fill = sex)) +
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
    aes(color = sex),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    alpha = 0.25,
    size = 0.6,
    stroke = 0,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(male = "#4F7CAC", female = "#E26D5A")) +
  scale_color_manual(values = c(male = "#4F7CAC", female = "#E26D5A")) +
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

ggsave("std_heterozygosity_by_chromosome_MF_violin.pdf", gg, width = 12, height = 3, dpi = 300)

```

# Plot manhattan plot of chr of interest
```
