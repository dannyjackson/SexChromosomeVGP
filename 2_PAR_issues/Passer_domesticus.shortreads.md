# Align short read genomes to the assembly and look for signals of recombination suppression
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread
```
## Accessions of short reads:
```
ERR2697464
ERR2697465
ERR2697466
ERR2697467
ERR2697468
ERR2697469
ERR2697470
ERR2697471
ERR2697472
ERR2697473
ERR2697477
ERR2697478
ERR2697479
ERR2697480
ERR2697481
ERR2697482
ERR2697483
ERR2697484
ERR2697485
ERR2697486
ERR2697495
```
## Download fastas
```
#!/bin/bash
#SBATCH --job-name=sra_fasta
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-21
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/logs/%x_%A_%a.err

module load sratoolkit/3.3.0

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/fastqs || exit 1

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
## Remove W from reference genome
```
#!/bin/bash
#SBATCH --job-name=index_genome
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --gres=lscratch:500
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/logs/index.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/logs/index.err

module load samtools
REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Passer_domesticus/ncbi_dataset/data/GCF_036417665.1/GCF_036417665.1_bPasDom1.hap1_genomic.fna"

samtools faidx "$REF_GENOME"
samtools faidx "$REF_GENOME" $(cut -f1 "$REF_GENOME.fai" | grep -v '^NC_087511\.1$') > ref_without_NC_087511.1.fa

```
## Index reference genome
```
#!/bin/bash
#SBATCH --job-name=index_genome
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --gres=lscratch:500
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/logs/index.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/logs/index.err

module load bwa
REF_GENOME="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/ref_without_NC_087511.1.fa"

bwa index $REF_GENOME
```
## align
```
#!/bin/bash
#SBATCH --job-name=align
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-35
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/logs/%x_%A_%a.err

module load bwa
module load samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread || exit 1

mkdir -p logs bam_files

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

REF_GENOME="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/ref_without_NC_087511.1.fa"
R1="fastqs/${run}_out/${run}_1.fastq"
R2="fastqs/${run}_out/${run}_2.fastq"

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
  echo "Missing FASTQ files for $run"
  exit 1
fi

bwa mem -t "${SLURM_CPUS_PER_TASK}" "$REF_GENOME" "$R1" "$R2" | \
  samtools view -@ "${SLURM_CPUS_PER_TASK}" -o "bam_files/${run}.bam" -S
```
## Sort bams
```
#!/bin/bash
#SBATCH --job-name=sortbams
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-20
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/logs/%x_%A_%a.err

module load samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread || exit 1

mkdir -p logs bam_files

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

REF_GENOME="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/ref_without_NC_087511.1.fa"

mkdir -p sorted_bam_files

samtools sort "bam_files/${run}.bam" -o "sorted_bam_files/${run}.sorted.bam"

samtools index "sorted_bam_files/${run}.sorted.bam"
```
## Compute depth of Z for each genome
```
#!/bin/bash
#SBATCH --job-name=depth
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=6:00:00
#SBATCH --gres=lscratch:500

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread

module load samtools

module load samtools
    samtools index $bam


for bam in sorted_bam_files/*.bam; do
    sample=$(basename "$bam" .bam)

    avg_depth=$(samtools depth -a -r NC_087512.1 "$bam" \
        | awk '{sum += $3; n++} END {if (n > 0) print sum/n; else print 0}')

    echo -e "${sample}\t${avg_depth}" 

    echo -e "${sample}\t${avg_depth}" >> avg_depth_Z.tsv
done 


for bam in sorted_bam_files/*.bam; do
    sample=$(basename "$bam" .bam)

    avg_depth=$(samtools depth -a -r NC_087474.1 "$bam" \
        | awk '{sum += $3; n++} END {if (n > 0) print sum/n; else print 0}')

    echo -e "${sample}\t${avg_depth}" 

    echo -e "${sample}\t${avg_depth}" >> avg_depth_chr1.tsv
done 

```
## Infer sex from depth
```
library(readr)
library(dplyr)

z_depth <- read_tsv(
  "avg_depth_Z.tsv",
  col_names = c("INDV", "Z_depth"),
  show_col_types = FALSE
)

chr1_depth <- read_tsv(
  "avg_depth_chr1.tsv",
  col_names = c("INDV", "chr1_depth"),
  show_col_types = FALSE
)

depth <- z_depth %>%
  inner_join(chr1_depth, by = "INDV") %>%
  mutate(
    INDV = sub("\\.sorted$", "", INDV),
    Z_chr1_ratio = Z_depth / chr1_depth
  )

# Inspect ratios
depth %>%
  arrange(Z_chr1_ratio) %>%
  select(INDV, Z_depth, chr1_depth, Z_chr1_ratio)

set.seed(1)

km <- kmeans(
  depth$Z_chr1_ratio,
  centers = 2,
  nstart = 100
)

depth$cluster <- km$cluster

# Higher Z:chr1 ratio = ZZ = Male
male_cluster <- which.max(km$centers)

depth <- depth %>%
  mutate(
    sex = if_else(
      cluster == male_cluster,
      "Male",
      "Female"
    )
  )

depth %>%
  arrange(Z_chr1_ratio) %>%
  select(
    INDV,
    Z_depth,
    chr1_depth,
    Z_chr1_ratio,
    sex
  )

table(depth$sex)

km$centers

# save file 
write_tsv(
  depth %>%
    select(
      INDV,
      Z_depth,
      chr1_depth,
      Z_chr1_ratio,
      sex
    ),
  "sex_from_Z_chr1_depth.tsv"
)

```
## Call variants
### First, generate list of chromosomes
```
module load samtools
REF_GENOME="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/ref_without_NC_087511.1.fa"
samtools faidx "$REF_GENOME"

cut -f1 "${REF_GENOME}.fai" > contigs.txt
```
### Then, call variants on the Z
```
#!/bin/bash
#SBATCH --job-name=callvariants
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=48:00:00
#SBATCH --gres=lscratch:500
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/logs/callvariants_%j.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/logs/callvariants_%j.err

set -euo pipefail

module load bcftools
module load samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread

bamdir="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/sorted_bam_files"
REF_GENOME="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/ref_without_NC_087511.1.fa"

contig="NC_087512.1"

bcftools mpileup \
  --threads "$SLURM_CPUS_PER_TASK" \
  -Ou \
  -f "$REF_GENOME" \
  -r "${contig}" \
  -a FORMAT/AD,DP,INFO/AD,SP \
  "${bamdir}"/*.bam \
| bcftools call \
  --threads "$SLURM_CPUS_PER_TASK" \
  -mv -V indels \
  -Ob \
  -o "${contig}.bcf"

bcftools index --threads "$SLURM_CPUS_PER_TASK" "${contig}.bcf"
```

## Compute average heterozygosity
```
mkdir -p heterozygosity
module load vcftools bcftools

bcftools view NC_087512.1.bcf | sed -E '/^#CHROM/ s#(/[^[:space:]]*/)(ERR[0-9]+)\.sorted\.bam#\2#g' | bgzip > NC_087512.1.renamed.vcf.gz

zcat NC_087512.1.renamed.vcf.gz | grep CHROM
VCF=NC_087512.1.renamed.vcf.gz
vcftools --gzvcf "$VCF" --het --out heterozygosity/NC_087512.1
```

## Compute average heterozygosity from these output files using this equation: Observed heterozygosity H_O = (N_Sites - O(HOM)) / N_Sites (from the .het file)
```
CHROM=NC_087512.1
infile=heterozygosity/"${CHROM}.het"
outfile=heterozygosity/"${CHROM}.with_HO.tsv"

echo "Processing ${infile}..."
awk 'BEGIN {OFS="\t"} 
    NR==1 {print $0, "H_O"; next} 
    {
        H_O = ($4 - $2) / $4;  # (N_SITES - O(HOM)) / N_SITES
        print $0, H_O
    }' "$infile" > "$outfile"

# assign sex
awk 'BEGIN{OFS="\t"}
NR==1 {print $0, "Sex"; next}
{$7 = ($6 < 0.1 ? "Female" : "Male"); print}
' "$outfile" > "${outfile%.txt}_sex.txt"
```
## compute windowed heterozygosity
```
### Heterozygosity
VCF=NC_087512.1.renamed.vcf.gz

tabix -p vcf "$VCF"


# Filter to biallelic SNPs on chromosome of interest
bcftools view \
    -v snps \
    -m2 -M2 \
    -O z \
    -o ${OUT}.snps.vcf.gz \
    "$VCF"

tabix -p vcf ${OUT}.snps.vcf.gz


# ------------------------------------------------------------
# Calculate per-individual heterozygosity in 50-kb windows
# ------------------------------------------------------------

WIN=50000

# Determine chromosome length from VCF header
CHR_LENGTH=$(bcftools query \
    -f '%POS\n' \
    ${OUT}.snps.vcf.gz | tail -1)

echo "Length represented in VCF: $CHR_LENGTH"
echo "Window size: $WIN"

HETOUT=heterozygosity/NC_087512.1.hets_50kb.tsv

CHR=NC_087512.1


printf "CHROM\tBIN_START\tBIN_END\tBIN_MID\tINDV\tO_HOM\tE_HOM\tN_SITES\tF\tH_O\n" \
    > "$HETOUT"

for START in $(seq 1 $WIN $CHR_LENGTH); do

    END=$((START + WIN - 1))
    MID=$((START + WIN / 2))

    bcftools view \
        -r "${CHR}:${START}-${END}" \
        -Ov \
        ${OUT}.snps.vcf.gz |
    vcftools \
        --vcf - \
        --het \
        --stdout 2>/dev/null |
    awk \
        -v OFS="\t" \
        -v chr="$CHR" \
        -v start="$START" \
        -v end="$END" \
        -v mid="$MID" '
        NR > 1 && $4 > 0 {
            HO = ($4 - $2) / $4
            print chr, start, end, mid,
                  $1, $2, $3, $4, $5, HO
        }
        ' >> "$HETOUT"

done

echo "Wrote: $HETOUT"

head $HETOUT
```
## compute windowed depth
```
# ============================================================
# Windowed depth on Z chromosome
# Per individual, 50-kb windows
# ============================================================

module load samtools

CHR="NC_087512.1"
WIN=50000

REF_GENOME="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/ref_without_NC_087511.1.fa"
BAMDIR="sorted_bam_files"

mkdir -p depth

DEPTHOUT="depth/${CHR}.depth_50kb.tsv"

# Chromosome length from reference index
CHR_LENGTH=$(awk -v chr="$CHR" '$1 == chr {print $2}' "${REF_GENOME}.fai")

printf "CHROM\tBIN_START\tBIN_END\tBIN_MID\tINDV\tMEAN_DEPTH\n" > "$DEPTHOUT"

for bam in "${BAMDIR}"/*.bam; do

    sample=$(basename "$bam" .bam)
    sample=${sample%.sorted}

    echo "Processing ${sample}..."

    samtools depth \
        -aa \
        -r "$CHR" \
        "$bam" |
    awk \
        -v OFS="\t" \
        -v chr="$CHR" \
        -v win="$WIN" \
        -v chrlen="$CHR_LENGTH" \
        -v ind="$sample" '
        {
            start = int(($2 - 1) / win) * win + 1
            sum[start] += $3
            n[start]++
        }

        END {
            for (start = 1; start <= chrlen; start += win) {

                end = start + win - 1
                if (end > chrlen)
                    end = chrlen

                mid = start + int((end - start) / 2)

                if (n[start] > 0)
                    mean_depth = sum[start] / n[start]
                else
                    mean_depth = 0

                print chr, start, end, mid, ind, mean_depth
            }
        }
        ' >> "$DEPTHOUT"

done

echo "Wrote: $DEPTHOUT"
head "$DEPTHOUT"
```
## Make plot of differences in heterozygosity between males and females
```
# ============================================================
# Female:Male heterozygosity ratio across Z chromosome
# Window-specific t-tests + BH FDR correction
# ============================================================


# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)


# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

OUT <- "heterozygosity/NC_087512.1"

het_file <- paste0(OUT, ".hets_50kb.tsv")

sex_file <- "sex_from_Z_chr1_depth.tsv"

depth_window_file <- "depth/NC_087512.1.depth_50kb.tsv"

# ------------------------------------------------------------
# Read windowed heterozygosity
# ------------------------------------------------------------

hets <- read_tsv(
  het_file,
  show_col_types = FALSE
)


# Inspect
print(head(hets))
print(names(hets))

# ============================================================
# Read windowed depth
# ============================================================


depth_windows <- read_tsv(
  depth_window_file,
  show_col_types = FALSE
)

# Read chromosome 1 mean depth used previously
chr1_depth <- read_tsv(
  "avg_depth_chr1.tsv",
  col_names = c("INDV", "chr1_depth"),
  show_col_types = FALSE
) %>%
  mutate(
    INDV = sub("\\.sorted$", "", INDV)
  )
  


# ------------------------------------------------------------
# Read sex assignments
# ------------------------------------------------------------

sex <- read_tsv(
  sex_file,
  show_col_types = FALSE
)


# Keep only ID and sex columns.
# Assumes sex file contains columns named INDV and Sex.

sex <- sex %>%
  select(INDV, sex) 


# Clean sex labels
sex <- sex %>%
  mutate(
    sex = str_squish(sex),
    sex = factor(
      sex,
      levels = c("Female", "Male")
    )
  )


# ------------------------------------------------------------
# Join sex assignments onto heterozygosity data
# ------------------------------------------------------------

hets_sex <- hets %>%
  left_join(
    sex,
    by = "INDV"
  )



# ------------------------------------------------------------
# Optional filtering
#
# Remove windows with very little genotype information.
#
# Since these are 50-kb windows, adjust this threshold
# depending on the SNP density in your data.
# ------------------------------------------------------------

min_sites <- 10

hets_sex <- hets_sex %>%
  filter(
    N_SITES >= min_sites,
    sex %in% c("Female", "Male")
  )



# ------------------------------------------------------------
# Join sex and autosomal depth
# ------------------------------------------------------------

depth_windows <- depth_windows %>%
  left_join(
    sex,
    by = "INDV"
  ) %>%
  left_join(
    chr1_depth,
    by = "INDV"
  ) %>%
  mutate(
    # normalize Z-window depth for overall sequencing depth
    depth_norm = MEAN_DEPTH / chr1_depth,

    species = substr(INDV, 1, 4)
  ) %>%
  filter(
    sex %in% c("Female", "Male"),
    !is.na(depth_norm),
    is.finite(depth_norm)
  )


# ============================================================
# Female : Male heterozygosity ratio
# ============================================================


# ------------------------------------------------------------
# Calculate mean female and male H_O per window
#
# Also perform a two-sample t-test within each window.
#
# This follows the same approach as the uploaded chr29 script.
# ------------------------------------------------------------

hets_ratio <- hets_sex %>%
  group_by(
    CHROM,
    BIN_START,
    BIN_END,
    BIN_MID
  ) %>%
  summarise(

    # Number of individuals contributing
    n_female = sum(
      sex == "Female" & !is.na(H_O)
    ),

    n_male = sum(
      sex == "Male" & !is.na(H_O)
    ),

    # Mean heterozygosity by sex
    mean_female = mean(
      H_O[sex == "Female"],
      na.rm = TRUE
    ),

    mean_male = mean(
      H_O[sex == "Male"],
      na.rm = TRUE
    ),

    # Female : Male heterozygosity ratio
    ratio = mean_female / mean_male,

    # Window-specific female vs male t-test
    pval = tryCatch(

      t.test(
        H_O[sex == "Female"],
        H_O[sex == "Male"],
        alternative = "two.sided"
      )$p.value,

      error = function(e) NA_real_
    ),

    .groups = "drop"
  )


# ------------------------------------------------------------
# Multiple-testing correction for heterozygosity
#
# Benjamini-Hochberg FDR correction across windows.
# ------------------------------------------------------------

hets_ratio <- hets_ratio %>%
  mutate(

    padj = p.adjust(
      pval,
      method = "BH"
    ),

    sig_flag = case_when(
      !is.na(padj) ~ padj < 0.001,
      TRUE ~ FALSE
    ),

    sig_label = if_else(
      sig_flag,
      "Significant",
      "Not significant"
    )
  )


# Inspect results
print(head(hets_ratio))

print(
  hets_ratio %>%
    select(
      CHROM,
      BIN_START,
      BIN_END,
      mean_female,
      mean_male,
      ratio,
      pval,
      padj,
      sig_flag
    )
)

# ============================================================
# Female vs Male depth statistics
# ============================================================

depth_stats <- depth_windows %>%
  group_by(
    CHROM,
    BIN_START,
    BIN_END,
    BIN_MID
  ) %>%
  summarise(

    n_female = sum(
      sex == "Female" & !is.na(depth_norm)
    ),

    n_male = sum(
      sex == "Male" & !is.na(depth_norm)
    ),

    mean_female = mean(
      depth_norm[sex == "Female"],
      na.rm = TRUE
    ),

    mean_male = mean(
      depth_norm[sex == "Male"],
      na.rm = TRUE
    ),

    ratio = mean_female / mean_male,

    pval = tryCatch(
      t.test(
        depth_norm[sex == "Female"],
        depth_norm[sex == "Male"],
        alternative = "two.sided"
      )$p.value,
      error = function(e) NA_real_
    ),

    .groups = "drop"
  ) %>%
  mutate(

    padj = p.adjust(
      pval,
      method = "BH"
    ),

    sig_flag = case_when(
      !is.na(padj) ~ padj < 0.001,
      TRUE ~ FALSE
    ),

    sig_label = if_else(
      sig_flag,
      "Significant",
      "Not significant"
    )
  )

write_tsv(
  depth_stats,
  "depth/NC_087512.1.depth_50kb.FM_stats.tsv"
)
# ------------------------------------------------------------
# Save statistics
# ------------------------------------------------------------

write_tsv(
  hets_ratio,
  paste0(
    OUT,
    ".hets_200kb.FM_ratio.tsv"
  )
)


# ============================================================
# Female vs Male depth statistics
# ============================================================

depth_stats <- depth_windows %>%
  group_by(
    CHROM,
    BIN_START,
    BIN_END,
    BIN_MID
  ) %>%
  summarise(

    n_female = sum(
      sex == "Female" & !is.na(depth_norm)
    ),

    n_male = sum(
      sex == "Male" & !is.na(depth_norm)
    ),

    mean_female = mean(
      depth_norm[sex == "Female"],
      na.rm = TRUE
    ),

    mean_male = mean(
      depth_norm[sex == "Male"],
      na.rm = TRUE
    ),

    ratio = mean_female / mean_male,

    pval = tryCatch(
      t.test(
        depth_norm[sex == "Female"],
        depth_norm[sex == "Male"],
        alternative = "two.sided"
      )$p.value,
      error = function(e) NA_real_
    ),

    .groups = "drop"
  ) %>%
  mutate(

    padj = p.adjust(
      pval,
      method = "BH"
    ),

    sig_flag = case_when(
      !is.na(padj) ~ padj < 0.001,
      TRUE ~ FALSE
    ),

    sig_label = if_else(
      sig_flag,
      "Significant",
      "Not significant"
    )
  )

write_tsv(
  depth_stats,
  "depth/NC_087512.1.depth_50kb.FM_stats.tsv"
)

# ============================================================
# Read W:Z PAF
# ============================================================
#
# PAF columns:
#
# 1  query name
# 2  query length
# 3  query start
# 4  query end
# 5  strand
# 6  target name
# 7  target length
# 8  target start
# 9  target end
# 10 number matching bases
# 11 alignment block length
# 12 MAPQ
#
# NC_087512.1 is the chromosome used for our x coordinate,
# so we use query_start / query_end.
# ============================================================

paf_file <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/Passer_domesticus_WtoZ.aln.paf"

CHR <- "NC_087512.1"
xmin <- 0
xmax <- 4000000

paf <- read_tsv(
  paf_file,
  col_names = FALSE,
  show_col_types = FALSE,
  progress = FALSE
)


# We only need the first 12 standard PAF columns
paf <- paf %>%
  select(1:12)

colnames(paf) <- c(
  "query",
  "query_length",
  "query_start",
  "query_end",
  "strand",
  "target",
  "target_length",
  "target_start",
  "target_end",
  "n_match",
  "aln_length",
  "mapq"
)


# ============================================================
# Calculate percent sequence identity
# ============================================================

paf_z <- paf %>%
  filter(
    query == CHR,
    aln_length > 0
  ) %>%
  mutate(

    pct_identity =
      100 * n_match / aln_length,

    aln_size =
      query_end - query_start
  )


# Check identity distribution
summary(paf_z$pct_identity)

quantile(
  paf_z$pct_identity,
  probs = c(
    0,
    0.01,
    0.05,
    0.25,
    0.5,
    0.75,
    0.95,
    0.99,
    1
  ),
  na.rm = TRUE
)


# ============================================================
# Identity bins
# ============================================================
#
# IMPORTANT:
#
# Replace these breaks and colors with the EXACT breaks/colors
# from your SV-by-eye plot.
#
# These values are placeholders until that SV plotting code
# is supplied.
# ============================================================

id_breaks <- c(
  -Inf,
  90,
  95,
  97,
  98.5,
  100
)

id_labels <- c(
  "<90%",
  "90-95%",
  "95-97%",
  "97-98.5%",
  "98.5-100%"
)


paf_z <- paf_z %>%
  mutate(
    identity_bin = cut(
      pct_identity,
      breaks = id_breaks,
      labels = id_labels,
      right = FALSE
    )
  )


# ------------------------------------------------------------
# SV-by-eye colors
# ------------------------------------------------------------

identity_colors <- c(
  "<90%"    = "#aacbd7",
  "90-95%"  = "#ece5b1",
  "95-97%"  = "#ece5b1",
  "97-98.5%"  = "#edc699",
  "98.5-100%" = "#ee9b90"
)


# ============================================================
# W:Z sequence identity track
# ============================================================

p_identity <- ggplot(
  paf_z
) +

  geom_rect(
    aes(
      xmin = query_start,
      xmax = query_end,
      ymin = 0,
      ymax = 1,
      fill = identity_bin
    ),
    color = NA
  ) +

  scale_fill_manual(
    values = identity_colors,
    drop = FALSE
  ) +

  scale_x_continuous(
    limits = c(0, xmax),
    expand = c(0, 0),
    labels = scales::label_number(
      scale = 1e-6,
      suffix = " Mb"
    )
  ) +

  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +

  labs(
    x = "Position on Z",
    y = "Z:W\nidentity",
    fill = "% identity"
  ) +

  annotate(
    "rect",
    xmin = 0,
    xmax = 667073,
    ymin = 0.1,
    ymax = 0.9,
    fill = NA,
    color = "#8a66ac",
    linewidth = 1
  ) +

  theme_bw(
    base_size = 15
  ) +

  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),

    legend.position = "right",

    plot.margin = margin(
      t = 0,
      r = 5.5,
      b = 5.5,
      l = 5.5
    )
  )

# ============================================================
# Three-track plot:
# 1. Windowed standardized depth by sex
# 2. Female vs male significance in depth
# 3. Windowed heterozygosity by sex
# 4. Female vs male significance in heterozygosity
# 5. W:Z PAF percent identity
# ============================================================

# ============================================================
# 1. Windowed normalized depth by individual
# ============================================================

depth_plot <- depth_windows %>%
  filter(
    BIN_START < xmax,
    BIN_END > xmin
  )

p_depth <- ggplot(
  depth_plot,
  aes(
    x = BIN_MID,
    y = depth_norm,
    color = sex,
    group = INDV
  )
) +
  geom_line(
    linewidth = 0.35,
    alpha = 0.75
  ) +
  scale_color_manual(
    values = c(
      Female = "#E26D5A",
      Male   = "#4F7CAC"
    )
  ) +
  scale_x_continuous(
    limits = c(xmin, xmax),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0, 2)) +
  labs(
    x = NULL,
    y = "Normalized\ndepth",
    color = "Sex"
  ) +
  theme_bw(
    base_size = 15
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.margin = margin(
      t = 5.5,
      r = 5.5,
      b = 0,
      l = 5.5
    )
  )

# ------------------------------------------------------------
# Restrict heterozygosity data to plotting region
# ------------------------------------------------------------

hets_plot <- hets_sex %>%
  filter(
    BIN_START < xmax,
    BIN_END > xmin
  )

# ============================================================
# 2. Significant female vs male depth difference
# ============================================================

depth_sig_plot <- depth_stats %>%
  filter(
    BIN_START < xmax,
    BIN_END > xmin
  )

p_depth_sig <- ggplot(
  depth_sig_plot
) +
  geom_rect(
    aes(
      xmin = pmax(BIN_START, xmin),
      xmax = pmin(BIN_END, xmax),
      ymin = 0,
      ymax = 1,
      fill = sig_flag
    ),
    color = NA
  ) +
  scale_fill_manual(
    values = c(
      `FALSE` = "grey70",
      `TRUE`  = "black"
    ),
    breaks = c(FALSE, TRUE),
    labels = c(
      "Not significant",
      "Significant"
    )
  ) +
  scale_x_continuous(
    limits = c(xmin, xmax),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "F vs M\ndepth",
    fill = "F vs M depth"
  ) +
  theme_bw(
    base_size = 15
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(
      t = 0,
      r = 5.5,
      b = 0,
      l = 5.5
    )
  )

# ============================================================
# 3. Windowed heterozygosity by sex
# ============================================================

p_het <- ggplot(
  hets_plot,
  aes(
    x = BIN_MID,
    y = H_O,
    color = sex,
    group = INDV
  )
) +

  geom_line(
    linewidth = 0.35,
    alpha = 0.75
  ) +

  scale_color_manual(
    values = c(
      Female = "#E26D5A",
      Male   = "#4F7CAC"
    )
  ) +

  scale_x_continuous(
    limits = c(xmin, xmax),
    expand = c(0, 0),
    labels = scales::label_number(
      scale = 1e-6,
      suffix = " Mb"
    )
  ) +

  labs(
    x = NULL,
    y = "Heterozygosity",
    color = "Sex"
  ) +

  theme_bw(
    base_size = 15
  ) +

  theme(
    panel.grid.minor = element_blank(),

    # x axis is shown only on bottom PAF panel
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),

    legend.position = "right",

    plot.margin = margin(
      t = 5.5,
      r = 5.5,
      b = 0,
      l = 5.5
    )
  )


# ============================================================
# 4. Significant female vs male heterozygosity difference
# ============================================================

sig_plot <- hets_ratio %>%
  filter(
    BIN_START < xmax,
    BIN_END > xmin
  )


p_sig <- ggplot(
  sig_plot
) +

  geom_rect(
    aes(
      xmin = pmax(BIN_START, xmin),
      xmax = pmin(BIN_END, xmax),
      ymin = 0,
      ymax = 1,
      fill = sig_flag
    ),
    color = NA
  ) +

  scale_fill_manual(
    values = c(
      `FALSE` = "grey70",
      `TRUE`  = "black"
    ),
    breaks = c(FALSE, TRUE),
    labels = c(
      "Not significant",
      "Significant"
    )
  ) +

  scale_x_continuous(
    limits = c(xmin, xmax),
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +

  labs(
    x = NULL,
    y = "F vs M",
    fill = "F vs M heterozygosity"
  ) +

  theme_bw(
    base_size = 15
  ) +

  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),

    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),

    panel.grid = element_blank(),

    legend.position = "right",

    plot.margin = margin(
      t = 0,
      r = 5.5,
      b = 0,
      l = 5.5
    )
  )


# ============================================================
# 5. W:Z percent identity
#
# Uses p_identity already constructed above
# ============================================================

p_identity <- p_identity +

  scale_x_continuous(
    limits = c(xmin, xmax),
    expand = c(0, 0),
    labels = scales::label_number(
      scale = 1e-6,
      suffix = " Mb"
    )
  ) +

  labs(
    x = "Position on Z",
    y = "Z:W\nidentity",
    fill = "% identity"
  )


# ============================================================
# Five-track plot
#
# 1. Windowed depth by individual
# 2. Female vs male depth significance
# 3. Windowed heterozygosity by individual
# 4. Female vs male heterozygosity significance
# 5. W:Z percent identity
# ============================================================

combined <- (
  p_depth /
  p_depth_sig /
  p_het /
  p_sig /
  p_identity
) +
  plot_layout(
    heights = c(
      3.0,   # windowed depth
      0.35,  # depth significance
      3.0,   # windowed heterozygosity
      0.35,  # heterozygosity significance
      0.70   # percent identity
    ),
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )


# ============================================================
# Save
# ============================================================

ggsave(
  paste0(
    OUT,
    ".hets_50kb.by_sex.significance.WZ_identity.pdf"
  ),
  combined,
  width = 12,
  height = 4,
  dpi = 300
)
```