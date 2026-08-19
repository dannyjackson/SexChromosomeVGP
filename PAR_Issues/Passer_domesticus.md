# Passer domesticus PAR inversion investigation

### Download the raw reads
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues
mkdir -p Passer_domesticus
cd Passer_domesticus

module load aws
aws s3 cp s3://genomeark/species/Passer_domesticus/bPasDom1/genomic_data/pacbio_hifi/ . --recursive --exclude '*' --include '*.ccsreads.fastq.gz' --no-sign-request

wget https://42basepairs.com/browse/s3/genomeark/species/Passer_domesticus/bPasDom1/genomic_data/pacbio_hifi/bPasDom1.cell1.ccsreads.fastq.gz
wget https://42basepairs.com/browse/s3/genomeark/species/Passer_domesticus/bPasDom1/genomic_data/pacbio_hifi/bPasDom1.cell2.ccsreads.fastq.gz
```

Breakpoint: NC_087512:667073

```
#!/bin/bash
#SBATCH --job-name=align_passer_domesticus
#SBATCH --output=slurm_output/%x_%A.out
#SBATCH --error=slurm_output/%x_%A.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

module load minimap2 samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus

FASTA=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Passer_domesticus/ncbi_dataset/data/GCF_036417665.1/GCF_036417665.1_bPasDom1.hap1_genomic.fna
FASTQ1=bPasDom1.cell1.ccsreads.fastq.gz
FASTQ2=bPasDom1.cell2.ccsreads.fastq.gz

minimap2 -ax map-hifi --secondary=no -t 24 "$FASTA" "$FASTQ1" "$FASTQ2" | samtools sort -@ 8 -m 4G -T "bPasDom1_sort" -o "bPasDom1.hifi_vs_hap1.bam"

samtools index -@ 8 bPasDom1.hifi_vs_hap1.bam
```

### Look for split alignments:

```
BAM=bPasDom1.hifi_vs_hap1.bam
samtools view -b "$BAM" "NC_087512.1:1-800000" > NC_087512_1-800kb.bam
samtools index NC_087512_1-800kb.bam
samtools view -c -f 2048 "$BAM" "NC_087512.1:1-800000"
samtools view -f 2048 "$BAM" "NC_087512.1:550000-750000"

samtools view "$BAM" "NC_087512.1:550000-750000" | grep 'SA:Z:' 


samtools view "$BAM" "NC_087512.1:550000-750000" | awk '$6 ~ /[0-9]+[SH]/ {print $1,$2,$3,$4,$5,$6}'

samtools view "$BAM" "NC_087512.1:550000-750000" | awk '$6 ~ /^[0-9]+[SH]/ {print $4,$1,$6}' | sort -n 

# Inspect candidate read
samtools view "$BAM" "NC_087512.1:650000-680000" | grep -F 'm64077_221211_010033/179046732/ccs'

printf '%s\n' 'm64077_221211_010033/179046732/ccs' > candidate_read.txt
samtools view -h -N candidate_read.txt "$BAM" > candidate_read.sam

grep -v '^@' candidate_read.sam

samtools depth -aa -r "NC_087512.1:600000-800000" "$BAM" > depth_600_800kb.txt
awk '{bin=int($2/1000)*1000; sum[bin]+=$3; n[bin]++} END {for (b in sum) print b,sum[b]/n[b]}' depth_600_800kb.txt | sort -n

samtools view "$BAM" "NC_087512.1:550000-950000" | awk '$0 ~ /SA:Z:/ {sa=""; for(i=12;i<=NF;i++) if($i ~ /^SA:Z:/) sa=$i; print $1,$2,$3,$4,$5,$6,sa}'


# check if reads span

samtools view "$BAM" "NC_087512.1:695000-705000" | awk '$4 < 700000 {print $1,$2,$4,$5}'

# Do any regions from 600-1000mb lack any reads spanning?

samtools view -F 2308 "$BAM" "NC_087512.1:600000-1000000" | \
    awk 'BEGIN{OFS="\t"} \
        {start=$4; cigar=$6; ref=0; \
        while(match(cigar,/^[0-9]+[MIDNSHP=X]/))\
        {op=substr(cigar,RSTART,RLENGTH); n=op; gsub(/[^0-9]/,"",n); t=op; gsub(/[0-9]/,"",t); if(t ~ /[MDN=X]/) ref+=n; cigar=substr(cigar,RLENGTH+1)} \
        end=start+ref-1; \
        if(end>=600000 && start<=1000000){if(start<600000) start=600000; if(end>1000000) end=1000000; print start,end}}' | \
        sort -k1,1n -k2,2n | \
        awk 'BEGIN{end=599999} {if($1>end+1) print end+1,$1-1; if($2>end) end=$2} END{if(end<1000000) print end+1,1000000}'

# none output.
# What is the minimum number of reads spanning this region?
samtools view -F 2308 "$BAM" "NC_087512.1:580000-1020000" | awk 'BEGIN{OFS="\t"} {start=$4; cigar=$6; ref=0; while(match(cigar,/^[0-9]+[MIDNSHP=X]/)){op=substr(cigar,RSTART,RLENGTH); n=op; gsub(/[^0-9]/,"",n); t=op; gsub(/[0-9]/,"",t); if(t ~ /[MDN=X]/) ref+=n; cigar=substr(cigar,RLENGTH+1)} print start,start+ref-1}' > read_spans.tsv

for POS in $(seq 600000 1000 1000000); do N=$(awk -v p="$POS" '$1 < p && $2 > p {n++} END{print n+0}' read_spans.tsv); echo "$POS $N"; done > spanning_reads_1kb.txt

sort -k2 -n spanning_reads_1kb.txt | head -n 50

# inspect reads aligning to this break point

samtools view "$BAM" "NC_087512.1:670000-682000" | awk '{print $1,$2,$4,$5,$6}'

awk '$1 < 676000 && $2 > 676000' read_spans.tsv
```
### make a plot of overlapping reads
```
library(Rsamtools)
library(GenomicAlignments)
library(ggplot2)

bam <- "bPasDom1.hifi_vs_hap1.bam"

which <- GRanges(
  "NC_087512.1",
  IRanges(650000, 700000)
)

param <- ScanBamParam(
  which = which,
  what = c("qname", "flag", "mapq", "cigar", "pos")
)

x <- scanBam(bam, param = param)[[1]]

df <- data.frame(
  read = x$qname,
  flag = x$flag,
  mapq = x$mapq,
  cigar = x$cigar,
  start = x$pos
)

cigar_ref_width <- function(cigar) {
  cigarWidthAlongReferenceSpace(cigar)
}

df$end <- df$start + cigar_ref_width(df$cigar) - 1

df$secondary <- bitwAnd(df$flag, 256) != 0
df$supplementary <- bitwAnd(df$flag, 2048) != 0
df$reverse <- bitwAnd(df$flag, 16) != 0

df_primary <- subset(df, !secondary & !supplementary)

df_primary <- df_primary[order(df_primary$start), ]
df_primary$y <- seq_len(nrow(df_primary))

p <- ggplot(df_primary) +
  geom_segment(
    aes(
      x = start,
      xend = end,
      y = y,
      yend = y,
      color = mapq
    ),
    linewidth = 1.5
  ) +
  coord_cartesian(xlim = c(650000, 700000)) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "NC_087512.1 position (bp)",
    y = "HiFi read",
    color = "MAPQ",
    title = "PacBio HiFi alignments across NC_087512.1:650-700 kb"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave("PAR_nonPAR_boundary.png", p)

```
### plot wider region
```
library(Rsamtools)
library(GenomicAlignments)
library(ggplot2)

bam <- "bPasDom1.hifi_vs_hap1.bam"

which <- GRanges(
  "NC_087512.1",
  IRanges(650000, 1000000)
)

param <- ScanBamParam(
  which = which,
  what = c("qname", "flag", "mapq", "cigar", "pos")
)

x <- scanBam(bam, param = param)[[1]]

df <- data.frame(
  read = x$qname,
  flag = x$flag,
  mapq = x$mapq,
  cigar = x$cigar,
  start = x$pos
)

cigar_ref_width <- function(cigar) {
  cigarWidthAlongReferenceSpace(cigar)
}

df$end <- df$start + cigar_ref_width(df$cigar) - 1

df$secondary <- bitwAnd(df$flag, 256) != 0
df$supplementary <- bitwAnd(df$flag, 2048) != 0
df$reverse <- bitwAnd(df$flag, 16) != 0

df_primary <- subset(df, !secondary & !supplementary)

df_primary <- df_primary[order(df_primary$start), ]
df_primary$y <- seq_len(nrow(df_primary))

p <- ggplot(df_primary) +
  geom_segment(
    aes(
      x = start,
      xend = end,
      y = y,
      yend = y,
      color = mapq
    ),
    linewidth = 1.5
  ) +
  coord_cartesian(xlim = c(650000, 1000000)) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "NC_087512.1 position (bp)",
    y = "HiFi read",
    color = "MAPQ",
    title = "PacBio HiFi alignments across NC_087512.1:650-1000 kb"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(
  "PAR_nonPAR_boundary_650_1000kb.png",
  p,
  width = 12,
  height = 8,
  dpi = 300
)

```
### Plot mean depth 650-1000kb
```
samtools depth -aa -r "NC_087512.1:650000-1000000" bPasDom1.hifi_vs_hap1.bam > depth_650_1000kb.txt

R 

library(ggplot2)

depth <- read.table(
  "depth_650_1000kb.txt",
  header = FALSE,
  col.names = c("chr", "pos", "depth")
)

p <- ggplot(depth, aes(x = pos, y = depth)) +
  geom_line() +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "NC_087512.1 position (bp)",
    y = "Read depth",
    title = "PacBio HiFi coverage across NC_087512.1:650-1000 kb"
  ) +
  theme_classic()

ggsave(
  "read_depth_650_1000kb.png",
  p,
  width = 10,
  height = 5,
  dpi = 300
)
```
### Plot mean depth 0-1000kb
```

samtools depth -aa -r "NC_087512.1:0-1000000" bPasDom1.hifi_vs_hap1.bam > depth_0_1000kb.txt

R 

library(ggplot2)

depth <- read.table(
  "depth_0_1000kb.txt",
  header = FALSE,
  col.names = c("chr", "pos", "depth")
)

p <- ggplot(depth, aes(x = pos, y = depth)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "NC_087512.1 position (bp)",
    y = "Read depth",
    title = "PacBio HiFi coverage across NC_087512.1:0-1000 kb"
  ) +
  theme_classic()

ggsave(
  "read_depth_0_1000kb.png",
  p,
  width = 10,
  height = 5,
  dpi = 300
)
```

### Identify the exact point with no depth of coverage
```
awk '$3 == 0 {print}' depth_650_1000kb.txt | head

NC_087512.1     676818  0

awk '$3 < 3 {print}' depth_650_1000kb.txt | head

NC_087512.1     675659  1
NC_087512.1     676216  1
NC_087512.1     676237  1
NC_087512.1     676469  1
NC_087512.1     676818  0
NC_087512.1     677549  1
```

# Align short read genomes to the assembly and look for signals of recombination suppression

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
# Remove W from reference genome
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
# Index reference genome
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
# align
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
# Sort bams
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
# Call variants
## First, generate list of chromosomes
```
module load samtools
REF_GENOME="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Passer_domesticus/shortread/ref_without_NC_087511.1.fa"
samtools faidx "$REF_GENOME"

cut -f1 "${REF_GENOME}.fai" > contigs.txt
```
## Then, call variants on the Z
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

# Compute average heterozygosity
```
mkdir -p heterozygosity
module load vcftools bcftools

bcftools view NC_087512.1.bcf | sed -E '/^#CHROM/ s#(/[^[:space:]]*/)(ERR[0-9]+)\.sorted\.bam#\2#g' | bgzip > NC_087512.1.renamed.vcf.gz

zcat NC_087512.1.renamed.vcf.gz | grep CHROM
VCF=NC_087512.1.renamed.vcf.gz
vcftools --gzvcf "$VCF" --het --out heterozygosity/NC_087512.1
```

# Compute average heterozygosity from these output files using this equation: Observed heterozygosity H_O = (N_Sites - O(HOM)) / N_Sites (from the .het file)
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
{$7 = ($6 < 0.1 ? "Male" : "Female"); print}
' "$outfile" > "${outfile%.txt}_sex.txt"
```
# compute windowed heterozygosity
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


# ------------------------------------------------------------
# Calculate per-individual heterozygosity in 200-kb windows
# ------------------------------------------------------------

WIN=200000

# Determine chromosome length from VCF header
CHR_LENGTH=$(bcftools query \
    -f '%POS\n' \
    ${OUT}.snps.vcf.gz | tail -1)

echo "Length represented in VCF: $CHR_LENGTH"
echo "Window size: $WIN"

HETOUT=heterozygosity/NC_087512.1.hets_200kb.tsv

CHR=NC_087512.1
WIN=200000

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
# Plot in R

```
# --- packages ---
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

OUT <- "heterozygosity/NC_087512.1"


# ------------------------------------------------------------
# Read 50-kb VCFtools heterozygosity
# ------------------------------------------------------------

hets <- read_tsv(
  paste0(OUT, ".hets_50kb.tsv"),
  show_col_types = FALSE
)


# Species assignment
hets <- hets %>%
  mutate(
    species = substr(INDV, 1, 4)
  )


# ------------------------------------------------------------
# Sex assignment
# ------------------------------------------------------------

sex_file <- "heterozygosity/NC_087512.1.with_HO.tsv_sex.txt"


sex <- read_tsv(
  sex_file,
  show_col_types = FALSE
) %>%
  select(INDV, Sex) %>%
  rename(sex = Sex)


hets_sex <- hets %>%
  left_join(sex, by = "INDV") %>%
  mutate(
    sex = str_squish(sex),
    sex = coalesce(sex, "Unknown"),
    sex = factor(
      sex,
      levels = c("Female", "Male", "Unknown")
    )
  )


# ------------------------------------------------------------
# Optional filtering
#
# 10-kb windows may contain very few SNPs, so retaining
# N_SITES lets us remove poorly estimated windows.
# ------------------------------------------------------------

min_sites <- 5

hets_plot <- hets_sex %>%
  filter(N_SITES >= min_sites)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p <- ggplot(
  hets_plot,
  aes(
    x = BIN_MID,
    y = H_O,
    color = sex,
    group = INDV
  )
) +
  geom_line(
    linewidth = 0.2,
    alpha = 0.95
  ) +
  scale_color_manual(
    values = c(
      Female = "#E26D5A",
      Male = "#4F7CAC",
      Unknown = "grey70"
    )
  ) +
  labs(
    title = "Z Chromosome (10 kb bins)",
    x = "Position (bp)",
    y = "Observed heterozygosity",
    color = "Sex"
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.35, "lines"),
    legend.text = element_text(size = 18)
  )


ggsave(
  paste0(OUT, ".hets_50kb.pdf"),
  p,
  width = 12,
  height = 3,
  dpi = 300
)

# plot first 1,000,000 sites
hets_plot <- hets_sex %>%
  filter(BIN_MID <= 1000000)

p <- ggplot(
  hets_plot,
  aes(
    x = BIN_MID,
    y = H_O,
    color = sex,
    group = INDV
  )
) +
  geom_line(linewidth = 0.4, alpha = 0.9) +
  geom_point(size = 1) +
  scale_color_manual(
    values = c(
      Female = "#E26D5A",
      Male = "#4F7CAC",
      Unknown = "grey70"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1000000),
    breaks = seq(0, 1000000, by = 200000)
  ) +
  labs(
    title = "Z Chromosome: first 1 Mb",
    x = "Position (bp)",
    y = "Observed heterozygosity",
    color = "Sex"
  ) +
  theme_bw(base_size = 18)


ggsave(
  paste0(OUT, ".hets_50kb.first1Mb.pdf"),
  p,
  width = 12,
  height = 3,
  dpi = 300
)



# plot first 2,000,000 sites
hets_plot <- hets_sex %>%
  filter(BIN_MID <= 2000000)

p <- ggplot(
  hets_plot,
  aes(
    x = BIN_MID,
    y = H_O,
    color = sex,
    group = INDV
  )
) +
  geom_line(linewidth = 0.4, alpha = 0.9) +
  geom_point(size = 1) +
  scale_color_manual(
    values = c(
      Female = "#E26D5A",
      Male = "#4F7CAC",
      Unknown = "grey70"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1000000),
    breaks = seq(0, 1000000, by = 200000)
  ) +
  labs(
    title = "Z Chromosome: first 2 Mb",
    x = "Position (bp)",
    y = "Observed heterozygosity",
    color = "Sex"
  ) +
  theme_bw(base_size = 18)


ggsave(
  paste0(OUT, ".hets_50kb.first2Mb.pdf"),
  p,
  width = 12,
  height = 3,
  dpi = 300
)



# plot first 3,000,000 sites
hets_plot <- hets_sex %>%
  filter(BIN_MID <= 3000000)

p <- ggplot(
  hets_plot,
  aes(
    x = BIN_MID,
    y = H_O,
    color = sex,
    group = INDV
  )
) +
  geom_line(linewidth = 0.4, alpha = 0.9) +
  geom_point(size = 1) +
  scale_color_manual(
    values = c(
      Female = "#E26D5A",
      Male = "#4F7CAC",
      Unknown = "grey70"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1000000),
    breaks = seq(0, 1000000, by = 200000)
  ) +
  labs(
    title = "Z Chromosome: first 3 Mb",
    x = "Position (bp)",
    y = "Observed heterozygosity",
    color = "Sex"
  ) +
  theme_bw(base_size = 18)


ggsave(
  paste0(OUT, ".hets_50kb.first3Mb.pdf"),
  p,
  width = 12,
  height = 3,
  dpi = 300
)
```
# Make stepwise plot of differences in heterozygosity between males and females
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


# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

OUT <- "heterozygosity/NC_087512.1"

het_file <- paste0(OUT, ".hets_50kb.tsv")

sex_file <- "heterozygosity/NC_087512.1.with_HO.tsv_sex.txt"


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


# ------------------------------------------------------------
# Species assignment
# ------------------------------------------------------------

hets <- hets %>%
  mutate(
    species = substr(INDV, 1, 4)
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
  select(INDV, Sex) %>%
  rename(sex = Sex)


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


# Check sex assignments
print(table(hets_sex$sex, useNA = "ifany"))


# ------------------------------------------------------------
# Optional filtering
#
# Remove windows with very little genotype information.
#
# Since these are 200-kb windows, adjust this threshold
# depending on the SNP density in your data.
# ------------------------------------------------------------

min_sites <- 5

hets_sex <- hets_sex %>%
  filter(
    N_SITES >= min_sites,
    sex %in% c("Female", "Male")
  )


# ------------------------------------------------------------
# Plot individual heterozygosity by sex
# ------------------------------------------------------------

p_individual <- ggplot(
  hets_sex,
  aes(
    x = BIN_MID,
    y = H_O,
    color = sex,
    group = INDV
  )
) +
  geom_line(
    linewidth = 0.3,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c(
      Female = "#E26D5A",
      Male = "#4F7CAC"
    )
  ) +
  labs(
    x = "Position (bp)",
    y = "Heterozygosity",
    color = "Sex"
  ) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )


ggsave(
  paste0(OUT, ".hets_200kb.by_sex.pdf"),
  p_individual,
  width = 12,
  height = 3,
  dpi = 300
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
        alternative = "greater"
      )$p.value,

      error = function(e) NA_real_
    ),

    .groups = "drop"
  )


# ------------------------------------------------------------
# Multiple-testing correction
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
# Stepwise F:M heterozygosity plot
# ============================================================


# ------------------------------------------------------------
# Plot each genomic window as a horizontal segment.
#
# Black = FDR-adjusted p < 0.001
# Grey  = not significant
#
# Using BIN_START / BIN_END makes each statistic span its
# actual genomic window.
# ------------------------------------------------------------

p_ratio <- ggplot() +

  geom_segment(
    data = hets_ratio,
    aes(
      x = BIN_START,
      xend = BIN_END,
      y = ratio,
      yend = ratio,
      color = sig_flag
    ),
    linewidth = 0.9,
    alpha = 0.9,
    lineend = "butt"
  ) +

  # Reference line:
  # F:M = 1 means equal female and male heterozygosity
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +

  scale_color_manual(
    values = c(
      `FALSE` = "grey70",
      `TRUE` = "black"
    ),
    breaks = c(FALSE, TRUE),
    labels = c(
      "Not significant",
      "Significant"
    )
  ) +

  labs(
    x = "Position (bp)",
    y = "H F:M",
    color = "Significance"
  ) +

  theme_bw(
    base_size = 15
  ) +

  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )


p_ratio


ggsave(
  paste0(
    OUT,
    ".hets_200kb.FM_ratio.pdf"
  ),
  p_ratio,
  width = 12,
  height = 3,
  dpi = 300
)


# ============================================================
# Optional: plot only first 3 Mb
# ============================================================

hets_ratio_1Mb <- hets_ratio %>%
  filter(
    BIN_START < 3000000
  )


p_het <- ggplot() +

  geom_segment(
    data = hets_ratio_1Mb,
    aes(
      x = BIN_START,
      xend = pmin(BIN_END, 3000000),
      y = ratio,
      yend = ratio,
      color = sig_flag
    ),
    linewidth = 1.1,
    alpha = 0.9,
    lineend = "butt"
  ) +

  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +

  scale_color_manual(
    values = c(
      `FALSE` = "grey70",
      `TRUE` = "black"
    ),
    breaks = c(FALSE, TRUE),
    labels = c(
      "Not significant",
      "Significant"
    )
  ) +

  scale_x_continuous(
    limits = c(0, 3000000),
    breaks = seq(
      0,
      3000000,
      by = 500000
    ),
    expand = c(0, 0)
  ) +

  labs(
    title = "Female : Male heterozygosity — first 3 Mb",
    x = "Position (bp)",
    y = "H F:M",
    color = "Significance"
  ) +

  annotate(
    "rect",
    xmin = 0,
    xmax = 667073,
    ymin = 0.5,
    ymax = 20,
    fill = NA,
    color = "#8a66ac",
    linewidth = 1
  ) +

  theme_bw(
    base_size = 15
  ) +

  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )



ggsave(
  paste0(
    OUT,
    ".hets_200kb.FM_ratio.first_1Mb.pdf"
  ),
  p_het,
  width = 8,
  height = 4,
  dpi = 300
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
xmax <- 3000000

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
# Combine plots
# ============================================================
library(patchwork)

combined <- p_het / p_identity +
  plot_layout(
    heights = c(4, 0.7),
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
    ".hets_50kb.FM_ratio.WZ_identity.pdf"
  ),
  combined,
  width = 12,
  height = 5,
  dpi = 300
)
```