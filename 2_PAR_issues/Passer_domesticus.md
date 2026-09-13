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
  IRanges(660000, 690000)
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
  geom_vline(
    xintercept = 667073,
    linetype = "dashed",
    linewidth = 0.7,
    color = "#8a65b9"
  ) +
  coord_cartesian(xlim = c(660000, 690000)) +
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

ggsave("PAR_nonPAR_boundary.pdf", p,
  width = 7,
  height = 7)

ggsave("PAR_nonPAR_boundary.png", p,
  width = 7,
  height = 7)

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

