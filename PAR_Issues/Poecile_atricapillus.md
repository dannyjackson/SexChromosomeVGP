# Poecile atricapillus PAR inversion investigation

### Set up the working environment
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues
mkdir -p Poecile_atricapillus
cd Poecile_atricapillus
```

minimap2 -ax map-hifi reference.fa reads.fastq.gz |
  samtools sort -o aligned.bam

samtools index aligned.bam

### Download the fastq file
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Poecile_atricapillus

module load aws
aws configure set default.s3.max_concurrent_requests 50`
aws s3 cp \
  s3://genomeark/species/Poecile_atricapillus/bPoeAtr1/genomic_data/pacbio_hifi/ \
  . \
  --recursive \
  --exclude "*" \
  --include "*.hifi_reads.fastq.gz" \
  --no-sign-request
```
### Align and index
```
#!/bin/bash
#SBATCH --job-name=align_poecile_atricapillus
#SBATCH --output=slurm_output/%x_%A.out
#SBATCH --error=slurm_output/%x_%A.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

module load minimap2 samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Poecile_atricapillus

FASTA=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Poecile_atricapillus/ncbi_dataset/data/GCF_030490865.1/GCF_030490865.1_bPoeAtr1.hap1_genomic.fna
FASTQ1=m54306Ue_220518_182133.demultiplex.bc1001--bc1001.hifi_reads.fastq.gz
FASTQ2=m64055e_220515_161029.demultiplex.bc1001--bc1001.hifi_reads.fastq.gz
FASTQ3=m64055e_220719_201855.demultiplex.bc1001--bc1001.hifi_reads.fastq.gz
FASTQ4=m64330e_220501_035120.demultiplex.bc1001--bc1001.hifi_reads.fastq.gz
FASTQ5=m64334e_220509_031744.demultiplex.bc1001--bc1001.hifi_reads.fastq.gz

minimap2 -ax map-hifi --secondary=no -t 24 "$FASTA" "$FASTQ1" "$FASTQ2" "$FASTQ3" "$FASTQ4" "$FASTQ5" | samtools sort -@ 8 -m 4G -T "Poecile_atricapillus_sort" -o "Poecile_atricapillus.hifi_vs_hap1.bam"

samtools index -@ 8 Poecile_atricapillus.hifi_vs_hap1.bam
```

PAR: 	NC_081289.1	143407603-143747811 
ROI: 	NC_081289.1	143000000-145000000
ROI: 	NC_081289.1	143500000-143407603

### make a plot of overlapping reads telomeric to the PAR
```
library(Rsamtools)
library(GenomicAlignments)
library(ggplot2)

bam <- "Poecile_atricapillus.hifi_vs_hap1.bam"

which <- GRanges(
  "NC_081289.1",
  IRanges(143380000, 143410000)
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
    xintercept = 143407603,
    linetype = "dashed",
    linewidth = 0.7,
    color = "#8a65b9"
  ) +
  coord_cartesian(xlim = c(143380000, 143410000)) +
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

ggsave("PAR_nonPAR_boundary.telomeric.pdf", p,
  width = 7,
  height = 7)

ggsave("PAR_nonPAR_boundary.telomeric.png", p,
  width = 7,
  height = 7)

```
# Make a plot of overlapping reads centromeric to the PAR
```
library(Rsamtools)
library(GenomicAlignments)
library(ggplot2)

bam <- "Poecile_atricapillus.hifi_vs_hap1.bam"

which <- GRanges(
  "NC_081289.1",
  IRanges(143740000, 143770000)
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
    xintercept = 143747811,
    linetype = "dashed",
    linewidth = 0.7,
    color = "#8a65b9"
  ) +
  coord_cartesian(xlim = c(143740000, 143770000)) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "NC_081289.1 position (bp)",
    y = "HiFi read",
    color = "MAPQ",
    title = "PacBio HiFi alignments across NC_081289.1:650-700 kb"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave("PAR_nonPAR_boundary.centromeric.pdf", p,
  width = 7,
  height = 7)

ggsave("PAR_nonPAR_boundary.centromeric.png", p,
  width = 7,
  height = 7)

```