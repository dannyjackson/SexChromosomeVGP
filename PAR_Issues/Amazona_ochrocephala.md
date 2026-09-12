# Amazona ochrocephala PAR inversion investigation

### Set up the working environment
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues
mkdir -p Amazona_ochrocephala
cd Amazona_ochrocephala
```
### Download the fastq file
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Amazona_ochrocephala

module load aws
aws configure set default.s3.max_concurrent_requests 50
aws s3 cp \
  s3://genomeark/species/Amazona_ochrocephala/bAmaOch1/genomic_data/pacbio_hifi/ \
  . \
  --recursive \
  --exclude "*" \
  --include "*.hifi_reads.fastq.gz" \
  --no-sign-request
```
### Align and index
```
#!/bin/bash
#SBATCH --job-name=align_Amazona_ochrocephala
#SBATCH --output=slurm_output/%x_%A.out
#SBATCH --error=slurm_output/%x_%A.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

module load minimap2 samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Amazona_ochrocephala

FASTA=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Amazona_ochrocephala/ncbi_dataset/data/GCA_039720435.1/GCA_039720435.1_bAmaOch1.hap1_genomic.fna
FASTQ1=m54306Ue_230330_053511.bc1002--bc1002.hifi_reads.fastq.gz
FASTQ2=m64055e_230122_092023.bc1002--bc1002.hifi_reads.fastq.gz
FASTQ3=m64334e_230205_093940.bc1002--bc1002.hifi_reads.fastq.gz

minimap2 -ax map-hifi --secondary=no -t 24 "$FASTA" "$FASTQ1" "$FASTQ2" "$FASTQ3" | samtools sort -@ 8 -m 4G -T "Amazona_ochrocephala_sort" -o "Amazona_ochrocephala.hifi_vs_hap1.bam"

samtools index -@ 8 Amazona_ochrocephala.hifi_vs_hap1.bam
```
# Filter for mapq > 30
```
samtools view -@ 8 -b -q 30 \
    Amazona_ochrocephala.hifi_vs_hap1.bam \
    -o Amazona_ochrocephala.hifi_vs_hap1.mapq30.bam

samtools index -@ 8 Amazona_ochrocephala.hifi_vs_hap1.mapq30.bam
```

PAR: 	CM077909 119902688-120227954
ROI: 	CM077909 116902688-119902688
ROI: 	NC_081289.1	143500000-143407603

### make a plot of overlapping reads
```
library(Rsamtools)
library(GenomicAlignments)
library(ggplot2)

bam <- "Amazona_ochrocephala.hifi_vs_hap1.mapq30.bam"

which <- GRanges(
  "CM077909.1",
  IRanges(118000000, 119950000)
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
    xintercept = 119902688,
    linetype = "dashed",
    linewidth = 0.2,
    color = "#8a65b9"
  ) +
  coord_cartesian(xlim = c(118000000, 119950000)) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "X/Z position (bp)",
    y = "HiFi read",
    color = "MAPQ",
    title = "PacBio HiFi alignments across X/Z kb"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave("PAR_nonPAR_boundary.pdf", p,
  width = 12,
  height = 5)

ggsave("PAR_nonPAR_boundary.png", p,
  width = 12,
  height = 5)

```