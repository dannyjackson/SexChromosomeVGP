# Podarcis neo sex
Is there a fusion on the Podarcis bocagei 

/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexlimited/usable_bed/Podarcis_bocagei.bed

source myconda
mamba activate genespace
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Podarcis

mkdir bed peptide
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexlimited/usable_bed/Podarcis_bocagei.bed Podarcis_bocagei/bed/

cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexlimited/usable_peptide/Podarcis_bocagei.fa Podarcis_bocagei/peptide/

cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexlimited/usable_bed/Podarcis_vaucheri.bed Podarcis_vaucheri/bed/

cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexlimited/usable_peptide/Podarcis_vaucheri.fa Podarcis_vaucheri/peptide/

```
#!/bin/bash
#SBATCH --job-name=genespace_podarcis
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

source myconda
mamba activate genespace_py3.10

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/scripts

Rscript genespace.podarcis.R
```
```
library(GENESPACE)
library(ggplot2)
library(data.table)

wd="/vf/users/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Podarcis"

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = "/vf/users/jacksondan/conda/envs/genespace/bin/"
)

out <- run_genespace(gpar, overwrite = FALSE)

```
#!/usr/bin/env Rscript

library(GENESPACE)
library(ggplot2)
library(data.table)


# Working directory: prefer env var WORKING_DIR, otherwise current dir (where you cd'd in SLURM)

f <- file.path(
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

# Genome IDs to plot
```
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
)

SPECIES_ID="Podarcis_bocagei"
SEXCHR="OZ076891.1"

roi <- data.frame( 
    genome = SPECIES_ID, 
    chr = SEXCHR, 
    color = c("#c23d3d")) 


ripDat <- plot_riparian(
    gsParam = gsParam,
    highlightBed = roi, 
    backgroundColor = NULL,
    minChrLen2plot = 1,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    reorderBySynteny = FALSE
)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Podarcis_bocagei_W.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)
```
# explore with lastz
```
samtools faidx /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/rPodBoc1/rPodBoc1.fa OZ076891.1 > PodBoc.W.fa

#!/bin/bash
#SBATCH --job-name=genespace_podarcis
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

module load lastz

cd /vf/users/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Podarcis/lastz

lastz PodBoc.W.fa /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/rPodVau1/rPodVau1.fa \
  --format=rdotplot --strand=plus \
  > PodBoc_W.PodVau_all.lastz.tsv

samtools faidx rPodBoc1.fa OZ076891.1 > PodBoc.W.fa

lastz PodBoc.W.fa[multiple] /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/rPodVau1/rPodVau1.fa[multiple] \
  --format=rdotplot \
  --strand=both \
  > PodBoc_W.PodVau_all.rdotplot
```
# Explore with minimap
```
cd /vf/users/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Podarcis/minimap

module load minimap2/2.30
module load seqkit

seqkit grep -r -v -p 'scaffold' \
  /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/rPodBoc1/rPodBoc1.fa \
  > rPodBoc1.no_scaffolds.fa

seqkit grep -r -v -p 'scaffold' \
  /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/rPodVau1/rPodVau1.fa \
  > rPodVau1.no_scaffolds.fa

seqkit grep -r -v -p 'scaffold' \
  /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Gallus_gallus/Gallus_gallus.fna \
  > Gallus_gallus.no_scaffolds.fa

seqkit grep -r -v -p 'scaffold' \
  /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Anolis_sagrei/ncbi_dataset/data/GCF_037176765.1/GCF_037176765.1_rAnoSag1.mat_genomic.fna \
  > Anolis_sagrei.no_scaffolds.fa

#!/bin/bash
#SBATCH --job-name=minimap_podarcis
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/minimap2

module load minimap2/2.30

minimap2 \
  -x asm20 \
  --eqx \
  -k 23 \
  -t 12 \
  -m 200 \
  -n 5 \
  -w 19 \
  rPodBoc1.no_scaffolds.fa \
  Anolis_sagrei.no_scaffolds.fa \
  > PodBoc_AnoSag.paf

minimap2 \
  -x asm20 \
  --eqx \
  -k 23 \
  -t 12 \
  -m 200 \
  -n 5 \
  -w 19 \
  rPodBoc1.no_scaffolds.fa \
  Gallus_gallus.no_scaffolds.fa \
  > PodBoc_GalGal.paf

minimap2 \
  -x asm20 \
  --eqx \
  -k 23 \
  -t 12 \
  -m 200 \
  -n 5 \
  -w 19 \
  rPodBoc1.no_scaffolds.fa \
  rPodVau1.no_scaffolds.fa \
  > PodBoc_PodVau.paf

git clone https://github.com/moold/paf2dotplot

./paf2dotplot/paf2dotplot.r \
  -o PodBoc_AnoSag \
  -q 1 \
  -m 1 \
  -r 1 \
  -p 10 \
  PodBoc_AnoSag.paf

./paf2dotplot/paf2dotplot.r \
  -o PodBoc_PodVau \
  -q 1 \
  -m 1 \
  -r 1 \
  -p 10 \
  PodBoc_PodVau.paf

./paf2dotplot/paf2dotplot.r \
  -o PodBoc_GalGal \
  -q 1 \
  -m 1 \
  -r 1 \
  -p 10 \
  PodBoc_GalGal.paf

# Filter and write BED:
# Output BED columns:
# ref_id  ref_start  ref_end  query_id  query_start  query_end  identity_percent  ref_span_len
awk -v thr=0.985 -v minlen=10000 -v OFS='\t' '
  function ident_from_de(de){ return (de=="" ? -1 : 1.0 - de) }
  function ident_from_core(m,a){ return (a>0 ? m/a : -1) }
  {
    # PAF fields:
    # 1 qname 2 qlen 3 qstart 4 qend 5 strand 6 tname 7 tlen 8 tstart 9 tend 10 nmatch 11 alnlen 12 mapq
    tspan = $9 - $8
    de=""; for(i=13;i<=NF;i++) if($i ~ /^de:f:/){ split($i,x,":"); de=x[3]; break }
    id = ident_from_de(de); if (id < 0) id = ident_from_core($10,$11)

    if (id >= thr && tspan >= minlen)
      printf "%s\t%d\t%d\t%s\t%d\t%d\t%.4f\t%d\n", $6, $8, $9, $1, $3, $4, id*100, tspan
  }
' PodBoc_W.PodVau_all.paf > PodBoc_W.PodVau_all.bed

```
# Explore with lastz properly
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/lastz

module load lastz/1.04.22

lastz \
  /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/rPodBoc1/rPodBoc1.fa[multiple] \
  /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/rPodVau1/rPodVau1.fa \
  --strand=both \
  --notransition \
  --step=50 \
  --seed=match12 \
  --filter=identity:70..100 \
  --format=maf \
  > PodBoc_vs_PodVau.maf

cp /vf/users/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Podarcis/lastz/PodBoc.W.fa .

lastz PodBoc.W.fa \
  --format=rdotplot \
  --strand=plus \
  --self \
  --nomirror \
  > PodBoc_W.self.lastz.tsv
```
# one chr plot
```
R

# Load libraries
library(ggplot2)
library(readr)
library(dplyr)
# Read & clean (robust to header/NA rows in rdotplot)
dots_raw <- read.table("PodBoc_W.self.lastz.tsv", header = TRUE, sep = "", stringsAsFactors = FALSE)
x <- suppressWarnings(as.numeric(dots_raw[[1]]))
y <- suppressWarnings(as.numeric(dots_raw[[2]]))

dots_raw <- read.table(
  "PodBoc_W.self.lastz.tsv",
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE,
  fill = TRUE
)

# Plot
png("self.long.lastz.png", width = 1800, height = 1800, res = 300)
par(mar = c(4,4,2,1))
xlim <- range(x, na.rm = TRUE)
ylim <- range(y, na.rm = TRUE)
plot(x, y, type = "l", xlab = "Reference", ylab = "Query", xlim = xlim, ylim = ylim, asp = 1)
dev.off()
```
# multiple chr plot
```
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

dat <- read_tsv(
  "PodBoc_W_vs_PodVau.tsv",
  col_names = c("V1", "V2"),
  col_types = cols(.default = col_character())
)

dat2 <- dat %>%
  mutate(
    is_header = !is.na(V1) & V1 == "seq1",
    chr = case_when(
      is_header & !is.na(V2) & !str_detect(V2, "^CB") ~ V2,
      is_header & !is.na(V2) &  str_detect(V2, "^CB") ~ "DROP",
      TRUE ~ NA_character_
    ),
    x = suppressWarnings(as.numeric(V1)),
    y = suppressWarnings(as.numeric(V2))
  ) %>%
  fill(chr, .direction = "down") %>%
  filter(!is.na(chr), chr != "DROP") %>%
  filter(!is.na(x), !is.na(y))

p <- ggplot(dat2, aes(x, y)) +
  geom_point(size = 0.001, alpha = 0.5) +
  facet_wrap(~ chr, scales = "free") +
  theme_bw()

ggsave(
  filename = "PodBoc_PodVau.lastz.png",
  plot = p,
  width = 10,
  height = 10,
  units = "in"
)

p <- ggplot(dat2, aes(x, y)) +
  geom_path(linewidth = 0.01) +
  facet_wrap(~ chr, scales = "free") +
  theme_bw()

ggsave(
  filename = "PodBoc_PodVau.lastz.path.png",
  plot = p,
  width = 10,
  height = 10,
  units = "in"
)
```
# Explore with cactus
```
module load cactus/3.0.0


cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/cactus

cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/minimap2/rPodVau1.no_scaffolds.fa .
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/minimap2/rPodBoc1.no_scaffolds.fa .
bgzip rPodVau1.no_scaffolds.fa
bgzip rPodBoc1.no_scaffolds.fa

echo '(rPodVau1, rPodBoc1);' > cactus.podarcis.txt
echo -e 'rPodVau1\t/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/cactus/rPodVau1.no_scaffolds.fa.gz' >> cactus.podarcis.txt
echo -e 'rPodBoc1\t/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/cactus/rPodBoc1.no_scaffolds.fa.gz' >> cactus.podarcis.txt

#!/bin/bash
#SBATCH --job-name=cactus_podarcis
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

module load cactus/3.0.0

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/cactus

cactus \
  /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/cactus/jobstore/ \
  /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/cactus/cactus.podarcis.txt \
  /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/cactus/cactus.podarcis.output.hal
