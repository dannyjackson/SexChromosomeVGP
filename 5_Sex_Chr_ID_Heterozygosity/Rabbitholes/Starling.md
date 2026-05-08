GCA_052056855.1

source myconda
mamba activate ncbi_datasets

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/VGP_Phase2

FILES_TO_DOWNLOAD="gff3,rna,cds,protein,genome,seq-report"
ACCESSION=GCA_052056855.1
datasets download genome accession "${ACCESSION}" \
  --include "${FILES_TO_DOWNLOAD}" 

unzip ncbi_dataset.zip

REF_GENOME="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/VGP_Phase2/ncbi_dataset/data/GCA_052056855.1/GCA_052056855.1_bStuVul1.hap1_genomic.fna"

module load samtools
bgzip $REF_GENOME
samtools faidx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/VGP_Phase2/ncbi_dataset/data/GCA_052056855.1/GCA_052056855.1_bStuVul1.hap1_genomic.fna.gz 

samtools faidx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/VGP_Phase2/ncbi_dataset/data/GCA_052056855.1/GCA_052056855.1_bStuVul1.hap1_genomic.fna.gz CM125341.1 > W.chrom.fa

samtools faidx ${REF_GENOME}.gz CM125342.1 > Z.chrom.fa

# Align W to Z

# Align (PAF) + include CIGAR in 'cg:Z:' tag via -c (important for SVbyEye readPaf/plotMiro)
# Here: query = Morphnus chr, target = Gallus chr

module load minimap2

THREADS=16

minimap2 -x asm5 -t "${THREADS}" --secondary=no -c --cs=long \
  W.chrom.fa Z.chrom.fa \
  > "starling.W_to_Z.paf"

Rscript plot_svbyeye.R \
  "starling.W_to_Z.paf" \
  "CM125341.1" \
  "CM125342.1" \
  "starling.CM125341_to_CM125342.pdf"

Rscript plot_svbyeye.inverted.R \
  "starling.W_to_Z.paf" \
  "CM125341.1" \
  "CM125342.1" \
  "starling.CM125341_to_CM125342.inverted.pdf"

```
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: plot_svbyeye.R <paf> <target_name> <query_name> <out_pdf>")
}

paf_file   <- args[1]
target_nm  <- args[2]  # e.g. NC_052577.1 (top)
query_nm   <- args[3]  # e.g. CM098430.1 (bottom)
out_pdf    <- args[4]

suppressPackageStartupMessages({
  library(SVbyEye)
  library(ggplot2)
})

# Read minimap2 PAF (includes tags like cg if present)
paf <- readPaf(paf.file = paf_file, include.paf.tags = FALSE)


# Miropeats-style visualization (polygons), colored by direction by default
p <- plotMiro(
  paf.table = paf,
  color.by = "direction",
  add.alignment.arrows = TRUE,
  outline.alignments = FALSE,
  offset.alignments = FALSE
) +
  ggtitle(paste0(query_nm, " (query) vs ", target_nm, " (target)")) +
  theme_bw()

ggsave(out_pdf, plot = p, width = 14, height = 6, units = "in")
```
# inverted R
```
#!/usr/bin/env Rscript

Rscript plot_svbyeye.inverted.R \
  "starling.W_to_Z.paf" \
  "CM125341.1" \
  "CM125342.1" \
  "starling.CM125341_to_CM125342.inverted.pdf"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: plot_svbyeye.R <paf> <target_name> <query_name> <out_pdf>")
}

paf_file   <- "starling.W_to_Z.paf"
target_nm  <- "CM125341.1"
query_nm   <- "CM125342.1"
out_pdf    <- "starling.CM125341_to_CM125342.inverted.pdf"

suppressPackageStartupMessages({
  library(SVbyEye)
  library(ggplot2)
})

# Read minimap2 PAF (includes tags like cg if present)
paf <- readPaf(paf.file = paf_file, include.paf.tags = FALSE)

paf$t.start <- paf$t.len - paf$t.start
paf$t.end   <- paf$t.len - paf$t.end

tmp <- paf$t.start
paf$t.start <- paf$t.end
paf$t.end   <- tmp

# Miropeats-style visualization (polygons), colored by direction by default
p <- plotMiro(
  paf.table = paf,
  color.by = "direction",
  add.alignment.arrows = TRUE,
  outline.alignments = FALSE,
  offset.alignments = FALSE
) +
  ggtitle(paste0(query_nm, " (query) vs ", target_nm, " (target)")) +
  theme_bw()

ggsave(out_pdf, plot = p, width = 14, height = 6, units = "in")
```