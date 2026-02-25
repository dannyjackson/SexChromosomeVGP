#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GENESPACE)
  library(ggplot2)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript genespace.synteny_to_chicken.wholegenome.sexchr.r <SPECIES_ID> <OUTDIR> <SEXCHR>\n",
       "Example: Rscript genespace.synteny_to_chicken.wholegenome.sexchr.r Bos_taurus /path/out/ X")
}

SPECIES_ID <- args[1]
OUTDIR <- args[2]
SEXCHR <- args[3]

# Working directory: prefer env var WORKING_DIR, otherwise current dir (where you cd'd in SLURM)

f <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace",
  SPECIES_ID,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

# Genome IDs to plot
genomeIDs <- c(SPECIES_ID, "Gallus_gallus_REF")

ordFun <- function(x) {
  data.table::frank(
    list(grepl("[a-zA-Z]", x),
         suppressWarnings(as.numeric(gsub("[a-zA-Z]", "", x)))),
    ties.method = "random"
  )
}

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
)

roi <- data.frame( 
    genome = SPECIES_ID, 
    chr = SEXCHR, 
    color = c("#c23d3d")) 


ripDat <- plot_riparian(
    gsParam = gsParam,
    highlightBed = roi, 
    backgroundColor = NULL,
    refGenome = "Gallus_gallus_REF",
    genomeIDs = genomeIDs,
    minChrLen2plot = 1,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    reorderBySynteny = FALSE
)

p_list <- ripDat$plot
p <- p_list[[1]]

# pdf_name <- file.path(OUTDIR, sprintf("%s.Gallus_gallus.5x3.%s.38.pdf", SPECIES_ID, SEXCHR))
pdf_name <- "test.pdf"

ggsave(
  filename = pdf_name,
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

message("Wrote: ", pdf_name)

png_name <- file.path(OUTDIR, sprintf("%s.Gallus_gallus.5x3.%s.png", SPECIES_ID, SEXCHR))

ggsave(
  filename = png_name,
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

message("Wrote: ", png_name)