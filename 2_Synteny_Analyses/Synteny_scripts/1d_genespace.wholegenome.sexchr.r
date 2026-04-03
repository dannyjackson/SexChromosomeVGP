#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GENESPACE)
  library(ggplot2)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript genespace.wholegenome.sexchr.r <SPECIES_ID> <OUTDIR> <SEXCHR>\n",
       "Example: Rscript genespace.wholegenome.sexchr.r Bos_taurus /path/out/ X")
}

SPECIES_ID <- args[1]
OUTDIR <- args[2]
SEXCHR <- args[3]
REF <- args[4]

# Working directory: prefer env var WORKING_DIR, otherwise current dir (where you cd'd in SLURM)

f <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  REF, "sexshared",
  SPECIES_ID,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

# Genome IDs to plot
genomeIDs <- c(SPECIES_ID, REF)

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
    refGenome = REF,
    genomeIDs = genomeIDs,
    minChrLen2plot = 1,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    reorderBySynteny = FALSE
)

p_list <- ripDat$plot
p <- p_list[[1]]

pdf_name <- file.path(OUTDIR, sprintf("%s.%s.5x3.%s.pdf", SPECIES_ID, REF, SEXCHR))

ggsave(
  filename = pdf_name,
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

message("Wrote: ", pdf_name)

png_name <- file.path(OUTDIR, sprintf("%s.%s.5x3.%s.png", SPECIES_ID, REF, SEXCHR))

ggsave(
  filename = png_name,
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

message("Wrote: ", png_name)