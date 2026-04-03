#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GENESPACE)
  library(ggplot2)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript genespace.wholegenome.r <SPECIES> <OUTDIR>\n",
       "Example: Rscript ... X Bos_taurus /path/out/")
}

SPECIES_ID <- args[1]
OUTDIR <- args[2]
REF <- args[3]
NUMCORES <- as.integer(Sys.getenv("OMP_NUM_THREADS", "16"))

# Working directory: prefer env var WORKING_DIR, otherwise current dir (where you cd'd in SLURM)
wd <- Sys.getenv("WORKING_DIR", unset = getwd())

message("SPECIES:  ", SPECIES_ID)
message("OUTDIR:   ", OUTDIR)
message("WD:       ", wd)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

gpar <- init_genespace(
  wd = wd,
  nCores = NUMCORES,
  nGaps = 100,
  path2mcscanx = "/vf/users/jacksondan/conda/envs/genespace/bin/"
)

out <- run_genespace(gpar, overwrite = FALSE)

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


ripDat <- plot_riparian(
  gsParam = out,
  refGenome = REF,
  genomeIDs = genomeIDs,
  refChrOrdFun = ordFun,
  addThemes = ggthemes,
  chrFill = "lightgrey",
  reorderBySynteny = FALSE
)

p_list <- ripDat$plot
p <- p_list[[1]]

pdf_name <- file.path(OUTDIR, sprintf("%s.%s.5x3.wholegenome.pdf", SPECIES_ID, REF))

ggsave(
  filename = pdf_name,
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

message("Wrote: ", pdf_name)

png_name <- file.path(OUTDIR, sprintf("%s.%s.5x3.wholegenome.png", SPECIES_ID, REF))

ggsave(
  filename = pdf_name,
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

message("Wrote: ", png_name)