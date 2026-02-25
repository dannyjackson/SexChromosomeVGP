#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GENESPACE)
  library(ggplot2)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript genespace.synteny_to_chicken.r <CHR_TYPE> <SPECIES> <OUTDIR>\n",
       "Example: Rscript ... X Bos_taurus /path/out/")
}

CHR_TYPE <- args[1]
SPECIES_ID <- args[2]
OUTDIR <- args[3]

# Working directory: prefer env var WORKING_DIR, otherwise current dir (where you cd'd in SLURM)
wd <- Sys.getenv("WORKING_DIR", unset = getwd())

message("CHR_TYPE: ", CHR_TYPE)
message("SPECIES:  ", SPECIES_ID)
message("OUTDIR:   ", OUTDIR)
message("WD:       ", wd)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = "/vf/users/jacksondan/conda/envs/genespace/bin/"
)

out <- run_genespace(gpar, overwrite = FALSE)

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


ripDat <- plot_riparian(
  gsParam = out,
  refGenome = "Gallus_gallus_REF",
  genomeIDs = genomeIDs,
  refChrOrdFun = ordFun,
  minChrLen2plot = 1,
  addThemes = ggthemes,
  chrFill = "lightgrey",
  reorderBySynteny = FALSE,

)

p_list <- ripDat$plot
p <- p_list[[1]]

pdf_name <- file.path(OUTDIR, sprintf("%s.Gallus_gallus.5x3.%s.pdf", SPECIES_ID, CHR_TYPE))

ggsave(
  filename = pdf_name,
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

message("Wrote: ", pdf_name)

pdf_name <- file.path(OUTDIR, sprintf("%s.Gallus_gallus.5x3.%s.png", SPECIES_ID, CHR_TYPE))

ggsave(
  filename = pdf_name,
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

message("Wrote: ", png_name)