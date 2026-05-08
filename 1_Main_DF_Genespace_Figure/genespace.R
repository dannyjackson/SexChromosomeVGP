#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GENESPACE)
  library(ggplot2)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

OUTDIR <- args[1]
NUMCORES <- as.integer(Sys.getenv("OMP_NUM_THREADS", "16"))

message("OUTDIR:   ", OUTDIR)


gpar <- init_genespace(
  wd = OUTDIR,
  nCores = NUMCORES,
  nGaps = 100,
  path2mcscanx = "/vf/users/jacksondan/conda/envs/genespace/bin/"
)

out <- run_genespace(gpar, overwrite = FALSE)