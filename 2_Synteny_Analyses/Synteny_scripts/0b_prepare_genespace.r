args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
repo <- args[2]
stopifnot(!is.na(outdir), nzchar(outdir))

library(GENESPACE)

SPECIES <- sort(basename(list.dirs(repo, full.names = TRUE, recursive = FALSE)))

parsedPaths <- parse_annotations(
  rawGenomeRepo = repo,
  genomeDirs    = SPECIES,
  genomeIDs     = SPECIES,
  gffString     = "gff",
  faString      = "translated.cds",
  presets       = "ncbi",
  genespaceWd   = outdir
)