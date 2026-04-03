args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
repo <- args[2]
stopifnot(!is.na(outdir), nzchar(outdir))

library(GENESPACE)

SPECIES <- sort(basename(list.dirs(repo, full.names = TRUE, recursive = FALSE)))

parsedPaths <- parse_annotations(
  rawGenomeRepo = repo,
  genomeDirs = SPECIES,
  genomeIDs = SPECIES,
  gffString = "gff$",
  faString = "translated.cds$",
  gffIdColumn = "ID",
  gffStripText = "^rna-",
  headerEntryIndex = 1,
  headerSep = " ",
  headerStripText = "^rna-|_[0-9]+$",
  genespaceWd = outdir,
)
