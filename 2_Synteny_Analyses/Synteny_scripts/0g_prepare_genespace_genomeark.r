args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
repo <- args[2]
stopifnot(!is.na(outdir), nzchar(outdir))

SPECIES <- sort(basename(list.dirs(repo, full.names = TRUE, recursive = FALSE)))

library(GENESPACE)

parsedPaths2 <- parse_annotations(
  rawGenomeRepo = repo,
  genomeDirs    = SPECIES,
  genomeIDs     = SPECIES,
  gffString     = "gff",
  faString      = "cds",
  presets       = "none",
  genespaceWd   = outdir,
  gffIdColumn        = "protein_id",
  headerSep          = " ",
  headerEntryIndex   = 1,
  headerStripText    = ".*\\|",   # strip everything up to last '|'
  gffStripText       = ".*\\|"    # same normalization on protein_id values
)