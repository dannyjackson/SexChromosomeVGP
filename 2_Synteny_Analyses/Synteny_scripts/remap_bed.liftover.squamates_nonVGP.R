#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
SPECIES <- args[1]

if (is.na(SPECIES) || !nzchar(SPECIES)) {
  stop("Usage: Rscript remap_bed.R <SPECIES>", call. = FALSE)
}

library(readr)

bed_file <- sprintf(
  "/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_nonVGP_lifted_gffs/Genespace_input_unfiltered/bed/%s.bed",
  SPECIES
)

map_file <- sprintf(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/%s.chromosome_mapping.tsv",
  SPECIES
)

out_file <- sprintf(
  "/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/squamate_nonVGP_lifted_gffs/Genespace_input_unfiltered/bed_reformatted/%s.remapped.bed",
  SPECIES
)

bed <- read_tsv(
  bed_file,
  col_names = FALSE,
  show_col_types = FALSE
)

map_tbl <- read_tsv(
  map_file,
  col_names = FALSE,
  show_col_types = FALSE
)

# names = accessions, values = chromosome labels
map_vec <- setNames(map_tbl$X2, map_tbl$X1)

# remap column 1, leave unmatched entries unchanged
bed$X1 <- ifelse(
  bed$X1 %in% names(map_vec),
  unname(map_vec[bed$X1]),
  bed$X1
)

write_tsv(
  bed,
  out_file,
  col_names = FALSE
)