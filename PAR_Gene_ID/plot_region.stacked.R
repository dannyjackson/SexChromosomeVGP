#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Gviz)
})

options(ucscChromosomeNames = FALSE)

base <- "/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

plot_types <- c("gene")

parse_region <- function(x) {
  x <- trimws(x)

  m <- regexec("^([^,]+),([^:]+):(\\d+)-(\\d+)$", x)
  hit <- regmatches(x, m)[[1]]

  if (length(hit) != 5) {
    stop(
      "Region must look like: Species,seqid:start-end\n",
      "Bad line: ", x
    )
  }

  list(
    original = x,
    species = hit[2],
    seqid = hit[3],
    start = as.integer(hit[4]),
    end = as.integer(hit[5])
  )
}

make_plot_label <- function(region) {
  paste0(
    region$species,
    "\n",
    region$seqid,
    ":",
    region$start,
    "-",
    region$end
  )
}

load_region_features <- function(region) {
  species <- region$species
  seqid <- region$seqid
  region_start <- region$start
  region_end <- region$end

  species_dir <- file.path(base, species)
  gff_file <- file.path(species_dir, paste0(species, ".gff"))
  fna_file <- file.path(species_dir, paste0(species, ".fna"))

  if (!dir.exists(species_dir)) {
    stop("Species directory not found: ", species_dir)
  }

  if (!file.exists(gff_file)) {
    stop("GFF file not found: ", gff_file)
  }

  if (!file.exists(fna_file)) {
    warning("FNA file not found: ", fna_file)
  }

  gr_region <- GRanges(
    seqnames = seqid,
    ranges = IRanges(start = region_start + 1, end = region_end)
  )

  message("Importing GFF region: ", region$original)
  message("  from: ", gff_file)

  features <- import(
    gff_file,
    format = "gff",
    which = gr_region
  )

  if (length(features) == 0) {
    warning("No GFF features found in region: ", region$original)
    return(NULL)
  }

  features <- GenomeInfoDb::keepSeqlevels(
    features,
    seqid,
    pruning.mode = "coarse"
  )

  features <- features[as.character(features$type) %in% plot_types]

  if (length(features) == 0) {
    warning("No selected feature types found in region: ", region$original)
    return(NULL)
  }

  features$type <- as.character(features$type)
  features$type[is.na(features$type) | features$type == ""] <- "feature"

  label_col <- NULL
  for (candidate in c("Name", "gene", "gene_name", "product", "ID")) {
    if (candidate %in% colnames(mcols(features))) {
      label_col <- candidate
      break
    }
  }

  if (!is.null(label_col)) {
    features$plot_id <- as.character(mcols(features)[[label_col]])
  } else {
    features$plot_id <- as.character(features$type)
  }

  features$plot_id[is.na(features$plot_id) | features$plot_id == ""] <- features$type
  features$plot_id[is.na(features$plot_id) | features$plot_id == ""] <- "feature"

  features
}

make_region_track <- function(region, features) {
  AnnotationTrack(
    range = features,
    name = make_plot_label(region),
    chromosome = region$seqid,
    genome = region$species,
    feature = features$type,
    id = features$plot_id,
    stacking = "dense",
    shape = "box",
    showFeatureId = TRUE,
    fontcolor.feature = "black",
    cex.feature = 0.5,
    background.title = "gray90",
    col.title = "black",
    col.axis = "black",
    legend = TRUE
  )
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    "Usage:\n",
    "  Rscript plot_stacked_regions.R regions.txt [output.pdf]\n"
  )
}

input_file <- args[1]

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

out_file <- if (length(args) >= 2) {
  args[2]
} else {
  "stacked_regions.pdf"
}

region_lines <- readLines(input_file)
region_lines <- trimws(region_lines)
region_lines <- region_lines[region_lines != ""]
region_lines <- region_lines[!grepl("^#", region_lines)]

if (length(region_lines) == 0) {
  stop("No regions found in input file.")
}

regions <- lapply(region_lines, parse_region)

tracks <- list()

for (region in regions) {
  features <- load_region_features(region)

  if (is.null(features)) {
    next
  }

  tracks[[length(tracks) + 1]] <- make_region_track(region, features)
}

 if (length(tracks) == 0) {
  stop("No tracks could be generated.")
}

# One axis at the top.
genome_axis <- GenomeAxisTrack(name = "Position")

all_tracks <- c(list(genome_axis), tracks)

# Height scales with number of regions.
pdf_height <- max(4, 1.2 + length(tracks) * 2.2)

pdf(out_file, width = 16, height = pdf_height)

plotTracks(
  all_tracks,
  main = "Stacked genome region annotations",
  cex.title = 0.85,
  cex.axis = 0.75,
  sizes = c(0.08, rep(0.92 / length(tracks), length(tracks)))
)

dev.off()

message("Wrote plot: ", out_file)