#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Gviz)
})

options(ucscChromosomeNames = FALSE)

base <- "/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

parse_region <- function(x) {
  m <- regexec("^([^,]+),([^:]+):(\\d+)-(\\d+)$", x)
  hit <- regmatches(x, m)[[1]]

  if (length(hit) != 5) {
    stop(
      "Region must look like: Species,seqid:start-end\n",
      "Example: Taeniopygia_guttata,NC_133063.1:0-2824462"
    )
  }

  list(
    species = hit[2],
    seqid = hit[3],
    start = as.integer(hit[4]),
    end = as.integer(hit[5])
  )
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    "Usage:\n",
    "  Rscript plot_region.R 'Taeniopygia_guttata,NC_133063.1:0-2824462'\n"
  )
}

region_string <- args[1]
region <- parse_region(region_string)

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

message("Importing GFF region from: ", gff_file)

features <- import(
  gff_file,
  format = "gff",
  which = gr_region
)

features <- GenomeInfoDb::keepSeqlevels(
  features,
  seqid,
  pruning.mode = "coarse"
)

if (length(features) == 0) {
  stop("No GFF features found in this region.")
}

# Choose feature types to show
plot_types <- c("gene")

features <- features[as.character(features$type) %in% plot_types]

if (length(features) == 0) {
  stop("No selected feature types found in this region.")
}

# Clean feature type
features$type <- as.character(features$type)
features$type[is.na(features$type) | features$type == ""] <- "feature"

# Build clean labels
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

genome_axis <- GenomeAxisTrack(name = "Position")

annotation_track <- AnnotationTrack(
  range = features,
  name = "Annotation",
  chromosome = seqid,
  genome = species,
  feature = features$type,
  id = features$plot_id,
  stacking = "dense",
  shape = "box",
  showFeatureId = TRUE,
  fontcolor.feature = "black",
  cex.feature = 0.55,
  background.title = "gray90",
  col.title = "black",
  col.axis = "black",
  legend = TRUE
)

out_file <- paste0(
  species, "_",
  gsub("[^A-Za-z0-9_.-]", "_", seqid), "_",
  region_start, "_", region_end,
  "_features.pdf"
)

pdf(out_file, width = 16, height = 3)

plotTracks(
  list(genome_axis, annotation_track),
  from = region_start + 1,
  to = region_end,
  chromosome = seqid,
  main = paste0(species, " | ", seqid, ":", region_start, "-", region_end),
  cex.title = 0.9,
  cex.axis = 0.8,
  sizes = c(0.12, 0.88)
)

dev.off()

message("Wrote plot: ", out_file)