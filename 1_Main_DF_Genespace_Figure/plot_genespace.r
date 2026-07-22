#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ape)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop(
    "Usage: Rscript plot_genespace.r <FOCAL_SPECIES> [species.txt]\n",
    "Example: Rscript plot_genespace.r Anolis_sagrei species.txt"
  )
}

FOCAL <- args[1]
species_file <- ifelse(length(args) >= 2, args[2], "species.txt")

BASE_DIR <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace"
IN_DIR   <- file.path(BASE_DIR, FOCAL, "combined_phaseblks")
OUT_DIR <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/VGP_DF_Figure"
infile <- file.path(IN_DIR, paste0(FOCAL, "_syntentic_chromosomes.sex_shared.csv"))
if (!file.exists(infile)) stop("Input CSV not found: ", infile)

if (!file.exists(species_file)) stop("Species file not found: ", species_file)

out_pdf <- file.path(OUT_DIR, paste0(FOCAL, "_tree.sexshared.RefChrs.species_of_interest.pdf"))

# --- read species of interest --------------------------------------------------

species_interest <- readLines(species_file, warn = FALSE)
species_interest <- str_trim(species_interest)
species_interest <- species_interest[species_interest != ""]
species_interest <- species_interest[!grepl("^Species:*$", species_interest, ignore.case = TRUE)]

# Convert binomials like "Carcharodon carcharias" to tree/CSV style:
# "Carcharodon_carcharias"
species_interest <- gsub("\\s+", "_", species_interest)

if (length(species_interest) == 0) {
  stop("No species found in species file: ", species_file)
}

# --- read input ----------------------------------------------------------------

out <- read.csv(infile, stringsAsFactors = FALSE)

if (!"Species" %in% names(out)) {
  stop("Input CSV must contain a column named Species")
}

# Normalize Species column just in case it contains spaces
out$Species <- gsub("\\s+", "_", str_trim(out$Species))

# Filter CSV to only species of interest
out <- out %>%
  filter(Species %in% species_interest)

if (nrow(out) == 0) {
  stop(
    "None of the species in ", species_file,
    " were found in the input CSV: ", infile
  )
}

tree_file <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/roadies_v1.1.16b.numbers.scientific.nwk"
if (!file.exists(tree_file)) stop("Tree file not found: ", tree_file)

tr_full <- ape::read.tree(tree_file)

# --- reconcile tip labels to match Species names in CSV ------------------------

# 1) rename Accipiter_gentilis -> Astur_gentilis
tr_full$tip.label <- ifelse(
  tr_full$tip.label == "Accipiter_gentilis",
  "Astur_gentilis",
  tr_full$tip.label
)

# 2) drop ONE trailing underscore, e.g. "Pongo_abelii_" -> "Pongo_abelii"
tr_full$tip.label <- sub("_$", "", tr_full$tip.label)

# 3) prune Canis lupus subspecies; retain only one node named "Canis_lupus"
canis_subspp <- grep("^Canis_lupus_", tr_full$tip.label, value = TRUE)

if (length(canis_subspp) > 0) {
  keep_tip <- if ("Canis_lupus_baileyi" %in% canis_subspp) {
    "Canis_lupus_baileyi"
  } else {
    canis_subspp[1]
  }

  drop_tips <- setdiff(canis_subspp, keep_tip)
  tr_full <- drop.tip(tr_full, drop_tips)

  tr_full$tip.label[tr_full$tip.label == keep_tip] <- "Canis_lupus"
}

# Safety: enforce unique tip labels after edits
if (anyDuplicated(tr_full$tip.label)) {
  dup <- tr_full$tip.label[duplicated(tr_full$tip.label)]
  stop("Duplicate tip labels after reconciliation: ", paste(unique(dup), collapse = ", "))
}

# --- prune tree to species of interest -----------------------------------------

species_in_tree <- intersect(species_interest, tr_full$tip.label)
species_missing_tree <- setdiff(species_interest, tr_full$tip.label)
species_missing_csv  <- setdiff(species_interest, unique(out$Species))

if (length(species_in_tree) == 0) {
  stop("None of the species in ", species_file, " were found in the tree.")
}

if (length(species_missing_tree) > 0) {
  warning(
    "Species in species file but missing from tree: ",
    paste(species_missing_tree, collapse = ", ")
  )
}

if (length(species_missing_csv) > 0) {
  warning(
    "Species in species file but missing from CSV: ",
    paste(species_missing_csv, collapse = ", ")
  )
}

tr_pruned <- keep.tip(tr_full, species_in_tree)

# Make ultrametric
tr_ultra <- chronos(tr_pruned, lambda = 1)

# ---- helper: split a cell like "22,Z" into character vector c("22","Z") -------

split_chr <- function(x) {
  if (is.null(x) || length(x) == 0) return(character())
  x <- as.character(x)
  if (is.na(x) || x == "") return(character())
  parts <- unlist(strsplit(x, ",", fixed = TRUE))
  parts <- str_trim(parts)
  parts[parts != ""]
}

# ---- Build long format: Species x RefChr --------------------------------------

out_long <- out %>%
  mutate(
    both_list = map(RefChr_both, split_chr),
    c1_list   = map(RefChr_chr1_only, split_chr),
    c2_list   = map(RefChr_chr2_only, split_chr),
    ref_all   = map2(map2(both_list, c1_list, union), c2_list, union)
  ) %>%
  select(Species, ref_all) %>%
  unnest(ref_all) %>%
  filter(!is.na(ref_all), ref_all != "") %>%
  distinct(Species, ref_all) %>%
  rename(RefChr = ref_all) %>%
  distinct(Species, RefChr)

if (nrow(out_long) == 0) {
  stop("No RefChr data remained after filtering to species of interest.")
}

# Keep only data for species also retained in the pruned tree
out_long <- out_long %>%
  filter(Species %in% tr_ultra$tip.label)

# Normalize chromosome labels to match c("1", ..., "39", "Z")
out_long <- out_long %>%
  mutate(
    RefChr = str_trim(as.character(RefChr)),
    RefChr = sub("^chr", "", RefChr, ignore.case = TRUE),
    RefChr = sub("^GGA", "", RefChr, ignore.case = TRUE)
  )

# Chicken reference chromosomes: autosomes 1-39 plus Z
chr_levels <- c(as.character(1:39), "Z")

# full grid: species of interest x all chromosomes
grid_df <- tidyr::expand_grid(
  Species = tr_ultra$tip.label,
  RefChr = chr_levels
) %>%
  left_join(out_long %>% mutate(present = 1L), by = c("Species", "RefChr")) %>%
  mutate(present = if_else(is.na(present), 0L, present))

# ---- Plot ---------------------------------------------------------------------

pdf(out_pdf, width = 6, height = 4)

# extra top/right margin for chromosome labels and grid
par(mar = c(2, 2, 7, 25), xpd = NA)

tmp <- plot(tr_ultra, plot = FALSE)
nchr <- length(chr_levels)

# Use a fixed grid spacing that is wide enough for all 40 chicken chromosomes
tree_width <- diff(range(tmp$x.lim))
dx <- tree_width * 0.035

x0 <- max(tmp$x.lim) + dx * 4
x_chr <- x0 + (seq_along(chr_levels) - 1) * dx
x_max <- max(x_chr) + dx * 4

plot(
  tr_ultra,
  cex = 0.35,
  edge.width = 0.5,
  no.margin = TRUE,
  x.lim = c(min(tmp$x.lim), x_max),
  y.lim = c(1, Ntip(tr_ultra) + 4)
)

lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_y <- lp$yy[1:Ntip(tr_ultra)]
tips  <- tr_ultra$tip.label

y_min <- min(tip_y)
y_max <- max(tip_y)

# Draw every chromosome column first, even if there are zero hits
for (j in seq_along(chr_levels)) {
  chrj <- chr_levels[j]
  xj <- x_chr[j]

  # faint vertical guide so empty columns remain visible
  segments(
    x0 = xj,
    y0 = y_min,
    x1 = xj,
    y1 = y_max,
    col = "gray90",
    lwd = 0.4
  )

  present_species <- grid_df %>%
    filter(RefChr == chrj, present == 1L) %>%
    pull(Species)

  is_present <- tips %in% present_species

  # empty/background dots for every species
  points(
    rep(xj, length(tip_y)),
    tip_y,
    pch = 21,
    cex = 0.45,
    col = "gray70",
    bg = "white",
    lwd = 0.4
  )

  # filled dots where chromosome is present
  if (any(is_present)) {
    points(
      rep(xj, sum(is_present)),
      tip_y[is_present],
      pch = 21,
      cex = 0.45,
      col = "black",
      bg = "black",
      lwd = 0.4
    )
  }
}

# Column labels above every grid column
y_lab <- Ntip(tr_ultra) + 2.2

text(
  x = x_chr,
  y = y_lab,
  labels = chr_levels,
  srt = 90,
  adj = c(0, 0.5),
  cex = 0.45,
  xpd = NA
)

# Optional title for the chromosome grid
text(
  x = mean(range(x_chr)),
  y = Ntip(tr_ultra) + 3.5,
  labels = "Chicken chromosome",
  cex = 0.45,
  xpd = NA
)

dev.off()

cat("FOCAL:", FOCAL, "\n")
cat("Species file:", species_file, "\n")
cat("Read:", infile, "\n")
cat("Wrote:", out_pdf, "\n")
cat("Species requested:", length(species_interest), "\n")
cat("Species retained in tree:", Ntip(tr_ultra), "\n")
cat("Species with data:", length(unique(out_long$Species)), "\n")
cat("Chicken chromosomes plotted:", length(chr_levels), "\n")