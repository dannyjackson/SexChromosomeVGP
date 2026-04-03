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
  stop("Usage: Rscript plot_genespace.r <FOCAL_SPECIES>\nExample: Rscript plot_genespace.r Anolis_sagrei")
}
FOCAL <- args[1]

BASE_DIR <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace"
IN_DIR   <- file.path(BASE_DIR, FOCAL, "combined_phaseblks")

infile <- file.path(IN_DIR, paste0(FOCAL, "_syntentic_chromosomes.sex_shared.csv"))
if (!file.exists(infile)) stop("Input CSV not found: ", infile)

out_pdf <- file.path(IN_DIR, paste0(FOCAL, "_tree.sexshared.RefChrs.pdf"))

# --- read input
out <- read.csv(infile, stringsAsFactors = FALSE)

# If these columns were written as comma-separated strings, keep them as-is.
# We'll union them by splitting inside the mutate below.

tree_file <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/roadies_v1.1.16b.numbers.scientific.nwk"
if (!file.exists(tree_file)) stop("Tree file not found: ", tree_file)

tr_full <- ape::read.tree(tree_file)

# --- reconcile tip labels to match Species names in CSV ------------------------

# 1) rename Accipiter_gentilis -> Astur_gentilis
tr_full$tip.label <- ifelse(tr_full$tip.label == "Accipiter_gentilis",
                            "Astur_gentilis",
                            tr_full$tip.label)

# 2) drop ONE trailing underscore (e.g., "Pongo_abelii_" -> "Pongo_abelii")
tr_full$tip.label <- sub("_$", "", tr_full$tip.label)

# 3) prune Canis lupus subspecies; retain only one node named "Canis_lupus"
canis_subspp <- grep("^Canis_lupus_", tr_full$tip.label, value = TRUE)

if (length(canis_subspp) > 0) {
  # choose which one to keep as the representative tip
  keep_tip <- if ("Canis_lupus_baileyi" %in% canis_subspp) {
    "Canis_lupus_baileyi"
  } else {
    canis_subspp[1]
  }

  drop_tips <- setdiff(canis_subspp, keep_tip)
  tr_full <- drop.tip(tr_full, drop_tips)

  # rename the kept subspecies tip to the species-level name
  tr_full$tip.label[tr_full$tip.label == keep_tip] <- "Canis_lupus"
}

# (optional safety) enforce unique tip labels after edits
if (anyDuplicated(tr_full$tip.label)) {
  dup <- tr_full$tip.label[duplicated(tr_full$tip.label)]
  stop("Duplicate tip labels after reconciliation: ", paste(unique(dup), collapse = ", "))
}

# Now do chronos on the reconciled tree
tr_ultra <- chronos(tr_full, lambda = 1)

# Make ultrametric (chronos)
tr_ultra <- chronos(tr_full, lambda = 1)

# ---- helper: split a cell like "22,Z" into character vector c("22","Z")
split_chr <- function(x) {
  if (is.null(x) || length(x) == 0) return(character())
  x <- as.character(x)
  if (is.na(x) || x == "") return(character())
  parts <- unlist(strsplit(x, ",", fixed = TRUE))
  parts <- str_trim(parts)
  parts[parts != ""]
}

# ---- Build long format: Species x RefChr
out_long <- out %>%
  mutate(
    # make each of the three columns a list column of tokens
    both_list  = map(RefChr_both,      split_chr),
    c1_list    = map(RefChr_chr1_only, split_chr),
    c2_list    = map(RefChr_chr2_only, split_chr),
    ref_all     = map2(map2(both_list, c1_list, union), c2_list, union)
  ) %>%
  select(Species, ref_all) %>%
  unnest(ref_all) %>%
  filter(!is.na(ref_all), ref_all != "") %>%
  distinct(Species, ref_all) %>%
  rename(RefChr = ref_all) %>%
  distinct(Species, RefChr)

# order reference chromosomes: numeric first, then others (Z/W/etc)
chr_levels <- out_long %>%
  distinct(RefChr) %>%
  mutate(
    is_num = suppressWarnings(!is.na(as.integer(RefChr))),
    numval = suppressWarnings(as.integer(RefChr))
  ) %>%
  arrange(desc(is_num), numval, RefChr) %>%
  pull(RefChr)

# full grid: all tips with data x all chrs
grid_df <- tidyr::expand_grid(
  Species = unique(out_long$Species),
  RefChr   = chr_levels
) %>%
  left_join(out_long %>% mutate(present = 1L), by = c("Species", "RefChr")) %>%
  mutate(present = if_else(is.na(present), 0L, present))

# ---- Plot
pdf(out_pdf, width = 3, height = 10)

# more right margin for dot grid + top margin for labels
par(mar = c(2, 2, 6, 12))

tmp <- plot(tr_ultra, plot = FALSE)
nchr <- length(chr_levels)

dx <- max(tmp$x.lim) * 0.02
grid_width <- dx * (nchr + 6)

plot(tr_ultra,
     cex = 0.1,
     edge.width = 0.25,
     no.margin = TRUE,
     x.lim = c(0, max(tmp$x.lim) + grid_width))

lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_y <- lp$yy[1:Ntip(tr_ultra)]
tips  <- tr_ultra$tip.label

x0 <- max(tmp$x.lim) + dx * 3

tips_with_data <- intersect(tips, unique(grid_df$Species))
idx_data <- match(tips_with_data, tips)
y_data <- tip_y[idx_data]

for (j in seq_along(chr_levels)) {
  chrj <- chr_levels[j]
  xj <- x0 + (j - 1) * dx

  present_species <- grid_df %>%
    filter(RefChr == chrj, present == 1L) %>%
    pull(Species)

  is_present_data <- tips_with_data %in% present_species

  points(rep(xj, length(y_data)), y_data,
         pch = 21, cex = 0.05, col = "lightgray", bg = "white")

  points(rep(xj, sum(is_present_data)), y_data[is_present_data],
         pch = 21, cex = 0.05, bg = "black")
}

# Column labels at the top (keep your fixed y, but fall back if it’s too small)
usr <- par("usr")
y_lab <- 585

text(
  x = x0 + (seq_along(chr_levels) - 1) * dx,
  y = y_lab,
  labels = chr_levels,
  srt = 90,
  adj = c(0, 0.5),
  xpd = NA,
  cex = 0.2
)

dev.off()

cat("FOCAL:", FOCAL, "\n")
cat("Read:", infile, "\n")
cat("Wrote:", out_pdf, "\n")
cat("Species with data:", length(unique(out_long$Species)), "\n")
cat("Chicken chromosomes plotted:", length(chr_levels), "\n")