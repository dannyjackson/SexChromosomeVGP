# Plot it on a phylogeny

## Blanks for empty?

library(data.table)
library(ggplot2)
library(patchwork)
library(ape)
library(dplyr)
library(tidyr)
library(purrr)


out <- read.csv("syntentic_chromosomes.sex_shared.Gallus_gallus_REF.csv")

tree_file <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/roadies_v1.1.16b.numbers.scientific.nwk"

tr_full <- ape::read.tree(tree_file)

tip_key <- tr_full$tip.label

tr_ultra <- chronos(tr_full, lambda = 1)


# out: your tibble with Species and list-cols
# tr_ultra: your ultrametric tree

# union the three list-cols

out_long <- out %>%
  mutate(
    gg_all = map2(
      map2(GgChr_both, GgChr_chr1_only, ~ union(.x, .y)),
      GgChr_chr2_only,
      ~ union(.x, .y)
    )
  ) %>%
  select(Species, gg_all) %>%
  unnest(gg_all) %>%
  filter(!is.na(gg_all), gg_all != "") %>%
  distinct(Species, gg_all) %>%
  rename(GgChr = gg_all) %>%
  # --- NEW: split combos like "22,Z" into separate rows "22" and "Z"
  separate_rows(GgChr, sep = ",") %>%
  mutate(GgChr = str_trim(GgChr)) %>%
  filter(GgChr != "") %>%
  distinct(Species, GgChr)

# order chicken chromosomes: numeric first, then others (Z/W/etc)
chr_levels <- out_long %>%
  distinct(GgChr) %>%
  mutate(
    is_num = suppressWarnings(!is.na(as.integer(GgChr))),
    numval = suppressWarnings(as.integer(GgChr))
  ) %>%
  arrange(desc(is_num), numval, GgChr) %>%
  pull(GgChr)

# full grid for all tips x all chrs
grid_df <- tidyr::expand_grid(
  Species = out_long$Species,
  GgChr   = chr_levels
) %>%
  left_join(out_long %>% mutate(present = 1L), by = c("Species", "GgChr")) %>%
  mutate(present = if_else(is.na(present), 0L, present))





# Plot tree

pdf("tree.pdf", width = 3, height = 10)

# give more right margin for the dot grid + top margin for labels
par(mar = c(2, 2, 6, 12))

# 1) Compute how much horizontal space you need for the dot grid
#    (dx * nchr) plus some padding.
tmp <- plot(tr_ultra, plot = FALSE)
nchr <- length(chr_levels)

# spacing between chromosome columns
dx <- max(tmp$x.lim) * 0.02

# extra horizontal room to the right of the tree (in tree x-units)
grid_width <- dx * (nchr + 6)

# 2) Plot tree but "squish it left" by expanding xlim to the right
plot(tr_ultra,
     cex = 0.1,
     edge.width = 0.25,
     no.margin = TRUE,
     x.lim = c(0, max(tmp$x.lim) + grid_width))

lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_y <- lp$yy[1:Ntip(tr_ultra)]
tips  <- tr_ultra$tip.label

# start x position for the grid (a bit to the right of the tree)
x0 <- max(tmp$x.lim) + dx * 3

# full tips/y from the plotted tree
lp   <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tips <- tr_ultra$tip.label
tip_y <- lp$yy[1:Ntip(tr_ultra)]

# subset: tips that appear in grid_df
tips_with_data <- intersect(tips, unique(grid_df$Species))
idx_data <- match(tips_with_data, tips)
y_data <- tip_y[idx_data]

# within each chromosome column:
for (j in seq_along(chr_levels)) {
  chrj <- chr_levels[j]
  xj <- x0 + (j - 1) * dx

  present_species <- grid_df %>%
    dplyr::filter(GgChr == chrj, present == 1L) %>%
    dplyr::pull(Species)

  is_present_data <- tips_with_data %in% present_species

  # light-gray circles ONLY for species present in grid_df
  points(rep(xj, length(y_data)), y_data,
         pch = 21, cex = 0.05, col = "lightgray", bg = "white")

  # black filled circles for present==1 among those species
  points(rep(xj, sum(is_present_data)), y_data[is_present_data],
         pch = 21, cex = 0.05, bg = "black")
}

# 2) Column labels at the top
usr <- par("usr")

# Put labels slightly ABOVE the top of the plotting region
y_lab <- 585
text(
  x = x0 + (seq_along(chr_levels) - 1) * dx,
  y = y_lab,
  labels = chr_levels,
  srt = 90,                 # rotate
  adj = c(0, 0.5),          # anchor at bottom of rotated text
  xpd = NA,                 # allow drawing into margin area
  cex = 0.2                 # tweak as needed
)

dev.off()