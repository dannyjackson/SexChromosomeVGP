#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript plot_gg_regions.R gg_regions.csv gg_chr_len.csv palette.txt output.png")
}

gg_regions_file <- args[1]
gg_chr_len_file <- args[2]
palette_file <- args[3]
output_file <- args[4]

# Read chromosome lengths
chr_len <- read_csv(
  gg_chr_len_file,
  col_names = c("chr", "length"),
  show_col_types = FALSE
) %>%
  mutate(
    chr = as.character(chr),
    length = as.numeric(length)
  )

# Preserve chromosome order from file
chr_order <- chr_len$chr

# Build cumulative offsets so chromosomes are laid out end-to-end
gap <- 0

chr_len <- chr_len %>%
  mutate(
    offset = cumsum(lag(length + gap, default = 0)),
    center = offset + length / 2
  )

# Make df of chromosome boundaries
chr_bounds <- chr_len %>%
  transmute(
    chr,
    start = offset,
    end = offset + length,
    center = offset + length / 2
  )

boundaries <- sort(unique(c(chr_len$offset, chr_len$offset + chr_len$length)))

# Read regions with clade column
regions <- read_csv(gg_regions_file, show_col_types = FALSE) %>%
  mutate(
    Gg_regions = ifelse(is.na(Gg_regions), "", Gg_regions),
    Species = as.character(Species),
    Clade = as.character(Clade)
  )

# Read palette
palette_df <- read_csv(
  palette_file,
  col_names = c("Clade", "Color"),
  show_col_types = FALSE
) %>%
  mutate(
    Clade = str_trim(as.character(Clade)),
    Color = str_trim(as.character(Color))
  )

palette_named <- setNames(palette_df$Color, palette_df$Clade)

# Expand semicolon-separated blocks
regions_long <- regions %>%
  separate_rows(Gg_regions, sep = ";") %>%
  mutate(Gg_regions = str_trim(Gg_regions)) %>%
  filter(Gg_regions != "")

# Parse chr:start-end
regions_long <- regions_long %>%
  extract(
    Gg_regions,
    into = c("chr", "start", "end"),
    regex = "^([^:]+):(\\d+)-(\\d+)$",
    remove = FALSE
  ) %>%
  mutate(
    chr = as.character(chr),
    start = as.numeric(start),
    end = as.numeric(end)
  ) %>%
  filter(!is.na(chr), !is.na(start), !is.na(end)) %>%
  mutate(
    xmin_local = pmin(start, end),
    xmax_local = pmax(start, end)
  ) %>%
  left_join(chr_len, by = "chr") %>%
  filter(!is.na(offset)) %>%
  mutate(
    xmin = offset + xmin_local,
    xmax = offset + xmax_local
  )

# Preserve species order from input file
species_order <- unique(regions$Species)

regions_long <- regions_long %>%
  mutate(Species = factor(Species, levels = rev(species_order)))

# Also set factor order for the full regions table
regions <- regions %>%
  mutate(Species = factor(Species, levels = rev(species_order)))

p <- ggplot() +
  geom_vline(
  xintercept = boundaries,
  linewidth = 0.1,
  color = "lightgray"
) +
  geom_segment(
    data = regions_long,
    aes(x = xmin, xend = xmax, y = Species, yend = Species, color = Clade),
    linewidth = 4,
    lineend = "butt"
  ) +
  scale_color_manual(
    values = palette_named,
    na.value = "grey50"
  ) +
  scale_x_continuous(
    breaks = chr_bounds$center,
    labels = chr_bounds$chr,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Gallus gallus chromosomes",
    y = "Unique SC systems",
    color = "Clade",
    title = "Chicken syntenic regions by SC system"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

height_inches <- max(4, 0.3 * length(species_order))

ggsave(
  output_file,
  plot = p,
  width = 16,
  height = height_inches,
  dpi = 300
)