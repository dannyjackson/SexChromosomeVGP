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
  stop("Usage: Rscript plot_gg_regions_by_clade.R gg_regions.csv gg_chr_len.csv palette.txt output.png")
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

gap <- 5000000

chr_len <- chr_len %>%
  mutate(
    offset = cumsum(lag(length + gap, default = 0)),
    center = offset + length / 2
  )

# Chromosome boundaries for thin vertical divider lines
boundaries <- sort(unique(c(chr_len$offset, chr_len$offset + chr_len$length)))

# Read regions
regions <- read_csv(gg_regions_file, show_col_types = FALSE) %>%
  mutate(
    Species = as.character(Species),
    Clade = as.character(Clade),
    Gg_regions = ifelse(is.na(Gg_regions), "", Gg_regions)
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

# Keep clade order based on first appearance in input
clade_order <- unique(regions$Clade)

# Expand semicolon-separated regions
regions_long <- regions %>%
  separate_rows(Gg_regions, sep = ";") %>%
  mutate(Gg_regions = str_trim(Gg_regions)) %>%
  filter(Gg_regions != "") %>%
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
  filter(!is.na(chr), chr != "", !is.na(start), !is.na(end)) %>%
  mutate(
    start = pmin(start, end),
    end = pmax(start, end)
  )

# Merge overlapping or touching intervals within each clade and chromosome
merge_intervals <- function(df) {
  df <- df %>%
    filter(!is.na(start), !is.na(end)) %>%
    arrange(start, end)

  if (nrow(df) == 0) {
    return(tibble(start = numeric(0), end = numeric(0)))
  }

  out_start <- c()
  out_end <- c()

  cur_start <- df$start[1]
  cur_end <- df$end[1]

  if (nrow(df) >= 2) {
    for (i in 2:nrow(df)) {
      s <- df$start[i]
      e <- df$end[i]

      if (s <= cur_end + 1) {
        cur_end <- max(cur_end, e)
      } else {
        out_start <- c(out_start, cur_start)
        out_end <- c(out_end, cur_end)
        cur_start <- s
        cur_end <- e
      }
    }
  }

  out_start <- c(out_start, cur_start)
  out_end <- c(out_end, cur_end)

  tibble(
    start = out_start,
    end = out_end
  )
}

clade_regions <- regions_long %>%
  group_by(Clade, chr) %>%
  group_modify(~ merge_intervals(.x)) %>%
  ungroup() %>%
  left_join(chr_len, by = "chr") %>%
  filter(!is.na(offset)) %>%
  mutate(
    xmin = offset + start,
    xmax = offset + end
  )

# Assign clades to touching rows with no gaps
clade_levels <- rev(clade_order)

clade_positions <- tibble(
  Clade = clade_levels,
  y = seq_along(clade_levels)
)

clade_regions <- clade_regions %>%
  left_join(clade_positions, by = "Clade") %>%
  mutate(
    ymin = y - 0.5,
    ymax = y + 0.5
  )

p <- ggplot() +
  geom_vline(
    xintercept = boundaries,
    linewidth = 0.1,
    color = "lightgray"
  ) +
  geom_rect(
    data = clade_regions,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = Clade
    ),
    color = NA
  ) +
  scale_fill_manual(
    values = palette_named,
    na.value = "grey50"
  ) +
  scale_x_continuous(
    breaks = chr_len$center,
    labels = chr_len$chr,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = clade_positions$y,
    labels = clade_positions$Clade,
    expand = c(0, 0)
  ) +
  coord_cartesian(
    ylim = c(0.5, length(clade_levels) + 0.5),
    expand = FALSE
  ) +
  labs(
    x = "Gallus gallus chromosomes",
    y = "Clade",
    fill = "Clade",
    title = "Chicken syntenic regions by clade"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

height_inches <- max(4, 0.22 * length(clade_order))

ggsave(
  output_file,
  plot = p,
  width = 16,
  height = height_inches,
  dpi = 300
)