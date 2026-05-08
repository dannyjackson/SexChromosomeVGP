#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
  stop("Usage: Rscript plot_gg_regions_by_clade.Neo_and_Ancestral.R gg_regions.Neo.csv gg_regions.Ancestral.csv gg_chr_len.csv palette.Neo.txt palette.Ancestral.txt output.png")
}

gg_regions_neo_file <- args[1]
gg_regions_ancestral_file <- args[2]
gg_chr_len_file <- args[3]
palette_neo_file <- args[4]
palette_ancestral_file <- args[5]
output_file <- args[6]

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

# Read palettes
palette_neo <- read_csv(
  palette_neo_file,
  col_names = c("Clade", "Color"),
  show_col_types = FALSE
) %>%
  mutate(
    Clade = str_trim(as.character(Clade)),
    Color = str_trim(as.character(Color)),
    RowID = paste0("Neo__", Clade)
  )

palette_ancestral <- read_csv(
  palette_ancestral_file,
  col_names = c("Clade", "Color"),
  show_col_types = FALSE
) %>%
  mutate(
    Clade = str_trim(as.character(Clade)),
    Color = str_trim(as.character(Color)),
    RowID = paste0("Ancestral__", Clade)
  )

palette_named <- c(
  setNames(palette_ancestral$Color, palette_ancestral$RowID),
  setNames(palette_neo$Color, palette_neo$RowID)
)

# Read region files
read_regions_file <- function(file, dataset_name) {
  read_csv(file, show_col_types = FALSE) %>%
    mutate(
      Species = as.character(Species),
      Clade = as.character(Clade),
      Gg_regions = ifelse(is.na(Gg_regions), "", Gg_regions),
      Dataset = dataset_name
    )
}

regions_neo <- read_regions_file(gg_regions_neo_file, "Neo")
regions_ancestral <- read_regions_file(gg_regions_ancestral_file, "Ancestral")

# Keep clade order based on first appearance
ancestral_clade_order <- unique(regions_ancestral$Clade)
neo_clade_order <- unique(regions_neo$Clade)

# Build row order:
# for each ancestral clade, put Ancestral row first and Neo row directly below if present
# then append neo-only clades at the bottom
row_order <- c()

for (cl in ancestral_clade_order) {
  row_order <- c(row_order, paste0("Ancestral__", cl))
  if (cl %in% neo_clade_order) {
    row_order <- c(row_order, paste0("Neo__", cl))
  }
}

neo_only_clades <- setdiff(neo_clade_order, ancestral_clade_order)
for (cl in neo_only_clades) {
  row_order <- c(row_order, paste0("Neo__", cl))
}

# Combine datasets
regions <- bind_rows(regions_ancestral, regions_neo) %>%
  mutate(
    RowID = paste0(Dataset, "__", Clade),
    RowLabel = paste0(Clade, ifelse(Dataset == "Neo", " (Neo)", " (Ancestral)"))
  )

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

# Merge overlapping or touching intervals within each dataset/clade/chromosome
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
  group_by(Dataset, Clade, RowID, RowLabel, chr) %>%
  group_modify(~ merge_intervals(.x)) %>%
  ungroup() %>%
  left_join(chr_len, by = "chr") %>%
  filter(!is.na(offset)) %>%
  mutate(
    xmin = offset + start,
    xmax = offset + end
  )

# Assign rows so they touch with no vertical gaps
row_levels <- rev(row_order)

row_positions <- tibble(
  RowID = row_levels,
  y = seq_along(row_levels)
) %>%
  mutate(
    RowLabel = case_when(
      str_starts(RowID, "Ancestral__") ~ paste0(str_remove(RowID, "^Ancestral__"), " (Ancestral)"),
      str_starts(RowID, "Neo__") ~ paste0(str_remove(RowID, "^Neo__"), " (Neo)")
    )
  )

clade_regions <- clade_regions %>%
  left_join(row_positions, by = c("RowID", "RowLabel")) %>%
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
      fill = RowID
    ),
    color = NA
  ) +
  scale_fill_manual(
    values = palette_named,
    breaks = row_positions$RowID,
    labels = row_positions$RowLabel,
    na.value = "grey50"
  ) +
  scale_x_continuous(
    breaks = chr_len$center,
    labels = chr_len$chr,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = row_positions$y,
    labels = row_positions$RowLabel,
    expand = c(0, 0)
  ) +
  coord_cartesian(
    ylim = c(0.5, length(row_levels) + 0.5),
    expand = FALSE
  ) +
  labs(
    x = "Gallus gallus chromosomes",
    y = "Clade",
    fill = "Dataset / Clade",
    title = "Chicken syntenic regions by clade"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

height_inches <- max(4, 0.22 * length(row_levels))

ggsave(
  output_file,
  plot = p,
  width = 16,
  height = height_inches,
  dpi = 300
)