#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(readr)
library(forcats)
library(scales)

outdir <- "PAR_gap_results"

summary <- read_tsv(
  file.path(outdir, "PAR.gaps.summary.tsv"),
  show_col_types = FALSE
)

gap_file <- file.path(outdir, "PAR.gaps.bed")

gap_cols <- c(
  "species",
  "chrom",
  "gap_start",
  "gap_end",
  "gap_len",
  "PAR_start",
  "PAR_end",
  "rel_gap_start",
  "rel_gap_end"
)

if (file.info(gap_file)$size == 0) {
  gaps <- tibble(
    species = character(),
    chrom = character(),
    gap_start = numeric(),
    gap_end = numeric(),
    gap_len = numeric(),
    PAR_start = numeric(),
    PAR_end = numeric(),
    rel_gap_start = numeric(),
    rel_gap_end = numeric()
  )
} else {
  gaps <- read_tsv(
    gap_file,
    col_names = gap_cols,
    show_col_types = FALSE
  )
}

summary <- summary %>%
  arrange(PAR_len) %>%
  mutate(
    species = factor(species, levels = species),
    PAR_start_rel = 0,
    PAR_end_rel = PAR_len
  )

gaps <- gaps %>%
  mutate(
    species = factor(species, levels = levels(summary$species))
  )

p1 <- ggplot() +
  geom_segment(
    data = summary,
    aes(
      x = PAR_start_rel,
      xend = PAR_end_rel,
      y = species,
      yend = species
    ),
    linewidth = 0.7,
    color = "grey70"
  ) +
  geom_segment(
    data = gaps,
    aes(
      x = rel_gap_start,
      xend = rel_gap_end,
      y = species,
      yend = species
    ),
    linewidth = 2.2,
    color = "black"
  ) +
  scale_x_continuous(labels = comma) +
  labs(
    title = "Assembly gaps within PAR regions",
    x = "Position within PAR region, bp",
    y = "Species"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7)
  )

ggsave(
  filename = file.path(outdir, "PAR_gaps_plot.png"),
  plot = p1,
  width = 12,
  height = max(6, 0.28 * nrow(summary)),
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "PAR_gaps_plot.pdf"),
  plot = p1,
  width = 12,
  height = max(6, 0.28 * nrow(summary))
)

summary2 <- summary %>%
  mutate(
    gap_fraction = total_gap_bp / PAR_len,
    species = fct_reorder(species, gap_fraction)
  )

p2 <- ggplot(summary2, aes(x = gap_fraction, y = species)) +
  geom_col() +
  scale_x_continuous(labels = percent_format(accuracy = 0.01)) +
  labs(
    title = "PAR gap fraction by species",
    x = "Fraction of PAR region that is gap sequence",
    y = "Species"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7)
  )

ggsave(
  filename = file.path(outdir, "PAR_gap_fraction_by_species.png"),
  plot = p2,
  width = 10,
  height = max(6, 0.28 * nrow(summary)),
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "PAR_gap_fraction_by_species.pdf"),
  plot = p2,
  width = 10,
  height = max(6, 0.28 * nrow(summary))
)

message("Wrote:")
message("  ", file.path(outdir, "PAR_gaps_plot.png"))
message("  ", file.path(outdir, "PAR_gaps_plot.pdf"))
message("  ", file.path(outdir, "PAR_gap_fraction_by_species.png"))
message("  ", file.path(outdir, "PAR_gap_fraction_by_species.pdf"))