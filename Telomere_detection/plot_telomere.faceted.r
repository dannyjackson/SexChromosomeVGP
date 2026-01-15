library(tidyverse)
library(ragg)
base_dir <- "/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/Telomere_detection/sexchrs/tidk_results"

cat("Starting telomere plot generation\n")
cat("Base directory:", base_dir, "\n\n")

all_data <- list()

species_dirs <- list.dirs(base_dir, recursive = FALSE)
cat("Found", length(species_dirs), "species directories\n\n")

for (i in seq_along(species_dirs)) {

  species_dir <- species_dirs[i]
  species <- basename(species_dir)

  cat(
    sprintf(
      "[%d/%d] Processing species: %s\n",
      i, length(species_dirs), species
    )
  )

  tsv_files <- list.files(
    species_dir,
    pattern = "\\.tsv$",
    full.names = TRUE
  )

  cat("  Found", length(tsv_files), "TSV files\n")

  for (j in seq_along(tsv_files)) {

    tsv <- tsv_files[j]
    fname <- basename(tsv)

    cat(
      sprintf(
        "    (%d/%d) Reading %s\n",
        j, length(tsv_files), fname
      )
    )

    sexchr <- str_split(fname, "\\.")[[1]][2]

    df <- read_tsv(tsv, show_col_types = FALSE) %>%
      mutate(
        species = species,
        sexchr = sexchr,
        chr_label = paste(sexchr, id, sep = ":")
      )

    all_data[[length(all_data) + 1]] <- df
  }

  cat("  Finished species:", species, "\n\n")
}

cat("Combining data frames...\n")
df_all <- bind_rows(all_data)

cat("Total rows:", nrow(df_all), "\n")
cat("Total species:", length(unique(df_all$species)), "\n")
cat("Total chromosomes:", length(unique(df_all$chr_label)), "\n\n")

cat("Reshaping data for plotting...\n")
df_long <- df_all %>%
  pivot_longer(
    cols = c(forward_repeat_number, reverse_repeat_number),
    names_to = "direction",
    values_to = "count"
  )

cat("Generating plot...\n")

p <- ggplot(df_long, aes(x = window, y = count, color = direction)) +
  geom_point(size = 0.6) +
  scale_color_manual(
    values = c(
      forward_repeat_number = "red",
      reverse_repeat_number = "blue"
    ),
    labels = c("Forward", "Reverse")
  ) +
  facet_grid(
    species ~ chr_label,
    scales = "free_x",
    space = "free_x"
  ) +
  labs(
    x = "Genomic window position (bp)",
    y = "Telomeric repeat count",
    color = "Repeat orientation"
  ) +
  theme_bw() +
  theme(
    strip.text.x = element_text(angle = 90),
    panel.spacing = unit(0.2, "lines")
  )

ggsave(
  filename = "telomere_counts.png",
  plot = p,
  device = ragg::agg_png,
  width = 18, height = 10, units = "in",
  dpi = 150   # start lower
)

out_file <- "all_species_all_chr_telomere_counts.png"

ggsave(
  out_file,
  plot = p,
  width = 16,
  height = 10,
  dpi = 300
)

cat("\nPlot saved to:", out_file, "\n")
cat("Done.\n")
