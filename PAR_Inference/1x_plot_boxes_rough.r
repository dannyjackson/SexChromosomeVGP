library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(forcats)

# /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/align_PAR_99thr
# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------
info <- read_tsv("bed_sexchr_info_with_matches.tsv",
                 show_col_types = FALSE)

bed_dir <- "sorted_beds"

# ------------------------------------------------------------
# Read BED regions into a long table
# ------------------------------------------------------------
read_bed <- function(file) {
    path <- file.path(bed_dir, file)
    if(!file.exists(path)) {
        return(NULL)
    }
    bed <- read_tsv(path, col_names = FALSE,
                    show_col_types = FALSE)
    if(nrow(bed) == 0) {
        return(NULL)
    }
    colnames(bed) <- c("chrom", "start", "end")
    bed$Filename <- file
    return(bed)
}

bed_regions <- bind_rows(lapply(info$Filename, read_bed))

# merge with metadata
plotdata <- info %>%
    left_join(bed_regions, by="Filename")

# Remove species whose bed was empty or missing
# They will still be plotted as empty chromosomes
plotdata$start[is.na(plotdata$start)] <- 0
plotdata$end[is.na(plotdata$end)] <- 0

# ------------------------------------------------------------
# SCALE WITHIN EACH SPECIES: every chromosome bar ranges 0–1
# ------------------------------------------------------------
plotdata <- plotdata %>%
    mutate(
        scaled_start = start / Seq_length,
        scaled_end   = end   / Seq_length
    )

# ------------------------------------------------------------
# Order species for plotting
# ------------------------------------------------------------
plotdata$Species <- factor(plotdata$Species,
                           levels = rev(unique(plotdata$Species)))



# ------------------------------------------------------------
# Build a "line id" that becomes Species (chrom) when needed
# ------------------------------------------------------------
species_levels <- rev(unique(plotdata$Species))  # your desired species order

plotdata2 <- plotdata %>%
  mutate(Species = factor(Species, levels = species_levels)) %>%
  group_by(Species) %>%
  mutate(
    n_chrom = n_distinct(chrom),
    chrom_rank = dense_rank(as.character(chrom)),
    Species_line = if_else(
      n_chrom > 1,
      paste0(as.character(Species), " (", chrom, ")"),
      as.character(Species)
    )
  ) %>%
  ungroup() %>%
  # order rows: by Species order, then chrom order within species
  arrange(Species, chrom_rank) %>%
  mutate(Species_line = factor(Species_line, levels = rev(unique(Species_line))))

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p <- ggplot() +

  # Chromosome background bars
  geom_rect(
    data = distinct(plotdata2, Species_line, Seq_length, Filename, chrom),
    aes(
      xmin = 0, xmax = 1,
      ymin = as.numeric(Species_line) - 0.35,
      ymax = as.numeric(Species_line) + 0.35
    ),
    fill = "white", color = "black"
  ) +

  # Red BED intervals
  geom_rect(
    data = plotdata2,
    aes(
      xmin = scaled_start,
      xmax = scaled_end,
      ymin = as.numeric(Species_line) - 0.35,
      ymax = as.numeric(Species_line) + 0.35
    ),
    fill = "red", alpha = 0.7
  ) +

  # Labels
  scale_y_continuous(
    breaks = seq_along(levels(plotdata2$Species_line)),
    labels = levels(plotdata2$Species_line),
    expand = expansion(mult = 0)
  ) +

  labs(
    title = "Chromosome Regions Across Species",
    x = "Scaled Chromosome Position (0–1)",
    y = "Species"
  ) +

  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10)
  )


# ------------------------------------------------------------
# Output PNG
# ------------------------------------------------------------
png("chromosome_regions.png", width = 2000, height = 9000, res = 200)
print(p)
dev.off()


min(plotdata %>%
  filter(grepl("Eudromia", Species)) %>%
  select(scaled_start))

min(plotdata %>%
  filter(grepl("Eudromia", Species)) %>%
  select(start))

# ℹ Use `print(n = ...)` to see more rows
> min(plotdata %>%
+   filter(grepl("Eudromia", Species)) %>%
+   select(start))
[1] 10842918
> 
> min(plotdata %>%
+   filter(grepl("Eudromia", Species)) %>%
+   select(scaled_start))
[1] 0.1475924
> q()


plotdata %>%
      dplyr::filter(grepl("Hoplias", Species))%>%
      select(GenBank_accession)%>%
      unique()
