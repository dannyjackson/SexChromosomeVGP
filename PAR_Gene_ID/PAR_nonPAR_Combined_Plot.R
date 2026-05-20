#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
  library(ggtree)
  library(patchwork)
  library(scales)
})

# ============================================================
# Command-line arguments
# Usage:
#   Rscript PAR_Combined_Plot.R birds
#   Rscript PAR_Combined_Plot.R mammals
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1 || !args[1] %in% c("birds", "mammals")) {
  stop(
    "Usage:\n",
    "  Rscript PAR_Combined_Plot.R birds\n",
    "  Rscript PAR_Combined_Plot.R mammals\n",
    call. = FALSE
  )
}

taxon <- args[1]

# ============================================================
# Inputs
# ============================================================

tree_file <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/roadies_v1.1.16b.numbers.scientific.nwk"
par_size_file <- "species_par.csv"

if (taxon == "birds") {
  gene_file <- "gene_locations_by_species.with_chr_label.all.In_PAR.Z_only.gapless_species.csv"
  telomere_file <- "gapless_species.telomeres.txt"
  output_file <- "combined_phylogeny_upset_gene_order_PAR_nonPAR_size.birds.pdf"
  plot_title <- "Phylogeny, PAR gene intersections, PAR gene order, and PAR size - birds"
}

if (taxon == "mammals") {
  gene_file <- "gene_locations_by_species.with_chr_label.all.In_PAR.X_only.gapless_species.csv"
  telomere_file <- "gapless_species.telomeres.txt"
  output_file <- "combined_phylogeny_upset_gene_order_PAR_nonPAR_size.mammals.pdf"
  plot_title <- "Phylogeny, PAR gene intersections, PAR gene order, and PAR size - mammals"
}

# ============================================================
# Read data
# ============================================================

df <- read_csv(gene_file, show_col_types = FALSE)
tree <- read.tree(tree_file)

# ============================================================
# Keep genes that are in the PAR in at least one species
# ============================================================

genes_in_PAR_any_species <- df %>%
  filter(In_PAR == "Y") %>%
  distinct(Gene) %>%
  pull(Gene)

df_par_relevant <- df %>%
  filter(Gene %in% genes_in_PAR_any_species)

# ============================================================
# Build PAR binary matrix for UpSet-style panel
# rows = genes, columns = species, TRUE/FALSE = gene in PAR
# ============================================================

par_binary <- df %>%
  filter(Gene %in% genes_in_PAR_any_species) %>%
  mutate(in_par_logical = In_PAR == "Y") %>%
  group_by(Gene, Species) %>%
  summarise(
    in_PAR = any(in_par_logical),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Species,
    values_from = in_PAR,
    values_fill = FALSE
  )

species_cols <- setdiff(names(par_binary), "Gene")

# ============================================================
# Tree pruning and species order
# ============================================================

tree_filtered <- keep.tip(
  tree,
  intersect(tree$tip.label, species_cols)
)

p_tree_tmp <- ggtree(tree_filtered)

tree_plot_order <- p_tree_tmp$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)

species_order <- tree_plot_order

# Restrict matrix and data to species in the tree
par_binary_tree <- par_binary %>%
  select(Gene, all_of(species_order))

df_par_relevant <- df_par_relevant %>%
  filter(Species %in% species_order)

# ============================================================
# 1. Phylogeny panel
# ============================================================

p_tree <- ggtree(tree_filtered) +
  geom_tiplab(size = 3, align = FALSE) +
  xlim_tree(0.25) +
  theme_tree2() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

# ============================================================
# 2. UpSet-style intersection panel
# ============================================================

gene_intersections <- par_binary_tree %>%
  rowwise() %>%
  mutate(
    intersection_id = paste(
      species_order[c_across(all_of(species_order))],
      collapse = "|"
    )
  ) %>%
  ungroup() %>%
  filter(intersection_id != "")

intersection_counts <- gene_intersections %>%
  count(intersection_id, name = "n_genes") %>%
  arrange(desc(n_genes)) %>%
  mutate(
    intersection_index = row_number(),
    intersection_index = factor(intersection_index, levels = intersection_index)
  )

intersection_matrix <- intersection_counts %>%
  separate_longer_delim(intersection_id, delim = "|") %>%
  rename(Species = intersection_id) %>%
  mutate(
    Species = factor(Species, levels = species_order),
    intersection_index = factor(
      intersection_index,
      levels = levels(intersection_counts$intersection_index)
    )
  )

p_bar <- ggplot(intersection_counts, aes(x = intersection_index, y = n_genes)) +
  geom_col() +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = NULL,
    y = NULL,
    title = "Genes in intersection"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 10, face = "bold"),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

p_matrix <- ggplot(intersection_matrix, aes(x = intersection_index, y = Species)) +
  geom_line(aes(group = intersection_index), linewidth = 0.3) +
  geom_point(size = 2.5) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(limits = species_order) +
  labs(
    x = "PAR species intersection",
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 5.5, 5.5, 5.5)
  )

# ============================================================
# 3. PAR gene order panel
# ============================================================

par_genes <- df %>%
  filter(In_PAR %in% c("Y", "Edge")) %>%
  filter(Species %in% species_order) %>%
  mutate(
    midpoint = (Start_pos + Stop_pos) / 2
  ) %>%
  group_by(Species, Chrom) %>%
  mutate(
    PAR_at_end = min(Start_pos, na.rm = TRUE) >= 10000000,
    adjusted_pos = if_else(PAR_at_end, -midpoint, midpoint)
  ) %>%
  arrange(Species, Chrom, adjusted_pos) %>%
  mutate(
    PAR_order = row_number()
  ) %>%
  ungroup()

gene_freq <- par_genes %>%
  distinct(Species, Gene) %>%
  count(Gene, name = "species_frequency")

par_genes <- par_genes %>%
  left_join(gene_freq, by = "Gene") %>%
  mutate(
    Species = factor(Species, levels = species_order)
  )

# ============================================================
# 3. PAR gene order panel
# Includes telomeres, ortholog segments, non-PAR genes, PAR genes, and Edge genes
# ============================================================

# Keep genes that are PAR in at least one species, but include their
# non-PAR occurrences in other species as light gray.
display_genes <- df_par_relevant %>%
  filter(Species %in% species_order) %>%
  mutate(
    midpoint = (Start_pos + Stop_pos) / 2,
    PAR_status = case_when(
      In_PAR == "Y" ~ "PAR",
      In_PAR == "Edge" ~ "Edge",
      TRUE ~ "Non-PAR"
    )
  ) %>%
  group_by(Species, Chrom) %>%
  mutate(
    PAR_at_end = min(Start_pos[PAR_status %in% c("PAR", "Edge")], na.rm = TRUE) >= 10000000,
    adjusted_pos = if_else(PAR_at_end, -midpoint, midpoint)
  ) %>%
  arrange(Species, Chrom, adjusted_pos) %>%
  mutate(
    PAR_order = row_number()
  ) %>%
  ungroup()

gene_freq <- display_genes %>%
  filter(PAR_status == "PAR") %>%
  distinct(Species, Gene) %>%
  count(Gene, name = "species_frequency")

display_genes <- display_genes %>%
  left_join(gene_freq, by = "Gene") %>%
  mutate(
    Species = factor(Species, levels = species_order),
    species_index = as.integer(Species)
  ) %>%
  group_by(Species, Gene) %>%
  mutate(
    Gene_label = if_else(
      n() > 1,
      paste0(Gene, "_hit", Hit_rank),
      Gene
    )
  ) %>%
  ungroup()

par_genes_y <- display_genes %>%
  filter(PAR_status == "PAR")

par_genes_edge <- display_genes %>%
  filter(PAR_status == "Edge")

non_par_genes <- display_genes %>%
  filter(PAR_status == "Non-PAR")

# Purple boundary line between PAR/Edge and non-PAR genes
par_nonpar_boundaries <- display_genes %>%
  arrange(Species, Chrom, PAR_order) %>%
  group_by(Species, Chrom) %>%
  mutate(
    next_PAR_order = lead(PAR_order),
    next_status = lead(PAR_status),
    boundary_x = (PAR_order + next_PAR_order) / 2,
    is_boundary = !is.na(next_status) &
      ((PAR_status %in% c("PAR", "Edge") & next_status == "Non-PAR") |
         (PAR_status == "Non-PAR" & next_status %in% c("PAR", "Edge")))
  ) %>%
  filter(is_boundary) %>%
  ungroup() %>%
  mutate(
    y = species_index - 0.35,
    yend = species_index + 0.35
  )

ortholog_segments <- display_genes %>%
  group_by(Species, Gene) %>%
  slice_min(PAR_order, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Gene, species_index) %>%
  group_by(Gene) %>%
  mutate(
    xend = lead(PAR_order),
    yend = lead(species_index),
    next_species_index = lead(species_index)
  ) %>%
  filter(!is.na(xend)) %>%
  filter(next_species_index == species_index + 1) %>%
  ungroup()

telomere_points <- read_csv(
  telomere_file,
  col_names = c("Species", "PAR_side", "Telomere_present"),
  show_col_types = FALSE
) %>%
  mutate(
    Species = str_trim(Species),
    PAR_side = str_trim(PAR_side),
    Telomere_present = str_trim(Telomere_present),
    Telomere_present = toupper(Telomere_present),
    Species = factor(Species, levels = species_order),
    species_index = as.integer(Species),
    PAR_order = 0,
    Gene_label = "telomere"
  ) %>%
  filter(!is.na(Species))

telomere_present <- telomere_points %>%
  filter(Telomere_present == "YES")

telomere_absent <- telomere_points %>%
  filter(Telomere_present == "NO")

p_gene_order <- ggplot(display_genes, aes(x = PAR_order, y = species_index)) +
  geom_line(
    aes(group = Species),
    linewidth = 0.4,
    color = "gray70"
  ) +

  # Ortholog segments between adjacent species in tree order
  geom_segment(
    data = ortholog_segments,
    aes(
      x = PAR_order,
      xend = xend,
      y = species_index,
      yend = yend
    ),
    color = "gray50",
    linewidth = 0.5,
    alpha = 0.25,
    inherit.aes = FALSE
  ) +

  # Purple boundaries between PAR/Edge and non-PAR genes
  geom_segment(
    data = par_nonpar_boundaries,
    aes(
      x = boundary_x,
      xend = boundary_x,
      y = y,
      yend = yend
    ),
    color = "purple",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +

  # Telomere present: filled red dot
  geom_point(
    data = telomere_present,
    aes(x = PAR_order, y = species_index),
    color = "red",
    fill = "red",
    shape = 21,
    size = 3,
    inherit.aes = FALSE
  ) +

  # Telomere absent: white circle with red outline
  geom_point(
    data = telomere_absent,
    aes(x = PAR_order, y = species_index),
    color = "red",
    fill = "white",
    shape = 21,
    size = 3,
    stroke = 1,
    inherit.aes = FALSE
  ) +

  # Telomere labels
  geom_text(
    data = telomere_points,
    aes(x = PAR_order, y = species_index, label = Gene_label),
    color = "red",
    angle = 45,
    hjust = 1,
    vjust = -0.5,
    size = 2.8,
    inherit.aes = FALSE
  ) +

  # Non-PAR genes: light gray
  geom_point(
    data = non_par_genes,
    color = "lightgray",
    size = 3
  ) +
  geom_text(
    data = non_par_genes,
    aes(label = Gene_label),
    color = "lightgray",
    angle = 45,
    hjust = 0,
    vjust = -0.5,
    size = 2.8
  ) +

  # Full PAR genes: light-to-dark blue gradient by species frequency
  geom_point(
    data = par_genes_y,
    aes(color = species_frequency),
    size = 3
  ) +
  geom_text(
    data = par_genes_y,
    aes(label = Gene_label, color = species_frequency),
    angle = 45,
    hjust = 0,
    vjust = -0.5,
    size = 2.8
  ) +

  # Edge genes: purple
  geom_point(
    data = par_genes_edge,
    color = "purple",
    size = 3
  ) +
  geom_text(
    data = par_genes_edge,
    aes(label = Gene_label),
    color = "purple",
    angle = 45,
    hjust = 0,
    vjust = -0.5,
    size = 2.8
  ) +

  scale_color_gradient(
    low = "lightblue",
    high = "darkblue",
    name = "PAR gene\nspecies frequency"
  ) +
  scale_x_continuous(
    breaks = sort(unique(c(0, display_genes$PAR_order))),
    expand = expansion(mult = c(0.06, 0.18))
  ) +
  scale_y_continuous(
    breaks = seq_along(species_order),
    labels = species_order,
    limits = c(0.5, length(species_order) + 0.5),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Gene order across PAR and neighboring non-PAR genes",
    y = NULL,
    title = "PAR and non-PAR gene order"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 10, face = "bold")
  )

# ============================================================
# 4. PAR size panel
# Read PAR sizes from species_par.csv
# Format: Species,start-stop
# ============================================================

par_size <- read_csv(
  par_size_file,
  col_names = c("Species", "PAR_interval"),
  show_col_types = FALSE
) %>%
  separate(
    PAR_interval,
    into = c("par_start", "par_stop"),
    sep = "-",
    convert = TRUE
  ) %>%
  mutate(
    par_size_bp = abs(par_stop - par_start)
  ) %>%
  filter(Species %in% species_order) %>%
  mutate(
    Species = factor(Species, levels = species_order)
  )

p_par_size <- ggplot(par_size, aes(x = par_size_bp, y = Species)) +
  geom_col() +
  scale_y_discrete(limits = species_order) +
  scale_x_continuous(
    labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "PAR size",
    y = NULL,
    title = "PAR size"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 10, face = "bold")
  )

# ============================================================
# Combine panels
# Layout:
# top row:    blank tree space | UpSet bars   | blank gene-order | blank PAR-size
# bottom row: tree             | UpSet matrix | gene order       | PAR size
# ============================================================

blank_tree_space <- ggplot() + theme_void()
blank_gene_space <- ggplot() + theme_void()
blank_size_space <- ggplot() + theme_void()

top_row <- blank_tree_space + p_bar + blank_gene_space + blank_size_space +
  plot_layout(widths = c(2.2, 2, 12, 2))

bottom_row <- p_tree + p_matrix + p_gene_order + p_par_size +
  plot_layout(widths = c(2.2, 2, 12, 2))

combined_plot <- top_row / bottom_row +
  plot_layout(heights = c(1.2, 4)) +
  plot_annotation(
    title = plot_title,
    theme = theme(
      plot.title = element_text(size = 14, face = "bold")
    )
  )

# ============================================================
# Save
# ============================================================

ggsave(
  output_file,
  combined_plot,
  width = 24,
  height = max(8, length(species_order) * 0.35),
  limitsize = FALSE
)

message("Wrote: ", output_file)