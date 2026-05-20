library(tidyverse)
library(ape)
library(ggtree)
library(patchwork)
library(scales)

# ============================================================
# Inputs
# ============================================================

gene_file <- "gene_locations_by_species.with_chr_label.all.In_PAR.Z_only.gapless_species.csv"
telomere_file <- "gapless_species.telomeres.txt"
tree_file <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/roadies_v1.1.16b.numbers.scientific.nwk"

df <- read_csv(gene_file)
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
    PAR_order = row_number(),
    Gene_label = if_else(
      duplicated(Gene) | duplicated(Gene, fromLast = TRUE),
      paste0(Gene, "_hit", Hit_rank),
      Gene
    )
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

# ------------------------------------------------------------
# Telomere points
# telomere file format:
# Species,Side,Telomere_present
# Example:
# Lathamus_discolor,R,NO
# ------------------------------------------------------------

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
    PAR_order = 0,
    Gene_label = "telomere"
  ) %>%
  filter(!is.na(Species))

telomere_present <- telomere_points %>%
  filter(Telomere_present == "YES")

telomere_absent <- telomere_points %>%
  filter(Telomere_present == "NO")

par_genes_y <- par_genes %>%
  filter(In_PAR == "Y")

par_genes_edge <- par_genes %>%
  filter(In_PAR == "Edge")
  
p_gene_order <- ggplot(par_genes, aes(x = PAR_order, y = Species)) +
  geom_line(aes(group = Species), linewidth = 0.4, color = "gray70") +

  # Telomere present: filled red dot
  geom_point(
    data = telomere_present,
    aes(x = PAR_order, y = Species),
    color = "red",
    fill = "red",
    shape = 21,
    size = 3,
    inherit.aes = FALSE
  ) +

  # Telomere absent: white circle with red outline
  geom_point(
    data = telomere_absent,
    aes(x = PAR_order, y = Species),
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
    aes(x = PAR_order, y = Species, label = Gene_label),
    color = "red",
    angle = 45,
    hjust = 1,
    vjust = -0.5,
    size = 2.8,
    inherit.aes = FALSE
  ) +

  # PAR genes
  geom_point(aes(color = species_frequency), size = 3) +
  geom_text(
    aes(label = Gene_label, color = species_frequency),
    angle = 45,
    hjust = 0,
    vjust = -0.5,
    size = 2.8
  ) +
  scale_color_gradient(
    low = "lightgray",
    high = "black",
    name = "PAR gene\nspecies frequency"
  ) +
  scale_x_continuous(
    breaks = sort(unique(c(0, par_genes$PAR_order))),
    expand = expansion(mult = c(0.06, 0.18))
  ) +
  scale_y_discrete(limits = species_order) +
  labs(
    x = "Gene order within PAR",
    y = NULL,
    title = "PAR gene order"
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
# 4. PAR size panel
# Read PAR sizes from species_par.csv
# Format: Species,start-stop
# ============================================================

par_size_file <- "species_par.csv"

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
# top row:    blank tree space | UpSet bars | blank gene-order | blank PAR-size
# bottom row: tree             | UpSet matrix | gene order      | PAR size
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
    title = "Phylogeny, PAR gene intersections, PAR gene order, and PAR size",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold")
    )
  )


ggsave(
  "combined_phylogeny_upset_gene_order_PAR_size.pdf",
  combined_plot,
  width = 24,
  height = max(8, length(species_order) * 0.35),
  limitsize = FALSE
)