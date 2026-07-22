library(tidyverse)
library(ape)
library(ggtree)
library(patchwork)
library(scales)

# ============================================================
# Inputs
# ============================================================

gene_file <- "birds_PAR_genes.all.tsv"

tree_file <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/roadies_v1.1.16b.numbers.scientific.nwk"

par_size_file <- "PAR.species_chr_region.gaps.txt"

telomere_file <- "species.telomeres.txt"

df <- read_tsv(
  gene_file,
  col_types = cols(
    Species = col_character(),
    TOGADir = col_character(),
    Chromosome = col_character(),
    StartPos = col_double(),
    StopPos = col_double(),
    GeneName = col_character(),
    InPAR = col_character()
  )
)
tree <- read.tree(tree_file)

species_par <- read_csv(
  par_size_file,
  col_names = c("Species", "PAR"),
  show_col_types = FALSE
) %>%
  mutate(
    Species = as.character(Species),
    PAR = str_trim(as.character(PAR))
  ) %>%
  separate(
    PAR,
    into = c("CHROM", "PAR"),
    sep = ":",
    convert = TRUE
  ) %>%
  separate(
    PAR,
    into = c("PARStart", "PARStop"),
    sep = "-",
    convert = TRUE
  ) %>%
  mutate(
    PARStart = as.numeric(PARStart),
    PARStop = as.numeric(PARStop),
    par_size_bp = abs(PARStop - PARStart)
  )

# ============================================================
# Normalize column names from birds_PAR_genes.tsv
# ============================================================

df <- df %>%
  rename(
    Gene = GeneName,
    Chrom = Chromosome,
    Start_pos = StartPos,
    Stop_pos = StopPos,
    In_PAR = InPAR
  ) %>%
  mutate(
    Start_pos = as.numeric(Start_pos),
    Stop_pos = as.numeric(Stop_pos),
    Gene = as.character(Gene),
    In_PAR = as.character(In_PAR)
  )

# ============================================================
# Drop genes excluded from all plots/analyses
# ============================================================

genes_to_drop <- c("LOC121469005")

df <- df %>%
  filter(!Gene %in% genes_to_drop)

# ============================================================
# Keep genes that are in the PAR in at least one species
# ============================================================

genes_in_PAR_any_species <- df %>%
  filter(In_PAR %in% c("Y", "Edge", "N")) %>%
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


tree_filtered <- ape::rotate(tree_filtered, node = 28)
tree_filtered <- ape::rotate(tree_filtered, node = 29)
tree_filtered <- ape::rotate(tree_filtered, node = 30)
tree_filtered <- ape::rotate(tree_filtered, node = 31)
tree_filtered <- ape::rotate(tree_filtered, node = 32)
tree_filtered <- ape::rotate(tree_filtered, node = 33)
tree_filtered <- ape::rotate(tree_filtered, node = 34)
tree_filtered <- ape::rotate(tree_filtered, node = 35)

p_tree_tmp <- ggtree(tree_filtered, ladderize = FALSE)

tree_plot_order <- p_tree_tmp$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)

species_order <- tree_plot_order

par_binary_tree <- par_binary %>%
  select(Gene, all_of(species_order))

df_par_relevant <- df_par_relevant %>%
  filter(Species %in% species_order)

species_par <- species_par %>%
  filter(Species %in% species_order)

# ============================================================
# 1. Phylogeny panel
# ============================================================

p_tree <- ggtree(tree_filtered, ladderize = FALSE) +
  geom_tiplab(size = 3, align = FALSE) +
  xlim_tree(0.4) +
  coord_cartesian(clip = "off") +
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
  mutate(
    n_species = if_else(
      intersection_id == "",
      0L,
      str_count(intersection_id, fixed("|")) + 1L
    )
  ) %>%
  count(intersection_id, n_species, name = "n_genes") %>%
  arrange(desc(n_species), desc(n_genes), intersection_id) %>%
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
  geom_point(size = 1) +
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
  filter(In_PAR %in% c("Y", "Edge", "N")) %>%
  filter(Species %in% species_order) %>%
  mutate(
    midpoint = (Start_pos + Stop_pos) / 2
  ) %>%
  group_by(Species, Chrom) %>%
  mutate(
    PAR_at_end = min(Start_pos[In_PAR %in% c("Y", "Edge")], na.rm = TRUE) >= 10000000,
    adjusted_pos = if_else(PAR_at_end, -midpoint, midpoint)
  ) %>%
  arrange(Species, Chrom, adjusted_pos) %>%
  mutate(
    PAR_order = row_number()
  ) %>%
  ungroup() %>%
  filter(PAR_order <= 90)

gene_freq <- par_genes %>%
  filter(In_PAR == "Y") %>%
  distinct(Species, Gene) %>%
  count(Gene, name = "species_frequency")

par_genes <- par_genes %>%
  left_join(gene_freq, by = "Gene") %>%
  mutate(
    Species = factor(Species, levels = species_order),
    species_index = as.integer(Species)
  ) %>%
  group_by(Species, Gene) %>%
  mutate(
    Hit_rank = row_number(),
    n_hits = n(),
    Gene_label = if_else(
      n_hits > 1,
      paste0(Gene, "_hit", Hit_rank),
      Gene
    )
  ) %>%
  ungroup()

par_genes <- par_genes %>%
  mutate(
    is_array_gene = str_detect(Gene, regex("array", ignore_case = TRUE))
  )

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
  
nonpar_genes <- par_genes %>%
  filter(In_PAR == "N", !is_array_gene)

nonpar_array_genes <- par_genes %>%
  filter(In_PAR == "N", is_array_gene)

par_genes_y <- par_genes %>%
  filter(In_PAR == "Y", !is_array_gene)

par_array_genes_y <- par_genes %>%
  filter(In_PAR == "Y", is_array_gene)

par_genes_edge <- par_genes %>%
  filter(In_PAR == "Edge", !is_array_gene)

par_array_genes_edge <- par_genes %>%
  filter(In_PAR == "Edge", is_array_gene)

ref_species_labels <- par_genes %>%
  filter(Species == "Taeniopygia_guttata")

ortholog_segments <- par_genes %>%

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

p_gene_order <- ggplot(par_genes, aes(x = PAR_order, y = species_index)) +
  geom_line(
    aes(group = Species),
    linewidth = 0.4,
    color = "gray70"
  ) +
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
  geom_point(
    data = nonpar_genes,
    color = "gray80",
    size = 1.2
  ) +
  geom_point(
    data = nonpar_array_genes,
    color = "gray80",
    shape = 24,
    size = 2
  ) +
  geom_point(
    data = par_genes_y,
    aes(color = species_frequency),
    size = 1.5
  ) +
  geom_point(
    data = par_array_genes_y,
    aes(color = species_frequency),
    shape = 24,
    size = 2
  ) +
  geom_point(
    data = par_genes_edge,
    color = "blue",
    size = 1.5
  ) +
  geom_point(
    data = par_array_genes_edge,
    color = "blue",
    shape = 24,
    size = 2
  ) +
  geom_text(
    data = ref_species_labels,
    aes(label = Gene_label),
    color = "black",
    angle = 45,
    hjust = 0,
    vjust = -0.5,
    nudge_y = 0.08,
    size = 2.5
  ) +
  scale_color_gradient(
    low = "#fcbba1",
    high = "#cb181d",
    name = "PAR gene\nspecies frequency"
  ) +
  geom_point(
  data = telomere_present,
  aes(x = PAR_order, y = species_index),
  inherit.aes = FALSE,
  shape = 23,
  fill = "black",
  color = "black",
  size = 2
) +
  geom_point(
    data = telomere_absent,
    aes(x = PAR_order, y = species_index),
    inherit.aes = FALSE,
    shape = 23,
    fill = "white",
    color = "black",
    size = 2
  ) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),
    limits = c(0, 90)
  ) +
  scale_y_continuous(
    breaks = seq_along(species_order),
    labels = species_order,
    limits = c(0.5, length(species_order) + 0.5),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
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
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 10, face = "bold")
  )

# ============================================================
# 4. PAR size panel
# Use PAR coordinates from species_par.csv
# ============================================================

par_size <- species_par %>%
  mutate(
    Species = factor(Species, levels = species_order)
  ) %>%
  filter(!is.na(Species))

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
# ============================================================

blank_tree_space <- ggplot() + theme_void()
blank_gene_space <- ggplot() + theme_void()
blank_size_space <- ggplot() + theme_void()

top_row <- blank_tree_space + p_bar + blank_gene_space + blank_size_space +
  plot_layout(widths = c(3, 2.5, 20, 2))

bottom_row <- p_tree + p_matrix + p_gene_order + p_par_size +
  plot_layout(widths = c(3, 2.5, 20, 2))

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
  height = 8,
  limitsize = FALSE
)

