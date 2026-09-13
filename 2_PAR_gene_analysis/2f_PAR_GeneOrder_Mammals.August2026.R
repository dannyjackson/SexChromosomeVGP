library(tidyverse)
library(ape)
library(ggtree)
library(patchwork)
library(scales)

# ============================================================
# Inputs
# ============================================================

# Gene to use as ordinal origin point
GENE <- "SHROOM2"

XG_ORIGIN_SPECIES <- c(
  "Pan_paniscus",
  "Pan_troglodytes",
  "Homo_sapiens",
  "Gorilla_gorilla",
  "Pongo_abelii",
  "Pongo_pygmaeus",
  "Symphalangus_syndactylus",
  "Macaca_nemestrina"
)

FLIP_NO_PAR_SPECIES <- c( 
  "Thomomys_bottae", 
  "Microtus_pennsylvanicus", 
  "Nyctalus_leisleri", 
  "Myotis_mystacinus")


XG_ORIGIN_SPECIES <- c(
)

# How many ordinal positions to show on either side of GENE
ORDINAL_MIN <- -45
ORDINAL_MAX <- 37

gene_file <- "mammals_PAR_genes.all.ZNF_arrays.tsv"

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
# Normalize column names from mammals_PAR_genes.tsv
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
# Keep genes that are in the PAR in at least one species
# ============================================================

genes_in_PAR_any_species <- df %>%
  filter(!str_detect(Gene, regex("array", ignore_case = TRUE))) %>%
  filter(Species %in% c("Homo_sapiens")) %>%
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

tree_filtered <- ape::rotate(tree_filtered, node = 42)
tree_filtered <- ape::rotate(tree_filtered, node = 71)
tree_filtered <- ape::rotate(tree_filtered, node = 72)
tree_filtered <- ape::rotate(tree_filtered, node = 73)
tree_filtered <- ape::rotate(tree_filtered, node = 74)
tree_filtered <- ape::rotate(tree_filtered, node = 75)


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
  geom_tiplab(size = 5, align = FALSE) +
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
# TEST: Diagnostic tree with node numbers for choosing rotations
# ============================================================

test_tree_png <- "tree_filtered_node_numbers.png"

p_tree_nodes <- ggtree(tree_filtered, ladderize = FALSE) +
  geom_tree() +
  
  # Tip labels
  geom_tiplab(size = 3, align = FALSE) +
  
  # Internal node numbers
  geom_text2(
    aes(label = node, subset = !isTip),
    hjust = -0.3,
    vjust = -0.3,
    size = 3,
    color = "red"
  ) +
  
  # Optional: tip node numbers too, useful for debugging
  geom_text2(
    aes(label = node, subset = isTip),
    hjust = 1.2,
    vjust = -0.4,
    size = 2.5,
    color = "blue"
  ) +
  
  xlim_tree(0.6) +
  coord_cartesian(clip = "off") +
  theme_tree2() +
  theme(
    plot.margin = margin(5.5, 80, 5.5, 5.5)
  )

ggsave(
  filename = test_tree_png,
  plot = p_tree_nodes,
  width = 10,
  height = max(6, length(tree_filtered$tip.label) * 0.22),
  dpi = 300
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
    title = NULL
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
# Gene order is centered on GENE, so GENE = ordinal position 0
# ============================================================


par_orientation <- df %>%
  group_by(Species, Chrom) %>%
  mutate(
    chr_start = 0,
    chr_end = max(Stop_pos, na.rm = TRUE)
  ) %>%
  filter(Gene == GENE) %>%
  summarise(
    chr_start = first(chr_start),
    chr_end = first(chr_end),

    origin_midpoint = first(
      (Start_pos + Stop_pos) / 2
    ),

    distance_to_start = origin_midpoint - chr_start,
    distance_to_end = chr_end - origin_midpoint,

    # Flip when SHROOM2 is closer to the chromosome end
    flip_orientation = distance_to_end < distance_to_start,

    .groups = "drop"
  )

par_genes <- df %>%
  filter(Gene %in% genes_in_PAR_any_species) %>%
  filter(!str_detect(Gene, regex("array", ignore_case = TRUE))) %>%
  filter(Species %in% species_order) %>%
  mutate(
    midpoint = (Start_pos + Stop_pos) / 2
  ) %>%
  left_join(
    par_orientation,
    by = c("Species", "Chrom")
  ) %>%
  mutate(
    flip_orientation = coalesce(flip_orientation, FALSE),

    adjusted_pos = if_else(
      flip_orientation,
      -midpoint,
      midpoint
    )
  ) %>%
  arrange(Species, Chrom, adjusted_pos) %>%
  group_by(Species, Chrom) %>%
  mutate(
    raw_PAR_order = row_number()
  ) %>%
  ungroup()

# ------------------------------------------------------------
# Find the origin gene in each species
# If there are multiple hits for GENE in a species, use the first
# ordinal occurrence after chromosome/PAR orientation adjustment.
# ------------------------------------------------------------

origin_gene_by_species <- tibble(
  Species = species_order,
  Origin_gene = if_else(
    Species %in% XG_ORIGIN_SPECIES,
    "XG",
    GENE
  )
)

gene_origin <- par_genes %>%
  inner_join(origin_gene_by_species, by = "Species") %>%
  filter(Gene == Origin_gene) %>%
  group_by(Species, Origin_gene) %>%
  slice_min(raw_PAR_order, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    origin_raw_PAR_order = if_else(
      Species %in% XG_ORIGIN_SPECIES,
      raw_PAR_order + 24,
      raw_PAR_order
    )
  ) %>%
  select(
    Species,
    Origin_gene,
    origin_raw_PAR_order
  )

missing_origin_species <- setdiff(
  as.character(species_order),
  as.character(gene_origin$Species)
)

if (length(missing_origin_species) > 0) {
  missing_origin_tbl <- origin_gene_by_species %>%
    filter(Species %in% missing_origin_species)

  warning(
    paste0(
      "Origin gene was not found in these species and they will be omitted from the gene-order panel: ",
      paste(
        paste0(
          missing_origin_tbl$Species,
          " expected ",
          missing_origin_tbl$Origin_gene
        ),
        collapse = ", "
      )
    )
  )
}

par_genes <- par_genes %>%
  inner_join(gene_origin, by = "Species") %>%
  mutate(
    PAR_order = raw_PAR_order - origin_raw_PAR_order
  ) %>%
  filter(
    PAR_order >= ORDINAL_MIN,
    PAR_order <= ORDINAL_MAX
  )

plot_x_min <- min(par_genes$PAR_order, na.rm = TRUE) - 2
plot_x_max <- plot_x_min + (ORDINAL_MAX - ORDINAL_MIN)

gene_freq <- par_genes %>%
  filter(In_PAR %in% c("Y", "Edge")) %>%
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
  ungroup() %>%
  mutate(
    is_array_gene = str_detect(Gene, regex("array", ignore_case = TRUE))
  )

# ------------------------------------------------------------
# Telomere points are also shifted relative to the origin gene.
# Previously telomere was plotted at raw ordinal position 0.
# Now telomere position = 0 - origin_raw_PAR_order.
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
    Telomere_present = toupper(Telomere_present)
  ) %>%
  inner_join(gene_origin, by = "Species") %>%
  mutate(
    Species = factor(Species, levels = species_order),
    species_index = as.integer(Species),
    PAR_order = 0 - origin_raw_PAR_order,
    Gene_label = "telomere"
  ) %>%
  filter(
    !is.na(Species),
    PAR_order >= ORDINAL_MIN,
    PAR_order <= ORDINAL_MAX
  )

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
  filter(Species == "Homo_sapiens")

# ------------------------------------------------------------
# Label the last PAR/Edge gene in each species
# "Last" means the greatest plotted PAR_order after any inversion.
# ------------------------------------------------------------

last_par_gene_labels <- par_genes %>%
  filter(
    In_PAR %in% c("Y", "Edge"),
    !is_array_gene
  ) %>%
  group_by(Species) %>%
  slice_max(
    order_by = PAR_order,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()
  
# ============================================================
# Gene-level PAR/Edge conservation panel
# Bar above each Homo sapiens gene = number of genomes where
# that gene is found in PAR or Edge
# ============================================================

gene_par_edge_counts <- df %>%
  filter(Gene %in% genes_in_PAR_any_species) %>%
  filter(!str_detect(Gene, regex("array", ignore_case = TRUE))) %>%
  filter(In_PAR %in% c("Y", "Edge")) %>%
  distinct(Species, Gene) %>%
  count(Gene, name = "n_genomes_PAR_or_Edge")

homo_gene_par_edge_counts <- ref_species_labels %>%
  distinct(Gene, Gene_label, PAR_order) %>%
  left_join(gene_par_edge_counts, by = "Gene") %>%
  mutate(
    n_genomes_PAR_or_Edge = replace_na(n_genomes_PAR_or_Edge, 0L)
  )

p_gene_par_edge_count <- ggplot(
  homo_gene_par_edge_counts,
  aes(x = PAR_order, y = n_genomes_PAR_or_Edge)
) +
  geom_col(width = 0.8) +
  scale_x_continuous(
    breaks = seq(
      ceiling(plot_x_min / 10) * 10,
      floor(plot_x_max / 10) * 10,
      by = 10
    ),
    limits = c(plot_x_min, plot_x_max),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 5.5, 0, 5.5)
  )
  
ortholog_segments <- par_genes %>%
  group_by(Species, Gene) %>%
  slice_min(abs(PAR_order), n = 1, with_ties = FALSE) %>%
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

highlight_segments <- par_genes %>%
  filter(Gene %in% c("SHROOM2", "XG")) %>%
  group_by(Species, Gene) %>%
  slice_min(abs(PAR_order), n = 1, with_ties = FALSE) %>%
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
  geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    linetype = "dashed",
    color = "gray40"
  ) +
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
  geom_segment(
    data = highlight_segments,
    aes(
      x = PAR_order,
      xend = xend,
      y = species_index,
      yend = yend
    ),
    color = "black",
    linewidth = 0.5,
    alpha = 1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = nonpar_genes,
    color = "#f1e1f9",
    size = 3
  ) +
  geom_point(
    data = nonpar_array_genes,
    color = "#f1e1f9",
    shape = 24,
    size = 2
  ) +
  geom_point(
    data = par_genes_y,
    aes(color = species_frequency),
    size = 3
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
    size = 3
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
    angle = 90,
    hjust = 0,
    vjust = 0.5,
    nudge_y = 0.15,
    size = 4
  ) +
  geom_text(
    data = last_par_gene_labels,
    aes(
      x = PAR_order,
      y = species_index,
      label = Gene_label
    ),
    inherit.aes = FALSE,
    color = "black",
    hjust = -0.15,
    vjust = 1.5,
    size = 4
  ) +
  scale_color_gradient(
    low = "#8a65b9",
    high = "#8a65b9",
    name = "PAR gene\nspecies frequency"
  ) +
  geom_point(
    data = telomere_present,
    aes(x = PAR_order, y = species_index),
    inherit.aes = FALSE,
    shape = 23,
    fill = "black",
    color = "black",
    size = 3
  ) +
  geom_point(
    data = telomere_absent,
    aes(x = PAR_order, y = species_index),
    inherit.aes = FALSE,
    shape = 23,
    fill = "white",
    color = "black",
    size = 3
  ) +
  scale_x_continuous(
    breaks = seq(
      ceiling(plot_x_min / 10) * 10,
      floor(plot_x_max / 10) * 10,
      by = 10
    ),
    limits = c(plot_x_min, plot_x_max),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq_along(species_order),
    labels = species_order,
    limits = c(0.5, length(species_order) + 0.5),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = paste0("Ordinal gene order"),
    y = NULL,
    title = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(
      t = 60,
      r = 10,
      b = 5.5,
      l = 5.5
    ),
    panel.border = element_blank()
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
    title = NULL
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

top_row <- blank_tree_space + blank_size_space + p_gene_par_edge_count +
  plot_layout(widths = c(3, 2, 20))

bottom_row <- p_tree + p_par_size + p_gene_order +
  plot_layout(widths = c(3, 2, 20))

combined_plot <- top_row / bottom_row +
  plot_layout(heights = c(0.5, 4)) +
  plot_annotation(
    title = "Phylogeny, PAR gene intersections, PAR gene order, and PAR size",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold")
    )
  )

ggsave(
  filename = paste0("combined_phylogeny_upset_gene_order_PAR_size.SHROOM_XG.conserved_genes.pdf"),
  plot = combined_plot,
  width = 24,
  height = 18,
  limitsize = FALSE
)










# ============================================================
# PAR size vs PAR gene count
# ALL MAMMALS -- ONE DATASET
#
# Color = taxonomic Order
# Shape = T2T status
#   triangle = T2T
#   circle   = non-T2T
#
# Models:
#   1. Ordinary LM
#   2. PGLS with estimated Pagel's lambda
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ape)
library(nlme)
library(ggplot2)
library(ggrepel)


# ============================================================
# 1. Define species metadata
# ============================================================

T2T <- c(
  "Gorilla_gorilla",
  "Pan_paniscus",
  "Pan_troglodytes",
  "Pongo_abelii",
  "Pongo_pygmaeus",
  "Symphalangus_syndactylus",
  "Ovis_canadensis",
  "Ovis_aries",
  "Capra_hircus",
  "Callithrix_jacchus"
)


species_orders <- tibble::tribble(
  ~Species,                     ~Order,
  "Balaenoptera_physalus",      "Cetartiodactyla",
  "Bos_taurus",                 "Cetartiodactyla",
  "Callithrix_jacchus",         "Primates",
  "Camelus_dromedarius",        "Cetartiodactyla",
  "Capra_hircus",               "Cetartiodactyla",
  "Eubalaena_glacialis",        "Cetartiodactyla",
  "Gorilla_gorilla",            "Primates",
  "Grampus_griseus",            "Cetartiodactyla",
  "Homo_sapiens",               "Primates",
  "Inia_geoffrensis",           "Cetartiodactyla",
  "Loxodonta_africana",         "Proboscidea",
  "Lycaon_pictus",              "Carnivora",
  "Macaca_nemestrina",          "Primates",
  "Manis_pentadactyla",         "Pholidota",
  "Marmota_flaviventris",       "Rodentia",
  "Meles_meles",                "Carnivora",
  "Mesoplodon_bidens",          "Cetartiodactyla",
  "Molossus_nigricans",         "Chiroptera",
  "Mustela_nivalis_vulgaris",   "Carnivora",
  "Myotis_nattereri",           "Chiroptera",
  "Ovis_aries",                 "Cetartiodactyla",
  "Ovis_canadensis",            "Cetartiodactyla",
  "Pan_paniscus",               "Primates",
  "Pan_troglodytes",            "Primates",
  "Panthera_onca",              "Carnivora",
  "Pongo_abelii",               "Primates",
  "Pongo_pygmaeus",             "Primates",
  "Pseudorca_crassidens",       "Cetartiodactyla",
  "Rhynchocyon_petersi",        "Macroscelidea",
  "Rhynchonycteris_naso",       "Chiroptera",
  "Symphalangus_syndactylus",   "Primates",
  "Trichechus_inunguis",        "Sirenia",
  "Urocitellus_parryii",        "Rodentia",
  "Canis_lupus_baileyi",        "Carnivora",
  "Nyctalus_leisleri",          "Chiroptera",
  "Dasypus_novemcinctus",       "Cingulata",
  "Artibeus_intermedius",       "Chiroptera",
  "Artibeus_lituratus",         "Chiroptera",
  "Myotis_mystacinus",          "Chiroptera",
  "Corynorhinus_townsendii",    "Chiroptera",
  "Miniopterus_schreibersii",   "Chiroptera"
)


# ============================================================
# Common plotting aesthetics
# ============================================================

point_shapes <- c(
  "Non-T2T" = 21,  # fillable circle
  "T2T"     = 24   # fillable triangle
)

# Fixed Order levels
order_levels <- c(
  "Cetartiodactyla",
  "Primates",
  "Carnivora",
  "Chiroptera",
  "Rodentia",
  "Proboscidea",
  "Pholidota",
  "Macroscelidea",
  "Sirenia",
  "Cingulata"
)


# ============================================================
# 2. Build ONE all-mammal dataset
# ============================================================

par_gene_counts <- df %>%
  filter(
    Gene %in% genes_in_PAR_any_species,
    In_PAR %in% c("Y", "Edge"),
    !str_detect(Gene, regex("array", ignore_case = TRUE))
  ) %>%
  distinct(Species, Gene) %>%
  count(Species, name = "par_gene_count")


model_data <- species_par %>%
  transmute(
    Species,
    par_size_bp = PARStop - PARStart,
    par_size_Mb = (PARStop - PARStart) / 1e6
  ) %>%
  left_join(par_gene_counts, by = "Species") %>%
  left_join(species_orders, by = "Species") %>%
  mutate(
    par_gene_count = replace_na(par_gene_count, 0L),

    T2T_status = if_else(
      Species %in% T2T,
      "T2T",
      "Non-T2T"
    ),

    # Log variables use the SAME species dataset
    log_par_size = log(par_size_Mb),
    log_par_gene_count = log1p(par_gene_count)
  ) %>%
  filter(
    is.finite(par_size_Mb),
    par_size_Mb > 0,
    is.finite(par_gene_count)
  ) %>%
  distinct(Species, .keep_all = TRUE)


# Make sure every species has an Order
if (any(is.na(model_data$Order))) {
  warning(
    "Missing Order for: ",
    paste(
      model_data$Species[is.na(model_data$Order)],
      collapse = ", "
    )
  )
}


# Use identical factor ordering everywhere
model_data$Order <- factor(
  model_data$Order,
  levels = order_levels
)

# ============================================================
# 3. Match ONE dataset to ONE phylogeny
# ============================================================

species_for_model <- intersect(
  tree_filtered$tip.label,
  model_data$Species
)

model_tree <- keep.tip(
  tree_filtered,
  species_for_model
)

model_data <- model_data %>%
  filter(Species %in% model_tree$tip.label) %>%
  arrange(match(Species, model_tree$tip.label))


stopifnot(
  identical(
    as.character(model_data$Species),
    model_tree$tip.label
  )
)


model_data <- as.data.frame(model_data)
rownames(model_data) <- model_data$Species


# ============================================================
# x = PAR gene count
# y = PAR size
#
# Color = Order
# Shape:
#   triangle = T2T
#   circle   = Non-T2T
# ============================================================


# ============================================================
# 4. RAW DATA MODELS
# ============================================================


# ----------------------------
# Ordinary linear model
# ----------------------------

lm_raw <- lm(
  par_size_Mb ~ par_gene_count,
  data = model_data
)

summary(lm_raw)


# ----------------------------
# PGLS
# ----------------------------
library(phylolm)

lambda_raw <- corPagel(
  value = 0.5,
  phy = model_tree,
  fixed = FALSE,
  form = ~ Species
)


pgls_raw <- phylolm(
  par_size_Mb ~ par_gene_count,
  data = model_data,
  phy = model_tree,
  model = "lambda"
)

summary(pgls_raw)



# ============================================================
# 5. Add residuals
# ============================================================

model_data$lm_raw_resid <- rstandard(
  lm_raw
)

model_data$pgls_raw_resid <- residuals(
  pgls_raw
  )


# ============================================================
# 6. Common plotting aesthetics
# ============================================================

point_shapes <- c(
  "Non-T2T" = 21,  # fillable circle
  "T2T"     = 24   # fillable triangle
)



# ============================================================
# 7. Define Order levels from the ALL-MAMMAL dataset
# ============================================================

mammal_order_levels <- model_data %>%
  distinct(Order) %>%
  filter(!is.na(Order)) %>%
  arrange(Order) %>%
  pull(Order) %>%
  as.character()


# ============================================================
# 5. Set identical factor levels in both datasets
# ============================================================

model_data <- model_data %>%
  mutate(
    Order = factor(
      Order,
      levels = mammal_order_levels
    )
  )


# ============================================================
# 6. Color-blind-friendly Order colors
# ============================================================

mammal_order_colors <- setNames(
  viridisLite::viridis(
    length(mammal_order_levels),
    option = "D",
    begin = 0.05,
    end = 0.95
  ),
  mammal_order_levels
)


# ============================================================
# 7. Alternate filled/open circles
# ============================================================

mammal_order_shapes <- setNames(
  rep(
    c(16, 1),
    length.out = length(mammal_order_levels)
  ),
  mammal_order_levels
)



mammal_order_fills <- mammal_order_colors

# Every second Order is open
mammal_order_fills[seq(2, length(mammal_order_fills), by = 2)] <- "white"

# ============================================================
# 8. Sanity checks
# ============================================================

table(model_data$Order, useNA = "ifany")


# ============================================================
# 8. RAW DATA -- ORDINARY LM
# ============================================================

p_lm_raw <- ggplot(
  model_data,
  aes(
    x = par_gene_count,
    y = par_size_Mb
  )
) +
  geom_point(
    aes(
      color = Order,
      fill = Order,
      shape = T2T_status
    ),
    size = 4,
    stroke = 1.2
  ) +

  scale_color_manual(
    values = mammal_order_colors,
    limits = mammal_order_levels,
    drop = FALSE
  ) +

  scale_fill_manual(
    values = mammal_order_fills,
    limits = mammal_order_levels,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  # LM line + 95% confidence ribbon
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black"
  ) +

  geom_text_repel(
    aes(label = Species),
    size = 3,
    color = "black",
    show.legend = FALSE
  ) +

  scale_color_manual(
    values = mammal_order_colors,
    limits = mammal_order_levels,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  labs(
    title = "Ordinary linear model",
    x = "Number of genes in the PAR",
    y = "PAR size (Mb)",
    color = "Order",
    shape = "Assembly"
  ) +


  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,
        shape = mammal_order_shapes
      )
    ),
    shape = "none"
  ) +

  theme_bw() +

  theme(
    legend.position = "right",
    legend.key.height = unit(0.45, "cm"),
    legend.text = element_text(size = 8)
  )






# ============================================================
# 9. RAW DATA -- PGLS
# ============================================================

# Add all-mammal PGLS residuals
model_data$pgls_resid <- residuals(
  pgls_raw
)


p_pgls_raw <- ggplot(
  model_data,
  aes(
    x = par_gene_count,
    y = pgls_resid
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.7
  ) +

  geom_point(
    aes(
      color = Order,
      fill = Order,
      shape = T2T_status
    ),
    size = 4,
    stroke = 1.2
  ) +

  scale_color_manual(
    values = mammal_order_colors,
    limits = mammal_order_levels,
    drop = FALSE
  ) +

  scale_fill_manual(
    values = mammal_order_fills,
    limits = mammal_order_levels,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  geom_text_repel(
    aes(label = Species),
    size = 3,
    color = "black",
    show.legend = FALSE
  ) +

  scale_color_manual(
    values = mammal_order_colors,
    limits = mammal_order_levels,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  labs(
    title = paste0(
      "PGLS residuals: All mammals; lambda = ",
      round(pgls_raw$optpar, 2)
    ),
    x = "PAR size (Mb)",
    y = "PGLS residual",
    color = "Order",
    shape = "Assembly"
  ) +

  theme_bw()



# ============================================================
# 16. Save plots
# ============================================================

# PNGs

ggsave(
  "LM.AllMammals.png",
  p_lm_raw,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "PGLS.AllMammals.png",
  p_pgls_raw,
  width = 8,
  height = 6,
  dpi = 300
)


# PDFs

ggsave(
  "LM.AllMammals.pdf",
  p_lm_raw,
  width = 8,
  height = 6
)

ggsave(
  "PGLS.AllMammals.pdf",
  p_pgls_raw,
  width = 8,
  height = 6
)

# ============================================================
# ANOVA of Order by size
# ============================================================

lm_order <- lm(
  par_size_Mb ~ Order,
  data = model_data
)

summary(lm_order)


# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================



# ============================================================
# Drop Primates and Rodentia
# ============================================================

model_data_reduced <- model_data %>%
  filter(
    !Order %in% c("Primates", "Rodentia")
  )


# ============================================================
# Re-match reduced dataset to phylogeny
# ============================================================

species_for_model <- intersect(
  tree_filtered$tip.label,
  model_data_reduced$Species
)

model_tree_reduced <- keep.tip(
  tree_filtered,
  species_for_model
)

model_data_reduced <- model_data_reduced %>%
  filter(Species %in% model_tree_reduced$tip.label) %>%
  arrange(match(Species, model_tree_reduced$tip.label))

stopifnot(
  identical(
    as.character(model_data_reduced$Species),
    model_tree_reduced$tip.label
  )
)

model_data_reduced <- as.data.frame(model_data_reduced)
rownames(model_data_reduced) <- model_data_reduced$Species


# ============================================================
# RAW MODELS
# ============================================================

# Ordinary LM
lm_raw_reduced <- lm(
  par_size_Mb ~ par_gene_count,
  data = model_data_reduced
)

summary(lm_raw_reduced)


library(phylolm)
# estimate λ with phylolm

pgls_raw_reduced <- phylolm(
  par_size_Mb ~ par_gene_count,
  data = model_data_reduced,
  phy = model_tree_reduced,
  model = "lambda"
)

summary(pgls_raw_reduced)


# ============================================================
# Add residuals
# ============================================================

model_data_reduced$lm_raw_resid <- rstandard(
  lm_raw_reduced
)

model_data_reduced$pgls_raw_resid <- residuals(
  pgls_raw_reduced
  )


# ============================================================
# 8. RAW DATA -- ORDINARY LM
# ============================================================

p_lm_reduced <- ggplot(
  model_data_reduced,
  aes(
    x = par_gene_count,
    y = par_size_Mb
  )
) +
  geom_point(
    aes(
      color = Order,
      fill = Order,
      shape = T2T_status
    ),
    size = 4,
    stroke = 1.2
  ) +

  scale_color_manual(
    values = mammal_order_colors,
    limits = mammal_order_levels,
    drop = FALSE
  ) +

  scale_fill_manual(
    values = mammal_order_fills,
    limits = mammal_order_levels,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  # LM line + 95% confidence ribbon
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black"
  ) +

  geom_text_repel(
    aes(label = Species),
    size = 3,
    color = "black",
    show.legend = FALSE
  ) +

  scale_color_manual(
    values = mammal_order_colors,
    limits = mammal_order_levels,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  labs(
    title = "Ordinary linear model",
    x = "Number of genes in the PAR",
    y = "PAR size (Mb)",
    color = "Order",
    shape = "Assembly"
  ) +


  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,
        shape = mammal_order_shapes
      )
    ),
    shape = "none"
  ) +

  theme_bw() +

  theme(
    legend.position = "right",
    legend.key.height = unit(0.45, "cm"),
    legend.text = element_text(size = 8)
  )






# ============================================================
# 9. RAW DATA -- PGLS
# ============================================================

# Add all-mammal PGLS residuals
model_data_reduced$pgls_resid <- residuals(
  pgls_raw_reduced
)


p_pgls_reduced <- ggplot(
  model_data_reduced,
  aes(
    x = par_gene_count,
    y = pgls_resid
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.7
  ) +

  geom_point(
    aes(
      color = Order,
      fill = Order,
      shape = T2T_status
    ),
    size = 4,
    stroke = 1.2
  ) +

  scale_color_manual(
    values = mammal_order_colors,
    limits = mammal_order_levels,
    drop = FALSE
  ) +

  scale_fill_manual(
    values = mammal_order_fills,
    limits = mammal_order_levels,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  geom_text_repel(
    aes(label = Species),
    size = 3,
    color = "black",
    show.legend = FALSE
  ) +

  scale_color_manual(
    values = mammal_order_colors,
    limits = mammal_order_levels,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  labs(
    title = paste0(
      "PGLS residuals: All mammals; lambda = ",
      round(pgls_raw$optpar, 2)
    ),
    x = "Number of genes in the PAR",
    y = "PGLS residual",
    color = "Order",
    shape = "Assembly"
  ) +

  theme_bw()



# ============================================================
# 16. Save plots
# ============================================================

# PNGs

ggsave(
  "LM.SHROOM2_Mammals.png",
  p_lm_reduced,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "PGLS.SHROOM2_Mammals.png",
  p_pgls_reduced,
  width = 8,
  height = 6,
  dpi = 300
)


# PDFs

ggsave(
  "LM.SHROOM2_Mammals.pdf",
  p_lm_reduced,
  width = 8,
  height = 6
)

ggsave(
  "PGLS.SHROOM2_Mammals.pdf",
  p_pgls_reduced,
  width = 8,
  height = 6
)


# ============================================================
# Write model summaries to text files
# ============================================================

capture.output(
  summary(lm_order),
  file = "ANOVA.AllMammals.summary.txt"
)

capture.output(
  summary(lm_order_reduced),
  file = "ANOVA.ReducedMammals.summary.txt"
)

capture.output(
  summary(lm_raw),
  file = "LM.AllMammals.summary.txt"
)

capture.output(
  summary(pgls_raw),
  file = "PGLS.AllMammals.summary.txt"
)

capture.output(
  summary(lm_raw_reduced),
  file = "LM.ReducedMammals.summary.txt"
)

capture.output(
  summary(pgls_raw_reduced),
  file = "PGLS.ReducedMammals.summary.txt"
)