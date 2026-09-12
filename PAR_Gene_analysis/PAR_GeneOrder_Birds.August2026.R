library(tidyverse)
library(ape)
library(ggtree)
library(patchwork)
library(scales)

# ============================================================
# Inputs
# ============================================================

# Gene to use as ordinal origin point
GENE <- "DYM"

INVERT_PAR_ORDER_SPECIES <- "Amazona_ochrocephala"

# How many ordinal positions to show on either side of GENE
ORDINAL_MIN <- -30
ORDINAL_MAX <- 35

gene_file <- "birds_PAR_genes.all.tsv"

tree_file <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/roadies_v1.1.16b.numbers.scientific.edited.nwk"

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
    par_size_Mb = abs(PARStop - PARStart)
  )

# Rename ostrich
species_par <- species_par %>%
  mutate(
    Species = recode(
      Species,
      "Struthio_camelus" = "Struthio_camelus_australis"
    )
  )

# ============================================================
# Normalize column names from birds_PAR_genes.all.tsv
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
  filter(Species %in% c("Taeniopygia_guttata")) %>%
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

# ============================================================
# Diagnostic tree with node numbers for choosing rotations
# ============================================================

tree_filtered <- ape::rotate(tree_filtered, node = 35)
tree_filtered <- ape::rotate(tree_filtered, node = 36)
tree_filtered <- ape::rotate(tree_filtered, node = 37)
tree_filtered <- ape::rotate(tree_filtered, node = 38)
tree_filtered <- ape::rotate(tree_filtered, node = 39)
tree_filtered <- ape::rotate(tree_filtered, node = 40)
tree_filtered <- ape::rotate(tree_filtered, node = 41)
tree_filtered <- ape::rotate(tree_filtered, node = 42)
tree_filtered <- ape::rotate(tree_filtered, node = 43)
tree_filtered <- ape::rotate(tree_filtered, node = 44)

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
  scale_y_continuous(limits = c(0, 10)) +
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
    raw_PAR_order = row_number()
  ) %>%
  ungroup()

# ------------------------------------------------------------
# Find the origin gene in each species
# If there are multiple hits for GENE in a species, use the first
# ordinal occurrence after chromosome/PAR orientation adjustment.
# ------------------------------------------------------------

gene_origin <- par_genes %>%
  filter(Gene == GENE) %>%
  group_by(Species) %>%
  slice_min(raw_PAR_order, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Species, origin_raw_PAR_order = raw_PAR_order)

missing_origin_species <- setdiff(
  as.character(species_order),
  as.character(gene_origin$Species)
)

if (length(missing_origin_species) > 0) {
  warning(
    paste0(
      "Origin gene '", GENE, "' was not found in these species and they will be omitted from the gene-order panel: ",
      paste(missing_origin_species, collapse = ", ")
    )
  )
}


par_genes <- par_genes %>%
  inner_join(gene_origin, by = "Species") %>%
  mutate(
    PAR_order = raw_PAR_order - origin_raw_PAR_order
  ) %>%
  group_by(Species) %>%
  mutate(
    par_min_order = min(PAR_order[In_PAR %in% c("Y", "Edge")], na.rm = TRUE),
    par_max_order = max(PAR_order[In_PAR %in% c("Y", "Edge")], na.rm = TRUE),
    PAR_order = if_else(
      Species == INVERT_PAR_ORDER_SPECIES &
        In_PAR %in% c("Y", "Edge"),
      par_min_order + par_max_order - PAR_order,
      PAR_order
    )
  ) %>%
  ungroup() %>%
  select(-par_min_order, -par_max_order) %>%
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
  filter(Species == "Taeniopygia_guttata")

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

finch_gene_par_edge_counts <- ref_species_labels %>%
  distinct(Gene, Gene_label, PAR_order) %>%
  left_join(gene_par_edge_counts, by = "Gene") %>%
  mutate(
    n_genomes_PAR_or_Edge = replace_na(n_genomes_PAR_or_Edge, 0L)
  )

p_gene_par_edge_count <- ggplot(
  finch_gene_par_edge_counts,
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
    expand = expansion(mult = c(0, 0.03))
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
  filter(Gene %in% c("DYM")) %>%
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
    x = paste0("Ordinal gene order relative to ", GENE),
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

p_par_size <- ggplot(par_size, aes(x = par_size_Mb, y = Species)) +
  geom_col() +
  scale_y_discrete(limits = species_order) +
  scale_x_continuous(
    labels = scales::label_number(scale = 1e-6),
    limits = c(0, 20e6),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
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
    title = "Phylogeny, PAR size, and PAR gene order",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold")
    )
  )

ggsave(
  filename = paste0("combined_phylogeny_upset_gene_order_PAR_size.", GENE, ".pdf"),
  plot = combined_plot,
  width = 24,
  height = 18,
  limitsize = FALSE
)















# ============================================================
# PAR size vs PAR gene count in birds
#
# Models:
#   1. lm -- all birds
#   2. PGLS -- all birds
#   3. lm -- ratites removed
#   4. PGLS -- ratites removed
#
# Response:  par_gene_count
# Predictor: par_size_Mb
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ape)
library(phylolm)
library(ggplot2)
library(ggrepel)


# ============================================================
# 1. Define species metadata
# ============================================================

T2T <- c(
  "Taeniopygia_guttata"
)


species_orders <- tibble::tribble(
  ~Species,                     ~Order,
  "Aegotheles_albertisi",       "Aegotheliformes",
  "Anas_platyrhynchos",         "Anseriformes",
  "Aythya_ferina",              "Anseriformes",
  "Aythya_marila",              "Anseriformes",
  "Calonectris_borealis",       "Procellariiformes",
  "Coloeus_monedula",           "Passeriformes",
  "Columba_livia",              "Columbiformes",
  "Cyanocitta_cristata",        "Passeriformes",
  "Cygnus_columbianus",         "Anseriformes",
  "Dixiphia_pipra",             "Passeriformes",
  "Heliangelus_exortis",        "Apodiformes",
  "Larus_argentatus",           "Charadriiformes",
  "Lathamus_discolor",          "Psittaciformes",
  "Mergus_octosetaceus",        "Anseriformes",
  "Morphnus_guianensis",        "Accipitriformes",
  "Numenius_arquata",           "Charadriiformes",
  "Opisthocomus_hoazin",        "Opisthocomiformes",
  "Passer_domesticus",          "Passeriformes",
  "Patagioenas_fasciata",       "Columbiformes",
  "Platalea_leucorodia",        "Pelecaniformes",
  "Phaethon_aethereus",         "Phaethontiformes",
  "Rissa_tridactyla",           "Charadriiformes",
  "Sarcoramphus_papa",          "Cathartiformes",
  "Strix_aluco",                "Strigiformes",
  "Struthio_camelus_australis",           "Struthioniformes",
  "Taeniopygia_guttata",        "Passeriformes",
  "Zosterops_lateralis",        "Passeriformes",
  "Poecile_atricapillus",       "Passeriformes",
  "Amazona_ochrocephala",       "Psittaciformes",
  "Anser_anser",                "Anseriformes",
  "Eudromia_elegans",           "Tinamiformes",
  "Dromaius_novaehollandiae",   "Casuariiformes"
)

# ============================================================
# Common plotting aesthetics
# ============================================================

# ============================================================
# Point shapes
# ============================================================

# Fillable shapes:
# circle   = Non-T2T
# triangle = T2T
point_shapes <- c(
  "Non-T2T" = 21,
  "T2T"     = 24
)


# ============================================================
# 1. Build ONE all-bird dataset
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
    par_size_Mb = (PARStop - PARStart) / 1e6
  ) %>%
  left_join(par_gene_counts, by = "Species") %>%
  mutate(
    par_gene_count = replace_na(par_gene_count, 0L)
  ) %>%
  filter(
    is.finite(par_size_Mb),
    par_size_Mb > 0,
    is.finite(par_gene_count),
    par_gene_count > 0
  ) %>%
  distinct(Species, .keep_all = TRUE)


# ============================================================
# 2. Match all-bird dataset to phylogeny
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
# 3. ALL BIRDS -- lm
# ============================================================

lm_all <- lm(
  par_size_Mb ~ par_gene_count,
  data = model_data
)

summary(lm_all)


# ============================================================
# 4. ALL BIRDS -- PGLS
# ============================================================

pgls_all <- phylolm(
  par_size_Mb ~ par_gene_count,
  data = model_data,
  phy = model_tree,
  model = "lambda"
)

summary(pgls_all)

# Estimated Pagel's lambda
pgls_all$optpar


# ============================================================
# 5. Plot -- lm, all birds
# ============================================================
# ============================================================
# Bird Order plotting aesthetics
# ============================================================

# Species with no PAR
no_par_species <- c(
  "Colius_striatus",
  "Falco_naumanni"
)


# ============================================================
# 1. Add taxonomic Order FIRST
# ============================================================

model_data <- model_data %>%
  select(-any_of("Order")) %>%
  left_join(
    species_orders,
    by = "Species"
  )

# ============================================================
# 2. Remove species with no PAR
# ============================================================

model_data <- model_data %>%
  filter(!Species %in% no_par_species)

model_data <- model_data %>%
  mutate(
    T2T_status = if_else(
      Species %in% T2T,
      "T2T",
      "Non-T2T"
    )
  )

# ============================================================
# 3. Check that every remaining species has an Order
# ============================================================

model_data %>%
  filter(is.na(Order)) %>%
  select(Species)

# ============================================================
# 4. Define Order levels from the remaining ALL-BIRD dataset
# ============================================================

bird_order_levels <- model_data %>%
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
      levels = bird_order_levels
    )
  )


# ============================================================
# 6. Color-blind-friendly Order colors
# ============================================================

bird_order_colors <- setNames(
  viridisLite::viridis(
    length(bird_order_levels),
    option = "D",
    begin = 0.05,
    end = 0.95
  ),
  bird_order_levels
)


# ============================================================
# 7. Alternate filled/open circles
# ============================================================

bird_order_shapes <- setNames(
  rep(
    c(16, 1),
    length.out = length(bird_order_levels)
  ),
  bird_order_levels
)


bird_order_fills <- bird_order_colors

# Every second Order is open
bird_order_fills[seq(2, length(bird_order_fills), by = 2)] <- "white"

# ============================================================
# 8. Sanity checks
# ============================================================

table(model_data$Order, useNA = "ifany")


# ============================================================
# 5. Plot -- LM, all birds
# ============================================================

p_lm_all <- ggplot(
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
    values = bird_order_colors,
    limits = bird_order_levels,
    drop = FALSE
  ) +

  scale_fill_manual(
    values = bird_order_fills,
    limits = bird_order_levels,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    linewidth = 0.8
  ) +

  geom_text_repel(
    aes(label = Species),
    size = 3,
    color = "black",
    show.legend = FALSE
  ) +

  scale_color_manual(
    values = bird_order_colors,
    limits = bird_order_levels,
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
        shape = bird_order_shapes
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
# Plot -- PGLS residuals, all birds
# ============================================================

# Add all-bird PGLS residuals
model_data$pgls_resid <- residuals(
  pgls_all
)

p_pgls_resid_all <- ggplot(
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
    values = bird_order_colors,
    limits = bird_order_levels,
    drop = FALSE
  ) +

  scale_fill_manual(
    values = bird_order_fills,
    limits = bird_order_levels,
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
    values = bird_order_colors,
    limits = bird_order_levels,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  labs(
    title = paste0(
      "PGLS residuals: All birds; lambda = ",
      round(pgls_all$optpar, 2)
    ),
    x = "Number of genes in the PAR",
    y = "PGLS residual",
    color = "Order",
    shape = "Assembly"
  ) +

  theme_bw()






# ============================================================
# 7. Remove ratites
#
# This reproduces your existing bird analysis, which used
# PAR size < 14 Mb to define the reduced dataset.
# ============================================================

model_data_reduced <- model_data %>%
  filter(par_size_Mb < 14)


model_data_reduced <- model_data_reduced %>%
  select(-any_of("Order")) %>%
  left_join(
    species_orders,
    by = "Species"
  )


model_data_reduced <- model_data_reduced %>%
  filter(!Species %in% no_par_species)


model_data_reduced <- model_data_reduced %>%
  mutate(
    T2T_status = if_else(
      Species %in% T2T,
      "T2T",
      "Non-T2T"
    )
  )


model_data_reduced %>%
  filter(is.na(Order)) %>%
  select(Species)


model_data_reduced <- model_data_reduced %>%
  mutate(
    Order = factor(
      Order,
      levels = bird_order_levels
    )
  )
table(model_data_reduced$Order, useNA = "ifany")

# ============================================================
# 8. Re-match reduced dataset to phylogeny
# ============================================================

species_for_model_reduced <- intersect(
  tree_filtered$tip.label,
  model_data_reduced$Species
)

model_tree_reduced <- keep.tip(
  tree_filtered,
  species_for_model_reduced
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
# 9. NO RATITES -- lm
# ============================================================

lm_no_ratites <- lm(
  par_size_Mb ~ par_gene_count,
  data = model_data_reduced
)

summary(lm_no_ratites)


# ============================================================
# 10. NO RATITES -- PGLS
# ============================================================

pgls_no_ratites <- phylolm(
  par_size_Mb ~ par_gene_count,
  data = model_data_reduced,
  phy = model_tree_reduced,
  model = "lambda"
)

summary(pgls_no_ratites)

# Estimated Pagel's lambda
pgls_no_ratites$optpar


# ============================================================
# Add PGLS residuals to model_data_reduced
# ============================================================

model_data_reduced$pgls_resid <- residuals(
  pgls_no_ratites
)


# ============================================================
# 11. Plot -- LM, ratites removed
# ============================================================

p_lm_no_ratites <- ggplot(
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
    values = bird_order_colors,
    limits = bird_order_levels,
    drop = FALSE
  ) +

  scale_fill_manual(
    values = bird_order_fills,
    limits = bird_order_levels,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    linewidth = 0.8
  ) +

  geom_text_repel(
    aes(label = Species),
    size = 3,
    color = "black",
    show.legend = FALSE
  ) +

  scale_color_manual(
    values = bird_order_colors,
    limits = bird_order_levels,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  labs(
    title = "LM: Ratites removed",
    x = "Number of genes in the PAR",
    y = "PAR size (Mb)",
    color = "Order",
    shape = "Assembly"
  ) +

  theme_bw()



# ============================================================
# 12. Plot -- PGLS residuals, ratites removed
# ============================================================

p_pgls_resid_no_ratites <- ggplot(
  model_data_reduced,
  aes(
    x = par_gene_count,
    y = pgls_resid
  )
) +

  # Expected residual = 0
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
    values = bird_order_colors,
    limits = bird_order_levels,
    drop = FALSE
  ) +

  scale_fill_manual(
    values = bird_order_fills,
    limits = bird_order_levels,
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
    values = bird_order_colors,
    limits = bird_order_levels,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = point_shapes
  ) +

  labs(
    title = paste0(
      "PGLS residuals: Ratites removed; lambda = ",
      round(pgls_no_ratites$optpar, 2)
    ),
    x = "Number of genes in the PAR",
    y = "PGLS residual",
    color = "Order",
    shape = "Assembly"
  ) +

  theme_bw()


# ============================================================
# 13. Save plots
# ============================================================

ggsave(
  "lm.AllBirds.pdf",
  p_lm_all,
  width = 8,
  height = 6
)

ggsave(
  "PGLS.residuals.AllBirds.pdf",
  p_pgls_resid_all,
  width = 8,
  height = 6
)

ggsave(
  "lm.noRatites.pdf",
  p_lm_no_ratites,
  width = 8,
  height = 6
)

ggsave(
  "PGLS.residuals.noRatites.pdf",
  p_pgls_resid_no_ratites,
  width = 8,
  height = 6
)



# ============================================================
# Write model summaries to text files
# ============================================================

capture.output(
  summary(lm_all),
  file = "LM.AllBirds.summary.txt"
)

capture.output(
  summary(pgls_all),
  file = "PGLS.AllBirds.summary.txt"
)

capture.output(
  summary(lm_no_ratites),
  file = "LM.ReducedBirds.summary.txt"
)

capture.output(
  summary(pgls_no_ratites),
  file = "PGLS.ReducedBirds.summary.txt"
)