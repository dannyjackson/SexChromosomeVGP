# ============================================================
# THREE-PANEL ORDER FIGURE
# Phylogeny | PAR gene count | PAR size
#
# ANOVA only for Orders with >5 taxa
# All Orders shown
# ============================================================

library(dplyr)
library(ggplot2)
library(emmeans)
library(multcompView)
library(ape)
library(ggtree)
library(patchwork)


# ============================================================
# 1. SET UP PHYLOGENY ONCE
# ============================================================

# Use all species represented in model_data_reduced
tree_data_all <- model_data_reduced %>%
  dplyr::filter(
    !is.na(Order),
    !is.na(Species)
  ) %>%
  dplyr::distinct(Species, Order)

species_keep <- intersect(
  tree_filtered$tip.label,
  tree_data_all$Species
)

tree_reduced <- ape::keep.tip(
  tree_filtered,
  species_keep
)

tree_data_all <- tree_data_all %>%
  dplyr::filter(
    Species %in% tree_reduced$tip.label
  )


# ------------------------------------------------------------
# Get phylogenetic ordering of Orders
# ------------------------------------------------------------

phylo_tip_order <- tree_reduced$tip.label

tree_data_all <- tree_data_all %>%
  dplyr::mutate(
    tip_position = match(
      Species,
      phylo_tip_order
    )
  )

order_phylo_order <- tree_data_all %>%
  dplyr::group_by(Order) %>%
  dplyr::summarise(
    phylo_position = mean(tip_position),
    .groups = "drop"
  ) %>%
  dplyr::arrange(phylo_position) %>%
  dplyr::pull(Order)


# ============================================================
# 2. CREATE COLLAPSED PHYLOGENY
# ============================================================

tree_plot <- ggtree(
  tree_reduced,
  ladderize = TRUE
)

collapse_ggtree <- getS3method(
  "collapse",
  "ggtree",
  envir = asNamespace("ggtree")
)

# Store Order names + nodes for labeling
order_nodes <- data.frame(
  Order = character(),
  node = integer()
)

for (ord in order_phylo_order) {

  spp <- tree_data_all %>%
    dplyr::filter(Order == ord) %>%
    dplyr::pull(Species)

  # More than one species: collapse Order
  if (length(spp) > 1) {

    node <- ape::getMRCA(
      tree_reduced,
      spp
    )

    if (!is.null(node)) {

      tree_plot <- collapse_ggtree(
        tree_plot,
        node = node,
        mode = "none"
      )

      order_nodes <- rbind(
        order_nodes,
        data.frame(
          Order = ord,
          node = node
        )
      )
    }

  # One species: label that species tip with its Order
  } else if (length(spp) == 1) {

    node <- which(
      tree_reduced$tip.label == spp
    )

    order_nodes <- rbind(
      order_nodes,
      data.frame(
        Order = ord,
        node = node
      )
    )
  }
}


# Get plotted x/y coordinates of the collapsed nodes
order_labels <- tree_plot$data %>%
  dplyr::filter(node %in% order_nodes$node) %>%
  dplyr::left_join(
    order_nodes,
    by = "node"
  )


# Add simple text labels -- nothing else
tree_plot <- tree_plot +
  geom_text(
    data = order_labels,
    aes(
      x = x,
      y = y,
      label = Order
    ),
    inherit.aes = FALSE,
    hjust = -0.15,
    size = 3
  ) +
  theme_tree2() +
  theme(
    plot.margin = margin(5.5, 15, 5.5, 5.5)
  )


# ============================================================
# 3. PAR SIZE
# ============================================================

plot_data <- model_data_reduced %>%
  dplyr::filter(
    !is.na(Order),
    !is.na(Species),
    !is.na(par_size_Mb)
  )


# ------------------------------------------------------------
# Count taxa and identify ANOVA Orders
# ------------------------------------------------------------

order_counts <- plot_data %>%
  dplyr::distinct(Order, Species) %>%
  dplyr::count(
    Order,
    name = "n_taxa"
  )

anova_orders <- order_counts %>%
  dplyr::filter(n_taxa > 4) %>%
  dplyr::pull(Order)


# ------------------------------------------------------------
# ANOVA
# ------------------------------------------------------------

model_data_order_anova <- plot_data %>%
  dplyr::filter(
    Order %in% anova_orders
  ) %>%
  droplevels()

lm_order_size <- lm(
  par_size_Mb ~ Order,
  data = model_data_order_anova
)

anova(lm_order_size)


# ------------------------------------------------------------
# Tukey pairwise comparisons
# ------------------------------------------------------------

emm_order <- emmeans(
  lm_order_size,
  ~ Order
)

pairwise_order <- pairs(
  emm_order,
  adjust = "tukey"
)

pairwise_df <- as.data.frame(
  summary(pairwise_order)
)


# ------------------------------------------------------------
# Compact letters from Tukey-adjusted p-values
# ------------------------------------------------------------

pvals <- pairwise_df$p.value

names(pvals) <- gsub(
  " - ",
  "-",
  pairwise_df$contrast,
  fixed = TRUE
)

letters_tukey <- multcompView::multcompLetters(
  pvals,
  threshold = 0.05
)$Letters

order_letters <- data.frame(
  Order = names(letters_tukey),
  .group = unname(letters_tukey)
)


# ------------------------------------------------------------
# Summarize ALL Orders
# ------------------------------------------------------------

order_summary <- plot_data %>%
  dplyr::group_by(Order) %>%
  dplyr::summarise(
    n_taxa = dplyr::n_distinct(Species),
    mean_size = mean(par_size_Mb),
    min_size = min(par_size_Mb),
    max_size = max(par_size_Mb),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    ANOVA_group = n_taxa > 4
  ) %>%
  dplyr::left_join(
    order_letters,
    by = "Order"
  )


# ------------------------------------------------------------
# Use SAME phylogenetic Order ordering
# ------------------------------------------------------------

order_summary <- order_summary %>%
  dplyr::filter(
    Order %in% order_phylo_order
  ) %>%
  dplyr::mutate(
    Order = factor(
      Order,
      levels = rev(order_phylo_order)
    )
  )


# ------------------------------------------------------------
# PAR size plot
# ------------------------------------------------------------

offset <- 0.03 * diff(
  range(
    plot_data$par_size_Mb,
    na.rm = TRUE
  )
)

par_size_plot <- ggplot(
  order_summary,
  aes(
    y = Order,
    x = mean_size
  )
) +

  geom_errorbar(
    aes(
      xmin = min_size,
      xmax = max_size,
      color = ANOVA_group
    ),
    orientation = "y",
    width = 0.20,
    linewidth = 0.7
  ) +

  geom_point(
    aes(
      color = ANOVA_group
    ),
    size = 3
  ) +

  geom_text(
    data = order_summary %>%
      dplyr::filter(!is.na(.group)),
    aes(
      x = max_size + offset,
      label = .group
    ),
    color = "black",
    hjust = 0,
    size = 4,
    fontface = "bold"
  ) +

  scale_color_manual(
    values = c(
      `TRUE` = "black",
      `FALSE` = "grey65"
    ),
    breaks = c(TRUE, FALSE),
    labels = c(
      "ANOVA (n >= 5)",
      "Not tested (n <= 5)"
    )
  ) +

  scale_x_continuous(
    expand = expansion(
      mult = c(0.03, 0.15)
    )
  ) +

  labs(
    x = "PAR size (Mb)",
    y = NULL,
    color = NULL
  ) +

  theme_classic(base_size = 12) +

  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )


# ============================================================
# 4. PAR GENE COUNT
# ============================================================

plot_data <- model_data_reduced %>%
  dplyr::filter(
    !is.na(Order),
    !is.na(Species),
    !is.na(par_gene_count)
  )


# ------------------------------------------------------------
# Count taxa and identify ANOVA Orders
# ------------------------------------------------------------

order_counts <- plot_data %>%
  dplyr::distinct(Order, Species) %>%
  dplyr::count(
    Order,
    name = "n_taxa"
  )

anova_orders <- order_counts %>%
  dplyr::filter(n_taxa > 4) %>%
  dplyr::pull(Order)


# ------------------------------------------------------------
# ANOVA
# ------------------------------------------------------------

model_data_order_anova <- plot_data %>%
  dplyr::filter(
    Order %in% anova_orders
  ) %>%
  droplevels()

lm_order_gene_count <- lm(
  par_gene_count ~ Order,
  data = model_data_order_anova
)

anova(lm_order_gene_count)


# ------------------------------------------------------------
# Tukey pairwise comparisons
# ------------------------------------------------------------

emm_order <- emmeans(
  lm_order_gene_count,
  ~ Order
)

pairwise_order <- pairs(
  emm_order,
  adjust = "tukey"
)

pairwise_df <- as.data.frame(
  summary(pairwise_order)
)


# ------------------------------------------------------------
# Compact letters
# ------------------------------------------------------------

pvals <- pairwise_df$p.value

names(pvals) <- gsub(
  " - ",
  "-",
  pairwise_df$contrast,
  fixed = TRUE
)

letters_tukey <- multcompView::multcompLetters(
  pvals,
  threshold = 0.05
)$Letters

order_letters <- data.frame(
  Order = names(letters_tukey),
  .group = unname(letters_tukey)
)


# ------------------------------------------------------------
# Summarize ALL Orders
# ------------------------------------------------------------

order_summary <- plot_data %>%
  dplyr::group_by(Order) %>%
  dplyr::summarise(
    n_taxa = dplyr::n_distinct(Species),
    mean_size = mean(par_gene_count),
    min_size = min(par_gene_count),
    max_size = max(par_gene_count),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    ANOVA_group = n_taxa > 4
  ) %>%
  dplyr::left_join(
    order_letters,
    by = "Order"
  )


# ------------------------------------------------------------
# SAME phylogenetic Order ordering as size panel
# ------------------------------------------------------------

order_summary <- order_summary %>%
  dplyr::filter(
    Order %in% order_phylo_order
  ) %>%
  dplyr::mutate(
    Order = factor(
      Order,
      levels = rev(order_phylo_order)
    )
  )


# ------------------------------------------------------------
# Gene count plot
# ------------------------------------------------------------

offset <- 0.03 * diff(
  range(
    plot_data$par_gene_count,
    na.rm = TRUE
  )
)

gene_count_plot <- ggplot(
  order_summary,
  aes(
    y = Order,
    x = mean_size
  )
) +

  geom_errorbar(
    aes(
      xmin = min_size,
      xmax = max_size,
      color = ANOVA_group
    ),
    orientation = "y",
    width = 0.20,
    linewidth = 0.7
  ) +

  geom_point(
    aes(
      color = ANOVA_group
    ),
    size = 3
  ) +

  geom_text(
    data = order_summary %>%
      dplyr::filter(!is.na(.group)),
    aes(
      x = max_size + offset,
      label = .group
    ),
    color = "black",
    hjust = 0,
    size = 4,
    fontface = "bold"
  ) +

  scale_color_manual(
    values = c(
      `TRUE` = "black",
      `FALSE` = "grey65"
    ),
    breaks = c(TRUE, FALSE),
    labels = c(
      "ANOVA (n >= 5)",
      "Not tested (n <= 5)"
    )
  ) +

  scale_x_continuous(
    expand = expansion(
      mult = c(0.03, 0.15)
    )
  ) +

  labs(
    x = "PAR gene count",
    y = NULL,
    color = NULL
  ) +

  theme_classic(base_size = 12) +

  theme(
    legend.position = "top",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )


# ============================================================
# 5. THREE-PANEL FIGURE
# ============================================================

mammal_order_plot <-
  tree_plot +
  gene_count_plot +
  par_size_plot +
  patchwork::plot_layout(
    widths = c(1.2, 2, 2)
  )

mammal_order_plot


# ============================================================
# 6. SAVE
# ============================================================

ggsave(
  "Orders.PAR_GeneCount_PARsize_Phylogeny.pdf",
  mammal_order_plot,
  width = 5,
  height = 8
)