#!/usr/bin/env Rscript

library(ape)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(ggtree)   # Bioconductor
library(dplyr)

## --------------------------
## Inputs
## --------------------------
tree_file <- Sys.getenv("TREE")
ann_file  <- Sys.getenv("ANN")
bed_dir   <- Sys.getenv("BEDDIR")
sex_chr_len   <- Sys.getenv("SEXCHRLEN")
args <- commandArgs(trailingOnly = TRUE)

if (tree_file == "" || ann_file == "" || bed_dir == "") {
  stop("Environment variables TREE, ANN, and BEDDIR must all be set.")
}

target_lineage <- args[1] 

## --------------------------
## Read data
## --------------------------
tree <- read.tree(tree_file)
ann  <- fread(ann_file)

# Read file from environment variable
sex_chr_len_file <- Sys.getenv("SEXCHRLEN")
sex_chr <- fread(sex_chr_len_file)

# Get names of the accession columns
col_sex <- names(sex_chr)[1]
col_ann <- names(ann)[1]

# Normalize accession IDs by removing leading GCA_ or GCF_
sex_chr[, acc_norm := sub("^GC[AF]_", "", get(col_sex))]
ann[,     acc_norm := sub("^GC[AF]_", "", get(col_ann))]

# Subset rows where chromosome = X or Z
# Assuming the chromosome identity is stored in column 1 of sex_chr
sex_chr_sub <- sex_chr[
  get("Chromosome name") %chin% c("X", "Z")
]

# Join to ann on the first column of both tables
joined_df <- ann %>%
  inner_join(
    sex_chr_sub,
    by = "acc_norm"
  )

ann <- joined_df

acc_col        <- "# accession"
name_col       <- "EnglishName"
lineage_col    <- "Lineage"
order_col      <- "Orders_English"
sciname_col    <- "ScientificName"
chrlen_col     <- "Seq length"
for (col in c(acc_col, name_col, lineage_col, order_col, sciname_col)) {
  if (!(col %in% names(ann))) {
    stop(paste("Column not found in annotation file:", col))
  }
}

## --------------------------
## Determine which ScientificNames have BED files
## --------------------------
bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = FALSE)
if (length(bed_files) == 0) {
  stop("No .bed files found in BEDDIR.")
}

# First two underscore-separated tokens => scientific name part
bed_prefix   <- sub("^([^_]+_[^_]+).*", "\\1", bed_files)
bed_sciname  <- gsub("_", " ", bed_prefix)  # "Acanthisitta chloris"
scientific_with_bed <- unique(bed_sciname)

## --------------------------
## Subset annotations: lineage AND has BED
## (if you want all mammals regardless of BED, drop the scientific_with_bed filter)
## --------------------------
ann_sub <- ann[
  get(lineage_col) == target_lineage &
    get(sciname_col) %in% scientific_with_bed
]

if (nrow(ann_sub) == 0) {
  stop("No individuals in the target lineage have matching BED files.")
}

## --------------------------
## Match to tree tips via accession, then prune
## --------------------------
tips_to_keep <- intersect(tree$tip.label, ann_sub[[acc_col]])
if (length(tips_to_keep) == 0) {
  stop("No matching accessions in the tree for this lineage with BEDs.")
}

subtree <- keep.tip(tree, tips_to_keep)

## --------------------------
## Build per-tip metadata in tree order
## --------------------------
ann_sub_pruned <- ann_sub[get(acc_col) %in% subtree$tip.label]

tip_meta <- data.frame(
  accession       = subtree$tip.label,
  EnglishName     = ann_sub_pruned[match(subtree$tip.label, ann_sub_pruned[[acc_col]]),
                                   get(name_col)],
  Orders_English  = ann_sub_pruned[match(subtree$tip.label, ann_sub_pruned[[acc_col]]),
                                   get(order_col)],
  ScientificName  = ann_sub_pruned[match(subtree$tip.label, ann_sub_pruned[[acc_col]]),
                                   get(sciname_col)],
  Chr_length      = ann_sub_pruned[match(subtree$tip.label, ann_sub_pruned[[acc_col]]),
                                   get(chrlen_col)],
  stringsAsFactors = FALSE
)

## --------------------------
## Tip colors by Orders_English
## --------------------------
orders <- sort(unique(tip_meta$Orders_English))
n_orders <- length(orders)

base_cols <- brewer.pal(min(max(3, n_orders), 12), "Set3")
if (n_orders > length(base_cols)) {
  cols <- colorRampPalette(base_cols)(n_orders)
} else {
  cols <- base_cols[seq_len(n_orders)]
}

order_pal  <- setNames(cols, orders)

## --------------------------
## Tree object with EnglishName labels
## --------------------------
subtree_plot <- subtree
subtree_plot$tip.label <- tip_meta$EnglishName

## --------------------------
## Helper: read PAR segments from ALL BEDs for a sciname
## Return one row per segment, with a bam_id
## --------------------------
get_par_segments_df <- function(sciname) {
  prefix <- paste0(gsub(" ", "_", sciname), "_")
  f <- bed_files[startsWith(bed_files, prefix)]
  
  if (length(f) == 0) {
    return(NULL)  # no BED for this species
  }
  
  seg_list <- lapply(seq_along(f), function(i) {
    bed_fname <- f[i]
    bed_path  <- file.path(bed_dir, bed_fname)
    
    # Zero-size BED → treat as a single 0-length "segment"
    if (file.info(bed_path)$size == 0) {
      return(data.frame(
        bam_id = paste0("bam", i),
        len_bp = 0,
        stringsAsFactors = FALSE
      ))
    }
    
    bed <- fread(bed_path, header = FALSE)
    seg_lengths <- bed[[3]] - bed[[2]]
    
    data.frame(
      bam_id = paste0("bam", i),
      len_bp = seg_lengths,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, seg_list)
}

## --------------------------
## Build long df of PAR segments per tip (with bam_id)
## --------------------------
par_list <- lapply(seq_len(nrow(tip_meta)), function(i) {
  sci   <- tip_meta$ScientificName[i]
  name  <- tip_meta$EnglishName[i]
  ord   <- tip_meta$Orders_English[i]
  chrlen   <- tip_meta$Chr_length[i]

  seg_df <- get_par_segments_df(sci)
  if (is.null(seg_df) || nrow(seg_df) == 0) {
    return(NULL)
  }
  
  cbind(
    EnglishName    = name,
    Orders_English = ord,
    seg_df,
    Chr_length = chrlen,
    stringsAsFactors = FALSE
  )
})

par_df <- bind_rows(par_list)

if (nrow(par_df) > 0) {
  par_df <- par_df %>%
    mutate(len_Mb = len_bp / 1e6)
}


par_df$ratio <- par_df$len_bp/par_df$Chr_length

library(dplyr)
library(ggtree)
library(ggplot2)
library(cowplot)

## ------------------------------------
## 1. Base ggtree with colored tips + labels
## ------------------------------------
p_tree <- ggtree(subtree_plot, size = 0.4)

p_tree$data <- p_tree$data %>%
  left_join(
    tip_meta %>%
      transmute(label = EnglishName,
                Orders_English = Orders_English),
    by = "label"
  )

p_tree <- p_tree +
  geom_tippoint(aes(color = Orders_English), size = 2, na.rm = TRUE) +
  geom_tiplab(size = 2.5, hjust = 0, offset = 0.002) +
  scale_color_manual(values = order_pal, name = "Order") +
  theme_tree2() +
  theme(
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )
max_x <- max(p_tree$data$x, na.rm = TRUE)
p_tree <- p_tree + xlim(0, max_x + 0.1)

## ------------------------------------
## 2 PAR totals & bar plot (unchanged)
## ------------------------------------
par_totals <- par_df %>%
  group_by(EnglishName, Orders_English, bam_id) %>%
  summarise(total_Mb = sum(len_Mb), .groups = "drop")

tip_order <- p_tree$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)

par_totals <- par_totals %>%
  mutate(EnglishName = factor(EnglishName, levels = tip_order))


## Small x offset for the gray X just right of the y-axis
max_x <- max(par_totals$total_Mb, na.rm = TRUE)
x_nudge <- if (is.finite(max_x) && max_x > 0) max_x * 0.02 else 0.1

library(ggpattern)

p_bar <- ggplot(par_totals,
                aes(y = EnglishName,
                    x = total_Mb,
                    fill = Orders_English,
                    pattern = bam_id)) +
  ggpattern::geom_col_pattern(
    position = "stack",
    colour  = NA,             # no border
    pattern_fill   = "white",      # let base fill (Order color) show through
    pattern_colour = "white", # color of pattern lines
    pattern_size = 0.05,
    pattern_spacing = 0.05,
    pattern_key_scale_factor = 0.1
  ) +
  scale_fill_manual(values = order_pal, guide = "none") +
  # use two valid pattern names; 'weave' may not exist in your ggpattern version
  scale_pattern_manual(values = c(bam1 = "none", bam2 = "stripe")) +
  theme_minimal(base_size = 9) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  labs(x = "PAR length (Mb)") +
  geom_text(
    data = subset(par_totals, total_Mb == 0 & !is.na(total_Mb)),
    aes(x = 0, y = EnglishName),
    label = "X",
    color = "grey40",
    size  = 2.8,
    nudge_x = x_nudge,
    inherit.aes = FALSE
  )


## ------------------------------------
## 3 PAR totals & bar plot (ratio)
## ------------------------------------
par_ratio <- par_df %>%
  group_by(EnglishName, Orders_English, bam_id) %>%
  summarise(ratio = sum(ratio), .groups = "drop")

tip_order <- p_tree$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)

par_ratio <- par_ratio %>%
  mutate(EnglishName = factor(EnglishName, levels = tip_order))


## Small x offset for the gray X just right of the y-axis
max_x <- max(par_ratio$ratio, na.rm = TRUE)
x_nudge <- if (is.finite(max_x) && max_x > 0) max_x * 0.02 else 0.1


p_ratio_bar <- ggplot(par_ratio,
                aes(y = EnglishName,
                    x = ratio,
                    fill = Orders_English,
                    pattern = bam_id)) +
  ggpattern::geom_col_pattern(
    position = "stack",
    colour  = NA,             # no border
    pattern_fill   = "white",      # let base fill (Order color) show through
    pattern_colour = "white", # color of pattern lines
    pattern_size = 0.05,
    pattern_spacing = 0.05,
    pattern_key_scale_factor = 0.1
  ) +
  scale_fill_manual(values = order_pal, guide = "none") +
  # use two valid pattern names; 'weave' may not exist in your ggpattern version
  scale_pattern_manual(values = c(bam1 = "none", bam2 = "stripe")) +
  theme_minimal(base_size = 9) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  labs(x = "PAR ratio") +
  geom_text(
    data = subset(par_ratio, ratio == 0 & !is.na(ratio)),
    aes(x = 0, y = EnglishName),
    label = "X",
    color = "grey40",
    size  = 2.8,
    nudge_x = x_nudge,
    inherit.aes = FALSE
  )

legend_pattern <- cowplot::get_legend(
  p_ratio_bar +
    guides(pattern = guide_legend(
      title = "X1 or X2",
      ncol  = 1,
      override.aes = list(
        fill          = "grey70",      # neutral base so pattern visible
        pattern_color = "black",       # pattern lines color
        size          = 3
      )
    )) +
    theme(
      legend.key.size = unit(0.3, "cm"),
      legend.text     = element_text(size = 7),
      legend.title    = element_text(size = 8)
    )
)


## ------------------------------------
## 4. Extract smaller legend from a COPY of p_tree
## ------------------------------------
p_tree_legend <- p_tree  # keep legend on this one only

legend <- cowplot::get_legend(
  p_tree_legend +
    guides(color = guide_legend(title = "Order",
                                ncol = 1,
                                override.aes = list(size = 3))) +
    theme(
      legend.key.size = unit(0.3, "cm"),
      legend.text     = element_text(size = 7),
      legend.title    = element_text(size = 8)
    )
)

## Make a version of the tree WITHOUT legend for the combined plot
p_tree_noleg <- p_tree + theme(legend.position = "none")
p_bar <- p_bar + theme(legend.position = "none")
p_ratio_bar <- p_ratio_bar + theme(legend.position = "none")

## ------------------------------------
## 5. Combine tree (no legend) + bar, then attach legend at far right
## ------------------------------------
combined_plot <- plot_grid(
  p_tree_noleg,
  p_bar,
  p_ratio_bar,
  ncol = 3,
  rel_widths = c(3, 1, 1),
  align = "h"
)

combined_legend <- plot_grid(
  legend,
  legend_pattern,
  ncol = 1,
  align = "v"
)

final_plot <- plot_grid(
  combined_plot,
  combined_legend,
  ncol = 2,
  rel_widths = c(1, 0.25)
)

## ------------------------------------
## 6. Save
## ------------------------------------
pdf_file <- paste0("subtree_", target_lineage,
                   "_tree_with_PARbars_byOrder.pdf")
ggsave(pdf_file, plot = final_plot, width = 10, height = 12, useDingbats = FALSE)


png_file <- paste0("subtree_", target_lineage,
                   "_tree_with_PARbars_byOrder.png")

ggsave(
  filename = png_file,
  plot     = final_plot,
  width    = 10,
  height   = 12,
  dpi      = 600,
  bg       = "white"  
)


cat("Wrote tree + PAR bar plot + small legend to:", pdf_file, " and ", png_file, "\n")
