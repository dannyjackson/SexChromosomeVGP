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

acc_col        <- "# accession"
name_col       <- "EnglishName"
lineage_col    <- "VGPLineage"
order_col      <- "Orders_English"
sciname_col    <- "ScientificName"

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
## Helper: read PAR segments from BED and get lengths
## If BED file is size 0 → return empty numeric vector
## --------------------------
get_par_segment_lengths <- function(sciname) {
  # BED filename prefix: e.g. "Nyctalus_leisleri_"
  prefix <- paste0(gsub(" ", "_", sciname), "_")
  f <- bed_files[startsWith(bed_files, prefix)]
  
  if (length(f) == 0)
    return(numeric(0))   # No BED file found → no segments
  
  if (length(f) > 1) {
    warning("Multiple BED files for ", sciname, "; using first: ", f[1])
    f <- f[1]
  }
  
  bed_path <- file.path(bed_dir, f)
  
  # Skip if BED file is size 0
  if (file.info(bed_path)$size == 0) {
    message("Skipping zero-size BED for: ", sciname)
    return(numeric(0))
  }
  
  # Snakemake-style BED, no header
  bed <- fread(bed_path, header = FALSE)
  seg_lengths <- bed[[3]] - bed[[2]]  # V3 - V2
  seg_lengths
}

## --------------------------
## Build long df of PAR segments per tip
## Carry Orders_English through explicitly
## --------------------------
par_list <- lapply(seq_len(nrow(tip_meta)), function(i) {
  sci   <- tip_meta$ScientificName[i]
  name  <- tip_meta$EnglishName[i]
  ord   <- tip_meta$Orders_English[i]
  
  segs <- get_par_segment_lengths(sci)
  if (length(segs) == 0) {
    return(NULL)   # keep tip, just no PAR rows
  }
  
  data.frame(
    EnglishName    = name,
    Orders_English = ord,
    segment_id     = paste0("seg", seq_along(segs)),
    len_bp         = segs,
    stringsAsFactors = FALSE
  )
})

par_df <- bind_rows(par_list)

# For plotting in Mb (optional)
if (nrow(par_df) > 0) {
  par_df <- par_df %>%
    mutate(len_Mb = len_bp / 1e6)
}

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

## ------------------------------------
## 2–3. PAR totals & bar plot (unchanged)
## ------------------------------------
par_totals <- par_df %>%
  group_by(EnglishName, Orders_English) %>%
  summarise(total_Mb = sum(len_Mb), .groups = "drop")

tip_order <- p_tree$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)

par_totals <- par_totals %>%
  mutate(EnglishName = factor(EnglishName, levels = tip_order))

p_bar <- ggplot(par_totals,
                aes(y = EnglishName,
                    x = total_Mb,
                    fill = Orders_English)) +
  geom_col() +
  scale_fill_manual(values = order_pal, guide = "none") +
  theme_minimal(base_size = 9) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  labs(x = "PAR length (Mb)")

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

## ------------------------------------
## 5. Combine tree (no legend) + bar, then attach legend at far right
## ------------------------------------
combined <- plot_grid(
  p_tree_noleg,
  p_bar,
  ncol = 2,
  rel_widths = c(3, 1),
  align = "h"
)

final_plot <- plot_grid(
  combined,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.25)
)

## ------------------------------------
## 6. Save
## ------------------------------------
pdf_file <- paste0("subtree_", target_lineage,
                   "_tree_with_PARbars_byOrder.pdf")
ggsave(pdf_file, plot = final_plot, width = 8, height = 12, useDingbats = FALSE)

cat("Wrote tree + PAR bar plot + small legend to:", pdf_file, "\n")
