PAR inference with tree

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_inference/alignment/continuous_percentID


# ============================================================
# Continuous percentID: ALL files stacked, with chromosome boxes
# PLUS: tree panel on the left, ordering species by the tree
# ============================================================

library(data.table)
library(ggplot2)
library(patchwork)

# Tree packages
library(ape)

# ----------------------------
# Inputs
# ----------------------------
tree_file <- "/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/roadies_v1.1.16b.numbers.scientific.nwk"


thr <- 98.5
min_aln <- 50000

thr_tag <- gsub("\\.", "p", sprintf("%.1f", thr))  # "98p5"

mb_formatter <- scales::label_number(
  scale = 1e-6,
  accuracy = 1e-6
)

mb_endpoints_only <- function(lims) {
  mids <- mean(lims)
  list(
    breaks = c(lims[1], mids, lims[2]),
    labels = c(
      sprintf("%.1f", lims[1] / 1e6),
      "",
      sprintf("%.1f", lims[2] / 1e6)
    )
  )
}

# ----------------------------
# Paths
# ----------------------------
setwd("/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_inference/alignment/continuous_percentID")

dir_in <- "/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/minimap2/continuous_percentID"
files  <- list.files(dir_in, pattern = "\\.refqry\\.csv$", full.names = TRUE)


# ----------------------------
# Read + standardize one file
# ----------------------------
read_one <- function(f) {
  raw <- fread(f)

  # file label early so we can return a placeholder if needed
  file_label <- sub("\\.aln\\.refqry\\.csv$", "", basename(f))

  # standardize
  dt <- raw[, list(
    chrom_qry = as.character(chrom_qry),
    len_qry   = as.numeric(len_qry),
    qs        = as.integer(bp_start_qry),
    qe        = as.integer(bp_end_qry),

    chrom_ref = as.character(chrom_ref),
    len_ref   = as.numeric(len_ref),
    rs        = as.integer(bp_start_ref),
    re        = as.integer(bp_end_ref),

    pid       = as.numeric(percent_identity_qry)
  )]

  # basic validity
  dt <- dt[!is.na(pid) & qe > qs & re > rs]

  # length + pid filters
  dt[, aln_qry := qe - qs]
  dt[, aln_ref := re - rs]
  dt[, aln_len := pmin(aln_qry, aln_ref)]
  dt <- dt[aln_len >= min_aln & pid >= thr]

  # If nothing passes, return a placeholder row for plotting row/labels/tree
  if (nrow(dt) == 0) {
    return(data.table(
      chrom_qry = NA_character_,
      len_qry   = NA_real_,
      qs = 0L, qe = 0L,
      chrom_ref = NA_character_,
      len_ref   = NA_real_,
      rs = 0L, re = 0L,
      pid = NA_real_,
      aln_qry = 0L, aln_ref = 0L, aln_len = 0L,
      file_label = file_label,
      chrom_pair = "NO_ALIGNMENTS"
    ))
  }

  dt[, file_label := file_label]
  dt[, chrom_pair := paste(chrom_qry, "→", chrom_ref)]
  dt[]
}


# ----------------------------
# Read all files
# ----------------------------
dt_list <- lapply(files, read_one)
dt_list <- Filter(Negate(is.null), dt_list)

if (length(dt_list) == 0) stop("No non-empty *.refqry.csv files found after filtering.")

dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# ============================================================
# TREE: read + define species order
# ============================================================
tr <- read.tree(tree_file)
tr <- ladderize(tr)

# Parse species ID from file_label:
# default: first two underscore-separated tokens => Genus_species
tipset <- unique(tr$tip.label)

get_species <- function(x, tipset) {
  parts <- strsplit(x, "_", fixed = TRUE)[[1]]

  # cumulative prefixes
  prefs <- vapply(seq_along(parts), function(k) paste(parts[1:k], collapse = "_"), character(1))

  # keep only those that exist in the tree tips
  hits <- prefs[prefs %in% tipset]

  if (length(hits) > 0) {
    hits[length(hits)]              # longest match
  } else if (length(parts) >= 2) {
    paste(parts[1], parts[2], sep = "_")  # fallback
  } else {
    x
  }
}

dt[, species := vapply(as.character(file_label), get_species, character(1), tipset = tipset)]


# --- after: dt[, species := ...]
library(ape)

tr <- read.tree(tree_file)
tr <- ladderize(tr)

# species present in the dataset
sp_data <- unique(dt$species)

# intersect with tree tips
keep_tips <- intersect(tr$tip.label, sp_data)

if (length(keep_tips) < 2) {
  stop("After pruning, tree has <2 tips. Check species parsing (get_species()) vs tip labels.")
}

# prune tree to dataset species
tr <- keep.tip(tr, keep_tips)
tr <- ladderize(tr)

# now use pruned tree tip order
tree_species <- tr$tip.label

# any species in dt but not in tree should be dropped or sent to end; since we pruned, drop them
dt <- dt[species %in% tree_species]

# factor species in tree order
dt[, species := factor(species, levels = tree_species)]

# Keep only species present in tree, but don't crash if some missing
tree_species <- tr$tip.label
dt[, in_tree := species %in% tree_species]

if (!any(dt$in_tree)) {
  stop("None of the species parsed from file_label match tree tip labels. Check get_species().")
}

# Order species by tree tip order; any non-tree species go to the end
dt[, species := factor(species, levels = tree_species)]

# ============================================================
# Build stacked y label and y geometry (NOW tree-ordered)
# ============================================================
# Stable within-species ordering: preserve original file order + first-seen chrom_pair
dt[, file_label := factor(file_label, levels = unique(file_label))]
setorder(dt, species, file_label)

dt[, y_label := paste0(as.character(file_label), "  |  ", chrom_pair)]

# Define y_label factor levels in the current dt order (tree-ordered)
dt[, y_label := factor(y_label, levels = unique(y_label))]

dt[, y := as.numeric(y_label)]
dt[, `:=`(ymin = y - 0.45, ymax = y + 0.45)]
n_rows <- length(levels(dt$y_label))

# mean y position for each species (used to align tree tips)
species_y <- dt[, .(y_tip = mean(y)), by = species]
species_y[, species := as.character(species)]

# ----------------------------
# Build chromosome outline boxes for Query / Reference
# ----------------------------
make_chr_boxes <- function(dt, which = c("qry", "ref")) {
  which <- match.arg(which)

  if (which == "qry") {
    tmp <- unique(dt[, list(y_label, y, ymin, ymax, chrom = chrom_qry, seqlen = len_qry)])
  } else {
    tmp <- unique(dt[, list(y_label, y, ymin, ymax, chrom = chrom_ref, seqlen = len_ref)])
  }

  # keep valid lengths only
  tmp <- tmp[!is.na(chrom) & !is.na(seqlen) & seqlen > 0]

  # if duplicates exist per (y_label, chrom), keep the max length
  tmp <- tmp[, .(
    y = y[1],
    ymin = ymin[1],
    ymax = ymax[1],
    seqlen = max(seqlen, na.rm = TRUE)
  ), by = .(y_label, chrom)]

  tmp[, `:=`(xmin = 0, xmax = seqlen)]
  tmp[]
}

boxes_qry <- make_chr_boxes(dt, which = "qry")
boxes_ref <- make_chr_boxes(dt, which = "ref")


# ----------------------------
# Query panel
# ----------------------------
p_qry <- ggplot(dt) +
  geom_rect(
    data = boxes_qry,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA,
    color = "grey80",
    linewidth = 0.3
  ) +
  geom_rect(aes(
    xmin = qs,
    xmax = qe,
    ymin = ymin,
    ymax = ymax,
    fill = "coral"
  )) +
  scale_x_continuous(
  limits = c(0, NA),
  breaks = function(lims) c(lims[1], mean(lims), lims[2]),
  labels = function(lims) c(
    "0",
    "",
    sprintf("%.1f", lims[2] / 1e6)
  )
  )+
  scale_y_continuous(
    breaks = NULL,
    labels = NULL,
    expand = expansion(mult = 0)
  ) +
  labs(
    x = "Position (Mb)",
    y = NULL,
    title = "Query"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(hjust = 1),  # right-justify tick text
    axis.text.x.top = element_text(hjust = 1),
    axis.text.x.bottom = element_text(hjust = 1)
  )

# ----------------------------
# Reference panel
# ----------------------------
p_ref <- ggplot(dt) +
  geom_rect(
    data = boxes_ref,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA,
    color = "grey80",
    linewidth = 0.3
  ) +
  geom_rect(aes(
    xmin = rs,
    xmax = re,
    ymin = ymin,
    ymax = ymax,
    fill = "coral"
  )) +
  scale_x_continuous(
  limits = c(0, NA),
  breaks = function(lims) c(lims[1], mean(lims), lims[2]),
  labels = function(lims) c(
    "0",
    "",
    sprintf("%.1f", lims[2] / 1e6)
  )
  ) +
  scale_y_continuous(
    breaks = seq_len(n_rows),
    labels = levels(dt$y_label),
    position = "right",   # <<< move axis to the right
    expand = expansion(mult = 0)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Reference"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y.left  = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.right = element_text(size = 8),
    axis.ticks.y.right = element_line(),
    legend.position = "none",
    axis.text.x = element_text(hjust = 1),  # right-justify tick text
    axis.text.x.top = element_text(hjust = 1),
    axis.text.x.bottom = element_text(hjust = 1)
  )

# ============================================================
# TREE PANEL (rectangular cladogram; tip labels on the right)
# ============================================================

tr2 <- reorder.phylo(tr, order = "postorder")

Ntip <- length(tr2$tip.label)
Nnode <- tr2$Nnode
total_nodes <- Ntip + Nnode

# y: align each tip to the mean y of its species' rows
y <- rep(NA_real_, total_nodes)
names(y) <- as.character(seq_len(total_nodes))

for (i in seq_len(Ntip)) {
  lab <- tr2$tip.label[i]
  hit <- species_y[species == lab, y_tip]
  if (length(hit) == 1 && !is.na(hit)) {
    y[i] <- hit
  } else {
    stop(paste0("Tip label missing from species_y after pruning: ", lab))
  }
}

# internal node y = mean(children y)
edge <- tr2$edge
parents <- unique(edge[, 1])
for (p in parents) {
  kids <- edge[edge[, 1] == p, 2]
  y[as.character(p)] <- mean(y[as.character(kids)], na.rm = TRUE)
}

# x: cladogram depth (ignore branch lengths)
# build x by accumulating +1 per edge from root
x <- rep(0, total_nodes)
names(x) <- as.character(seq_len(total_nodes))

tr3 <- reorder.phylo(tr2, order = "cladewise")
edge3 <- tr3$edge

for (i in seq_len(nrow(edge3))) {
  p <- edge3[i, 1]
  c <- edge3[i, 2]
  x[as.character(c)] <- x[as.character(p)] + 1
}

# Build *rectangular* segments:
# 1) horizontal: parent (x_p, y_child) -> child (x_c, y_child)
# 2) vertical:   parent (x_p, y_parent) -> parent (x_p, y_child)
seg_h <- data.table(
  x    = x[as.character(edge[, 1])],
  y    = y[as.character(edge[, 2])],
  xend = x[as.character(edge[, 2])],
  yend = y[as.character(edge[, 2])]
)

seg_v <- data.table(
  x    = x[as.character(edge[, 1])],
  y    = y[as.character(edge[, 1])],
  xend = x[as.character(edge[, 1])],
  yend = y[as.character(edge[, 2])]
)

seg_dt <- rbind(seg_h, seg_v)

tips_dt <- data.table(
  label = tr2$tip.label,
  node  = seq_len(Ntip),
  x     = x[as.character(seq_len(Ntip))],
  y     = y[as.character(seq_len(Ntip))]
)

xmax <- max(tips_dt$x)

# Tip labels on the RIGHT of the tree
p_tree <- ggplot() +
  geom_segment(data = seg_dt, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.3) +
  geom_text(
    data = tips_dt,
    aes(x = x, y = y, label = label),
    hjust = 0,                       # left-justify so text extends rightward
    nudge_x = 0.15,                  # push labels to the right of tip
    size = 2.4
  ) +
  scale_x_continuous(limits = c(0, xmax + 1), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0.5, n_rows + 0.5), expand = expansion(mult = 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 12) +
  theme(
    plot.margin = margin(5.5, 80, 5.5, 5.5, "pt")  # extra RIGHT room for labels
  )


# ----------------------------
# Combine + save
# ----------------------------
p <- (p_tree | p_qry | p_ref) +
  plot_layout(widths = c(0.9, 1, 1)) +
  plot_annotation(
    title = "Continuous percent identity: query vs reference (all files) + species tree",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

height_in <- max(20, n_rows * 0.09)

out_png <- sprintf(
  "continuous_percentID.ALLFILES.with_tree_and_chr_boxes.thr%s.len%d.png",
  thr_tag,
  min_aln
)
ggsave(
  filename = out_png,
  plot = p,
  width = 10,          # slightly wider to fit tree labels
  height = height_in,
  units = "in",
  dpi = 200
)

message("Wrote: ", out_png)
message("Rows plotted: ", n_rows,
        " | Boxes (qry/ref): ", nrow(boxes_qry), "/", nrow(boxes_ref))


