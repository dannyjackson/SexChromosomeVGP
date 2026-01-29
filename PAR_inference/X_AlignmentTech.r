# Robust ALL-by-AssemblyTech plotting script (phylo-ordered)
library(data.table)
library(ggplot2)
library(patchwork)
library(ape)

# ----------------------------
# Params / inputs
# ----------------------------
tree_file <- "/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/roadies_v1.1.16b.numbers.scientific.nwk"
thr <- 98.5
min_aln <- 10000
thr_tag <- gsub("\\.", "p", sprintf("%.1f", thr))
setwd("/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_inference/alignment/AlignmentTech")
dir_in <- "/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/minimap2/continuous_percentID"
files  <- list.files(dir_in, pattern = "\\.refqry\\.csv$", full.names = TRUE)

# ----------------------------
# Read + ladderize tree once (global)
# ----------------------------
tr_full <- ape::read.tree(tree_file)
tr_full <- ape::ladderize(tr_full)
global_tipset <- tr_full$tip.label

# species parsing helper (same heuristic as original)
get_species_from_label <- function(x, tipset = global_tipset) {
  parts <- strsplit(x, "_", fixed = TRUE)[[1]]
  prefs <- vapply(seq_along(parts), function(k) paste(parts[1:k], collapse = "_"), character(1))
  hits <- prefs[prefs %in% tipset]
  if (length(hits) > 0) {
    hits[length(hits)]
  } else if (length(parts) >= 2) {
    paste(parts[1], parts[2], sep = "_")
  } else {
    x
  }
}

# ----------------------------
# read_one (unchanged)
# ----------------------------
read_one <- function(f) {
  raw <- fread(f)
  file_label <- sub("\\.aln\\.refqry\\.csv$", "", basename(f))

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

  dt <- dt[!is.na(pid) & qe > qs & re > rs]
  dt[, aln_qry := qe - qs]
  dt[, aln_ref := re - rs]
  dt[, aln_len := pmin(aln_qry, aln_ref)]
  dt <- dt[aln_len >= min_aln & pid >= thr]

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
# Read + merge ref file
# ----------------------------
dt_list <- lapply(files, read_one)
dt_list <- Filter(Negate(is.null), dt_list)
if (length(dt_list) == 0) stop("No non-empty *.refqry.csv files found.")
dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

ref_csv <- "/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/VGP_OrdinalList_Phase1Freeze_v1.2_Sept.30.2025_sex_chrs_HalfDeep_SCINKD.csv"
ref <- fread(ref_csv)

dt[, Genus_species := sub("^([A-Za-z]+_[A-Za-z]+).*", "\\1", file_label)]
ref2 <- unique(ref[, .(Genus_species, Assembly.tech)])
dt <- merge(dt, ref2, by = "Genus_species", all.x = TRUE)

# report unmapped
unmapped <- unique(dt[is.na(Assembly.tech), .(file_label, Genus_species)])
if (nrow(unmapped) > 0) {
  message("WARNING: some files have no Assembly.tech mapping (they will be grouped under 'NA').")
  print(head(unmapped, 20))
}

# ----------------------------
# normalize Assembly.tech (you had this)
# ----------------------------
clean_tech_string <- function(x) {
  x2 <- ifelse(is.na(x), "NA", as.character(x))
  x2 <- tolower(x2)
  x2 <- trimws(x2)
  x2 <- gsub("[[:space:]]+", " ", x2)
  x2 <- gsub("\\s*,\\s*", ", ", x2)
  x2 <- gsub("\\+", " + ", x2)
  x2 <- gsub("\\band\\b", "and", x2)
  x2 <- gsub("[^[:alnum:],+\\s\\-\\.]", "", x2)
  x2 <- trimws(x2)
  x2 <- gsub("\\bdamar\\b", "damar", x2)
  x2 <- gsub("hifiams", "hifiasm", x2)
  x2
}
dt[, AssemblyTech_clean := clean_tech_string(Assembly.tech)]

# ----------------------------
# Simplified grouping (CLR anywhere, hifiasm-solo, NA, other)
# ----------------------------
is_clr    <- grepl("clr", dt$Assembly.tech, ignore.case = TRUE)
is_hifiasm_solo <- grepl("hifiasm", dt$Assembly.tech, ignore.case = TRUE) & grepl("solo", dt$Assembly.tech, ignore.case = TRUE)
is_na     <- is.na(dt$Assembly.tech)
is_other  <- !(is_clr | is_hifiasm_solo | is_na)

group_map <- list(
  CLR_anywhere     = which(is_clr),
  Hifiasm_solo     = which(is_hifiasm_solo),
  Assembly_NA      = which(is_na),
  Everything_else  = which(is_other)
)

# ----------------------------
# safe chr boxes (same as you had)
# ----------------------------
make_chr_boxes <- function(dt_sub, which = c("qry", "ref")) {
  which <- match.arg(which)
  if (which == "qry") {
    tmp <- unique(dt_sub[, .(y_label, y, ymin, ymax, chrom = chrom_qry, seqlen = len_qry)])
  } else {
    tmp <- unique(dt_sub[, .(y_label, y, ymin, ymax, chrom = chrom_ref, seqlen = len_ref)])
  }
  if (nrow(tmp) == 0) return(NULL)
  tmp <- tmp[!is.na(chrom) & !is.na(seqlen) & seqlen > 0]
  if (nrow(tmp) == 0) return(NULL)
  tmp <- tmp[, .(y = y[1], ymin = ymin[1], ymax = ymax[1], seqlen = max(seqlen, na.rm = TRUE)), by = .(y_label, chrom)]
  tmp[, `:=`(xmin = 0, xmax = seqlen)]
  tmp[]
}

# ----------------------------
# Plot function: uses phylogenetic ordering for samples
# ----------------------------
plot_simple_group <- function(dt_sub, label, out_dir,
                              per_row = 0.09, max_height = 40, tr_global = tr_full) {

  # annotate species (using global tipset)
  dt_sub[, species := vapply(as.character(file_label),
                             function(x) get_species_from_label(x, tr_global$tip.label),
                             character(1))]

  # determine species present in the tree and those not in tree
  species_present <- intersect(unique(dt_sub$species), tr_global$tip.label)
  species_not_in_tree <- setdiff(unique(dt_sub$species), tr_global$tip.label)

  # build species order: tree-order first (via pruned tree), then any non-tree species in their original appearance order
  species_order_tree <- character(0)
  if (length(species_present) > 0) {
    trp <- ape::keep.tip(tr_global, species_present)
    trp <- ape::ladderize(trp)
    species_order_tree <- trp$tip.label
  }
  species_final <- c(species_order_tree, species_not_in_tree)

  # set species factor to the phylogenetic order
  dt_sub[, species := factor(as.character(species), levels = species_final)]

  # Now order rows by species (phylo) then by file_label to preserve within-species file order
  # But DO NOT set file_label factor levels yet (we want appearance order after this ordering)
  data.table::setorder(dt_sub, species, file_label)

  # Now set file_label factor levels to the appearance order (so y_label preserves this)
  dt_sub[, file_label := factor(file_label, levels = unique(as.character(dt_sub$file_label)))]

  # create y_label and geometry
  dt_sub[, y_label := paste0(as.character(file_label), "  |  ", chrom_pair)]
  dt_sub[, y_label := factor(y_label, levels = unique(y_label))]
  dt_sub[, y := as.numeric(y_label)]
  dt_sub[, `:=`(ymin = y - 0.45, ymax = y + 0.45)]
  n_rows <- length(levels(dt_sub$y_label))
  if (n_rows == 0) {
    message("No rows to plot for group: ", label, " — skipping.")
    return(NULL)
  }

  # chr boxes + x-limits
  boxes_qry <- make_chr_boxes(dt_sub, "qry")
  boxes_ref <- make_chr_boxes(dt_sub, "ref")

  all_x <- c(
    if (!is.null(boxes_qry)) boxes_qry$xmax else NA_real_,
    if (!is.null(boxes_ref)) boxes_ref$xmax else NA_real_,
    as.numeric(dt_sub$qs), as.numeric(dt_sub$qe), as.numeric(dt_sub$rs), as.numeric(dt_sub$re)
  )
  all_x <- all_x[is.finite(all_x) & !is.na(all_x)]
  max_x <- if (length(all_x) > 0) max(all_x) else 1
  xlim_left <- 0; xlim_right <- max(1, max_x)
  mid_x <- (xlim_left + xlim_right) / 2

  # Query and Ref panels
  p_qry <- ggplot() +
    { if (!is.null(boxes_qry)) geom_rect(data = boxes_qry, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = NA, color = "grey80", linewidth = 0.3) } +
    { if (nrow(dt_sub) > 0) geom_rect(data = dt_sub, aes(xmin = qs, xmax = qe, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "coral") } +
    annotate("text", x = xlim_left, y = max(dt_sub$ymax, na.rm = TRUE) + 0.3, label = label, hjust = 0, vjust = 1, size = 4, fontface = "bold") +
    scale_x_continuous(limits = c(xlim_left, xlim_right), breaks = c(xlim_left, mid_x, xlim_right), labels = c("0", "", sprintf("%.1f", xlim_right / 1e6))) +
    scale_y_continuous(expand = expansion(mult = 0)) +
    labs(x = "Query position (Mb)", y = NULL, title = "Query") +
    theme_bw() + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

  p_ref <- ggplot() +
    { if (!is.null(boxes_ref)) geom_rect(data = boxes_ref, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = NA, color = "grey80", linewidth = 0.3) } +
    { if (nrow(dt_sub) > 0) geom_rect(data = dt_sub, aes(xmin = rs, xmax = re, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "coral") } +
    scale_x_continuous(limits = c(xlim_left, xlim_right), breaks = c(xlim_left, mid_x, xlim_right), labels = c("0", "", sprintf("%.1f", xlim_right / 1e6))) +
    scale_y_continuous(breaks = seq_len(n_rows), labels = levels(dt_sub$y_label), position = "right", expand = expansion(mult = 0)) +
    labs(title = "Reference") +
    theme_bw() + theme(panel.grid = element_blank(), axis.text.y.left = element_blank(), axis.text.y.right = element_text(size = 8))

  # Tree panel construction: align tip y positions to species mean y
  species_y <- dt_sub[, .(y_tip = mean(y, na.rm = TRUE)), by = species]
  species_y[, species := as.character(species)]

  present_species <- intersect(species_final, tr_global$tip.label)
  if (length(present_species) > 0) {
    trp <- ape::keep.tip(tr_global, present_species)
    trp <- ape::ladderize(trp)

    tr2 <- ape::reorder.phylo(trp, "postorder")
    Ntip <- length(tr2$tip.label)
    total_nodes <- Ntip + tr2$Nnode

    ycoords <- rep(NA_real_, total_nodes)
    for (i in seq_len(Ntip)) {
      lab <- tr2$tip.label[i]
      hit <- species_y[species == lab, y_tip]
      if (length(hit) != 1 || is.na(hit)) hit <- i
      ycoords[i] <- hit
    }

    edge <- tr2$edge
    if (nrow(edge) > 0) {
      for (p in unique(edge[,1])) {
        kids <- edge[edge[,1] == p, 2]
        ycoords[p] <- mean(ycoords[kids], na.rm = TRUE)
      }
      tr3 <- ape::reorder.phylo(tr2, "cladewise")
      xcoords <- rep(0, total_nodes)
      for (i in seq_len(nrow(tr3$edge))) xcoords[tr3$edge[i,2]] <- xcoords[tr3$edge[i,1]] + 1

      seg_dt <- rbind(
        data.table(x = xcoords[edge[,1]], y = ycoords[edge[,2]], xend = xcoords[edge[,2]], yend = ycoords[edge[,2]]),
        data.table(x = xcoords[edge[,1]], y = ycoords[edge[,1]], xend = xcoords[edge[,1]], yend = ycoords[edge[,2]])
      )
    } else {
      seg_dt <- data.table(x=numeric(), y=numeric(), xend=numeric(), yend=numeric())
      xcoords <- rep(0, total_nodes)
    }

    tips_dt <- data.table(label = tr2$tip.label, x = xcoords[seq_len(Ntip)], y = ycoords[seq_len(Ntip)])
    xmax <- if (nrow(tips_dt) > 0) max(tips_dt$x) else 0

    p_tree <- ggplot() +
      { if (nrow(seg_dt) > 0) geom_segment(data = seg_dt, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.3) } +
      { if (nrow(tips_dt) > 0) geom_text(data = tips_dt, aes(x = x, y = y, label = label), hjust = 0, nudge_x = 0.15, size = 2.4) } +
      scale_x_continuous(limits = c(0, xmax + 1), expand = expansion(mult = 0)) +
      scale_y_continuous(limits = c(0.5, max(1, n_rows) + 0.5), expand = expansion(mult = 0)) +
      coord_cartesian(clip = "off") +
      theme_void(base_size = 12) +
      theme(plot.margin = margin(5.5, 40, 5.5, 5.5, "pt"))
  } else {
    # fallback empty tree panel
    p_tree <- ggplot() + theme_void()
  }

  combined <- (p_tree | p_qry | p_ref) + plot_layout(widths = c(0.9, 1, 1))
  height_in <- max(3.5, n_rows * per_row)
  height_in <- min(height_in, max_height)

  out_file <- file.path(out_dir, paste0("continuous_percentID.", gsub("[^A-Za-z0-9_]", "_", label), ".png"))
  ggsave(out_file, combined, width = 12, height = height_in, units = "in", dpi = 200, limitsize = FALSE)
  message("Wrote: ", out_file, " | rows = ", n_rows, " | height_in = ", sprintf("%.2f", height_in))
  return(out_file)
}

# ----------------------------
# Run for your 4 groups
# ----------------------------
results <- list()
for (grp_name in names(group_map)) {
  idx <- group_map[[grp_name]]
  if (length(idx) == 0) {
    message("Group ", grp_name, " has 0 files: skipping.")
    next
  }
  dt_sub <- dt[idx]
  results[[grp_name]] <- tryCatch({
    plot_simple_group(dt_sub, grp_name, out_dir = getwd(), tr_global = tr_full)
  }, error = function(e) {
    message("Error plotting group ", grp_name, " : ", e$message)
    NULL
  })
}

print(results)
