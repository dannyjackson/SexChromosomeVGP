# Continuous alignment plotting

setwd("/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_inference/alignment/continuous_percentID")

dir_in <- "/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/minimap2/continuous_percentID"
files <- list.files(dir_in, pattern = "\\.refqry\\.csv$", full.names = TRUE)

library(data.table)
library(ggplot2)
library(patchwork)
library(readr)

infile <- "/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/minimap2/continuous_percentID/Homo_sapiens_YtoX.aln.refqry.csv"

dt <- fread(infile)

info <- read_tsv("/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_inference/alignment/bed_sexchr_info_with_matches.tsv",
                 show_col_types = FALSE)

# Coerce + keep only needed columns
dt <- dt[, .(
  chrom_qry = as.character(chrom_qry),
  qs = as.integer(bp_start_qry),
  qe = as.integer(bp_end_qry),
  chrom_ref = as.character(chrom_ref),
  rs = as.integer(bp_start_ref),
  re = as.integer(bp_end_ref),
  pid = as.numeric(percent_identity_qry)
)]

# Basic cleaning
dt <- dt[!is.na(pid) & qe > qs & re > rs]

# Factor ordering (keeps rows aligned between panels)
dt[, chrom_pair := factor(
  paste(chrom_qry, "→", chrom_ref),
  levels = unique(paste(chrom_qry, "→", chrom_ref))
)]

# Numeric y for rectangles
dt[, y := as.numeric(chrom_pair)]
dt[, `:=`(ymin = y - 0.45, ymax = y + 0.45)]

# ----------------------------
# Query panel
# ----------------------------
p_qry <- ggplot(dt) +
  geom_rect(aes(
    xmin = qs,
    xmax = qe,
    ymin = ymin,
    ymax = ymax,
    fill = pid
  )) +
  scale_y_continuous(
    breaks = seq_along(levels(dt$chrom_pair)),
    labels = levels(dt$chrom_pair),
    expand = expansion(mult = 0)
  ) +
  labs(
    x = "Query position (bp)",
    y = NULL,
    title = "Query"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 9)
  )

# ----------------------------
# Reference panel
# ----------------------------
p_ref <- ggplot(dt) +
  geom_rect(aes(
    xmin = rs,
    xmax = re,
    ymin = ymin,
    ymax = ymax,
    fill = pid
  )) +
  scale_y_continuous(
    breaks = seq_along(levels(dt$chrom_pair)),
    labels = NULL,
    expand = expansion(mult = 0)
  ) +
  labs(
    x = "Reference position (bp)",
    y = NULL,
    title = "Reference"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank()
  )

# ----------------------------
# Combine panels
# ----------------------------
p <- (p_qry | p_ref) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Continuous percent identity: query vs reference",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )


png("continuous_percentID.Homo_sapiens.png", width = 2000, height = 1000, res = 200)
print(p)
dev.off()

png("continuous_percentID.png", width = 2000, height = 9000, res = 200)
print(p)
dev.off()



















setwd("/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_inference/alignment/continuous_percentID")

dir_in <- "/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/minimap2/continuous_percentID"
files <- list.files(dir_in, pattern = "\\.refqry\\.csv$", full.names = TRUE)

library(data.table)
library(ggplot2)
library(patchwork)

# ----------------------------
# Read + standardize all files
# ----------------------------
read_one <- function(f) {
  dt <- fread(f)

  # Coerce + keep only needed columns
  dt <- dt[, .(
    chrom_qry = as.character(chrom_qry),
    qs        = as.integer(bp_start_qry),
    qe        = as.integer(bp_end_qry),
    chrom_ref = as.character(chrom_ref),
    rs        = as.integer(bp_start_ref),
    re        = as.integer(bp_end_ref),
    pid       = as.numeric(percent_identity_qry)
  )]

  # Basic cleaning
  dt <- dt[!is.na(pid) & qe > qs & re > rs]
  if (nrow(dt) == 0) return(NULL)

  # File label (strip extension)
  dt[, file_label := sub("\\.aln\\.refqry\\.csv$", "", basename(f))]

  # Within-file chromosome pair label
  dt[, chrom_pair := paste(chrom_qry, "→", chrom_ref)]

  dt[]
}

dt_list <- lapply(files, read_one)
dt_list <- Filter(Negate(is.null), dt_list)

if (length(dt_list) == 0) stop("No non-empty *.refqry.csv files found after filtering.")

dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# ----------------------------
# Build stacked y labels: file then chrom_pair
# ----------------------------
dt[, y_label := paste0(file_label, "  |  ", chrom_pair)]

# Keep stable ordering: by file, then first-seen chrom_pair order within each file
# (this avoids alphabetic shuffling across files)
dt[, file_label := factor(file_label, levels = unique(file_label))]
setorder(dt, file_label)

# make y_label a factor in appearance order
dt[, y_label := factor(y_label, levels = unique(y_label))]

# Numeric y for rectangles
dt[, y := as.numeric(y_label)]
dt[, `:=`(ymin = y - 0.45, ymax = y + 0.45)]

n_rows <- length(levels(dt$y_label))

# ----------------------------
# Query panel
# ----------------------------
p_qry <- ggplot(dt) +
  geom_rect(aes(
    xmin = qs,
    xmax = qe,
    ymin = ymin,
    ymax = ymax,
    fill = pid
  )) +
  scale_y_continuous(
    breaks = seq_len(n_rows),
    labels = levels(dt$y_label),
    expand = expansion(mult = 0)
  ) +
  labs(
    x = "Query position (bp)",
    y = NULL,
    title = "Query"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8)
  )

# ----------------------------
# Reference panel
# ----------------------------
p_ref <- ggplot(dt) +
  geom_rect(aes(
    xmin = rs,
    xmax = re,
    ymin = ymin,
    ymax = ymax,
    fill = pid
  )) +
  scale_y_continuous(
    breaks = seq_len(n_rows),
    labels = NULL,
    expand = expansion(mult = 0)
  ) +
  labs(
    x = "Reference position (bp)",
    y = NULL,
    title = "Reference"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank()
  )

# ----------------------------
# Combine panels
# ----------------------------
p <- (p_qry | p_ref) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Continuous percent identity: query vs reference (all files)",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# ----------------------------
# Save: height scales with number of rows
# ----------------------------
# ~18 px per row at res=200 => 0.09 inches per row; tweak if you want denser/taller
height_in <- max(8, n_rows * 0.09)

ggsave(
  filename = "continuous_percentID.ALLFILES.png",
  plot = p,
  width = 12,
  height = height_in,
  units = "in",
  dpi = 200
)
