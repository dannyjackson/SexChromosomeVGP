#!/usr/bin/env Rscript

# Stacked 3-track BED plot (shared X axis)
#
# Usage:
#   Rscript 3x_plot_chr_regions_stack.R bed1 bed2 bed3 chr_len out.png
#
# Example (your call):
#   Rscript 3x_plot_chr_regions_stack.R par.Y.bed $Y_TE_BED ${Y_ALIGN_BED} 62460029 Y.TE.png

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Need 5 args: <bed1> <bed2> <bed3> <chr_len> <out_png>\n")
}

bed1 <- args[1]
bed2 <- args[2]
bed3 <- args[3]
chr_len <- as.numeric(args[4])
out_png <- args[5]

read_bed3 <- function(path, chr_len, force_chr = NULL) {
  x <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")
  if (ncol(x) < 3) stop(sprintf("BED file has < 3 columns: %s", path))
  x <- x[, 1:3]
  colnames(x) <- c("chr", "start", "end")

  x <- x[!is.na(x$start) & !is.na(x$end), , drop = FALSE]
  x$start <- as.numeric(x$start)
  x$end   <- as.numeric(x$end)

  if (!is.null(force_chr)) {
    x <- x[x$chr == force_chr, , drop = FALSE]
  }

  # Clamp and filter
  x$start <- pmax(0, pmin(x$start, chr_len))
  x$end   <- pmax(0, pmin(x$end, chr_len))
  x <- x[x$end > x$start, , drop = FALSE]

  x
}

# Read first bed to get chromosome label; then filter others to same chr (like prior script)
b1 <- read_bed3(bed1, chr_len, force_chr = NULL)
chr_name <- if (nrow(b1) > 0) b1$chr[1] else NA

b2 <- read_bed3(bed2, chr_len, force_chr = if (!is.na(chr_name)) chr_name else NULL)
b3 <- read_bed3(bed3, chr_len, force_chr = if (!is.na(chr_name)) chr_name else NULL)

track_titles <- c("PAR regions", "TEs", "Alignments")
track_colors <- c(
  adjustcolor("dodgerblue", alpha.f = 0.35),
  adjustcolor("red",        alpha.f = 0.85),
  adjustcolor("forestgreen",alpha.f = 0.60)
)

beds <- list(b1, b2, b3)

# Plot
png(out_png, width = 2400, height = 900, res = 300)

# More left margin so track labels fit; only bottom panel gets x-axis labels
par(mfrow = c(3, 1), mar = c(1.2, 1.0, 2.0, 1.0), oma = c(3.0, 0, 2.0, 0))

# Shared geometry for all panels
y0 <- 0.30
y1 <- 0.70

for (i in 1:3) {
  xa <- "n"
  plot(NA,
     xlim = c(0, chr_len),
     ylim = c(0, 1),
     xaxs = "i", yaxs = "i",
     xaxt = "n",
     yaxt = "n",   # <-- this turns off y-axis + labels
     xlab = "",
     ylab = "",
     main = track_titles[i])


  # Backbone per track
  rect(0, y0, chr_len, y1, border = "black", col = "white", lwd = 2)

  # Draw regions IN ORDER OF ENTRY (important): rect is vectorized and preserves row order
  if (nrow(beds[[i]]) > 0) {
    rect(beds[[i]]$start, y0, beds[[i]]$end, y1,
         col = track_colors[i], border = NA)
  }

}

# Shared x-axis (bottom panel already has axis ticks)
tick_step <- 10e6
ticks <- seq(0, chr_len, by = tick_step)
axis(1, at = ticks, labels = format(ticks / 1e6, trim = TRUE), lwd = 0, lwd.ticks = 1)
mtext("Mb", side = 1, line = 1.8, adj = 1, outer = TRUE)

# Overall title (outer margin)
overall_chr <- if (!is.na(chr_name)) chr_name else "chr"
mtext(sprintf("Regions on %s (length = %d bp)", overall_chr, chr_len),
      side = 3, line = 0.5, outer = TRUE, cex = 1.1)

dev.off()
