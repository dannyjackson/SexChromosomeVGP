#!/usr/bin/env Rscript

# Usage:
#   Rscript 2x_plot_Z_repeats_bar.R Z_TE_BED 87715099 out.png
#
# Example:
#   Rscript plot_Z_repeats_bar.R "$Z_TE_BED" 87715099 Z_repeats_bar.png

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Need 3 args: <Z_TE_BED> <chr_len> <out_png>\n")
}
bedfile <- args[1]
chr_len <- as.numeric(args[2])
out_png <- args[3]

# Read the file (tab-delimited, no header)
# Columns: 1=chr, 2=start, 3=end, ...
z <- read.table(bedfile, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")

colnames(z)[1:3] <- c("chr", "start", "end")

# Basic sanity
z <- z[!is.na(z$start) & !is.na(z$end), ]
z$start <- as.numeric(z$start)
z$end   <- as.numeric(z$end)

# Clamp to [0, chr_len]
z$start <- pmax(0, pmin(z$start, chr_len))
z$end   <- pmax(0, pmin(z$end, chr_len))
z <- z[z$end > z$start, ]

# PNG device
png(out_png, width = 2400, height = 300, res = 300)
par(mar = c(2, 1, 2, 1))


# Set up blank plotting region
plot(NA,
     xlim = c(0, chr_len),
     ylim = c(0, 1),
     xaxs = "i", yaxs = "i",
     xaxt = "n", yaxt = "n",
     xlab = "Position on Z (bp)",
     ylab = "",
     main = sprintf("Repeats on %s (length = %d bp)", z$chr[1], chr_len))

# Chromosome "backbone" rectangle (light outline)
y0 <- 0.40
y1 <- 0.60
rect(0, y0, chr_len, y1, border = "black", col = "white", lwd = 2)

# Draw repeats as red blocks inside the backbone
# (No border so dense regions look continuous)
rect(z$start, y0, z$end, y1, col = "red", border = NA)

# Axis: ticks every 10 Mb (adjust if you want)
tick_step <- 10e6
ticks <- seq(0, chr_len, by = tick_step)
axis(1, at = ticks, labels = format(ticks/1e6, trim = TRUE), lwd = 0, lwd.ticks = 1)
mtext("Mb", side = 1, line = 2.2, adj = 1)

# Optional: label PAR-ish telomeric zone (if useful)
# text(chr_len * 0.98, 0.75, "telomere", cex = 0.8, adj = 1)

dev.off()
