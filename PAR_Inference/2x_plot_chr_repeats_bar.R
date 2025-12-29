#!/usr/bin/env Rscript

# Usage (old):
#   Rscript 2x_plot_chr_repeats_bar.R TE.bed chr_len out.png
#
# Usage (new, with highlight regions):
#   Rscript 2x_plot_chr_repeats_bar.R TE.bed par.bed chr_len out.png

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 3) {
  bedfile <- args[1]
  parfile <- NA
  chr_len <- as.numeric(args[2])
  out_png <- args[3]
} else if (length(args) == 4) {
  bedfile <- args[1]
  parfile <- args[2]
  chr_len <- as.numeric(args[3])
  out_png <- args[4]
} else {
  stop("Need 3 or 4 args:\n  <TE_BED> [par.bed] <chr_len> <out_png>\n")
}

# Read TE bed
z <- read.table(bedfile, sep = "\t", header = FALSE,
                stringsAsFactors = FALSE, quote = "")
colnames(z)[1:3] <- c("chr", "start", "end")

# Basic sanity
z <- z[!is.na(z$start) & !is.na(z$end), ]
z$start <- as.numeric(z$start)
z$end   <- as.numeric(z$end)

# Clamp to [0, chr_len]
z$start <- pmax(0, pmin(z$start, chr_len))
z$end   <- pmax(0, pmin(z$end, chr_len))
z <- z[z$end > z$start, ]

# Optional: read PAR/regions-to-highlight bed
par_regions <- NULL
if (!is.na(parfile) && file.exists(parfile)) {
  par_regions <- read.table(parfile, sep = "\t", header = FALSE,
                            stringsAsFactors = FALSE, quote = "")
  colnames(par_regions)[1:3] <- c("chr", "start", "end")
  par_regions <- par_regions[!is.na(par_regions$start) & !is.na(par_regions$end), ]
  par_regions$start <- as.numeric(par_regions$start)
  par_regions$end   <- as.numeric(par_regions$end)

  # If TE bed might contain multiple chr, only keep matching chr
  if (nrow(z) > 0 && !is.na(z$chr[1])) {
    par_regions <- par_regions[par_regions$chr == z$chr[1], , drop = FALSE]
  }

  # Clamp + filter
  par_regions$start <- pmax(0, pmin(par_regions$start, chr_len))
  par_regions$end   <- pmax(0, pmin(par_regions$end, chr_len))
  par_regions <- par_regions[par_regions$end > par_regions$start, , drop = FALSE]

  if (nrow(par_regions) == 0) par_regions <- NULL
}

# PNG device
png(out_png, width = 2400, height = 300, res = 300)
par(mar = c(2, 1, 2, 1))

plot(NA,
     xlim = c(0, chr_len),
     ylim = c(0, 1),
     xaxs = "i", yaxs = "i",
     xaxt = "n", yaxt = "n",
     xlab = sprintf("Position on %s (bp)", ifelse(nrow(z) > 0, z$chr[1], "chr")),
     ylab = "",
     main = sprintf("Repeats on %s (length = %d bp)",
                    ifelse(nrow(z) > 0, z$chr[1], "chr"), chr_len))

# Chromosome "backbone"
y0 <- 0.40
y1 <- 0.60
rect(0, y0, chr_len, y1, border = "black", col = "white", lwd = 2)

# Draw PAR/highlight regions first (behind repeats)
if (!is.null(par_regions)) {
  rect(par_regions$start, y0, par_regions$end, y1,
       col = adjustcolor("darkorange", alpha.f = 0.5),
       border = NA)
  # optional outline:
  rect(par_regions$start, y0, par_regions$end, y1,
       border = adjustcolor("darkorange", alpha.f = 1),
       col = NA, lwd = 1)
}

# Draw repeats on top
rect(z$start, y0, z$end, y1, col = "red", border = NA)

# Axis ticks every 10 Mb
tick_step <- 10e6
ticks <- seq(0, chr_len, by = tick_step)
axis(1, at = ticks, labels = format(ticks/1e6, trim = TRUE), lwd = 0, lwd.ticks = 1)
mtext("Mb", side = 1, line = 2.2, adj = 1)

dev.off()
