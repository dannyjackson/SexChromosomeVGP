#!/usr/bin/env Rscript

alignment_file <- paste0(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/",
  "jacksondan/analyses/PAR_inference/alignment/continuous_percentID/",
  "surviving_regions.qry.thr98p5.len10000.ALLFILES.csv"
)

par_file <- "PAR.species_chr_region.txt"
output_file <- "PAR_alignment_coverage.csv"

# Read alignment intervals.
# sep = "" handles whitespace-delimited input; use sep = "," if truly comma-separated.
alignments <- read.csv(
  alignment_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Read inferred PAR definitions.
par <- read.table(
  par_file,
  sep = ",",
  header = FALSE,
  col.names = c("Species", "Region"),
  stringsAsFactors = FALSE
)

# Parse chromosome, start, and stop.
par$Chromosome <- sub(":.*$", "", par$Region)
coordinates <- sub("^.*:", "", par$Region)

par$Start <- as.numeric(sub("-.*$", "", coordinates))
par$Stop <- as.numeric(sub("^.*-", "", coordinates))


merge_and_sum <- function(starts, stops) {
  if (length(starts) == 0) {
    return(0)
  }

  intervals <- data.frame(Start = starts, Stop = stops)
  intervals <- intervals[order(intervals$Start, intervals$Stop), ]

  current_start <- intervals$Start[1]
  current_stop <- intervals$Stop[1]
  total_bp <- 0

  if (nrow(intervals) > 1) {
    for (i in 2:nrow(intervals)) {
      if (intervals$Start[i] <= current_stop) {
        current_stop <- max(current_stop, intervals$Stop[i])
      } else {
        total_bp <- total_bp + (current_stop - current_start)
        current_start <- intervals$Start[i]
        current_stop <- intervals$Stop[i]
      }
    }
  }

  total_bp + (current_stop - current_start)
}


results <- lapply(seq_len(nrow(par)), function(i) {
  species <- par$Species[i]
  par_start <- par$Start[i]
  par_stop <- par$Stop[i]
  par_length <- par_stop - par_start

  species_alignments <- alignments[alignments$Species == species, ]

  # Clip alignment intervals to the inferred PAR boundaries.
  clipped_start <- pmax(species_alignments$Start, par_start)
  clipped_stop <- pmin(species_alignments$Stop, par_stop)

  keep <- clipped_start < clipped_stop

  aligned_bp <- merge_and_sum(
    clipped_start[keep],
    clipped_stop[keep]
  )

  proportion <- if (par_length > 0) {
    aligned_bp / par_length
  } else {
    NA_real_
  }

  data.frame(
    Species = species,
    Region = paste0(
      par$Chromosome[i], ":",
      par_start, "-", par_stop
    ),
    Aligned_bp = aligned_bp,
    Proportion_aligned = proportion,
    stringsAsFactors = FALSE
  )
})

results <- do.call(rbind, results)

write.table(
  results,
  file = output_file,
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

png(
  filename = "PAR_proportion_aligned_histogram.png",
  width = 1200,
  height = 900,
  res = 150
)

hist(
  results$Proportion_aligned,
  breaks = 20,
  main = "Distribution of PAR Alignment Proportions",
  xlab = "Proportion aligned",
  ylab = "Number of species",
  xlim = c(0, 1)
)

dev.off()