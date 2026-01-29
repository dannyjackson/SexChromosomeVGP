suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-s", "--species"), type="character",
              help="Species prefix, e.g. Eubalaena_glacialis"),
  make_option(c("-o", "--outdir"), type="character",
              help="Species output directory, e.g. /.../cetaceans/Eubalaena_glacialis"),
  make_option(c("-n", "--top_n"), type="integer", default=0,
              help="If >0, restrict analysis to top N scaffolds by length (from .fai). Default: 0 = use all.")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$species) || opt$species == "") stop("ERROR: must provide -s SPECIES", call.=FALSE)
if (is.null(opt$outdir)  || opt$outdir  == "") stop("ERROR: must provide -o OUTDIR",  call.=FALSE)

SPECIES <- opt$species
OUTDIR  <- opt$outdir
TOP_N   <- opt$top_n

if (!dir.exists(OUTDIR)) stop(paste("ERROR: outdir does not exist:", OUTDIR), call.=FALSE)

# Input file paths (in OUTDIR)
hap1_kmer_file <- file.path(OUTDIR, paste0(SPECIES, ".hap1-minus-hap2.results"))
hap2_kmer_file <- file.path(OUTDIR, paste0(SPECIES, ".hap2-minus-hap1.results"))

hap1_fai_file  <- file.path(OUTDIR, paste0(SPECIES, ".hap1.fasta.gz.fai"))
hap2_fai_file  <- file.path(OUTDIR, paste0(SPECIES, ".hap2.fasta.gz.fai"))

hap1_bed_file  <- file.path(OUTDIR, paste0(SPECIES, ".hap1-minus-hap2.bed"))
hap2_bed_file  <- file.path(OUTDIR, paste0(SPECIES, ".hap2-minus-hap1.bed"))

needed <- c(hap1_kmer_file, hap2_kmer_file, hap1_fai_file, hap2_fai_file, hap1_bed_file, hap2_bed_file)
missing <- needed[!file.exists(needed)]
if (length(missing) > 0) {
  stop(paste("ERROR: missing required input files:\n", paste(missing, collapse="\n")), call.=FALSE)
}

# Read inputs
hap1_kmer <- read.delim(hap1_kmer_file, header=FALSE)
hap2_kmer <- read.delim(hap2_kmer_file, header=FALSE)

hap1_idx  <- read.delim(hap1_fai_file, header=FALSE)
hap2_idx  <- read.delim(hap2_fai_file, header=FALSE)

hap1_bed  <- read.delim(hap1_bed_file, header=FALSE)
hap2_bed  <- read.delim(hap2_bed_file, header=FALSE)
colnames(hap1_bed) <- c("scaffold","start","end")
colnames(hap2_bed) <- c("scaffold","start","end")

# Keep only scaffolds present in kmer results
hap1_tmp <- subset(hap1_idx, V1 %in% hap1_kmer$V1)
hap2_tmp <- subset(hap2_idx, V1 %in% hap2_kmer$V1)

hap1_merge <- merge(hap1_tmp, hap1_kmer, by="V1", sort=FALSE)
hap2_merge <- merge(hap2_tmp, hap2_kmer, by="V1", sort=FALSE)

# Density = kmers / length (your original V2.y/V2.x)
hap1_merge$density <- with(hap1_merge, V2.y / V2.x)
hap2_merge$density <- with(hap2_merge, V2.y / V2.x)

# Optional top-N by length (V2.x)
if (TOP_N > 0) {
  hap1_merge <- hap1_merge[order(hap1_merge$V2.x, decreasing=TRUE), ]
  hap2_merge <- hap2_merge[order(hap2_merge$V2.x, decreasing=TRUE), ]
  hap1_merge <- head(hap1_merge, TOP_N)
  hap2_merge <- head(hap2_merge, TOP_N)
}

hap1_final <- hap1_merge[, c("V1","V2.x","V2.y","density")]
hap2_final <- hap2_merge[, c("V1","V2.x","V2.y","density")]
hap1_final$Dataset <- "hap1"
hap2_final$Dataset <- "hap2"

# Save a boxplot (non-interactive)
png(file.path(OUTDIR, paste0(SPECIES, ".kmer_density.boxplot.png")), width=1600, height=800, res=200)
boxplot(hap1_final$density, hap2_final$density, names=c("hap1","hap2"), pch=19,
        main=paste0("Kmer density: ", SPECIES), ylab="kmers / bp")
dev.off()

# Hap1-only histograms
xchr_1m <- ggplot(hap1_bed, aes(x=start)) +
  geom_histogram(aes(y=after_stat(density)), color="#5F00F6", fill="black",
                 alpha=0.8, binwidth=1e6) +
  ggtitle(paste0("Kmer Density hap1: ", SPECIES, " (1 Mb bins)")) +
  xlab("Position on chromosome") +
  ylab("Density") +
  theme_bw()

ggsave(
  filename = file.path(OUTDIR, paste0(SPECIES, ".kmer_masked.hap1.1M.png")),
  plot = xchr_1m,
  width = 8, height = 4, dpi = 300
)

xchr_500k <- ggplot(hap1_bed, aes(x=start)) +
  geom_histogram(aes(y=after_stat(density)), color="#5F00F6", fill="black",
                 alpha=0.8, binwidth=500000) +
  ggtitle(paste0("Kmer Density hap1: ", SPECIES, " (500 kb bins)")) +
  xlab("Position on chromosome") +
  ylab("Density") +
  theme_bw()

ggsave(
  filename = file.path(OUTDIR, paste0(SPECIES, ".kmer_masked.hap1.500k.png")),
  plot = xchr_500k,
  width = 8, height = 4, dpi = 300
)


xchr_100k <- ggplot(hap1_bed, aes(x=start)) +
  geom_histogram(aes(y=after_stat(density)), color="#5F00F6", fill="black",
                 alpha=0.8, binwidth=100000) +
  ggtitle(paste0("Kmer Density hap1: ", SPECIES, " (100 kb bins)")) +
  xlab("Position on chromosome") +
  ylab("Density") +
  theme_bw()

ggsave(
  filename = file.path(OUTDIR, paste0(SPECIES, ".kmer_masked.hap1.100k.png")),
  plot = xchr_100k,
  width = 8, height = 4, dpi = 300
)

xchr_10k <- ggplot(hap1_bed, aes(x=start)) +
  geom_histogram(aes(y=after_stat(density)), color="#5F00F6", fill="black",
                 alpha=0.8, binwidth=10000) +
  ggtitle(paste0("Kmer Density hap1: ", SPECIES, " (10 kb bins)")) +
  xlab("Position on chromosome") +
  ylab("Density") +
  theme_bw()

ggsave(
  filename = file.path(OUTDIR, paste0(SPECIES, ".kmer_masked.hap1.10k.png")),
  plot = xchr_10k,
  width = 8, height = 4, dpi = 300
)


xchr_50k <- ggplot(hap1_bed, aes(x=start)) +
  geom_histogram(aes(y=after_stat(density)), color="#5F00F6", fill="black",
                 alpha=0.8, binwidth=50000) +
  ggtitle(paste0("Kmer Density hap1: ", SPECIES, " (50 kb bins)")) +
  xlab("Position on chromosome") +
  ylab("Density") +
  theme_bw()

ggsave(
  filename = file.path(OUTDIR, paste0(SPECIES, ".kmer_masked.hap1.50k.png")),
  plot = xchr_50k,
  width = 8, height = 4, dpi = 300
)

message("Done: wrote plots to ", OUTDIR)
