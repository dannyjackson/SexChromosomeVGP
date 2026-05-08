# Pop structure Gar

mamba create -n popgen -c conda-forge -c bioconda r-base r-essentials bioconductor-gdsfmt bioconductor-snprelate -y
mamba activate popgen
R -q -e "library(gdsfmt); library(SNPRelate); sessionInfo()"
R -q -e "library(data.table); library(ggplot2); sessionInfo()"

```
#!/usr/bin/env Rscript

library(gdsfmt)
library(SNPRelate)

vcf.fn <- "Lepisosteus_oculatus_snps_multiallelic.renamed.contigs.vcf.gz"
sex.file <- "accession_sex.csv"
outDir <- "."
name <- "Lepisosteus_oculatus_PCA"

gds.fn <- file.path(outDir, paste0(name, ".gds"))
summary.fn <- file.path(outDir, paste0(name, "_summary.txt"))
var.fn <- file.path(outDir, paste0(name, "_percentvariation.txt"))
table.fn <- file.path(outDir, paste0(name, ".txt"))
pdf.fn <- file.path(outDir, paste0(name, "_pca_by_sex.pdf"))

snpgdsVCF2GDS(vcf.fn, gds.fn, method = "biallelic.only")

sink(summary.fn)
snpgdsSummary(gds.fn)
sink()

genofile <- snpgdsOpen(gds.fn)

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
snp.id <- read.gdsn(index.gdsn(genofile, "snp.id"))
snp.chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
snp.pos <- read.gdsn(index.gdsn(genofile, "snp.position"))


sex.df <- read.csv(sex.file, header = TRUE, stringsAsFactors = FALSE)


colnames(sex.df) <- c("sample.id", "sex")

sex.df$sample.id <- as.character(sex.df$sample.id)
sex.df$sex <- trimws(as.character(sex.df$sex))

sample.df <- read.csv(sex.file, header = TRUE, stringsAsFactors = FALSE)

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

run_and_plot_pca <- function(genofile, sample.ids, metadata, prefix, main.title) {
  if (length(sample.ids) < 2) {
    cat("Skipping", prefix, "- fewer than 2 samples\n")
    return(NULL)
  }

  pca <- snpgdsPCA(
    genofile,
    sample.id = sample.ids,
    num.thread = 2,
    autosome.only = FALSE
  )

  pc.percent <- pca$varprop * 100

  tab <- data.frame(
    sample.id = pca$sample.id,
    EV1 = pca$eigenvect[, 1],
    EV2 = pca$eigenvect[, 2],
    stringsAsFactors = FALSE
  )

  tab <- merge(tab, metadata, by = "sample.id", all.x = TRUE, sort = FALSE)
  tab <- tab[match(pca$sample.id, tab$sample.id), ]
  tab$sex[is.na(tab$sex) | tab$sex == ""] <- "Unknown"
  tab$sex <- factor(tab$sex)

  write.table(
    tab,
    file = file.path(outDir, paste0(prefix, ".txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    data.frame(PC = seq_along(pc.percent), PercentVariance = round(pc.percent, 2)),
    file = file.path(outDir, paste0(prefix, "_percentvariation.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cols <- seq_along(levels(tab$sex))
  names(cols) <- levels(tab$sex)

  pdf(file = file.path(outDir, paste0(prefix, ".pdf")), useDingbats = FALSE)
  plot(
    tab$EV1, tab$EV2,
    col = cols[as.character(tab$sex)],
    pch = 19,
    xlab = paste0("PC1 (", round(pc.percent[1], 2), "%)"),
    ylab = paste0("PC2 (", round(pc.percent[2], 2), "%)"),
    main = main.title
  )
  legend(
    "topright",
    legend = levels(tab$sex),
    col = cols,
    pch = 19,
    bty = "n"
  )
  dev.off()

  return(tab)
}

all.meta <- sex.df[sex.df$sample.id %in% sample.id, ]

run_and_plot_pca(
  genofile = genofile,
  sample.ids = intersect(sample.id, all.meta$sample.id),
  metadata = all.meta,
  prefix = paste0(name, "_all"),
  main.title = "PCA of all samples colored by sex"
)

male.ids <- all.meta$sample.id[tolower(all.meta$sex) %in% c("male", "m")]
female.ids <- all.meta$sample.id[tolower(all.meta$sex) %in% c("female", "f")]

run_and_plot_pca(
  genofile = genofile,
  sample.ids = intersect(sample.id, male.ids),
  metadata = all.meta[all.meta$sample.id %in% male.ids, ],
  prefix = paste0(name, "_males"),
  main.title = "PCA of male samples"
)

run_and_plot_pca(
  genofile = genofile,
  sample.ids = intersect(sample.id, female.ids),
  metadata = all.meta[all.meta$sample.id %in% female.ids, ],
  prefix = paste0(name, "_females"),
  main.title = "PCA of female samples"
)

cat("Done.\n")
```
# Make a PCA for every chromosome
```
library(data.table)
library(ggplot2)

run_pca_by_chromosome <- function(genofile, sample.ids, metadata, snp.id, snp.chr, chrom.file, outDir, prefix) {
  chrom_map <- fread(
    chrom.file,
    header = FALSE,
    col.names = c("accession", "chrom_plot")
  )

  dir.create(file.path(outDir, "pca_by_chromosome"), showWarnings = FALSE)

  for (i in seq_len(nrow(chrom_map))) {
    chr_acc <- chrom_map$accession[i]
    chr_plot <- chrom_map$chrom_plot[i]

    snp.sel <- which(as.character(snp.chr) == chr_acc)

    if (length(snp.sel) < 2) {
      cat("Skipping chromosome", chr_acc, "- fewer than 2 SNPs\n")
      next
    }

    pca <- snpgdsPCA(
      genofile,
      sample.id = sample.ids,
      snp.id = snp.id[snp.sel],
      num.thread = 2,
      autosome.only = FALSE
    )

    pc.percent <- pca$varprop * 100

    tab <- data.frame(
      sample.id = pca$sample.id,
      EV1 = pca$eigenvect[, 1],
      EV2 = pca$eigenvect[, 2],
      stringsAsFactors = FALSE
    )

    tab <- merge(tab, metadata, by = "sample.id", all.x = TRUE, sort = FALSE)
    tab <- tab[match(pca$sample.id, tab$sample.id), ]
    tab$sex[is.na(tab$sex) | tab$sex == ""] <- "Unknown"

    sex.plot <- tolower(as.character(tab$sex))
    sex.plot[!(sex.plot %in% c("male", "female"))] <- "unknown"

    cols <- c(
      male = "blue",
      female = "black",
      unknown = "grey70"
    )

    write.table(
      tab,
      file = file.path(outDir, "pca_by_chromosome", paste0(prefix, "_chr_", chr_plot, ".txt")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    write.table(
      data.frame(PC = seq_along(pc.percent), PercentVariance = round(pc.percent, 2)),
      file = file.path(outDir, "pca_by_chromosome", paste0(prefix, "_chr_", chr_plot, "_percentvariation.txt")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    pdf(
      file = file.path(outDir, "pca_by_chromosome", paste0(prefix, "_chr_", chr_plot, ".pdf")),
      useDingbats = FALSE
    )
    plot(
      tab$EV1, tab$EV2,
      col = cols[sex.plot],
      pch = 19,
      xlab = paste0("PC1 (", round(pc.percent[1], 2), "%)"),
      ylab = paste0("PC2 (", round(pc.percent[2], 2), "%)"),
      main = paste("Chromosome", chr_plot, "PCA by sex")
    )
    legend(
      "topright",
      legend = c("male", "female", "unknown"),
      col = cols[c("male", "female", "unknown")],
      pch = 19,
      bty = "n"
    )
    dev.off()
  }
}

run_pca_by_chromosome(
  genofile = genofile,
  sample.ids = intersect(sample.id, all.meta$sample.id),
  metadata = all.meta,
  snp.id = snp.id,
  snp.chr = snp.chr,
  chrom.file = "chrom_conversion.txt",
  outDir = outDir,
  prefix = paste0(name, "_bychr")
)
```