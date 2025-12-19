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
vgp_file <- Sys.getenv("VGP")
df <- read.csv(vgp_file)       

# 1. "Genome": evaluate Accession column, rename as "Genome"
library(dplyr)

df <- df %>%
  mutate(
    Genome = coalesce(
      Accession_num_main_haplotype,
      RefSeq.annotation.main.haplotype,
      Scientific.Name
    )
  )

# Identify rows where the two accession columns differ
diff_rows <- df[df$Genome != df$Accession_num_main_haplotype, ]

# Select and display desired columns
result <- diff_rows[, c("Genome",
                        "Accession_num_main_haplotype",
                        "Scientific.Name",
                        "RefSeq.annotation.main.haplotype")]

print(result)

# No results found -- all things are the same luckily.

# 2. Genetic_Sex -- recode column "Sex_Chromosome_System_Karyotpe" (which contains more information than necessary) into just "sex is or is not genetic"
df$Genetic_Sex <- with(df, ifelse(
  Sex_Chromosome_System_Karyotype %in% c("ZW", "XY", "complex XY", "XO", "complex ZW", "homomorphic"),
  "Y",
  ifelse(
    Sex_Chromosome_System_Karyotype %in% c("?", "UN", "", " "),
    "U",
    ifelse(
      Sex_Chromosome_System_Karyotype %in% c("None", NA, "<NA>"),
      "N",
      NA  # fallback if an unexpected value appears
    )
  )
))

# 3. All_Sex_Chr -- Use column "Hetergametic_Sequenced" to infer whether all of the sex chromosomes were assembled or not
df$All_Sex_Chr <- df$Hetergametic_Sequenced
df$All_Sex_Chr[df$All_Sex_Chr %in% c("", " ")] <- NA

# Write it out to a new DF
df_new <- df[, c("Genome", "Genetic_Sex", "All_Sex_Chr")]

write.csv(df_new, "VGP_PAR_datafile.csv", row.names = FALSE)
