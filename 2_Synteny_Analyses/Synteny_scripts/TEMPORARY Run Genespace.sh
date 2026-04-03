# Create directory structure
```
SCRIPTS=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/

mkdir -p $SCRIPTS

${SCRIPTS}/1_Create_Directory_Structure.sh

```
## Submit batch job
```
# create output directory for plots
mkdir -p ${OUTDIR}/${SPECIES}/plots

N=$(wc -l < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv)
sbatch --array=1-${N} ${SCRIPTS}/2_genespace_array.sh
```
## Plot just synteny with sex chrom
```
mkdir -p ${OUTDIR}/${SPECIES}/plots_sexchr

N=$(wc -l < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv)
sbatch --array=1-${N} ${SCRIPTS}/3_genespace_array.plot_sex_chr.sh
```

## Create a table of syntenic regions to chicken
### Sex-shared to chicken WG
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/

```
# Combine all *_phasedBlks.csv files into one CSV (one header), adding provenance columns.






































df <- df3
gallus <- "Gallus_gallus_REF"

# Helper to safely compute bp length
bp_len <- function(start, end) {
  s <- suppressWarnings(as.double(start))
  e <- suppressWarnings(as.double(end))
  abs(e - s) + 1
}



bp_len <- function(start, end) {
  s <- suppressWarnings(as.double(start))
  e <- suppressWarnings(as.double(end))
  abs(e - s) + 1
}

# VGP lookup table: "Genus species" -> "Genus_species"
vgp_lut <- VGP_DATAFRAME %>%
  transmute(
    TaxonKey = str_replace_all(str_squish(Scientific.Name), " ", "_"),
    Lineage,
    Superorder
  ) %>%
  distinct(TaxonKey, .keep_all = TRUE)


summary_df <- df %>%
  # only Gallus vs non-Gallus comparisons
  filter((genome1 == gallus & genome2 != gallus) | (genome2 == gallus & genome1 != gallus)) %>%
  mutate(
    Species = if_else(genome1 == gallus, genome2, genome1),
    Chromosome = if_else(genome1 == gallus, chr1, chr2),

    BP_Gallus = if_else(genome1 == gallus,
                        bp_len(startBp1, endBp1),
                        bp_len(startBp2, endBp2)),
    BP_Query = if_else(genome1 == gallus,
                       bp_len(startBp2, endBp2),
                       bp_len(startBp1, endBp1)),

    nHits_Gallus = suppressWarnings(as.double(if_else(genome1 == gallus, nHits1, nHits2))),
    nHits_Query  = suppressWarnings(as.double(if_else(genome1 == gallus, nHits2, nHits1))),

    # join key for VGP (underscored binomial)
    TaxonKey = Species
  ) %>%
  # drop reciprocal duplicates (optional; remove if you want both directions counted)
  distinct(Species, Chromosome, blkID, startBp1, endBp1, startBp2, endBp2, nHits1, nHits2, .keep_all = TRUE) %>%
  group_by(Species, Chromosome, TaxonKey) %>%
  summarise(
    TotalBP_Gallus = sum(BP_Gallus, na.rm = TRUE),
    nHits_Gallus   = sum(nHits_Gallus, na.rm = TRUE),
    TotalBP_Query  = sum(BP_Query, na.rm = TRUE),
    nHits_Query    = sum(nHits_Query, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(vgp_lut, by = "TaxonKey") %>%
  relocate(Lineage, Superorder, .after = Species) %>%
  select(Species, Lineage, Superorder, Chromosome,
         TotalBP_Gallus, nHits_Gallus, TotalBP_Query, nHits_Query) %>%
  arrange(Lineage, Superorder, Species, Chromosome)


# Optional: show which Species didn't match VGP after stripping _REF
unmatched <- summary_df %>%
  filter(is.na(Lineage) | is.na(Superorder)) %>%
  distinct(Species)
if (nrow(unmatched) > 0) {
  cat("\nUnmatched Species (after stripping _REF for lookup):\n")
  print(unmatched)
}
# Umatched:
1 Canis_lupus    
2 Guaruba_guaruba

summary_df <- summary_df %>%
  mutate(
    Lineage = case_when(
      Species == "Canis_lupus" ~ "Mammal",
      Species == "Guaruba_guaruba" ~ "Bird",
      TRUE ~ Lineage
    )
  )

summary_df <- summary_df %>%
  mutate(
    Superorder = case_when(
      Species == "Canis_lupus" ~ "Laurasiatheria",
      Species == "Guaruba_guaruba" ~ "Psittaciformes",
      TRUE ~ Superorder
    )
  )

write_csv(summary_df, "Gallus_vs_species_window_sums_with_lineage.csv")
cat("Wrote Gallus_vs_species_window_sums_with_lineage.csv with", nrow(summary_df), "rows\n")

























```
# Make a plot NOT WORKING YET
```
# Plot synteny blocks on Gallus Z (0..86,044,486 bp) with species on Y

library(readr)
library(dplyr)
library(purrr)
library(stringr)

library(ggplot2)

gallus <- "Gallus_gallus_REF"
chr_target <- "Z"
Z_len <- 86044486

to_num <- function(x) suppressWarnings(as.numeric(x))

synteny_plot_df <- df %>%
  # keep Gallus vs non-Gallus
  filter((genome1 == gallus & genome2 != gallus) | (genome2 == gallus & genome1 != gallus)) %>%
  mutate(
    Species = if_else(genome1 == gallus, genome2, genome1),

    # Gallus chromosome for this row
    GallusChr = if_else(genome1 == gallus, chr1, chr2),

    # Gallus interval for this row
    g_start = if_else(genome1 == gallus, startBp1, startBp2),
    g_end   = if_else(genome1 == gallus, endBp1,   endBp2),

    g_start = to_num(g_start),
    g_end   = to_num(g_end),

    x_start = pmin(g_start, g_end, na.rm = TRUE),
    x_end   = pmax(g_start, g_end, na.rm = TRUE)
  ) %>%
  # only Gallus Z
  filter(as.character(GallusChr) == chr_target) %>%
  # drop reciprocal duplicates (optional but usually desirable)
  distinct(Species, GallusChr, blkID, x_start, x_end, .keep_all = TRUE) %>%
  filter(!is.na(x_start), !is.na(x_end)) %>%
  mutate(Species = factor(Species, levels = sort(unique(as.character(Species)))))

# Ensure A->Z (top to bottom depends on your preference; this is the common "A at top")
synteny_plot_df <- synteny_plot_df %>%
  mutate(Species = factor(as.character(Species), levels = sort(unique(as.character(Species)))))

p <- ggplot(synteny_plot_df, aes(y = Species)) +
  geom_segment(aes(x = x_start, xend = x_end, yend = Species), linewidth = 0.4, alpha = 0.8) +
  coord_cartesian(xlim = c(0, Z_len), clip = "off") +
  scale_y_discrete(limits = levels(synteny_plot_df$Species)) +  # enforce the order
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "Gallus_gallus Z position (bp)",
    y = "Species",
    title = "Synteny blocks to Gallus_gallus Z"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7),
    plot.margin = margin(t = 15, r = 10, b = 15, l = 10)  # <- THIS fixes top/bottom cutoff
  )

p <- p +
  geom_vline(xintercept = 0,  linewidth = 0.2, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = Z_len, linewidth = 0.2, linetype = "dashed", colour = "grey60")

# give the output a touch more vertical room too (optional but helps)
ggsave(
  "synteny_to_gallus_Z.png",
  p,
  width = 10,
  height = max(4, 0.22 * nlevels(synteny_plot_df$Species)), # slightly larger multiplier
  dpi = 300
)


ggsave("synteny_to_gallus_Z.png", p, width = 10, height = max(4, 0.2 * nlevels(synteny_plot_df$Species)), dpi = 300)
```