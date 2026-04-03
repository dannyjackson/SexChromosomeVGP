# Run Genespace
## Define reference genomes
save as /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/reference_lists/ref_species.txt
```
Homo_sapiens
Gallus_gallus
Taeniopygia_guttata
Anolis_sagrei
Pristis_pectinata
Gasterosteus_aculeatus
```
# Move old dirs out of the way
```
while read -r REF; do
    echo ${REF}
    mv /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${REF} /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/trash
done < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/reference_lists/ref_species.txt

sbatch --job-name=rm_trash --time=01:00:00 --mem=1G --cpus-per-task=1 --wrap="rm -rf /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/trash"
```
# Create requisite genespace directory structure for all NCBI-annotated genomes
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts

source myconda
mamba activate genespace_py3.10

chmod +x 1a_make_genespace_dirs.sh

# test on Gallus
./1a_make_genespace_dirs.sh Gallus_gallus

# make sure all chrs are in numeric format
awk '{print $1}' ../Gallus_gallus/sexshared/*/bed/*bed | sort -u

# run on the rest
while IFS= read -r REF; do
  ./1a_make_genespace_dirs.sh "$REF"
done < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/reference_lists/ref_species.noGal.txt
```
## Submit batch job
```
awk -F',' 'NR>1 && $2 ~ /^(X|Z)/ {print $1 "\t" $2}' \
/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv

# test on Gallus
N=$(wc -l < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv)
sbatch --array=1-${N} 1b_genespace_array.sh Gallus_gallus

# submit for the rest
while IFS= read -r REF; do
    N=$(wc -l < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv)
    sbatch --array=1-${N} 1b_genespace_array.sh ${REF}
done < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/reference_lists/ref_species.noGal.txt
```

## Create a table of syntenic regions to each reference genome
```
Rscript 1e_combine_phsblks.r Gallus_gallus
while read -r REF; do
    echo ${REF}
    Rscript 1e_combine_phsblks.r ${REF}
done < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/reference_lists/ref_species.noGal.txt
```
# Make adjustments for species with newly identified sex chromosomes
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/combined_phaseblks/
FILE=Gallus_gallus_syntentic_chromosomes.sex_shared.csv

sed -i 's/Accipiter_gentilis/Astur_gentilis/g' $FILE
sed -i 's/Lagenorhynchus_acutus/Leucopleurus_acutus/g' $FILE
sed -i 's/Eptesicus_nilssonii/Cnephaeus_nilssonii/g' $FILE
echo 'Narcine_bancroftii,X,"14,27,34",,' >> $FILE

## Accipiter gentilis
## Lagenorhynchus acutus > Leucopleurus acutus
## Eptesicus nilssonii > Cnephaeus nilssonii
## Narcine bancroftii -- chr 12

# Anniella stebbinsi
```
# Plot it on the phylogeny
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts
# Run this interactively on one species to double check that all species are being matched on the phylogeny
# This checks: What species have RefChr data but are not in the tree?
setdiff(unique(grid_df$Species), tr_ultra$tip.label) 
# [1] "Accipiter_gentilis"   "Amazona_ochrocephala" "Canis_lupus"         
# [4] "Pongo_abelii"         "Sciurus_carolinensis"

"Sciurus_carolinensis" %in% tr_ultra$tip.label
grep("Sciurus", tr_ultra$tip.label, value = TRUE)

"Pongo_abelii" %in% tr_ultra$tip.label
grep("Pongo", tr_ultra$tip.label, value = TRUE)

# test on Gallus
Rscript 1f_plot_genespace.r Gallus_gallus

# run on other refs
while read -r REF; do
    echo ${REF}
    Rscript 1f_plot_genespace.r ${REF}
done < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/reference_lists/ref_species.noGal.txt
```
# Identify why genespace did not run on some species
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/

GAL=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/


## NCBI annotations
### No genespace dir at all
Panthera onca
### Has genespace dir
#### Just bed and peptide
Pipistrellus nathusii # no bed file
Lemur catta

## Genomeark annotations
### No genespace dir at all
Anniella stebbinsi
Monodon monocero
Lycocorax pyrrhopterus
Corynorhinus townsendii
Doryrhina cyclops
Phascolarctos cinereus
Tupaia tana

### Has genespace dir
#### Just bed and peptide
Willisornis vidua
Larus fuscus
Anas platyrhynchos
Inia geoffrensis
Artibeus lituratus
Myotis emarginatus
Myotis mystacinus
Pipistrellus pygmaeus
Vespertilio murinus
Podarcis siculus
#### All subdirs
Echiichthys vipera


# Rerun:
PAIRFILE=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv
PAIRFILE_TRBL=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.troubleshoot.csv
grep 'Echiichthys_vipera' $PAIRFILE > ${PAIRFILE_TRBL}
grep 'Chlamydotis_macqueenii' $PAIRFILE >> ${PAIRFILE_TRBL}
grep 'Girardinichthys_multiradiatus' $PAIRFILE >> ${PAIRFILE_TRBL}
grep 'Amblyraja_radiata' $PAIRFILE >> ${PAIRFILE_TRBL}
grep 'Tursiops_truncatus' $PAIRFILE >> ${PAIRFILE_TRBL}

N=$(wc -l < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.troubleshoot.csv)
sbatch --array=1-${N} 1b_genespace_array.troubleshoot.sh Gallus_gallus

## Neo sex
grep 'Hoplias_malabaricus' $PAIRFILE | grep 'X1' >> ${PAIRFILE_TRBL}



























# test script
```
library(data.table)
library(ggplot2)
library(patchwork)
library(ape)
library(dplyr)
library(tidyr)
library(purrr)


out <- read.csv("/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Anolis_sagrei/sexshared/combined_phaseblks/Anolis_sagrei_syntentic_chromosomes.sex_shared.Gallus_gallus_REF.csv")

tree_file <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/roadies_v1.1.16b.numbers.scientific.nwk"

tr_full <- ape::read.tree(tree_file)

tip_key <- tr_full$tip.label

tr_ultra <- chronos(tr_full, lambda = 1)


# out: your tibble with Species and list-cols
# tr_ultra: your ultrametric tree

# union the three list-cols

out_long <- out %>%
  mutate(
    gg_all = map2(
      map2(GgChr_both, GgChr_chr1_only, ~ union(.x, .y)),
      GgChr_chr2_only,
      ~ union(.x, .y)
    )
  ) %>%
  select(Species, gg_all) %>%
  unnest(gg_all) %>%
  filter(!is.na(gg_all), gg_all != "") %>%
  distinct(Species, gg_all) %>%
  rename(GgChr = gg_all) %>%
  # --- NEW: split combos like "22,Z" into separate rows "22" and "Z"
  separate_rows(GgChr, sep = ",") %>%
  mutate(GgChr = str_trim(GgChr)) %>%
  filter(GgChr != "") %>%
  distinct(Species, GgChr)

# order chicken chromosomes: numeric first, then others (Z/W/etc)
chr_levels <- out_long %>%
  distinct(GgChr) %>%
  mutate(
    is_num = suppressWarnings(!is.na(as.integer(GgChr))),
    numval = suppressWarnings(as.integer(GgChr))
  ) %>%
  arrange(desc(is_num), numval, GgChr) %>%
  pull(GgChr)

# full grid for all tips x all chrs
grid_df <- tidyr::expand_grid(
  Species = out_long$Species,
  GgChr   = chr_levels
) %>%
  left_join(out_long %>% mutate(present = 1L), by = c("Species", "GgChr")) %>%
  mutate(present = if_else(is.na(present), 0L, present))





# Plot tree

pdf("tree.pdf", width = 3, height = 10)

# give more right margin for the dot grid + top margin for labels
par(mar = c(2, 2, 6, 12))

# 1) Compute how much horizontal space you need for the dot grid
#    (dx * nchr) plus some padding.
tmp <- plot(tr_ultra, plot = FALSE)
nchr <- length(chr_levels)

# spacing between chromosome columns
dx <- max(tmp$x.lim) * 0.02

# extra horizontal room to the right of the tree (in tree x-units)
grid_width <- dx * (nchr + 6)

# 2) Plot tree but "squish it left" by expanding xlim to the right
plot(tr_ultra,
     cex = 0.1,
     edge.width = 0.25,
     no.margin = TRUE,
     x.lim = c(0, max(tmp$x.lim) + grid_width))

lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_y <- lp$yy[1:Ntip(tr_ultra)]
tips  <- tr_ultra$tip.label

# start x position for the grid (a bit to the right of the tree)
x0 <- max(tmp$x.lim) + dx * 3

# full tips/y from the plotted tree
lp   <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tips <- tr_ultra$tip.label
tip_y <- lp$yy[1:Ntip(tr_ultra)]

# subset: tips that appear in grid_df
tips_with_data <- intersect(tips, unique(grid_df$Species))
idx_data <- match(tips_with_data, tips)
y_data <- tip_y[idx_data]

# within each chromosome column:
for (j in seq_along(chr_levels)) {
  chrj <- chr_levels[j]
  xj <- x0 + (j - 1) * dx

  present_species <- grid_df %>%
    dplyr::filter(GgChr == chrj, present == 1L) %>%
    dplyr::pull(Species)

  is_present_data <- tips_with_data %in% present_species

  # light-gray circles ONLY for species present in grid_df
  points(rep(xj, length(y_data)), y_data,
         pch = 21, cex = 0.05, col = "lightgray", bg = "white")

  # black filled circles for present==1 among those species
  points(rep(xj, sum(is_present_data)), y_data[is_present_data],
         pch = 21, cex = 0.05, bg = "black")
}

# 2) Column labels at the top
usr <- par("usr")

# Put labels slightly ABOVE the top of the plotting region
y_lab <- 585
text(
  x = x0 + (seq_along(chr_levels) - 1) * dx,
  y = y_lab,
  labels = chr_levels,
  srt = 90,                 # rotate
  adj = c(0, 0.5),          # anchor at bottom of rotated text
  xpd = NA,                 # allow drawing into margin area
  cex = 0.2                 # tweak as needed
)

dev.off()
