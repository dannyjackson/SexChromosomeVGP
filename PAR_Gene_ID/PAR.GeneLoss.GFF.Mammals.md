# Identify all genes found in any avian PAR, curate fastas for blast analysis, then identify genes within PARs of all genomes
## 0. Quantify gappiness of each PAR
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/mammals

## Excluding:
Microtus_pennsylvanicus,CHROM:0-116367
Ochotona_princeps,CHROM:107724012-107801096 



# Confirm that accessions in sexchrfile still match correctly (genomes have been updated)
SPECIES=Urocitellus_parryii
grep '>' ${BASE}/${SPECIES}/${SPECIES}.fna | grep 'chromosome'

# Revise those that need revision
SEXCHR_FILE="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv"

grep 'Urocitellus_parryii' $SEXCHR_FILE
sed -i 's/CM099876/NC_135547/g' $SEXCHR_FILE
sed -i 's/CM099877/NC_135548/g' $SEXCHR_FILE

Capra_hircus # no matches at all


cat > PAR.species_chr_region.txt <<'EOF'
Eubalaena_glacialis,NC_083736.1:0-7069966
Macaca_nemestrina,NC_092145.1:158149255-159757195 
Callithrix_jacchus,NC_133525.1.1:0-1757279
Mesoplodon_bidens,OZ073217.1:135117924-142816029
Inia_geoffrensis,CM070920.1:0-7142434
Camelus_dromedarius,NC_087472.1:108540129-114202744
Ovis_canadensis,NC_091727.1:0-7072606
Manis_pentadactyla,NC_080038.1:0-4881534
Pan_paniscus,NC_073272.2:0-2524164                      
Marmota_flaviventris,NC_092518.1:521224-10347709        
Balaenoptera_physalus,OZ239531.1:0-7257636
Pan_troglodytes,NC_072421.2:0-3170188          
Rhynchonycteris_naso,CM073052.1:0-10832599
Ovis_aries,CP162266.1:0-7021881
Trichechus_inunguis,CM102173.1:0-9558467
Pseudorca_crassidens,NC_090317.1:127639906-136059651 
Myotis_nattereri,OZ125678.2:0-1449976
Rhynchocyon_petersi,CM091802.1:0-20182068
Grampus_griseus,OZ206318.1:0-7148355
Mustela_nivalis_vulgaris,OZ211688.1:127035994-138477583
Capra_hircus,CP168640.1:0-7040093
Loxodonta_africana,NC_087369.1:0-10314587
Urocitellus_parryii,NC_135547.1:122998411-131433435 
Meles_meles,NC_060087.1:0-6387864
Homo_sapiens,NC_060947.1:0-2394410
EOF

chmod +x find_PAR_gaps.sh 
./find_PAR_gaps.sh

Rscript plot_PAR_gaps.R
```
### Species without gaps
Loxodonta africana has four gaps each of 13 base pairs... this feels a little suspicious to me? Since assemblers can put in any random number of bp for a gap? Excluding it for now.
```
cat > species_par.csv <<'EOF'
Pseudorca_crassidens,127639906-136059651 
Inia_geoffrensis,0-7142434
Ovis_canadensis,0-7072606
Capra_hircus,0-7040093
Ovis_aries,0-7021881
Meles_meles,0-6387864
Pan_troglodytes,0-3170188          
Pan_paniscus,0-2524164                      
Homo_sapiens,0-2394410
Callithrix_jacchus,0-1757279
EOF
```
## 1. Curate annotated genes found in the PARs
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/GFFs/mammals

chmod +x extract_x_par_genes.sh
./extract_x_par_genes.sh species_par.csv > x_par_genes.tsv

# reannotate all genes with only "-like" gene names

chmod +x infer_gene_names_from_descriptions.py

./infer_gene_names_from_descriptions.py \
  -i x_par_genes.tsv \
  -o x_par_genes.with_inferred_names.tsv \
  -f x_par_genes.failed_to_rename.tsv

```
# Compute freq of loc genes
```
#!/usr/bin/env Rscript

library(dplyr)
library(readr)

# Input file from your AWK script
infile <- "x_par_genes.with_inferred_names.tsv"

df <- read_tsv(infile, show_col_types = FALSE)

df2 <- df %>%
  mutate(
    Region = if_else(PAR_status == "Y", "PAR", "nonPAR"),
    Is_LOC = grepl("^LOC", GeneName),
    Is_Uncharacterized_LOC = grepl("uncharacterized", GeneDescription, ignore.case = TRUE)
  )

loc_freq <- df2 %>%
  group_by(Region) %>%
  summarise(
    TotalGenes = n(),
    LOCGenes = sum(Is_LOC, na.rm = TRUE),
    NonLOCGenes = sum(!Is_LOC, na.rm = TRUE),
    LOCFrequency = LOCGenes / TotalGenes,
    LOCPercent = 100 * LOCFrequency,
    UncharacterizedLOCGenes = sum(Is_Uncharacterized_LOC, na.rm = TRUE),
    UncharacterizedLOCFrequency = UncharacterizedLOCGenes / TotalGenes,
    UncharacterizedLOCPercent = 100 * UncharacterizedLOCFrequency,
    .groups = "drop"
  )

print(loc_freq, width = Inf)

write_tsv(loc_freq, "loc_gene_frequency_PAR_vs_nonPAR.tsv")
```
















### Note if telomere is ID'd on the end of the PAR
```
telomere_file="gapless_species.telomeres.txt"

cat > "$telomere_file" <<'EOF'
Pseudorca_crassidens,R,YES
Inia_geoffrensis,L,YES
Ovis_canadensis,L,YES
Capra_hircus,L,YES
Ovis_aries,L,YES
Meles_meles,L,YES
Pan_troglodytes,L,
Pan_paniscus,L,YES
Homo_sapiens,L,YES
Callithrix_jacchus,L,YES
EOF
```
## 3. Plot the data
### Plot a PCA of PAR status by gene and species
```
module load R

R

library(tidyverse)
library(ggrepel)

infile <- "gene_locations_by_species.with_chr_label.all.In_PAR.X_only.gapless_species.csv"

df <- read_csv(infile, show_col_types = FALSE)

# Collapse multiple Z hits per Species/Gene.
# Priority:
#   Y/Edge > N
z_states <- df %>%
  mutate(
    state = case_when(
      In_PAR %in% c("Y", "Edge") ~ 2,
      In_PAR == "N" ~ 1,
      In_PAR == "U" ~ 0,
      is.na(In_PAR) ~ 0,
      TRUE ~ 0
    )
  ) %>%
  group_by(Species, Gene) %>%
  summarise(
    state = max(state),
    .groups = "drop"
  )

# Make the full Species x Gene matrix.
# Missing Species/Gene combinations are interpreted as not found on Z.
mat <- z_states %>%
  pivot_wider(
    names_from = Gene,
    values_from = state,
    values_fill = 0
  )

species <- mat$Species

pca_mat <- mat %>%
  select(-Species) %>%
  as.data.frame()

rownames(pca_mat) <- species

# Remove genes with no variation across species
pca_mat <- pca_mat[, apply(pca_mat, 2, var) > 0]

pca <- prcomp(pca_mat, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("Species")

var_explained <- summary(pca)$importance[2, ] * 100

ggplot(pca_df, aes(PC1, PC2, label = Species)) +
  geom_point(size = 3) +
  geom_text_repel(
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.1,
    point.padding = 0.3,
    segment.alpha = 0.5
  ) +
  xlab(paste0("PC1 (", round(var_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  theme_bw()

ggsave("In_PAR_PCA_from_X_only.pdf", width = 7, height = 6)

write_csv(
  mat,
  "species_by_gene.In_PAR_numeric_matrix.from_X_only.csv"
)
```
### Plot a PCA of PAR status by gene
```
library(ggrepel)

# mat is Species x Gene:
# first column = Species
# remaining columns = numeric states 0/1/2

gene_pca_mat <- mat %>%
  column_to_rownames("Species") %>%
  t() %>%
  as.data.frame()

# Remove species columns with no variation across genes
gene_pca_mat <- gene_pca_mat[, apply(gene_pca_mat, 2, var) > 0]

# Remove genes with no variation across species after transpose, if any
gene_pca_mat <- gene_pca_mat[apply(gene_pca_mat, 1, var) > 0, ]

gene_pca <- prcomp(gene_pca_mat, center = TRUE, scale. = TRUE)

gene_pca_df <- as.data.frame(gene_pca$x) %>%
  rownames_to_column("Gene")

gene_var_explained <- summary(gene_pca)$importance[2, ] * 100

ggplot(gene_pca_df, aes(PC1, PC2)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = Gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.alpha = 0.5
  ) +
  xlab(paste0("PC1 (", round(gene_var_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(gene_var_explained[2], 1), "%)")) +
  theme_bw()

ggsave("In_PAR_PCA_genes_labeled.pdf", width = 10, height = 8)

write_csv(gene_pca_df, "In_PAR_PCA_gene_coordinates.csv")
```
# Make combined plot
```
Rscript ../PAR_Combined_Plot.R mammals
```
# Permissions modification
chmod -R g+rwx /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR*
chmod -R g+rwx /data/Wilson_Lab/data/VGP_genomes_phase1/
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles