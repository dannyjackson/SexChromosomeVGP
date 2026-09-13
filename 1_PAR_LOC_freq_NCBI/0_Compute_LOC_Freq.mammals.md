# Identify all genes found in any avian PAR, curate fastas for blast analysis, then identify genes within PARs of all genomes
## 0. Quantify gappiness of each PAR
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Compute_LOC_Freq/mammals

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Compute_LOC_Freq/mammals

```
# Create list of PARs
Dropped Ovis aries and Macaca nemestrina from this analysis because all of the NCBI annotations are LOC or missing.

```
cat > PAR.species_chr_region.ncbi.txt <<'EOF'
Artibeus_lituratus,CM076392.3:0-4918541
Balaenoptera_physalus,OZ239531.1:0-7257636
Bos_taurus,NC_040105.1:0-6843483
Callithrix_jacchus,NC_133524.1:0-1757279
Camelus_dromedarius,NC_087472.1:108540129-114202744
Capra_hircus,CP168640.1:0-7202443
Eubalaena_glacialis,NC_083736.1:0-7069966
Gorilla_gorilla,NC_073247.2:0-12039748
Grampus_griseus,OZ206318.1:0-7148355
Homo_sapiens,NC_060947.1:0-2394410
Inia_geoffrensis,CM070920.1:0-7142434
Loxodonta_africana,NC_087369.1:0-10314587
Lycaon_pictus,CM082710.1:0-6575129
Manis_pentadactyla,NC_080038.1:0-5004179
Marmota_flaviventris,NC_092518.1:0-10556221
Meles_meles,NC_060087.1:0-6399546
Mesoplodon_bidens,OZ073217.1:135117924-142816029
Molossus_nigricans,CM078089.1:0-4370896
Mustela_nivalis_vulgaris,OZ211688.1:131948904-138477583
Myotis_mystacinus,OZ075425.2:125458419-127971689
Myotis_nattereri,OZ125678.2:0-2687171
Ovis_canadensis,NC_091727.1:0-7072606
Pan_paniscus,NC_073272.2:0-2524164
Pan_troglodytes,NC_072421.2:0-3170188
Panthera_onca,CM102116.1:122788602-130846311
Pongo_abelii,NC_072008.2:0-2382235
Pongo_pygmaeus,NC_072396.2:0-2356740
Pseudorca_crassidens,NC_090317.1:127639906-136059651
Rhynchocyon_petersi,CM091802.1:0-20182068
Rhynchonycteris_naso,CM073052.1:138579092-142670400
Symphalangus_syndactylus,NC_072447.2:0-16994066
Trichechus_inunguis,CM102173.1:0-9624845
Urocitellus_parryii,NC_135547.1:122998411-131433435
Artibeus_intermedius,CM076326.2:147359214-151570710
Miniopterus_schreibersii,OZ071095.2:103813303-106815783
Canis_lupus,NC_132876.1:118611595-125236116
Nyctalus_leisleri,OZ183628.2:0-1665649
Dasypus_novemcinctus,NC_080704.1:0-9679476
Corynorhinus_townsendii,CM133721.1:0-1829836
EOF
```
## 1. Curate annotated genes found in the PARs
```
chmod +x 0a_extract_x_par_genes.sh
./0a_extract_x_par_genes.sh PAR.species_chr_region.ncbi.txt > x_par_genes.tsv
```
# Compute freq of loc genes 
```
#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(ggrepel)
library(ggrepel)

# Input file from your AWK script
infile <- "x_par_genes.tsv"

df <- read_tsv(infile, show_col_types = FALSE) %>%
  filter(Species != "Homo_sapiens")
  

df <- read_tsv(infile, show_col_types = FALSE) 


df2 <- df %>%
  mutate(
    Region = if_else(PAR_status == "Y", "PAR", "nonPAR"),

    Is_LOC =
    grepl("^LOC", coalesce(GeneName, ""), ignore.case = TRUE) |
    grepl("-like$", coalesce(GeneName, ""), ignore.case = TRUE) |
    grepl("-like$", coalesce(GeneDescription, ""), ignore.case = TRUE),

    Is_Uncharacterized_LOC =
      Is_LOC &
      grepl(
        "uncharacterized",
        coalesce(GeneDescription, ""),
        ignore.case = TRUE
      )
  )


# LOC frequency by species and region
loc_freq_by_species <- df2 %>%
  group_by(Species, Region) %>%
  summarise(
    TotalGenes = n(),
    LOCGenes = sum(Is_LOC, na.rm = TRUE),
    NonLOCGenes = sum(!Is_LOC, na.rm = TRUE),
    LOCFrequency = LOCGenes / TotalGenes,
    LOCPercent = 100 * LOCFrequency,

    UncharacterizedLOCGenes =
      sum(Is_Uncharacterized_LOC, na.rm = TRUE),

    UncharacterizedLOCFrequency =
      UncharacterizedLOCGenes / TotalGenes,

    UncharacterizedLOCPercent =
      100 * UncharacterizedLOCFrequency,

    .groups = "drop"
  )

print(loc_freq_by_species, n = Inf, width = Inf)

write_tsv(
  loc_freq_by_species,
  "LOC_frequency_by_species_and_region.tsv"
)

# Put PAR and nonPAR LOC percentages side by side
loc_variation <- loc_freq_by_species %>%
  select(
    Species,
    Region,
    TotalGenes,
    LOCGenes,
    LOCPercent
  ) %>%
  pivot_wider(
    names_from = Region,
    values_from = c(TotalGenes, LOCGenes, LOCPercent)
  ) %>%
  mutate(
    LOCPercent_difference_PAR_minus_nonPAR =
      LOCPercent_PAR - LOCPercent_nonPAR
  ) %>%
  arrange(desc(LOCPercent_difference_PAR_minus_nonPAR))

print(loc_variation, n = Inf, width = Inf)

write_tsv(
  loc_variation,
  "LOC_percent_variation_by_species.tsv"
)

plot_df <- plot_df %>%
  mutate(
    Hap2_Y = if_else(
      Species %in% c("Bos_taurus", "Lycaon_pictus"),
      "Y",
      "N"
    )
  )

p <- ggplot(
  plot_df,
  aes(
    x = Region,
    y = LOCPercent,
    group = Species
  )
) +
  geom_line(
    aes(color = Hap2_Y),
    alpha = 0.5
  ) +
  geom_point(
    aes(color = Hap2_Y),
    size = 2.5
  ) +
  geom_text_repel(
    data = plot_df %>% filter(Region == "PAR"),
    aes(
      label = Species,
      color = Hap2_Y
    ),
    nudge_x = 0.15,
    direction = "y",
    hjust = 0,
    segment.alpha = 0.4,
    size = 3,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "N" = "black",
      "Y" = "blue"
    ),
    breaks = "Y",
    labels = "Y in hap2",
    name = NULL
  ) +
  scale_x_discrete(
    limits = c("nonPAR", "PAR"),
    expand = expansion(mult = c(0.1, 0.45))
  ) +
  labs(
    x = NULL,
    y = "% LOC-like genes"
  ) +
  theme_bw()

print(p)

ggsave(
  "LOC_percent_PAR_vs_nonPAR_by_species.labeled.pdf",
  p,
  width = 8,
  height = 7
)

average_difference <- loc_variation %>%
  filter(
    !Species %in% c(
      "Lycaon_pictus",
      "Bos_taurus",
      "Homo_sapiens"
    )
  ) %>%
  summarise(
    n_species = sum(!is.na(LOCPercent_difference_PAR_minus_nonPAR)),
    mean_PAR_minus_nonPAR =
      mean(
        LOCPercent_difference_PAR_minus_nonPAR,
        na.rm = TRUE
      ),
    sd_PAR_minus_nonPAR =
      sd(
        LOCPercent_difference_PAR_minus_nonPAR,
        na.rm = TRUE
      )
  )

print(average_difference)

lycaon_bos_difference <- loc_variation %>%
  filter(Species %in% c("Lycaon_pictus", "Bos_taurus")) %>%
  select(
    Species,
    LOCPercent_nonPAR,
    LOCPercent_PAR,
    LOCPercent_difference_PAR_minus_nonPAR
  )

print(lycaon_bos_difference, n = Inf, width = Inf)