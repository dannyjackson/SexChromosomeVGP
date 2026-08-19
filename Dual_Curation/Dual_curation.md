cd /Users/jacksondan/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/Documents/VGP_Phase1/PARs/Curation_Effects

R

```
library(dplyr)

df <- read.csv("VGPPhase1-haplotype-comparison-complete.csv")

# subset to mammals and birds


# subset to mammals
df_m <- df %>% 
  filter(grepl("m", Order))

# subset to birds
df_b <- df %>% 
  filter(grepl("b", Order))

# combine mammals and birds
df_mb <- bind_rows(df_m, df_b)

df_mb_simple <- unique(df_mb[,c("Scientific.Name","curation_method")])

write.csv(df_mb_simple, file = "species_curationmethod.csv")

```
# Combine with Sex Chromosome data
```
library(dplyr)

df_curation <- read.csv("species_curationmethod.csv")

df_sexchr <- read.csv("VGP_Phase1plus_Freeze_v1.3_sex_chr_supplement.csv")

df_sexchr_mb <- df_sexchr %>% 
  filter(Lineage %in% c("Mammals", "Birds"))


df_sexchr_mb_simple <- unique(df_sexchr_mb[,c("Scientific.Name","Lineage", "sex.v1.3", "Hetergametic_Sequenced")])

write.csv(df_sexchr_mb_simple, "df_sexchr_mb_simple.csv")

df_all <- left_join(df_curation, df_sexchr_mb_simple, by = "Scientific.Name")

write.csv(df_all, "species_curationmethod_genomesex.csv")

```
# Saved and manually edited any instances that didn't join properly
```
library(dplyr)
library(ggplot2)

df <- read.csv("species_curationmethod_genomesex.edited.csv")

df_nodata <- df %>% filter(is.na(Hetergametic_Sequenced))

df_filtered <- df %>% filter(Hetergametic_Sequenced == "Y")

par_species <- c("Aegotheles_albertisi",
"Anas_platyrhynchos",
"Aythya_ferina",
"Aythya_marila",
"Calonectris_borealis",
"Coloeus_monedula",
"Columba_livia",
"Cyanocitta_cristata",
"Cygnus_columbianus",
"Dixiphia_pipra",
"Heliangelus_exortis",
"Larus_argentatus",
"Lathamus_discolor",
"Mergus_octosetaceus",
"Morphnus_guianensis",
"Numenius_arquata",
"Opisthocomus_hoazin",
"Passer_domesticus",
"Patagioenas_fasciata",
"Platalea_leucorodia",
"Phaethon_aethereus",
"Rissa_tridactyla",
"Sarcoramphus_papa",
"Strix_aluco",
"Struthio_camelus_australis",
"Taeniopygia_guttata",
"Zosterops_lateralis",
"Balaenoptera_physalus",
"Bos_taurus",
"Callithrix_jacchus",
"Camelus_dromedarius",
"Capra_hircus",
"Eubalaena_glacialis",
"Gorilla_gorilla",
"Grampus_griseus",
"Homo_sapiens",
"Inia_geoffrensis",
"Loxodonta_africana",
"Lycaon_pictus",
"Macaca_nemestrina",
"Manis_pentadactyla",
"Marmota_flaviventris",
"Meles_meles",
"Mesoplodon_bidens",
"Molossus_nigricans",
"Mustela_nivalis_vulgaris",
"Myotis_nattereri",
"Ovis_aries",
"Ovis_canadensis",
"Pan_paniscus",
"Pan_troglodytes",
"Panthera_onca",
"Pongo_abelii",
"Pongo_pygmaeus",
"Pseudorca_crassidens",
"Rhynchocyon_petersi",
"Rhynchonycteris_naso",
"Symphalangus_syndactylus",
"Trichechus_inunguis",
"Urocitellus_parryii"
)

# Match the space-separated names in Scientific.Name
par_species <- gsub("_", " ", par_species)

df_filtered <- df_filtered %>%
  mutate(
    PAR_status = if_else(
      Scientific.Name %in% par_species,
      "Y",
      "N"
    )
  )

write.csv(df_filtered, "species_curationmethod_genomesex_PARstatus.csv")
```
### R script
```
library(dplyr)
library(ggplot2)

df_filtered <- read.csv(
  "species_curationmethod_genomesex_PARstatus.edited.csv"
)

df_filtered <- df_filtered[
  df_filtered$PAR_status_revised != "Biological" |
    is.na(df_filtered$PAR_status_revised),
]

df_filtered <- df_filtered %>%
  mutate(
    PAR_status_revised = if_else(
      PAR_status_revised == "T2T",
      "Y",
      PAR_status_revised
    )
  )

df_filtered$PAR_status_revised <- factor(
  df_filtered$PAR_status_revised,
  levels = c(
    "TBD",
    "N",
    "Misassembly",
    "Biological (putative)",
    "Y"
  )
)


# Sample size for each bar
df_sample_size <- df_filtered %>%
  count(curation_method, name = "sample_size")

p_proportion <- ggplot(
  df_filtered,
  aes(x = curation_method, fill = PAR_status_revised)
) +
  geom_bar(position = "fill") +
  # Add sample size above each bar
  geom_text(
    data = df_sample_size,
    aes(
      x = curation_method,
      y = 1.03,
      label = paste0("N=", sample_size)
    ),
    inherit.aes = FALSE,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "Biological" = "#858585",
      "Biological (putative)" = "#c7c7c7",
      "Misassembly" = "#008b31",
      "N" = "#3c3c3c",
      "TBD" = "#f9a7a7",
      "Y" = "#762A83"
    ),
    labels = c(
      "Biological" = "Biological loss of PAR (known)",
      "Biological (putative)" = "Putative biological loss of PAR",
      "Misassembly" = "Misassembly of PAR",
      "N" = "PAR Collapse",
      "TBD" = "TBD",
      "Y" = "PAR Assembled"
    ),
    na.value = "grey80"
  ) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.10))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Curation method",
    y = "Proportion",
    fill = "PAR status"
  ) +
  theme_bw() +
  theme( 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank() )


ggsave(
  "proportion.pdf",
  plot = p_proportion,
  width = 5,
  height = 7
)

ggsave(
  "proportion.png",
  plot = p_proportion,
  width = 5,
  height = 7,
  dpi = 300
)

ggsave(
  "proportion.svg",
  plot = p_proportion,
  width = 5,
  height = 7,
  dpi = 300
)


df_chisq <- df_filtered

df_chisq <- df_chisq %>%
  mutate(
    PAR_status_binary = case_when(
      PAR_status_revised %in% c("Misassembly", "N", "Biological (putative)") ~ "N",
      PAR_status_revised %in% c("Y") ~ "Y",
      TRUE ~ NA_character_
    )
  )

# Contingency table
chisq_tab <- table(
  df_chisq$curation_method,
  df_chisq$PAR_status_binary
)

chisq_tab

# Chi-squared test
chisq_result <- chisq.test(chisq_tab)

chisq_result