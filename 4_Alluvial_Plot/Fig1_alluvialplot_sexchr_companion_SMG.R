## Alluvial Plot Summary Fig. 1 Companion Paper
setwd("/Users/jacksondan/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/Documents/VGP_Phase1/Alluvial")

# Load (run every session you need them)
invisible(lapply(packages, library, character.only = TRUE))

library(ggalluvial)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(dplyr)
library(viridis)
library(readr)
library(purrr)
library(data.table)
library(patchwork)
library(cowplot)
library(gridExtra)
library(svglite)

VGP_sex_df_curated <- as.data.frame(read.csv("VGP_list_sex_chroms_curated.ForAlluvial.Csv", header = TRUE, sep = "\t"))

VGP_sex_df.plot <- VGP_sex_df_curated %>%
  filter(Sex.chromosomes.main.haploptype != "") %>%
  filter(Accession_num_main_haplotype != "") %>%
  ## There shouldn't be any blanks in the sex column anymore, but just in case I rename them explicitly
  dplyr::mutate(
    sex = case_when(
      sex == "" | is.na(sex)           ~ "Blank",
      TRUE                             ~ sex
    )
  )

sex_label_mapping <- c(
  "M" = "Male",
  "F" = "Female",
  "UN" = "Unknown",
  "Blank" = "Blank",
  "H" = "Hermaphrodite"
)

# 
VGP_sex_df.plot_clean <- VGP_sex_df.plot %>%
  dplyr::mutate(sex = dplyr::recode_factor(
    sex, "M" = "Male", "F" = "Female", "UN" = "Unknown", "Blank" = "Blank", "H" = "Hermaphrodite"
  ))

#view(VGP_sex_df.plot_clean)

length(VGP_sex_df.plot_clean %>% filter(Sex.chromosomes.main.haploptype.sex == "No", 
                                        Environmental %in% c("No", "UN")))
no_sex_chr_ID_VGP_sex_df.plot_clean <- VGP_sex_df.plot_clean %>% filter(Sex.chromosomes.main.haploptype.sex == "No", 
                                                                        Environmental %in% c("No", "UN"), Accession_num_main_haplotype != "", sex != "Hermaphrodite")

# sex count 
VGP_sex_df.plot_sexCount <- VGP_sex_df.plot_clean %>%
  dplyr::count(sex) %>%
  dplyr::mutate(percentage = n / sum(n) * 100)
VGP_sex_df.plot_sexCount
sum(VGP_sex_df.plot_sexCount$n)

## % species with heterogametic sex sequenced

heterogam_label_mapping <- c(
  "Y" = "Yes",
  "N" = "No",
  "UN" = "Unknown",
  "NA" = "Not Applicable"
)
# Replace NA values with "NA" and map the labels
VGP_sex_df.plot_heterog <- VGP_sex_df.plot %>%
  dplyr::mutate(Hetergametic_Sequenced = ifelse(is.na(Hetergametic_Sequenced), "NA", Hetergametic_Sequenced)) %>%
  dplyr::mutate(Hetergametic_Sequenced = factor(Hetergametic_Sequenced, levels = names(heterogam_label_mapping), labels = heterogam_label_mapping))

heterogametic_data <- VGP_sex_df.plot_heterog %>%
  dplyr::count(Hetergametic_Sequenced) %>%
  dplyr::mutate(percentage = n / sum(n) * 100)

# % with both sex chroms identified

category_order <- c("Yes", "X only (Female)", "Z only (Male)", 
                    "Y is missing (Male)", "W is missing (Female)", "X only (Sex Unknown)",
                    "Z is missing (Female)",
                    "None found", "Unknown SD",
                    "Needs chr named", "Needs curation", "Hermaphrodite", "ESD")

unique(VGP_sex_df.plot$Sex.chromosomes.main.haploptype)
unique(VGP_sex_df.plot$sex)
unique(VGP_sex_df.plot$Environmental)

# Categorize rows
library(dplyr)
library(stringr)

# --- Define label groups in one place ----------------------------------------

yes_vals <- c(
  "Has X and Y", "Has Z and W", "Has 5X and 4Y", "Has 5X and 5Y",
  "Has X and Y1, Y2", "Has X and Y, Y2", "Has X1, X2, and Y",
  "Has Z1, Z2, and W",
  "has Z and W", "Has Z and W (remove identical ones from broiler)", "Has sex chroms",
  "Has X1, X2, Y1, Y2", "Has Y1, Y2",
  "Has sex chromosomes"
)

x_only_female_vals <- c("Has X only")

x_only_unknown_vals <- c("Has X", "Has X (sex unknown)", "Has X, sex unknown")

z_only_male_vals <- c("Has Z only")

y_missing_male_vals <- c(
  "Has X (and no Y, but a male)", "Has X (need to move Y from pat)",
  "Has X (where is Y? male)", "Has X (Y not labelled?)",
  "Has X only (move Y from pat)",
  "Has X only (create new assembly combo and move Y from pat)",
  "Has X (need Y labelled)", "Has X (need Y named?)", "Has X (need Y named)",
  "Has X (need to move Y into maternal haplotype)", "Has X (missing y?)",
  "Has X only (Y missing)"
)

w_missing_female_vals <- c(
  "Has Z only (where is W, female)", "Has Z only (Where is W)",
  "Has Z only (W missing)", "Has Z only (where is W, male)"
)

z_missing_female_vals <- c("Has W (Chr9)")

needs_chr_named_vals <- c(
  "(needs chr named)", "(needs sex chr named)", "(chromosomes need naming)", "(sex chromosomes not named)",
  "Needs sex chr named", "Need sex chr named", "(need chr named)",
  "(need sex chr named)", "Needs curation GCA_026018925.1 high priority", "Needs curation GCA_026229955.1"
)

needs_curation_vals <- c("needs curation", "Needs curation", "(needs curation)")

unknown_sd_vals <- c("No", "NO", "PM_G_Scaffold labels are germline chromosomes")

esd_environment_vals <- c("TSD II", "ESD_other", "size", "TSD Ia", "TSD")

none_found_vals <- c("No", "NO")

# --- Categorize rows ---------------------------------------------------------

VGP_sex_df.plot_chrs.found <- VGP_sex_df.plot %>%
  mutate(
    sch = `Sex.chromosomes.main.haploptype`,
    Category = case_when(
      # Yes categories
      sch %in% yes_vals ~ "Yes",

      # X only categories
      sch %in% x_only_female_vals  ~ "X only (Female)",
      sch %in% x_only_unknown_vals ~ "X only (Sex Unknown)",

      # Z only categories
      sch %in% z_only_male_vals ~ "Z only (Male)",

      # Missing sex chromosome categories
      sch %in% y_missing_male_vals   ~ "Y is missing (Male)",
      sch %in% w_missing_female_vals ~ "W is missing (Female)",
      sch %in% z_missing_female_vals ~ "Z is missing (Female)",

      # Needs naming / curation
      sch %in% needs_chr_named_vals ~ "Needs chr named",
      sch %in% needs_curation_vals  ~ "Needs curation",

      # Split "No"
      sex == "Hermaphrodite" ~ "Hermaphrodite",
      sch == "No" & Environmental %in% esd_environment_vals ~ "ESD",
      sch %in% unknown_sd_vals & Genotypic == "UN" ~ "Unknown SD",
      sch %in% none_found_vals & Genotypic %in% c("male heterogametic", "female heterogametic") ~ "None found",

      TRUE ~ NA_character_
    ),
    Category = factor(Category, levels = category_order)
  ) %>%
  select(-sch)  # remove helper column


# Filter out rows with NA categories
#VGP_sex_df.plot <- VGP_sex_df.plot %>% filter(!is.na(Category))
VGP_sex_df.plot_chrs.found
homogametic_VGP_sex_df.plot_clean <- VGP_sex_df.plot_chrs.found %>% 
  filter(Category %in% c("X only (Female)", "X only (Sex Unknown)", "Y is missing (Male)", "Z only (Male)", "W is missing (Female)"))
# Count and calculate percentages
haplotype_data <- VGP_sex_df.plot_chrs.found %>%
  dplyr::count(Category) %>%
  dplyr::mutate(percentage = n / sum(n) * 100)

##### New sex pie chart with heterogamety:

# Clean sex column
sex_heterogam_df <- VGP_sex_df.plot_chrs.found %>%
  mutate(
    sex_clean = case_when(
      sex == "M" ~ "Male",
      sex == "F" ~ "Female",
      sex == "H" ~ "Hermaphrodite",
      is.na(sex) | sex %in% c("", "UN") ~ "Unknown",
      TRUE ~ sex
    ),
    Genotypic = tolower(Genotypic),
    Environmental = tolower(Environmental)
  )

# Create new pie chart grouping
sex_heterogam_df_grouped <- sex_heterogam_df %>%
  mutate(SexGroup = case_when(
    # Male heterogamety
    Genotypic == "male heterogametic" & sex_clean %in% c("Male", "Female","Unknown") ~ paste0("Male heterogametic: ", sex_clean),
    # Female heterogamety
    Genotypic == "female heterogametic" & sex_clean %in% c("Male", "Female","Unknown") ~ paste0("Female heterogametic: ", sex_clean),
    # ESD
    Environmental %in% c("tsd", "tsd ia", "tsd ib", "tsd ii", "esd_other", "size", "density", "ph") ~ "ESD",
    # Hermaphrodites
    sex_clean == "Hermaphrodite" ~ "Hermaphrodite",
    # Unknown heterogamety
    sex_clean %in% c("Male", "Female", "Unknown") ~ paste0("Unknown SD: ", sex_clean),
    TRUE ~ "Unknown"
  ))
#view(sex_heterogam_df_grouped)
# Step 3: Count and reorder based on your desired levels
custom_order <- c(
  "Male heterogametic: Male", "Male heterogametic: Female", "Male heterogametic: Unknown",
  "Female heterogametic: Male", "Female heterogametic: Female", "Female heterogametic: Unknown",
  "Unknown SD: Male", "Unknown SD: Female", "Unknown SD: Unknown",
  "ESD", "Hermaphrodite"
)

sexgroup_counts <- sex_heterogam_df_grouped %>%
  count(SexGroup) %>%
  filter(SexGroup %in% custom_order) %>%
  mutate(
    SexGroup = factor(SexGroup, levels = custom_order),
    percentage = round(n / sum(n) * 100, 1)
  )
sum(sexgroup_counts$n)

gsd_counts <- sex_heterogam_df_grouped %>%
  count(Environmental == "no")
gsd_counts
# Define custom colors (from viridis palette)
custom_colors_sex <- c(
  "Male heterogametic: Male" = "#440154",     # light purple
  "Male heterogametic: Female" =  "#482173",   # dark purple
  "Male heterogametic: Unknown" = "#433e85", 
  "Female heterogametic: Male" = "#86d549",   # light blue
  "Female heterogametic: Female" = "#c2df23", # dark blue
  "Female heterogametic: Unknown" = "#fde725",
  "Unknown SD: Male" = "#2d708e",  # light green
  "Unknown SD: Female" = "#25858e",# dark green
  "Unknown SD: Unknown" = "#1e9b8a", # lighter green for unknown
  "ESD" = "darkgrey",                # 
  "Hermaphrodite" = "lightgrey"       # 
)

# Recategorize for the bar chart
category_order_bar <- c(
  "Both sex chromosomes annotated", "Sex-limited chromosome not sequenced", 
  "Sex-limited chromosome missing", "Unidentified", "Unknown sex determination", 
  "Environmental sex determination (ESD)", "Needs curation"
)
# Define custom colors 
custom_colors_barp <- c(
  "Both sex chromosomes annotated" = "#0d0887",     # light purple
  "Sex-limited chromosome not sequenced" =  "#b12a90",   # dark purple
  "Sex-limited chromosome missing" = "#6a00a8", 
  "Unidentified" = "#fca636",   # light blue
  "Unknown sex determination" = "#e16462", # dark blue
  "Environmental sex determination (ESD)" = "#f0f921",                # 
  "Needs curation" = "darkgrey"       # 
)

VGP_sex_df.plot_chrs.found_noHerm <- VGP_sex_df.plot_chrs.found %>% 
  dplyr::filter(!Category %in% c("Hermaphrodite")) %>% dplyr::filter(!Lineage %in% c("Invertebrates"))

sex_alluvial_df <- VGP_sex_df.plot_chrs.found_noHerm %>%
  mutate(
    Category_bar = case_when(
      Category %in% c("Yes", "Has sex chromosomes") ~ "Both sex chromosomes annotated",
      Category %in% c("X only (Female)", "Z only (Male)", "X only (Sex Unknown)") ~ "Sex-limited chromosome not sequenced",
      Category %in% c("Y is missing (Male)", "W is missing (Female)", "Z is missing (Female)") ~ "Sex-limited chromosome missing",
      Category %in% c("Needs curation", "Needs chr named") ~ "Needs curation",
      Category %in% c("Unknown SD") ~ "Unknown sex determination",
      Category %in% c("ESD") ~ "Environmental sex determination (ESD)",
      Category %in% c("None found" ) ~ "Unidentified",
      TRUE ~ as.character(Category)
    ),
    Category_bar = factor(Category_bar, levels = category_order_bar),
    Group = case_when(
      Lineage %in% c("Mammals") ~ "Mammals",
      Lineage %in% c("Birds") ~ "Birds",
      Lineage %in% c("Amphibians") ~ "Amphibians",
      #Lineage %in% c("Invertebrates") ~ "Other\nDeuterostomes",
      Superorder %in% c("Testudines-Cryptodira", "Testudines-Pleurodira", "Crocodylia") ~ "Turtles &\nCrocodilians",
      Superorder %in% c("Actinopterygii-Ray-finned fishes", "Sarcopterygii-Lobe-finned fishes", "Agnatha-Jawless fishes") ~ "Ray-finned,\nLobe-finned Fishes\n& Cyclostomes",
      Superorder %in% c("Chondrichthyes-Cartilaginous fishes") ~ "Cartilaginous\nFishes",
      Superorder %in% c("Squamata") ~ "Lepidosaurs",
      TRUE ~ "Other"
    )
  )

group_order <- c("Mammals", "Birds", "Turtles &\nCrocodilians", "Lepidosaurs",  
                 "Amphibians", "Ray-finned,\nLobe-finned Fishes\n& Cyclostomes", "Cartilaginous\nFishes")

sex_alluvial_df <- sex_alluvial_df %>%
  mutate(Group = factor(Group, levels = group_order))

# Make frequency df for alluvial
alluvial_freq_df <- sex_alluvial_df %>%
  transmute(
    Group = Group,
    sex = na_if(as.character(sex), ""),
    Sexual.system = na_if(as.character(Sexual.system), ""),
    Hetergametic_Sequenced = na_if(as.character(Hetergametic_Sequenced), "ESD/H"),
    Genotypic = na_if(as.character(Genotypic), "ESD/H"),
    Category_bar = na_if(as.character(Category_bar), "")
  ) %>%
  count(Group, sex, Sexual.system, Hetergametic_Sequenced, Genotypic, Category_bar, 
        name = "Freq")

fig1_alluvial <- ggplot(data = alluvial_freq_df,
       aes(axis1 = sex,
           axis2 = Sexual.system,
           axis3 = Hetergametic_Sequenced,
           axis4 = Genotypic,
           axis5 = Category_bar,
           y = Freq)) +
  scale_x_discrete(
    limits = c("Sex", "Sexual system", "Heterogametic sequenced", "Genotypic", "Sex Chr Assembly"),
    expand = c(.05, .05)
  ) +
  scale_fill_viridis_d(option = "viridis")+
  xlab("Metadata") +
  geom_alluvium(aes(fill = Group), alpha = 0.85, width = 1/12, na.rm = F) +
  geom_stratum(width = 1/4, color = "grey40", na.rm = F) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 2.6
  ) +
  theme_minimal() +
  ggtitle("Overview of VGP Data Freeze Sex Chromosome Assemblies")
fig1_alluvial

fig1_alluvial_vertical <- fig1_alluvial + coord_flip()
fig1_alluvial_vertical
