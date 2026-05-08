# ------------------------------------------------------------------
# VGP sex chromosome alluvial workflow
# 1) Subset to necessary columns
# 2) Recode relevant columns
# 3) Save essential/reformatted dataframe as CSV
# 4) Plot as an alluvial plot
# ------------------------------------------------------------------

setwd("/Users/jacksondan/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/Documents/Alluvial")

# ---- Packages ----
packages <- c("tidyverse", "ggalluvial", "viridis", "readr", "forcats")

install.packages(setdiff(packages, rownames(installed.packages())))
invisible(lapply(packages, library, character.only = TRUE))

# ---- 0) Read data ----
infile <- "VGP_OrdinalList_Phase1Freeze_v1.2_Sept.30.2025_sex_chrs_HalfDeep_SCINKD_v2.tsv"

VGP_raw <- readr::read_tsv(infile, show_col_types = FALSE)

# ---- 1) Subset to necessary columns ----
# Keep ONLY what is used downstream (filters, recodes, grouping, and alluvial axes)
VGP <- VGP_raw %>%
  dplyr::select(
    Accession_num_main_haplotype,
    Sex.chromosomes.main.haploptype,
    Sex.chromosomes.main.haploptype.sex,
    sex,
    Sexual.system,
    Hetergametic_Sequenced,
    Genotypic,
    Environmental,
    Lineage,
    Superorder
  )

# ---- 2) Recode values of relevant columns ----
# Helper vectors used for categorization
yes_vals <- c(
  "Has X and Y", "Has Z and W", "Has 5X and 4Y", "Has 5X and 5Y",
  "Has X and Y1, Y2", "Has X and Y, Y2", "Has X1, X2, and Y",
  "Has Z1, Z2, and W",
  "has Z and W", "Has Z and W (remove identical ones from broiler)", "Has sex chroms",
  "Has X1, X2, Y1, Y2", "Has Y1, Y2",
  "Has sex chromosomes"
)

x_only_female_vals  <- c("Has X only")
x_only_unknown_vals <- c("Has X", "Has X (sex unknown)", "Has X, sex unknown")
z_only_male_vals    <- c("Has Z only")

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
unknown_sd_vals     <- c("No", "NO", "PM_G_Scaffold labels are germline chromosomes")
none_found_vals     <- c("No", "NO")
esd_environment_vals <- c("TSD II", "ESD_other", "size", "TSD Ia", "TSD", "TSD Ib", "density", "ph")

# Desired orders
category_order <- c(
  "Yes", "X only (Female)", "Z only (Male)",
  "Y is missing (Male)", "W is missing (Female)", "X only (Sex Unknown)",
  "Z is missing (Female)",
  "None found", "Unknown SD",
  "Needs chr named", "Needs curation", "Hermaphrodite", "ESD"
)

category_order_bar <- c(
  "Both sex chromosomes annotated", "Sex-limited chromosome not sequenced",
  "Sex-limited chromosome missing", "Unidentified", "Unknown sex determination",
  "Environmental sex determination (ESD)", "Needs curation"
)

group_order <- c(
  "Mammals", "Birds", "Turtles &\nCrocodilians", "Lepidosaurs",
  "Amphibians", "Ray-finned,\nLobe-finned Fishes\n& Cyclostomes", "Cartilaginous\nFishes"
)

VGP_clean <- VGP %>%
  # Keep rows needed for analyses/plot
  dplyr::filter(!is.na(Sex.chromosomes.main.haploptype), Sex.chromosomes.main.haploptype != "") %>%
  dplyr::filter(!is.na(Accession_num_main_haplotype), Accession_num_main_haplotype != "") %>%
  # Standardize/canonicalize key columns
  dplyr::mutate(
  sex = dplyr::case_when(
    is.na(sex) | sex == "" ~ "Blank",
    sex == "M"             ~ "Male",
    sex == "F"             ~ "Female",
    sex == "UN"            ~ "Unknown",
    sex == "H"             ~ "Hermaphrodite",
    TRUE                   ~ as.character(sex)
  ),
  sex = factor(sex, levels = c("Male", "Female", "Unknown", "Hermaphrodite", "Blank")),
  Genotypic = tolower(as.character(Genotypic)),
  Environmental = tolower(as.character(Environmental)),
  Sex.chromosomes.main.haploptype = as.character(Sex.chromosomes.main.haploptype)
) %>%
  # Derive Category (sex chromosome annotation status)
  dplyr::mutate(
    Category = dplyr::case_when(
      Sex.chromosomes.main.haploptype %in% yes_vals ~ "Yes",

      Sex.chromosomes.main.haploptype %in% x_only_female_vals  ~ "X only (Female)",
      Sex.chromosomes.main.haploptype %in% x_only_unknown_vals ~ "X only (Sex Unknown)",

      Sex.chromosomes.main.haploptype %in% z_only_male_vals ~ "Z only (Male)",

      Sex.chromosomes.main.haploptype %in% y_missing_male_vals   ~ "Y is missing (Male)",
      Sex.chromosomes.main.haploptype %in% w_missing_female_vals ~ "W is missing (Female)",
      Sex.chromosomes.main.haploptype %in% z_missing_female_vals ~ "Z is missing (Female)",

      Sex.chromosomes.main.haploptype %in% needs_chr_named_vals ~ "Needs chr named",
      Sex.chromosomes.main.haploptype %in% needs_curation_vals  ~ "Needs curation",

      as.character(sex) == "Hermaphrodite" ~ "Hermaphrodite",

      Sex.chromosomes.main.haploptype %in% c("No", "NO") &
        Environmental %in% tolower(esd_environment_vals) ~ "ESD",

      Sex.chromosomes.main.haploptype %in% unknown_sd_vals &
        Genotypic == "un" ~ "Unknown SD",

      Sex.chromosomes.main.haploptype %in% none_found_vals &
        Genotypic %in% c("male heterogametic", "female heterogametic") ~ "None found",

      TRUE ~ NA_character_
    ),
    Category = factor(Category, levels = category_order)
  ) %>%
  # Recategorize to "Category_bar" used on final axis
  dplyr::mutate(
    Category_bar = dplyr::case_when(
      Category %in% c("Yes") ~ "Both sex chromosomes annotated",
      Category %in% c("X only (Female)", "Z only (Male)", "X only (Sex Unknown)") ~ "Sex-limited chromosome not sequenced",
      Category %in% c("Y is missing (Male)", "W is missing (Female)", "Z is missing (Female)") ~ "Sex-limited chromosome missing",
      Category %in% c("Needs curation", "Needs chr named") ~ "Needs curation",
      Category %in% c("Unknown SD") ~ "Unknown sex determination",
      Category %in% c("ESD") ~ "Environmental sex determination (ESD)",
      Category %in% c("None found") ~ "Unidentified",
      TRUE ~ NA_character_
    ),
    Category_bar = factor(Category_bar, levels = category_order_bar)
  ) %>%
  # Grouping variable for fill (your original taxonomic bucketing)
  dplyr::mutate(
    Group = dplyr::case_when(
      Lineage %in% c("Mammals") ~ "Mammals",
      Lineage %in% c("Birds") ~ "Birds",
      Lineage %in% c("Amphibians") ~ "Amphibians",
      Superorder %in% c("Testudines-Cryptodira", "Testudines-Pleurodira", "Crocodylia") ~ "Turtles &\nCrocodilians",
      Superorder %in% c("Actinopterygii-Ray-finned fishes", "Sarcopterygii-Lobe-finned fishes", "Agnatha-Jawless fishes") ~
        "Ray-finned,\nLobe-finned Fishes\n& Cyclostomes",
      Superorder %in% c("Chondrichthyes-Cartilaginous fishes") ~ "Cartilaginous\nFishes",
      Superorder %in% c("Squamata") ~ "Lepidosaurs",
      TRUE ~ "Other"
    ),
    Group = factor(Group, levels = group_order)
  ) %>%
  # Optional: drop rows that still can't be categorized for plotting
  dplyr::filter(!is.na(Category_bar), !is.na(Group))

# ---- Optional: apply your prior exclusion rules for the alluvial dataset ----
VGP_alluvial_input <- VGP_clean %>%
  dplyr::filter(as.character(Category) != "Hermaphrodite") %>%
  dplyr::filter(as.character(Lineage) != "Invertebrates")

# ---- 3) Save essential and reformatted dataframe as CSV ----
# Save the cleaned, analysis-ready data (one row per accession)
readr::write_csv(VGP_alluvial_input, "VGP_sex_df_essential_reformatted.csv")

# Also save the frequency table used for ggalluvial (useful for reproducibility)
alluvial_freq_df <- VGP_alluvial_input %>%
  dplyr::transmute(
    Group = Group,
    sex = na_if(as.character(sex), ""),
    Sexual.system = na_if(as.character(Sexual.system), ""),
    Hetergametic_Sequenced = na_if(as.character(Hetergametic_Sequenced), ""),
    Genotypic = na_if(as.character(Genotypic), ""),
    Category_bar = na_if(as.character(Category_bar), "")
  ) %>%
  dplyr::count(Group, sex, Sexual.system, Hetergametic_Sequenced, Genotypic, Category_bar, name = "Freq")

readr::write_csv(alluvial_freq_df, "VGP_alluvial_frequency_table.csv")

group_colors <- c(
  "Mammals" = "#E69F00",
  "Birds" = "#00796B",
  "Crocodiles" = "#009688",
  "Turtles" = "#4DB6AC",
  "Lepidosauria" = "#80CBC4",
  "Amphibians" = "#984EA3",
  "Lobe-finned fishes" = "#A6761D",
  "Ray-finned fishes" = "#56B4E9",
  "Cartilaginous fishes" = "#0072B2",
  "Cyclostomes" = "#CC79A7",
  "Other Deurostomes" = "#999999"
)


# ---- 4) Plot the data as an alluvial plot ----
fig1_alluvial <- ggplot(
  data = alluvial_freq_df,
  aes(
    axis1 = sex,
    axis2 = Sexual.system,
    axis3 = Hetergametic_Sequenced,
    axis4 = Genotypic,
    axis5 = Category_bar,
    y = Freq
  )
) +
  scale_x_discrete(
    limits = c("Sex", "Sexual system", "Heterogametic sequenced", "Genotypic", "Sex Chr Assembly"),
    expand = c(.05, .05)
  ) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  xlab("Metadata") +
  geom_alluvium(aes(fill = Group), alpha = 0.85, width = 1/12, na.rm = FALSE) +
  geom_stratum(width = 1/4, color = "grey40", na.rm = FALSE) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 2.6
  ) +
  theme_minimal() +
  ggtitle("Overview of VGP Data Freeze Sex Chromosome Assemblies")

fig1_alluvial

# Vertical version (optional)
fig1_alluvial_vertical <- fig1_alluvial + coord_flip()
fig1_alluvial_vertical

> colnames(alluvial_freq_df)
[1] "Group"                  "sex"                   
[3] "Sexual.system"          "Hetergametic_Sequenced"
[5] "Genotypic"              "Category_bar"          
[7] "Freq"  

# Is the data useful?
1. All samples spatially separated by "group" (alluvial_freq_df$group)
2. Genetic sex determination system y/n (alluvial_freq_df$Genetic_SD)
4. Both sex chromosomes sequenced? (alluvial_freq_df$Hetergametic_Sequenced)

Recode alluvial_freq_df$Category_bar as a new column:
Column: Genetic_SD
"Needs curation"     > U                   
"Sex-limited chromosome not sequenced" > Y
"Both sex chromosomes annotated"     > Y
"Sex-limited chromosome missing"       > Y
"Environmental sex determination (ESD)" > ESD
"Unknown sex determination"  > U         
"Unidentified" > U



alluvial_freq_df <- alluvial_freq_df %>%
  dplyr::mutate(
    Genetic_SD = dplyr::case_when(
      Category_bar %in% c(
        "Both sex chromosomes annotated",
        "Sex-limited chromosome not sequenced",
        "Sex-limited chromosome missing"
      ) ~ "Y",
      Category_bar == "Environmental sex determination (ESD)" ~ "ESD",
      Category_bar %in% c(
        "Needs curation",
        "Unknown sex determination",
        "Unidentified"
      ) ~ "U",
      TRUE ~ NA_character_
    ),
    Genetic_SD = factor(Genetic_SD, levels = c("Y", "ESD", "U"))
  )

# If alluvial_freq_df$Genetic_SD is ESD or U, then convert alluvial_freq_df$Hetergametic_Sequenced to NA
alluvial_freq_df <- alluvial_freq_df %>%
  dplyr::mutate(
    Hetergametic_Sequenced = dplyr::if_else(
      Genetic_SD %in% c("ESD", "U"),
      "UN",
      as.character(Hetergametic_Sequenced)
    )
  )



plot_df <- alluvial_freq_df %>%
  dplyr::mutate(
    Genetic_SD = factor(Genetic_SD, levels = c("Y", "U", "ESD")),
    Hetergametic_Sequenced = factor(Hetergametic_Sequenced, levels = c("Y", "N", "UN")),
    Group = factor(Group)  # set levels too if you want Group order fixed
  )

fig <- ggplot(
  plot_df,
  aes(
    axis1 = Group,
    axis2 = Genetic_SD,
    axis3 = Hetergametic_Sequenced,
    y = Freq
  )
) +
  scale_x_discrete(
    limits = c("Group", "Genetic SD", "Heterogametic sequenced"),
    expand = c(.05, .05)
  ) +
  geom_alluvium(aes(fill = Group), alpha = 0.85, width = 1/12,
                na.rm = FALSE, decreasing = FALSE) +
  geom_stratum(width = 1/16, color = "grey40",
               na.rm = FALSE, decreasing = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2) +
  theme_minimal() 




ggsave("VGP_alluvial_plot_horizontal.pdf", fig, width = 12, height = 6, dpi = 300)
























alluvial_freq_df <- alluvial_freq_df %>%
  dplyr::mutate(
    Genetic_SD = factor(Genetic_SD, levels = c("Y", "U", "ESD"))
  )

alluvial_freq_df <- alluvial_freq_df %>%
  dplyr::mutate(
    Hetergametic_Sequenced = factor(Hetergametic_Sequenced, levels = c("Y", "N", "UN"))
  )











# Put Group on the TOP by making it the LAST axis, then flip
fig_alluvial_vertical_group_on_top <- ggplot(
  alluvial_freq_df,
  aes(
    axis1 = Group,
    axis2 = Genetic_SD,
    axis3 = Hetergametic_Sequenced,
    y = Freq
  )
) +
  scale_x_discrete(
    limits = c("Group", "Genetic SD", "Heterogametic sequenced"),
    expand = c(.05, .05)
  ) +
  geom_alluvium(aes(fill = Group), alpha = 0.85, width = 1/12, na.rm = FALSE) +
  geom_stratum(width = 1/16, color = "grey40", na.rm = FALSE) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 2
  ) +
  theme_minimal() +
  labs(
    title = "Genetic SD → Heterogametic Sequenced → Group",
    x = NULL,
    y = "Count"
  )

fig_alluvial_vertical_group_on_top


5. Is there a par?


# Save plots (optional)
ggsave("VGP_alluvial_plot_horizontal.pdf", fig_alluvial_vertical_group_on_top, width = 12, height = 6, dpi = 300)
ggsave("VGP_alluvial_plot_vertical.png", fig1_alluvial_vertical, width = 8, height = 10, dpi = 300)
