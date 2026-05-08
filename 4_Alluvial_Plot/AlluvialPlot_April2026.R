library(ggalluvial)
library(ggplot2)
library(tidyverse)
library(readr)
library(dplyr)
library(viridis)

# Read file
df <- read_csv("VGP_list_sex_chroms_curated.ForAlluvial.csv", show_col_types = FALSE)

lineage_colors <- c(
  "Mammals" = "#E69F00",              # Okabe-Ito Orange
  "Birds" = "#00796B",                # Teal 700
  "Crocodiles" = "#009688",           # Teal 500
  "Turtles" = "#4DB6AC",              # Teal 300
  "Lepidosauria" = "#80CBC4",         # Teal 200
  "Amphibians" = "#984EA3",           # Purple
  "Lobe-finned fishes" = "#A6761D",   # Brown
  "Ray-finned fishes" = "#56B4E9",    # Okabe-Ito Sky Blue
  "Cartilaginous fishes" = "#0072B2", # Okabe-Ito Blue
  "Cyclostomes" = "#CC79A7",          # Okabe-Ito Magenta
  "Other Deurostomes" = "#999999"     # Neutral Gray
)

# Clean and standardize values
plot_df <- df %>%
  mutate(
    Extended_lineage = if_else(
      is.na(Extended_lineage) | Extended_lineage == "",
      "Unknown lineage",
      Extended_lineage
    ),

    GSD = case_when(
      is.na(GSD) | GSD == "" ~ "Unknown",
      GSD == "Y" ~ "Yes",
      GSD == "N" ~ "No",
      GSD == "UN" ~ "Unknown",
      TRUE ~ GSD
    ),

    XZ_Idd = case_when(
      is.na(XZ_Idd) | XZ_Idd == "" ~ "Unknown",
      XZ_Idd == "Y" ~ "Yes",
      XZ_Idd == "N" ~ "No",
      XZ_Idd == "UN" ~ "Unknown",
      TRUE ~ XZ_Idd
    ),

    Heterogametic_sequenced = case_when(
      is.na(Heterogametic_sequenced) | Heterogametic_sequenced == "" ~ "Unknown",
      Heterogametic_sequenced == "Y" ~ "Yes",
      Heterogametic_sequenced == "N" ~ "No",
      Heterogametic_sequenced == "UN" ~ "Unknown",
      TRUE ~ Heterogametic_sequenced
    ),

    YW_Idd = case_when(
      is.na(YW_Idd) | YW_Idd == "" ~ "Unknown",
      YW_Idd == "Y" ~ "Yes",
      YW_Idd == "N" ~ "No",
      YW_Idd == "UN" ~ "Unknown",
      TRUE ~ YW_Idd
    ),

    PAR_assembled = case_when(
      PAR_assembled == "Y" ~ "Yes",
      PAR_assembled == "N" ~ "No",
      PAR_assembled == "NA" ~ "NA",
      TRUE ~ PAR_assembled
    ),

    Total_sampled = Extended_lineage
  )

# Optional: set axis orders explicitly
plot_df <- plot_df %>%
  mutate(
    Total_sampled = factor(
      Total_sampled, 
      levels = c("Other Deurostomes",
        "Cyclostomes",
        "Cartilaginous fishes",
        "Ray-finned fishes",
        "Lobe-finned fishes",
        "Amphibians",
        "Mammals",
        "Lepidosauria",
        "Turtles",
        "Crocodiles",
        "Birds")
    ),
    GSD = factor(
      GSD,
      levels = c("Yes", "No", "Unknown", "NA")
    ),
    XZ_Idd = factor(
      XZ_Idd,
      levels = c("Yes", "No", "Unknown", "NA")
    ),
    Heterogametic_sequenced = factor(
      Heterogametic_sequenced,
      levels = c("Yes", "No", "Unknown", "NA")
    ),
    YW_Idd = factor(
      YW_Idd,
      levels = c("Yes", "No", "Unknown", "NA")
    ),
    PAR_assembled = case_when(
      is.na(PAR_assembled) ~ "NA",
      TRUE ~ as.character(PAR_assembled)
    ),
    PAR_assembled = factor(
      PAR_assembled,
      levels = c("Yes", "No", "Unknown", "NA")
    )
  )

lineage_levels <- c(
  "Other Deurostomes",
  "Cyclostomes",
  "Cartilaginous fishes",
  "Ray-finned fishes",
  "Lobe-finned fishes",
  "Amphibians",
  "Mammals",
  "Lepidosauria",
  "Turtles",
  "Crocodiles",
  "Birds"
)

plot_df <- plot_df %>%
  mutate(
    Extended_lineage = factor(Extended_lineage, levels = lineage_levels),
    lineage_order = match(Extended_lineage, lineage_levels)
  )

# Build frequency table for alluvial
alluvial_freq_df <- plot_df %>%
  count(
    Extended_lineage,
    Total_sampled,
    GSD,
    XZ_Idd,
    Heterogametic_sequenced,
    YW_Idd,
    PAR_assembled,
    name = "Freq"
  )

# Alluvial plot
fig_alluvial <- ggplot(
  alluvial_freq_df,
  aes(
    axis1 = Total_sampled,
    axis2 = GSD,
    axis3 = XZ_Idd,
    axis4 = Heterogametic_sequenced,
    axis5 = YW_Idd,
    axis6 = PAR_assembled,
    y = Freq
  )
) +
  scale_x_discrete(
    limits = c(
      "Total sampled",
      "GSD",
      "XZ identified",
      "Heterogametic sequenced",
      "YW identified",
      "PAR assembled"
    ),
    expand = c(.06, .06)
  ) +
  geom_alluvium(
    aes(fill = Extended_lineage),
    width = 1/10,
    alpha = 0.9,
    aes.bind = "alluvia"
  ) +
  geom_stratum(
    width = 1/4,
    fill = "grey90",
    color = "grey40"
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3
  ) +
  scale_fill_manual(values = lineage_colors, drop = FALSE) +
  labs(
    title = "Overview of sampled species, sex chromosome identification, and PAR assembly",
    x = NULL,
    y = "Number of species",
    fill = "Extended lineage"
  ) +
  theme_minimal(base_size = 12)

fig_alluvial

output_file <- "Horizontal_alluvial.pdf"
ggsave(
  output_file,
  plot = fig_alluvial,
  width = 16,
  height = 8,
  dpi = 300
)




















# Optional vertical version
fig_alluvial_vertical <- ggplot(
  alluvial_freq_df,
  aes(
    axis1 = PAR_assembled,
    axis2 = Heterogametic_sequenced,
    axis3 = XZ_Idd,
    axis4 = GSD,
    axis5 = Total_sampled,
    y = Freq
  )
) +
  scale_x_discrete(
    limits = c(
      "PAR assembled",
      "Heterogametic sequenced",
      "Sex chr identified",
      "GSD",
      "Total sampled"
    ),
    expand = c(.06, .06)
  ) +
  geom_alluvium(aes(fill = Extended_lineage), width = 1/10, alpha = 0.9) +
  geom_stratum(width = 1/4, fill = "grey90", color = "grey40") +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3
  ) +
  scale_fill_manual(values = lineage_colors, drop = FALSE) +
  coord_flip() +
  labs(
    title = "Overview of sampled species, sex chromosome identification, and PAR assembly",
    x = NULL,
    y = "Number of species",
    fill = "Extended lineage"
  ) +
  theme_minimal(base_size = 12)

fig_alluvial_vertical