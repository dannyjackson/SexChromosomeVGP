library(dplyr)
library(ggplot2)

# Read in the csv
df_filtered <- read.csv(
  "species_curationmethod_PAR_status.csv"
)

# Remove murid species from analyses
df_filtered <- df_filtered[
  df_filtered$PAR_status_revised != "Biological" |
    is.na(df_filtered$PAR_status_revised),
]

# Recode "T2T" genomes to Y
# T2T was shorthand for PAR assembled + T2T + single curated
df_filtered <- df_filtered %>%
  mutate(
    PAR_status_revised = if_else(
      PAR_status_revised == "T2T",
      "Y",
      PAR_status_revised
    )
  )

# Set factor levels to determine plotting order
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


# Sample size of each curation method
df_sample_size <- df_filtered %>%
  count(curation_method, name = "sample_size")

# Make a stacked bar plot showing the data
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

# Save in various file types
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


# Perform a chi-square test of independence
df_chisq <- df_filtered

# Combine various information categories into either:
# "Y" -- the ancestral PAR was assembled to both sex chromosomes
# "N" -- the ancestral PAR was not assembled

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