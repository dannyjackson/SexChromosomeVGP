#!/usr/bin/env Rscript

# Filter VGP CSV into good / poor / uncertain assembly-tech files.
#
# Rules (case-insensitive):
#   - POOR: any Assembly.tech containing "clr"
#   - POOR: any Assembly.tech that is "hifi solo" (interpreted as any HiFi/Hifiasm + solo)
#           e.g., "HiFi, Hifiasm solo", "HiFi hifiasm, solo", "ONT duplex, Hifiasm solo", etc.
#   - GOOD: everything else with a non-empty Assembly.tech
#   - UNCERTAIN: missing/blank Assembly.tech
#
# Outputs:
#   good:      ..._goodAssemblytech.csv   (exact path you provided)
#   poor:      same name with _poorAssemblytech.csv
#   uncertain: same name with _uncertainAssemblytech.csv

library(readr)
library(dplyr)
library(stringr)
library(tidyr)


infile <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/VGP_OrdinalList_Phase1Freeze_v1.1_sex_chroms_seqCov_HalfDeep_SCINKD_May8.25_DwnldDec16.25.csv"

good_out <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/VGP_OrdinalList_Phase1Freeze_v1.1_goodAssemblytech.csv"
poor_out <- sub("_goodAssemblytech\\.csv$", "_poorAssemblytech.csv", good_out)
uncertain_out <- sub("_goodAssemblytech\\.csv$", "_uncertainAssemblytech.csv", good_out)

colname <- "Assembly.tech"

# Read as all character columns to avoid surprises
df <- read_csv(infile, col_types = cols(.default = col_character()), na = c("", "NA", "NaN"))

if (!(colname %in% names(df))) {
  stop(sprintf('Column "%s" not found. Columns are: %s', colname, paste(names(df), collapse = ", ")))
}

# Normalize for matching
tech_norm <- df[[colname]] %>%
  replace_na("") %>%
  str_replace_all("[\u2013\u2014]", "-") %>%   # en/em dash
  str_squish() %>%
  str_to_lower()

uncertain_mask <- tech_norm == ""

# POOR if contains "clr"
poor_clr <- str_detect(tech_norm, "clr")

# "HiFi solo" exclusion:
# interpret as any line that contains "solo" AND indicates HiFi/Hifiasm.
# This catches:
#  - "HiFi, Hifiasm solo"
#  - "HiFi hifiasm, solo"
#  - "ONT duplex, Hifiasm solo"
poor_hifi_solo <- str_detect(tech_norm, "solo") & (str_detect(tech_norm, "hifi") | str_detect(tech_norm, "hifiasm"))

poor_mask <- (!uncertain_mask) & (poor_clr | poor_hifi_solo)
good_mask <- (!uncertain_mask) & (!poor_mask)

good_df <- df[good_mask, , drop = FALSE]
poor_df <- df[poor_mask, , drop = FALSE]
uncertain_df <- df[uncertain_mask, , drop = FALSE]

# Ensure output directory exists
out_dir <- dirname(good_out)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(good_df, good_out, na = "")
write_csv(poor_df, poor_out, na = "")
write_csv(uncertain_df, uncertain_out, na = "")

cat("Wrote:\n")
cat(sprintf("  GOOD      %6d  -> %s\n", nrow(good_df), good_out))
cat(sprintf("  POOR      %6d  -> %s\n", nrow(poor_df), poor_out))
cat(sprintf("  UNCERTAIN %6d  -> %s\n", nrow(uncertain_df), uncertain_out))
