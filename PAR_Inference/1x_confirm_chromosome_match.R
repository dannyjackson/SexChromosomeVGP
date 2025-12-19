
library(readr)
library(dplyr)
library(stringr)

# Load TSV
df <- read_tsv("bed_sexchr_info.tsv", show_col_types = FALSE)

# Function to get chromosome names from a BED file
get_bed_chroms <- function(filepath) {

  # prepend directory
  filepath <- file.path("sorted_beds", filepath)

  # if file doesn't exist
  if (!file.exists(filepath)) return("NOT_FOUND")

  # Read BED allowing blank lines
  bed <- read_tsv(
    filepath,
    col_names = FALSE,
    col_types = cols(.default = "c"),
    comment = "#",
    progress = FALSE
  )

  # If file is completely empty (no rows or no columns)
  if (ncol(bed) == 0 || nrow(bed) == 0) {
    return("EMPTY")
  }

  # Extract first column (chrom names)
  chroms <- bed[[1]]

  # Remove NA, empty, whitespace-only chromosome entries
  chroms <- chroms[!is.na(chroms) & str_trim(chroms) != ""]

  # If after cleaning there is nothing informative, treat as empty
  if (length(chroms) == 0) {
    return("EMPTY")   # empty or whitespace-only file
  }

  return(unique(chroms))
}

# Loop through rows
df_out <- df %>%
  rowwise() %>%
  mutate(
    bed_chroms = list(get_bed_chroms(Filename)),
    Match_type = case_when(
      identical(bed_chroms, "NOT_FOUND") ~ "bed_not_found",
      identical(bed_chroms, "EMPTY")     ~ "empty_bed",
      GenBank_accession %in% bed_chroms &
        RefSeq_accession %in% bed_chroms ~ "GenBank+RefSeq",
      GenBank_accession %in% bed_chroms  ~ "GenBank",
      RefSeq_accession %in% bed_chroms   ~ "RefSeq",
      TRUE ~ "neither"
    )
  ) %>%
  ungroup()

# Print output

df_neither <- df_out %>% filter(Match_type == "neither")
print(df_neither)

print(df_neither$Species)
# Optional: write result to file
write_tsv(df_out, "bed_sexchr_info_with_matches.tsv")
