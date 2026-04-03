#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript combine_phaseblks.r <FOCAL_SPECIES>\nExample: Rscript combine_phaseblks.r Anolis_sagrei")
}
FOCAL <- args[1]

BASE_DIR <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace"
ROOT <- file.path(BASE_DIR, FOCAL, "sexshared")
ROOT_OUT <- file.path(BASE_DIR, FOCAL)

if (!dir.exists(ROOT)) {
  stop("Root directory does not exist: ", ROOT)
}

OUT_DIR <- file.path(ROOT_OUT, "combined_phaseblks")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

files <- list.files(
  ROOT,
  pattern = "_phasedBlks\\.csv$",
  recursive = TRUE,
  full.names = TRUE
) %>% sort()

if (length(files) == 0) {
  stop("No *_phasedBlks.csv files found under: ", ROOT)
}

# Parse metadata from full path:
# .../<FOCAL>/sexshared/<PARTNER>/riparian/<FILE_SPECIES>_phasedBlks.csv
parse_meta <- function(path, focal = FOCAL) {
  partner <- str_match(path, paste0("/", focal, "/sexshared/([^/]+)/"))[, 2]
  file_species <- str_match(path, "/riparian/([^/]+)_phasedBlks\\.csv$")[, 2]

  tibble(
    focal = focal,
    partner = partner %||% NA_character_,
    file_species = file_species %||% NA_character_,
    sourceFile = path
  )
}

`%||%` <- function(x, y) {
  if (!is.null(x) && length(x) > 0 && !is.na(x)) x else y
}

# Read all phasedBlks, keep union of columns, add provenance cols
df <- map_dfr(files, function(f) {
  dat <- read_csv(
    f,
    col_types = cols(.default = col_character()),
    progress = FALSE,
    show_col_types = FALSE
  )
  bind_cols(parse_meta(f), dat)
})

# Normalize chromosome label
if ("chr1" %in% names(df)) df$chr1[df$chr1 == "CP162266.1"] <- "X"
if ("chr2" %in% names(df)) df$chr2[df$chr2 == "CP162266.1"] <- "X"

# Filter to rows where a sex chromosome is in either chr1 or chr2
df2 <- df %>%
  filter(grepl("[X-Zx-z]", chr1) | grepl("[X-Zx-z]", chr2))

# Remove rows that are reference Z/X > other species (same logic as your script)
df3 <- df2 %>%
  filter(!(genome1 == FOCAL & grepl("^(Z|X)", chr1) & !grepl("^(Z|X)", chr2))) %>%
  filter(!(genome2 == FOCAL & grepl("^(Z|X)", chr2) & !grepl("^(Z|X)", chr1)))

# Load VGP dataframe
VGP_DATAFRAME <- read_csv(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/VGP_OrdinalList_Phase1Freeze_v1.1_sex_chroms_seqCov_HalfDeep_SCINKD_May8.25_DwnldDec16.25.csv",
  show_col_types = FALSE
)

# Rename Canis lupus baileyi -> Canis lupus
VGP_DATAFRAME$Scientific.Name[VGP_DATAFRAME$Scientific.Name == "Canis lupus baileyi"] <- "Canis lupus"

# Rename Guaruba_guaruba -> Guaruba_guarouba in df (your original fix)
df3$file_species[df3$file_species == "Guaruba_guaruba"] <- "Guaruba_guarouba"

# Join on species (file_species with "_" -> " ")
df4 <- df3 %>%
  mutate(species_join = gsub("_", " ", file_species)) %>%
  left_join(
    VGP_DATAFRAME %>%
      transmute(
        Scientific_join = trimws(Scientific.Name),
        Lineage,
        Superorder
      ),
    by = c("species_join" = "Scientific_join")
  ) %>%
  select(-species_join)

unmatched <- df4 %>%
  filter(is.na(Lineage) & is.na(Superorder)) %>%
  distinct(file_species)

# Write combined outputs
all_out <- file.path(OUT_DIR, paste0(FOCAL, "_all_phasedBlks_combined.csv"))
sex_out <- file.path(OUT_DIR, paste0(FOCAL, "_sexchrs.phasedBlks.csv"))
uniq_out <- file.path(OUT_DIR, paste0(FOCAL, "_unique_chrs.csv"))
syn_out  <- file.path(OUT_DIR, paste0(FOCAL, "_syntentic_chromosomes.sex_shared.csv"))
unm_out  <- file.path(OUT_DIR, paste0(FOCAL, "_unmatched_species.csv"))

write_csv(df,  all_out)
write_csv(df4, sex_out)
write_csv(unmatched, unm_out)

# Unique chromosome pairs table
unique_chrs <- unique(df4[c("Lineage", "Superorder", "genome1","genome2","chr1","chr2")])
write_csv(unique_chrs, uniq_out)

# Condense to just one instance of matching against the reference
ref_name <- FOCAL

out <- unique_chrs %>%
  filter(xor(genome1 == ref_name, genome2 == ref_name)) %>%
  mutate(
    Species = if_else(genome1 == ref_name, genome2, genome1),
    SpChr   = if_else(genome1 == ref_name, chr2,   chr1),
    Ref_Chr1 = if_else(genome1 == ref_name, chr1, NA_character_),
    Ref_Chr2 = if_else(genome2 == ref_name, chr2, NA_character_)
  ) %>%
  group_by(Species, SpChr) %>%
  summarise(
    Ref_Chr1 = list(sort(unique(na.omit(Ref_Chr1)))),
    Ref_Chr2 = list(sort(unique(na.omit(Ref_Chr2)))),
    RefChr_both      = list(intersect(Ref_Chr1[[1]], Ref_Chr2[[1]])),
    RefChr_chr1_only = list(setdiff(Ref_Chr1[[1]], Ref_Chr2[[1]])),
    RefChr_chr2_only = list(setdiff(Ref_Chr2[[1]], Ref_Chr1[[1]])),
    .groups = "drop"
  ) %>%
  select(Species, SpChr, RefChr_both, RefChr_chr1_only, RefChr_chr2_only)

out_csv <- out %>%
  mutate(across(starts_with("RefChr"),
                ~ vapply(.x, \(v) paste(v, collapse = ","), character(1))))

write_csv(out_csv, syn_out)

cat("FOCAL:", FOCAL, "\n")
cat("Searched:", ROOT, "\n")
cat("Found", length(files), "phasedBlks files\n")
cat("Wrote:\n  ", all_out, "\n  ", sex_out, "\n  ", uniq_out, "\n  ", syn_out, "\n  ", unm_out, "\n", sep = "")
cat("Rows combined:", nrow(df), "\n")
cat("Rows after sex-chr filters + join:", nrow(df4), "\n")
if (nrow(unmatched) > 0) {
  cat("Unmatched species written to:", unm_out, "\n")
}