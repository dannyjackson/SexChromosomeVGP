library(readr)
library(dplyr)
library(purrr)
library(stringr)

roots <- c("genespace")

files <- unlist(lapply(roots, function(r) {
  if (dir.exists(r)) list.files(r, pattern = "_phasedBlks\\.csv$", recursive = TRUE, full.names = TRUE) else character()
}), use.names = FALSE) %>% sort()

stopifnot(length(files) > 0)

parse_meta <- function(path) {
  # Example: Z/Accipiter_gentilis/riparian/Accipiter_gentilis_phasedBlks.csv
  parts <- strsplit(path, .Platform$file.sep, fixed = TRUE)[[1]]
  tibble(
    topDir = parts[1] %||% NA_character_,
    species = parts[2] %||% NA_character_,
    subdir = parts[3] %||% NA_character_,
    sourceFile = path
  )
}

# Read all, keep union of columns (fills missing with NA), add provenance cols

df <- map_dfr(files, function(f) {
  dat <- read_csv(
    f,
    col_types = cols(.default = col_character()),
    progress = FALSE,
    show_col_types = FALSE
  )
  bind_cols(parse_meta(f), dat)
})

write_csv(df, "all_phasedBlks_combined.csv")
cat("Wrote all_phasedBlks_combined.csv with", nrow(df), "rows from", length(files), "files\n")

# rename instances of "CP162266.1" to X
df$chr1[df$chr1 == "CP162266.1"] <- "X"
df$chr2[df$chr2 == "CP162266.1"] <- "X"
  
# Filter to just rows where a sex chromosome is in either chr 1 or chr 2
df2 <- df %>%
  filter(grepl("[X-Zx-z]", chr1) | grepl("[X-Zx-z]", chr2))

# remove rows that are chicken Z > other species 
df3 <- df2 %>%
  filter(!(genome1 == "Gallus_gallus_REF" & chr1 == "Z" & !grepl("^(Z|X)", chr2))) %>%
  filter(!(genome2 == "Gallus_gallus_REF" & chr2 == "Z" & !grepl("^(Z|X)", chr1)))

VGP_DATAFRAME <- read_csv("/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/VGP_OrdinalList_Phase1Freeze_v1.1_sex_chroms_seqCov_HalfDeep_SCINKD_May8.25_DwnldDec16.25.csv")
# rename Canis lupus baileyi to Canis lupus
VGP_DATAFRAME$Scientific.Name[VGP_DATAFRAME$Scientific.Name == "Canis lupus baileyi"] <- "Canis lupus"

# rename to Guaruba_guaruba to Guaruba guarouba in df
df3$species[df3$species == "Guaruba_guaruba"] <- "Guaruba_guarouba"

# Join df3 and VGP_DATAFRAME on df3$species and VGP_DATAFRAME$Scientific.Name for columns Lineage, Superorder


df4 <- df3 %>%
  mutate(species_join = gsub("_", " ", species)) %>%
  left_join(
    VGP_DATAFRAME %>% 
      transmute(Scientific.Name, Lineage, Superorder,
                Scientific_join = trimws(Scientific.Name)),
    by = c("species_join" = "Scientific_join")
  ) %>%
  select(-species_join)

unmatched <- df4 %>% filter(is.na(Lineage) & is.na(Superorder)) %>% distinct(species)
unmatched


write_csv(df4, "sexchrs.phasedBlks.csv")
unique_chrs <- unique(df4[c("Lineage", "Superorder", "genome1","genome2","chr1","chr2")])
write_csv(unique_chrs, "unique_chrs.csv")

unique_chrs <- read.csv("unique_chrs.csv")

# condense to just one instance of matching, retaining info if a chr is only matched once not twice
gg_name <- "Gallus_gallus_REF"   # change if your column uses a different exact string

out <- unique_chrs %>%
  # keep rows where chicken is in exactly one of the genome columns
  filter(xor(genome1 == gg_name, genome2 == gg_name)) %>%
  mutate(
    Species = if_else(genome1 == gg_name, genome2, genome1),
    SpChr   = if_else(genome1 == gg_name, chr2,   chr1),

    # chicken chromosome when chicken is genome1 -> chicken chr is chr1
    gg_chr1 = if_else(genome1 == gg_name, chr1, NA_character_),

    # chicken chromosome when chicken is genome2 -> chicken chr is chr2
    gg_chr2 = if_else(genome2 == gg_name, chr2, NA_character_)
  ) %>%
  group_by(Species, SpChr) %>%
  summarise(
    gg_chr1 = list(sort(unique(na.omit(gg_chr1)))),
    gg_chr2 = list(sort(unique(na.omit(gg_chr2)))),

    GgChr_both      = list(intersect(gg_chr1[[1]], gg_chr2[[1]])),
    GgChr_chr1_only = list(setdiff(gg_chr1[[1]], gg_chr2[[1]])),
    GgChr_chr2_only = list(setdiff(gg_chr2[[1]], gg_chr1[[1]])),
    .groups = "drop"
  ) %>%
  select(Species, SpChr, GgChr_both, GgChr_chr1_only, GgChr_chr2_only)


out_csv <- out %>%
  mutate(across(starts_with("GgChr"),
                ~ vapply(.x, \(v) paste(v, collapse = ","), character(1))))

readr::write_csv(out_csv, "syntentic_chromosomes.sex_shared.Gallus_gallus_REF.csv")