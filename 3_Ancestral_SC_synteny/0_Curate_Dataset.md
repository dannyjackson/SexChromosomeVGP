# Goal:
Create a plot that has chicken chromosomes on the X axis and unique SC systems on the Y. For each unique SC on the Y, there will be a horizontal line from start-stop of each syntenic block for a given clade. This will be min and max of all species that have that SC.

# Dataset:
## Sharks and Rays
```
Carcharodon_carcharias
```
## Ray finned fish
```
Colia_mystus
Argentina_silus
Gasterosteus_aculeatus
```
## Mammals
```
Sarcophilus_harrisii
```
## Lepidosaurs
```
Cyclura_pinguis
Vipera_latastei
Podarcis_raffonei
```
## Turtles
```
```
## Frogs
```
Pseudacris_triseriata
```
## Birds
```
Gallus_gallus
```
module load samtools
samtools faidx /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Gallus_gallus/ncbi_dataset/data/GCF_016700215.2/GCF_016700215.2_bGalGal1.pat.whiteleghornlayer.GRCg7w_genomic.fna
../Gallus_gallus/combined_phaseblks/

```
1,196449156
2,149539284
3,110642502
4,90861225
5,59506338
6,36220557
7,36382834
8,29578256
9,23733309
10,20453248
11,19638187
12,20119077
13,17905061
14,15331188
15,12703657
16,2706039
17,11092391
18,11623896
19,10455293
20,14265659
21,6970754
22,4686657
23,6253421
24,6478339
25,3067737
26,5349051
27,5228753
28,5437364
29,726478
30,755666
31,2457334
32,125424
33,3839931
34,3469343
35,554126
36,358375
37,157853
38,667312
39,177356
Z,86044486
```
# Save as test_SC.list.txt
```
Carcharodon_carcharias
Colia_mystus
Argentina_silus
Gasterosteus_aculeatus
Sarcophilus_harrisii
Cyclura_pinguis
Vipera_latastei
Podarcis_raffonei
```
# Subset combined_phased_blks
```
awk -F',' '
BEGIN {
    OFS=","
}

# species list
FNR==NR {
    gsub(/\r/, "", $1)
    if ($1 != "") {
        want[$1] = 1
        order[++n] = $1
    }
    next
}

# header row of CSV
FNR==1 {
    for (i=1; i<=NF; i++) {
        col[$i] = i
    }
    next
}

{
    sp = $(col["partner"])

    if (!(sp in want)) next
    if ($(col["genome1"]) == "Gallus_gallus") next
    if ($(col["genome2"]) != "Gallus_gallus") next
    if ($(col["chr1"]) != "X" && $(col["chr1"]) != "Z") next

    region = $(col["chr2"]) ":" $(col["startBp2"]) "-" $(col["endBp2"])

    key = sp SUBSEP region
    if (!(key in seen)) {
        seen[key] = 1
        if (sp in regions) {
            regions[sp] = regions[sp] ";" region
        } else {
            regions[sp] = region
        }
    }
}

END {
    print "Species","Gg_regions"
    for (i=1; i<=n; i++) {
        sp = order[i]
        print sp, regions[sp]
    }
}
' test_SC.list.txt ../Gallus_gallus/combined_phaseblks/Gallus_gallus_sexchrs.phasedBlks.csv > gg_regions.csv
```
# Plot the data
```
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript plot_gg_regions.R gg_regions.csv gg_chr_len.csv output.png")
}

gg_regions_file <- "gg_regions.csv"
gg_chr_len_file <- "gg_chr_len.csv"
output_file <- "gg_regions_plot.png"


gg_regions_file <- args[1]
gg_chr_len_file <- args[2]
output_file <- args[3]

# Read chromosome lengths
chr_len <- read_csv(
  gg_chr_len_file,
  col_names = c("chr", "length"),
  show_col_types = FALSE
) %>%
  mutate(chr = as.character(chr),
         length = as.numeric(length))

# Preserve chromosome order from file
chr_order <- chr_len$chr

# Build cumulative offsets so chromosomes are laid out end-to-end
gap <- 5000000

chr_len <- chr_len %>%
  mutate(offset = cumsum(lag(length + gap, default = 0)),
         center = offset + length / 2)

# Read regions
regions <- read_csv(gg_regions_file, show_col_types = FALSE) %>%
  mutate(Gg_regions = ifelse(is.na(Gg_regions), "", Gg_regions))

# Expand semicolon-separated blocks
regions_long <- regions %>%
  separate_rows(Gg_regions, sep = ";") %>%
  mutate(Gg_regions = str_trim(Gg_regions)) %>%
  filter(Gg_regions != "")

# Parse chr:start-end
regions_long <- regions_long %>%
  extract(
    Gg_regions,
    into = c("chr", "start", "end"),
    regex = "^([^:]+):(\\d+)-(\\d+)$",
    remove = FALSE
  ) %>%
  mutate(
    chr = as.character(chr),
    start = as.numeric(start),
    end = as.numeric(end)
  ) %>%
  filter(!is.na(chr), !is.na(start), !is.na(end)) %>%
  mutate(
    xmin_local = pmin(start, end),
    xmax_local = pmax(start, end)
  ) %>%
  left_join(chr_len, by = "chr") %>%
  filter(!is.na(offset)) %>%
  mutate(
    xmin = offset + xmin_local,
    xmax = offset + xmax_local
  )

# Preserve species order from input file
species_order <- unique(regions$Species)

regions_long <- regions_long %>%
  mutate(Species = factor(Species, levels = rev(species_order)))

# Background rectangles for chromosomes
chr_bg <- chr_len %>%
  mutate(
    xmin = offset,
    xmax = offset + length
  )

p <- ggplot() +
  geom_rect(
    data = chr_bg,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.08
  ) +
  geom_segment(
    data = regions_long,
    aes(x = xmin, xend = xmax, y = Species, yend = Species),
    linewidth = 1.2,
    lineend = "butt"
  ) +
  scale_x_continuous(
    breaks = chr_len$center,
    labels = chr_len$chr,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Gallus gallus chromosomes",
    y = "Unique SC systems",
    title = "Chicken syntenic regions by SC system"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

height_inches <- max(4, 0.3 * length(species_order))

ggsave(
  output_file,
  plot = p,
  width = 16,
  height = 3,
  dpi = 300
)
```