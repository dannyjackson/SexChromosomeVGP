# Plot Subclades

# Chondricthyes
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Chondrichthyes
source myconda
mamba activate genespace

  Narcine_bancroftii \
  Pristis_pectinata \
  Hypanus_sabinus \
  Mobula_birostris \
  Heptranchias_perlo \
  Pristiophorus_japonicus \
  Heterodontus_francisci \
  Stegostoma_tigrinum \
  Hemiscyllium_ocellatum \
  Mustelus_asterias \
  Scyliorhinus_canicula \
  Cetorhinus_maximus \
  Carcharodon_carcharias


cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Chondrichthyes

15006456
R

library(GENESPACE)
library(ggplot2)
CLADE <- "Chondrichthyes"
wd <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  CLADE
)

f <- file.path(
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

SPECIES=c("Heptranchias_perlo", "Pristiophorus_japonicus", "Heterodontus_francisci", "Stegostoma_tigrinum", "Hemiscyllium_ocellatum", "Mustelus_asterias", "Scyliorhinus_canicula", "Cetorhinus_maximus", "Carcharodon_carcharias")

roi <- data.frame( 
    genome = c("Carcharodon_carcharias"),
    chr = c("X"), 
    color = c("#0072B2")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x9.pdf"), 
    plot = p, 
    width = 10, 
    height = 9, 
    units = "in" )

# Option 2

SPECIES=c("Narcine_bancroftii", "Pristis_pectinata", "Hypanus_sabinus", "Mobula_birostris", "Heptranchias_perlo", "Pristiophorus_japonicus", "Heterodontus_francisci", "Stegostoma_tigrinum", "Hemiscyllium_ocellatum", "Mustelus_asterias", "Scyliorhinus_canicula", "Cetorhinus_maximus", "Carcharodon_carcharias")

roi <- data.frame( 
    genome = c("Carcharodon_carcharias"),
    chr = c("X"), 
    color = c("#0072B2")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x13.Scyliorhinus_canicula.pdf"), 
    plot = p, 
    width = 10, 
    height = 13, 
    units = "in" )
```
# Shrews
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Shrews

R

library(GENESPACE)
library(ggplot2)
CLADE <- "Shrews"
wd <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  CLADE
)

f <- file.path(
  wd,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

SPECIES=c("Suncus_etruscus", "Sorex_araneus")

roi <- data.frame( 
    genome = "Sorex_araneus", 
    chr = c("X"), 
    color = c("#E69F00")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Shrews.10x6.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )

ggsave( 
    filename = "Shrews.10x4.pdf", 
    plot = p, 
    width = 10, 
    height = 4, 
    units = "in" )


ggsave( 
    filename = "Shrews.10x2.pdf", 
    plot = p, 
    width = 10, 
    height = 2, 
    units = "in" )
```
# Bat neo sex
```

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Bat_Neo


R

library(GENESPACE)
library(ggplot2)
CLADE <- "Bat_Neo"
wd <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  CLADE
)

f <- file.path(
  wd,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

SPECIES=c("Glossophaga_mutica", "Artibeus_intermedius")

roi <- data.frame( 
    genome = "Artibeus_intermedius", 
    chr = c("X"), 
    color = c("#E69F00")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x2.pdf"), 
    plot = p, 
    width = 10, 
    height = 2, 
    units = "in" )

```
# Snakes
```
CLADE="Snakes"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Vipera_latastei \
  Vipera_berus \
  Vipera_ursinii \
  Erythrolamprus_reginae \
  Thamnophis_elegans \
  Natrix_helvetica
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Squamates
```
CLADE="Squamates"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Rhineura_floridana \
  Lacerta_agilis \
  Zootoca_vivipara \
  Podarcis_siculus \
  Podarcis_muralis \
  Podarcis_liolepis \
  Podarcis_bocagei \
  Podarcis_vaucheri \
  Podarcis_pityusensis \
  Podarcis_tiliguerta \
  Podarcis_raffonei \
  Podarcis_filfolensis \
  Podarcis_melisellensis \
  Podarcis_gaigeae \
  Podarcis_cretensis \
  Podarcis_erhardii
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Strisores
```
for SPECIES in \
  Caprimulgus_europaeus \
  Nyctibius_grandis \
  Podargus_strigoides \
  Aegotheles_albertisi \
  Hemiprocne_comata \
  Apus_apus \
  Phaethornis_superciliosus \
  Heliangelus_exortis \
  Calypte_anna


cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Strisores


R

library(GENESPACE)
library(ggplot2)
CLADE <- "Strisores"
wd <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  CLADE
)

f <- file.path(
  wd,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

SPECIES=c("Caprimulgus_europaeus", "Nyctibius_grandis", "Podargus_strigoides")

roi <- data.frame( 
    genome = c("Nyctibius_grandis","Podargus_strigoides"), 
    chr = c("Z", "Z"), 
    color = c("#00796B")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    minChrLen2plot = 1,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x3.pdf"), 
    plot = p, 
    width = 10, 
    height = 2, 
    units = "in" )
```
# Pigeons
```
for SPECIES in \
  Pterocles_gutturalis \
  Caloenas_nicobarica \
  Patagioenas_fasciata \
  Columba_livia \
  Nesoenas_mayeri \
  Streptopelia_turtur \
  Streptopelia_decaocto


cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Pigeons


R

library(GENESPACE)
library(ggplot2)
CLADE <- "Pigeons"
wd <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  CLADE
)

f <- file.path(
  wd,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

SPECIES=c("Pterocles_gutturalis", "Caloenas_nicobarica", "Patagioenas_fasciata")


roi <- data.frame( 
    genome = c("Pterocles_gutturalis", "Caloenas_nicobarica"), 
    chr = c("Z", "Z"), 
    color = c("#00796B", "#00796B")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    minChrLen2plot = 1,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x3.pdf"), 
    plot = p, 
    width = 10, 
    height = 3, 
    units = "in" )
```
# Spoonbills
```
CLADE="Spoonbills"

  Pelecanus_crispus \
  Theristicus_caerulescens \
  Platalea_leucorodia \
  Morus_bassanus

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Spoonbills


R

library(GENESPACE)
library(ggplot2)
CLADE <- "Spoonbills"
wd <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  CLADE
)

f <- file.path(
  wd,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

SPECIES=c("Pelecanus_crispus", "Theristicus_caerulescens", "Platalea_leucorodia", "Morus_bassanus")


roi <- data.frame( 
    genome = c("Pelecanus_crispus", "Theristicus_caerulescens", "Platalea_leucorodia"), 
    chr = c("Z", "Z", "Z"), 
    color = c("#00796B", "#00796B", "#00796B")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    minChrLen2plot = 1,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x4.pdf"), 
    plot = p, 
    width = 10, 
    height = 2, 
    units = "in" )
```
# manually run chondrichthyes and snakes
```
CLADE="Chondrichthyes"
./1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/

CLADE="Snakes"
./1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Squamates_broad
```
CLADE="Squamates"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Dibamus_smithi \
  Anolis_sagrei \
  Cyclura_pinguis \
  Vipera_latastei
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Hawks
```

for SPECIES in \
  Athene_noctua \
  Strix_aluco \
  Sarcoramphus_papa \
  Gypaetus_barbatus \
  Aquila_chrysaetos \
  Morphnus_guianensis \
  Harpia_harpyja \
  Accipiter_gentilis \
  Buteo_buteo \
  Haliaeetus_albicilla



cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Hawks


R

library(GENESPACE)
library(ggplot2)
CLADE <- "Hawks"
wd <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  CLADE
)

f <- file.path(
  wd,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

SPECIES=c("Athene_noctua", "Strix_aluco", "Sarcoramphus_papa", "Gypaetus_barbatus", "Aquila_chrysaetos", "Morphnus_guianensis", "Harpia_harpyja", "Accipiter_gentilis", "Buteo_buteo", "Haliaeetus_albicilla")


roi <- data.frame( 
    genome = c("Sarcoramphus_papa", "Gypaetus_barbatus", "Aquila_chrysaetos", "Morphnus_guianensis", "Harpia_harpyja", "Accipiter_gentilis", "Buteo_buteo", "Haliaeetus_albicilla"), 
    chr = c("Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z"), 
    color = c("#00796B", "#00796B", "#00796B", "#00796B", "#00796B", "#00796B", "#00796B", "#00796B")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    minChrLen2plot = 10,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x8.pdf"), 
    plot = p, 
    width = 10, 
    height = 8, 
    units = "in" )
```
# Coraciimorphs
```
for SPECIES in \
  Leptosomus_discolor \
  Pogoniulus_pusillus \
  Dryobates_pubescens \
  Colius_striatus

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Coraciimorphs


R

library(GENESPACE)
library(ggplot2)
CLADE <- "Coraciimorphs"
wd <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  CLADE
)

f <- file.path(
  wd,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

SPECIES=c("Colius_striatus", "Dryobates_pubescens", "Pogoniulus_pusillus", "Leptosomus_discolor")


roi <- data.frame( 
    genome = c("Leptosomus_discolor", "Dryobates_pubescens", "Colius_striatus"), 
    chr = c("Z", "Z", "Z"), 
    color = c("#00796B", "#00796B", "#00796B")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    minChrLen2plot = 10,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x4.pdf"), 
    plot = p, 
    width = 10, 
    height = 4, 
    units = "in" )
```
# Parrots
```
for SPECIES in \
  Strigops_habroptilus \
  Amazona_ochrocephala \
  Guaruba_guaruba \
  Ara_ararauna \
  Psittacula_echo \
  Lathamus_discolor \
  Melopsittacus_undulatus


cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Parrots


R

library(GENESPACE)
library(ggplot2)
CLADE <- "Parrots"
wd <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  CLADE
)

f <- file.path(
  wd,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

SPECIES=c("Strigops_habroptilus", "Amazona_ochrocephala", "Guaruba_guaruba", "Ara_ararauna", "Psittacula_echo", "Lathamus_discolor", "Melopsittacus_undulatus", "Acanthisitta_chloris")


roi <- data.frame( 
    genome = c("Strigops_habroptilus", "Amazona_ochrocephala", "Guaruba_guaruba", "Ara_ararauna", "Psittacula_echo", "Lathamus_discolor", "Melopsittacus_undulatus"), 
    chr = c("Z", "Z", "Z", "Z", "Z", "Z", "Z"), 
    color = c("#00796B", "#00796B", "#00796B", "#00796B", "#00796B", "#00796B", "#00796B")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    minChrLen2plot = 10,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x8.pdf"), 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )
```
# Sylviida
```
CLADE="Sylviida"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Poecile_atricapillus \
  Hirundo_rustica \
  Zosterops_lateralis \
  Sylvia_atricapilla \
  Sylvia_borin



cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Sylviida


R

library(GENESPACE)
library(ggplot2)
CLADE <- "Sylviida"
wd <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  CLADE
)

f <- file.path(
  wd,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

SPECIES=c("Poecile_atricapillus", "Hirundo_rustica")


roi <- data.frame( 
    genome = c("Poecile_atricapillus", "Hirundo_rustica"), 
    chr = c("Z", "Z"), 
    color = c("#00796B", "#00796B")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    minChrLen2plot = 10,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x2.pdf"), 
    plot = p, 
    width = 10, 
    height = 2, 
    units = "in" )

```
# Muscicapoids
```
CLADE="Muscicapoids"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Acridotheres_tristis \
  Cinclus_cinclus \
  Catharus_ustulatus \
  Erithacus_rubecula \
  Oenanthe_melanoleuca
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Monotremes
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Monotremes


R

library(GENESPACE)
library(ggplot2)
CLADE <- "Monotremes"
wd <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  CLADE
)

f <- file.path(
  wd,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

SPECIES=c("Tachyglossus_aculeatus", "Ornithorhynchus_anatinus")

roi <- data.frame( 
    genome = c("Tachyglossus_aculeatus", "Tachyglossus_aculeatus", "Tachyglossus_aculeatus", "Tachyglossus_aculeatus", "Tachyglossus_aculeatus", "Ornithorhynchus_anatinus", "Ornithorhynchus_anatinus", "Ornithorhynchus_anatinus", "Ornithorhynchus_anatinus", "Ornithorhynchus_anatinus"),
    chr = c("X1", "X2", "X3", "X4", "X5", "X1", "X2", "X3", "X4", "X5"), 
    color = c("#E69F00", "#E69F00", "#E69F00", "#E69F00", "#e53600", "#E69F00", "#E69F00", "#E69F00", "#e53600", "#E69F00")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x2.pdf"), 
    plot = p, 
    width = 10, 
    height = 2, 
    units = "in" )
```
# Marsupials
```
CLADE="Marsupials"
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}

R

library(GENESPACE)
library(ggplot2)
CLADE <- "Marsupials"
wd <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  CLADE
)

f <- file.path(
  wd,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

SPECIES=c("Macrotis_lagotis",
    "Trichosurus_vulpecula")

roi <- data.frame( 
    genome = c("Sminthopsis_crassicaudata", "Macrotis_lagotis"),
    chr = c("X", "X"), 
    color = c("#E69F00", "#E69F00")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x2.pdf"), 
    plot = p, 
    width = 10, 
    height = 2, 
    units = "in" )
```
# Mammals_broad
```
CLADE="Mammals_broad"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Ornithorhynchus_anatinus \
  Monodelphis_domestica \
  Macrotis_lagotis \
  Homo_sapiens


CLADE="Mammals_broad"
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}

R

library(GENESPACE)
library(ggplot2)
CLADE <- "Mammals_broad"
wd <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  CLADE
)

f <- file.path(
  wd,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

SPECIES=c("Ornithorhynchus_anatinus", "Monodelphis_domestica", "Homo_sapiens")

roi <- data.frame( 
    genome = c("Ornithorhynchus_anatinus", "Ornithorhynchus_anatinus", "Ornithorhynchus_anatinus", "Ornithorhynchus_anatinus", "Ornithorhynchus_anatinus", "Monodelphis_domestica", "Homo_sapiens"),
    chr = c("X1", "X2", "X3", "X4", "X5", "X", "X"), 
    color = c("#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    minChrLen2plot = 1,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x6.pdf"), 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )


SPECIES=c("Homo_sapiens", "Monodelphis_domestica")

roi <- data.frame( 
    genome = c("Monodelphis_domestica", "Homo_sapiens"),
    chr = c("X", "X"), 
    color = c("#E69F00", "#E69F00")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    minChrLen2plot = 1,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x4.pdf"), 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )
```
# Eutherian_fusion
```
for SPECIES in \
  Monodelphis_domestica \
  Dromiciops_gliroides \
  Sminthopsis_crassicaudata \
  Sarcophilus_harrisii \
  Trichosurus_vulpecula \
  Macropus_eugenii \
  Homo_sapiens


CLADE="Eutherian_fusion"
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}

R

library(GENESPACE)
library(ggplot2)
CLADE <- "Eutherian_fusion"
wd <- file.path(
  "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace",
  CLADE
)

f <- file.path(
  wd,
  "results",
  "gsParams.rda"
)

load(f, verbose = TRUE)

SPECIES=c("Monodelphis_domestica", "Dromiciops_gliroides", "Sminthopsis_crassicaudata", "Sarcophilus_harrisii", "Trichosurus_vulpecula", "Macropus_eugenii", "Homo_sapiens")

roi <- data.frame( 
    genome = c("Homo_sapiens"),
    chr = c("X"), 
    color = c("#E69F00")
    )

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    minChrLen2plot = 1,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = paste0(CLADE, ".10x7.pdf"), 
    plot = p, 
    width = 10, 
    height = 7, 
    units = "in" )
```