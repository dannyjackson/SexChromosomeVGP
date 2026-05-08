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
```
## Nyctibius_grandis
```
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

SPECIES=c("Caprimulgus_europaeus", "Nyctibius_grandis")

roi <- data.frame( 
    genome = c("Caprimulgus_europaeus","Nyctibius_grandis"), 
    chr = c("16", "Z"), 
    color = c("#00796B","#d3d3d3")
    )
invert <- data.frame(genome = "Nyctibius_grandis", chr = "Z")

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
    backgroundColor = NA, 
    invertTheseChrs = invert,
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Nyctibus.2x4.pdf", 
    plot = p, 
    width = 4, 
    height = 2, 
    units = "in" )
```
## Podargus_strigoides
```
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

SPECIES=c("Caprimulgus_europaeus", "Podargus_strigoides")

roi <- data.frame( 
    genome = c("Caprimulgus_europaeus","Caprimulgus_europaeus","Podargus_strigoides"), 
    chr = c("5", "6", "Z"), 
    color = c("#00796B","#00796B","#d3d3d3")
    )
    
invert <- data.frame(
    genome = c("Podargus_strigoides","Caprimulgus_europaeus"),
    chr = c("Z", "5"))

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
    backgroundColor = NA, 
    invertTheseChrs = invert,
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Podargus.2x4.pdf", 
    plot = p, 
    width = 4, 
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

```
## Pterocles_gutturalis
```
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

SPECIES=c("Patagioenas_fasciata", "Pterocles_gutturalis")


roi <- data.frame( 
    genome = c("Patagioenas_fasciata","Pterocles_gutturalis"), 
    chr = c("8","Z"), 
    color = c("#00796B","#d3d3d3")
    )
    
invert <- data.frame(
    genome = c("Pterocles_gutturalis","Patagioenas_fasciata"),
    chr = c("Z","8"))

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
    backgroundColor = NA, 
    invertTheseChrs = invert,
    reorderBySynteny = FALSE)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Pterocles.2x4.pdf", 
    plot = p, 
    width = 4, 
    height = 2, 
    units = "in" )
```
## Caloenas_nicobarica
```
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

SPECIES=c("Patagioenas_fasciata", "Caloenas_nicobarica")


roi <- data.frame( 
    genome = c("Patagioenas_fasciata","Caloenas_nicobarica"), 
    chr = c("6","Z"), 
    color = c("#00796B","#000000")
    )
    
invert <- data.frame(
    genome = c("Caloenas_nicobarica"),
    chr = c("Z"))

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
    backgroundColor = NA, 
    invertTheseChrs = invert,
    reorderBySynteny = FALSE)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Caloenas.2x4.pdf", 
    plot = p, 
    width = 4, 
    height = 2, 
    units = "in" )
```
# Platalea_leucorodia
```
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

SPECIES=c("Morus_bassanus", "Platalea_leucorodia")
SPECIES=c("Platalea_leucorodia", "Pelecanus_crispus")


roi <- data.frame( 
    genome = c("Pelecanus_crispus","Platalea_leucorodia"), 
    chr = c("Z","Z"), 
    color = c("#00796B","#00796B")
    )

invert <- data.frame(
    genome = c("Pelecanus_crispus", "Platalea_leucorodia"),
    chr = c("Z", "Z"))

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
    backgroundColor = NA, 
    invertTheseChrs = invert,
    reorderBySynteny = FALSE)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Pelecanus.2x4.pdf", 
    plot = p, 
    width = 4, 
    height = 2, 
    units = "in" )
```
# Hawks and eagles

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

## Accipiter_gentilis
```
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

SPECIES=c("Sarcoramphus_papa", "Accipiter_gentilis")

roi <- data.frame( 
    genome = c("Sarcoramphus_papa","Accipiter_gentilis"), 
    chr = c("30","Z"), 
    color = c("#00796B","#000000")
    )
    
invert <- data.frame(
    genome = c("Sarcoramphus_papa"),
    chr = c("30"))

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
    backgroundColor = NA, 
    invertTheseChrs = invert,
    reorderBySynteny = FALSE)



p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Accipiter.2x4.pdf", 
    plot = p, 
    width = 4, 
    height = 2, 
    units = "in" )
```
## Harpy eagle
```
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

SPECIES=c("Accipiter_gentilis", "Harpia_harpyja")

roi <- data.frame( 
    genome = c("Accipiter_gentilis","Harpia_harpyja"), 
    chr = c("23","Z"), 
    color = c("#00796B","#000000")
    )
    
invert <- data.frame(
    genome = c("Caloenas_nicobarica"),
    chr = c("Z"))

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
    backgroundColor = NA, 
    reorderBySynteny = FALSE)

    invertTheseChrs = invert,

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Harpia.2x4.pdf", 
    plot = p, 
    width = 4, 
    height = 2, 
    units = "in" )
```
## All hawks
```
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

SPECIES=c("Accipiter_gentilis", "Harpia_harpyja")

roi <- data.frame( 
    genome = c("Accipiter_gentilis","Harpia_harpyja"), 
    chr = c("23","Z"), 
    color = c("#00796B","#000000")
    )
    
invert <- data.frame(
    genome = c("Caloenas_nicobarica"),
    chr = c("Z"))

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
    backgroundColor = NA, 
    reorderBySynteny = FALSE)

    invertTheseChrs = invert,

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Hawks.sexchrs.6x4.pdf", 
    plot = p, 
    width = 4, 
    height = 2, 
    units = "in" )
```
## All hawks
```
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

SPECIES=c("Athene_noctua", "Strix_aluco", "Sarcoramphus_papa", "Gypaetus_barbatus", "Morphnus_guianensis", "Harpia_harpyja", "Accipiter_gentilis", "Buteo_buteo", "Haliaeetus_albicilla")

roi <- data.frame( 
    genome = c("Accipiter_gentilis","Sarcoramphus_papa","Athene_noctua", "Strix_aluco","Gypaetus_barbatus", "Morphnus_guianensis", "Harpia_harpyja", "Accipiter_gentilis", "Buteo_buteo", "Haliaeetus_albicilla","Gypaetus_barbatus","Morphnus_guianensis"), 
    chr = c("23","30","Z","Z","Z","Z","Z","Z","Z","Z","6","2"), 
    color = c("#00796B","#00796B","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#d3d3d3","#d3d3d3")
    )

invert <- data.frame(
    genome = c("Buteo_buteo","Buteo_buteo","Haliaeetus_albicilla","Haliaeetus_albicilla","Morphnus_guianensis","Gypaetus_barbatus","Sarcoramphus_papa","Athene_noctua"),
    chr = c("21","Z","24","Z","2","6","30","Z"))

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
    backgroundColor = NA, 
    invertTheseChrs = invert,
    reorderBySynteny = FALSE)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Hawks.sexchrs.6x5.pdf", 
    plot = p, 
    width = 5, 
    height = 6, 
    units = "in" )
```
# Coraciimorphs
## All
```

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
    genome = c("Colius_striatus", "Colius_striatus", "Colius_striatus", "Dryobates_pubescens", "Dryobates_pubescens", "Dryobates_pubescens", "Pogoniulus_pusillus", "Pogoniulus_pusillus", "Pogoniulus_pusillus", "Pogoniulus_pusillus", "Pogoniulus_pusillus", "Leptosomus_discolor", "Leptosomus_discolor", "Leptosomus_discolor"), 
    chr = c("1", "24", "Z", "20", "36", "Z", "4", "27", "37", "40", "Z", "1", "25", "Z"), 
    color = c("#00796B", "#00796B","#000000","#00796B", "#00796B","#000000","#00796B", "#00796B","#00796B", "#00796B","#000000","#00796B", "#00796B","#000000")
    )
    
invert <- data.frame(
    genome = c("Caloenas_nicobarica"),
    chr = c("Z"))

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
    backgroundColor = NA, 
    reorderBySynteny = FALSE)

    invertTheseChrs = invert,

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "All.4x4.pdf", 
    plot = p, 
    width = 4, 
    height = 4, 
    units = "in" )
```

## Leptostomus discolor
```
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

SPECIES=c( "Trogon_surrucura", "Leptosomus_discolor")

roi <- data.frame( 
    genome = c("Trogon_surrucura", "Leptosomus_discolor"), 
    chr = c("24", "Z"), 
    color = c("#00796B","#000000")
    )
    
invert <- data.frame(
    genome = c("Trogon_surrucura"),
    chr = c("24"))

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
    backgroundColor = NA, 
    invertTheseChrs = invert,
    reorderBySynteny = FALSE)

    

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Leptostomus.2x4.pdf", 
    plot = p, 
    width = 4, 
    height = 2, 
    units = "in" )
```
## Dryobates pubescens
```
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

SPECIES=c( "Trogon_surrucura", "Dryobates_pubescens")

roi <- data.frame( 
    genome = c("Trogon_surrucura", "Dryobates_pubescens"), 
    chr = c("1", "Z"), 
    color = c("#00796B","#000000")
    )
    
invert <- data.frame(
    genome = c("Dryobates_pubescens", "Dryobates_pubescens", "Dryobates_pubescens", "Trogon_surrucura"),
    chr = c("7", "15", "Z", "Z"))

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
    backgroundColor = NA, 
    invertTheseChrs = invert,
    reorderBySynteny = FALSE)
    

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Dryobates.2x4.pdf", 
    plot = p, 
    width = 4, 
    height = 2, 
    units = "in" )
```
## Colius striatus
```
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

SPECIES=c( "Trogon_surrucura", "Colius_striatus")

roi <- data.frame( 
    genome = c("Trogon_surrucura", "Colius_striatus"), 
    chr = c("23", "Z"), 
    color = c("#00796B","#000000")
    )
    
invert <- data.frame(
    genome = c("Dryobates_pubescens"),
    chr = c("Z"))

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
    backgroundColor = NA, 
    reorderBySynteny = FALSE)
    

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Colius.2x4.pdf", 
    plot = p, 
    width = 4, 
    height = 2, 
    units = "in" )
```
# Poecile
```
R

library(GENESPACE)
library(ggplot2)
CLADE <- "Poecile"
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

SPECIES=c("Taeniopygia_guttata", "Poecile_atricapillus")

roi <- data.frame( 
    genome = c("Poecile_atricapillus", "Taeniopygia_guttata"), 
    chr = c("Z", "1A"), 
    color = c("#000000","#00796B")
    )
    
invert <- data.frame(
    genome = c("Taeniopygia_guttata"),
    chr = c("Z"))

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
    backgroundColor = NA, 
    invertTheseChrs = invert,
    reorderBySynteny = FALSE)
    
p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Peocile.2x4.pdf", 
    plot = p, 
    width = 4, 
    height = 2, 
    units = "in" )
```
# Acridotheres

```
R

library(GENESPACE)
library(ggplot2)
CLADE <- "Acridotheres"
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

SPECIES=c("Taeniopygia_guttata", "Acridotheres_tristis")

roi <- data.frame( 
    genome = c("Acridotheres_tristis", "Taeniopygia_guttata", "Taeniopygia_guttata"), 
    chr = c("Z", "5", "28"), 
    color = c("#000000","#00796B", "#00796B")
    )
    
invert <- data.frame(
    genome = c("Taeniopygia_guttata", "Acridotheres_tristis", "Acridotheres_tristis", "Taeniopygia_guttata"),
    chr = c("Z", "12", "19", "28"))

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
    backgroundColor = NA, 
    invertTheseChrs = invert,
    reorderBySynteny = FALSE)
    
p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Acridotheres.2x4.pdf", 
    plot = p, 
    width = 4, 
    height = 2, 
    units = "in" )
```
# Mammals -- shrews
```
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
    genome = c("Sorex_araneus", "Suncus_etruscus", "Suncus_etruscus", "Suncus_etruscus", "Suncus_etruscus", "Suncus_etruscus", "Suncus_etruscus", "Suncus_etruscus"), 
    chr = c("X", "1", "2", "5", "7", "9", "10", "12"), 
    color = c("#000000","#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00")
    )
    
invert <- data.frame(
    genome = c("Taeniopygia_guttata", "Acridotheres_tristis", "Acridotheres_tristis", "Taeniopygia_guttata"),
    chr = c("Z", "12", "19", "28"))

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
    backgroundColor = NA, 
    reorderBySynteny = FALSE)
    
    invertTheseChrs = invert,
p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Shrews.2x4.pdf", 
    plot = p, 
    width = 4, 
    height = 2, 
    units = "in" )
```
# Mammals -- bats
```
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
    genome = c("Artibeus_intermedius", "Glossophaga_mutica"), 
    chr = c("X", "12"), 
    color = c("#000000","#E69F00")
    )
    
invert <- data.frame(
    genome = c("Artibeus_intermedius"),
    chr = c("6"))

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
    backgroundColor = NA, 
    invertTheseChrs = invert,
    reorderBySynteny = FALSE)
    
p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Bats.2x4.pdf", 
    plot = p, 
    width = 4, 
    height = 2, 
    units = "in" )
```
# Sharks and Rays
```
R

library(GENESPACE)
library(ggplot2)
CLADE <- "Chondrichthyes"
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

SPECIES=c(
  "Mobula_birostris",
  "Hypanus_sabinus", 
  "Pristis_pectinata", 
  "Narcine_bancroftii", 
  "Amblyraja_radiata", 
  "Raja_brachyura"
  )

roi <- data.frame( 
    genome = c("Amblyraja_radiata", "Narcine_bancroftii", "Amblyraja_radiata", "Hypanus_sabinus"), 
    chr = c("X", "X", "16", "X2"), 
    color = c("#000000", "#0072B2", "#0072B2", "#0072B2")
    )
  
invert <- data.frame(
    genome = c("Mobula_birostris", "Amblyraja_radiata", "Amblyraja_radiata"),
    chr = c("X", "22", "33"))

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

chrLabFun <- function(x) {
  x2 <- gsub("^0", "",
        gsub("chr|scaf|chromosome|scaffold|^lg|_", "", tolower(x)))

  x2[x2 == "20"] <- "x2"
  return(x2)
}

ordFun <- function(x) {
  lab <- chrLabFun(x)

  is_x2 <- lab == "x2"

  chr_num <- suppressWarnings(as.numeric(gsub("[a-zA-Z]", "", lab)))

  data.table::frank(
    list(
      is_x2,      # FALSE first, TRUE last
      chr_num
    ),
    ties.method = "random"
  )
}

ripDat <- plot_riparian( 
    gsParam = gsParam, 
    highlightBed = roi, 
    genomeIDs = SPECIES, 
    minChrLen2plot = 1,
    refChrOrdFun = ordFun,
    chrLabFun = chrLabFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = NA, 
    invertTheseChrs = invert,
    reorderBySynteny = FALSE)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Sharks.4x4.pdf", 
    plot = p, 
    width = 4, 
    height = 4, 
    units = "in" )
```