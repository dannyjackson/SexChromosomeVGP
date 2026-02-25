# Investigations of putative neo sex chromosomes
Set up environment and directories
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/neo_sex/Morphnus_etc/genespace

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/neo_sex/Morphnus_etc/genespace

cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Morphnus_guianensis/bed .
cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Morphnus_guianensis/peptide .


cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Accipiter_gentilis/bed .
cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Accipiter_gentilis/peptide .


cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Strix_aluco/bed .
cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Strix_aluco/peptide .

```


Morphnus_guianensis
Accipiter_gentilis
Strix_alco
Gallus_gallus_REF

#!/bin/bash

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace

source myconda

mamba activate genespace

Rscript /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/neo_sex/Morphnus_etc/Morphnus_etc.genespace.r
```
sbatch \
  -c 16 \
  -t 1:00:00 \
  --mem-per-cpu=20G \
  --mail-user=jacksondan@nih.gov \
  --mail-type=ALL \
  --gres=lscratch:500 \
  --output=slurm_output/genespace_morphnus.%j \
  submit_morphnus.sh
```
library(GENESPACE)
library(ggplot2)
library(data.table)

setwd("/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/neo_sex/Morphnus_etc/")

gpar <- init_genespace(
  wd = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/neo_sex/Morphnus_etc/genespace",
  path2mcscanx = "/vf/users/jacksondan/conda/envs/genespace/bin/"
)

out <- run_genespace(gpar, overwrite = FALSE)

# Genome IDs to plot
genomeIDs <- c("Morphnus_guianensis", "Accipiter_gentilis", "Strix_aluco", "Gallus_gallus_REF")

ordFun <- function(x) {
  data.table::frank(
    list(grepl("[a-zA-Z]", x),
         suppressWarnings(as.numeric(gsub("[a-zA-Z]", "", x)))),
    ties.method = "random"
  )
}

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
)


ripDat <- plot_riparian(
  gsParam = out,
  refGenome = "Gallus_gallus_REF",
  genomeIDs = genomeIDs,
  refChrOrdFun = ordFun,
  addThemes = ggthemes,
  chrFill = "lightgrey",
  reorderBySynteny = FALSE
)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

```
# Redoing
```
library(GENESPACE)
library(ggplot2)
library(data.table)

load('/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/neo_sex/Morphnus_etc/genespace/results/gsParams.rda',
        verbose = TRUE)


# Genome IDs to plot
```
genomeIDs <- c("Morphnus_guianensis", "Accipiter_gentilis", "Strix_aluco")

ordFun <- function(x) {
  data.table::frank(
    list(grepl("[a-zA-Z]", x),
         suppressWarnings(as.numeric(gsub("[a-zA-Z]", "", x)))),
    ties.method = "random"
  )
}

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
)


ripDat <- plot_riparian(
  gsParam = gsParam,
  refGenome = "Gallus_gallus_REF",
  genomeIDs = genomeIDs,
  refChrOrdFun = ordFun,
  addThemes = ggthemes,
  chrFill = "lightgrey",
  reorderBySynteny = FALSE
)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

roi <- data.frame( 
    genome = c("Gallus_gallus_REF", "Strix_aluco", "Accipiter_gentilis", "Morphnus_guianensis"), 
    chr = c("Z", "Z", "Z", "CM098430.1"), 
    color = c("#c23d3d")) 

roi <- data.frame( 
    genome = c("Morphnus_guianensis"), 
    chr = c("CM098430.1"), 
    color = c("#c23d3d")) 

ripDat <- plot_riparian(
    gsParam = gsParam,
    highlightBed = roi, 
    backgroundColor = NULL,
    refGenome = "Gallus_gallus_REF",
    genomeIDs = genomeIDs,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    reorderBySynteny = FALSE
)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc.roi.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)


roi <- data.frame( 
    genome = c("Accipiter_gentilis"), 
    chr = c("Z"), 
    color = c("#c23d3d")) 


ripDat <- plot_riparian(
    gsParam = gsParam,
    highlightBed = roi, 
    backgroundColor = NULL,
    refGenome = "Gallus_gallus_REF",
    genomeIDs = genomeIDs,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    reorderBySynteny = FALSE
)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc.roi.2.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)




roi <- data.frame( 
    genome = c("Strix_aluco"), 
    chr = c("Z"), 
    color = c("#c23d3d")) 


ripDat <- plot_riparian(
    gsParam = gsParam,
    highlightBed = roi, 
    backgroundColor = NULL,
    refGenome = "Gallus_gallus_REF",
    genomeIDs = genomeIDs,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    reorderBySynteny = FALSE
)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc.roi.3.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)
```
# Analyze within clade, excluding chicken

Set up environment and directories
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/neo_sex/Morphnus_etc_2/genespace

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/neo_sex/Morphnus_etc_2/genespace

cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Morphnus_guianensis/bed .
cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Morphnus_guianensis/peptide .

cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Accipiter_gentilis/bed .
cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Accipiter_gentilis/peptide .

cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Strix_aluco/bed .
cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Strix_aluco/peptide .

cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Haliaeetus_albicilla/bed .
cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Haliaeetus_albicilla/peptide .

cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Harpia_harpyja/bed .
cp -r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/Harpia_harpyja/peptide .

```


Morphnus_guianensis
Accipiter_gentilis
Haliaeetus_albicilla
Harpia_harpyja
Strix_alco
```
#!/bin/bash

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace

source myconda

mamba activate genespace

Rscript /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/neo_sex/Morphnus_etc_2/Morphnus_etc_2.genespace.r
```
## Rscript:
```
library(GENESPACE)
library(ggplot2)
library(data.table)

setwd("/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/neo_sex/Morphnus_etc_2/")

gpar <- init_genespace(
  wd = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/neo_sex/Morphnus_etc_2/genespace",
  path2mcscanx = "/vf/users/jacksondan/conda/envs/genespace/bin/"
)

out <- run_genespace(gpar, overwrite = FALSE)

# Genome IDs to plot
genomeIDs <- c("Morphnus_guianensis", "Accipiter_gentilis", "Strix_aluco", "Haliaeetus_albicilla", "Harpia_harpyja")

ordFun <- function(x) {
  data.table::frank(
    list(grepl("[a-zA-Z]", x),
         suppressWarnings(as.numeric(gsub("[a-zA-Z]", "", x)))),
    ties.method = "random"
  )
}

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
)


ripDat <- plot_riparian(
  gsParam = out,
  refGenome = "Strix_aluco",
  genomeIDs = genomeIDs,
  refChrOrdFun = ordFun,
  addThemes = ggthemes,
  chrFill = "lightgrey",
  reorderBySynteny = FALSE
)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc_2.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)
```
## Sbatch submission:
```
sbatch \
  -c 16 \
  -t 1:00:00 \
  --mem-per-cpu=20G \
  --mail-user=jacksondan@nih.gov \
  --mail-type=ALL \
  --gres=lscratch:500 \
  --output=slurm_output/genespace_morphnus_2.%j \
  submit_morphnus_2.sh
```
## Redoing the plot

```
library(GENESPACE)
library(ggplot2)
library(data.table)

load('/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/neo_sex/Morphnus_etc_2/genespace/results/gsParams.rda',
        verbose = TRUE)

# Genome IDs to plot

genomeIDs <- c("Morphnus_guianensis", "Accipiter_gentilis", "Strix_aluco", "Haliaeetus_albicilla", "Harpia_harpyja")

ordFun <- function(x) {
  data.table::frank(
    list(grepl("[a-zA-Z]", x),
         suppressWarnings(as.numeric(gsub("[a-zA-Z]", "", x)))),
    ties.method = "random"
  )
}

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
)


ripDat <- plot_riparian(
  gsParam = gsParam,
  refGenome = "Strix_aluco",
  genomeIDs = genomeIDs,
  refChrOrdFun = ordFun,
  addThemes = ggthemes,
  chrFill = "lightgrey",
  reorderBySynteny = FALSE
)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc_2.redo.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

roi <- data.frame( 
    genome = c("Accipiter_gentilis"), 
    chr = c("Z"), 
    color = c("#c23d3d")) 

ripDat <- plot_riparian(
    gsParam = gsParam,
    highlightBed = roi, 
    backgroundColor = NULL,
    refGenome = "Strix_aluco",
    genomeIDs = genomeIDs,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    reorderBySynteny = FALSE
)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc.Accipiter.roi.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)


roi <- data.frame( 
    genome = c("Morphnus_guianensis"), 
    chr = c("CM098430.1"), 
    color = c("#c23d3d")) 


ripDat <- plot_riparian(
    gsParam = gsParam,
    highlightBed = roi, 
    backgroundColor = NULL,
    refGenome = "Strix_aluco",
    genomeIDs = genomeIDs,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    reorderBySynteny = FALSE
)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc.roi.Morphnus.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)


roi <- data.frame( 
    genome = c("Strix_aluco"), 
    chr = c("Z"), 
    color = c("#c23d3d")) 


ripDat <- plot_riparian(
    gsParam = gsParam,
    highlightBed = roi, 
    backgroundColor = NULL,
    refGenome = "Strix_aluco",
    genomeIDs = genomeIDs,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    reorderBySynteny = FALSE
)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc.roi.Strix_aluco.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

# Add Gallus gallus

genomeIDs <- c("Morphnus_guianensis", "Accipiter_gentilis", "Strix_aluco", "Haliaeetus_albicilla", "Harpia_harpyja", "Gallus_gallus_REF")

ordFun <- function(x) {
  data.table::frank(
    list(grepl("[a-zA-Z]", x),
         suppressWarnings(as.numeric(gsub("[a-zA-Z]", "", x)))),
    ties.method = "random"
  )
}

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
)


ripDat <- plot_riparian(
  gsParam = gsParam,
  refGenome = "Gallus_gallus_REF",
  genomeIDs = genomeIDs,
  refChrOrdFun = ordFun,
  addThemes = ggthemes,
  chrFill = "lightgrey",
  reorderBySynteny = FALSE
)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc_2.redo.Gallus_gallus_REF.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

roi <- data.frame( 
    genome = c("Accipiter_gentilis"), 
    chr = c("Z"), 
    color = c("#c23d3d")) 

ripDat <- plot_riparian(
    gsParam = gsParam,
    highlightBed = roi, 
    backgroundColor = NULL,
    refGenome = "Strix_aluco",
    genomeIDs = genomeIDs,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    reorderBySynteny = FALSE
)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc.Accipiter.Gallus_gallus_REF.roi.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)


roi <- data.frame( 
    genome = c("Morphnus_guianensis"), 
    chr = c("CM098430.1"), 
    color = c("#c23d3d")) 


ripDat <- plot_riparian(
    gsParam = gsParam,
    highlightBed = roi, 
    backgroundColor = NULL,
    refGenome = "Strix_aluco",
    genomeIDs = genomeIDs,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    reorderBySynteny = FALSE
)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc.roi.Morphnus.Gallus_gallus_REF.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)


roi <- data.frame( 
    genome = c("Strix_aluco"), 
    chr = c("Z"), 
    color = c("#c23d3d")) 


ripDat <- plot_riparian(
    gsParam = gsParam,
    highlightBed = roi, 
    backgroundColor = NULL,
    refGenome = "Strix_aluco",
    genomeIDs = genomeIDs,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    reorderBySynteny = FALSE
)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc.roi.Strix_aluco.Gallus_gallus_REF.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

# Harpia_harpyja

roi <- data.frame( 
    genome = c("Harpia_harpyja"), 
    chr = c("Z"), 
    color = c("#c23d3d")) 


ripDat <- plot_riparian(
    gsParam = gsParam,
    highlightBed = roi, 
    backgroundColor = NULL,
    refGenome = "Strix_aluco",
    genomeIDs = genomeIDs,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    reorderBySynteny = FALSE
)


p_list <- ripDat$plot
p <- p_list[[1]]

ggsave(
  filename = "Morphnus_etc.roi.Harpia_harpyja.Gallus_gallus_REF.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)


# Paired plots
genomeIDs <- c("Morphnus_guianensis")


ordFun <- function(x) {
  data.table::frank(
    list(grepl("[a-zA-Z]", x),
         suppressWarnings(as.numeric(gsub("[a-zA-Z]", "", x)))),
    ties.method = "random"
  )
}

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
)


ripDat <- plot_riparian(
  gsParam = gsParam,
  refGenome = "Gallus_gallus_REF",
  genomeIDs = genomeIDs,
  refChrOrdFun = ordFun,
  addThemes = ggthemes,
  chrFill = "lightgrey",
  reorderBySynteny = FALSE
)

ggsave(
  filename = "Morphnus_guianensis.VS.Gallus_gallus_REF.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

genomeIDs <- c("Accipiter_gentilis", "Gallus_gallus_REF")

ordFun <- function(x) {
  data.table::frank(
    list(grepl("[a-zA-Z]", x),
         suppressWarnings(as.numeric(gsub("[a-zA-Z]", "", x)))),
    ties.method = "random"
  )
}

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
)


ripDat <- plot_riparian(
  gsParam = gsParam,
  refGenome = "Gallus_gallus_REF",
  genomeIDs = genomeIDs,
  refChrOrdFun = ordFun,
  addThemes = ggthemes,
  chrFill = "lightgrey",
  reorderBySynteny = FALSE
)

ggsave(
  filename = "Accipiter_gentilis.VS.Gallus_gallus_REF.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

genomeIDs <- c("Strix_aluco", "Gallus_gallus_REF")

ordFun <- function(x) {
  data.table::frank(
    list(grepl("[a-zA-Z]", x),
         suppressWarnings(as.numeric(gsub("[a-zA-Z]", "", x)))),
    ties.method = "random"
  )
}

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
)


ripDat <- plot_riparian(
  gsParam = gsParam,
  refGenome = "Gallus_gallus_REF",
  genomeIDs = genomeIDs,
  refChrOrdFun = ordFun,
  addThemes = ggthemes,
  chrFill = "lightgrey",
  reorderBySynteny = FALSE
)

ggsave(
  filename = "Strix_aluco.VS.Gallus_gallus_REF.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

genomeIDs <- c("Haliaeetus_albicilla", "Gallus_gallus_REF")

ordFun <- function(x) {
  data.table::frank(
    list(grepl("[a-zA-Z]", x),
         suppressWarnings(as.numeric(gsub("[a-zA-Z]", "", x)))),
    ties.method = "random"
  )
}

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
)


ripDat <- plot_riparian(
  gsParam = gsParam,
  refGenome = "Gallus_gallus_REF",
  genomeIDs = genomeIDs,
  refChrOrdFun = ordFun,
  addThemes = ggthemes,
  chrFill = "lightgrey",
  reorderBySynteny = FALSE
)

ggsave(
  filename = "Haliaeetus_albicilla.VS.Gallus_gallus_REF.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)

genomeIDs <- c("Harpia_harpyja", "Gallus_gallus_REF")

ordFun <- function(x) {
  data.table::frank(
    list(grepl("[a-zA-Z]", x),
         suppressWarnings(as.numeric(gsub("[a-zA-Z]", "", x)))),
    ties.method = "random"
  )
}

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white")
)


ripDat <- plot_riparian(
  gsParam = gsParam,
  refGenome = "Gallus_gallus_REF",
  genomeIDs = genomeIDs,
  refChrOrdFun = ordFun,
  addThemes = ggthemes,
  chrFill = "lightgrey",
  reorderBySynteny = FALSE
)

ggsave(
  filename = "Harpia_harpyja.VS.Gallus_gallus_REF.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in"
)


```