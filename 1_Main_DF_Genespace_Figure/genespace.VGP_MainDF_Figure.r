library(GENESPACE)
library(ggplot2)

wd <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure"

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = "/vf/users/jacksondan/conda/envs/genespace/bin/")

out <- run_genespace(gpar, overwrite = F)



SPECIES=c("Homo_sapiens", "Gallus_gallus", "Anolis_sagrei", "Podarcis_raffonei", "Pseudacris_triseriata", "Hoplias_malabaricus", "Narcine_bancroftii") 

SPECIES=c("Narcine_bancroftii", "Hoplias_malabaricus", "Pseudacris_triseriata", "Podarcis_raffonei", "Anolis_sagrei", "Gallus_gallus", "Homo_sapiens") 

roi <- data.frame( 
    genome = c("Homo_sapiens", "Gallus_gallus", "Anolis_sagrei", "Podarcis_raffonei", "Pseudacris_triseriata", "Hoplias_malabaricus", "Hoplias_malabaricus", "Narcine_bancroftii"), 
    chr = c("X", "Z", "X", "Z", "X", "X1", "X2", "X"), 
    color = c("#332288", "#88CCEE", "#44AA99", "#999933", "#000000", "#CC6677", "#CC6677", "#882255")) 

ordFun <- function(x) 
  data.table::frank(
    list(grepl("[a-zA-Z]", x), 
         as.numeric(gsub("[a-zA-Z]", "", x))), 
    ties.method = "random")

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Homo_sapiens", 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ffffff", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "VGP_MainDF_Figure.10x6.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )

ggsave( 
    filename = "VGP_MainDF_Figure.10x4.pdf", 
    plot = p, 
    width = 10, 
    height = 4, 
    units = "in" )

ggsave( 
    filename = "VGP_MainDF_Figure.10x3.pdf", 
    plot = p, 
    width = 10, 
    height = 3, 
    units = "in" )


ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Homo_sapiens", 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ffffff", 
    scaleGapSize = .75,
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "VGP_MainDF_Figure.10x6.gapProp.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )

ggsave( 
    filename = "VGP_MainDF_Figure.10x4.gapProp.pdf", 
    plot = p, 
    width = 10, 
    height = 4, 
    units = "in" )

ggsave( 
    filename = "VGP_MainDF_Figure.10x3.gapProp.pdf", 
    plot = p, 
    width = 10, 
    height = 3, 
    units = "in" )

## light background synteny

ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Homo_sapiens", 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    scaleGapSize = .75,
    backgroundColor = "#ececec", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "VGP_MainDF_Figure.10x6.gapProp.lightBackground.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )

## reversed order of coloring

roi <- data.frame( 
    genome = c("Homo_sapiens", "Gallus_gallus", "Anolis_sagrei", "Podarcis_raffonei", "Pseudacris_triseriata", "Hoplias_malabaricus", "Hoplias_malabaricus", "Narcine_bancroftii"), 
    chr = c("X", "Z", "X", "Z", "X", "X1", "X2", "X"), 
    color = c("#332288", "#88CCEE", "#44AA99", "#999933", "#000000", "#CC6677", "#CC6677", "#882255")) 


roi <- roi[rev(seq_len(nrow(roi))), ]


ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Homo_sapiens", 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    scaleGapSize = .75,
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "VGP_MainDF_Figure.10x6.gapProp.reverseColors.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )

# Highlight sex chr of only one species at a time
## Pseudacris_triseriata
roi <- data.frame( 
    genome = c("Pseudacris_triseriata"), 
    chr = c("1"), 
    color = c("#000000")) 


ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Homo_sapiens", 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    scaleGapSize = .75,
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "VGP_MainDF_Figure.10x6.gapProp.Pseudacris_triseriata.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )


roi <- data.frame( 
    genome = c("Homo_sapiens"), 
    chr = c("X"), 
    color = c("#332288")) 

ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Homo_sapiens", 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    scaleGapSize = .75,
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "VGP_MainDF_Figure.10x6.gapProp.Homo_sapiens.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )


# Gallus gallus

roi <- data.frame( 
    genome = c("Gallus_gallus"), 
    chr = c("Z"), 
    color = c("#88CCEE")) 

ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Gallus_gallus", 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    scaleGapSize = .75,
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "VGP_MainDF_Figure.10x6.gapProp.Gallus_gallus.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )

# Anolis_sagrei

roi <- data.frame( 
    genome = c("Anolis_sagrei"), 
    chr = c("X"), 
    color = c("#44AA99")) 

ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Anolis_sagrei", 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    scaleGapSize = .75,
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "VGP_MainDF_Figure.10x6.gapProp.Anolis_sagrei.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )

# Podarcis_raffonei

roi <- data.frame( 
    genome = c("Podarcis_raffonei"), 
    chr = c("Z"), 
    color = c("#999933")) 

ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Podarcis_raffonei", 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    scaleGapSize = .75,
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "VGP_MainDF_Figure.10x6.gapProp.Podarcis_raffonei.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )

# Hoplias_malabaricus

roi <- data.frame( 
    genome = c("Hoplias_malabaricus"), 
    chr = c( "X1", "X2"), 
    color = c("#CC6677", "#7329c6")) 

ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Hoplias_malabaricus", 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    scaleGapSize = .75,
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "VGP_MainDF_Figure.10x6.gapProp.Hoplias_malabaricus.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )

# Hoplias X1

roi <- data.frame( 
    genome = c("Hoplias_malabaricus"), 
    chr = c( "X1"), 
    color = c("#CC6677")) 

ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Hoplias_malabaricus", 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    scaleGapSize = .75,
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "VGP_MainDF_Figure.10x6.gapProp.Hoplias_malabaricus.X1.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )

# Hoplias X2

roi <- data.frame( 
    genome = c("Hoplias_malabaricus"), 
    chr = c( "X2"), 
    color = c("#CC6677")) 

ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Hoplias_malabaricus", 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    scaleGapSize = .75,
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "VGP_MainDF_Figure.10x6.gapProp.Hoplias_malabaricus.X2.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )

# Narcine_bancroftii

roi <- data.frame( 
    genome = c("Narcine_bancroftii"), 
    chr = c( "12"), 
    color = c("#882255")) 

ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Narcine_bancroftii", 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ececec", 
    scaleGapSize = .75,
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "VGP_MainDF_Figure.10x6.gapProp.Narcine_bancroftii.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )

