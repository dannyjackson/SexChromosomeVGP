library(GENESPACE)
library(ggplot2)

load('/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure_GWS//results/gsParams.rda')

wd <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure_GWS/"

SPECIES=c("Carcharodon_carcharias", "Gasterosteus_aculeatus", "Pseudacris_triseriata", "Podarcis_raffonei", "Gallus_gallus", "Homo_sapiens") 

roi <- data.frame( 
    genome = c("Homo_sapiens", "Gallus_gallus", "Podarcis_raffonei", "Pseudacris_triseriata", "Gasterosteus_aculeatus", "Carcharodon_carcharias"), 
    chr = c("X", "Z", "Z", "X", "X", "X"), 
    color = c("#E69F00", "#00796B", "#80CBC4", "#984EA3", "#56B4E9", "#0072B2")) 

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
    refGenome = "Homo_sapiens", 
    genomeIDs = SPECIES, 
    minChrLen2plot = 1,
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#ffffff", 
    scaleGapSize = .75,
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "VGP_MainDF_Figure.10x6.gapProp.GWS.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )