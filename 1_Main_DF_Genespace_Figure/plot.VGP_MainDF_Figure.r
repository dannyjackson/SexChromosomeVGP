# Plotting Genespace

```
R

library(GENESPACE)
library(ggplot2)
wd <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure"

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = "/vf/users/jacksondan/conda/envs/genespace/bin/")

out <- run_genespace(gpar, overwrite = F)



SPECIES=c("Homo_sapiens", "Gallus_gallus", "Anolis_sagrei", "Podarcis_raffonei", "Pseudacris_triseriata", "Hoplias_malabaricus", "Narcine_bancroftii") 

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