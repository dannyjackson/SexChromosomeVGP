library(GENESPACE)
library(ggplot2)

wd <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/anuran_plot"

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = "/vf/users/jacksondan/conda/envs/genespace/bin/")

out <- run_genespace(gpar, overwrite = F)


SPECIES=c("Pseudacris_triseriata", "Hyla_sarda") 

roi <- data.frame( 
    genome = c("Pseudacris_triseriata"),
    chr = c("X"), 
    color = c("#000000")) 

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
    refGenome = "Pseudacris_triseriata", 
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
    filename = "Anurans.10x6.gapProp.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )
