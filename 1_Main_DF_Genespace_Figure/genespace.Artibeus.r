library(GENESPACE)
library(ggplot2)

wd <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/figure"

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = "/vf/users/jacksondan/conda/envs/genespace/bin/")

out <- run_genespace(gpar, overwrite = F)

# load('/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Artibeus/figure/results/gsParams.rda', verbose = TRUE)


SPECIES=c("Homo_sapiens", "Artibeus_lituratus", "Artibeus_intermedius") 

roi <- data.frame( 
    genome = c("Homo_sapiens", "Homo_sapiens"), 
    chr = c("X", "Y"), 
    color = c("#7f00e6", "#045b1b")) 

ordFun <- function(x) {
  x2 <- gsub("^chr", "", x, ignore.case=TRUE)
  sexmap <- c(X=23, Y=24, M=25, MT=25)
  num <- suppressWarnings(as.numeric(x2))
  num[is.na(num) & x2 %in% names(sexmap)] <- sexmap[x2[is.na(num) & x2 %in% names(sexmap)]]
  num[is.na(num)] <- 1e9
  rank(num, ties.method="first")
}


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
    backgroundColor = "#7b7b7b", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Artibeus.v4.10x6.pdf", 
    plot = p, 
    width = 10, 
    height = 6, 
    units = "in" )



# Plot 2

SPECIES=c("Homo_sapiens", "Artibeus_lituratus") 

roi <- data.frame( 
    genome = c("Artibeus_lituratus", "Artibeus_lituratus", "Artibeus_lituratus"), 
    chr = c("X", "Y", "Y2"), 
    color = c("#7f00e6", "#045b1b", "#c23d3d")) 

ordFun <- function(x) {
  x2 <- gsub("^chr", "", x, ignore.case=TRUE)
  sexmap <- c(X=23, Y=24, M=25, MT=25)
  num <- suppressWarnings(as.numeric(x2))
  num[is.na(num) & x2 %in% names(sexmap)] <- sexmap[x2[is.na(num) & x2 %in% names(sexmap)]]
  num[is.na(num)] <- 1e9
  rank(num, ties.method="first")
}


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
    backgroundColor = "#7b7b7b", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Artibeus2.v4.10x3.pdf", 
    plot = p, 
    width = 10, 
    height = 3, 
    units = "in" )


# Plot 3

SPECIES=c("Homo_sapiens", "Artibeus_intermedius") 

roi <- data.frame( 
    genome = c("Artibeus_intermedius", "Artibeus_intermedius", "Artibeus_intermedius"), 
    chr = c("X", "Y1", "Y2"), 
    color = c("#7f00e6", "#045b1b", "#c23d3d")) 

ordFun <- function(x) {
  x2 <- gsub("^chr", "", x, ignore.case=TRUE)
  sexmap <- c(X=23, Y=24, M=25, MT=25)
  num <- suppressWarnings(as.numeric(x2))
  num[is.na(num) & x2 %in% names(sexmap)] <- sexmap[x2[is.na(num) & x2 %in% names(sexmap)]]
  num[is.na(num)] <- 1e9
  rank(num, ties.method="first")
}


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
    backgroundColor = "#7b7b7b", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Artibeus3.v4.10x3.pdf", 
    plot = p, 
    width = 10, 
    height = 3, 
    units = "in" )



# Plot 4

SPECIES=c("Artibeus_lituratus", "Artibeus_intermedius") 

roi <- data.frame( 
    genome = c("Artibeus_lituratus", "Artibeus_lituratus", "Artibeus_lituratus"), 
    chr = c("X", "Y", "Y2"), 
    color = c("#7f00e6", "#045b1b", "#c23d3d")) 

ordFun <- function(x) {
  x2 <- gsub("^chr", "", x, ignore.case=TRUE)
  sexmap <- c(X=23, Y=24, M=25, MT=25)
  num <- suppressWarnings(as.numeric(x2))
  num[is.na(num) & x2 %in% names(sexmap)] <- sexmap[x2[is.na(num) & x2 %in% names(sexmap)]]
  num[is.na(num)] <- 1e9
  rank(num, ties.method="first")
}


ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Artibeus_lituratus", 
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = "#7b7b7b", 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Artibeus_spp.v4.10x3.pdf", 
    plot = p, 
    width = 10, 
    height = 3, 
    units = "in" )


invchr <- data.frame(
  genome = c("Artibeus_lituratus", "Artibeus_intermedius"), 
  chr = c("X", "Y2"))

ripDat <- plot_riparian( 
    gsParam = out, 
    highlightBed = roi, 
    refGenome = "Artibeus_lituratus", 
    invertTheseChrs = invchr,
    genomeIDs = SPECIES, 
    refChrOrdFun = ordFun,
    addThemes = ggthemes,
    chrFill = "lightgrey",
    backgroundColor = NULL, 
    reorderBySynteny = FALSE)

p_list <- ripDat$plot
p <- p_list[[1]]

ggsave( 
    filename = "Artibeus_spp.v4.sex_chrs.10x3.pdf", 
    plot = p, 
    width = 5, 
    height = 3, 
    units = "in" )