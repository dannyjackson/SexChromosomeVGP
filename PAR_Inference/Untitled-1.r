cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_inference/alignment/SV_by_eye

## Install required packages
install.packages("devtools")
library(devtools)

## Install from GitHub repository
devtools::install_github("daewoooo/SVbyEye", branch="master")

/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/minimap2

## Load the SVbyEye package
library(SVbyEye)
## Get PAF file to read
paf.file <- "/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/minimap2/Eudromia_elegans_WtoZ.aln.paf"
## Read in PAF
paf.table <- readPaf(
    paf.file = paf.file,
    include.paf.tags = TRUE, restrict.paf.tags = "cg"
)
paf.table

filterPaf(paf.table = paf.table, min.align.len = 100000)

flipPaf(paf.table = paf.table, force = TRUE)

paf.file <- system.file("extdata", "test1.paf",
    package = "SVbyEye"
)
## Read in PAF
paf.table <- readPaf(
    paf.file = paf.file,
    include.paf.tags = TRUE, restrict.paf.tags = "cg"
)
## Make a plot colored by alignment directionality
png(
  filename = "plotMiro_direction.png",
  width = 2000,
  height = 2000,
  res = 300
)

plotMiro(
  paf.table = paf.table,
  color.by = "direction"
)

dev.off()
