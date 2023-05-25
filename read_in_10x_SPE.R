#
library(SpatialExperiment)
library(ggspavis)
# path must contain a 10x "/outs/" folder
sample <- file.path("/Users/cathal.king/Documents/SAHMRI/Visium_Analysis/Outs/C1_outs/outs/")
# read in
C1_sce <- read10xVisium(samples = sample, type = "sparse", images = "lowres")
C1_sce
# View image in SPE object
plotVisium(C1_sce)