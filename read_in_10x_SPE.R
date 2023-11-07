### to a SPE (SpatialExperiment) object
library(SpatialExperiment)
library(ggspavis)
# path must contain a 10x "/outs/" folder
sample <- file.path("/PATH/outs/")
# read in
C1_sce <- read10xVisium(samples = sample, type = "sparse", images = "lowres")
C1_sce
# View image in SPE object
plotVisium(C1_sce)


### to a Seurat object
# Load a 10x Genomics Visium Spatial Experiment into a Seurat object
seu_pros_data <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir = "/PATH/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  slice = "slice1", # specify name of the stored image
  filter.matrix = TRUE, 
  to.upper = FALSE
)
#
seu_pros_data
# plot
SpatialPlot(object = seu_pros_data)

