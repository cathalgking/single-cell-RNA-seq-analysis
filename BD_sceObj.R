# Create a SingleCellExperiment (sce) object from a BD Rhapsody generated dataset
# Use the RSEC corrected reads per cell counts data table outputted from the BD Targeted analysis pipeline on SevenBridges Genomics
# The sce object is of S4 class and can be converted to a Seurat object and used with a variety of packages

library(SingleCellExperiment)

# read in data table
counts_matrix <- read.table("/home/cathal/Main_RSEC_and_SampleTagReads_perCell_merged.txt", header = T, sep = '\t')

# now construct our first SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = counts_matrix))

# access the count data assay in the sce object 
assay(sce, "counts")
# assay(test.fromseurat2SCE, "counts")

# see what assays are in the object
assays(sce)
