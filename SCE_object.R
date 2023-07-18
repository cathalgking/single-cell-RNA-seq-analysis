## convert SCE object back to a Seurat object or use the SCE object for analysis.

## load the following packages
library(Seurat)
library(SingleCellExperiment)

## read in SCE object with rds
stromal_sce <- readRDS(file = "/Users/cathal.king/Desktop/Stromal_cell_cluster_SCE.rds")

## convert the SCE object to a Seurat object
stromal_cluster_seurat <- as.Seurat(stromal_sce)

## if that works then great. Use the Seurat object.

## double check that both have the same number of cells and genes
dim(stromal_cluster_seurat) # seurat
dim(stromal_sce)

## dimplot with the Seurat object. It just seems to have capitalised the names of the dims.
DimPlot(stromal_seurat, reduction = "UMAP.HARMONY")




## The SCE object is good too tho and you could do some analysis or plots something like the following.

## first you will need these packages
library(scater)
library(scuttle)

## check number of genes and cells in SCE object
dim(stromal_sce)

## the meta-data is in the coldata slot of a SCE object
head(colData(stromal_sce), 2)

## dimensionality reduction plot with a SCE object
plotReducedDim(object = stromal_sce, dimred = "UMAP.HARMONY")




