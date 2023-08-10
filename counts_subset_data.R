## subset object to only contain cells from a certain Ident in the Seurat object. For example below, subset the Seurat object to only contain cells labelled as Monocytes.

#
library(Seurat)
library(SeuratObject)

## check dimplot
# DimPlot(seu_no_PCs, reduction = "umap.harmony", group.by = "seurat_clusters", label = T)

# extract assay from Seurat object
counts <- as.data.frame(seu_no_PCs@assays$RNA$counts)
# View
View(counts)

# check current Idents
table(Idents(seu_no_PCs))
levels(seu_no_PCs)

## set Idents to the variable that you want to subset on e.g cluster number, cell label etc.
# set Idents to clusters number in meta-data
Idents(seu_no_PCs) <- "seurat_clusters"

## use subset to reduce object based on a condition
seu_cluster5 <- subset(seu_no_PCs, idents = 5) # subset object to only contain cells from cluster 5

# set Idents to cell annotation variable in meta-data
Idents(seu_no_PCs) <- "predicted.id"
seu_monocytes <- subset(seu_no_PCs, idents = "monocyte") # subset object to only contain cells labelled as monocyte

###

## split object

annotations <- SplitObject(object = All_merged_noPCs_none, split.by = "predicted.id")


