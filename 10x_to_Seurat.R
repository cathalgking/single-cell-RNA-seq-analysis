# import 10x single-cell data to a Seurat object
# some other functions to check data


library(Seurat)
library(tidyverse)
library(MAST)
library(RColorBrewer)
library(data.table)
library(reticulate)

# single-cell RNA-seq for the enterocytes https://zenodo.org/record/3403670#.XrF8S9NKjyU
# /home/cathal_king/Downloads/raw_data/GSM2644349_Lgr5eGFP_neg_1

# Read in 10x data folder
ent.data <- Read10X(data.dir = "/home/cathal_king/Downloads/raw_data/GSM2644349_Lgr5eGFP_neg_1")
# Initialize the Seurat object with the raw (non-normalized data).
ent <- CreateSeuratObject(counts = ent.data, project = "enterocytes_continuum", min.cells = 3, min.features = 200)
ent

dim(ent)
ent.data[1:5, 1:7]


dense.size <- object.size(x = as.matrix(x = ent.data))
dense.size
sparse.size <- object.size(x = ent.data)
sparse.size
dense.size / sparse.size
