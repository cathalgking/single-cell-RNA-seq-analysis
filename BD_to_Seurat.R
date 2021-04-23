#########################################################################################################################
# Take a BD Rhapsody generated single-cell dataset and preprocess it to a Seurat object for analysis with Seurat
# This assumes the counts file was generated with SevenBridges and was aligned with the STAR aligner
# read in the RSEC_MolsPerCell.csv file from SBG 
# Find a sample rhapsody dataset here  https://drive.google.com/file/d/1q8KGAey6tajeJPZdfryGVVmcW-Fwq14U/view?usp=sharing
# read in genes, barcodes, matrix files if analysing 10x data
# Seurat manual https://cran.r-project.org/web/packages/Seurat/Seurat.pdf
#########################################################################################################################


library(Seurat)
library(tidyverse)
library(MAST)
library(RColorBrewer)
library(doMC)
library(data.table)
library(reticulate)


# Use the RSEC Mol per cell csv file from SBG output. Read it in but notice 'skip=8' as the first 8 lines just contain info about the run.
Abseq_1 <- read.csv(file = '~/Downloads/Combined_BD-Demo-WTA-AbSeq-SMK_RSEC_MolsPerCell.csv', sep = ',', header = TRUE, row.names = 1, check.names = F, skip = 8)
# transpose count matrix 
Abseq_1P <- t(Abseq_1[, str_detect(string = colnames(Abseq_1), pattern = 'pAbO')])

# look at part of the matrix 
Abseq_1P[1:10, 1:15]
# get dimensions - cols and rows
dim(Abseq_1P)

# use regex to seperate out identifiers (genes) contained in same column based on the pAb0 identifier
Abseq_1RNA <- (Abseq_1[, !str_detect(string = colnames(Abseq_1), pattern = 'pAbO')])
# transpose to be left with a gene by cell count matrix
Abseq_1RNA <- t(Abseq_1RNA[, !str_detect(string = colnames(Abseq_1RNA), pattern = 'Sample')])

# look again; sanity check
Abseq_1RNA[1:10, 1:15]
dim(Abseq_1RNA)

# sample tag identifiers
Sample_Tag1 <- t(Abseq_1[, str_detect(string = colnames(Abseq_1), pattern = 'Sample_Name')]) # 

#- Features - Proteins. Extract protein names, again by pattern recognition.
P_names <- sapply(X = str_split(string = rownames(Abseq_1P), pattern = '\\|'), 
                  FUN = function(x) ifelse(x[1] == 'CD197', 
                                           paste(paste(x[1], x[2], sep = '|'), x[3], sep = '_'),
                                           paste(x[1], x[2], sep = '|')))

rownames(Abseq_1P) <- paste0(P_names, '_P')

RNA_names <- str_replace(string = rownames(Abseq_1RNA), pattern = '\\|[^(PolyA)]*',  replacement = '_')
RNA_names <- str_replace(string = RNA_names, pattern = '_?$', replacement = '_RNA')
rownames(Abseq_1RNA) <- RNA_names

#- Have a look at dimensions and the rownames of the RNA and antibody matrix
dim(Abseq_1RNA)
dim(Abseq_1P)

# Finally, create Seurat object
Cartridge1 <- CreateSeuratObject(counts = Abseq_1RNA, min.cells = 1, project = "BD_SBG_data")

### merge multiple experiments with MergeSeurat()
### VlnPlot(object = Cartridge1, features = RNA_names)

## Add protein expression data to the Seurat object
# Make sure that the Cell IDs have the same names
colnames(Abseq_1P) <- paste("Cartridge1", colnames(Abseq_1P), sep="_")


#########################################################################################################################################################
# From here, the rhapsody dataset is in a Seurat Object, contains the relevant information and is ready to be analysed in Seurat in whatever way you wish.
# Begin with usual steps; Pre-processing, QC, normalisation etc.
#########################################################################################################################################################



PBMC <- NormalizeData(object = Cartridge1, assay.type = 'RNA', normalization.method = 'LogNormalize', scale.factor = 10000)

PBMC[["percent.mt"]] <- PercentageFeatureSet(object = PBMC, pattern = "^MT-")

head(x = PBMC@meta.data, 5)


#Visualize QC metrics as a violin plot
VlnPlot(object = PBMC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualize feature-feature relationships
plot1 <- FeatureScatter(object = PBMC, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(object = PBMC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
CombinePlots(plots = list(plot1,plot2))
PBMC <- subset(x = PBMC, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


## Identification of highly variable features
PBMC <- FindVariableFeatures(object = PBMC,selection.method = 'vst', nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = PBMC), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = PBMC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))





