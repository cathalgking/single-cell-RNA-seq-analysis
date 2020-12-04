###########################################################
# pre-processing of BD Rhapsody single-cell data (.st file)
###########################################################

library(Seurat)
library(deMULTIplex)
library(DoubletFinder)

## Step 1: Load raw gene expression (long-format .st file from SevenBridgesGenomics output) -----------------------------------------------------------------------------------------------------------
rna.BD.raw <- read.table("/home/cathal_king/Desktop/BD-Demo-WTA-AbSeq-SMK/Combined_BD-Demo-WTA-AbSeq-SMK_Expression_Data.st", header=T, stringsAsFactors = F)
rna.BD.raw <- rna.BD.raw[,c(1,2,7)]
colnames(rna.BD.raw) <- c("Cell","Gene","Count")
cellIDs <- unique(rna.BD.raw$Cell)
geneIDs <- unique(rna.BD.raw$Gene)

## Step 2: Convert to sparse, wide matrix) ------------------------------------------------------------------------------------------------------------------------------------------------------------
library(Matrix)
raw.mat.BD <- matrix(0L, nrow=length(geneIDs), ncol=length(cellIDs))
raw.mat.BD <- Matrix(raw.mat.BD, sparse = TRUE)
rownames(raw.mat.BD) <- geneIDs; colnames(raw.mat.BD) <- cellIDs
x <- 0
for (i in cellIDs) {
  ind <- which(rna.BD.raw$Cell == i)
  x <- x + 1
  print(paste(x, "out of 5665", sep=" "))
  raw.mat.BD[rna.BD.raw$Gene[ind],x] <- rna.BD.raw$Count[ind]
  rna.BD.raw <- rna.BD.raw[-ind, ]
}


## Step 3: Remove uninformative genes and CITE-seq Ab counts) -----------------------------------------------------------------------------------------------------------------------------------------
gene.counts <- Matrix::rowSums(raw.mat.BD)
raw.mat.BD <- raw.mat.BD[which(gene.counts > 3), ]
ind <- grep("pAb",rownames(raw.mat.BD))
raw.mat.BD <- raw.mat.BD[-ind, ]
seu_bd <- CreateSeuratObject(raw.mat.BD)

# convert Seurat obj to SCE object
test.fromseurat2SCE <- as.SingleCellExperiment(seu_bd)
# or continue with seu_bd object with Seurat functions


# If you are starting here with another count matrix, just construct the SCE object (SingleCellExperiment) and proceed
# sce <- SingleCellExperiment(assays = list(counts = counts_matrix))

# access the count data assay in the sce object 
assay(sce, "counts")
assay(test.fromseurat2SCE, "counts")



# log normalise using scater
sce <- scater::logNormCounts(sce)
test.fromseurat2SCE <- scater::logNormCounts(test.fromseurat2SCE)


# look at the logcount matrix
logcounts(sce)

# see what assays are in the object
assays(sce)


# add some meta-data or information about the columns or samples or cells
cell_metadata <- data.frame(batch = c(1, 1, 2))
rownames(cell_metadata) <- paste0("cell_", 1:3)

# now make a new sce object with the meta data in it
sce <- SingleCellExperiment(assays = list(counts = counts_matrix),
                            colData = cell_metadata)

# access the col data 
colData(sce)
# alternatively
sce$batch


# if we only want cells in batch 1 for example
sce[, sce$batch == 1]

# rowRanges tells you the genomic intervals or ranges
rowRanges(sce)


# add QC metrics for rows (Genes)
sce <- scater::addPerFeatureQC(sce)
rowData(sce)

# subset certain genes into their own sce object
gene1_4 <- sce[c("gene_1", "gene_4"), ]
# alternatively
#gene1_4 <- sce[c(1, 4), ]


## IF you have a certain list of genes that you want to look at then input them like so
# they will be put into the metadata slot here
my_genes <- c("gene_1", "gene_5")
metadata(sce) <- list(favorite_genes = my_genes)
metadata(sce)


# run PCA and add PCA to the reducedDim slot
sce <- scater::logNormCounts(sce)
sce <- scater::runPCA(sce)
reducedDim(sce, "PCA")


# run a tNSE analysis
sce <- scater::runTSNE(sce, perplexity = 0.1)
# look at the results
reducedDim(sce, "TSNE")

# see all entries in 
reducedDims(sce)

# run UMAP using the uwot package and add it to the reducedDim slot. 
u <- uwot::umap(t(logcounts(sce)), n_neighbors = 2)
reducedDim(sce, "UMAP_uwot") <- u
reducedDims(sce) # Now stored in the object.
# can also do it with the scater package, and it will add automatically to the slot







