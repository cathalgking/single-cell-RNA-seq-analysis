# Create a SingleCellExperiment (sce) object from a BD Rhapsody generated dataset
# Use the RSEC corrected reads per cell counts data table outputted from the BD Targeted analysis pipeline on SevenBridges Genomics
# The sce object is of S4 class and can be converted to a Seurat object and used with a variety of packages


### Either read in data table from SevenBridges and transpose, then construct the SCE object

### OR 

### Read in an already transposed df and process from line 26


library(SingleCellExperiment)
library(data.table)

# read in data table
#df <- read.table("/home/cathal/Main_RSEC_and_SampleTagReads_perCell_merged.txt", header = T, sep = '\t')

# transpose
# tt <- t(df)

## write to file
## write.table(tt, file = "transposed_df4.txt", sep = '\t', col.names = FALSE, row.names = FALSE)

# read in data table
counts_matrix <- read.table("/home/cathal/Downloads/transposed_df3.txt", header = T, sep = '\t')

rownames(counts_matrix) <- counts_matrix$Cell_Index

counts_matrix$Cell_Index <- NULL

counts_matrix <- as.matrix(counts_matrix)

sce <- SingleCellExperiment(assays = list(counts = counts_matrix))


# access the count data assay in the sce object 
assay(sce, "counts")
# assay(test.fromseurat2SCE, "counts")

# see what assays are in the object
assays(sce)
