## subset single-cell object to only contain certain genes

# define genes
genes <- c("CCL14", "MRC1", "TGM2", "CXADR", "VWF", "SEMA3A", "NAMPT", "A2M", "CD59", "HSPG2")

# first check if genes are in the rownames
genes %in% rownames(sce)

## gene names for rows and column names for cell/barcode

# Seurat object
seu_filt <- seu[genes,]

# SCE object
sce_filt <- sce[genes,]
