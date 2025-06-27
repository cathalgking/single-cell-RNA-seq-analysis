library(SingleCellExperiment)
library(ggplot2)
library(scater)
library(scran)
library(bluster)
library(BiocNeighbors)

# read in SCE object
sce <- readRDS(file = "other_data/Song_dge_E.rds")

# function to run clustering and dim redux pipeline on SCE object
run_leiden_bluster <- function(sce, n_pcs = 10, resolution = 0.5, objective = "modularity",
                              umap_n_neighbors = 15, umap_min_dist = 0.1) {
  # Compute log-normalized counts if missing
  if (!"logcounts" %in% assayNames(sce)) {
    sce <- logNormCounts(sce)
  }

  # Run PCA if not already done, using "logcounts"
  if (!"PCA" %in% reducedDimNames(sce)) {
    sce <- runPCA(sce, ncomponents = n_pcs, exprs_values = "logcounts")
  }

  # Leiden clustering on first n_pcs PCs
  clusters <- clusterRows(
    reducedDim(sce, "PCA")[, 1:n_pcs],
    BLUSPARAM = SNNGraphParam(
      cluster.fun = "leiden",
      cluster.args = list(
        resolution_parameter = resolution,
        objective_function = objective
      )
    )
  )

  colData(sce)$cluster <- clusters

  # Run UMAP on PCA results and store in reducedDims
  sce <- runUMAP(sce, dimred = "PCA", ncomponents = 2,
                 n_neighbors = umap_n_neighbors, 
                 min_dist = umap_min_dist)

  return(sce)
}

# execute pipeline
sce <- run_leiden_bluster(sce)
