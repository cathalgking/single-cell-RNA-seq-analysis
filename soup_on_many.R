###
library(SoupX)
library(Seurat)

# Setup paths
starting_dir <- '/homes/cathal.king/Myeloma/fastqs'
first_level_subdirs <- list.dirs(starting_dir, full.names = TRUE, recursive = FALSE)
subdirs_two_levels_down <- c()
# Loop through the first-level subdirectories
for (first_level_subdir in first_level_subdirs) {
  # Get subdirectories directly under the first-level subdirectory
  second_level_subdirs <- list.dirs(first_level_subdir, full.names = TRUE, recursive = FALSE)
  # Add the second-level subdirectories to the list
  subdirs_two_levels_down <- c(subdirs_two_levels_down, second_level_subdirs)
}
subdirs_two_levels_down <- subdirs_two_levels_down[c(4:9,13:15,20:22,29:31,35:37)]

output.dir <- '/homes/cathal.king/Myeloma/soup_corrected_matrices_6p/'

# start the loop
for (file.dir in subdirs_two_levels_down){
  # Extract sample name from path
  sample.name <- basename(file.dir)
  
  # Set up file dir for the 10x base folder
  file.dir_10x <- paste0(file.dir,"/outs")
  
  # Run SoupX
  sc = load10X(file.dir_10x)
  sc = setContaminationFraction(sc, 0.06)
  out = adjustCounts(sc, roundToInt=TRUE)
  
  # Create and save Seurat object
  srat = CreateSeuratObject(out)
  saveRDS(srat, paste0(output.dir,sample.name,'.rds'))
}