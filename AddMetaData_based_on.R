## Add a new column to a Seurat objects meta-data which contains entries based on another column in the meta-data.
## Steps:
# Extract the meta-data from the Seurat object to a data-frame
# Add a column to that data-frame
# Add the new meta-data back to the Seurat object

# extract df
new_patient_PCs <- as.data.frame(x = PC_Seu@meta.data$patient)

# Multiple conditions when adding new column to data-frame:
# In the new column to add, if P1, input ND, if P2, call it Normal etc. 
new_patient_PCs_diagnosis <- new_patient_PCs %>% mutate(Diagnosis =
                                                          case_when(p_name == 'P1' ~ "ND", 
                                                                    p_name == 'P3' ~ "ND",
                                                                    p_name == 'P2' ~ "Normal",
                                                                    p_name == 'P5' ~ "Tx",
                                                                    p_name == 'P6' ~ "Tx",
                                                                    p_name == 'P7' ~ "Tx")
)

# Add new metadata to Seurat object.
PC_Seu <- AddMetaData(object = PC_Seu, metadata = new_patient_PCs_diagnosis, col.name = "Diagnosis")