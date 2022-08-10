# Seurat Pre-process
library(Seurat)
library(tidyverse)
library(magrittr)

# Create each individual Seurat object for every sample
for (sample in c("DMP", "BAd1", "BAd3", "BAd6", "SKMd1", "SKMd3", "SKMd6")){
  seurat_data <- Read10X(data.dir = paste0("/Users/tianming1027/Dropbox/Mac/Desktop/Seq_Files/SBS_scRNAseq/cellranger_out/", sample, "/filtered_feature_bc_matrix/"))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100,
                                   min.cells = 3,
                                   project = sample)
  assign(sample, seurat_obj)
}

# Check the metadata in the new Seurat objects
head(DMP@meta.data)
head(BAd3@meta.data)

# Create a merged Seurat object
sample.list <- list(BAd1, BAd3, BAd6, SKMd1, SKMd3, SKMd6)
combined <- merge(x = DMP, 
                  y = sample.list, 
                  add.cell.id = c("DMP", "BAd1", "BAd3", "BAd6", "SKMd1", "SKMd3", "SKMd6"))

# Check that the merged object has the appropriate sample-specific prefixes
head(combined@meta.data)
tail(combined@meta.data)
