# Seurat Integrate samples using shared highly variable genes
## Reference1: https://satijalab.org/seurat/archive/v3.1/sctransform_vignette.html
#@ Reference2: https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

## Specifically, this integration method expects “correspondences” or shared biological states among at least 
## a subset of single cells across the groups. 

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

################# memory exhausted #########################
################# remote HPC instead #########################

# Integrate across conditions
options(future.globals.maxSize = 12000 * 1024^2) #use 12Gb memory
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# Save integrated seurat object
saveRDS(seurat_integrated, "data/integrated_seurat.rds")

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP                             
DimPlot(seurat_integrated)                             

DimPlot(seurat_integrated,
        split.by = "sample") 








