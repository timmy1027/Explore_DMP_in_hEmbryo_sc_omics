# Seurat Integrate samples using shared highly variable genes
## Reference1: https://satijalab.org/seurat/archive/v3.1/sctransform_vignette.html
## Reference2: https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html
## Reference3: https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/stacked-vlnplot-for-given-features-sets.html
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
################# remote HPC instead ######################

# Integrate across conditions
options(future.globals.maxSize = 100000 * 1024^2) #use 12Gb memory
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# Save integrated seurat object
saveRDS(seurat_integrated, "data/integrated_seurat.RData")

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
PCAPlot(seurat_integrated,
        split.by = "sample")

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP                             
tiff(filename = "figure/s4_Integrated_UMAP.tiff", width = 1500, height = 1300, units = "px", res = 300)
DimPlot(seurat_integrated)
dev.off()

tiff(filename = "figure/s4_Integrated_UMAP_bySample.tiff", width = 3000, height = 3500, units = "px", res = 300)
DimPlot(seurat_integrated,
        split.by = "sample",
        ncol = 3) 
dev.off()








