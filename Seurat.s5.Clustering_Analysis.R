# Seurat Clustering Analysis
## Reference1: https://satijalab.org/seurat/archive/v3.1/sctransform_vignette.html
## Reference2: https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html
## Reference3: https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/stacked-vlnplot-for-given-features-sets.html
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)


## Seurat assigns cells to clusters based on their PCA scores 
## derived from the expression of the integrated most variable genes.

## Determining how many PCs to include in the clustering step is 
## therefore important to ensure that we are capturing the 
## majority of the variation, or cell types, present in our dataset.

# Explore heatmap of PCs
tiff(filename = "figure/s5_Determine_PCs.tiff", 
     width = 3500, height = 3000, units = "px", res = 300)
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
dev.off()


# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
# Plot the elbow plot
tiff(filename = "figure/s5_Elbow_Plot_PC40.tiff", 
     width = 1500, height = 1200, units = "px", res = 300)
ElbowPlot(object = seurat_integrated, 
          ndims = 40)
dev.off()


# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
# Explore resolutions
seurat_integrated@meta.data %>% 
  View()


# Assign identity of clusters
## To choose a resolution to start with, 
## pick something in the middle of the range like 0.6 or 0.8. 
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

# Plot the UMAP
tiff(filename = "figure/s5_Cell_Cluster_UMAP_res_0.4.tiff", 
         width = 1500, height = 1200, units = "px", res = 300)
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 3)
dev.off()





























