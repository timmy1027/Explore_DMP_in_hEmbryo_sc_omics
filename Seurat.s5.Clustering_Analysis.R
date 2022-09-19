# Seurat Clustering Analysis
## Reference1: https://satijalab.org/seurat/archive/v3.1/sctransform_vignette.html
## Reference2: https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html
## Reference3: https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/stacked-vlnplot-for-given-features-sets.html
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

seurat_integrated <- readRDS("data/s4_integrated_seurat.RData")
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
tiff(filename = "figure/s5_Cell_Cluster_UMAP_res_0.4_bySample.tiff", 
         width = 3500, height = 1200, units = "px", res = 300)
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 3,
        split.by = "sample")
dev.off()

# Goals:To generate cell type-specific clusters and use known markers to determine the identities of the clusters.

# Recommendations:
## Have a good idea of your expectations for the cell types to be present and a handful of marker genes for these cell types. 
## Know whether you expect cell types of low complexity or higher mitochondrial content AND whether the cells are differentiating
## Identify any junk clusters for removal. Possible junk clusters could include those with high mitochondrial content and 
## low UMIs/genes. If not detecting all cell types as separate clusters, try changing the UMAP resolution, the number of 
## PCs used for clustering, or the number of variable genes used

# Segregation of clusters by sample
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)

# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()

# Explore whether clusters segregate by cell cycle phase
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

tiff(filename = "figure/s5_Integrated_Seurat_Metrics.tiff", 
     width = 2000, height = 2500, units = "px", res = 300)
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
dev.off()


seurat.DMP <- seurat_integrated[, seurat_integrated$sample == "DMP"]
seurat.SKMd1 <- seurat_integrated[, seurat_integrated$sample == "SKMd1"]
seurat.SKMd3 <- seurat_integrated[, seurat_integrated$sample == "SKMd3"]
seurat.SKMd6 <- seurat_integrated[, seurat_integrated$sample == "SKMd6"]
seurat.BAd1 <- seurat_integrated[, seurat_integrated$sample == "BAd1"]
seurat.BAd3 <- seurat_integrated[, seurat_integrated$sample == "BAd3"]
seurat.BAd6 <- seurat_integrated[, seurat_integrated$sample == "BAd6"]

FeaturePlot(seurat.DMP, 
            reduction = "umap", 
            features = c("PAX3", "FOXC1", "MYF5", "PRDM16", "NR2A4"),
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
tiff(filename = "figure/s5_PC1_genes_SKMd3.tiff", 
     width = 2000, height = 2500, units = "px", res = 300)
FeaturePlot(seurat.SKMd1, 
            reduction = "umap", 
            features = c("COL1A1", "FN1", "TGFBI", "NRG3", "PAX3", "BMP5"),
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            )
dev.off()



# Exploration of the PCs driving the different clusters
# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)
# Check the top 16 PCs
# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
tiff(filename = "figure/s5_PC_weighting.tiff", 
     width = 2500, height = 2500, units = "px", res = 300)
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = "none", 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)
dev.off()

# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:16, nfeatures = 20)

















