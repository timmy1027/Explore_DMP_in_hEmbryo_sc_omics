# scRNA-seq analysis, modified from Kagawa's Nautre 2022 paper
## standard normalization method, not SCTransformation

# Seurat Sample Loading & Merging
library(Seurat)
library(tidyverse)
library(magrittr)

# Create each individual Seurat object for every sample
for (sample in c("DMP", "BAd1", "BAd3", "BAd6", "SKMd1", "SKMd3", "SKMd6")){
  seurat_data <- Read10X(data.dir = paste0("/Users/tianming1027/Dropbox/Mac/Desktop/Seq_Files/SBS_scRNAseq/cellranger_out/", sample, "/filtered_feature_bc_matrix/"))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100,
                                   min.cells = 5,
                                   project = sample)
  assign(sample, seurat_obj)
}

# Create a merged Seurat object
sample.list <- list(BAd1, BAd3, BAd6, SKMd1, SKMd3, SKMd6)
combined <- merge(x = DMP, 
                  y = sample.list, 
                  add.cell.id = c("DMP", "BAd1", "BAd3", "BAd6", "SKMd1", "SKMd3", "SKMd6"))
rm(list = c("BAd1", "BAd3", "BAd6", "DMP", "SKMd1", "SKMd3", "SKMd6", "sample.list", "seurat_data", "seurat_obj"))

# Add number of genes per UMI for each cell to metadata
combined$log10GenesPerUMI <- log10(combined$nFeature_RNA) / log10(combined$nCount_RNA)


# Compute percent mito ratio
combined$mitoRatio <- PercentageFeatureSet(object = combined, pattern = "^MT-")
combined$mitoRatio <- combined@meta.data$mitoRatio / 100

# create a metadata data.frame for ease.
metadata <- combined@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^DMP_"))] <- "DMP"
metadata$sample[which(str_detect(metadata$cells, "^SKMd1_"))] <- "SKMd1"
metadata$sample[which(str_detect(metadata$cells, "^SKMd3_"))] <- "SKMd3"
metadata$sample[which(str_detect(metadata$cells, "^SKMd6_"))] <- "SKMd6"
metadata$sample[which(str_detect(metadata$cells, "^BAd1_"))] <- "BAd1"
metadata$sample[which(str_detect(metadata$cells, "^BAd3_"))] <- "BAd3"
metadata$sample[which(str_detect(metadata$cells, "^BAd6_"))] <- "BAd6"


# Add metadata back to Seurat object
combined@meta.data <- metadata

# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_combined_seurat <- subset(x = combined, 
                                   subset= (nUMI >= 500) & 
                                     (nUMI <= 60000) &
                                     (nGene >= 250) & 
                                     (log10GenesPerUMI > 0.80) & 
                                     (mitoRatio < 0.125))

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_combined_seurat, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_combined_seurat@meta.data)

rm(list = c("combined", "counts", "filtered_combined_seurat", "nonzero", "sample", "keep_genes"))

saveRDS(filtered_seurat, file = "data/modified_from_Kagawa/filtered_seurat.RDS")

# Count-data were log-normalized, top 3000 highly variable were selected, and standardization of per gene expression values across cells 
# was performed using NormalizeData, FindVariableFeatures and ScaleData data functions in Seurat. Principal component analysis (PCA) based 
# on the standardized highly variable features was used for linear dimension reduction, a shared nearest neighbor (SNN) graph was constructed 
# on the dimensionally reduced data, and the graph was partitioned using a SNN modularity optimization based clustering algorithm at a range 
# of resolutions using RunPCA, FindNeighbors and FindClusters from Seurat with default settings. Cluster marker genes were identified with the 
# Wilcox likelihood-ratio test using the FindAllMarkers function. Uniform Manifold Approximation and Projection (UMAP) was used for visualization.
set.seed(1234)
varGene = 3000
pc = 20
resolutionwanted=c(0.02,0.1, 0.2, 0.4, 0.6, 0.8, 1)
filtered_seurat <- filtered_seurat %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = varGene) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(verbose = FALSE) %>%
  FindClusters(resolution = resolutionwanted, verbose=FALSE) %>%
  RunUMAP(dims = 1:pc, min.dist = 0.5, verbose=FALSE)

# UMAP projection colored according to sample time
idents <- c("DMP","BAd1","BAd3","BAd6","SKMd1","SKMd3", "SKMd6")
colors <- c("#ff7f00","#4d076a","#ea7bc0","#fa6e6e","#fdbf6f","#33a02c", "#215d6e")
Idents(filtered_seurat) <- "sample"
p <- DimPlot(filtered_seurat, order = idents, cols = colors, pt.size=0.4, split.by = "sample", ncol = 4)
p[[1]]$layers[[1]]$aes_params$alpha = .7
p <- p + theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=3)))

tiff(filename = "figure/figures_Kagawa_pipeline/UMAP_split.tiff", width = 700, height = 500, units = "px")
print(p)
dev.off()

p <- DimPlot(filtered_seurat, order = idents, cols = colors, pt.size=0.4)
p[[1]]$layers[[1]]$aes_params$alpha = .7
p <- p + theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=3)))

tiff(filename = "figure/figures_Kagawa_pipeline/UMAP.tiff", width = 700, height = 500, units = "px")
print(p)
dev.off()

# UMAP projection colored according to cluster identity at various resolutions
ResolutionList <- grep("_snn_res", colnames(filtered_seurat@meta.data), value = TRUE)

for(Resolution in ResolutionList){
  p <- DimPlot(filtered_seurat, label = TRUE, group.by = Resolution, order = T, pt.size = 0.4) +
          theme(aspect.ratio=1) +
          theme(legend.title=element_blank())
  
  ggsave(p, file=paste0("figure/figures_Kagawa_pipeline/UMAP_clusters_", Resolution,".png"), width = 2500, height = 2200, units = "px")
}



counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Visualize selected gene expression
markers.list <- c(
                  "FABP3", "TOP2A", "UBE2C", "BMP7",
                  "TGFBI", "NFIA", "NR2F2", "RUNX1", "RUNX2",
                  "FOXC1", "PAX3", "PDGFRA","BMP5","MEIS1",
                  "NRG3", "ERBB4","SIX1", "EBF1", "HIF3A", "LAMA2", 
                  "COL1A2", "IGFBP2", "IGFBP5", "HOXB3")

tiff(filename = "figure/figures_Kagawa_pipeline/selected_gene_exp_umap_reduced_extended.mincutq10.tiff", width = 1000, height = 1100, units = "px")
FeaturePlot(object = filtered_seurat, reduction = "umap", features = markers.list, order = T, pt.size = 0.4, 
            max.cutoff = "q99", min.cutoff = "q10",cols = c("#FEE0D2","#67000D"))
dev.off()

tiff(filename = "figure/figures_Kagawa_pipeline/selected_gene_exp_umap_reduced_extended.default.tiff", width = 1000, height = 900, units = "px")
FeaturePlot(object = filtered_seurat, reduction = "umap", features = markers.list, order = T, pt.size = 0.4, 
            max.cutoff = "q99",cols = c("#FEE0D2","#67000D"))
dev.off()

# ScreePlot, PC heatmap, Represented PC genes

tiff(filename = "figure/figures_Kagawa_pipeline/Elbow_Plot_PC20.tiff", 
     width = 1500, height = 1200, units = "px", res = 300)
ElbowPlot(object = filtered_seurat, 
          ndims = 20)
dev.off()

tiff(filename = "figure/figures_Kagawa_pipeline/Heatmapt_PC9.tiff", width = 1800, height = 1800, units = "px", res = 300)
DimHeatmap(filtered_seurat, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
dev.off()

print(x = filtered_seurat[["pca"]], 
      dims = 1:10, 
      nfeatures = 20)


VlnPlot(filtered_seurat, features = c("BMP5", "BMP7", "NRG3", "FABP3", "PAX3"), slot = "counts", log = TRUE)


# Find markers and marker heatmap

#wanted.resolution <- "RNA_snn_res.0.2" 
wanted.resolution <- "RNA_snn_res.0.2"
## Identify top markers at wanted resolution

Idents(filtered_seurat) <- wanted.resolution
top.markers <- FindAllMarkers(filtered_seurat, only.pos = TRUE, verbose = FALSE) %>%
  dplyr::filter(p_val_adj < 0.01) %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) 

write.csv(top.markers, file = "FindMarkers.ssn0.2.csv")




library(dplyr)
library(scran)
library(batchelor)
library(Seurat)
library(scales)
library(kableExtra)
















