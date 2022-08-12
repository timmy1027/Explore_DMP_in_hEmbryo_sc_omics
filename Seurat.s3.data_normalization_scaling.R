# Seurat data Normalization and Scaling

library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

filtered_seurat <- readRDS("data/seurat_filtered.RData")
seurat_phase <- NormalizeData(filtered_seurat)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m.genes, 
                                 s.features = s.genes)
# View cell cycle scores and phases assigned to cells
View(seurat_phase@meta.data) 
table(seurat_phase@meta.data$Phase)


# Assess whether cell cycle is a major source of variation in merged dataset using PCA
# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
tiff(filename = "figure/s3_cellcycle_PCA_splitbysample.tiff", width = 2500, height = 2000, units = "px", res = 300)
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "sample",
        ncol = 3)
dev.off()
## Different samples are heavily biased towards different cc stage.

# Regress out cell cycle scores during data scaling
seurat_phase_regressed <- ScaleData(seurat_phase,
                          vars.to.regress = c("S.Score", "G2M.Score"))
seurat_phase_regressed <- RunPCA(seurat_phase_regressed)

tiff(filename = "figure/s3_cellcycle_PCA_cc_regressed.tiff", width = 2500, height = 2000, units = "px", res = 300)
DimPlot(seurat_phase_regressed,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
dev.off()
tiff(filename = "figure/s3_cellcycle_PCA_splitbysample_cc_regressed.tiff", width = 2500, height = 2000, units = "px", res = 300)
DimPlot(seurat_phase_regressed,
        reduction = "pca",
        group.by= "Phase",
        split.by = "sample",
        ncol = 3)
dev.off()
## cell cycle variants are correctly regressed out. Next, perform more accurate SCTransform. By default, 
## SCTransform accounts for cellular sequencing depth, or nUMIs.

# Before we run this for loop, we know that the output can generate large R objects/variables in terms of memory. 
# If we have a large dataset, then we might need to adjust the limit for allowable object sizes 
# within R (Default is 500 * 1024 ^ 2 = 500 Mb) using the following code:
options(future.globals.maxSize = 10000 * 1024^2)

# Now Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

split_seurat <- split_seurat[c("DMP", "BAd1", "BAd3", "BAd6", "SKMd1", "SKMd3", "SKMd6")]

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m.genes, s.features=s.genes)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("S.Score", "G2M.Score", "mitoRatio"))
}

# Check which assays are stored in objects
split_seurat$DMP@assays

# Check a few marker genes in individual sample
VlnPlot(split_seurat$BAd6, features = c("FOXC1", "PAX3", "UCP1", "NR4A2", "PRDM16"), 
        pt.size = 0.2, ncol = 4)




























