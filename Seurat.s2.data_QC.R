# Seurat data QC

## reference: 
### https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
### https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/seurat-pre-process.html#load-count-matrix-from-cellranger
### https://hbctraining.github.io/scRNA-seq/lessons/03_SC_quality_control-setup.html

# Explore merged metadata
View(combined@meta.data)

## orig.ident: project as we had assigned it
## nCount_RNA: number of UMIs per cell
## nFeature_RNA: number of genes detected per cell

## number of genes detected per UMI: this metric with give us an idea of the complexity of our dataset (more genes detected per UMI, more complex our data)
## mitochondrial ratio: this metric will give us a percentage of cell reads originating from the mitochondrial genes

# Add number of genes per UMI for each cell to metadata
combined$log10GenesPerUMI <- log10(combined$nFeature_RNA) / log10(combined$nCount_RNA)


# Compute percent mito ratio
combined$mitoRatio <- PercentageFeatureSet(object = combined, pattern = "^MT-")
combined$mitoRatio <- combined@meta.data$mitoRatio / 100


# Visualize QC metrics as a violin plot
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), split.by = "orig.ident")
# Print out  QC metrics as a violin plot
tiff(filename = "DMP_diff_scRNA_QC.tiff", width = 2000, height = 1200, units = "px", res = 300)
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), split.by = "orig.ident")
dev.off()






















