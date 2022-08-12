# Seurat data QC and filtering

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
tiff(filename = "figure/DMP_diff_scRNA_QC.tiff", width = 2000, height = 1200, units = "px", res = 300)
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), split.by = "orig.ident")
dev.off()


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

# Create .RData object to load at any time
saveRDS(combined, file="data/combined_raw_seurat.RData")


# Assessing the quality metrics
## Cell counts
## UMI counts per cell
##Genes detected per cell
## UMIs vs. genes detected
## Mitochondrial counts ratio
## Novelty

# Visualize the number of cell counts per sample
tiff(filename = "figure/QC_nCells.tiff", width = 1400, height = 1200, units = "px", res = 300)
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells")
dev.off()


# Visualize the number UMIs/transcripts per cell
tiff(filename = "figure/QC_nUMIs.tiff", width = 1400, height = 1200, units = "px", res = 300)
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500) +
  geom_vline(xintercept = 60000)
dev.off()
## We can see that majority of our cells in both samples have 1000 UMIs or greater, which is great.

# Visualize the distribution of genes detected per cell via histogram
tiff(filename = "figure/QC_nGenes.tiff", width = 1400, height = 1200, units = "px", res = 300)
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
dev.off()

# Visualize the distribution of genes detected per cell via boxplot
tiff(filename = "figure/QC_nCells_vs_nGenes.tiff", width = 1400, height = 1200, units = "px", res = 300)
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nGenes Detected per Cell")
dev.off()

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
tiff(filename = "figure/QC_nUMIs_vs_nCells.tiff", width = 3600, height = 4000, units = "px", res = 300)
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per cell
tiff(filename = "figure/QC_MitoRatio.tiff", width = 1400, height = 1200, units = "px", res = 300)
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
dev.off()
## This metric can identify whether there is a large amount of mitochondrial contamination from dead or dying cells.
## We define poor quality samples for mitochondrial counts as cells which surpass the 0.2 mitochondrial ratio mark

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
tiff(filename = "figure/QC_Complexity.tiff", width = 1400, height = 1200, units = "px", res = 300)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
dev.off()


# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_combined_seurat <- subset(x = combined, 
                          subset= (nUMI >= 500) & 
                            (nUMI <= 60000) &
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

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
# Save filtered Seurat project
save(filtered_seurat, file="data/seurat_filtered.RData")

saveRDS(filtered_seurat, file="data/seurat_filtered.RData")
test.obj <- readRDS("data/seurat_filtered.RData")


# Extract the new metadata from the filtered Seurat object to go through the same plots as with the unfiltered data
metadata_clean <- filtered_seurat@meta.data

# Visualize the number of cell counts per sample
tiff(filename = "figure/cleaned_QC/QC_nCells.tiff", width = 1400, height = 1200, units = "px", res = 300)
metadata_clean %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells")
dev.off()


# Visualize the number UMIs/transcripts per cell
tiff(filename = "figure/cleaned_QC/QC_nUMIs.tiff", width = 1400, height = 1200, units = "px", res = 300)
metadata_clean %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500) +
  geom_vline(xintercept = 60000)
dev.off()
## We can see that majority of our cells in both samples have 1000 UMIs or greater, which is great.

# Visualize the distribution of genes detected per cell via histogram
tiff(filename = "figure/cleaned_QC/QC_nGenes.tiff", width = 1400, height = 1200, units = "px", res = 300)
metadata_clean %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
dev.off()

# Visualize the distribution of genes detected per cell via boxplot
tiff(filename = "figure/cleaned_QC/QC_nCells_vs_nGenes.tiff", width = 1400, height = 1200, units = "px", res = 300)
metadata_clean %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nGenes Detected per Cell")
dev.off()

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
tiff(filename = "figure/cleaned_QC/QC_nUMIs_vs_nCells.tiff", width = 3600, height = 4000, units = "px", res = 300)
metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per cell
tiff(filename = "figure/cleaned_QC/QC_MitoRatio.tiff", width = 1400, height = 1200, units = "px", res = 300)
metadata_clean %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
dev.off()
## This metric can identify whether there is a large amount of mitochondrial contamination from dead or dying cells.
## We define poor quality samples for mitochondrial counts as cells which surpass the 0.2 mitochondrial ratio mark

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
tiff(filename = "figure/cleaned_QC/QC_Complexity.tiff", width = 1400, height = 1200, units = "px", res = 300)
metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
dev.off()

