# GSE157329 human CS12-CS16 embryo data integration into DMP scRNA-seq data

path <- "/Users/tianming1027/Desktop/Seq_Files/SBS_scRNAseq/GSE157329_Xu_Human _CS12_to_CS16"

library(Seurat)
library(Matrix)
library(dplyr)
matrix_dir = "/Users/tianming1027/Desktop/Seq_Files/SBS_scRNAseq/GSE157329_Xu_Human _CS12_to_CS16/So/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
raw_mx <- mat
rm(mat)
nFeatures <- read.delim(features.path,
                           header = FALSE, sep = "")
nBarcodes <- read.delim(barcode.path,
                           header = FALSE)

colnames(raw_mx) = nBarcodes$V1
rownames(raw_mx) = nFeatures$V2

seurat_obj_emb <- CreateSeuratObject(counts = raw_mx, 
                                 min.features = 100,
                                 min.cells = 5,
                                 project = all_sam)
# Export mata.data as csv file
library(dplyr)
library(tibble)
library(Seurat)

export_df <- obj@meta.data %>% 
  rownames_to_column("barcodes") %>%
  select(barcodes, genotype, Clustername)

write.csv(export_df, "/path/to/file.csv")

####################
## Load functions
####################

source("~/Documents/GitHub repository/heoa/src/function.r") # all functions used in iterative clustering


####################
# Preprocessing
####################z 

# Red blood cells have much fewer expressed genes than other cells. To avoid the potential distortion on clustering by the large population of red blood cells, we identified red blood cells by the expression of HBA1, excluded them from clustering.

# Cells with extremely high expression of hemoglobin were considered as erythroids
all_sam <- c('h0','h22','h5','h9a','h9b','ht7','l0','l21','l22','l9','lv5','lv6','t0','t21','t22','t5','t9','tv7','v0','v21','v22','v6','vl5') # all samples
hbg <- c("HBA1","HBA2","HBE1","HBG1","HBG2","HBZ") # Signatures of erythroid
hbg_sum <- apply( raw_mx[hbg,],2,sum)

meta.data <- seurat_obj_emb@meta.data
meta.data$cells <- nBarcodes[1:100,]

# Add number of genes per UMI for each cell to metadata
seurat_obj_emb$log10GenesPerUMI <- log10(seurat_obj_emb$nFeature_RNA) / log10(seurat_obj_emb$nCount_RNA)


# Compute percent mito ratio
seurat_obj_emb$mitoRatio <- PercentageFeatureSet(object = seurat_obj_emb, pattern = "^MT-")
seurat_obj_emb$mitoRatio <- seurat_obj_emb@meta.data$mitoRatio / 100

sum(seurat_obj_emb@meta.data$mitoRatio)

meta.data <- meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
seurat_obj_emb@meta.data <- meta.data

# Filter out low quality reads using selected thresholds - these will change with experiment
options(future.globals.maxSize = 13000 * 1024^2)
sobj.emb.filtered <- subset(x = seurat_obj_emb, 
                                   subset= (nUMI >= 500) & 
                                     (nUMI <= 60000) &
                                     (nGene >= 250) & 
                                     (log10GenesPerUMI > 0.80) & 
                                     (mitoRatio < 0.125))


















