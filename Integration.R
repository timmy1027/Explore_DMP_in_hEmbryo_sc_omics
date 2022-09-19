# GSE157329 human CS12-CS16 embryo data integration into DMP scRNA-seq data

path <- "/Users/tianming1027/Desktop/Seq_Files/SBS_scRNAseq/GSE157329_Xu_Human _CS12_to_CS16"

library(dplyr)
library(scran)
library(batchelor)
library(Seurat)
library(scales)
library(kableExtra)
library(Matrix)

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
table(nBarcodes$V2)
table(nBarcodes$V4)


colnames(raw_mx) = nBarcodes$V1
rownames(raw_mx) = nFeatures$V2

seurat_obj_emb <- CreateSeuratObject(counts = raw_mx, 
                                 min.features = 100,
                                 min.cells = 5
                                 )


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
sum(rownames(meta.data) == nBarcodes$V1) == nrow(nBarcodes)
meta.data$hoea.tissue.cluster <- nBarcodes$V2
meta.data$hoea.cell.type <- nBarcodes$V3
meta.data$hoea.cell.subtype <- nBarcodes$V4
meta.data$hoea.embryo.no <- nBarcodes$V5
meta.data$hoea.embryo.stage <- nBarcodes$V6
meta.data$hoea.tissue.region <- nBarcodes$V7


# Add number of genes per UMI for each cell to metadata
meta.data$log10GenesPerUMI <- log10(seurat_obj_emb$nFeature_RNA) / log10(seurat_obj_emb$nCount_RNA)


# Compute percent mito ratio
meta.data$mitoRatio <- PercentageFeatureSet(object = seurat_obj_emb, pattern = "^MT-")
meta.data$mitoRatio <- meta.data$mitoRatio / 100


meta.data <- meta.data %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Export mata.data as csv file
library(dplyr)
library(tibble)

# not used
export_df <- obj@meta.data %>% 
  rownames_to_column("barcodes") %>%
  select(barcodes, genotype, Clustername)
#not used

#use this. path="/hoea/data/So"
write.csv(meta.data, "TMW.hoea.all.metadata.csv")

# Add metadata back to Seurat object
seurat_obj_emb@meta.data <- meta.data



# Filter out low quality reads using selected thresholds - these will change with experiment
options(future.globals.maxSize = 120000 * 1024^2) #120 Gb
sobj.emb.filtered <- subset(x = seurat_obj_emb, 
                                   subset= (nUMI >= 500) & 
                                     (nUMI <= 60000) &
                                     (nGene >= 250) & 
                                     (log10GenesPerUMI > 0.80) & 
                                     (mitoRatio < 0.125))
sobj.emb.filtered@meta.data$project <- rep("hoea", nrow(sobj.emb.filtered@meta.data))
sobj.emb.filtered@meta.data <- sobj.emb.filtered@meta.data[, -c(12,13)]
Idents(sobj.emb.filtered) <- sobj.emb.filtered@meta.data$hoea.cell.type

#merge DMP object with hoea seurat object
sobj.DMP.filtered <- readRDS("filtered_seurat.RDS")
sobj.DMP.filtered@meta.data$project <- rep("DMP", nrow(sobj.DMP.filtered@meta.data))
##re-level idents based on snn resolution = 0.1
Idents(sobj.DMP.filtered) <- "RNA_snn_res.0.1"
levels(Idents(sobj.DMP.filtered))
DMP.cluster.ids <- c("SKMd6_uniq", "DMP", "BA_PAX3+", 
                     "SKM_BALike", "SKM_shared", "SKMd3_uniq",
                     "BA_SIX1+")
names(DMP.cluster.ids) <- levels(sobj.DMP.filtered)
sobj.DMP.filtered <- RenameIdents(sobj.DMP.filtered, DMP.cluster.ids)

DMP.metadata <- sobj.DMP.filtered@meta.data
colnames(DMP.metadata)
colnames(meta.data)
DMP.metadata <- DMP.metadata[, -c(1:3)]
DMP.metadata$seq_folder  <- as.factor(DMP.metadata$seq_folder)
DMP.metadata <- DMP.metadata[, -c(6,8:15)]
sobj.DMP.filtered@meta.data <- DMP.metadata
unique(Idents(sobj.DMP.filtered))
colnames(sobj.DMP.filtered@meta.data) <- c("seq_folder",
                                           "nUMI", "nGene",
                                           "log10GenesPerUMI",
                                           "mitoRatio",
                                           "DMP_sample",
                                           "project")

#merge the two objects
combined <- merge(sobj.emb.filtered, y = sobj.DMP.filtered, 
                  add.cell.ids = c("hoea", "DMP"))

colnames(combined@meta.data)
dim(combined)

combined[["old.ident"]] <- Idents(object = combined)

GetAssayData(combined)["PITX1", 130000:130025]
library(tidyr)
meta1 <- combined@meta.data
meta1$DMP_sample <- replace_na(meta1$DMP_sample, "other")
unique(meta1$DMP_sample)
meta1$hoea.cell.type <- replace_na(meta1$hoea.cell.type, "other")
meta1$hoea.cell.subtype <- replace_na(meta1$hoea.cell.subtype, "other")
meta1$hoea.embryo.no <- replace_na(meta1$hoea.embryo.no, "other")
meta1$hoea.embryo.stage <- replace_na(meta1$hoea.embryo.stage, "other")
meta1$hoea.tissue.region <- replace_na(meta1$hoea.tissue.region, "other")
meta1$hoea.tissue.cluster <- replace_na(meta1$hoea.tissue.cluster, "other")

combined@meta.data <- meta1 

write.csv(meta1, "combined.raw.metadata.csv")
saveRDS(combined, file = "combined.raw.RDS")


#Data Integration
combined.normalized <- NormalizeData(combined, verbose = FALSE)
combined.normalized@assays$RNA@data[1:10, 1:20]

experiments <- c("DMP", "hoea")
varGene = 2000
pc = 20
sobj.list <- SplitObject(combined.normalized, split.by = "project") %>%
  lapply(function(x){ x = FindVariableFeatures(x, verbose=F, nfeatures = varGene)})
sobj.list <- sobj.list[experiments]

set.seed(123)
combined.normalized <- SeuratWrappers::RunFastMNN(sobj.list, verbose = F, features = varGene) %>%
  RunUMAP(reduction = "mnn", dims = 1:pc, verbose = F) %>%
  FindNeighbors(reduction = "mnn", dims = 1:pc, verbose = F)

combined.normalized[["old.ident"]] <- Idents(object = combined.normalized)

idents <- as.character(unique(Idents(combined.normalized)))
DimPlot(combined.normalized, order = idents, 
        reduction = "mnn",raster = FALSE)

tiff(filename = test, width = 800, height = 1000)
FeaturePlot(combined.normalized, features = "PITX1", split.by = "sample", max.cutoff = 3, 
            cols = c("grey", "red"), raster=FALSE)
dev.off()


#https://satijalab.org/seurat/articles/essential_commands.html
