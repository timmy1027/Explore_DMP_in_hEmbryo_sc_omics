# GSE157329 human CS12-CS16 embryo data integration into DMP scRNA-seq data

path <- "/Users/tianming1027/Desktop/Seq_Files/SBS_scRNAseq/GSE157329_Xu_Human _CS12_to_CS16"

library(dplyr)
library(scran)
library(batchelor)
library(Seurat)
library(scales)
library(kableExtra)
library(Matrix)
library(scater)
library(cowplot)
library(BiocParallel)
library(BiocNeighbors)


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





#https://satijalab.org/seurat/articles/essential_commands.html

##
meta.data <- combined@meta.data
#unique sample names
samples <- as.character(unique(meta.data$seq_folder))
#number of cells per sample
table(meta.data$seq_folder)
#BAd1  BAd3  BAd6   DMP    h0    h5   h9a   h9b   ht7    l0    l9   lv5   lv6 SKMd1 SKMd3 SKMd6    t0    t5    t6    t9   tv7    v0    v6   vl5 
#1841  2347  1841  3735  3527 10081  2462  2389  9300  8414  2759 12859  5239  4091  4017 12068  3380  7224     1  2503  8693  3457  6301 13781 
BAd1  BAd3  BAd6   DMP    h0    h5   h9a   h9b   ht7    l0    l9   lv5   lv6 SKMd1 SKMd3 SKMd6 t0    t5    t6    t9   tv7    v0    v6   vl5

#Data Integration, https://nbisweden.github.io/excelerate-scRNAseq/session-integration/Data_Integration.html
options(future.globals.maxSize = 120000 * 1024^2) #120 Gb
sample.list <- SplitObject(combined, split.by = "seq_folder")

#filtered h0
counts <- GetAssayData(object = sample.list$t6, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered.t6 <- CreateSeuratObject(filtered_counts,
                                  meta.data = sample.list$t6@meta.data)

keep_genes <- Reduce(intersect, list(rownames(filtered.BAd1),rownames(filtered.BAd3),
                                     rownames(filtered.BAd6),rownames(filtered.DMP),
                                     rownames(filtered.SKMd1),rownames(filtered.SKMd3),
                                     rownames(filtered.SKMd6),rownames(filtered.t0),
                                     rownames(filtered.t5),
                                     rownames(filtered.t9),rownames(filtered.l9),
                                     rownames(filtered.tv7),rownames(filtered.h0),
                                     rownames(filtered.h5),rownames(filtered.h9a),
                                     rownames(filtered.ht7),rownames(filtered.v0),
                                     rownames(filtered.v6),rownames(filtered.lv5),
                                     rownames(filtered.lv6),rownames(filtered.vl5),
                                     rownames(filtered.l0)
                                     ))


data.l0 <- as.SingleCellExperiment(sample.list$l0)
data.l0 <- data.l0[match(keep_genes, rownames(data.l0)), ]
data.l0 <- computeSumFactors(data.l0)
data.l0 <- logNormCounts(data.l0)

# total raw counts
counts_samples <- cbind(counts(data.BAd1),counts(data.BAd3),
                        counts(data.BAd6),counts(data.DMP),
                        counts(data.SKMd1),counts(data.SKMd3),
                        counts(data.SKMd6),counts(data.t0),
                        counts(data.t5),
                        counts(data.t9),counts(data.l9),
                        counts(data.tv7),counts(data.h0),
                        counts(data.h5),counts(data.h9a),
                        counts(data.ht7),counts(data.v0),
                        counts(data.v6),counts(data.lv5),
                        counts(data.lv6),counts(data.vl5),
                        counts(data.l0)
                        )

mBN.counts.samples <- multiBatchNorm(data.BAd1, data.BAd3, data.BAd6,data.DMP,
                                     data.SKMd1, data.SKMd3,data.SKMd6,
                                     data.t0, data.t5, data.t9, data.l9,
                                     )

# total normalized counts (with multibatch normalization)
logcounts_samples <- cbind(logcounts(data.BAd1),logcounts(data.BAd3),
                           logcounts(data.BAd6),logcounts(data.DMP),
                           logcounts(data.SKMd1),logcounts(data.SKMd3),
                           logcounts(data.SKMd6),logcounts(data.t0),
                           logcounts(data.t5),
                           logcounts(data.t9),logcounts(data.l9),
                           logcounts(data.tv7),logcounts(data.h0),
                           logcounts(data.h5),logcounts(data.h9a),
                           logcounts(data.ht7),logcounts(data.v0),
                           logcounts(data.v6),logcounts(data.lv5),
                           logcounts(data.lv6),logcounts(data.vl5),
                           logcounts(data.l0)
)

# sce object of the combined data 
sce <- SingleCellExperiment( 
  assays = list(counts = counts_samples, logcounts = logcounts_samples),  
  rowData = rowData(data.l0), 
  colData = rbind(colData(data.BAd1),colData(data.BAd3),
                  colData(data.BAd6),colData(data.DMP),
                  colData(data.SKMd1),colData(data.SKMd3),
                  colData(data.SKMd6),colData(data.t0),
                  colData(data.t5),
                  colData(data.t9),colData(data.l9),
                  colData(data.tv7),colData(data.h0),
                  colData(data.h5),colData(data.h9a),
                  colData(data.ht7),colData(data.v0),
                  colData(data.v6),colData(data.lv5),
                  colData(data.lv6),colData(data.vl5),
                  colData(data.l0)) 
)

saveRDS(object = sce, file = "sce.RDS")

#################################################
combined <- readRDS("DMP.hoea.combined.filtered.RDS")
sce.ob <- readRDS("sce.ob.RDS")
meta.data <- combined@meta.data
combined.counts <- GetAssayData(object = combined, slot = "counts")
meta.data$cell <- rownames(meta.data)
experiments <- c("hoea", "DMP")
sce.ob <- list()
for (b in experiments){ 
  temp.M <- meta.data %>% filter(project == b) 
  temp.sce <-  SingleCellExperiment(list(counts = as.matrix(combined.counts[keep_genes,temp.M$cell])), colData = temp.M) %>%
    computeSumFactors()
  sce.ob[[b]] <- temp.sce
}


mBN.sce.ob <- multiBatchNorm(sce.ob$hoea, sce.ob$DMP)
names(mBN.sce.ob) <- experiments

lognormExp.mBN <- mBN.sce.ob %>%
  lapply(function(x){logcounts(x) %>% as.data.frame() 
    %>% return()}) %>%
  do.call("bind_cols",.)

varGene = 2000
pc = 20

sobj <-
  CreateSeuratObject(
    combined.counts[rownames(lognormExp.mBN),meta.data$cell], meta.data = meta.data
  ) %>%
  NormalizeData(verbose = FALSE)

sobj@assays$RNA@data <- as.matrix(lognormExp.mBN[rownames(lognormExp.mBN),colnames(sobj)])

sobj.list <- SplitObject(sobj, split.by = "project") %>%
  lapply(function(x){ x = FindVariableFeatures(x, verbose=F, nfeatures = varGene)})
sobj.list <- sobj.list[experiments]

rm(combiend.counts, combined, lognormExp.mBN, mBN.sce.ob, sce.ob)

set.seed(123)
sobj <- SeuratWrappers::RunFastMNN(sobj.list, verbose = F, features = varGene) %>%
  RunUMAP(reduction = "mnn", dims = 1:pc, verbose = F) %>%
  FindNeighbors(reduction = "mnn", dims = 1:pc, verbose = F)

for(plotby in c("hoea.cell.type","DMP_sample","hoea.embryo.stage")){
  Idents(sobj) <- plotby
  
  if(plotby == "hoea.embryo.stage"){
    orderidents <- unique(sobj@meta.data$hoea.embryo.stage)
    cols<-rev(c("#f0590e","#cca95e","#960000","gray90"))
  } else if(plotby == "DMP_sample"){
    orderidents <- c("DMP", "BAd1", "BAd3", "BAd6", "SKMd1", "SKMd3", "SKMd6", "other")
    cols<-rev(c("#ff6361","#58508d","#5ecc74","#bc5090","#cca95e","#5e69cc","#b81212","gray90"))
  }else{
    orderidents <- c(setdiff(levels(Idents(sobj)),"other"),"other")
    cols <- c("gray90", hue_pal()(length(orderidents)-1))
  }
  
  sobj[[plotby]] <- factor(Idents(sobj), levels = orderidents)
  ncol <- max(3, floor(length(orderidents)/6))
  
  print(
    DimPlot(sobj, order=orderidents, cols=cols, pt.size=0.3, raster = FALSE) +
      theme(aspect.ratio=1) +
      theme(legend.position="bottom", legend.title=element_blank()) +
      guides(color = guide_legend(override.aes=list(size=3), ncol=ncol))
  )
  
}
















