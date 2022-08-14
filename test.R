library(Seurat)
library(data.table)
setwd("~")
setwd("/mnt/bioinfo/longqy/coad/single_cell/result/")
rm(list = ls())


da1 = fread("GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt",header = T,sep = "\t")
da2 = fread("GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt",header = T,sep = "\t")
all = cbind(da1,da2)
all = data.frame(all)
rownames(all) = all[,1]
all2 = all[,c(2:63690,63692:ncol(all))]

sample_name = as.character(unique(sapply(colnames(all2),function(x){strsplit(x,split = "_")[[1]][1]})))

sample_each = list()
for (i in 1:length(sample_list)) {
  idx = sample_list[i]
  sample_da = all2[,which(as.character(sapply(colnames(all2),function(x){strsplit(x,split = "_")[[1]][1]})) == idx)]
  sample_seurat = CreateSeuratObject(counts = sample_da,project = "coad",assay = "RNA")
  sample_seurat$stim = idx
  sample_seurat = NormalizeData(sample_seurat,normalization.method = "LogNormalize",scale.factor = 10000)
  sample_seurat = FindVariableFeatures(sample_seurat,selection.method = "vst",nfeatures=2000)
  sample_each[[i]] = sample_seurat
}

save(sample_each,file="sample_each.rData")

load("sample_each.rData",verbose = T)
sample_list = sample_each


idx = sapply(sample_list, function(x){
  num = nrow(x@meta.data)
})

remain_idx = which(idx > 1000)

ss = as.character(sapply(sample_name,function(x){strsplit(x,split = "[.]")[[1]][1]}))
select = names(which(table(ss) == 2))
select2 = sample_name[which(ss %in% select)]

names(sample_list) = sample_name

used_list = sample_list[select2]

save(used_list,file = "used_sample_seurat.rData")

########## 3. 样本整合

anchors = FindIntegrationAnchors(object.list = used_list,dims=1:20)
integrated = IntegrateData(anchorset = anchors,dims = 1:20)


#integrated = FindVariableFeatures(integrated)
integrated = ScaleData(integrated,features = VariableFeatures(integrated))

integrated <- RunPCA(integrated, verbose = TRUE,features = VariableFeatures(object = integrated))
ElbowPlot(integrated,ndims = 50)

integrated<- FindNeighbors(integrated, reduction = "pca", dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.2)


#integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:20)
integrated <- RunTSNE(integrated, dims = 1:20,do.fast = TRUE)

object = integrated

#DimPlot(object, reduction = "umap", label = TRUE)

TSNEPlot(anno,label=TRUE)
FeaturePlot(object = anno,features = c("EPCAM","DCN","CD68","CD3D","CD79A"),reduction = "tsne",max.cutoff = "q80",min.cutoff = "q10",label=T)

save(object,file="SMC_cluster.rData")


anno = object
tt = as.character(anno@active.ident)
tt[tt=="0"] = "T cells"
tt[tt=="2"] = "T cells"
tt[tt=="7"] = "T cells"
tt[tt=="11"] = "T cells"
tt[tt=="1"] = "Epithelial"
tt[tt=="9"] = "Epithelial"
tt[tt=="6"] = "Stromal"
tt[tt=="8"] = "Stromal"
tt[tt=="10"] = "Stromal"
tt[tt=="3"] = "Myeloid"
tt[tt=="4"] = "B cells"
tt[tt=="5"] = "B cells"
tt[tt=="3"] = "B cells"
Idents(anno) = tt
TSNEPlot(anno,label=T)

save(anno,file = "SMC_anno.rData")

############ 4. T cell 分亚群
setwd("~")
setwd("/mnt/bioinfo/longqy/coad/single_cell/result/")

load("SMC_anno.rData",verbose = T)

tcell = subset(anno,idents = "T cells")
object = tcell

object <- RunPCA(object, verbose = TRUE,features = VariableFeatures(object = object))
ElbowPlot(object,ndims = 50)

object <- FindNeighbors(object, reduction = "pca", dims = 1:15)
object <- FindClusters(object, resolution = 0.8)
object <- RunTSNE(object, dims = 1:15,do.fast = TRUE)

TSNEPlot(object = object,label=T)

marker = read.csv("tcell_marker.csv",header = T)
marker2 = split.data.frame(marker,f = marker$type)

FeaturePlot(object,features = tmp,label = T,min.cutoff = "q10",max.cutoff = "q95",reduction = "tsne")



#############注释
anno_t = object
tt = as.character(anno_t@active.ident)
tt[tt=="4"] = "CD4+"
tt[tt=="7"] = "CD4+"
tt[tt=="10"] = "CD4+"
tt[tt=="1"] = "CD8+"
tt[tt=="3"] = "CD8+"
tt[tt=="6"] = "CD8+"
tt[tt=="9"] = "CD8+"
tt[tt=="2"] = "Regulatory T"
tt[tt=="8"] = "Regulatory T"
tt[tt=="0"] = "Th17"
tt[tt=="13"] = "Tfh"
tt[tt=="5"] = "NK"
tt[tt=="12"] = "NK"
tt[tt=="11"] = "Proliferation"

Idents(anno_t) = tt
TSNEPlot(anno_t,label=T)

save(anno_t,file = "tcell_anno.rData")

############ 5.ATPase 表达分布情况
load("tcell/tcell_anno.rData",verbose = T)
atp = read.csv("atp.csv",header = F)
rownames(atp) = atp[,1]

mettl = read.csv("mettl.csv",header = F)
rownames(mettl) = mettl[,1]

ss = anno_t$stim
n = c(paste("SMC0",1:9,".N",sep = ""),"SMC10.N")
t = c(paste("SMC0",1:9,".T",sep = ""),"SMC10.T")
for (i in 1:length(n)) {
  ni = n[i]
  ss[ss==ni] = "Normal"
}

for (i in 1:length(t)) {
  ti = t[i]
  ss[ss==ti] = "Tumor"
}

anno_t$stim2 = ss

VlnPlot(object = anno_t,features = rownames(atp),split.by  = "stim2",pt.size = 0)
VlnPlot(object = anno_t,features = rownames(mettl),split.by  = "stim2",pt.size = 0)

########### 6. stromal分亚群
setwd("~")
setwd("/mnt/bioinfo/longqy/coad/single_cell/result/")

load("SMC_anno.rData",verbose = T)
stromal = subset(anno,idents="Stromal")

object = stromal

object <- RunPCA(object, verbose = TRUE,features = VariableFeatures(object = object))
ElbowPlot(object,ndims = 50)

object <- FindNeighbors(object, reduction = "pca", dims = 1:10)
object <- FindClusters(object, resolution = 0.8)
object <- RunTSNE(object, dims = 1:10,do.fast = TRUE)

TSNEPlot(object = object,label=T)

marker = read.csv("stromal_marker.csv",header = T)
marker2 = as.character(marker[,2])
FeaturePlot(object = object,features = marker2,reduction = "tsne",label = T,min.cutoff = "q20",max.cutoff = "q80")
VlnPlot(object = object,features = marker2,pt.size = 0,sort = "increasing")

anno_stromal = object

tt = as.character(anno_stromal@active.ident)
tt[tt=="6"] = "Myofibroblasts"
tt[tt=="11"] = "Myofibroblasts"
tt[tt=="15"] = "Proliferative endothelial cells"
tt[tt=="0"] = "Stromal 1"
tt[tt=="7"] = "Stromal 1"
tt[tt=="9"] = "Stromal 1"
tt[tt=="10"] = "Stromal 1"
tt[tt=="1"] = "Stromal 2"
tt[tt=="3"] = "Stromal 3"
tt[tt=="4"] = "Stromal 3"
tt[tt=="5"] = "Tip-like endothelial cells"
tt[tt=="14"] = "Tip-like endothelial cells"
tt[tt=="13"] = "Tip-like endothelial cells"
tt[tt=="12"] = "Tip-like endothelial cells"
tt[tt=="2"] = "Stalk-like endothelial cells"
tt[tt=="8"] = "Lymphatic endothelial cells"


Idents(anno_stromal) = tt
TSNEPlot(anno_stromal,label=T)
save(anno_stromal,file = "stromal_sub_anno.rData")

load("stromal/stromal_sub_anno.rData",verbose = T)


atp = read.csv("atp.csv",header = F)
rownames(atp) = atp[,1]



ss = anno_stromal$stim
n = c(paste("SMC0",1:9,".N",sep = ""),"SMC10.N")
t = c(paste("SMC0",1:9,".T",sep = ""),"SMC10.T")
for (i in 1:length(n)) {
  ni = n[i]
  ss[ss==ni] = "Normal"
}

for (i in 1:length(t)) {
  ti = t[i]
  ss[ss==ti] = "Tumor"
}

anno_stromal$stim2 = ss

VlnPlot(object = anno_stromal,features = rownames(atp),split.by  = "stim2",pt.size = 0)
VlnPlot(object = anno_stromal,features = rownames(mettl),split.by  = "stim2",pt.size = 0)

########### 7. myeloid 分亚群
load("SMC_anno.rData",verbose = T)
TSNEPlot(anno,label=T)

myeloid = subset(anno,idents="Myeloid")

object = myeloid

object <- RunPCA(object, verbose = TRUE,features = VariableFeatures(object = object))
ElbowPlot(object,ndims = 50)

object <- FindNeighbors(object, reduction = "pca", dims = 1:10)
object <- FindClusters(object, resolution = 0.5)
object <- RunTSNE(object, dims = 1:10,do.fast = TRUE)

TSNEPlot(object = object,label=T)

marker = read.csv("myeloid/myeloid_marker.csv",header = T)
marker2 = as.character(marker[,2])
FeaturePlot(object = object,features = marker2,reduction = "tsne",label = T,min.cutoff = "q20",max.cutoff = "q80")
VlnPlot(object = object,features = marker2,pt.size = 0,sort = "increasing")


anno_myeloid = object

tt = as.character(anno_myeloid@active.ident)
tt[tt=="2"] = "Proinflammatory macrophages"
tt[tt=="3"] = "Proinflammatory macrophages"
tt[tt=="0"] = "Anti-Proinflammatory macrophages"
tt[tt=="6"] = "Anti-Proinflammatory macrophages"
tt[tt=="4"] = "SPP1+ macrophages"
tt[tt=="1"] = "Conventional dendritic cells"
tt[tt=="7"] = "Proliferating macrophages"
tt[tt=="5"] = "Unknown"


Idents(anno_myeloid) = tt
TSNEPlot(anno_myeloid,label=T)
save(anno_myeloid,file = "myeloid/myeloid_sub_anno.rData")

load("myeloid/myeloid_sub_anno.rData",verbose = T)



atp = read.csv("atp.csv",header = F)
rownames(atp) = atp[,1]

mettl = read.csv("mettl.csv",header = F)
rownames(mettl) = mettl[,1]

ss = anno_myeloid$stim
n = c(paste("SMC0",1:9,".N",sep = ""),"SMC10.N")
t = c(paste("SMC0",1:9,".T",sep = ""),"SMC10.T")
for (i in 1:length(n)) {
  ni = n[i]
  ss[ss==ni] = "Normal"
}

for (i in 1:length(t)) {
  ti = t[i]
  ss[ss==ti] = "Tumor"
}

anno_myeloid$stim2 = ss

VlnPlot(object = anno_myeloid,features = rownames(atp),split.by  = "stim2",pt.size = 0)
VlnPlot(object = anno_myeloid,features = rownames(mettl),split.by  = "stim2",pt.size = 0)

anno_myeloid$stim2

dd = as.matrix(anno_myeloid@assays$RNA@data)
head(dd[1:10,1:10])

yyn = names(which(anno_myeloid$stim2 == "Normal" & as.character(anno_myeloid@active.ident) == "SPP1"))
yyt = names(which(anno_myeloid$stim2 == "Tumor"))


t1 = dd["ATP6V0B",yyn]
t2 = dd["ATP6V0B",yyt]

########## 8. B cell分亚群

load("SMC_anno.rData",verbose = T)
TSNEPlot(anno,label=T)

bcell = subset(anno,idents="B cells")

object = bcell

object <- RunPCA(object, verbose = TRUE,features = VariableFeatures(object = object))
ElbowPlot(object,ndims = 50)

object <- FindNeighbors(object, reduction = "pca", dims = 1:10)
object <- FindClusters(object, resolution = 0.4)
object <- RunTSNE(object, dims = 1:10,do.fast = TRUE)

TSNEPlot(object = object,label=T)

marker = read.csv("bcell/bcell_marker.csv",header = T)
marker2 = as.character(marker[,2])
marker2 = as.character(na.omit(marker2))
FeaturePlot(object = object,features = marker2,reduction = "tsne",label = T,min.cutoff = "q20",max.cutoff = "q80")
VlnPlot(object = object,features = marker2,pt.size = 0,sort = "increasing")


anno_bcell = object

tt = as.character(anno_bcell@active.ident)
tt[tt=="4"] = "IgG+ plasma"
tt[tt=="7"] = "IgG+ plasma"
tt[tt=="0"] = "IgA+ plasma"
tt[tt=="2"] = "IgA+ plasma"
tt[tt=="6"] = "IgA+ plasma"
tt[tt=="8"] = "Unspecified plasma"
tt[tt=="1"] = "CD19+ CD20+ B"
tt[tt=="3"] = "CD19+ CD20+ B"
tt[tt=="5"] = "CD19+ CD20+ B"
tt[tt=="9"] = "Unknown"


Idents(anno_bcell) = tt
TSNEPlot(anno_bcell,label=T)
save(anno_bcell,file = "bcell//bcell_anno.rData")


load("bcell/bcell_anno.rData",verbose = T)


atp = read.csv("atp.csv",header = F)
rownames(atp) = atp[,1]
ss = anno_bcell$stim
n = c(paste("SMC0",1:9,".N",sep = ""),"SMC10.N")
t = c(paste("SMC0",1:9,".T",sep = ""),"SMC10.T")
for (i in 1:length(n)) {
  ni = n[i]
  ss[ss==ni] = "Normal"
}

for (i in 1:length(t)) {
  ti = t[i]
  ss[ss==ti] = "Tumor"
}

anno_bcell$stim2 = ss

VlnPlot(object = anno_bcell,features = rownames(atp),split.by  = "stim2",pt.size = 0)
VlnPlot(object = anno_bcell,features = rownames(mettl),split.by  = "stim2",pt.size = 0)
########### 9. epithelial 分亚群
load("SMC_anno.rData",verbose = T)
TSNEPlot(anno,label=T)

epith = subset(anno,idents="Epithelial")

object = epith

object <- RunPCA(object, verbose = TRUE,features = VariableFeatures(object = object))
ElbowPlot(object,ndims = 50)

object <- FindNeighbors(object, reduction = "pca", dims = 1:15)
object <- FindClusters(object, resolution = 0.5)
object <- RunTSNE(object, dims = 1:15,do.fast = TRUE)

TSNEPlot(object = object,label=T)

marker = read.csv("epith/epith_marker.csv",header = T)
marker2 = as.character(marker[,2])
marker2 = as.character(na.omit(marker2))
FeaturePlot(object = object,features = marker2,reduction = "tsne",label = T,min.cutoff = "q20",max.cutoff = "q80")
VlnPlot(object = object,features = marker2,pt.size = 0,sort = "increasing")


anno_epith = object

tt = as.character(anno_epith@active.ident)
tt[tt=="7"] = "Globlet cells"
tt[tt=="2"] = "Mature colonocytes type 2"
tt[tt=="6"] = "Mature colonocytes type 1"
tt[tt=="5"] = "Mature colonocytes type 1"
tt[tt=="1"] = "Intermediate"
tt[tt=="3"] = "Intermediate"
tt[tt=="0"] = "Stem-like/TA"
tt[tt=="4"] = "Stem-like/TA"


Idents(anno_epith) = tt
TSNEPlot(anno_epith,label=T)
save(anno_epith,file = "epith/epith_anno.rData")

load("epith/epith_anno.rData",verbose = T)

atp = read.csv("atp.csv",header = F)
rownames(atp) = atp[,1]
ss = anno_epith$stim
n = c(paste("SMC0",1:9,".N",sep = ""),"SMC10.N")
t = c(paste("SMC0",1:9,".T",sep = ""),"SMC10.T")
for (i in 1:length(n)) {
  ni = n[i]
  ss[ss==ni] = "Normal"
}

for (i in 1:length(t)) {
  ti = t[i]
  ss[ss==ti] = "Tumor"
}

anno_epith$stim2 = ss

VlnPlot(object = anno_epith,features = rownames(atp),split.by  = "stim2",pt.size = 0)
VlnPlot(object = anno_epith,features = rownames(mettl),split.by  = "stim2",pt.size = 0)


############## Tcell METTL9
load("tcell/tcell_anno.rData",verbose = T)
TSNEPlot(anno_t,label=T)
cd4 = subset(anno_t,idents="NK")
nor_idx = names(cd4$stim)[which(as.character(sapply(as.character(cd4$stim),function(x){strsplit(x,split = "[.]")[[1]][2]}))=="N")]
cd4_norm = subset(cd4,cell = nor_idx)
obj = cd4_norm

degs_function = function(da,idx){
  dause = da
  tmp=idx
  condition <- factor(tmp)
  colData <- data.frame(row.names=colnames(dause), condition)
  dds <- DESeqDataSetFromMatrix(dause, DataFrame(condition), design= ~ condition )
  degs = DESeq(dds)
  res <- as.matrix(results(degs))
  out = data.frame(dause,res[,c("log2FoldChange","pvalue","padj")])
  sig = rep("no",nrow(out))
  sig[which(out$log2FoldChange > 0.25 & out$pvalue < 0.05)] = "up"
  sig[which(out$log2FoldChange < -0.25 & out$pvalue < 0.05)] = "down"
  out$sig = sig
  return(out)
}


diff_obj = function(obj){
  da = as.matrix(obj@assays$RNA@counts)
  met9 = as.matrix(da["METTL9",])
  met9_pos = met9[which(met9>0),]
  met9_neg = met9[which(met9 ==0),]
  posidx = names(met9_pos)
  negidx = names(met9_neg)
  ss = as.matrix(obj@active.ident)
  ss[posidx,] = "pos"  
  ss[negidx,] = "neg"
  degs = degs_function(da = da,idx = ss)
  output = list(degs=degs,idx = ss)
}

obj2 = diff_obj(obj = obj)

degs = obj2$degs
used = degs[degs$sig != "no",646:649] 
write.csv(used,file = "degs/NK_degs.csv")

######## function
rm(list = ls())

degs = read.csv("degs/overlap/nk/nk_overlap_degs.csv")

x<- as.character(degs[,1])
eg <- bitr(x, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")

id_trans = function(da,eg){
  ss = da[,"geneID"]
  dd = lapply(ss, function(y){
    y2 = strsplit(y,split = "/")[[1]]
    ge = unique(eg[which(eg[,2] %in% y2),1])
    ge2 = paste(ge,sep = "/",collapse = "/")
    return(ge2)
  })
  dd2 = unlist(dd)
  da$geneID = dd2
  return(da)
}


ego_BP <- enrichGO(OrgDb="org.Hs.eg.db", gene = eg[,2], pvalueCutoff = 0.05, ont = "BP", readable=TRUE)
res = ego_BP@result
res2 = res[order(res[,"Count"],decreasing = T),]
write.csv(res2,file = "degs/overlap/nk/nk_overlap_degs_bp.csv")

kk <- enrichKEGG(gene = eg[,2], organism = 'hsa', pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
out = kk@result
out2 = out[order(out[,"Count"],decreasing = T),]
out3 = id_trans(da = out2,eg = eg)
write.csv(out3,file="degs/overlap/nk/nk_overlap_degs_kk.csv")


#######################normal and tumor
rm(list = ls())
load("tcell/tcell_anno.rData",verbose = T)
TSNEPlot(anno_t,label=T)
cd4 = subset(anno_t,idents="CD4+")
nor_idx = names(cd4$stim)[which(as.character(sapply(as.character(cd4$stim),function(x){strsplit(x,split = "[.]")[[1]][2]}))=="N")]
tum_idx = names(cd4$stim)[which(as.character(sapply(as.character(cd4$stim),function(x){strsplit(x,split = "[.]")[[1]][2]}))=="T")]

diff_obj2 = function(obj,nor_idx,tum_idx){
  da = as.matrix(obj@assays$RNA@counts)
  ss = as.matrix(obj@active.ident)
  ss[nor_idx,] = "normal"  
  ss[tum_idx,] = "tumor"
  degs = degs_function(da = da,idx = ss)
  output = list(degs=degs,idx = ss)
}

deg_norm = diff_obj2(obj = cd4,nor_idx = nor_idx,tum_idx = tum_idx)

save(deg_norm,file="degs/nt/cd4/cd4_deseq2.rData")

deg_norm2 = deg_norm$degs
deg_used = deg_norm2[deg_norm2$sig != "no",5278:5281] 
write.csv(deg_used,file = "degs/cd8_degs_nt.csv")

mett9_degs = read.csv(file = "degs/mettl9/cd8/cd8_degs.csv",header = T)
rownames(mett9_degs) = mett9_degs[,1]
ovge = intersect(rownames(mett9_degs),rownames(deg_used))

ovge2 = data.frame(mett9_degs[ovge,],deg_used[ovge,])
colnames(ovge2) = c("gene","log2FoldChange_mettl9","pvalue_mettl9","padj_mettl9","sig_mettl9","log2FoldChange_nt","pvalue_nt","padj_nt","sig_nt")
write.csv(ovge2,file = "nk_overlap_degs.csv")

################################################
degs = read.csv("degs/nt/cd4/cd4_degs_nt.csv",header = T)
head(degs)
id = "HSPB1"
nor_da = da[id,nor_idx]
t_da = da[id,tum_idx]
vioplot(as.numeric(nor_da),as.numeric(t_da),main = id,names = c("normal","tumor"))

load("tcell/tcell_anno.rData",verbose = T)
TSNEPlot(anno_t,label=T)

ss = anno_t$stim
n = c(paste("SMC0",1:9,".N",sep = ""),"SMC10.N")
t = c(paste("SMC0",1:9,".T",sep = ""),"SMC10.T")
for (i in 1:length(n)) {
  ni = n[i]
  ss[ss==ni] = "Normal"
}

for (i in 1:length(t)) {
  ti = t[i]
  ss[ss==ti] = "Tumor"
}

anno_t$stim2 = ss
FeaturePlot(anno_t,features = "METTL9",reduction = "tsne",label = T,split.by = "stim2")


cd4 = subset(anno_t,idents="CD4+")
nor_idx = names(cd4$stim)[which(as.character(sapply(as.character(cd4$stim),function(x){strsplit(x,split = "[.]")[[1]][2]}))=="N")]
tum_idx = names(cd4$stim)[which(as.character(sapply(as.character(cd4$stim),function(x){strsplit(x,split = "[.]")[[1]][2]}))=="T")]

####
ss = as.matrix(cd4@active.ident)
ss[nor_idx,] = "normal"
ss[tum_idx,] = "tumor"
ss = as.factor(ss)
names(ss) = names(cd4@active.ident)
cd4@active.ident = ss

DefaultAssay(cd4) = "RNA"
cd4 = ScaleData(cd4)
marker = FindAllMarkers(object = cd4,only.pos = F,min.pct = 0.5,logfc.threshold = 0,slot = "scale.data")
VlnPlot(cd4,features = "METTL9")
