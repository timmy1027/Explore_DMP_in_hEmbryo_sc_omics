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


























