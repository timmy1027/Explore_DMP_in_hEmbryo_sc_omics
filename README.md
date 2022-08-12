# Explore DMP in hEmbryo single-cell omics
CS12 to CS16 human embryo scRNA-seq data is from GSE157329

##Step1
Loop loading samples: DMP, BAd1, BAd3, BAd6, SKMd1, SKMd3, SKMd6

##Step2
Rename colname of metadata as:
seq_folder | orig.ident,
nUMI | nCount_RNA,
nGene | nFeature_RNA

save merged Seurat project in "data/combined_raw_seurat.RData"  

To exclude low-quality cells, cells with fewer than 500 and over 60,000 UMIs were discarded. Cells with total nGenes less than 250 were
excluded. RNA complexity with log10GenesPerUMI > 0.8 were saved. Cells with high mitochondral DNA contamination, mitoRatio >= 20%, 
were excluded.  

Non-zero genes were defined as that a gene must be expressed in at least 10 cells  

Filtered Seurat was created by "filtered_counts"" and "filtered_seurat@meta.data"  

save merged Seurat project in "data/seurat_filtered.RData"
