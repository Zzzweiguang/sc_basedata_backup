library(Seurat) 
library(anndata)
library(reticulate)
library(dplyr)
use_python("/home/lixiangyu/anaconda3/envs/zjx_r433/bin/python")
ad <- import("anndata")
getwd()##
 
#-----------------------------------------------CBT----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/CBT/CBT-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/CBT/CBT-pre_annt.rds")