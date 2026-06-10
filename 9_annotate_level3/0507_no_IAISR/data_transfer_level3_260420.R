library(Seurat) 
library(anndata)
library(reticulate)
library(dplyr)
# use_python("/home/lixiangyu/anaconda3/envs/zr_r433/bin/python")
use_python("/home/lixiangyu/anaconda3/envs/zjx_r433/bin/python")
ad <- import("anndata")
setwd('/home/lixiangyu/zr/Annotate/ANNOTATE_new/9_annotate_level3/0507_no_IAISR')
getwd()##
 
#-----------------------------------------------EC_Arterial----------------------------------------------------

final_ad <- ad$read_h5ad("./output_human/EC_Arterial/ec_Arterial_level3_harm-pre_annot.h5ad")
counts_matrix <- t(final_ad$X)
# colnames(counts_matrix) <- final_ad$obs_names 
# rownames(counts_matrix) <- final_ad$var_names  
colnames(counts_matrix) <- make.unique(as.character(final_ad$obs_names))##去重步骤
rownames(counts_matrix) <- make.unique(gsub("_", "-", as.character(final_ad$var_names)))##去重步骤
metadata <- final_ad$obs##去重步骤
rownames(metadata) <- colnames(counts_matrix)##去重步骤
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_human/EC_Arterial/ec_Arterial_level3_harm-pre_annot.rds")


 
#-----------------------------------------------Mac_Inflammatory----------------------------------------------------

final_ad <- ad$read_h5ad("./output_human/Mac_Inflammatory/macros_infla_level3_pre_annot.h5ad")
counts_matrix <- t(final_ad$X)
# colnames(counts_matrix) <- final_ad$obs_names 
# rownames(counts_matrix) <- final_ad$var_names  
colnames(counts_matrix) <- make.unique(as.character(final_ad$obs_names))##去重步骤
rownames(counts_matrix) <- make.unique(gsub("_", "-", as.character(final_ad$var_names)))##去重步骤
metadata <- final_ad$obs##去重步骤
rownames(metadata) <- colnames(counts_matrix)##去重步骤
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_human/Mac_Inflammatory/macros_infla_level3_pre_annot.rds")



 
#-----------------------------------------------Mac_Foamy----------------------------------------------------

final_ad <- ad$read_h5ad("./output_human/Mac_Foamy/macros_foamy_level3_pre_annot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_human/Mac_Foamy/macros_foamy_level3_pre_annot.rds")



 
#-----------------------------------------------Mac_Homeostatic----------------------------------------------------

final_ad <- ad$read_h5ad("./output_human/Mac_Homeostatic/macros_homeo_level3_pre_annot.h5ad")
counts_matrix <- t(final_ad$X)
# colnames(counts_matrix) <- final_ad$obs_names 
# rownames(counts_matrix) <- final_ad$var_names  
colnames(counts_matrix) <- make.unique(as.character(final_ad$obs_names))##去重步骤
rownames(counts_matrix) <- make.unique(gsub("_", "-", as.character(final_ad$var_names)))##去重步骤
metadata <- final_ad$obs##去重步骤
rownames(metadata) <- colnames(counts_matrix)##去重步骤
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_human/Mac_Homeostatic/macros_homeo_level3_pre_annot.rds")





