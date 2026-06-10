library(Seurat) 
library(anndata)
library(reticulate)
library(dplyr)
# use_python("/home/lixiangyu/anaconda3/envs/zr_r433/bin/python")
use_python("/home/lixiangyu/anaconda3/envs/zjx_r433/bin/python")
ad <- import("anndata")
# outdir = '/home/lixiangyu/zr/Annotate/ANNOTATE_new/2_annotation/annot_mouse_0518/output/'
outdir = '/home/lixiangyu/zr/Annotate/ANNOTATE_new/2_annotation/annot_mouse_0518/output_new_marker/'
getwd()##
#-----------------------------------------------GSE141031----------------------------------------------------
final_ad <- ad$read_h5ad(file.path(outdir,"GSE141031/GSE141031_pre_annt.h5ad"))
counts_matrix <- t(final_ad$X)
# colnames(counts_matrix) <- final_ad$obs_names 
# rownames(counts_matrix) <- final_ad$var_names  
colnames(counts_matrix) <- make.unique(as.character(final_ad$obs_names))##去重步骤
rownames(counts_matrix) <- make.unique(gsub("_", "-", as.character(final_ad$var_names)))##去重步骤
metadata <- final_ad$obs##去重步骤
rownames(metadata) <- colnames(counts_matrix)##去重步骤
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = metadata )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(seurat_obj,file = file.path(outdir,"GSE141031/GSE141031_pre_annt.rds"))


#-----------------------------------------------GSE233625----------------------------------------------------
final_ad <- ad$read_h5ad(file.path(outdir,"GSE233625/GSE233625_pre_annt.h5ad"))
counts_matrix <- t(final_ad$X)
# colnames(counts_matrix) <- final_ad$obs_names 
# rownames(counts_matrix) <- final_ad$var_names  
colnames(counts_matrix) <- make.unique(as.character(final_ad$obs_names))##去重步骤
rownames(counts_matrix) <- make.unique(gsub("_", "-", as.character(final_ad$var_names)))##去重步骤
metadata <- final_ad$obs##去重步骤
rownames(metadata) <- colnames(counts_matrix)##去重步骤
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = metadata )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(seurat_obj,file = file.path(outdir,"GSE233625/GSE233625_pre_annt.rds"))

#-----------------------------------------------GSE237067----------------------------------------------------
final_ad <- ad$read_h5ad(file.path(outdir,"GSE237067/GSE237067_pre_annt.h5ad"))
counts_matrix <- t(final_ad$X)
# colnames(counts_matrix) <- final_ad$obs_names 
# rownames(counts_matrix) <- final_ad$var_names  
colnames(counts_matrix) <- make.unique(as.character(final_ad$obs_names))##去重步骤
rownames(counts_matrix) <- make.unique(gsub("_", "-", as.character(final_ad$var_names)))##去重步骤
metadata <- final_ad$obs##去重步骤
rownames(metadata) <- colnames(counts_matrix)##去重步骤
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = metadata )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(seurat_obj,file = file.path(outdir,"GSE237067/GSE237067_pre_annt.rds"))

#-----------------------------------------------GSE269449----------------------------------------------------
final_ad <- ad$read_h5ad(file.path(outdir,"GSE269449/GSE269449_pre_annt.h5ad"))
counts_matrix <- t(final_ad$X)
# colnames(counts_matrix) <- final_ad$obs_names 
# rownames(counts_matrix) <- final_ad$var_names 
colnames(counts_matrix) <- make.unique(as.character(final_ad$obs_names))##去重步骤
rownames(counts_matrix) <- make.unique(gsub("_", "-", as.character(final_ad$var_names)))##去重步骤
metadata <- final_ad$obs##去重步骤
rownames(metadata) <- colnames(counts_matrix)##去重步骤 
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = metadata )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(seurat_obj,file = file.path(outdir,"GSE269449/GSE269449_pre_annt.rds"))


#-----------------------------------------------GSE284253----------------------------------------------------
final_ad <- ad$read_h5ad(file.path(outdir,"GSE284253/GSE284253_pre_annt.h5ad"))
counts_matrix <- t(final_ad$X)
# colnames(counts_matrix) <- final_ad$obs_names 
# rownames(counts_matrix) <- final_ad$var_names 
colnames(counts_matrix) <- make.unique(as.character(final_ad$obs_names))##去重步骤
rownames(counts_matrix) <- make.unique(gsub("_", "-", as.character(final_ad$var_names)))##去重步骤
metadata <- final_ad$obs##去重步骤
rownames(metadata) <- colnames(counts_matrix)##去重步骤 
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = metadata )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(seurat_obj,file = file.path(outdir,"GSE284253/GSE284253_pre_annt.rds"))


#-----------------------------------------------GSE279601----------------------------------------------------
final_ad <- ad$read_h5ad(file.path(outdir,"GSE279601/GSE279601_pre_annt.h5ad"))
counts_matrix <- t(final_ad$X)
# colnames(counts_matrix) <- final_ad$obs_names 
# rownames(counts_matrix) <- final_ad$var_names  
colnames(counts_matrix) <- make.unique(as.character(final_ad$obs_names))##去重步骤
rownames(counts_matrix) <- make.unique(gsub("_", "-", as.character(final_ad$var_names)))##去重步骤
metadata <- final_ad$obs##去重步骤
rownames(metadata) <- colnames(counts_matrix)##去重步骤 
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = metadata )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(seurat_obj,file = file.path(outdir,"GSE279601/GSE279601_pre_annt.rds"))

#-----------------------------------------------GSE155513----------------------------------------------------
final_ad <- ad$read_h5ad(file.path(outdir,"GSE155513/GSE155513_pre_annt.h5ad"))
counts_matrix <- t(final_ad$X)
# colnames(counts_matrix) <- final_ad$obs_names 
# rownames(counts_matrix) <- final_ad$var_names  
colnames(counts_matrix) <- make.unique(as.character(final_ad$obs_names))##去重步骤
rownames(counts_matrix) <- make.unique(gsub("_", "-", as.character(final_ad$var_names)))##去重步骤
metadata <- final_ad$obs##去重步骤
rownames(metadata) <- colnames(counts_matrix)##去重步骤 
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = metadata )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(seurat_obj,file = file.path(outdir,"GSE155513/GSE155513_pre_annt.rds"))


#-----------------------------------------------GSE239591----------------------------------------------------
final_ad <- ad$read_h5ad(file.path(outdir,"GSE239591/GSE239591_pre_annt.h5ad"))
counts_matrix <- t(final_ad$X)
# colnames(counts_matrix) <- final_ad$obs_names 
# rownames(counts_matrix) <- final_ad$var_names  
colnames(counts_matrix) <- make.unique(as.character(final_ad$obs_names))##去重步骤
rownames(counts_matrix) <- make.unique(gsub("_", "-", as.character(final_ad$var_names)))##去重步骤
metadata <- final_ad$obs##去重步骤
rownames(metadata) <- colnames(counts_matrix)##去重步骤 
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = metadata )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(seurat_obj,file = file.path(outdir,"GSE239591/GSE239591_pre_annt.rds"))
