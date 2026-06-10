library(Seurat) 
library(anndata)
library(reticulate)
library(dplyr)
# use_python("/home/lixiangyu/anaconda3/envs/zr_r433/bin/python")
use_python("/home/lixiangyu/anaconda3/envs/zjx_r433/bin/python")
ad <- import("anndata")
setwd('/home/lixiangyu/zr/Annotate/ANNOTATE_new/8_annotate_level2/0507_no_IAISR/')
getwd()##
 
#-----------------------------------------------ec----------------------------------------------------

final_ad <- ad$read_h5ad("./output_human/Endothelial_cell/level2_ec_pre_annot.h5ad")
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
saveRDS(object = seurat_obj, file = "./output_human/Endothelial_cell/level2_ec_pre_annot.rds")


#-----------------------------------------------mac----------------------------------------------------
final_ad <- ad$read_h5ad("./output_human/Macrophage/level2_macros_pre_annot.h5ad")
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
saveRDS(object = seurat_obj, file = "./output_human/Macrophage/level2_macros_pre_annot.rds")

#-----------------------------------------------dc----------------------------------------------------

final_ad <- ad$read_h5ad("./output_human/Dendritic_cell/level2_dc_pre_annot.h5ad")
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
saveRDS(object = seurat_obj, file = "./output_human/Dendritic_cell/level2_dc_pre_annot.rds")



#-----------------------------------------------tcell----------------------------------------------------

final_ad <- ad$read_h5ad("./output_human/T_cell/level2_tc_pre_annot.h5ad")
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
saveRDS(object = seurat_obj, file = "./output_human/T_cell/level2_tc_pre_annot.rds")



#-----------------------------------------------bcell----------------------------------------------------

final_ad <- ad$read_h5ad("./output_human/B_cell/level2_bc_pre_annot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_human/B_cell/level2_bc_pre_annot.rds")


#-----------------------------------------------Monocyte----------------------------------------------------

final_ad <- ad$read_h5ad("./output_human/Monocyte/level2_mono_pre_annot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_human/Monocyte/level2_mono_pre_annot.rds")



#---------------------------------------------------smc------------------------------------------------

final_ad <- ad$read_h5ad("./output_human/SMC/level2_smc_pre_annot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_human/SMC/level2_smc_pre_annot.rds")
