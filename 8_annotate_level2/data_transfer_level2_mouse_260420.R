library(Seurat) 
library(anndata)
library(reticulate)
library(dplyr)
# use_python("/home/lixiangyu/anaconda3/envs/zr_r433/bin/python")
use_python("/home/lixiangyu/anaconda3/envs/zjx_r433/bin/python")
ad <- import("anndata")
setwd('/home/lixiangyu/zr/Annotate/ANNOTATE_new/8_annotate_level2/')
getwd()##
 

#-----------------------------------------------EC----------------------------------------------------
final_ad <- ad$read_h5ad("./output_260420/Mouse_level2/Endothelial_cell/level2_ec_pre_annot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_260420/Mouse_level2/Endothelial_cell/level2_ec_pre_annot.rds")


#-----------------------------------------------mac----------------------------------------------------
final_ad <- ad$read_h5ad("./output_260420/Mouse_level2/Macrophage/level2_macros_pre_annot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_260420/Mouse_level2/Macrophage/level2_macros_pre_annot.rds")


#-----------------------------------------------DC----------------------------------------------------
final_ad <- ad$read_h5ad("./output_260420/Mouse_level2/Dendritic_cell/level2_dc_pre_annot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_260420/Mouse_level2/Dendritic_cell/level2_dc_pre_annot.rds")

#-----------------------------------------------tcell----------------------------------------------------
final_ad <- ad$read_h5ad("./output_260420/Mouse_level2/T_cell/level2_tc_pre_annot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_260420/Mouse_level2/T_cell/level2_tc_pre_annot.rds")


#-----------------------------------------------bcell----------------------------------------------------
final_ad <- ad$read_h5ad("./output_260420/Mouse_level2/B_cell/level2_bc_pre_annot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_260420/Mouse_level2/B_cell/level2_bc_pre_annot.rds")


#-----------------------------------------------mono----------------------------------------------------
final_ad <- ad$read_h5ad("./output_260420/Mouse_level2/Monocyte/level2_mono_pre_annot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_260420/Mouse_level2/Monocyte/level2_mono_pre_annot.rds")


#-----------------------------------------------smc----------------------------------------------------
final_ad <- ad$read_h5ad("./output_260420/Mouse_level2/SMC/level2_smc_pre_annot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_260420/Mouse_level2/SMC/level2_smc_pre_annot.rds")




