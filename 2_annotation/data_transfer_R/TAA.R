library(Seurat) 
library(anndata)
library(reticulate)
library(dplyr)
use_python("/home/lixiangyu/anaconda3/envs/zjx_r433/bin/python")
ad <- import("anndata")
getwd()##
 
#-----------------------------------------------TAA-A1-3-3LIB----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/TAA/TAA_A1_3_3LIB-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/TAA/TAA_A1_3_3LIB-pre_annt.rds")



#-----------------------------------------------TAA-Z1----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/TAA/TAA_Z1-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/TAA/TAA_Z1-pre_annt.rds")



#-----------------------------------------------TAD_AD1_2_3LIB----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/TAA/TAD_AD1_2_3LIB-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/TAA/TAD_AD1_2_3LIB-pre_annt.rds")



#-----------------------------------------------TAD1_Z3_3LIB----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/TAA/TAD1_Z3_3LIB-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/TAA/TAD1_Z3_3LIB-pre_annt.rds")



#-----------------------------------------------TAA_A1_5_3LIB-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/TAA/TAA_A1_5_3LIB-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/TAA/TAA_A1_5_3LIB-pre_annt.rds")



#-----------------------------------------------TAA_B1_5_3LIB-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/TAA/TAA_B1_5_3LIB-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/TAA/TAA_B1_5_3LIB-pre_annt.rds")



#-----------------------------------------------TAA_Z3-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/TAA/TAA_Z3-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/TAA/TAA_Z3-pre_annt.rds")



#-----------------------------------------------TAD2_Z1_3LIB-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/TAA/TAD2_Z1_3LIB-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/TAA/TAD2_Z1_3LIB-pre_annt.rds")


#-----------------------------------------------TAD2_Z2_3LIB-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/TAA/TAD2_Z2_3LIB-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/TAA/TAD2_Z2_3LIB-pre_annt.rds")




#-----------------------------------------------TAD2_Z3_3LIB-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/TAA/TAD2_Z3_3LIB-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/TAA/TAD2_Z3_3LIB-pre_annt.rds")