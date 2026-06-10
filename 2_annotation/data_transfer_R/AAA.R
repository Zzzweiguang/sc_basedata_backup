library(Seurat) 
library(anndata)
library(reticulate)
library(dplyr)
use_python("/home/lixiangyu/anaconda3/envs/zjx_r433/bin/python")
ad <- import("anndata")
getwd()##
 
#-----------------------------------------------1-JD----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/AAA/1_JD-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/AAA/1_JD-pre_annt.rds")

#-----------------------------------------------2-ZDZJ----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/AAA/2_ZDZJ-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/AAA/2_ZDZJ-pre_annt.rds")


#-----------------------------------------------AAA----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/AAA/3_AAA-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/AAA/3_AAA-pre_annt.rds")


#-----------------------------------------------AAA-1-3LIB----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/AAA/AAA_1_3LIB-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/AAA/AAA_1_3LIB-pre_annt.rds")



#-----------------------------------------------AAA-2-3LIB_0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/AAA/AAA_2_3LIB-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/AAA/AAA_2_3LIB-pre_annt.rds")



#-----------------------------------------------AAA-3-3LIB_0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/AAA/AAA_3_3LIB-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/AAA/AAA_3_3LIB-pre_annt.rds")


#-----------------------------------------------AAA-4-3LIB_0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/AAA/AAA_4_3LIB-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/AAA/AAA_4_3LIB-pre_annt.rds")


#-----------------------------------------------AAA-9_0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/AAA/AAA_9-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/AAA/AAA_9-pre_annt.rds")

#-----------------------------------------------AAA-8----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/AAA/AAA_8-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/AAA/AAA_8-pre_annt.rds")

#-----------------------------------------------AAA-D----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/AAA/AAA_D-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/AAA/AAA_D-pre_annt.rds")


#-----------------------------------------------AAA-MAX-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/AAA/AAA_max-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/AAA/AAA_MAX-pre_annt.rds")


#-----------------------------------------------AAA-P-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/AAA/AAA_P-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/AAA/AAA_P-pre_annt.rds")


#-----------------------------------------------AAA-PRO-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/AAA/AAA_PRO-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/AAA/AAA_PR0-pre_annt.rds")



#-----------------------------------------------RRAA-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/AAA/RRAA-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/AAA/RRAA-pre_annt.rds")