library(Seurat) 
library(anndata)
library(reticulate)
library(dplyr)
use_python("/home/lixiangyu/anaconda3/envs/zjx_r433/bin/python")
ad <- import("anndata")
getwd()##
 
#-----------------------------------------------ISR-7-1----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/Atherosclerosis_ISR/ISR_7_1-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/Atherosclerosis_ISR/ISR_7_1-pre_annt.rds")



#-----------------------------------------------ISR-8----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/Atherosclerosis_ISR/ISR_8-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/Atherosclerosis_ISR/ISR_8-pre_annt.rds")




#-----------------------------------------------POP-ISR-A----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/Atherosclerosis_ISR/POP_ISR_A-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/Atherosclerosis_ISR/POP_ISR_A-pre_annt.rds")


#-----------------------------------------------LW_matrix----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/Atherosclerosis_ISR/LW_matrix-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/Atherosclerosis_ISR/LW_matrix-pre_annt.rds")



#-----------------------------------------------ZNB_matrix----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/Atherosclerosis_ISR/ZNB-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/Atherosclerosis_ISR/ZNB-pre_annt.rds")



#-----------------------------------------------IAISR-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/Atherosclerosis_ISR/IAISR-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/Atherosclerosis_ISR/IAISR-pre_annt.rds")



#-----------------------------------------------ISR_7_2-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/Atherosclerosis_ISR/ISR_7_2-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/Atherosclerosis_ISR/ISR_7_2-pre_annt.rds")



#-----------------------------------------------ISR_9-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/Atherosclerosis_ISR/ISR_9-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/Atherosclerosis_ISR/ISR_9-pre_annt.rds")


#-----------------------------------------------LCH_matrix-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/Atherosclerosis_ISR/LCH_matrix-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/Atherosclerosis_ISR/LCH_matrix-pre_annt.rds")


#-----------------------------------------------LHJ20230414_matrix-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/Atherosclerosis_ISR/LHJ20230414_matrix-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/Atherosclerosis_ISR/LHJ20230414_matrix-pre_annt.rds")


#-----------------------------------------------POP_ISR2_B-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/Atherosclerosis_ISR/POP_ISR2_B-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/Atherosclerosis_ISR/POP_ISR2_B-pre_annt.rds")



#-----------------------------------------------S221031_matrix-0420----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/Atherosclerosis_ISR/S221031_matrix-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/Atherosclerosis_ISR/S221031_matrix-pre_annt.rds")

