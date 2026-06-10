library(Seurat) 
library(anndata)
library(reticulate)
library(dplyr)
# use_python("/home/lixiangyu/anaconda3/envs/zr_r433/bin/python")
use_python("/home/lixiangyu/anaconda3/envs/zjx_r433/bin/python")
ad <- import("anndata")
setwd('/home/lixiangyu/zr/Annotate/ANNOTATE_new/7.5_correct/')
getwd()##
#-----------------------------------------------ec----------------------------------------------------
final_ad <- ad$read_h5ad("./output_260420/ECs/EC-reload-harm-preannot.h5ad")
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
saveRDS(object = seurat_obj, file = "./output_260420/ECs/EC-reload-harm-preannot.rds")


#-----------------------------------------------mac----------------------------------------------------
final_ad <- ad$read_h5ad("./output_260420/Macrophages/Macrophage-reload-harm-preannot.h5ad")
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
saveRDS(object = seurat_obj, file = "./output_260420/Macrophages/Macrophage-reload-harm-preannot.rds")

#-----------------------------------------------dc----------------------------------------------------

final_ad <- ad$read_h5ad("./output_260420/DCs/DC-reload-harm-preannot.h5ad")
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
saveRDS(object = seurat_obj, file = "./output_260420/DCs/DC-reload-harm-preannot.rds")


#-----------------------------------------------nk----------------------------------------------------

final_ad <- ad$read_h5ad("./output_260420/NK/NK-reload-harm-preannot.h5ad")
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
saveRDS(object = seurat_obj, file = "./output_260420/NK/NK-reload-harm-preannot.rds")


#-----------------------------------------------tcell----------------------------------------------------

final_ad <- ad$read_h5ad("./output_260420/T_cells/T_cells-reload-harm-preannot.h5ad")
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
saveRDS(object = seurat_obj, file = "./output_260420/T_cells/T_cells-reload-harm-preannot.rds")



#-----------------------------------------------bcell----------------------------------------------------

final_ad <- ad$read_h5ad("./output_260420/B_cells/B_cell-reload-harm-preannot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_260420/B_cells/B_cell-reload-harm-preannot.rds")


#-----------------------------------------------Monocyte----------------------------------------------------

final_ad <- ad$read_h5ad("./output_260420/Monocytes/Monocyte-reload-harm-preannot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_260420/Monocytes/Monocyte-reload-harm-preannot.rds")


#---------------------------------------------------neu------------------------------------------------

final_ad <- ad$read_h5ad("./output_260420/Neutrophil/Neutrophil-reload-harm-preannot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_260420/Neutrophil/Neutrophil-reload-harm-preannot.rds")


#---------------------------------------------------fibro_smc------------------------------------------------

final_ad <- ad$read_h5ad("./output_260420/Fibro_smc/fibrossmc_adata_reload-harm-preannot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_260420/Fibro_smc/fibrossmc_adata_reload-harm-preannot.rds")



#---------------------------------------------------mast------------------------------------------------

final_ad <- ad$read_h5ad("./output_260420/Mast_cells/mast_adata_reload-harm-preannot.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./output_260420/Mast_cells/mast_adata_reload-harm-preannot.rds")
