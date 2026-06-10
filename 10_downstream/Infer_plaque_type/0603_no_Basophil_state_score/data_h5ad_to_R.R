library(Seurat) 
library(anndata)
library(reticulate)
library(dplyr)
use_python("/home/lixiangyu/anaconda3/envs/zjx_r433/bin/python")
ad <- import("anndata")
getwd()##
 ########################################################固定marker#####################################################
final_ad <- ad$read_h5ad("/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0603_no_Basophil_state_score/human_plaque_type_pred_ml_marker_score_no_threshold.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
# saveRDS(object = seurat_obj, file = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0603_no_Basophil_state_score/human_plaque_type_pred_state_score_no_threshold.rds")
print(seurat_obj)

var_meta <- as.data.frame(final_ad$var)
rownames(var_meta) <- rownames(seurat_obj[["RNA"]])
seurat_obj[["RNA"]] <- SeuratObject::AddMetaData(object = seurat_obj[["RNA"]],metadata = var_meta)
head(seurat_obj[["RNA"]]@meta.data)
dim(seurat_obj@meta.data)
saveRDS(object = seurat_obj, file = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0603_no_Basophil_state_score/human_plaque_type_pred_ml_marker_score_no_threshold_var.rds")

##############################################mouse############################################
final_ad <- ad$read_h5ad("//home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0603_no_Basophil_state_score/mouse_plaque_type_pred_ml_marker_score_no_threshold.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
# saveRDS(object = seurat_obj, file = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0603_no_Basophil_state_score/mouse_plaque_type_pred_no_threshold.rds")
print(seurat_obj)
var_meta <- as.data.frame(final_ad$var)
rownames(var_meta) <- rownames(seurat_obj[["RNA"]])
seurat_obj[["RNA"]] <- SeuratObject::AddMetaData(object = seurat_obj[["RNA"]],metadata = var_meta)
head(seurat_obj[["RNA"]]@meta.data)
dim(seurat_obj@meta.data)
saveRDS(object = seurat_obj, file = "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0603_no_Basophil_state_score/mouse_plaque_type_pred_ml_marker_score_no_threshold_var.rds")


