library(Seurat) 
library(anndata)
library(reticulate)
library(dplyr)
use_python("/home/lixiangyu/anaconda3/envs/zjx_r433/bin/python")
ad <- import("anndata")
getwd()##
 
#-----------------------------------------------cellxgene----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/Healthy/cellxgene-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/Healthy/cellxgene-pre_annt.rds")


#-----------------------------------------------Bleckwehl_et_al----------------------------------------------------
final_ad <- ad$read_h5ad("./data_transfer_R/output_data/Healthy/Bleckwehl_et_al-pre_annt.h5ad")
counts_matrix <- t(final_ad$X)
colnames(counts_matrix) <- final_ad$obs_names 
rownames(counts_matrix) <- final_ad$var_names  
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix, assay = "RNA",meta.data = final_ad$obs )
print(head(rownames(seurat_obj), 10))
print(head(colnames(seurat_obj), 10))
saveRDS(object = seurat_obj, file = "./data_transfer_R/output_data/Healthy/Bleckwehl_et_al-pre_annt.rds")
adata <- readRDS("./data_transfer_R/output_data/Healthy/Bleckwehl_et_al-pre_annt.rds")
level1_marker <- list(
    "B cell" = c('CD79A', 'CD79B', 'MS4A1', 'IGKC', 'CD22', 'FCER2'),
    'T cell'= c('CD2', 'TRAC', 'CD69', 'CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'EOMES', 'LAG3'),
    'Natural killer cell'= c('NKG7', 'XCL1', 'CTSW', 'XCL2', 'CD160', 'FCGR3A', 'PRF1', 'GNLY'),

    'Dendritic cell'= c('CLEC10A', 'FCER1A', 'CD1C', 'HLA-DRA', 'HLA-DRB1'),
    'Macrophage'= c('C1QA', 'C1QB', 'C1QC', 'CD74', 'CXCL8', 'AIF1', 'CD68', 'ITGAM', 'CSF1R', 'HLA-DRA', 'LGALS3','CD163'),##CD163 from cc,del CD14
    'Monocyte'= c('FCN1', 'S100A8', 'S100A9', 'S100A12', 'VCAN', 'CD52', 'LYZ', 'CTSS','CD14'),## CD14 from CC
    'Mast cell'= c('TPSAB1', 'TPSB2',  'HDC', 'CMA1') ,
    'Erythrocyte/Erythroid'= c('HBB', 'HBA1', 'HBA2', 'ALAS2', 'AHSP', 'SLC4A1', 'GYPA', 'KLF1', 'TMCC2'), 
    'Neutrophil' = c('NAMPT','IFITM2','G0S2','CXCL8','NEAT1','SRGN','AQP9','SOD2','FCGR3B','IVNS1ABP'),
    # 'Basophil'= c('IL3RA', 'FCER1A', 'MS4A2', 'HDC', 'GATA1', 'PRSS33', 'ENPP3'),
    'Basophil' =  c('TPSB2', 'CPA3', 'SLC24A3', 'FER', 'RP11-779O18.3', 'KIT', 'HPGDS', 'SYTL3', 'MAML3', 'ELL2', 'RP11-139E24.1', 'CCDC200', 'AKAP13', 'AREG', 'RHOH', 'LRMDA', 'ARID1B', 'IRAK3', 'TEX14', 'HPGD', 'ERCC1', 'CTNNBL1', 'LINC00486_ENSG00000230876', 'ZBTB20'),

    'Endothelial cell'= c('PECAM1', 'VWF', 'FABP4', 'CLDN5', 'IFI27', 'ECSCR', 'DYSF', 'CD34', 'COL4A1', 'COL4A2', 'SPARCL1', 'PLVAP', 'MPZL2', 'SULF1', 'EDN1'),

    'Fibroblast'= c('LUM', 'DCN', 'COL1A1', 'COL1A2', 'FBLN1', 'THY1'), 
    'Smooth muscle cell'= c('ACTA2', 'MYH11', 'MYL9', 'TPM2', 'CALD1', 'TAGLN', 'TNFRSF11B', 'LUM', 'APOE', 'APOC1', 'AGT', 'NOTCH3', 'PDGFRB', 'MFAP4'), 
    
    # 'Pericyte'= c('RGS5', 'PDGFRB', 'CSPG4', 'KCNJ8', 'ABCC9', 'NOTCH3', 'MCAM', 'DES', 'MYLK', 'ACTN1')
    'Pericyte'= c('TAGLN', 'LPP', 'CALD1', 'TPM2', 'MYL9', 'ACTA2', 'MAP1B', 'PRKG1', 'IGFBP5', 'SYNPO2', 'EPS8', 'TIMP3', 'LMOD1', 'C11orf96', 'INPP4B', 'NOTCH3', 'EBF1', 'STEAP4', 'MT-RNR1', 'CRISPLD2', 'SOX5', 'PPP1R14A', 'FILIP1L', 'LHFPL6', 'PTPRG')
)
str(level1_marker)

sc <- AddModuleScore(object = adata,features = level1_marker,ctrl = 100,name = "ModuleScore_",assay = "RNA",slot = "counts")
module_score_cols <- grep("ModuleScore_", colnames(sc@meta.data), value = TRUE)
#按聚类分组，计算每个模块的平均得分
cluster_module_score_table <- sc@meta.data %>%group_by(leiden) %>%summarise(across(all_of(module_score_cols), mean), .groups = "drop") %>%
  rename(Cluster = leiden)  
cluster_top_module <- apply(cluster_module_score_table[, -1], 1, function(x) colnames(cluster_module_score_table[, -1])[which.max(x)])
cluster_top_celltype <- gsub("ModuleScore_", "", cluster_top_module)
cluster_module_score_table$Top_Matched_CellType <- cluster_top_celltype
write.csv(cluster_module_score_table,file = "./data_transfer_R/output_data/Healthy/Bleckwehl_et_al-pre_annt.csv",row.names = FALSE, )

