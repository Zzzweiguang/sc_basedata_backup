library(Seurat)
library(Matrix)
library(jsonlite)

#.libPaths("/work/cxl230019/software/envs/my_R/lib/R/library/")
.libPaths()
library(Matrix)
library(Seurat)
library(magrittr)
library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(grid)
require(dendextend)
require(ggthemes)
library(tidyverse)
.libPaths()
#.libPaths("/home2/s231294/.conda/envs/my_R/lib/R/library")
.libPaths()
library(Matrix)
library(magrittr)
library(dplyr)
library(ggplot2)
library(grid)
require(dendextend)
require(ggthemes)
library(tidyverse)
library(Matrix)
library(magrittr)
library(dplyr)
library(ggplot2)
library(grid)
require(dendextend)
require(ggthemes)
library(tidyverse)
library(Seurat)
library(patchwork)
library(gridExtra)
library(cowplot)
library(monocle3)

library(Seurat)
library(Matrix)
library(jsonlite)

# ─────────────────────────────────────────────
# 路径配置
# ─────────────────────────────────────────────
export_dir <- "/groups/xuan/cxl230019/Data/Public_scGC/Normal/scpoli/ref_query_label_transfer/test_round5_eta20_ep150/R4_plus_CancerDiscov2022/seurat_export_R4_plus_CancerDiscov2022"

run_tag <- "R4_plus_CancerDiscov2022_eta20_ep150"

# ─────────────────────────────────────────────
# 1. 读取 counts(全基因)和 log-normalized 矩阵
# ─────────────────────────────────────────────
cat("Loading count matrix (all genes)...\n")
counts_mat <- readMM(file.path(export_dir, "data_counts_allgenes.mtx"))
barcodes   <- readLines(file.path(export_dir, "barcodes.tsv"))
genes_all  <- readLines(file.path(export_dir, "genes_all.tsv"))
rownames(counts_mat) <- genes_all
colnames(counts_mat) <- barcodes
counts_mat <- as(counts_mat, "CsparseMatrix")
cat("  counts_mat:", nrow(counts_mat), "genes ×", ncol(counts_mat), "cells\n")

cat("Loading log-normalized matrix (data_X = adata.X)...\n")
data_mat <- readMM(file.path(export_dir, "data_X.mtx"))
if (nrow(data_mat) == length(genes_all)) {
  rownames(data_mat) <- genes_all
  cat("  data_X uses full gene set\n")
} else {
  genes_hvg <- readLines(file.path(export_dir, "genes_hvg.tsv"))
  stopifnot(nrow(data_mat) == length(genes_hvg))
  rownames(data_mat) <- genes_hvg
  cat("  data_X uses HVG subset:", length(genes_hvg), "genes\n")
}
colnames(data_mat) <- barcodes
data_mat <- as(data_mat, "CsparseMatrix")

# ─────────────────────────────────────────────
# 2. meta.data
# ─────────────────────────────────────────────
cat("Loading meta.data...\n")
meta <- read.csv(file.path(export_dir, "meta_data.csv"), row.names = 1)
cat("  meta.data:", nrow(meta), "cells ×", ncol(meta), "columns\n")

# ─────────────────────────────────────────────
# 3. 创建 Seurat 对象
#    ⚠️ data slot 已经是 scanpy log-normalized,不要再跑 NormalizeData()
# ─────────────────────────────────────────────
seu <- CreateSeuratObject(
  counts    = counts_mat,
  meta.data = meta,
  assay     = "RNA",
  min.cells = 0, min.features = 0
)

# 把 log-data 塞进 data slot(若是 HVG 子集则先扩展到全基因)
if (nrow(data_mat) != nrow(counts_mat)) {
  cat("data_X is HVG subset; expanding to full gene set for data slot...\n")
  full_data <- counts_mat
  full_data[] <- 0
  full_data[rownames(data_mat), ] <- data_mat
  seu <- SetAssayData(seu, layer = "data", new.data = full_data, assay = "RNA")
} else {
  seu <- SetAssayData(seu, layer = "data", new.data = data_mat, assay = "RNA")
}

cat("Seurat object created:", ncol(seu), "cells ×", nrow(seu), "genes\n")

# ─────────────────────────────────────────────
# 4. HVG
# ─────────────────────────────────────────────
hvg <- readLines(file.path(export_dir, "hvg_list.txt"))
hvg <- intersect(hvg, rownames(seu))
VariableFeatures(seu) <- hvg
cat("HVG set:", length(hvg), "features\n")

# ─────────────────────────────────────────────
# 5. 添加降维结果
# ─────────────────────────────────────────────
add_reduction <- function(seu, csv_path, key, assay = "RNA",
                          loadings_path = NULL, stdev = NULL) {
  if (!file.exists(csv_path)) {
    cat("  [skip]", key, "-- file not found:", basename(csv_path), "\n")
    return(seu)
  }
  emb <- as.matrix(read.csv(csv_path, row.names = 1))
  colnames(emb) <- paste0(key, "_", seq_len(ncol(emb)))

  args <- list(embeddings = emb, key = paste0(key, "_"), assay = assay)
  if (!is.null(loadings_path) && file.exists(loadings_path)) {
    loadings <- as.matrix(read.csv(loadings_path, row.names = 1))
    colnames(loadings) <- paste0(key, "_", seq_len(ncol(loadings)))
    args$loadings <- loadings
  }
  if (!is.null(stdev)) args$stdev <- stdev

  seu[[key]] <- do.call(CreateDimReducObject, args)
  cat("  +", key, ":", ncol(emb), "dims\n")
  seu
}

# PCA(去批次前)
pca_var <- read.csv(file.path(export_dir, "pca_variance.csv"))
seu <- add_reduction(
  seu, file.path(export_dir, "obsm_X_pca.csv"), "pca",
  loadings_path = file.path(export_dir, "varm_PCs.csv"),
  stdev = pca_var$stdev
)

# scPoli latent(去批次后)
seu <- add_reduction(seu, file.path(export_dir, "obsm_X_scPoli.csv"), "scpoli")

# UMAP(基于 scPoli,去批次后)
seu <- add_reduction(seu, file.path(export_dir, "obsm_X_umap.csv"), "umap")

# UMAP before(基于 PCA,去批次前)
seu <- add_reduction(seu, file.path(export_dir, "obsm_X_umap_before.csv"), "umapbefore")

cat("Reductions added:", paste(names(seu@reductions), collapse = ", "), "\n")

# ─────────────────────────────────────────────
# 6. 配色(可选)
# ─────────────────────────────────────────────
colors_path <- file.path(export_dir, "colors.json")
if (file.exists(colors_path)) {
  colors <- fromJSON(colors_path)
  seu@misc$colors <- colors
  cat("Colors loaded into seu@misc$colors\n")
}

# ─────────────────────────────────────────────
# 7. 验证
# ─────────────────────────────────────────────
print(seu)
cat("\nMeta columns:\n")
print(colnames(seu@meta.data))

# 快速可视化(可选,如果在交互式环境)
# DimPlot(seu, reduction = "umap", group.by = "cell_type_pred", label = TRUE) + NoLegend()
# DimPlot(seu, reduction = "umapbefore", group.by = "journal")
# FeaturePlot(seu, features = c("LGR5", "OLFM4", "CCL25"), reduction = "umap")

# ─────────────────────────────────────────────
# 8. 保存
# ─────────────────────────────────────────────
out_file <- file.path(export_dir, paste0("seurat_scpoli_integrated_", run_tag, ".rds"))
saveRDS(seu, out_file)
cat("✓ Saved:", out_file, "\n")

# 测试 DimPlot
DimPlot(seu, reduction = "umap", group.by = "cell_type_pred", label = TRUE) 

setwd("/groups/xuan/cxl230019/Data/Public_scGC/Normal/scpoli/ref_query_label_transfer/test_round5_eta20_ep150/R4_plus_CancerDiscov2022/")
ggsave(filename = "umap_cell_type_pred.pdf", width = 14, height = 9)
# 测试 FeaturePlot(现在 Seurat 5.4.0 应该能用了)
#FeaturePlot(seu, features = c("LGR5", "OLFM4", "CCL25"), reduction = "umap")

# 用保存的 scanpy 配色画图
#DimPlot(seu, reduction = "umap", group.by = "cell_type_pred") +
#  scale_color_manual(values = colors$subtype_pred)

# 同时看两个 layer 的数值特征
cat("=== counts layer ===\n")
summary(LayerData(epi, assay = "RNA", layer = "counts")@x)

cat("\n=== data layer ===\n")
summary(LayerData(epi, assay = "RNA", layer = "data")@x)

# 是否相同(值上)?
cnt <- LayerData(epi, assay = "RNA", layer = "counts")
dat <- LayerData(epi, assay = "RNA", layer = "data")
cat("\nDimensions equal:", all(dim(cnt) == dim(dat)), "\n")
cat("Same nnz:", length(cnt@x) == length(dat@x), "\n")
if (length(cnt@x) == length(dat@x)) {
  cat("Max abs diff:", max(abs(cnt@x - dat@x)), "\n")
}

# 是不是整数?
cat("\ndata layer is integer-valued:",
    all(dat@x == round(dat@x)), "\n")

# 1. 看 counts 和 data 的维度是否一致
dim(LayerData(epi, assay = "RNA", layer = "counts"))
dim(LayerData(epi, assay = "RNA", layer = "data"))

# 2. 看两个 layer 的 feature 名是否一致
length(rownames(LayerData(epi, assay = "RNA", layer = "counts")))
length(rownames(LayerData(epi, assay = "RNA", layer = "data")))

# 3. 看 nnz 差异在哪 —— 如果 data 只是 counts 的子集(HVG),
#    那 data 的非零元素应该是 counts 的子集
cnt <- LayerData(epi, assay = "RNA", layer = "counts")
dat <- LayerData(epi, assay = "RNA", layer = "data")

# 检查 data 的非零位置是否都在 counts 的非零位置内
# (如果 data 是 counts 子集筛选后的结果,这应该成立)

epithelial_cells <- c("Chief cell", "PMC","Enteroendocrine cell","Gastric progenitor cell",
                     "Neck-chief hybrid cell","Parietal cell","Neck cell","Metaplastic stem cell")

epi <- subset(seu, idents = epithelial_cells)

table(Idents(epi))

# 先备份(可选,方便对比)
# epi[["RNA"]]$data_old <- LayerData(epi, "RNA", "data")

# 标准归一化 - 这会用 counts layer 重新生成正确的 data layer
epi <- NormalizeData(epi)

# 验证修复成功
cat("=== After normalization ===\n")
summary(LayerData(epi, "RNA", layer = "data")@x)
cat("\nIs integer-valued (should be FALSE):",
    all(LayerData(epi, "RNA", layer = "data")@x == 
        round(LayerData(epi, "RNA", layer = "data")@x)), "\n")
cat("\nRange (should be ~0 to ~8):",
    range(LayerData(epi, "RNA", layer = "data")@x), "\n")

# 定义要去除的细胞类型
remove_types <- c("Fibroblasts","Neutrophil","Plasma cell",  "Mast cell", "T cell")


# 方法2：使用Seurat的subset（推荐）
epi <- subset(epi, subset = cell_type_ori %in% remove_types, invert = TRUE)

# 验证结果
table(epi$cell_type_ori)

epi

# 提取UMAP坐标
umap_coords <- as.data.frame(epi@reductions$umap@cell.embeddings)

# 定义要保留的细胞（取反即去除条件）
# 去除条件1：umap_2 > 4
# 去除条件2：umap_2 > 10 且 umap_1 > 0
keep_cells <- rownames(umap_coords)[
  !(umap_coords$umap_2 < 3 )
]

# 子集过滤
epi <- subset(epi, cells = keep_cells)

# 验证结果
DimPlot(epi, raster=FALSE, reduction = "umap", label = TRUE, pt.size = 0.3,
        label.size = 5.5)

DimPlot(epi, raster=FALSE, reduction = "umap", pt.size = 0.3,
        label.size = 5.5,group.by="journal")

ggsave(filename = "umap_cell_type_pred_epi_journal.pdf", width = 9, height = 8)

# 查看scpoli维度数
ncol(epi@reductions$scpoli)

# 查看umap和umapbefore的维度
ncol(epi@reductions$umap)
ncol(epi@reductions$umapbefore)

# 查看各reduction的前几个值
head(epi@reductions$scpoli@cell.embeddings[, 1:3])
head(epi@reductions$umap@cell.embeddings[, 1:2])
head(epi@reductions$umapbefore@cell.embeddings[, 1:2])

df=epi
#df <- FindNeighbors(df, dims = 1:30)
df <- FindNeighbors(df, reduction = "scpoli", dims = 1:10)
df <- FindClusters(df, resolution = 3)
head(Idents(df), 5)

df <- RunUMAP(
  df,reduction = "scpoli", dims = 1:10,
  min.dist = 1.1,
  n.neighbors = 150
)
DimPlot(df, raster=FALSE,reduction = "umap", label = TRUE, pt.size = 0.4,label.size = 5.5)+
  theme(axis.text= element_text(colour= 'black',size=14),
        legend.text=element_text(colour= 'black',size=16),
        
  ) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))+
  theme(plot.title = element_text(size = 40, face = "bold"))

ggsave(filename = "umap_epi_res3_scpoli.pdf" ,width = 12,height = 9)

DimPlot(df, raster=FALSE,reduction = "umap", label = TRUE, pt.size = 0.4,label.size = 5.5, group.by = 'journal')+
  theme(axis.text= element_text(colour= 'black',size=14),
        legend.text=element_text(colour= 'black',size=16),
        
  ) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))+
  theme(plot.title = element_text(size = 40, face = "bold"))

ggsave(filename = "journal_epi_res3_scpoli.pdf" ,width = 12,height = 9)


library(Seurat)
library(FNN)
library(ggplot2)
library(patchwork)

# ===== Step 1: 用 scPoli embedding 计算每个细胞的局部纯度 =====
emb_coords <- Embeddings(epi, reduction = "scpoli")
current_idents <- as.character(Idents(epi))

K <- 30
nn_idx <- get.knn(emb_coords, k = K)$nn.index

purity <- sapply(seq_len(nrow(nn_idx)), function(i) {
  mean(current_idents[nn_idx[i, ]] == current_idents[i])
})
epi$local_purity <- purity

# ===== Step 2: 为每种细胞类型设置过滤阈值 =====
# 阈值越高 = 过滤越严格;NA = 不过滤(全部保留)
purity_thresh <- c(
  "PMC"                     = 0.3,   # 严格过滤
  "Parietal cell"           = 0.2,   # 严格过滤
  "Chief cell"              = 0.1,   # 温和过滤
  "Enteroendocrine cell"    = 0.1,
  "Gastric progenitor cell" = 0.15,
  "Neck-chief hybrid cell"  = 0.1,
  "Neck cell"               = 0.1,
  "Metaplastic stem cell"   = 0.1
)

# 决定每个细胞是否保留
keep_flag <- mapply(function(ct, p) {
  thr <- purity_thresh[ct]
  if (is.na(thr)) return(TRUE)        # 没有阈值的类型全保留
  return(p >= thr)
}, current_idents, purity)

# ===== Step 3: 报告过滤效果 =====
removed_summary <- data.frame(
  celltype = names(table(current_idents)),
  threshold = purity_thresh[names(table(current_idents))],
  before = as.numeric(table(current_idents))
)
removed_summary$after <- as.numeric(
  table(factor(current_idents[keep_flag], levels = removed_summary$celltype))
)
removed_summary$removed <- removed_summary$before - removed_summary$after
removed_summary$pct_removed <- round(
  removed_summary$removed / removed_summary$before * 100, 2
)
print(removed_summary)

cat("\nBefore:", ncol(epi), 
    "| After:", sum(keep_flag), 
    "| Removed:", sum(!keep_flag),
    sprintf("(%.2f%%)\n", sum(!keep_flag) / ncol(epi) * 100))

# ===== Step 4: subset =====
epi_clean <- subset(epi, cells = colnames(epi)[keep_flag])

# ===== Step 5: 重画 UMAP(8 张高亮图)=====
epithelial_cells <- c("Chief cell", "PMC", "Enteroendocrine cell", 
                      "Gastric progenitor cell", "Neck-chief hybrid cell",
                      "Parietal cell", "Neck cell", "Metaplastic stem cell")

plot_list <- list()
for (ct in epithelial_cells) {
  cells_highlight <- WhichCells(epi_clean, idents = ct)
  
  p <- DimPlot(epi_clean,
               reduction = "umap",
               cells.highlight = cells_highlight,
               cols.highlight = "red",
               cols = "lightgrey",
               sizes.highlight = 0.3,
               pt.size = 0.2,
               order = TRUE) +
       ggtitle(paste0(ct, " (n=", length(cells_highlight), ")")) +
       theme(legend.position = "none",
             plot.title = element_text(size = 11, hjust = 0.5),
             axis.text = element_blank(),
             axis.ticks = element_blank())
  plot_list[[ct]] <- p
}

combined <- wrap_plots(plot_list, ncol = 4)
print(combined)

ggsave("epi_umap_cleaned.pdf", combined, width = 16, height = 8)
ggsave("epi_umap_cleaned.png", combined, width = 16, height = 8, dpi = 300)


EC <- c("RFX6","SCG5","SCGN","CHGA")
msc <- c("OLFM4","LGR5","SOX9")
par <- c("ATP4A","ATP4B","GIF")
gpc <- c("TOP2A","MKI67","CENPF","CDK1")
pmc <- c("MUC5AC","TFF1")
cc <- c("CEACAM5","CEACAM6","CEACAM7")
zc_neck_hybrid_spem <- c("PGA3","PGA4","PGA5","MUC6","TFF2","CD44","LGR5","ANXA10","BHLHA15","TNFRSF19","MKI67","TOP2A")


# 合并所有基因列表，并去重
features <- unique(c(zc_neck_hybrid_spem, EC,  msc, par, gpc, pmc, cc))
df=epi
# 检查基因是否都在 Seurat 对象中
features <- features[features %in% rownames(df)]

DotPlot(df, features = features, dot.scale = 6) +
  RotatedAxis() +
  scale_size_continuous(range = c(2, 8)) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold")
  ) +
  theme(plot.title = element_text(size = 50, face = "bold"))

ggsave(filename = "dot_marker_transfer2.pdf", width = 18, height = 8)

EC <- c("RFX6","SCG5","SCGN","CHGA")
msc <- c("OLFM4","LGR5","SOX9")
par <- c("ATP4A","ATP4B","GIF")
gpc <- c("TOP2A","MKI67","CENPF","CDK1")
pmc <- c("MUC5AC","TFF1")
cc <- c("CEACAM5","CEACAM6","CEACAM7")
zc_neck_hybrid_spem <- c("PGA3","PGA4","PGA5","MUC6","TFF2","CD44","LGR5","ANXA10","BHLHA15","TNFRSF19","MKI67","TOP2A")


# 合并所有基因列表，并去重
features <- unique(c(zc_neck_hybrid_spem, EC,  msc, par, gpc, pmc, cc))

# 检查基因是否都在 Seurat 对象中
features <- features[features %in% rownames(df)]

DotPlot(df, features = features, dot.scale = 6) +
  RotatedAxis() +
  scale_size_continuous(range = c(2, 8)) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold")
  ) +
  theme(plot.title = element_text(size = 50, face = "bold"))

ggsave(filename = "dot_marker_transfer3.pdf", width = 18, height = 24)

EC <- c("RFX6","SCG5","SCGN","CHGA")
ent <- c("FABP1","FABP2","CES2","VIL1")
gob <- c("TFF3","SPINK4","MUC2")
msc <- c("OLFM4","LGR5","SOX9")
par <- c("ATP4A","ATP4B","GIF")
gpc <- c("TOP2A","MKI67","CENPF","CDK1")
pmc <- c("MUC5AC","TFF1")
cc <- c("CEACAM5","CEACAM6","CEACAM7")
zc_neck_hybrid_spem <- c("PGA3","PGA4","PGA5","MUC6","TFF2","CD44","LGR5","ANXA10","BHLHA15","TNFRSF19","MKI67","TOP2A")
TC <- c("CD2","CD3D","CD3E","CD3G","TRBC2")
SMC <- c("ACTA2","ACTG2","MYH11","SYNPO2","LMOD1")   # 把大写 C 改成 c
plasma <- c("MZB1","CD79A","ENAM","IGHA2")
BC <- c("MS4A1","BANK1")
mast <- c("TPSAB1","TPSB2","CPA3","KIT")
mac <- c("CD14","CD163","CSF1R","C1QC")
fib <- c("COL1A2","DCN","COL3A1","LUM","COL6A1","CXCL14")
endo <- c("PECAM1","VWF","ENG","PLVAP")

# 合并所有基因列表，并去重
features <- unique(c(zc_neck_hybrid_spem, EC, ent,  gob, msc, par, gpc, pmc, cc,
                     TC, SMC, plasma, BC, mast, mac, fib, endo))

# 检查基因是否都在 Seurat 对象中
features <- features[features %in% rownames(df)]

DotPlot(df, features = features, dot.scale = 6) +
  RotatedAxis() +
  scale_size_continuous(range = c(2, 8)) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold")
  ) +
  theme(plot.title = element_text(size = 50, face = "bold"))

ggsave(filename = "dot_marker_transfer1.pdf", width = 30, height = 20)

# 基本统计量
summary(df@meta.data$cell_type_uncert)

# 标准差
sd(df@meta.data$cell_type_uncert, na.rm = TRUE)

# 分位数(更细)
quantile(df@meta.data$cell_type_uncert, 
         probs = c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1), 
         na.rm = TRUE)

# NA 数量
cat("NA count:", sum(is.na(df@meta.data$cell_type_uncert)), "\n")

# 直方图
hist(df@meta.data$cell_type_uncert, 
     breaks = 50, 
     main = "Distribution of subtype_uncert",
     xlab = "Uncertainty", col = "steelblue", border = "white")


TC <- c("CD2","CD3D","CD3E","CD3G","TRBC2")
SMC <- c("ACTA2","ACTG2","MYH11","SYNPO2","LMOD1")   # 把大写 C 改成 c
plasma <- c("MZB1","CD79A","ENAM","IGHA2")
BC <- c("MS4A1","BANK1")
mast <- c("TPSAB1","TPSB2","CPA3","KIT")
mac <- c("CD14","CD163","CSF1R","C1QC")
fib <- c("COL1A2","DCN","COL3A1","LUM","COL6A1","CXCL14")
endo <- c("PECAM1","VWF","ENG","PLVAP")

# 合并所有基因列表，并去重
features <- unique(c(
                     TC, SMC, plasma, BC, mast, mac, fib, endo))

# 检查基因是否都在 Seurat 对象中
features <- features[features %in% rownames(df)]

DotPlot(df, features = features, dot.scale = 6) +
  RotatedAxis() +
  scale_size_continuous(range = c(2, 8)) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold")
  ) +
  theme(plot.title = element_text(size = 50, face = "bold"))

ggsave(filename = "dot_marker_transfer12.pdf", width = 15, height = 20)

library(dplyr)

summary_long <- df@meta.data %>%
  mutate(ident = as.character(Idents(df))) %>%
  group_by(ident, cell_type_pred) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(as.numeric(ident), desc(n_cells))

print(summary_long, n = Inf)

table(df$dataset_id[Idents(df) == 42])