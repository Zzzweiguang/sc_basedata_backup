# ============================================================
# pseudobulk + DESeq2: unstable vs stable
# 输入：Seurat RDS
# 输出：DESeq2 差异基因结果
# gene name 使用 RNA assay meta.data 中的 original_gene_name
# ============================================================

library(Seurat)
library(DESeq2)
library(dplyr)
library(tidyr)
library(openxlsx)

# =========================
# 0. 读取 RDS 文件
# =========================

rds_file <- "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0523_no_Basophil/output_auto/human/human_plaque_type_pred_var.rds"

seu <- readRDS(rds_file)

# 查看对象
print(seu)

# 查看 metadata 列名
print(colnames(seu@meta.data))
head(seu@meta.data)

# 查看 assay
print(Assays(seu))
DefaultAssay(seu) <- "RNA"

# =========================
# 0.1 提取 original_gene_name
# =========================
# original_gene_name 在 seu[["RNA"]]@meta.data 里

gene_meta <- seu[["RNA"]]@meta.data
table(is.na(gene_meta$original_gene_names))
sum(gene_meta$original_gene_names == "", na.rm = TRUE)

head(gene_meta$original_gene_names)
head(gene_meta)
# 检查 original_gene_name 是否存在
stopifnot("original_gene_names" %in% colnames(gene_meta))

gene_map <- gene_meta %>%
  dplyr::mutate(feature_id = ensembl_id) %>%
  dplyr::select(feature_id, original_gene_names) %>%
  dplyr::distinct(feature_id, .keep_all = TRUE)

print(head(gene_map))

# =========================
# 1. 确认 metadata
# =========================
# 如果你的列名不是 sample / Plaque_type_pred，需要在这里改
# 例如：
# seu$sample <- seu$orig.ident
# seu$Plaque_type_pred <- seu$plaque_type

table(seu$Plaque_type_pred)
table(seu$sample, seu$Plaque_type_pred)

# 只保留 non-unstable 和 unstable-like
seu <- subset(
  seu,
  subset = Plaque_type_pred %in% c("non-unstable", "unstable-like")
)

# 避免下划线影响后面 separate 拆分
seu$sample <- gsub("_", "-", seu$sample)
seu$Plaque_type_pred <- gsub("_", "-", seu$Plaque_type_pred)

# 再次确认
table(seu$Plaque_type_pred)
table(seu$sample, seu$Plaque_type_pred)

# =========================
# 2. pseudobulk 聚合
# =========================
# 每个 sample 按 Plaque_type_pred 聚合
# unstable vs stable 是基于样本做比较，不是基于单细胞直接比较

pb <- AggregateExpression(
  object = seu,
  assays = "RNA",
  slot = "counts",
  group.by = c("sample", "Plaque_type_pred"),
  return.seurat = FALSE
)

pb_counts <- pb$RNA
############################

print(head(rownames(pb_counts)))
print(head(rownames(gene_meta)))
print(head(gene_meta$original_gene_names))
print(head(gene_meta$ensembl_id))

##########################


# 查看 pseudobulk 矩阵
print(dim(pb_counts))
print(head(colnames(pb_counts)))

# =========================
# 3. 构建样本信息
# =========================

sample_info <- data.frame(
  pseudobulk_id = colnames(pb_counts)
) %>%
  separate(
    col = pseudobulk_id,
    into = c("sample", "Plaque_type_pred"),
    sep = "_",
    remove = FALSE
  )

rownames(sample_info) <- sample_info$pseudobulk_id

sample_info$Plaque_type_pred <- factor(
  sample_info$Plaque_type_pred,
  levels = c("non-unstable", "unstable-like")
)

# 检查
print(sample_info)
table(sample_info$Plaque_type_pred)

# =========================
# 4. count matrix
# =========================

count_matrix <- round(pb_counts)
count_matrix <- as.matrix(count_matrix)

# 保证 count 矩阵列顺序和 sample_info 行顺序一致
count_matrix <- count_matrix[, rownames(sample_info)]

stopifnot(all(colnames(count_matrix) == rownames(sample_info)))

# =========================
# 5. DESeq2 差异分析
# =========================

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ Plaque_type_pred
)

# 过滤低表达基因
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("Plaque_type_pred", "unstable-like", "non-unstable")
)

res <- as.data.frame(res)
res$feature_id <- rownames(res)

# 加入 original_gene_name
res <- res %>%
  left_join(gene_map, by = "feature_id") %>%
  mutate(
    gene_symbol = ifelse(
      is.na(original_gene_names) | original_gene_names == "",
      feature_id,
      original_gene_names
    )
  ) %>%
  relocate(gene_symbol, original_gene_names, feature_id) %>%
  arrange(padj)

# =========================
# 6. 筛选 DEG
# =========================

deg <- res %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

unstable_up <- deg %>%
  filter(log2FoldChange > 0)

stable_up <- deg %>%
  filter(log2FoldChange < 0)

cat("总差异基因数:", nrow(deg), "\n")
cat("unstable 上调基因数:", nrow(unstable_up), "\n")
cat("stable 上调基因数:", nrow(stable_up), "\n")

# =========================
# 7. 保存结果
# =========================

outdir <- "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0523_no_Basophil/output_DESeq2/human_auto"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

write.xlsx(
  res,
  file.path(outdir, "all_genes_unstable_vs_stable.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  deg,
  file.path(outdir, "DEG_padj0.05_log2FC1.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  unstable_up,
  file.path(outdir, "unstable_up.xlsx"),
  rowNames = FALSE
)

write.xlsx(
  stable_up,
  file.path(outdir, "stable_up.xlsx"),
  rowNames = FALSE
)

cat("DESeq2 分析完成，结果已保存到:", outdir, "\n")
