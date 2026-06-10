library(Seurat)
library(Scissor)
library(Matrix)
library(preprocessCore)
library(readxl)

set.seed(123)

work_dir <- "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0603_no_Basophil_state_score"
sc_rds <- file.path(work_dir, "human_plaque_type_pred_state_score_no_threshold_var.rds")
bulk_xlsx <- file.path(work_dir, "Scissors", "NC_atlas_bulk.xlsx")
out_dir <- file.path(work_dir, "Scissors", "output_human_NC_atlas_known_sampled_original")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 默认抽 10000 个 Plaque_type 已知细胞；运行前可用 SCISSOR_MAX_CELLS 修改。
# 例如：SCISSOR_MAX_CELLS=20000 Rscript Scissors/ScissorsTest_original_sampled.R
max_cells <- as.integer(Sys.getenv("SCISSOR_MAX_CELLS", "10000"))
if (is.na(max_cells) || max_cells <= 0) {
  stop("SCISSOR_MAX_CELLS must be a positive integer.")
}

get_counts <- function(seu) {
  tryCatch(
    GetAssayData(seu, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(seu, assay = "RNA", slot = "counts")
  )
}

filter_known_plaque_type <- function(meta) {
  if (!"Plaque_type" %in% colnames(meta)) {
    stop("Column 'Plaque_type' was not found in single-cell metadata.")
  }

  plaque_type <- trimws(as.character(meta$Plaque_type))
  is_known <- !is.na(plaque_type) &
    plaque_type != "" &
    !tolower(plaque_type) %in% c("unknown", "unkonwn")

  message("Plaque_type before filtering:")
  print(table(plaque_type, useNA = "ifany"))
  message("Known Plaque_type cells: ", sum(is_known), " / ", nrow(meta))

  if (sum(is_known) == 0) {
    stop("No cells remain after filtering Plaque_type != Unknown/Unkonwn.")
  }

  rownames(meta)[is_known]
}

sample_cells_by_plaque_type <- function(meta, cells, max_cells) {
  if (length(cells) <= max_cells) {
    message("No sampling needed: ", length(cells), " <= ", max_cells)
    return(cells)
  }

  plaque_type <- trimws(as.character(meta[cells, "Plaque_type"]))
  names(plaque_type) <- cells

  group_sizes <- table(plaque_type)
  raw_n <- as.numeric(group_sizes) / sum(group_sizes) * max_cells
  names(raw_n) <- names(group_sizes)
  n_take <- floor(raw_n)

  # 按小数余数补齐，保证总抽样数等于 max_cells，同时保持 Stable/Unstable 比例。
  while (sum(n_take) < max_cells) {
    candidates <- names(group_sizes)[n_take < as.integer(group_sizes)]
    if (length(candidates) == 0) break
    add_one <- candidates[which.max(raw_n[candidates] - n_take[candidates])]
    n_take[add_one] <- n_take[add_one] + 1
  }

  sampled <- unlist(
    lapply(names(n_take), function(group) {
      group_cells <- names(plaque_type)[plaque_type == group]
      sample(group_cells, size = n_take[group], replace = FALSE)
    }),
    use.names = FALSE
  )

  message("Plaque_type after stratified sampling:")
  print(table(meta[sampled, "Plaque_type"], useNA = "ifany"))
  message("Sampled cells: ", length(sampled), " / ", length(cells))

  sampled
}

compact_seurat_for_save <- function(seu) {
  assay <- seu[["RNA"]]
  if (!inherits(assay@counts, "dgCMatrix")) {
    assay@counts <- as(assay@counts, "dgCMatrix")
  }
  if (!inherits(assay@data, "dgCMatrix")) {
    assay@data <- as(assay@data, "dgCMatrix")
  }
  assay@scale.data <- matrix(nrow = 0, ncol = 0)
  seu[["RNA"]] <- assay
  seu
}

# 1. 读取单细胞数据，只保留 Plaque_type 已知的细胞，并按 Plaque_type 分层抽样。
seu <- readRDS(sc_rds)
DefaultAssay(seu) <- "RNA"

known_plaque_cells <- filter_known_plaque_type(seu@meta.data)
sampled_plaque_cells <- sample_cells_by_plaque_type(
  meta = seu@meta.data,
  cells = known_plaque_cells,
  max_cells = max_cells
)

sc_counts <- get_counts(seu)[, sampled_plaque_cells, drop = FALSE]
selected_meta <- seu@meta.data[sampled_plaque_cells, , drop = FALSE]

sc_dataset <- CreateSeuratObject(
  counts = sc_counts,
  project = "Scissor_Human_Plaque",
  min.cells = 0,
  min.features = 0,
  meta.data = selected_meta
)

rm(seu, sc_counts)
gc()

sc_dataset <- NormalizeData(sc_dataset, verbose = TRUE)
sc_dataset <- FindVariableFeatures(sc_dataset, selection.method = "vst", verbose = TRUE)
sc_dataset <- ScaleData(sc_dataset, verbose = TRUE)
sc_dataset <- RunPCA(sc_dataset, features = VariableFeatures(sc_dataset), verbose = TRUE)
sc_dataset <- FindNeighbors(sc_dataset, dims = 1:10, verbose = TRUE)
sc_dataset <- FindClusters(sc_dataset, resolution = 0.6, verbose = TRUE)
sc_dataset <- RunUMAP(sc_dataset, dims = 1:10, verbose = TRUE)

# 2. 读取 bulk 数据。
# NC_atlas_bulk.xlsx: 第 1 行是样本名，第 2 行是分组，后面是 gene expression。
bulk_raw <- as.data.frame(
  read_excel(bulk_xlsx, sheet = "bulk", col_names = FALSE, .name_repair = "minimal")
)
sample_ids <- as.character(unlist(bulk_raw[1, -1], use.names = FALSE))
sample_groups <- as.character(unlist(bulk_raw[2, -1], use.names = FALSE))

bulk_dataset <- bulk_raw[-c(1, 2), ]
gene_ids <- as.character(bulk_dataset[[1]])
bulk_dataset <- bulk_dataset[, -1, drop = FALSE]
bulk_dataset <- as.matrix(sapply(bulk_dataset, as.numeric))

rownames(bulk_dataset) <- gene_ids
colnames(bulk_dataset) <- sample_ids

keep_samples <- !is.na(sample_groups) & sample_groups %in% c("Early Lesion", "Late Lesion")
bulk_dataset <- bulk_dataset[, keep_samples, drop = FALSE]
sample_groups <- sample_groups[keep_samples]

phenotype <- ifelse(sample_groups == "Late Lesion", 1, 0)
names(phenotype) <- colnames(bulk_dataset)

common_genes <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
cat("bulk samples:", ncol(bulk_dataset), "\n")
cat("single cells used:", ncol(sc_dataset), "\n")
cat("common genes:", length(common_genes), "\n")
cat("phenotype table:\n")
print(table(sample_groups))

if (length(common_genes) == 0) {
  stop("No common genes between bulk data and single-cell data.")
}

# 3. 按原版 Scissor 调用：bulk expression + Seurat object + binary phenotype。
infos1 <- Scissor(
  bulk_dataset = bulk_dataset,
  sc_dataset = sc_dataset,
  phenotype = phenotype,
  tag = c("Early Lesion", "Late Lesion"),
  alpha = 0.05,
  family = "binomial",
  Save_file = file.path(tempdir(), "Scissor_inputs_human_NC_atlas_known_sampled.RData")
)

saveRDS(infos1, file.path(out_dir, "Scissor_result_human_NC_atlas_known_sampled.rds"))

# 4. 写回单细胞 metadata 并保存对象/图。
Scissor_select <- rep("Background", ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos1$Scissor_pos] <- "Scissor+"
Scissor_select[infos1$Scissor_neg] <- "Scissor-"

Scissor_coef <- rep(0, ncol(sc_dataset))
names(Scissor_coef) <- colnames(sc_dataset)
if (!is.null(infos1$Coefs) && length(infos1$Coefs) == ncol(sc_dataset)) {
  Scissor_coef <- as.numeric(infos1$Coefs)
  names(Scissor_coef) <- colnames(sc_dataset)
} else {
  warning("infos1$Coefs is missing or length does not match cells; scissor_coef will be 0.")
}

sc_dataset <- AddMetaData(
  sc_dataset,
  metadata = data.frame(
    scissor = Scissor_select,
    scissor_coef = Scissor_coef,
    row.names = names(Scissor_select)
  )
)
sc_dataset <- compact_seurat_for_save(sc_dataset)

saveRDS(sc_dataset, file.path(out_dir, "human_NC_atlas_known_sampled_scissor_seurat.rds"))

if ("umap" %in% Reductions(sc_dataset)) {
  pdf(file.path(out_dir, "human_NC_atlas_known_sampled_scissor_umap.pdf"), width = 8, height = 6)
  print(
    DimPlot(
      sc_dataset,
      reduction = "umap",
      group.by = "scissor",
      cols = c("Background" = "grey80", "Scissor+" = "indianred1", "Scissor-" = "royalblue"),
      pt.size = 0.2,
      order = c("Scissor-", "Scissor+")
    )
  )
  dev.off()
} else {
  message("No UMAP reduction found in the input object; skip UMAP plot.")
}

cell_meta_cols <- intersect(
  c(
    "sample",
    "dataset",
    "Plaque_type",
    "Plaque_type_prob_stable",
    "Plaque_type_prob_unstable",
    "Plaque_type_pred",
    "cell_type_level1",
    "cell_type_level1_human",
    "cell_type_level1_corrected",
    "cell_type_level2",
    "cell_type_level3",
    "seurat_clusters"
  ),
  colnames(sc_dataset@meta.data)
)

cell_table <- data.frame(
  cell = names(Scissor_select),
  scissor = Scissor_select,
  scissor_coef = Scissor_coef,
  stringsAsFactors = FALSE
)
if (length(cell_meta_cols) > 0) {
  cell_table <- cbind(
    cell_table,
    sc_dataset@meta.data[names(Scissor_select), cell_meta_cols, drop = FALSE]
  )
}

write.csv(
  cell_table,
  file.path(out_dir, "human_NC_atlas_known_sampled_scissor_cells.csv"),
  row.names = FALSE
)

cat("Scissor+ cells:", length(infos1$Scissor_pos), "\n")
cat("Scissor- cells:", length(infos1$Scissor_neg), "\n")
cat("Output:", out_dir, "\n")
