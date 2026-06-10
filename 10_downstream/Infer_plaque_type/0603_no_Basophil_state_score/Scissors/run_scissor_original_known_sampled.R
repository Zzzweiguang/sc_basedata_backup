library(Seurat)
library(Scissor)
library(Matrix)
library(preprocessCore)
library(readxl)

set.seed(as.integer(Sys.getenv("SCISSOR_SEED", "123")))
Sys.setenv(
  R_THREADS = Sys.getenv("R_THREADS", "1"),
  OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS", "1"),
  OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS", "1"),
  MKL_NUM_THREADS = Sys.getenv("MKL_NUM_THREADS", "1"),
  VECLIB_MAXIMUM_THREADS = Sys.getenv("VECLIB_MAXIMUM_THREADS", "1"),
  NUMEXPR_NUM_THREADS = Sys.getenv("NUMEXPR_NUM_THREADS", "1")
)

work_dir <- "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0603_no_Basophil_state_score"
sc_rds <- file.path(work_dir, "human_plaque_type_pred_state_score_no_threshold_var.rds")
bulk_xlsx <- file.path(work_dir, "Scissors", "NC_atlas_bulk.xlsx")
out_dir <- file.path(work_dir, "Scissors", "output_original_known_balanced_author_qc_logcpm")
max_cells <- as.integer(Sys.getenv("SCISSOR_MAX_CELLS", "10000"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

get_counts <- function(seu) {
  tryCatch(
    GetAssayData(seu, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(seu, assay = "RNA", slot = "counts")
  )
}

filter_known_plaque_type <- function(meta) {
  plaque_type <- trimws(as.character(meta$Plaque_type))
  is_known <- !is.na(plaque_type) &
    plaque_type != "" &
    !tolower(plaque_type) %in% c("unknown", "unkonwn")
  print(table(plaque_type, useNA = "ifany"))
  message("Known Plaque_type cells: ", sum(is_known), " / ", nrow(meta))
  rownames(meta)[is_known]
}

sample_cells <- function(meta, max_cells) {
  # 为了公平比较 Stable/Unstable，这里按 Plaque_type 平衡抽样。
  # 默认 SCISSOR_MAX_CELLS=10000，即 Stable 5000 + Unstable 5000。
  if (is.na(max_cells) || max_cells <= 0) {
    stop("SCISSOR_MAX_CELLS must be a positive integer.")
  }
  if (!"Plaque_type" %in% colnames(meta)) {
    stop("Column 'Plaque_type' was not found in single-cell metadata.")
  }

  plaque_type <- trimws(as.character(meta$Plaque_type))
  names(plaque_type) <- rownames(meta)
  cells_by_plaque <- split(rownames(meta), plaque_type, drop = TRUE)
  target_types <- intersect(c("Stable", "Unstable"), names(cells_by_plaque))
  if (length(target_types) != 2) {
    stop("Balanced sampling expects both Stable and Unstable cells in Plaque_type.")
  }

  per_type <- floor(max_cells / length(target_types))
  selected <- unlist(lapply(target_types, function(type) {
    cells <- cells_by_plaque[[type]]
    sample(cells, min(length(cells), per_type), replace = FALSE)
  }), use.names = FALSE)

  message("Balanced subsampled cells by Plaque_type: ", length(selected), " / ", nrow(meta))
  print(table(meta[selected, "Plaque_type"], useNA = "ifany"))
  selected
}

author_qc_removed_samples <- function() {
  # 作者代码里的 bulk QC 删除列表：
  # PCA outliers、mapped reads percentage < 75%、duplication > 60%、mapped reads < 10,000,000。
  samples_rm_pca <- c(
    "22L000719", "22L000728", "22L000872", "22L000884",
    "22L000886", "22L000889", "22L000893"
  )
  samples_rm_perc_mapped <- c(
    "120658", "120659", "120661", "120663", "120664", "120665", "120666", "120680",
    "120681", "120682", "120683", "120684", "120685", "130004", "22L000714", "22L000719",
    "22L000722", "22L000728", "22L000736", "22L000817", "22L000826", "22L000837",
    "22L000852", "22L000872", "22L000876", "22L000884", "22L000886", "22L000889",
    "22L000893", "22L001767", "22L001770", "22L001785"
  )
  samples_rm_pct_duplication <- "120655"
  samples_rm_num_mapped <- c(
    "22L000714", "22L000719", "22L000722", "22L000728", "22L000736",
    "22L000817", "22L000826", "22L000837", "22L000852", "22L000860",
    "22L000872", "22L000876", "22L000884", "22L000886", "22L000889",
    "22L000893", "22L001767", "22L001770", "22L001785"
  )

  unique(c(
    samples_rm_pca,
    samples_rm_perc_mapped,
    samples_rm_pct_duplication,
    samples_rm_num_mapped
  ))
}

read_bulk <- function(path) {
  bulk_raw <- as.data.frame(
    read_excel(path, sheet = "bulk", col_names = FALSE, .name_repair = "minimal")
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

  # 这些 bulk 值像 raw/estimated counts。先按作者 QC 列表删除样本，再转成 logCPM。
  lib_size <- colSums(bulk_dataset, na.rm = TRUE)
  samples_rm <- author_qc_removed_samples()
  keep_qc <- !colnames(bulk_dataset) %in% samples_rm
  message("Bulk library size summary before filtering:")
  print(summary(lib_size))
  message("Author QC removed samples listed: ", length(samples_rm))
  message("Author QC removed samples found in this bulk matrix: ", sum(!keep_qc))
  if (sum(keep_qc) == 0) {
    stop("No bulk samples remain after author QC filtering.")
  }
  if (any(!keep_qc)) {
    message("Removed author-QC-failed bulk samples:")
    print(data.frame(
      sample = colnames(bulk_dataset)[!keep_qc],
      group = sample_groups[!keep_qc],
      lib_size = lib_size[!keep_qc],
      row.names = NULL
    ))
  }
  bulk_dataset <- bulk_dataset[, keep_qc, drop = FALSE]
  sample_groups <- sample_groups[keep_qc]
  lib_size <- lib_size[keep_qc]

  # CPM 先校正样本测序深度，log2(x+1) 再压缩极高表达基因的影响。
  bulk_dataset <- t(t(bulk_dataset) / lib_size * 1e6)
  bulk_dataset <- log2(bulk_dataset + 1)

  phenotype <- ifelse(sample_groups == "Late Lesion", 1, 0)
  names(phenotype) <- colnames(bulk_dataset)

  list(expression = bulk_dataset, phenotype = phenotype, groups = sample_groups)
}

normalize_quantiles_single_thread <- function(x) {
  if (!is.matrix(x)) {
    stop("Matrix expected in normalize_quantiles_single_thread")
  }

  rows <- nrow(x)
  cols <- ncol(x)
  sorted <- apply(x, 2, sort, na.last = TRUE)
  target <- rowMeans(sorted, na.rm = TRUE)
  rm(sorted)
  gc()

  out <- matrix(NA_real_, nrow = rows, ncol = cols)
  rownames(out) <- rownames(x)
  colnames(out) <- colnames(x)

  for (j in seq_len(cols)) {
    if (j %% 1000 == 0) {
      message("Single-thread quantile normalization column: ", j, " / ", cols)
      gc()
    }

    values <- x[, j]
    ord <- order(values, na.last = NA)
    if (length(ord) == 0) {
      next
    }

    sorted_values <- values[ord]
    run_lengths <- rle(sorted_values)$lengths
    group_id <- rep(seq_along(run_lengths), run_lengths)
    group_target_sum <- rowsum(target[seq_along(ord)], group_id, reorder = FALSE)
    group_target_mean <- as.numeric(group_target_sum) / run_lengths
    out[ord, j] <- group_target_mean[group_id]
  }

  out
}

Scissor_original_binomial_qnorm_r <- function(
    bulk_dataset,
    sc_dataset,
    phenotype,
    tag,
    alpha = 0.05,
    cutoff = 0.2,
    Save_file = file.path(tempdir(), "Scissor_inputs_original_known_sampled.RData")
) {
  common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
  if (length(common) == 0) {
    stop("There is no common genes between the given single-cell and bulk samples.")
  }

  # 这里保持 Scissor 原版思路：构建 single-cell 表达矩阵和 SNN 网络。
  # 原包 Scissor() 会调用 preprocessCore::normalize.quantiles()，当前环境会 pthread 报错；
  # 所以下面只把 quantile normalization 换成单线程 R 实现，其余 APML/Net 流程保持一致。
  sc_exprs <- as.matrix(sc_dataset@assays$RNA@data)
  network <- as.matrix(sc_dataset@graphs$RNA_snn)
  diag(network) <- 0
  network[which(network != 0)] <- 1

  # bulk 和单细胞拼在一起做 quantile normalization 后，再计算 bulk-cell 相关性矩阵 X。
  dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, ])
  dataset1 <- normalize_quantiles_single_thread(dataset0)
  rownames(dataset1) <- rownames(dataset0)
  colnames(dataset1) <- colnames(dataset0)

  Expression_bulk <- dataset1[, seq_len(ncol(bulk_dataset)), drop = FALSE]
  Expression_cell <- dataset1[, (ncol(bulk_dataset) + 1):ncol(dataset1), drop = FALSE]
  X <- cor(Expression_bulk, Expression_cell)

  quality_check <- quantile(X)
  print("|**************************************************|")
  print("Performing quality-check for the correlations")
  print("The five-number summary of correlations:")
  print(quality_check)
  print("|**************************************************|")
  if (quality_check[3] < 0.01) {
    warning("The median correlation between the single-cell and bulk samples is relatively low.")
  }

  Y <- as.numeric(phenotype)
  z <- table(Y)
  if (length(z) != length(tag)) {
    stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
  }
  print(sprintf("Current phenotype contains %d %s and %d %s samples.", z[1], tag[1], z[2], tag[2]))
  print("Perform logistic regression on the given phenotypes:")

  save(X, Y, network, Expression_bulk, Expression_cell, file = Save_file)

  for (i in seq_along(alpha)) {
    set.seed(123)
    fit0 <- APML1(
      X, Y,
      family = "binomial",
      penalty = "Net",
      alpha = alpha[i],
      Omega = network,
      nlambda = 100,
      nfolds = min(10, nrow(X))
    )
    fit1 <- APML1(
      X, Y,
      family = "binomial",
      penalty = "Net",
      alpha = alpha[i],
      Omega = network,
      lambda = fit0$lambda.min
    )
    coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])
    cell1 <- colnames(X)[which(coefs > 0)]
    cell2 <- colnames(X)[which(coefs < 0)]
    percentage <- (length(cell1) + length(cell2)) / ncol(X)

    print(sprintf("alpha = %s", alpha[i]))
    print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", length(cell1), length(cell2)))
    print(sprintf("The percentage of selected cell is: %s%%", formatC(percentage * 100, format = "f", digits = 3)))
    if (percentage < cutoff) {
      break
    }
  }
  print("|**************************************************|")

  list(
    para = list(alpha = alpha[i], lambda = fit0$lambda.min, family = "binomial"),
    Coefs = coefs,
    Scissor_pos = cell1,
    Scissor_neg = cell2
  )
}

compact_seurat_for_save <- function(seu) {
  assay <- seu[["RNA"]]
  assay@counts <- as(assay@counts, "dgCMatrix")
  assay@data <- as(assay@data, "dgCMatrix")
  assay@scale.data <- matrix(nrow = 0, ncol = 0)
  seu[["RNA"]] <- assay
  seu
}

seu <- readRDS(sc_rds)
DefaultAssay(seu) <- "RNA"
known_cells <- filter_known_plaque_type(seu@meta.data)
known_meta <- seu@meta.data[known_cells, , drop = FALSE]
selected_cells <- sample_cells(known_meta, max_cells)

sc_dataset <- CreateSeuratObject(
  counts = get_counts(seu)[, selected_cells, drop = FALSE],
  project = "Scissor_Human_Plaque_Original_Sampled",
  min.cells = 0,
  min.features = 0,
  meta.data = known_meta[selected_cells, , drop = FALSE]
)
rm(seu)
gc()

sc_dataset <- NormalizeData(sc_dataset, verbose = TRUE)
sc_dataset <- FindVariableFeatures(sc_dataset, selection.method = "vst", verbose = TRUE)
sc_dataset <- ScaleData(sc_dataset, verbose = TRUE)
sc_dataset <- RunPCA(sc_dataset, features = VariableFeatures(sc_dataset), verbose = TRUE)
sc_dataset <- FindNeighbors(sc_dataset, dims = 1:10, verbose = TRUE)
sc_dataset <- FindClusters(sc_dataset, resolution = 0.6, verbose = TRUE)
sc_dataset <- RunUMAP(sc_dataset, dims = 1:10, verbose = TRUE)

bulk <- read_bulk(bulk_xlsx)
common_genes <- intersect(rownames(bulk$expression), rownames(sc_dataset))
cat("bulk samples:", ncol(bulk$expression), "\n")
cat("single cells used:", ncol(sc_dataset), "\n")
cat("common genes:", length(common_genes), "\n")
print(table(bulk$groups))
if (length(common_genes) == 0) stop("No common genes between bulk and single-cell data.")

infos1 <- Scissor_original_binomial_qnorm_r(
  bulk_dataset = bulk$expression,
  sc_dataset = sc_dataset,
  phenotype = bulk$phenotype,
  tag = c("Early Lesion", "Late Lesion"),
  alpha = 0.05,
  Save_file = file.path(tempdir(), "Scissor_inputs_original_known_sampled.RData")
)

saveRDS(infos1, file.path(out_dir, "Scissor_result_original_known_sampled.rds"))

scissor_select <- rep("Background", ncol(sc_dataset))
names(scissor_select) <- colnames(sc_dataset)
scissor_select[infos1$Scissor_pos] <- "Scissor+"
scissor_select[infos1$Scissor_neg] <- "Scissor-"

# Coefs 是每个细胞的模型系数：正值偏 Late Lesion，负值偏 Early Lesion，0 为 Background。
scissor_coef <- rep(0, ncol(sc_dataset))
names(scissor_coef) <- colnames(sc_dataset)
if (!is.null(infos1$Coefs) && length(infos1$Coefs) == ncol(sc_dataset)) {
  scissor_coef <- as.numeric(infos1$Coefs)
  names(scissor_coef) <- colnames(sc_dataset)
} else {
  warning("infos1$Coefs is missing or length does not match cells; scissor_coef will be 0.")
}

sc_dataset <- AddMetaData(
  sc_dataset,
  metadata = data.frame(
    scissor = scissor_select,
    scissor_coef = scissor_coef,
    row.names = names(scissor_select)
  )
)

write.csv(
  data.frame(
    cell = names(scissor_select),
    scissor = scissor_select,
    scissor_coef = scissor_coef,
    Plaque_type = sc_dataset@meta.data[names(scissor_select), "Plaque_type"],
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "scissor_cells_original_known_sampled.csv"),
  row.names = FALSE
)
saveRDS(compact_seurat_for_save(sc_dataset), file.path(out_dir, "scissor_seurat_original_known_sampled.rds"))

cat("Scissor+ cells:", length(infos1$Scissor_pos), "\n")
cat("Scissor- cells:", length(infos1$Scissor_neg), "\n")
cat("Output:", out_dir, "\n")
