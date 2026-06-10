library(Seurat)
library(Scissor)
library(Matrix)
library(readxl)

set.seed(as.integer(Sys.getenv("SCISSOR_SEED", "123")))

work_dir <- "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0603_no_Basophil_state_score"
sc_rds <- file.path(work_dir, "human_plaque_type_pred_state_score_no_threshold_var.rds")
bulk_xlsx <- file.path(work_dir, "Scissors", "NC_atlas_bulk.xlsx")
out_dir <- file.path(work_dir, "Scissors", "output_sparse_known_full")

alpha <- as.numeric(Sys.getenv("SCISSOR_ALPHA", "0.05"))
block_size <- as.integer(Sys.getenv("SCISSOR_BLOCK_SIZE", "20000"))
nlambda <- as.integer(Sys.getenv("SCISSOR_NLAMBDA", "50"))
nfolds <- as.integer(Sys.getenv("SCISSOR_NFOLDS", "5"))
use_network <- !identical(toupper(Sys.getenv("SCISSOR_USE_NETWORK", "TRUE")), "FALSE")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

get_assay_matrix <- function(seu, slot_or_layer) {
  tryCatch(
    GetAssayData(seu, assay = "RNA", layer = slot_or_layer),
    error = function(e) GetAssayData(seu, assay = "RNA", slot = slot_or_layer)
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
  phenotype <- ifelse(sample_groups == "Late Lesion", 1, 0)
  names(phenotype) <- colnames(bulk_dataset)

  list(expression = bulk_dataset, phenotype = phenotype, groups = sample_groups)
}

standardize_bulk <- function(bulk_matrix) {
  bulk_matrix <- log1p(bulk_matrix)
  bulk_z <- scale(bulk_matrix)
  bulk_z[is.na(bulk_z)] <- 0
  as.matrix(bulk_z)
}

compute_cor_block <- function(bulk_z, cell_block) {
  n_genes <- nrow(cell_block)
  cell_means <- Matrix::colMeans(cell_block)
  cell_sq_sums <- Matrix::colSums(cell_block * cell_block)
  cell_sds <- sqrt(pmax(cell_sq_sums - n_genes * cell_means^2, 0) / (n_genes - 1))
  numerator <- crossprod(bulk_z, cell_block)
  numerator <- numerator - tcrossprod(colSums(bulk_z), cell_means)
  denominator <- (n_genes - 1) * cell_sds
  cor_block <- sweep(as.matrix(numerator), 2, denominator, "/")
  cor_block[, !is.finite(denominator) | denominator == 0] <- 0
  cor_block[!is.finite(cor_block)] <- 0
  cor_block
}

compute_bulk_cell_cor <- function(bulk_matrix, cell_matrix, block_size) {
  n_cells <- ncol(cell_matrix)
  x <- matrix(0, nrow = ncol(bulk_matrix), ncol = n_cells)
  rownames(x) <- colnames(bulk_matrix)
  colnames(x) <- colnames(cell_matrix)
  bulk_z <- standardize_bulk(bulk_matrix)
  for (start in seq(1, n_cells, by = block_size)) {
    end <- min(start + block_size - 1, n_cells)
    message("Computing correlation block: ", start, "-", end, " / ", n_cells)
    x[, start:end] <- compute_cor_block(bulk_z, cell_matrix[, start:end, drop = FALSE])
    gc()
  }
  x
}

prepare_sparse_network <- function(seu, cells) {
  if (!"RNA_snn" %in% names(seu@graphs)) {
    message("RNA_snn missing; building graph on known Plaque_type cells.")
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = TRUE)
    seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = TRUE)
    n_pcs <- min(30, ncol(seu) - 1, length(VariableFeatures(seu)))
    seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = n_pcs, verbose = TRUE)
    seu <- FindNeighbors(seu, dims = seq_len(min(20, n_pcs)), verbose = TRUE)
  }
  network <- as(seu@graphs$RNA_snn, "dgCMatrix")
  network <- network[cells, cells]
  diag(network) <- 0
  network <- drop0(network)
  network@x <- rep(1, length(network@x))
  network
}

bulk <- read_bulk(bulk_xlsx)
seu <- readRDS(sc_rds)
DefaultAssay(seu) <- "RNA"
known_cells <- filter_known_plaque_type(seu@meta.data)
seu <- subset(seu, cells = known_cells)

data_matrix <- get_assay_matrix(seu, "data")
if (nrow(data_matrix) == 0 || ncol(data_matrix) == 0 || Matrix::nnzero(data_matrix) == 0) {
  message("RNA data is empty; running NormalizeData.")
  seu <- NormalizeData(seu, verbose = TRUE)
  data_matrix <- get_assay_matrix(seu, "data")
}
data_matrix <- as(data_matrix, "dgCMatrix")

common_genes <- intersect(rownames(bulk$expression), rownames(data_matrix))
if (length(common_genes) == 0) stop("No common genes between bulk and single-cell data.")

bulk_expression <- bulk$expression[common_genes, , drop = FALSE]
cell_expression <- data_matrix[common_genes, , drop = FALSE]
rm(data_matrix)
gc()

network <- NULL
if (use_network) {
  network <- prepare_sparse_network(seu, colnames(cell_expression))
}

cat("bulk samples:", ncol(bulk_expression), "\n")
cat("single cells used:", ncol(cell_expression), "\n")
cat("common genes:", length(common_genes), "\n")
print(table(bulk$groups))
cat("network mode:", ifelse(use_network, "sparse Net", "Enet without SNN"), "\n")

X <- compute_bulk_cell_cor(bulk_expression, cell_expression, block_size)
Y <- as.numeric(bulk$phenotype)
penalty <- ifelse(use_network, "Net", "Enet")

fit0 <- APML1(
  X, Y,
  family = "binomial",
  penalty = penalty,
  alpha = alpha,
  Omega = network,
  nlambda = nlambda,
  nfolds = min(nfolds, nrow(X))
)
fit1 <- APML1(
  X, Y,
  family = "binomial",
  penalty = penalty,
  alpha = alpha,
  Omega = network,
  lambda = fit0$lambda.min
)

coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])
scissor_pos <- colnames(X)[which(coefs > 0)]
scissor_neg <- colnames(X)[which(coefs < 0)]

infos1 <- list(
  para = list(alpha = alpha, lambda = fit0$lambda.min, family = "binomial", penalty = penalty),
  Coefs = coefs,
  Scissor_pos = scissor_pos,
  Scissor_neg = scissor_neg,
  fit_cv = fit0$fit
)

saveRDS(infos1, file.path(out_dir, "Scissor_result_sparse_known_full.rds"))
write.csv(
  data.frame(
    cell = colnames(X),
    scissor = ifelse(colnames(X) %in% scissor_pos, "Scissor+", ifelse(colnames(X) %in% scissor_neg, "Scissor-", "Background")),
    coefficient = coefs,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "scissor_cells_sparse_known_full.csv"),
  row.names = FALSE
)

cat("Scissor+ cells:", length(scissor_pos), "\n")
cat("Scissor- cells:", length(scissor_neg), "\n")
cat("Output:", out_dir, "\n")
