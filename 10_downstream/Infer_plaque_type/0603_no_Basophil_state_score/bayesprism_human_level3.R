#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

default_config <- function() {
  work_dir <- "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0603_no_Basophil_state_score"

  list(
    work_dir = work_dir,
    sc_rds = file.path(work_dir, "human_plaque_type_pred_state_score_no_threshold_var.rds"),
    bulk_xlsx = file.path(work_dir, "Scissors", "NC_atlas_bulk.xlsx"),
    output_dir = file.path(work_dir, "BayesPrism", "output_human_NC_atlas_level3"),

    assay = "RNA",
    count_layer = "counts",
    cell_type_label_col = "cell_type_level1_corrected",
    cell_state_label_col = "cell_type_level3",
    sc_gene_id_col = "ensembl_id",

    bulk_sheet = "bulk",
    bulk_sample_row = 1L,
    bulk_group_row = 2L,
    bulk_gene_col = 1L,
    group_keep = c("Early Lesion", "Late Lesion"),

    # Keep this moderate: the source Seurat object has >1M cells.
    max_cells_per_type = 1000L,
    seed = 123L,
    n_cores = 8L,

    input_type = "count.matrix",
    round_bulk_counts = TRUE,
    force_dense_reference = FALSE,
    force_dense_mixture = TRUE,

    key = NULL,
    outlier_cut = 0.01,
    outlier_fraction = 0.10,
    min_common_genes = 1000L,
    clr_pseudocount = 1e-6,
    exclude_author_qc_samples = TRUE,

    dry_run = FALSE
  )
}

author_qc_samples <- function() {
  samples_rm_pca <- c("22L000719", "22L000728", "22L000872", "22L000884", "22L000886", "22L000889", "22L000893")
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

  unique(c(samples_rm_pca, samples_rm_perc_mapped, samples_rm_pct_duplication, samples_rm_num_mapped))
}

parse_bool <- function(value) {
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

parse_args <- function(args) {
  cfg <- default_config()

  for (arg in args) {
    if (identical(arg, "--dry-run")) {
      cfg$dry_run <- TRUE
    } else if (grepl("^--", arg) && grepl("=", arg, fixed = TRUE)) {
      key <- sub("^--([^=]+)=.*$", "\\1", arg)
      value <- sub("^--[^=]+=", "", arg)

      if (key == "sc-rds") cfg$sc_rds <- value
      else if (key == "bulk-xlsx") cfg$bulk_xlsx <- value
      else if (key == "output-dir") cfg$output_dir <- value
      else if (key == "cell-type-col") {
        cfg$cell_type_label_col <- value
        cfg$cell_state_label_col <- value
      }
      else if (key == "cell-type-label-col") cfg$cell_type_label_col <- value
      else if (key == "cell-state-label-col") cfg$cell_state_label_col <- value
      else if (key == "sc-gene-id-col") cfg$sc_gene_id_col <- value
      else if (key == "bulk-sheet") cfg$bulk_sheet <- value
      else if (key == "group-keep") cfg$group_keep <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
      else if (key == "max-cells-per-type") cfg$max_cells_per_type <- as.integer(value)
      else if (key == "n-cores") cfg$n_cores <- as.integer(value)
      else if (key == "seed") cfg$seed <- as.integer(value)
      else if (key == "round-bulk-counts") cfg$round_bulk_counts <- parse_bool(value)
      else if (key == "force-dense-reference") cfg$force_dense_reference <- parse_bool(value)
      else if (key == "force-dense-mixture") cfg$force_dense_mixture <- parse_bool(value)
      else if (key == "exclude-author-qc-samples") cfg$exclude_author_qc_samples <- parse_bool(value)
      else if (key == "key") cfg$key <- if (nzchar(value)) value else NULL
      else stop("Unknown argument: ", arg)
    } else {
      stop("Unknown argument: ", arg)
    }
  }

  cfg
}

check_packages <- function(packages) {
  installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  missing <- packages[!installed]

  if (length(missing) > 0L) {
    msg <- paste0("Missing R package(s): ", paste(missing, collapse = ", "))
    if ("BayesPrism" %in% missing) {
      msg <- paste0(
        msg,
        "\nInstall BayesPrism in the R environment you will use, for example:\n",
        "  install.packages('remotes')\n",
        "  remotes::install_github('Danko-Lab/BayesPrism')"
      )
    }
    stop(msg, call. = FALSE)
  }
}

clean_gene_ids <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.[0-9]+$", "", x)
  x
}

get_assay_counts <- function(seu, assay, count_layer) {
  counts <- tryCatch(
    Seurat::GetAssayData(seu, assay = assay, layer = count_layer),
    error = function(e) Seurat::GetAssayData(seu, assay = assay, slot = count_layer)
  )

  if (!inherits(counts, "Matrix")) {
    counts <- Matrix::Matrix(as.matrix(counts), sparse = TRUE)
  }
  counts
}

read_bulk_matrix <- function(cfg) {
  message("Reading bulk matrix: ", cfg$bulk_xlsx)

  bulk_raw <- as.data.frame(
    readxl::read_excel(
      cfg$bulk_xlsx,
      sheet = cfg$bulk_sheet,
      col_names = FALSE,
      .name_repair = "minimal"
    ),
    stringsAsFactors = FALSE
  )

  if (nrow(bulk_raw) <= max(cfg$bulk_sample_row, cfg$bulk_group_row) ||
      length(setdiff(seq_len(ncol(bulk_raw)), cfg$bulk_gene_col)) == 0L) {
    stop("Bulk sheet is too small for configured sample/group/gene rows.")
  }

  sample_cols <- setdiff(seq_len(ncol(bulk_raw)), cfg$bulk_gene_col)
  sample_ids <- trimws(as.character(unlist(bulk_raw[cfg$bulk_sample_row, sample_cols], use.names = FALSE)))
  sample_groups <- trimws(as.character(unlist(bulk_raw[cfg$bulk_group_row, sample_cols], use.names = FALSE)))

  if (any(is.na(sample_ids) | sample_ids == "")) {
    stop("Bulk sample row contains empty sample IDs.")
  }
  if (anyDuplicated(sample_ids)) {
    stop("Bulk sample IDs are duplicated; please resolve duplicates before BayesPrism.")
  }

  expr_rows <- setdiff(seq_len(nrow(bulk_raw)), c(cfg$bulk_sample_row, cfg$bulk_group_row))
  expr_df <- bulk_raw[expr_rows, , drop = FALSE]
  gene_ids <- clean_gene_ids(expr_df[[cfg$bulk_gene_col]])
  expr_values <- expr_df[, sample_cols, drop = FALSE]
  expr_values <- as.data.frame(
    lapply(expr_values, function(x) suppressWarnings(as.numeric(as.character(x)))),
    check.names = FALSE
  )
  expr_mat <- as.matrix(expr_values)

  valid_gene <- !is.na(gene_ids) & gene_ids != ""
  expr_mat <- expr_mat[valid_gene, , drop = FALSE]
  gene_ids <- gene_ids[valid_gene]
  rownames(expr_mat) <- gene_ids
  colnames(expr_mat) <- sample_ids

  if (anyNA(expr_mat)) {
    message("Bulk matrix contains NA values; setting them to 0.")
    expr_mat[is.na(expr_mat)] <- 0
  }
  if (any(expr_mat < 0, na.rm = TRUE)) {
    stop("Bulk matrix contains negative values, which BayesPrism cannot use as counts.")
  }

  if (anyDuplicated(rownames(expr_mat))) {
    message("Collapsing duplicated bulk gene IDs by summing expression values.")
    expr_mat <- rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE)
  }

  keep_samples <- !is.na(sample_groups) & sample_groups != ""
  if (length(cfg$group_keep) > 0L) {
    keep_samples <- keep_samples & sample_groups %in% cfg$group_keep
  }
  expr_mat <- expr_mat[, keep_samples, drop = FALSE]
  sample_groups <- sample_groups[keep_samples]

  if (ncol(expr_mat) == 0L) {
    stop("No bulk samples remain after group filtering.")
  }

  if (cfg$round_bulk_counts) {
    expr_mat <- round(expr_mat)
  }

  keep_gene <- rowSums(expr_mat) > 0
  expr_mat <- expr_mat[keep_gene, , drop = FALSE]

  mixture <- t(expr_mat)
  storage.mode(mixture) <- "numeric"

  sample_info <- data.frame(
    sample = rownames(mixture),
    group = factor(sample_groups, levels = unique(c(cfg$group_keep, sample_groups))),
    stringsAsFactors = FALSE
  )

  list(mixture = mixture, sample_info = sample_info)
}

downsample_cells <- function(cell_types, max_cells_per_type, seed) {
  set.seed(seed)
  idx_by_type <- split(seq_along(cell_types), cell_types)
  selected <- unlist(
    lapply(idx_by_type, function(idx) {
      if (length(idx) <= max_cells_per_type) {
        idx
      } else {
        sort(sample(idx, max_cells_per_type))
      }
    }),
    use.names = FALSE
  )
  sort(selected)
}

clean_label_values <- function(x) {
  labels <- trimws(as.character(x))
  labels[labels %in% c("", "NA", "NaN", "nan")] <- NA_character_
  labels
}

make_reference_counts <- function(cell_type_labels, cell_state_labels, selected_idx) {
  input_df <- data.frame(
    cell_type = cell_type_labels,
    cell_state = cell_state_labels,
    stringsAsFactors = FALSE
  )
  selected_df <- input_df[selected_idx, , drop = FALSE]

  n_input <- aggregate(
    rep(1L, nrow(input_df)),
    by = list(cell_type = input_df$cell_type, cell_state = input_df$cell_state),
    FUN = length
  )
  colnames(n_input)[3] <- "n_input"

  n_selected <- aggregate(
    rep(1L, nrow(selected_df)),
    by = list(cell_type = selected_df$cell_type, cell_state = selected_df$cell_state),
    FUN = length
  )
  colnames(n_selected)[3] <- "n_selected"

  reference_counts <- merge(n_input, n_selected, by = c("cell_type", "cell_state"), all.x = TRUE, sort = TRUE)
  reference_counts$n_selected[is.na(reference_counts$n_selected)] <- 0L
  reference_counts
}

prepare_sc_reference <- function(cfg) {
  message("Reading single-cell reference: ", cfg$sc_rds)

  suppressPackageStartupMessages(library(Seurat))
  seu <- readRDS(cfg$sc_rds)
  Seurat::DefaultAssay(seu) <- cfg$assay

  if (!cfg$cell_type_label_col %in% colnames(seu@meta.data)) {
    stop("Cell type label column not found in Seurat metadata: ", cfg$cell_type_label_col)
  }
  if (!cfg$cell_state_label_col %in% colnames(seu@meta.data)) {
    stop("Cell state label column not found in Seurat metadata: ", cfg$cell_state_label_col)
  }

  counts <- get_assay_counts(seu, cfg$assay, cfg$count_layer)
  meta <- seu@meta.data

  cell_type_labels <- clean_label_values(meta[[cfg$cell_type_label_col]])
  cell_state_labels <- clean_label_values(meta[[cfg$cell_state_label_col]])

  valid_cells <- !is.na(cell_type_labels) & !is.na(cell_state_labels) & rownames(meta) %in% colnames(counts)
  meta <- meta[valid_cells, , drop = FALSE]
  cell_type_labels <- cell_type_labels[valid_cells]
  cell_state_labels <- cell_state_labels[valid_cells]

  selected_idx <- downsample_cells(cell_state_labels, cfg$max_cells_per_type, cfg$seed)
  selected_cells <- rownames(meta)[selected_idx]
  selected_type_labels <- cell_type_labels[selected_idx]
  selected_state_labels <- cell_state_labels[selected_idx]

  gene_meta <- seu[[cfg$assay]]@meta.data
  if (cfg$sc_gene_id_col %in% colnames(gene_meta)) {
    gene_ids <- clean_gene_ids(gene_meta[[cfg$sc_gene_id_col]])
  } else {
    message("sc_gene_id_col not found; using assay rownames as gene IDs: ", cfg$sc_gene_id_col)
    gene_ids <- clean_gene_ids(rownames(counts))
  }

  if (length(gene_ids) != nrow(counts)) {
    stop("Gene metadata length does not match count matrix rows.")
  }

  valid_genes <- !is.na(gene_ids) & gene_ids != ""
  if (anyDuplicated(gene_ids[valid_genes])) {
    stop("Single-cell gene IDs are duplicated after cleaning. Aggregate the reference first.")
  }

  counts_sel <- counts[valid_genes, selected_cells, drop = FALSE]
  rownames(counts_sel) <- gene_ids[valid_genes]
  sc_reference <- Matrix::t(counts_sel)

  cell_type_labels_selected <- selected_type_labels
  names(cell_type_labels_selected) <- rownames(sc_reference)
  cell_state_labels_selected <- selected_state_labels
  names(cell_state_labels_selected) <- rownames(sc_reference)

  reference_counts <- make_reference_counts(cell_type_labels, cell_state_labels, selected_idx)

  rm(seu, counts, counts_sel)
  gc(verbose = FALSE)

  list(
    reference = sc_reference,
    cell_type_labels = cell_type_labels_selected,
    cell_state_labels = cell_state_labels_selected,
    reference_counts = reference_counts
  )
}

align_inputs <- function(sc, bulk, cfg) {
  common_genes <- intersect(colnames(sc$reference), colnames(bulk$mixture))
  if (length(common_genes) < cfg$min_common_genes) {
    stop("Too few common genes between sc reference and bulk: ", length(common_genes))
  }

  sc$reference <- sc$reference[, common_genes, drop = FALSE]
  bulk$mixture <- bulk$mixture[, common_genes, drop = FALSE]

  keep_genes <- Matrix::colSums(sc$reference) > 0 & colSums(bulk$mixture) > 0
  sc$reference <- sc$reference[, keep_genes, drop = FALSE]
  bulk$mixture <- bulk$mixture[, keep_genes, drop = FALSE]

  if (ncol(sc$reference) < cfg$min_common_genes) {
    stop("Too few nonzero common genes after filtering: ", ncol(sc$reference))
  }

  if (cfg$force_dense_reference) {
    message("Converting single-cell reference to dense matrix.")
    sc$reference <- as.matrix(sc$reference)
  }
  if (cfg$force_dense_mixture) {
    bulk$mixture <- as.matrix(bulk$mixture)
  }

  list(sc = sc, bulk = bulk)
}

write_input_summaries <- function(cfg, sc, bulk) {
  dir.create(cfg$output_dir, showWarnings = FALSE, recursive = TRUE)

  reference_type_counts <- aggregate(
    sc$reference_counts[, c("n_input", "n_selected")],
    by = list(cell_type = sc$reference_counts$cell_type),
    FUN = sum
  )
  reference_state_counts <- aggregate(
    sc$reference_counts[, c("n_input", "n_selected")],
    by = list(cell_type = sc$reference_counts$cell_state),
    FUN = sum
  )

  write.csv(reference_type_counts, file.path(cfg$output_dir, "reference_cell_type_level1_counts.csv"), row.names = FALSE)
  write.csv(sc$reference_counts, file.path(cfg$output_dir, "reference_cell_state_level3_counts.csv"), row.names = FALSE)
  write.csv(reference_state_counts, file.path(cfg$output_dir, "reference_cell_type_level3_counts.csv"), row.names = FALSE)
  write.csv(
    bulk$sample_info,
    file.path(cfg$output_dir, "bulk_sample_groups.csv"),
    row.names = FALSE
  )
  writeLines(
    colnames(sc$reference),
    file.path(cfg$output_dir, "common_genes_used.txt")
  )

  manifest <- c(
    paste("sc_rds:", cfg$sc_rds),
    paste("bulk_xlsx:", cfg$bulk_xlsx),
    paste("output_dir:", cfg$output_dir),
    paste("cell_type_label_col:", cfg$cell_type_label_col),
    paste("cell_state_label_col:", cfg$cell_state_label_col),
    paste("sc_gene_id_col:", cfg$sc_gene_id_col),
    paste("max_cells_per_type:", cfg$max_cells_per_type),
    paste("n_reference_cells:", nrow(sc$reference)),
    paste("n_bulk_samples:", nrow(bulk$mixture)),
    paste("n_common_genes_used:", ncol(sc$reference)),
    paste("bulk_groups:", paste(levels(bulk$sample_info$group), collapse = ", ")),
    paste("round_bulk_counts:", cfg$round_bulk_counts),
    paste("force_dense_reference:", cfg$force_dense_reference),
    paste("force_dense_mixture:", cfg$force_dense_mixture),
    paste("exclude_author_qc_samples:", cfg$exclude_author_qc_samples),
    paste("dry_run:", cfg$dry_run)
  )
  writeLines(manifest, file.path(cfg$output_dir, "input_manifest.txt"))
}

run_bayesprism <- function(cfg, sc, bulk) {
  check_packages("BayesPrism")

  new_prism <- getExportedValue("BayesPrism", "new.prism")
  run_prism <- getExportedValue("BayesPrism", "run.prism")

  prism_arg_aliases <- list(
    reference = list(reference = sc$reference, reference.counts = sc$reference),
    mixture = list(mixture = bulk$mixture, mixture.input = bulk$mixture),
    input_type = list(input.type = cfg$input_type, input_type = cfg$input_type),
    cell_type = list(cell.type.labels = sc$cell_type_labels, cell_type.labels = sc$cell_type_labels),
    cell_state = list(cell.state.labels = sc$cell_state_labels, cell_state.labels = sc$cell_state_labels),
    key = list(key = cfg$key),
    outlier_cut = list(outlier.cut = cfg$outlier_cut, outlier_cut = cfg$outlier_cut),
    outlier_fraction = list(outlier.fraction = cfg$outlier_fraction, outlier_fraction = cfg$outlier_fraction)
  )

  new_formals <- names(formals(new_prism))
  accepts_dots <- "..." %in% new_formals
  prism_args <- list()

  for (arg_group in prism_arg_aliases) {
    matched_name <- intersect(names(arg_group), new_formals)
    if (length(matched_name) > 0L) {
      prism_args[[matched_name[1]]] <- arg_group[[matched_name[1]]]
    } else if (accepts_dots) {
      prism_args[[names(arg_group)[1]]] <- arg_group[[1]]
    }
  }

  message("Creating BayesPrism object.")
  prism <- do.call(new_prism, prism_args)

  message("Running BayesPrism with n.cores = ", cfg$n_cores)
  run_formals <- names(formals(run_prism))
  run_args <- list()
  if ("prism" %in% run_formals || "..." %in% run_formals) {
    run_args$prism <- prism
  } else {
    run_args[[run_formals[1]]] <- prism
  }
  if ("n.cores" %in% run_formals || "..." %in% run_formals) {
    run_args[["n.cores"]] <- cfg$n_cores
  } else if ("n_cores" %in% run_formals) {
    run_args$n_cores <- cfg$n_cores
  } else if ("n.core" %in% run_formals) {
    run_args[["n.core"]] <- cfg$n_cores
  }
  bp_res <- do.call(run_prism, run_args)

  saveRDS(bp_res, file.path(cfg$output_dir, "bayesprism_result.rds"))
  bp_res
}

try_get_fraction <- function(bp_res, which_theta, state_or_type) {
  get_fraction <- tryCatch(
    getExportedValue("BayesPrism", "get.fraction"),
    error = function(e) NULL
  )

  if (is.null(get_fraction)) {
    return(NULL)
  }

  first_arg_names <- c("bp", "object", "prism")
  for (first_arg in first_arg_names) {
    args <- list(bp_res, which.theta = which_theta, state.or.type = state_or_type)
    names(args)[1] <- first_arg
    theta <- tryCatch(do.call(get_fraction, args), error = function(e) NULL)
    if (!is.null(theta)) {
      return(theta)
    }
  }

  NULL
}

coerce_theta_matrix <- function(theta, sample_ids) {
  if (is.null(theta)) {
    return(NULL)
  }
  if (is.list(theta) && !is.data.frame(theta)) {
    if ("theta" %in% names(theta)) {
      theta <- theta$theta
    } else {
      theta <- theta[[1]]
    }
  }

  theta <- as.matrix(theta)
  row_overlap <- length(intersect(rownames(theta), sample_ids))
  col_overlap <- length(intersect(colnames(theta), sample_ids))
  if (col_overlap > row_overlap) {
    theta <- t(theta)
  }
  theta
}

clr_transform <- function(theta, pseudocount) {
  theta <- as.matrix(theta)
  theta <- theta + pseudocount
  log_theta <- log(theta)
  sweep(log_theta, 1L, rowMeans(log_theta), FUN = "-")
}

theta_to_long <- function(theta, sample_info, value_col) {
  theta_df <- as.data.frame(as.table(theta), stringsAsFactors = FALSE)
  colnames(theta_df) <- c("sample", "cell_type", value_col)
  theta_df <- merge(theta_df, sample_info, by = "sample", all.x = TRUE, sort = FALSE)
  theta_df
}

filter_author_qc_samples <- function(theta, sample_info, cfg) {
  if (!cfg$exclude_author_qc_samples) {
    return(list(theta = theta, sample_info = sample_info, removed_samples = character(0)))
  }

  removed_samples <- intersect(author_qc_samples(), rownames(theta))
  keep_samples <- setdiff(rownames(theta), removed_samples)
  theta <- theta[keep_samples, , drop = FALSE]
  sample_info <- sample_info[sample_info$sample %in% keep_samples, , drop = FALSE]
  sample_info <- sample_info[match(keep_samples, sample_info$sample), , drop = FALSE]
  sample_info$group <- factor(sample_info$group, levels = levels(sample_info$group))

  list(theta = theta, sample_info = sample_info, removed_samples = removed_samples)
}

order_cell_types_by_mean <- function(theta_long, value_col) {
  means <- aggregate(theta_long[[value_col]], by = list(cell_type = theta_long$cell_type), FUN = mean, na.rm = TRUE)
  means$cell_type[order(means$x, means$cell_type)]
}

compute_group_stats <- function(theta_long, value_col) {
  groups <- levels(theta_long$group)
  groups <- groups[groups %in% unique(as.character(theta_long$group))]
  if (length(groups) < 2L) {
    return(data.frame())
  }

  group_a <- groups[1]
  group_b <- groups[2]
  cell_types <- sort(unique(theta_long$cell_type))

  stats <- lapply(cell_types, function(cell_type) {
    dat <- theta_long[theta_long$cell_type == cell_type, , drop = FALSE]
    x <- dat[as.character(dat$group) == group_a, value_col]
    y <- dat[as.character(dat$group) == group_b, value_col]
    p_value <- NA_real_
    if (length(x) > 0L && length(y) > 0L) {
      p_value <- tryCatch(
        stats::t.test(x, y)$p.value,
        error = function(e) NA_real_
      )
    }
    data.frame(
      cell_type = cell_type,
      group_a = group_a,
      group_b = group_b,
      n_group_a = sum(!is.na(x)),
      n_group_b = sum(!is.na(y)),
      mean_group_a = mean(x, na.rm = TRUE),
      mean_group_b = mean(y, na.rm = TRUE),
      delta_mean_b_minus_a = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
      median_group_a = stats::median(x, na.rm = TRUE),
      median_group_b = stats::median(y, na.rm = TRUE),
      delta_median_b_minus_a = stats::median(y, na.rm = TRUE) - stats::median(x, na.rm = TRUE),
      t_p = p_value,
      stringsAsFactors = FALSE
    )
  })

  stats_df <- do.call(rbind, stats)
  stats_df$t_fdr <- stats::p.adjust(stats_df$t_p, method = "BH")
  stats_df$p_label <- ifelse(
    is.na(stats_df$t_p),
    "NA",
    ifelse(stats_df$t_fdr < 0.05, paste0("p = ", format(stats_df$t_p, digits = 2)), "ns")
  )
  stats_df[order(stats_df$t_fdr, stats_df$t_p), , drop = FALSE]
}

make_group_palette <- function(groups) {
  groups <- as.character(groups)
  cols <- grDevices::hcl.colors(length(groups), palette = "Dark 3")
  names(cols) <- groups
  if ("Early Lesion" %in% groups) cols["Early Lesion"] <- "#0072B2"
  if ("Late Lesion" %in% groups) cols["Late Lesion"] <- "#D55E00"
  cols
}

save_boxplot <- function(
  theta_long,
  value_col,
  out_prefix,
  y_label,
  percent_axis = FALSE,
  p_values = NULL,
  add_p_labels = FALSE
) {
  suppressPackageStartupMessages(library(ggplot2))

  theta_long$group <- factor(theta_long$group, levels = levels(theta_long$group))
  groups <- levels(theta_long$group)
  groups <- groups[groups %in% unique(as.character(theta_long$group))]
  palette <- make_group_palette(groups)
  label_position <- max(theta_long[[value_col]], na.rm = TRUE)
  label_position <- label_position + max(abs(label_position), 1) * 0.25

  p <- ggplot2::ggplot(
    theta_long,
    ggplot2::aes(
      x = cell_type,
      y = .data[[value_col]],
      fill = group
    )
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.25) +
    ggplot2::geom_point(
      ggplot2::aes(color = group),
      position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.65),
      size = 0.45,
      alpha = 0.45,
      stroke = 0
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::scale_color_manual(values = palette, drop = FALSE) +
    ggplot2::labs(x = NULL, y = y_label, fill = NULL, color = NULL) +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 6),
      axis.text.x = ggplot2::element_text(size = 7),
      axis.title.x = ggplot2::element_text(size = 8, face = "bold"),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 7),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )

  if (percent_axis) {
    p <- p + ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%"))
  }
  if (add_p_labels && !is.null(p_values) && nrow(p_values) > 0L) {
    p_values$cell_type <- factor(p_values$cell_type, levels = levels(theta_long$cell_type))
    p_values$label_position <- label_position
    p <- p +
      ggplot2::geom_text(
        data = p_values,
        ggplot2::aes(x = cell_type, y = label_position, label = p_label),
        inherit.aes = FALSE,
        hjust = 1,
        vjust = 0,
        size = 2.1
      ) +
      ggplot2::expand_limits(y = label_position)
  }

  ggplot2::ggsave(paste0(out_prefix, ".pdf"), p, width = 7.0, height = 7.5, units = "in")
  ggplot2::ggsave(paste0(out_prefix, ".png"), p, width = 7.0, height = 7.5, units = "in", dpi = 600)
}

write_theta_outputs <- function(cfg, theta, sample_info, tag, label_tag) {
  sample_ids <- sample_info$sample
  theta <- coerce_theta_matrix(theta, sample_ids)
  if (is.null(theta)) {
    return(FALSE)
  }

  rownames(theta) <- as.character(rownames(theta))
  theta <- theta[sample_ids[sample_ids %in% rownames(theta)], , drop = FALSE]

  theta_all_path <- file.path(cfg$output_dir, paste0("theta_", tag, "_", label_tag, "_all_samples.csv"))
  write.csv(
    data.frame(sample = rownames(theta), theta, check.names = FALSE),
    theta_all_path,
    row.names = FALSE
  )

  filtered <- filter_author_qc_samples(theta, sample_info, cfg)
  theta <- filtered$theta
  sample_info <- filtered$sample_info

  if (length(filtered$removed_samples) > 0L) {
    writeLines(
      filtered$removed_samples,
      file.path(cfg$output_dir, paste0("theta_", tag, "_", label_tag, "_excluded_author_qc_samples.txt"))
    )
  }

  theta_path <- file.path(cfg$output_dir, paste0("theta_", tag, "_", label_tag, ".csv"))
  write.csv(
    data.frame(sample = rownames(theta), theta, check.names = FALSE),
    theta_path,
    row.names = FALSE
  )

  theta_long <- theta_to_long(theta, sample_info, "fraction")
  cell_type_levels <- order_cell_types_by_mean(theta_long, "fraction")
  theta_long$cell_type <- factor(theta_long$cell_type, levels = cell_type_levels)
  write.csv(
    theta_long,
    file.path(cfg$output_dir, paste0("theta_", tag, "_", label_tag, "_long.csv")),
    row.names = FALSE
  )

  fraction_stats <- compute_group_stats(theta_long, "fraction")
  write.csv(
    fraction_stats,
    file.path(cfg$output_dir, paste0("theta_", tag, "_", label_tag, "_group_stats.csv")),
    row.names = FALSE
  )

  save_boxplot(
    theta_long,
    "fraction",
    file.path(cfg$output_dir, paste0("theta_", tag, "_", label_tag, "_fraction_boxplot")),
    "BayesPrism fraction",
    percent_axis = TRUE
  )

  theta_clr <- clr_transform(theta, cfg$clr_pseudocount)
  write.csv(
    data.frame(sample = rownames(theta_clr), theta_clr, check.names = FALSE),
    file.path(cfg$output_dir, paste0("theta_", tag, "_", label_tag, "_clr.csv")),
    row.names = FALSE
  )

  theta_clr_long <- theta_to_long(theta_clr, sample_info, "clr")
  theta_clr_long$cell_type <- factor(theta_clr_long$cell_type, levels = cell_type_levels)
  write.csv(
    theta_clr_long,
    file.path(cfg$output_dir, paste0("theta_", tag, "_", label_tag, "_clr_long.csv")),
    row.names = FALSE
  )

  clr_stats <- compute_group_stats(theta_clr_long, "clr")
  write.csv(
    clr_stats,
    file.path(cfg$output_dir, paste0("theta_", tag, "_", label_tag, "_clr_group_stats.csv")),
    row.names = FALSE
  )

  save_boxplot(
    theta_clr_long,
    "clr",
    file.path(cfg$output_dir, paste0("theta_", tag, "_", label_tag, "_clr_boxplot")),
    "CLR fraction",
    p_values = clr_stats,
    add_p_labels = TRUE
  )

  TRUE
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(cfg$output_dir, showWarnings = FALSE, recursive = TRUE)

  check_packages(c("Seurat", "Matrix", "readxl", "ggplot2"))
  if (!cfg$dry_run) {
    check_packages("BayesPrism")
  }

  bulk <- read_bulk_matrix(cfg)
  sc <- prepare_sc_reference(cfg)
  aligned <- align_inputs(sc, bulk, cfg)
  sc <- aligned$sc
  bulk <- aligned$bulk

  write_input_summaries(cfg, sc, bulk)

  message("Reference cells: ", nrow(sc$reference))
  message("Bulk samples: ", nrow(bulk$mixture))
  message("Common nonzero genes: ", ncol(sc$reference))
  print(table(bulk$sample_info$group))

  if (cfg$dry_run) {
    message("Dry run complete. BayesPrism was not executed.")
    return(invisible(NULL))
  }

  bp_res <- run_bayesprism(cfg, sc, bulk)
  sample_info <- bulk$sample_info

  theta_final_type <- try_get_fraction(bp_res, which_theta = "final", state_or_type = "type")
  wrote_final_type <- write_theta_outputs(cfg, theta_final_type, sample_info, "final", "cell_type_level1")

  theta_final_state <- try_get_fraction(bp_res, which_theta = "final", state_or_type = "state")
  wrote_final_state <- write_theta_outputs(cfg, theta_final_state, sample_info, "final", "cell_state_level3")

  theta_initial_type <- try_get_fraction(bp_res, which_theta = "first.gibbs", state_or_type = "type")
  wrote_initial_type <- write_theta_outputs(cfg, theta_initial_type, sample_info, "first_gibbs", "cell_type_level1")

  theta_initial_state <- try_get_fraction(bp_res, which_theta = "first.gibbs", state_or_type = "state")
  wrote_initial_state <- write_theta_outputs(cfg, theta_initial_state, sample_info, "first_gibbs", "cell_state_level3")

  if (!wrote_final_type && !wrote_final_state && !wrote_initial_type && !wrote_initial_state) {
    stop("BayesPrism finished, but theta fractions could not be extracted with get.fraction().")
  }

  message("BayesPrism output saved to: ", cfg$output_dir)
}

main()
