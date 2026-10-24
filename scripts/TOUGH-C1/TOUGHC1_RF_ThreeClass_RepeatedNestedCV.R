# =====================================================================
# FUSED-RF on the exact DeepDrug3D three-class TOUGH-C1 outer folds
#   - 5 repeated stratified 5-fold outer CV (25 evaluations)
#   - One three-class model: CONTROL, HEME, NUC
#   - Inner 5-fold joint selection of t_end and RF min.node.size
#   - t_end grid: 6, 7, ..., 20 A; raw descriptor grid remains 0.1 A
#   - Fold-specific MFPCA retaining at least 95% of training variation
#   - NUC and HEME reported using one-versus-rest metrics
# =====================================================================

rm(list = ls())
invisible(gc())
set.seed(1)

suppressPackageStartupMessages({
  library(fda)
  library(parallel)
  library(ranger)
})

descriptor_file <- Sys.getenv("FUSED_DESCRIPTOR_FILE", unset = file.path("data", "cache", "TOUGH_C1_LBS_ligand_ILR_CDPA_to20A.rds"))
fold_file <- Sys.getenv("FUSED_FOLD_FILE", unset = file.path(
  "scripts", "TOUGH-C1", "folds",
  "TOUGHC1_threeclass_5x5_outer_fold_assignments.csv"
))
allow_subset <- identical(Sys.getenv("FUSED_ALLOW_SUBSET"), "1")
cohort_name <- if (allow_subset) "threeclass_common4000" else "threeclass_full4002"
out_dir <- Sys.getenv("FUSED_OUTPUT_DIR", unset = file.path(
  "results", "TOUGH-C1", cohort_name, "RF"
))
checkpoint_dir <- file.path(out_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

T_END_GRID <- seq(6, 20, by = 1)
K_INNER <- 5L
VAR_TARGET <- 0.95
MAX_PC <- 20L
SMOOTH_LAMBDA <- 1e-3
NORDER <- 4L
BASIS_DIVS <- c(2, 3, 4, 5)
RF_NUM_TREES <- 400L
RF_MIN_NODE_SIZE_GRID <- c(2L, 5L, 7L)
N_CORES <- as.integer(Sys.getenv("FUSED_N_CORES", unset = "4"))
if (is.na(N_CORES) || N_CORES < 1L) N_CORES <- 1L

required_files <- c(descriptor_file, fold_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required file(s): ", paste(missing_files, collapse = ", "))
}

obj <- readRDS(descriptor_file)
thr_full <- as.numeric(obj$thr_full)
ilr_arr <- obj$ilr_arr
cov_arr <- obj$cov_arr
labels <- factor(obj$labels, levels = c("CONTROL", "HEME", "NUC"))
t_start <- as.numeric(obj$d_min_sel)
folds <- read.csv(fold_file, stringsAsFactors = FALSE)

retained_rows <- sort(unique(folds$descriptor_row_index))
if (length(labels) != 4002L || nrow(folds) != 5L * length(retained_rows) ||
    any(!retained_rows %in% seq_along(labels)) ||
    anyDuplicated(folds[c("repetition", "descriptor_row_index")])) {
  stop("Unexpected descriptor or fold-assignment dimensions.")
}
if (!allow_subset && !identical(retained_rows, seq_along(labels))) {
  stop("Fold descriptor_row_index values do not match the descriptor rows.")
}
fold_label <- folds$class[match(retained_rows, folds$descriptor_row_index)]
if (!identical(as.character(labels[retained_rows]), fold_label)) {
  stop("Class labels in the fold file do not match the descriptor object.")
}
if (!identical(sort(unique(folds$repetition)), 1:5) ||
    !identical(sort(unique(folds$outer_test_fold)), 1:5)) stop("Expected 5 repeats and 5 folds.")
for (repeat_id in 1:5) {
  if (!identical(sort(folds$descriptor_row_index[folds$repetition == repeat_id]), retained_rows)) {
    stop("Repetitions have different pocket cohorts.")
  }
}
if (is.null(obj$pairs_df) ||
    !identical(as.character(obj$pairs_df$p_base[folds$descriptor_row_index]), folds$p_base) ||
    !identical(as.character(obj$pairs_df$l_base[folds$descriptor_row_index]), folds$l_base)) {
  stop("Fold identities do not match cached protein/ligand identities.")
}
if (allow_subset) {
  if (length(retained_rows) != 4000L || is.null(obj$pairs_df)) stop("Expected common 4000-pocket cohort.")
  excluded <- obj$pairs_df$p_base[setdiff(seq_along(labels), retained_rows)]
  if (!setequal(excluded, c("1ukmB", "2o8cD"))) stop("Unexpected excluded pockets.")
  if (!identical(as.character(obj$pairs_df$p_base[folds$descriptor_row_index]), folds$p_base) ||
      !identical(as.character(obj$pairs_df$l_base[folds$descriptor_row_index]), folds$l_base)) {
    stop("Fold protein/ligand identifiers do not match descriptor rows.")
  }
}
input_signature <- paste(unname(tools::md5sum(c(descriptor_file, fold_file,
  file.path("scripts", "TOUGH-C1", "TOUGHC1_RF_ThreeClass_RepeatedNestedCV.R")))), collapse = ":")
cat("Verified", length(retained_rows), "retained pockets; training excludes all other descriptor rows.\n")
if (identical(Sys.getenv("FUSED_AUDIT_ONLY"), "1")) quit(save = "no", status = 0)
requested_repeat <- Sys.getenv("FUSED_REPETITION")
requested_fold <- Sys.getenv("FUSED_OUTER_FOLD")
repeat_ids <- if (nzchar(requested_repeat)) as.integer(requested_repeat) else 1:5
fold_ids <- if (nzchar(requested_fold)) as.integer(requested_fold) else 1:5
if (anyNA(c(repeat_ids, fold_ids)) || any(!c(repeat_ids, fold_ids) %in% 1:5)) stop("Invalid requested fold.")

make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  output <- vector("list", k)
  for (class_indices in split(seq_along(y), y)) {
    shuffled <- sample(class_indices)
    assignment <- rep(seq_len(k), length.out = length(shuffled))
    for (fold in seq_len(k)) {
      output[[fold]] <- c(output[[fold]], shuffled[assignment == fold])
    }
  }
  lapply(output, sort)
}

choose_nbasis <- function(tvec) {
  candidates <- sort(unique(pmax(NORDER, round(length(tvec) / BASIS_DIVS))))
  candidates[1]
}

smooth_channel <- function(mat, tvec, nbasis) {
  basis <- create.bspline.basis(range(tvec), nbasis, NORDER)
  fd_parameter <- fdPar(basis, int2Lfd(2), SMOOTH_LAMBDA)
  fit <- smooth.basis(tvec, t(mat), fd_parameter)
  t(eval.fd(tvec, fit$fd))
}

smooth_all_channels <- function(t_end) {
  threshold_indices <- which(thr_full >= t_start & thr_full <= t_end)
  tvec <- thr_full[threshold_indices]
  nbasis <- choose_nbasis(tvec)
  list(
    ilr1 = smooth_channel(ilr_arr[, threshold_indices, 1], tvec, nbasis),
    ilr2 = smooth_channel(ilr_arr[, threshold_indices, 2], tvec, nbasis),
    c12 = smooth_channel(cov_arr[, threshold_indices, 1], tvec, nbasis),
    c13 = smooth_channel(cov_arr[, threshold_indices, 2], tvec, nbasis),
    c23 = smooth_channel(cov_arr[, threshold_indices, 3], tvec, nbasis),
    v1 = smooth_channel(cov_arr[, threshold_indices, 4], tvec, nbasis),
    v2 = smooth_channel(cov_arr[, threshold_indices, 5], tvec, nbasis),
    v3 = smooth_channel(cov_arr[, threshold_indices, 6], tvec, nbasis)
  )
}

fit_pca_scores <- function(x_train, x_test) {
  center <- colMeans(x_train)
  train_centered <- sweep(x_train, 2, center, "-")
  test_centered <- sweep(x_test, 2, center, "-")
  eig <- eigen(cov(train_centered), symmetric = TRUE)
  keep <- which(eig$values > .Machine$double.eps)
  list(
    train = train_centered %*% eig$vectors[, keep, drop = FALSE],
    test = test_centered %*% eig$vectors[, keep, drop = FALSE],
    eigenvalues = pmax(eig$values[keep], 0)
  )
}

build_mfpca_scores <- function(channels, train_indices, test_indices) {
  ilr_energy <- mean(abs(c(
    channels$ilr1[train_indices, ], channels$ilr2[train_indices, ]
  )))
  cdpa_energy <- mean(abs(c(
    channels$c12[train_indices, ], channels$c13[train_indices, ],
    channels$c23[train_indices, ], channels$v1[train_indices, ],
    channels$v2[train_indices, ], channels$v3[train_indices, ]
  )))
  ilr_scale <- if (ilr_energy == 0) 1 else cdpa_energy / ilr_energy

  make_features <- function(indices) {
    cbind(
      channels$ilr1[indices, ] * ilr_scale,
      channels$ilr2[indices, ] * ilr_scale,
      channels$c12[indices, ], channels$c13[indices, ],
      channels$c23[indices, ], channels$v1[indices, ],
      channels$v2[indices, ], channels$v3[indices, ]
    )
  }

  pca <- fit_pca_scores(make_features(train_indices), make_features(test_indices))
  cumulative <- cumsum(pca$eigenvalues / sum(pca$eigenvalues))
  k_use <- which(cumulative >= VAR_TARGET)[1]
  if (is.na(k_use)) k_use <- length(cumulative)
  k_use <- min(k_use, MAX_PC)
  list(
    train = pca$train[, seq_len(k_use), drop = FALSE],
    test = pca$test[, seq_len(k_use), drop = FALSE],
    k_use = k_use
  )
}

fit_rf <- function(x_train, y_train, min_node_size, seed) {
  ranger(
    dependent.variable.name = "y",
    data = data.frame(y = y_train, x_train, check.names = FALSE),
    probability = TRUE,
    num.trees = RF_NUM_TREES,
    mtry = max(1L, floor(sqrt(ncol(x_train)))),
    min.node.size = min_node_size,
    oob.error = FALSE,
    num.threads = 1L,
    seed = seed
  )
}

predict_rf <- function(model, x_test) {
  probabilities <- predict(model, data.frame(x_test, check.names = FALSE))$predictions
  probabilities <- probabilities[, levels(labels), drop = FALSE]
  predictions <- factor(
    levels(labels)[max.col(probabilities, ties.method = "first")],
    levels = levels(labels)
  )
  list(predictions = predictions, probabilities = probabilities)
}

auc_rank <- function(is_positive, scores) {
  n_positive <- sum(is_positive)
  n_negative <- sum(!is_positive)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[is_positive]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * n_negative)
}

one_vs_rest_metrics <- function(y_true, predicted, probabilities, positive_label) {
  actual_positive <- y_true == positive_label
  predicted_positive <- predicted == positive_label
  tp <- sum(actual_positive & predicted_positive)
  fp <- sum(!actual_positive & predicted_positive)
  tn <- sum(!actual_positive & !predicted_positive)
  fn <- sum(actual_positive & !predicted_positive)
  safe_divide <- function(x, y) if (y == 0) NA_real_ else x / y
  data.frame(
    comparison = paste0(positive_label, "_vs_REST"),
    accuracy = safe_divide(tp + tn, tp + fp + tn + fn),
    ppv = safe_divide(tp, tp + fp),
    tpr = safe_divide(tp, tp + fn),
    tnr = safe_divide(tn, tn + fp),
    auc = auc_rank(actual_positive, probabilities[, positive_label]),
    tp = tp, fp = fp, tn = tn, fn = fn
  )
}

evaluate_t_end <- function(t_end, outer_train, inner_folds_relative, repetition, outer_fold) {
  channels <- smooth_all_channels(t_end)
  records <- list()
  for (inner_fold in seq_len(K_INNER)) {
    inner_validation <- outer_train[inner_folds_relative[[inner_fold]]]
    inner_training <- setdiff(outer_train, inner_validation)
    scores <- build_mfpca_scores(channels, inner_training, inner_validation)
    for (min_node_size in RF_MIN_NODE_SIZE_GRID) {
      seed <- 300000 + 10000 * repetition + 1000 * outer_fold +
        10 * inner_fold + min_node_size
      model <- fit_rf(
        scores$train, droplevels(labels[inner_training]), min_node_size, seed
      )
      prediction <- predict_rf(model, scores$test)
      records[[length(records) + 1L]] <- data.frame(
        repetition = repetition,
        outer_fold = outer_fold,
        inner_fold = inner_fold,
        t_end = t_end,
        min_node_size = min_node_size,
        retained_mfpcs = scores$k_use,
        accuracy = mean(prediction$predictions == labels[inner_validation])
      )
    }
  }
  do.call(rbind, records)
}

outer_records <- list()
inner_records <- list()
prediction_records <- list()
start_time <- Sys.time()

for (repetition in repeat_ids) {
  repetition_folds <- folds[folds$repetition == repetition, ]
  for (outer_fold in fold_ids) {
    checkpoint_file <- file.path(
      checkpoint_dir,
      sprintf("FUSED_threeclass_repeat%02d_outer%02d.rds", repetition, outer_fold)
    )
    if (file.exists(checkpoint_file)) {
      cat(sprintf("Loading completed repetition %d, outer fold %d.\n", repetition, outer_fold))
      completed <- readRDS(checkpoint_file)
      if (!identical(completed$input_signature, input_signature)) {
        stop("Checkpoint does not match common-cohort inputs.")
      }
      outer_records[[length(outer_records) + 1L]] <- completed$outer
      inner_records[[length(inner_records) + 1L]] <- completed$inner
      prediction_records[[length(prediction_records) + 1L]] <- completed$predictions
      next
    }

    cat(sprintf("\n===== Repetition %d, outer fold %d =====\n", repetition, outer_fold))
    outer_test <- repetition_folds$descriptor_row_index[
      repetition_folds$outer_test_fold == outer_fold
    ]
    outer_train <- setdiff(retained_rows, outer_test)
    inner_folds_relative <- make_stratified_folds(
      droplevels(labels[outer_train]), K_INNER,
      seed = 200000 + 100 * repetition + outer_fold
    )

    candidate_results <- mclapply(
      T_END_GRID,
      evaluate_t_end,
      outer_train = outer_train,
      inner_folds_relative = inner_folds_relative,
      repetition = repetition,
      outer_fold = outer_fold,
      mc.cores = min(N_CORES, length(T_END_GRID)),
      mc.preschedule = TRUE
    )
    inner_table <- do.call(rbind, candidate_results)
    aggregate_table <- aggregate(
      accuracy ~ t_end + min_node_size,
      data = inner_table,
      FUN = mean
    )
    aggregate_table <- aggregate_table[
      order(-aggregate_table$accuracy, aggregate_table$t_end, aggregate_table$min_node_size),
    ]
    selected <- aggregate_table[1, ]

    final_channels <- smooth_all_channels(selected$t_end)
    final_scores <- build_mfpca_scores(final_channels, outer_train, outer_test)
    final_model <- fit_rf(
      final_scores$train,
      droplevels(labels[outer_train]),
      selected$min_node_size,
      seed = 900000 + 100 * repetition + outer_fold
    )
    final_prediction <- predict_rf(final_model, final_scores$test)
    multiclass_accuracy <- mean(final_prediction$predictions == labels[outer_test])

    outer_table <- do.call(rbind, lapply(c("HEME", "NUC"), function(positive_label) {
      metrics <- one_vs_rest_metrics(
        labels[outer_test], final_prediction$predictions,
        final_prediction$probabilities, positive_label
      )
      metrics$repetition <- repetition
      metrics$outer_fold <- outer_fold
      metrics$multiclass_accuracy <- multiclass_accuracy
      metrics$selected_t_end <- selected$t_end
      metrics$selected_min_node_size <- selected$min_node_size
      metrics$retained_mfpcs <- final_scores$k_use
      metrics$mean_inner_accuracy <- selected$accuracy
      metrics
    }))

    prediction_table <- data.frame(
      repetition = repetition,
      outer_fold = outer_fold,
      descriptor_row_index = outer_test,
      true_class = as.character(labels[outer_test]),
      predicted_class = as.character(final_prediction$predictions),
      probability_CONTROL = final_prediction$probabilities[, "CONTROL"],
      probability_HEME = final_prediction$probabilities[, "HEME"],
      probability_NUC = final_prediction$probabilities[, "NUC"]
    )
    if (allow_subset) {
      write.csv(data.frame(descriptor_row_index = outer_train,
                          p_base = obj$pairs_df$p_base[outer_train],
                          class = as.character(labels[outer_train])),
                file.path(out_dir, "training_pockets.csv"), row.names = FALSE)
    }

    completed <- list(
      input_signature = input_signature,
      training_descriptor_rows = outer_train,
      test_descriptor_rows = outer_test,
      outer = outer_table,
      inner = inner_table,
      predictions = prediction_table
    )
    saveRDS(completed, checkpoint_file)
    outer_records[[length(outer_records) + 1L]] <- outer_table
    inner_records[[length(inner_records) + 1L]] <- inner_table
    prediction_records[[length(prediction_records) + 1L]] <- prediction_table

    cat(sprintf(
      "Selected t_end=%.1f, min.node.size=%d, MFPCs=%d, multiclass ACC=%.4f\n",
      selected$t_end, selected$min_node_size, final_scores$k_use, multiclass_accuracy
    ))
    rm(candidate_results, final_channels, final_scores, final_model)
    invisible(gc())
  }
}

outer_table <- do.call(rbind, outer_records)
inner_table <- do.call(rbind, inner_records)
prediction_table <- do.call(rbind, prediction_records)
outer_table <- outer_table[order(outer_table$comparison, outer_table$repetition, outer_table$outer_fold), ]

summary_table <- do.call(rbind, lapply(split(outer_table, outer_table$comparison), function(x) {
  data.frame(
    method = "FUSED_RF_three_class",
    comparison = x$comparison[1],
    n_outer_evaluations = nrow(x),
    mean_accuracy = mean(x$accuracy), sd_accuracy = sd(x$accuracy),
    mean_ppv = mean(x$ppv), sd_ppv = sd(x$ppv),
    mean_tpr = mean(x$tpr), sd_tpr = sd(x$tpr),
    mean_tnr = mean(x$tnr), sd_tnr = sd(x$tnr),
    mean_auc = mean(x$auc), sd_auc = sd(x$auc),
    mean_multiclass_accuracy = mean(x$multiclass_accuracy),
    sd_multiclass_accuracy = sd(x$multiclass_accuracy),
    mean_selected_t_end = mean(x$selected_t_end),
    sd_selected_t_end = sd(x$selected_t_end),
    mean_selected_min_node_size = mean(x$selected_min_node_size),
    sd_selected_min_node_size = sd(x$selected_min_node_size),
    mean_retained_mfpcs = mean(x$retained_mfpcs),
    sd_retained_mfpcs = sd(x$retained_mfpcs)
  )
}))
rownames(summary_table) <- NULL

write.csv(outer_table, file.path(out_dir, "FUSED_threeclass_outer_metrics.csv"), row.names = FALSE)
write.csv(inner_table, file.path(out_dir, "FUSED_threeclass_inner_metrics.csv"), row.names = FALSE)
write.csv(prediction_table, file.path(out_dir, "FUSED_threeclass_outer_predictions.csv"), row.names = FALSE)
write.csv(summary_table, file.path(out_dir, "FUSED_threeclass_summary.csv"), row.names = FALSE)

elapsed_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("\n===== FUSED three-class same-fold summary =====\n")
print(summary_table, row.names = FALSE)
cat(sprintf("\nElapsed script-session time: %.2f minutes\n", elapsed_minutes))
cat("Saved outputs under:", out_dir, "\n")
