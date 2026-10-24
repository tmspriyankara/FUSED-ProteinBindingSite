# =====================================================================
# TOUGH-C1 | Three-class FUSED models with one-versus-rest evaluation
#   - Exact 5 repeated 5-fold outer splits used for DeepDrug3D.
#   - One three-class model is fitted to CONTROL, HEME, and NUC.
#   - Inner 5-fold CV selects t_end and model hyperparameters using
#     three-class accuracy; the outer fold is used only for evaluation.
#   - RF results are imported from the completed exact-fold analysis.
#   - Multinomial lasso and nearest Mahalanobis mean are fitted here.
# =====================================================================

rm(list = ls()); invisible(gc())
set.seed(1)

suppressPackageStartupMessages({
  library(fda)
  library(glmnet)
  library(parallel)
})

descriptor_file <- Sys.getenv("FUSED_DESCRIPTOR_FILE", unset = file.path("data", "cache", "TOUGH_C1_LBS_ligand_ILR_CDPA_to20A.rds"))
fold_file <- Sys.getenv("FUSED_FOLD_FILE", unset = file.path(
  "scripts", "TOUGH-C1", "folds",
  "TOUGHC1_threeclass_5x5_outer_fold_assignments.csv"
))
allow_subset <- identical(Sys.getenv("FUSED_ALLOW_SUBSET"), "1")
cohort_name <- if (allow_subset) "threeclass_common4000" else "threeclass_full4002"
out_dir <- Sys.getenv("FUSED_OUTPUT_DIR", unset = file.path(
  "results", "TOUGH-C1", cohort_name, "all_classifiers"
))
skip_rf <- identical(Sys.getenv("FUSED_SKIP_RF_IMPORT"), "1")
rf_outer_file <- Sys.getenv("FUSED_RF_RESULTS", unset = file.path(
  "results", "TOUGH-C1", cohort_name, "RF", "FUSED_threeclass_outer_metrics.csv"
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
LASSO_LAMBDA_GRID <- 10^seq(-2, -4, length.out = 10)
MAH_EPS <- 1e-6
N_CORES <- as.integer(Sys.getenv("FUSED_N_CORES", unset = "4"))
if (is.na(N_CORES) || N_CORES < 1L) N_CORES <- 1L

required_files <- c(descriptor_file, fold_file, if (!skip_rf) rf_outer_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing required file(s): ", paste(missing_files, collapse = ", "))

obj <- readRDS(descriptor_file)
thr_full <- as.numeric(obj$thr_full)
ilr_arr <- obj$ilr_arr
cov_arr <- obj$cov_arr
labels <- factor(obj$labels, levels = c("CONTROL", "HEME", "NUC"))
t_start <- as.numeric(obj$d_min_sel)
folds <- read.csv(fold_file, stringsAsFactors = FALSE)

retained_rows <- sort(unique(folds$descriptor_row_index))
if (length(labels) != 4002L || nrow(folds) != 5L * length(retained_rows) ||
    any(!retained_rows %in% seq_along(labels)) || anyDuplicated(folds[c("repetition", "descriptor_row_index")])) {
  stop("Unexpected descriptor or fold-assignment dimensions.")
}
if (!allow_subset && !identical(retained_rows, seq_along(labels))) {
  stop("Fold descriptor_row_index values do not match descriptor rows.")
}
fold_label <- folds$class[match(retained_rows, folds$descriptor_row_index)]
if (!identical(as.character(labels[retained_rows]), fold_label)) {
  stop("Class labels in the fold file do not match the descriptor object.")
}
if (!identical(sort(unique(folds$repetition)), 1:5) ||
    !identical(sort(unique(folds$outer_test_fold)), 1:5)) stop("Expected five repeats and five folds.")
for (repeat_id in 1:5) {
  if (!identical(sort(folds$descriptor_row_index[folds$repetition == repeat_id]), retained_rows)) {
    stop("Repetitions have different cohorts.")
  }
}
if (is.null(obj$pairs_df) ||
    !identical(as.character(obj$pairs_df$p_base[folds$descriptor_row_index]), folds$p_base) ||
    !identical(as.character(obj$pairs_df$l_base[folds$descriptor_row_index]), folds$l_base)) {
  stop("Fold identities do not match cached protein/ligand identities.")
}
if (allow_subset && (length(retained_rows) != 4000L ||
    !setequal(obj$pairs_df$p_base[setdiff(seq_along(labels), retained_rows)], c("1ukmB", "2o8cD")))) {
  stop("Unexpected excluded pockets.")
}
input_signature <- paste(unname(tools::md5sum(c(descriptor_file, fold_file,
  "scripts/TOUGH-C1/TOUGHC1_MLR_NMM_ThreeClass_RepeatedNestedCV.R"))), collapse = ":")
cat("Verified", length(retained_rows), "retained pockets for MLR-l1 and NMM.\n")
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
  ilr_energy <- mean(abs(c(channels$ilr1[train_indices, ], channels$ilr2[train_indices, ])))
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
      channels$c12[indices, ], channels$c13[indices, ], channels$c23[indices, ],
      channels$v1[indices, ], channels$v2[indices, ], channels$v3[indices, ]
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

fit_lasso_path <- function(x_train, y_train) {
  fit_lambda_values <- sort(unique(c(1, 1e-1, LASSO_LAMBDA_GRID)), decreasing = TRUE)
  glmnet(
    x = as.matrix(x_train), y = y_train,
    family = "multinomial", alpha = 1,
    lambda = fit_lambda_values, standardize = TRUE,
    control = list(maxit = 1000000)
  )
}

lasso_grid_accuracy <- function(model, x_test, y_test) {
  prediction <- predict(
    model, newx = as.matrix(x_test), s = LASSO_LAMBDA_GRID, type = "class"
  )
  prediction <- as.matrix(prediction)
  if (ncol(prediction) != length(LASSO_LAMBDA_GRID)) {
    prediction <- matrix(as.vector(prediction), ncol = length(LASSO_LAMBDA_GRID))
  }
  vapply(seq_along(LASSO_LAMBDA_GRID), function(j) {
    mean(factor(prediction[, j], levels = levels(labels)) == y_test)
  }, numeric(1))
}

predict_lasso <- function(model, x_test, lambda_value) {
  predicted <- as.vector(predict(
    model, newx = as.matrix(x_test), s = lambda_value, type = "class"
  ))
  probability_array <- predict(
    model, newx = as.matrix(x_test), s = lambda_value, type = "response"
  )
  probabilities <- probability_array[, , 1, drop = FALSE][, , 1]
  probabilities <- probabilities[, levels(labels), drop = FALSE]
  list(
    predictions = factor(predicted, levels = levels(labels)),
    probabilities = probabilities
  )
}

fit_nmm <- function(x_train, y_train, eps_value = MAH_EPS) {
  x_train <- as.matrix(x_train)
  classes <- levels(droplevels(y_train))
  means <- lapply(classes, function(class_name) {
    colMeans(x_train[y_train == class_name, , drop = FALSE])
  })
  names(means) <- classes
  covariance <- cov(x_train) + diag(eps_value, ncol(x_train))
  inverse <- tryCatch(solve(covariance), error = function(e) MASS::ginv(covariance))
  list(classes = classes, means = means, inverse = inverse)
}

predict_nmm <- function(model, x_test) {
  x_test <- as.matrix(x_test)
  distances <- matrix(NA_real_, nrow(x_test), length(model$classes))
  colnames(distances) <- model$classes
  for (j in seq_along(model$classes)) {
    difference <- sweep(x_test, 2, model$means[[model$classes[j]]], "-")
    distances[, j] <- rowSums((difference %*% model$inverse) * difference)
  }
  predictions <- factor(
    model$classes[max.col(-distances, ties.method = "first")],
    levels = levels(labels)
  )
  list(predictions = predictions, distances = distances)
}

auc_rank <- function(is_positive, scores) {
  n_positive <- sum(is_positive)
  n_negative <- sum(!is_positive)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[is_positive]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * n_negative)
}

one_vs_rest_metrics <- function(y_true, predicted, scores, positive_label) {
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
    auc = auc_rank(actual_positive, scores),
    tp = tp, fp = fp, tn = tn, fn = fn
  )
}

evaluate_t_end <- function(t_end, outer_train, inner_folds_relative, repetition, outer_fold) {
  channels <- smooth_all_channels(t_end)
  lasso_records <- list()
  nmm_records <- list()

  for (inner_fold in seq_len(K_INNER)) {
    inner_validation <- outer_train[inner_folds_relative[[inner_fold]]]
    inner_training <- setdiff(outer_train, inner_validation)
    scores <- build_mfpca_scores(channels, inner_training, inner_validation)
    y_training <- droplevels(labels[inner_training])
    y_validation <- factor(labels[inner_validation], levels = levels(labels))

    lasso_model <- fit_lasso_path(scores$train, y_training)
    lasso_accuracy <- lasso_grid_accuracy(lasso_model, scores$test, y_validation)
    lasso_records[[inner_fold]] <- data.frame(
      model = "MLR_l1", repetition = repetition, outer_fold = outer_fold,
      inner_fold = inner_fold, t_end = t_end,
      lambda = LASSO_LAMBDA_GRID,
      eps = NA_real_,
      retained_mfpcs = scores$k_use,
      accuracy = lasso_accuracy
    )

    nmm_model <- fit_nmm(scores$train, y_training)
    nmm_prediction <- predict_nmm(nmm_model, scores$test)
    nmm_records[[inner_fold]] <- data.frame(
      model = "NMM", repetition = repetition, outer_fold = outer_fold,
      inner_fold = inner_fold, t_end = t_end,
      lambda = NA_real_, eps = MAH_EPS, retained_mfpcs = scores$k_use,
      accuracy = mean(nmm_prediction$predictions == y_validation)
    )
  }

  list(lasso = do.call(rbind, lasso_records), nmm = do.call(rbind, nmm_records))
}

select_best <- function(inner_table, hyperparameters) {
  formula <- as.formula(paste("accuracy ~", paste(hyperparameters, collapse = " + ")))
  aggregate_table <- aggregate(formula, data = inner_table, FUN = mean)
  ordering <- do.call(order, c(list(-aggregate_table$accuracy), aggregate_table[hyperparameters]))
  aggregate_table[ordering[1], , drop = FALSE]
}

build_outer_rows <- function(model_name, prediction, scores_by_class, y_test,
                             repetition, outer_fold, selected_t_end,
                             retained_mfpcs, mean_inner_accuracy,
                             selected_lambda = NA_real_) {
  multiclass_accuracy <- mean(prediction == y_test)
  do.call(rbind, lapply(c("HEME", "NUC"), function(positive_label) {
    metrics <- one_vs_rest_metrics(
      y_test, prediction, scores_by_class[, positive_label], positive_label
    )
    metrics$model <- model_name
    metrics$repetition <- repetition
    metrics$outer_fold <- outer_fold
    metrics$multiclass_accuracy <- multiclass_accuracy
    metrics$selected_t_end <- selected_t_end
    metrics$selected_lambda <- selected_lambda
    metrics$retained_mfpcs <- retained_mfpcs
    metrics$mean_inner_accuracy <- mean_inner_accuracy
    metrics
  }))
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
      sprintf("FUSED_MLR_NMM_repeat%02d_outer%02d.rds", repetition, outer_fold)
    )
    if (file.exists(checkpoint_file)) {
      cat(sprintf("Loading repetition %d, outer fold %d.\n", repetition, outer_fold))
      completed <- readRDS(checkpoint_file)
      if (!identical(completed$input_signature, input_signature)) stop("Checkpoint input/source mismatch.")
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
      T_END_GRID, evaluate_t_end,
      outer_train = outer_train,
      inner_folds_relative = inner_folds_relative,
      repetition = repetition,
      outer_fold = outer_fold,
      mc.cores = min(N_CORES, length(T_END_GRID)),
      mc.preschedule = TRUE
    )
    lasso_inner <- do.call(rbind, lapply(candidate_results, `[[`, "lasso"))
    nmm_inner <- do.call(rbind, lapply(candidate_results, `[[`, "nmm"))
    lasso_selected <- select_best(lasso_inner, c("t_end", "lambda"))
    nmm_selected <- select_best(nmm_inner, "t_end")

    y_train <- droplevels(labels[outer_train])
    y_test <- factor(labels[outer_test], levels = levels(labels))

    lasso_channels <- smooth_all_channels(lasso_selected$t_end)
    lasso_scores <- build_mfpca_scores(lasso_channels, outer_train, outer_test)
    lasso_model <- fit_lasso_path(lasso_scores$train, y_train)
    lasso_prediction <- predict_lasso(
      lasso_model, lasso_scores$test, lasso_selected$lambda
    )
    lasso_outer <- build_outer_rows(
      "MLR_l1", lasso_prediction$predictions, lasso_prediction$probabilities,
      y_test, repetition, outer_fold, lasso_selected$t_end,
      lasso_scores$k_use, lasso_selected$accuracy, lasso_selected$lambda
    )

    nmm_channels <- smooth_all_channels(nmm_selected$t_end)
    nmm_scores <- build_mfpca_scores(nmm_channels, outer_train, outer_test)
    nmm_model <- fit_nmm(nmm_scores$train, y_train)
    nmm_prediction <- predict_nmm(nmm_model, nmm_scores$test)
    nmm_outer <- build_outer_rows(
      "NMM", nmm_prediction$predictions, -nmm_prediction$distances,
      y_test, repetition, outer_fold, nmm_selected$t_end,
      nmm_scores$k_use, nmm_selected$accuracy
    )

    prediction_table <- rbind(
      data.frame(
        model = "MLR_l1", repetition = repetition, outer_fold = outer_fold,
        descriptor_row_index = outer_test,
        true_class = as.character(y_test),
        predicted_class = as.character(lasso_prediction$predictions),
        score_CONTROL = lasso_prediction$probabilities[, "CONTROL"],
        score_HEME = lasso_prediction$probabilities[, "HEME"],
        score_NUC = lasso_prediction$probabilities[, "NUC"]
      ),
      data.frame(
        model = "NMM", repetition = repetition, outer_fold = outer_fold,
        descriptor_row_index = outer_test,
        true_class = as.character(y_test),
        predicted_class = as.character(nmm_prediction$predictions),
        score_CONTROL = -nmm_prediction$distances[, "CONTROL"],
        score_HEME = -nmm_prediction$distances[, "HEME"],
        score_NUC = -nmm_prediction$distances[, "NUC"]
      )
    )

    completed <- list(
      input_signature = input_signature,
      training_descriptor_rows = outer_train,
      test_descriptor_rows = outer_test,
      outer = rbind(lasso_outer, nmm_outer),
      inner = rbind(lasso_inner, nmm_inner),
      predictions = prediction_table
    )
    saveRDS(completed, checkpoint_file)
    if (skip_rf) {
      write.csv(data.frame(descriptor_row_index = outer_train,
                          p_base = obj$pairs_df$p_base[outer_train],
                          class = as.character(labels[outer_train])),
                file.path(out_dir, "training_pockets.csv"), row.names = FALSE)
    }
    outer_records[[length(outer_records) + 1L]] <- completed$outer
    inner_records[[length(inner_records) + 1L]] <- completed$inner
    prediction_records[[length(prediction_records) + 1L]] <- completed$predictions

    cat(sprintf(
      "MLR-l1: t_end=%.1f, lambda=%.6g | NMM: t_end=%.1f\n",
      lasso_selected$t_end, lasso_selected$lambda, nmm_selected$t_end
    ))
    rm(candidate_results, lasso_channels, nmm_channels, lasso_scores, nmm_scores)
    invisible(gc())
  }
}

new_outer <- do.call(rbind, outer_records)
new_inner <- do.call(rbind, inner_records)
new_predictions <- do.call(rbind, prediction_records)

if (skip_rf) {
  all_outer <- new_outer
} else {
rf_outer <- read.csv(rf_outer_file, stringsAsFactors = FALSE)
rf_outer$model <- "RF"
rf_outer$selected_lambda <- NA_real_
common_columns <- intersect(names(rf_outer), names(new_outer))
all_outer <- rbind(
  rf_outer[, common_columns, drop = FALSE],
  new_outer[, common_columns, drop = FALSE]
)
}
all_outer <- all_outer[order(all_outer$model, all_outer$comparison,
                             all_outer$repetition, all_outer$outer_fold), ]

summary_table <- do.call(rbind, lapply(
  split(all_outer, interaction(all_outer$model, all_outer$comparison, drop = TRUE)),
  function(x) data.frame(
    model = x$model[1],
    comparison = x$comparison[1],
    n_outer_evaluations = nrow(x),
    mean_accuracy = mean(x$accuracy), sd_accuracy = sd(x$accuracy),
    mean_ppv = mean(x$ppv), sd_ppv = sd(x$ppv),
    mean_tpr = mean(x$tpr), sd_tpr = sd(x$tpr),
    mean_tnr = mean(x$tnr), sd_tnr = sd(x$tnr),
    mean_auc = mean(x$auc), sd_auc = sd(x$auc),
    mean_selected_t_end = mean(x$selected_t_end),
    sd_selected_t_end = sd(x$selected_t_end),
    mean_multiclass_accuracy = mean(x$multiclass_accuracy),
    sd_multiclass_accuracy = sd(x$multiclass_accuracy),
    mean_retained_mfpcs = mean(x$retained_mfpcs),
    sd_retained_mfpcs = sd(x$retained_mfpcs)
  )
))
rownames(summary_table) <- NULL
summary_table$model <- factor(summary_table$model, levels = c("RF", "MLR_l1", "NMM"))
summary_table$comparison <- factor(
  summary_table$comparison, levels = c("NUC_vs_REST", "HEME_vs_REST")
)
summary_table <- summary_table[order(summary_table$comparison, summary_table$model), ]
summary_table$model <- as.character(summary_table$model)
summary_table$comparison <- as.character(summary_table$comparison)

format_pm <- function(mean_value, sd_value) sprintf("%.4f \\pm %.4f", mean_value, sd_value)
manuscript_table <- data.frame(
  comparison = summary_table$comparison,
  classifier = c(RF = "RF", MLR_l1 = "MLR-l1", NMM = "NMM")[summary_table$model],
  ACC = mapply(format_pm, summary_table$mean_accuracy, summary_table$sd_accuracy),
  PPV = mapply(format_pm, summary_table$mean_ppv, summary_table$sd_ppv),
  TPR = mapply(format_pm, summary_table$mean_tpr, summary_table$sd_tpr),
  TNR = mapply(format_pm, summary_table$mean_tnr, summary_table$sd_tnr),
  AUC = mapply(format_pm, summary_table$mean_auc, summary_table$sd_auc),
  selected_t_max_A = mapply(
    format_pm, summary_table$mean_selected_t_end, summary_table$sd_selected_t_end
  ),
  stringsAsFactors = FALSE
)

write.csv(all_outer, file.path(out_dir, "TOUGHC1_FUSED_threeclass_all_models_outer_metrics.csv"), row.names = FALSE)
write.csv(new_inner, file.path(out_dir, "TOUGHC1_FUSED_threeclass_MLR_NMM_inner_metrics.csv"), row.names = FALSE)
write.csv(new_predictions, file.path(out_dir, "TOUGHC1_FUSED_threeclass_MLR_NMM_predictions.csv"), row.names = FALSE)
write.csv(summary_table, file.path(out_dir, "TOUGHC1_FUSED_threeclass_all_models_summary.csv"), row.names = FALSE)
write.csv(manuscript_table, file.path(out_dir, "TOUGHC1_FUSED_threeclass_all_models_manuscript_table.csv"), row.names = FALSE)

cat("\n===== TOUGH-C1 FUSED three-class one-versus-rest summary =====\n")
print(summary_table, row.names = FALSE)
cat(sprintf("\nElapsed script-session time: %.2f minutes\n",
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
cat("Saved outputs under:", out_dir, "\n")
