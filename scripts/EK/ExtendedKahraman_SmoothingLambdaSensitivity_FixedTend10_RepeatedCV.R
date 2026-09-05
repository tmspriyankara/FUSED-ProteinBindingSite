# =====================================================================
# Extended Kahraman | B-spline smoothing lambda sensitivity
#   - Fixed t_start = 4.8 A and fixed t_end = 10.0 A.
#   - Compares smoothing lambda values: 1e-4, 1e-3, 1e-2, 1e-1.
#   - Reference smoothing lambda is 1e-3.
#   - Repeated stratified 5-fold CV: 20 repeats x 5 folds = 100 paired folds.
#   - Same folds are used for all smoothing lambda values.
# =====================================================================

rm(list = ls()); invisible(gc())
set.seed(1)
options(repos = c(CRAN = "https://cran.rstudio.com/"))

suppressPackageStartupMessages({
  library(fda)
  library(ranger)
  library(glmnet)
})

# ------------------- paths -------------------
cache_dir <- file.path("data", "cache")
out_dir   <- file.path("results", "EK")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw_cache_file <- file.path(cache_dir, "EK_raw_ILR_CDPA_thr4p8_20p0_by0p1.rds")
sulfur_cache_file <- file.path(cache_dir, "EK_sulfur_present_thr4p8_20p0_by0p1.rds")

if (!file.exists(raw_cache_file)) {
  stop("Missing raw descriptor cache: ", raw_cache_file,
       "\nRun the EK repeated nested CV script first to build it.")
}

# ------------------- analysis settings -------------------
t_start <- 4.8
t_end_fixed <- 10.0

N_REPEATS <- 20
K_FOLDS <- 5

MAX_PC     <- 20
VAR_TARGET <- 0.95

SPLINE_DEGREE <- 3
SMOOTH_LAMBDA_GRID <- c(`1e-4` = 1e-4, `1e-3` = 1e-3, `1e-2` = 1e-2, `1e-1` = 1e-1)
REFERENCE_LAMBDA_NAME <- "1e-3"

BASIS_DIVS <- c(2, 3, 4, 5)
NORDER     <- SPLINE_DEGREE + 1

RF_NUM_TREES <- 400
RF_MIN_NODE_SIZE <- 5
LASSO_LAMBDA <- 1e-3
MAH_EPS <- 1e-6

N_PERM <- 10000
N_BOOT <- 10000

# =====================================================================
# helpers
# =====================================================================
make_strat_folds <- function(y, K, seed = 1) {
  set.seed(seed)
  idx_by_class <- split(seq_along(y), y)
  folds <- vector("list", K)

  for (cl in names(idx_by_class)) {
    idx <- sample(idx_by_class[[cl]])
    grp <- rep(seq_len(K), length.out = length(idx))
    for (k in seq_len(K)) {
      folds[[k]] <- c(folds[[k]], idx[grp == k])
    }
  }

  lapply(folds, sort)
}

choose_nbasis <- function(tvec) {
  n_t <- length(tvec)
  cand_nb <- sort(unique(pmax(NORDER, round(n_t / BASIS_DIVS))))
  cand_nb[1]
}

smooth_channel <- function(mat, tvec, nb, norder, lambda) {
  b <- create.bspline.basis(range(tvec), nb, norder)
  fdP <- fdPar(b, int2Lfd(2), lambda)
  sm <- smooth.basis(tvec, t(mat), fdP)
  t(eval.fd(tvec, sm$fd))
}

smooth_all_channels <- function(ilr_arr, cov_arr, idx_use, thr_full, smooth_lambda) {
  t_vec <- thr_full[idx_use]
  nb_use <- choose_nbasis(t_vec)

  list(
    t_vec = t_vec,
    nbasis = nb_use,
    smooth_lambda = smooth_lambda,
    channels = list(
      ilr1 = smooth_channel(ilr_arr[, idx_use, 1], t_vec, nb_use, NORDER, smooth_lambda),
      ilr2 = smooth_channel(ilr_arr[, idx_use, 2], t_vec, nb_use, NORDER, smooth_lambda),
      c12  = smooth_channel(cov_arr[, idx_use, 1], t_vec, nb_use, NORDER, smooth_lambda),
      c13  = smooth_channel(cov_arr[, idx_use, 2], t_vec, nb_use, NORDER, smooth_lambda),
      c23  = smooth_channel(cov_arr[, idx_use, 3], t_vec, nb_use, NORDER, smooth_lambda),
      v1   = smooth_channel(cov_arr[, idx_use, 4], t_vec, nb_use, NORDER, smooth_lambda),
      v2   = smooth_channel(cov_arr[, idx_use, 5], t_vec, nb_use, NORDER, smooth_lambda),
      v3   = smooth_channel(cov_arr[, idx_use, 6], t_vec, nb_use, NORDER, smooth_lambda)
    )
  )
}

fit_pca_scores <- function(Xtr, Xte) {
  center <- colMeans(Xtr)
  Xtr_c <- sweep(Xtr, 2, center, "-")
  Xte_c <- sweep(Xte, 2, center, "-")

  Sigma <- cov(Xtr_c)
  eg <- eigen(Sigma, symmetric = TRUE)
  keep <- which(eg$values > .Machine$double.eps)
  values <- pmax(eg$values[keep], 0)
  rotation <- eg$vectors[, keep, drop = FALSE]

  list(
    Ztr_full = Xtr_c %*% rotation,
    Zte_full = Xte_c %*% rotation,
    eigenvalues = values
  )
}

build_mfpca_scores <- function(smoothed, train_idx, test_idx) {
  ch <- smoothed$channels

  ilr_energy <- mean(abs(c(ch$ilr1[train_idx, ], ch$ilr2[train_idx, ])))
  cdpa_energy <- mean(abs(c(
    ch$c12[train_idx, ], ch$c13[train_idx, ], ch$c23[train_idx, ],
    ch$v1[train_idx, ], ch$v2[train_idx, ], ch$v3[train_idx, ]
  )))
  scale_ilr <- if (ilr_energy == 0) 1 else cdpa_energy / ilr_energy

  Xtr_fun <- cbind(
    ch$ilr1[train_idx, ] * scale_ilr,
    ch$ilr2[train_idx, ] * scale_ilr,
    ch$c12[train_idx, ], ch$c13[train_idx, ], ch$c23[train_idx, ],
    ch$v1[train_idx, ], ch$v2[train_idx, ], ch$v3[train_idx, ]
  )
  Xte_fun <- cbind(
    ch$ilr1[test_idx, ] * scale_ilr,
    ch$ilr2[test_idx, ] * scale_ilr,
    ch$c12[test_idx, ], ch$c13[test_idx, ], ch$c23[test_idx, ],
    ch$v1[test_idx, ], ch$v2[test_idx, ], ch$v3[test_idx, ]
  )

  if (!all(is.finite(Xtr_fun)) || !all(is.finite(Xte_fun))) {
    stop("Non-finite functional features encountered before PCA")
  }

  pca_fit <- fit_pca_scores(Xtr_fun, Xte_fun)
  eigenvalues <- pca_fit$eigenvalues
  prop_var <- eigenvalues / sum(eigenvalues)
  cum_var <- cumsum(prop_var)
  k_use <- which(cum_var >= VAR_TARGET)[1]
  if (is.na(k_use)) k_use <- length(cum_var)
  k_use <- min(k_use, MAX_PC)

  Ztr <- pca_fit$Ztr_full[, seq_len(k_use), drop = FALSE]
  Zte <- pca_fit$Zte_full[, seq_len(k_use), drop = FALSE]
  colnames(Ztr) <- paste0("PC", seq_len(k_use))
  colnames(Zte) <- paste0("PC", seq_len(k_use))

  list(Ztr = Ztr, Zte = Zte, k_use = k_use)
}

rf_predict_accuracy <- function(Ztr, y_tr, Zte, y_te, seed = 1) {
  p <- ncol(Ztr)
  mtry_use <- max(1, floor(sqrt(p)))
  df_tr <- data.frame(y = y_tr, Ztr, check.names = FALSE)

  fit <- ranger(
    dependent.variable.name = "y",
    data = df_tr,
    num.trees = RF_NUM_TREES,
    mtry = mtry_use,
    min.node.size = RF_MIN_NODE_SIZE,
    probability = TRUE,
    oob.error = FALSE,
    seed = seed
  )

  pp <- predict(fit, data = data.frame(Zte, check.names = FALSE))$predictions
  pred <- factor(colnames(pp)[max.col(pp)], levels = levels(y_te))
  mean(pred == y_te)
}

lasso_predict_accuracy <- function(Ztr, y_tr, Zte, y_te) {
  lambda_path <- sort(unique(c(1, 1e-1, 1e-2, LASSO_LAMBDA)), decreasing = TRUE)
  fit <- glmnet(
    x = as.matrix(Ztr),
    y = y_tr,
    family = "multinomial",
    alpha = 1,
    lambda = lambda_path,
    standardize = TRUE,
    control = list(maxit = 1000000)
  )

  pred <- predict(
    fit,
    newx = as.matrix(Zte),
    s = LASSO_LAMBDA,
    type = "class"
  )
  pred <- factor(as.vector(pred), levels = levels(y_te))
  mean(pred == y_te)
}

mah_predict_accuracy <- function(Ztr, y_tr, Zte, y_te) {
  Ztr_mat <- as.matrix(Ztr)
  Zte_mat <- as.matrix(Zte)
  cls_tr <- levels(droplevels(y_tr))

  mu_list <- lapply(cls_tr, function(cl) {
    colMeans(Ztr_mat[y_tr == cl, , drop = FALSE])
  })
  names(mu_list) <- cls_tr

  Sigma <- cov(Ztr_mat)
  Sigma_reg <- Sigma + diag(MAH_EPS, ncol(Sigma))
  Sigma_inv <- tryCatch(
    solve(Sigma_reg),
    error = function(e) MASS::ginv(Sigma_reg)
  )

  Dmat <- matrix(NA_real_, nrow = nrow(Zte_mat), ncol = length(cls_tr))
  colnames(Dmat) <- cls_tr
  for (j in seq_along(cls_tr)) {
    diff <- sweep(Zte_mat, 2, mu_list[[cls_tr[j]]], "-")
    Dmat[, j] <- rowSums((diff %*% Sigma_inv) * diff)
  }

  pred <- factor(cls_tr[max.col(-Dmat)], levels = levels(y_te))
  mean(pred == y_te)
}

paired_test_summary <- function(df_model) {
  ref <- df_model[df_model$smooth_lambda_name == REFERENCE_LAMBDA_NAME, ]
  ref <- ref[order(ref$repeat_id, ref$fold_id), ]

  out <- list()
  for (lambda_name in setdiff(names(SMOOTH_LAMBDA_GRID), REFERENCE_LAMBDA_NAME)) {
    alt <- df_model[df_model$smooth_lambda_name == lambda_name, ]
    alt <- alt[order(alt$repeat_id, alt$fold_id), ]
    stopifnot(all(ref$repeat_id == alt$repeat_id))
    stopifnot(all(ref$fold_id == alt$fold_id))

    d <- alt$accuracy - ref$accuracy
    obs <- mean(d)

    set.seed(9401)
    perm_mean <- replicate(N_PERM, mean(d * sample(c(-1, 1), length(d), replace = TRUE)))
    p_two_sided <- (sum(abs(perm_mean) >= abs(obs)) + 1) / (N_PERM + 1)

    set.seed(9402)
    boot_mean <- replicate(N_BOOT, mean(sample(d, length(d), replace = TRUE)))
    ci <- as.numeric(quantile(boot_mean, c(0.025, 0.975), names = FALSE))

    out[[length(out) + 1L]] <- data.frame(
      model = df_model$model[1],
      reference_lambda = SMOOTH_LAMBDA_GRID[[REFERENCE_LAMBDA_NAME]],
      comparison_lambda = SMOOTH_LAMBDA_GRID[[lambda_name]],
      n_repeats = N_REPEATS,
      n_folds = K_FOLDS,
      n_paired_evals = length(d),
      reference_mean_accuracy = mean(ref$accuracy),
      reference_sd_accuracy = sd(ref$accuracy),
      reference_mean_k = mean(ref$k_use),
      reference_sd_k = sd(ref$k_use),
      comparison_mean_accuracy = mean(alt$accuracy),
      comparison_sd_accuracy = sd(alt$accuracy),
      comparison_mean_k = mean(alt$k_use),
      comparison_sd_k = sd(alt$k_use),
      mean_paired_difference = obs,
      sd_paired_difference = sd(d),
      paired_difference_ci_lower = ci[1],
      paired_difference_ci_upper = ci[2],
      signflip_p_two_sided = p_two_sided,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}

# =====================================================================
# load cache and smooth for each lambda
# =====================================================================
raw <- readRDS(raw_cache_file)
ilr_arr <- raw$ilr_arr
cov_arr <- raw$cov_arr
labels <- droplevels(raw$labels)

if (!is.null(raw$thr_full)) {
  thr_full <- as.numeric(raw$thr_full)
} else if (file.exists(sulfur_cache_file)) {
  sulfur_raw <- readRDS(sulfur_cache_file)
  thr_full <- as.numeric(sulfur_raw$thr_full)
} else {
  thr_full <- seq(4.8, 20.0, by = 0.1)
}

if (length(thr_full) != dim(ilr_arr)[2] || length(thr_full) != dim(cov_arr)[2]) {
  stop("Threshold grid length does not match descriptor array dimensions.")
}

t_end_idx <- which.min(abs(thr_full - t_end_fixed))
if (abs(thr_full[t_end_idx] - t_end_fixed) > 1e-8) {
  stop("Could not find t_end = 10.0 A in threshold grid.")
}
idx_use <- which(thr_full >= t_start & thr_full <= t_end_fixed)

cat("Loaded EK cache with", length(labels), "structures.\n")
cat("Fixed interval:", t_start, "to", t_end_fixed, "Angstrom\n")
cat("Smoothing lambda grid:", paste(names(SMOOTH_LAMBDA_GRID), collapse = ", "), "\n")
print(table(labels))

smoothed_by_lambda <- list()
for (lambda_name in names(SMOOTH_LAMBDA_GRID)) {
  cat("Smoothing fixed-threshold functional descriptors for lambda =", lambda_name, "\n")
  smoothed_by_lambda[[lambda_name]] <- smooth_all_channels(
    ilr_arr, cov_arr, idx_use, thr_full, SMOOTH_LAMBDA_GRID[[lambda_name]]
  )
}

# =====================================================================
# repeated paired CV
# =====================================================================
fold_results <- list()

for (rep_id in seq_len(N_REPEATS)) {
  cat("\nRepeat", rep_id, "of", N_REPEATS, "\n")
  folds <- make_strat_folds(labels, K_FOLDS, seed = 16000 + rep_id)

  for (fold_id in seq_len(K_FOLDS)) {
    test_idx <- folds[[fold_id]]
    train_idx <- setdiff(seq_along(labels), test_idx)
    y_tr <- droplevels(labels[train_idx])
    y_te <- factor(labels[test_idx], levels = levels(y_tr))

    for (lambda_name in names(SMOOTH_LAMBDA_GRID)) {
      scores <- build_mfpca_scores(smoothed_by_lambda[[lambda_name]], train_idx, test_idx)

      seed_use <- 17000 + 100 * rep_id + fold_id
      rf_acc <- rf_predict_accuracy(scores$Ztr, y_tr, scores$Zte, y_te, seed = seed_use)
      lasso_acc <- lasso_predict_accuracy(scores$Ztr, y_tr, scores$Zte, y_te)
      mah_acc <- mah_predict_accuracy(scores$Ztr, y_tr, scores$Zte, y_te)

      fold_results[[length(fold_results) + 1L]] <- data.frame(
        repeat_id = rep_id,
        fold_id = fold_id,
        model = c("RF", "MultinomialLasso", "Mahalanobis"),
        t_end = t_end_fixed,
        smooth_lambda_name = lambda_name,
        smooth_lambda = SMOOTH_LAMBDA_GRID[[lambda_name]],
        accuracy = c(rf_acc, lasso_acc, mah_acc),
        k_use = scores$k_use,
        stringsAsFactors = FALSE
      )

      cat(sprintf(
        "Repeat %d | Fold %d | lambda=%s | k=%d | RF %.3f | Lasso %.3f | Mah %.3f\n",
        rep_id, fold_id, lambda_name, scores$k_use, rf_acc, lasso_acc, mah_acc
      ))

      rm(scores); invisible(gc())
    }
  }
}

fold_df <- do.call(rbind, fold_results)

accuracy_summary <- aggregate(
  cbind(accuracy, k_use) ~ model + smooth_lambda_name + smooth_lambda,
  data = fold_df,
  FUN = function(z) c(mean = mean(z), sd = sd(z), n = length(z))
)
acc_stats <- as.data.frame(accuracy_summary$accuracy)
k_stats <- as.data.frame(accuracy_summary$k_use)
names(acc_stats) <- c("mean_accuracy", "sd_accuracy", "n_evals")
names(k_stats) <- c("mean_k", "sd_k", "n_k")
accuracy_summary$accuracy <- NULL
accuracy_summary$k_use <- NULL
accuracy_summary <- cbind(accuracy_summary, acc_stats, k_stats[, c("mean_k", "sd_k")])

pairwise_summary <- do.call(
  rbind,
  lapply(split(fold_df, fold_df$model), paired_test_summary)
)
pairwise_summary <- pairwise_summary[
  order(match(pairwise_summary$model, c("RF", "MultinomialLasso", "Mahalanobis")),
        pairwise_summary$comparison_lambda),
]
rownames(pairwise_summary) <- NULL

write.csv(
  fold_df,
  file.path(out_dir, "EK_smoothingLambdaSensitivity_fixedTend10_repeated20x5_fold_results.csv"),
  row.names = FALSE
)
write.csv(
  accuracy_summary,
  file.path(out_dir, "EK_smoothingLambdaSensitivity_fixedTend10_repeated20x5_accuracy_summary.csv"),
  row.names = FALSE
)
write.csv(
  pairwise_summary,
  file.path(out_dir, "EK_smoothingLambdaSensitivity_fixedTend10_repeated20x5_pairwise_tests.csv"),
  row.names = FALSE
)

cat("\n===== EK smoothing lambda sensitivity at fixed t_end = 10 A =====\n")
print(accuracy_summary)

cat("\n===== Pairwise tests versus lambda = 1e-3 =====\n")
print(pairwise_summary)

cat("\nSaved outputs under:", out_dir, "\n")
