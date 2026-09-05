# =====================================================================
# Extended Kahraman | Sulfur sensitivity analysis at fixed t_end = 10 A
#   - Fixed t_start = 4.8 A and fixed t_end = 10.0 A.
#   - Repeated stratified 5-fold CV: 20 repeats x 5 folds = 100 pairs.
#   - Same folds are used for baseline FUSED and FUSED + S_present_at_10A.
#   - Paired sign-flip permutation tests assess whether sulfur changes accuracy.
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
if (!file.exists(sulfur_cache_file)) {
  stop("Missing sulfur indicator cache: ", sulfur_cache_file,
       "\nRun the EK sulfur indicator script first to build it.")
}

# ------------------- analysis settings -------------------
t_start <- 4.8
t_end_fixed <- 10.0

N_REPEATS <- 20
K_FOLDS <- 5

MAX_PC     <- 20
VAR_TARGET <- 0.95

SPLINE_DEGREE <- 3
SMOOTH_LAMBDA <- 1e-3
BASIS_DIVS    <- c(2, 3, 4, 5)
NORDER        <- SPLINE_DEGREE + 1

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

smooth_all_channels <- function(ilr_arr, cov_arr, idx_use, thr_full) {
  t_vec <- thr_full[idx_use]
  nb_use <- choose_nbasis(t_vec)

  list(
    t_vec = t_vec,
    nbasis = nb_use,
    channels = list(
      ilr1 = smooth_channel(ilr_arr[, idx_use, 1], t_vec, nb_use, NORDER, SMOOTH_LAMBDA),
      ilr2 = smooth_channel(ilr_arr[, idx_use, 2], t_vec, nb_use, NORDER, SMOOTH_LAMBDA),
      c12  = smooth_channel(cov_arr[, idx_use, 1], t_vec, nb_use, NORDER, SMOOTH_LAMBDA),
      c13  = smooth_channel(cov_arr[, idx_use, 2], t_vec, nb_use, NORDER, SMOOTH_LAMBDA),
      c23  = smooth_channel(cov_arr[, idx_use, 3], t_vec, nb_use, NORDER, SMOOTH_LAMBDA),
      v1   = smooth_channel(cov_arr[, idx_use, 4], t_vec, nb_use, NORDER, SMOOTH_LAMBDA),
      v2   = smooth_channel(cov_arr[, idx_use, 5], t_vec, nb_use, NORDER, SMOOTH_LAMBDA),
      v3   = smooth_channel(cov_arr[, idx_use, 6], t_vec, nb_use, NORDER, SMOOTH_LAMBDA)
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

add_sulfur_indicator <- function(scores, sulfur_train, sulfur_test) {
  Ztr <- cbind(scores$Ztr, S_present_at_10A = as.numeric(sulfur_train))
  Zte <- cbind(scores$Zte, S_present_at_10A = as.numeric(sulfur_test))
  list(Ztr = Ztr, Zte = Zte, k_use = scores$k_use)
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
  d <- df_model$sulfur_accuracy - df_model$baseline_accuracy
  obs <- mean(d)

  set.seed(9001)
  perm_mean <- replicate(N_PERM, mean(d * sample(c(-1, 1), length(d), replace = TRUE)))
  p_two_sided <- (sum(abs(perm_mean) >= abs(obs)) + 1) / (N_PERM + 1)

  set.seed(9002)
  boot_mean <- replicate(N_BOOT, mean(sample(d, length(d), replace = TRUE)))
  ci <- as.numeric(quantile(boot_mean, c(0.025, 0.975), names = FALSE))

  data.frame(
    model = df_model$model[1],
    n_repeats = N_REPEATS,
    n_folds = K_FOLDS,
    n_paired_evals = length(d),
    t_end = t_end_fixed,
    baseline_mean_accuracy = mean(df_model$baseline_accuracy),
    baseline_sd_accuracy = sd(df_model$baseline_accuracy),
    sulfur_mean_accuracy = mean(df_model$sulfur_accuracy),
    sulfur_sd_accuracy = sd(df_model$sulfur_accuracy),
    mean_paired_difference = obs,
    sd_paired_difference = sd(d),
    paired_difference_ci_lower = ci[1],
    paired_difference_ci_upper = ci[2],
    signflip_p_two_sided = p_two_sided,
    stringsAsFactors = FALSE
  )
}

# =====================================================================
# load caches and build fixed-threshold representation
# =====================================================================
raw <- readRDS(raw_cache_file)
ilr_arr <- raw$ilr_arr
cov_arr <- raw$cov_arr
labels <- droplevels(raw$labels)

sulfur_raw <- readRDS(sulfur_cache_file)
sulfur_present_arr <- sulfur_raw$sulfur_present_arr
thr_full <- if (!is.null(raw$thr_full)) as.numeric(raw$thr_full) else as.numeric(sulfur_raw$thr_full)

if (length(thr_full) != dim(ilr_arr)[2] || length(thr_full) != dim(cov_arr)[2]) {
  stop("Threshold grid length does not match descriptor array dimensions.")
}
if (length(thr_full) != ncol(sulfur_present_arr)) {
  stop("Threshold grid length does not match sulfur indicator cache dimensions.")
}

t_end_idx <- which.min(abs(thr_full - t_end_fixed))
if (abs(thr_full[t_end_idx] - t_end_fixed) > 1e-8) {
  stop("Could not find t_end = 10.0 A in threshold grid.")
}
idx_use <- which(thr_full >= t_start & thr_full <= t_end_fixed)
sulfur_at_10 <- sulfur_present_arr[, t_end_idx]

cat("Loaded EK cache with", length(labels), "structures.\n")
cat("Fixed interval:", t_start, "to", t_end_fixed, "Angstrom\n")
cat("Sulfur present at 10 A:", sum(sulfur_at_10), "of", length(sulfur_at_10), "\n")
print(table(labels))

cat("Smoothing fixed-threshold functional descriptors...\n")
smoothed <- smooth_all_channels(ilr_arr, cov_arr, idx_use, thr_full)

# =====================================================================
# repeated paired CV
# =====================================================================
fold_results <- list()

for (rep_id in seq_len(N_REPEATS)) {
  cat("\nRepeat", rep_id, "of", N_REPEATS, "\n")
  folds <- make_strat_folds(labels, K_FOLDS, seed = 5000 + rep_id)

  for (fold_id in seq_len(K_FOLDS)) {
    test_idx <- folds[[fold_id]]
    train_idx <- setdiff(seq_along(labels), test_idx)
    y_tr <- droplevels(labels[train_idx])
    y_te <- factor(labels[test_idx], levels = levels(y_tr))

    scores <- build_mfpca_scores(smoothed, train_idx, test_idx)
    scores_s <- add_sulfur_indicator(
      scores,
      sulfur_train = sulfur_at_10[train_idx],
      sulfur_test = sulfur_at_10[test_idx]
    )

    rf_base <- rf_predict_accuracy(
      scores$Ztr, y_tr, scores$Zte, y_te,
      seed = 6000 + 100 * rep_id + fold_id
    )
    rf_s <- rf_predict_accuracy(
      scores_s$Ztr, y_tr, scores_s$Zte, y_te,
      seed = 6000 + 100 * rep_id + fold_id
    )

    lasso_base <- lasso_predict_accuracy(scores$Ztr, y_tr, scores$Zte, y_te)
    lasso_s <- lasso_predict_accuracy(scores_s$Ztr, y_tr, scores_s$Zte, y_te)

    mah_base <- mah_predict_accuracy(scores$Ztr, y_tr, scores$Zte, y_te)
    mah_s <- mah_predict_accuracy(scores_s$Ztr, y_tr, scores_s$Zte, y_te)

    fold_results[[length(fold_results) + 1L]] <- data.frame(
      repeat_id = rep_id,
      fold_id = fold_id,
      model = c("RF", "MultinomialLasso", "Mahalanobis"),
      t_end = t_end_fixed,
      baseline_accuracy = c(rf_base, lasso_base, mah_base),
      sulfur_accuracy = c(rf_s, lasso_s, mah_s),
      paired_difference = c(rf_s - rf_base, lasso_s - lasso_base, mah_s - mah_base),
      k_use = scores$k_use,
      stringsAsFactors = FALSE
    )

    cat(sprintf(
      "Repeat %d | Fold %d | RF %.3f -> %.3f | Lasso %.3f -> %.3f | Mah %.3f -> %.3f\n",
      rep_id, fold_id, rf_base, rf_s, lasso_base, lasso_s, mah_base, mah_s
    ))

    rm(scores, scores_s); invisible(gc())
  }
}

fold_df <- do.call(rbind, fold_results)
summary_df <- do.call(
  rbind,
  lapply(split(fold_df, fold_df$model), paired_test_summary)
)
summary_df <- summary_df[match(c("RF", "MultinomialLasso", "Mahalanobis"), summary_df$model), ]
rownames(summary_df) <- NULL

write.csv(
  fold_df,
  file.path(out_dir, "EK_sulfurSensitivity_fixedTend10_repeated20x5_fold_results.csv"),
  row.names = FALSE
)
write.csv(
  summary_df,
  file.path(out_dir, "EK_sulfurSensitivity_fixedTend10_repeated20x5_summary.csv"),
  row.names = FALSE
)

cat("\n===== EK sulfur sensitivity at fixed t_end = 10 A =====\n")
print(summary_df)
cat("\nSaved outputs under:", out_dir, "\n")
