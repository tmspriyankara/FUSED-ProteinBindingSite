# =====================================================================
# Extended Kahraman | Single-threshold feature-set comparison
#   - Compares ILR + CDPA, ILR only, and CDPA only at fixed thresholds.
#   - Repeated stratified 5-fold CV: 20 repeats x 5 folds = 100 paired folds.
#   - Same folds are used for all feature sets within each classifier.
#   - Pairwise sign-flip tests compare fold-wise accuracy differences.
# =====================================================================

rm(list = ls()); invisible(gc())
set.seed(1)
options(repos = c(CRAN = "https://cran.rstudio.com/"))

suppressPackageStartupMessages({
  library(ranger)
  library(glmnet)
  library(MASS)
})

# ------------------- paths -------------------
cache_dir <- file.path("data", "cache")
out_dir   <- file.path("results", "EK")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw_cache_file <- file.path(cache_dir, "EK_raw_ILR_CDPA_thr4p8_20p0_by0p1.rds")
sulfur_cache_file <- file.path(cache_dir, "EK_sulfur_present_thr4p8_20p0_by0p1.rds")

if (!file.exists(raw_cache_file)) {
  stop("Missing raw descriptor cache: ", raw_cache_file)
}

# ------------------- analysis settings -------------------
thresholds_to_run <- c(5.3, 12.0)

N_REPEATS <- 20
K_FOLDS <- 5

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

standardize_train_test <- function(Xtr, Xte) {
  mu <- colMeans(Xtr)
  sdv <- apply(Xtr, 2, sd)
  sdv[sdv == 0] <- 1
  list(
    Xtr = scale(Xtr, center = mu, scale = sdv),
    Xte = scale(Xte, center = mu, scale = sdv)
  )
}

rf_predict_accuracy <- function(Xtr, y_tr, Xte, y_te, seed = 1) {
  p <- ncol(Xtr)
  mtry_use <- max(1, floor(sqrt(p)))

  fit <- ranger(
    dependent.variable.name = "y",
    data = data.frame(y = y_tr, Xtr, check.names = FALSE),
    num.trees = RF_NUM_TREES,
    mtry = mtry_use,
    min.node.size = RF_MIN_NODE_SIZE,
    probability = TRUE,
    oob.error = FALSE,
    seed = seed
  )

  pp <- predict(fit, data = data.frame(Xte, check.names = FALSE))$predictions
  pred <- factor(colnames(pp)[max.col(pp)], levels = levels(y_te))
  mean(pred == y_te)
}

lasso_predict_accuracy <- function(Xtr, y_tr, Xte, y_te) {
  lambda_path <- sort(unique(c(1, 1e-1, 1e-2, LASSO_LAMBDA)), decreasing = TRUE)

  fit <- glmnet(
    x = as.matrix(Xtr),
    y = y_tr,
    family = "multinomial",
    alpha = 1,
    lambda = lambda_path,
    standardize = FALSE,
    control = list(maxit = 1000000)
  )

  pred <- predict(
    fit,
    newx = as.matrix(Xte),
    s = LASSO_LAMBDA,
    type = "class"
  )
  pred <- factor(as.vector(pred), levels = levels(y_te))
  mean(pred == y_te)
}

mah_predict_accuracy <- function(Xtr, y_tr, Xte, y_te) {
  cls_tr <- levels(droplevels(y_tr))

  Sigma <- cov(Xtr)
  if (is.null(dim(Sigma))) Sigma <- matrix(Sigma, nrow = 1, ncol = 1)
  Sigma_reg <- Sigma + diag(MAH_EPS, ncol(Sigma))
  Sigma_inv <- tryCatch(
    solve(Sigma_reg),
    error = function(e) MASS::ginv(Sigma_reg)
  )

  mu_list <- lapply(cls_tr, function(cl) {
    colMeans(Xtr[y_tr == cl, , drop = FALSE])
  })
  names(mu_list) <- cls_tr

  Dmat <- matrix(NA_real_, nrow = nrow(Xte), ncol = length(cls_tr))
  colnames(Dmat) <- cls_tr
  for (j in seq_along(cls_tr)) {
    diff <- sweep(Xte, 2, mu_list[[cls_tr[j]]], "-")
    Dmat[, j] <- rowSums((diff %*% Sigma_inv) * diff)
  }

  pred <- factor(cls_tr[max.col(-Dmat)], levels = levels(y_te))
  mean(pred == y_te)
}

extract_features_at_threshold <- function(ilr_arr, cov_arr, thr_full, target_r) {
  idx <- which.min(abs(thr_full - target_r))
  actual_r <- thr_full[idx]

  list(
    radius = actual_r,
    features = list(
      ILR_CDPA = cbind(
        ILR1 = ilr_arr[, idx, 1],
        ILR2 = ilr_arr[, idx, 2],
        c12  = cov_arr[, idx, 1],
        c13  = cov_arr[, idx, 2],
        c23  = cov_arr[, idx, 3],
        v1   = cov_arr[, idx, 4],
        v2   = cov_arr[, idx, 5],
        v3   = cov_arr[, idx, 6]
      ),
      ILR_only = cbind(
        ILR1 = ilr_arr[, idx, 1],
        ILR2 = ilr_arr[, idx, 2]
      ),
      CDPA_only = cbind(
        c12 = cov_arr[, idx, 1],
        c13 = cov_arr[, idx, 2],
        c23 = cov_arr[, idx, 3],
        v1  = cov_arr[, idx, 4],
        v2  = cov_arr[, idx, 5],
        v3  = cov_arr[, idx, 6]
      )
    )
  )
}

paired_test_summary <- function(df_pair) {
  d <- df_pair$accuracy_a - df_pair$accuracy_b
  obs <- mean(d)

  set.seed(9101)
  perm_mean <- replicate(N_PERM, mean(d * sample(c(-1, 1), length(d), replace = TRUE)))
  p_two_sided <- (sum(abs(perm_mean) >= abs(obs)) + 1) / (N_PERM + 1)

  set.seed(9102)
  boot_mean <- replicate(N_BOOT, mean(sample(d, length(d), replace = TRUE)))
  ci <- as.numeric(quantile(boot_mean, c(0.025, 0.975), names = FALSE))

  data.frame(
    threshold = df_pair$threshold[1],
    classifier = df_pair$classifier[1],
    comparison = df_pair$comparison[1],
    n_repeats = N_REPEATS,
    n_folds = K_FOLDS,
    n_paired_evals = length(d),
    mean_accuracy_a = mean(df_pair$accuracy_a),
    sd_accuracy_a = sd(df_pair$accuracy_a),
    mean_accuracy_b = mean(df_pair$accuracy_b),
    sd_accuracy_b = sd(df_pair$accuracy_b),
    mean_paired_difference = obs,
    sd_paired_difference = sd(d),
    paired_difference_ci_lower = ci[1],
    paired_difference_ci_upper = ci[2],
    signflip_p_two_sided = p_two_sided,
    stringsAsFactors = FALSE
  )
}

make_pairwise_results <- function(fold_df) {
  pairs <- list(
    c("ILR_CDPA", "ILR_only"),
    c("ILR_CDPA", "CDPA_only"),
    c("CDPA_only", "ILR_only")
  )

  out <- list()
  split_keys <- split(fold_df, list(fold_df$threshold, fold_df$classifier), drop = TRUE)

  for (df_sub in split_keys) {
    for (pair in pairs) {
      a <- pair[1]
      b <- pair[2]
      df_a <- df_sub[df_sub$feature_set == a, ]
      df_b <- df_sub[df_sub$feature_set == b, ]
      df_a <- df_a[order(df_a$repeat_id, df_a$fold_id), ]
      df_b <- df_b[order(df_b$repeat_id, df_b$fold_id), ]

      stopifnot(all(df_a$repeat_id == df_b$repeat_id))
      stopifnot(all(df_a$fold_id == df_b$fold_id))

      df_pair <- data.frame(
        threshold = df_a$threshold,
        classifier = df_a$classifier,
        comparison = paste(a, "minus", b),
        repeat_id = df_a$repeat_id,
        fold_id = df_a$fold_id,
        accuracy_a = df_a$accuracy,
        accuracy_b = df_b$accuracy,
        stringsAsFactors = FALSE
      )
      out[[length(out) + 1L]] <- paired_test_summary(df_pair)
    }
  }

  do.call(rbind, out)
}

# =====================================================================
# load cached descriptors
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

cat("Loaded EK cache with", length(labels), "structures.\n")
print(table(labels))

# =====================================================================
# repeated paired CV
# =====================================================================
fold_results <- list()

for (target_r in thresholds_to_run) {
  feat <- extract_features_at_threshold(ilr_arr, cov_arr, thr_full, target_r)
  actual_r <- feat$radius

  cat("\n############################################################\n")
  cat("Single-threshold comparison at", actual_r, "Angstrom\n")
  cat("############################################################\n")

  ok_all <- Reduce(`&`, lapply(feat$features, function(X) {
    apply(X, 1, function(z) all(is.finite(z)))
  }))

  y_use <- droplevels(labels[ok_all])
  X_list <- lapply(feat$features, function(X) X[ok_all, , drop = FALSE])

  cat("Complete-case sample size:", length(y_use), "\n")
  print(table(y_use))

  for (rep_id in seq_len(N_REPEATS)) {
    cat("\nThreshold", actual_r, "| Repeat", rep_id, "of", N_REPEATS, "\n")
    folds <- make_strat_folds(y_use, K_FOLDS, seed = 10000 + round(actual_r * 10) + rep_id)

    for (fold_id in seq_len(K_FOLDS)) {
      test_idx <- folds[[fold_id]]
      train_idx <- setdiff(seq_along(y_use), test_idx)
      y_tr <- droplevels(y_use[train_idx])
      y_te <- factor(y_use[test_idx], levels = levels(y_tr))

      for (feature_set in names(X_list)) {
        X <- X_list[[feature_set]]
        sc <- standardize_train_test(X[train_idx, , drop = FALSE], X[test_idx, , drop = FALSE])

        seed_use <- 11000 + round(actual_r * 10) * 1000 + rep_id * 100 + fold_id
        rf_acc <- rf_predict_accuracy(sc$Xtr, y_tr, sc$Xte, y_te, seed = seed_use)
        lasso_acc <- lasso_predict_accuracy(sc$Xtr, y_tr, sc$Xte, y_te)
        mah_acc <- mah_predict_accuracy(sc$Xtr, y_tr, sc$Xte, y_te)

        fold_results[[length(fold_results) + 1L]] <- data.frame(
          threshold = actual_r,
          repeat_id = rep_id,
          fold_id = fold_id,
          feature_set = feature_set,
          classifier = c("RF", "MultinomialLasso", "Mahalanobis"),
          accuracy = c(rf_acc, lasso_acc, mah_acc),
          stringsAsFactors = FALSE
        )
      }

      cat(sprintf("Threshold %.1f | Repeat %d | Fold %d complete\n", actual_r, rep_id, fold_id))
      invisible(gc())
    }
  }
}

fold_df <- do.call(rbind, fold_results)

accuracy_summary <- aggregate(
  accuracy ~ threshold + feature_set + classifier,
  data = fold_df,
  FUN = function(z) c(mean = mean(z), sd = sd(z), n = length(z))
)
stats <- as.data.frame(accuracy_summary$accuracy)
names(stats) <- c("mean_accuracy", "sd_accuracy", "n_evals")
accuracy_summary$accuracy <- NULL
accuracy_summary <- cbind(accuracy_summary, stats)

pairwise_summary <- make_pairwise_results(fold_df)

write.csv(
  fold_df,
  file.path(out_dir, "EK_singleThreshold_featureSetComparison_repeated20x5_fold_results.csv"),
  row.names = FALSE
)
write.csv(
  accuracy_summary,
  file.path(out_dir, "EK_singleThreshold_featureSetComparison_repeated20x5_accuracy_summary.csv"),
  row.names = FALSE
)
write.csv(
  pairwise_summary,
  file.path(out_dir, "EK_singleThreshold_featureSetComparison_repeated20x5_pairwise_tests.csv"),
  row.names = FALSE
)

cat("\n===== Single-threshold feature-set accuracy summary =====\n")
print(accuracy_summary)

cat("\n===== Pairwise sign-flip tests =====\n")
print(pairwise_summary)

cat("\nSaved outputs under:", out_dir, "\n")
