# =====================================================================
# Extended Kahraman | Exploratory accuracy vs t_end plot
#   - t_start is fixed at 4.8 A.
#   - t_end is evaluated from 6 to 20 A in 0.2 A increments.
#   - One stratified 5-fold CV is used at each t_end.
#   - No nested t_end selection.
# =====================================================================

rm(list = ls()); invisible(gc())
set.seed(1)

suppressPackageStartupMessages({
  library(fda)
  library(ranger)
  library(glmnet)
})

cache_file <- file.path("data", "cache", "EK_raw_ILR_CDPA_thr4p8_20p0_by0p1.rds")
out_dir <- file.path("results", "EK")
fig_dir <- file.path("figures", "EK")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cache_file)) {
  stop("Missing descriptor cache: ", cache_file,
       "\nRun the EK nested-CV script once to create this cache.")
}

raw <- readRDS(cache_file)
ilr_arr <- raw$ilr_arr
cov_arr <- raw$cov_arr
labels <- droplevels(raw$labels)

thr_full <- seq(4.8, 20.0, by = 0.1)
t_start <- 4.8
t_end_grid <- seq(6.0, 20.0, by = 0.2)

K_FOLD <- 5
MAX_PC <- 20
VAR_TARGET <- 0.95

SPLINE_DEGREE <- 3
SMOOTH_LAMBDA <- 1e-3
BASIS_DIVS <- c(2, 3, 4, 5)
NORDER <- SPLINE_DEGREE + 1

RF_NUM_TREES <- 400
RF_MIN_NODE_SIZE <- 5
LASSO_LAMBDA <- 1e-3
LASSO_LAMBDA_GRID <- 10^seq(-2, -4, length.out = 10)
MAH_EPS <- 1e-6

make_strat_folds <- function(y, K, seed = 1) {
  set.seed(seed)
  idx_by_class <- split(seq_along(y), y)
  folds <- vector("list", K)
  for (cl in names(idx_by_class)) {
    idx <- sample(idx_by_class[[cl]])
    grp <- rep(seq_len(K), length.out = length(idx))
    for (k in seq_len(K)) folds[[k]] <- c(folds[[k]], idx[grp == k])
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

smooth_all_channels <- function(ilr_arr, cov_arr, idx_use) {
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

rf_accuracy <- function(Ztr, y_tr, Zte, y_te, seed = 1) {
  fit <- ranger(
    dependent.variable.name = "y",
    data = data.frame(y = y_tr, Ztr, check.names = FALSE),
    num.trees = RF_NUM_TREES,
    mtry = max(1, floor(sqrt(ncol(Ztr)))),
    min.node.size = RF_MIN_NODE_SIZE,
    probability = TRUE,
    oob.error = FALSE,
    num.threads = 1,
    seed = seed
  )
  pp <- predict(fit, data = data.frame(Zte, check.names = FALSE))$predictions
  pred <- factor(colnames(pp)[max.col(pp)], levels = levels(y_te))
  mean(pred == y_te)
}

lasso_accuracy <- function(Ztr, y_tr, Zte, y_te) {
  lambda_path <- sort(unique(c(1, 0.1, LASSO_LAMBDA_GRID, LASSO_LAMBDA)), decreasing = TRUE)
  fit <- glmnet(
    x = as.matrix(Ztr),
    y = y_tr,
    family = "multinomial",
    alpha = 1,
    lambda = lambda_path,
    standardize = TRUE,
    control = list(maxit = 1e6)
  )
  pred <- predict(fit, newx = as.matrix(Zte), s = LASSO_LAMBDA, type = "class")
  pred <- factor(as.vector(pred), levels = levels(y_te))
  mean(pred == y_te)
}

mah_accuracy <- function(Ztr, y_tr, Zte, y_te) {
  Ztr_mat <- as.matrix(Ztr)
  Zte_mat <- as.matrix(Zte)
  cls_tr <- levels(droplevels(y_tr))
  mu_list <- lapply(cls_tr, function(cl) colMeans(Ztr_mat[y_tr == cl, , drop = FALSE]))
  names(mu_list) <- cls_tr
  Sigma <- cov(Ztr_mat)
  Sigma_reg <- Sigma + diag(MAH_EPS, ncol(Sigma))
  Sigma_inv <- tryCatch(solve(Sigma_reg), error = function(e) MASS::ginv(Sigma_reg))
  Dmat <- matrix(NA_real_, nrow = nrow(Zte_mat), ncol = length(cls_tr))
  colnames(Dmat) <- cls_tr
  for (j in seq_along(cls_tr)) {
    diff <- sweep(Zte_mat, 2, mu_list[[cls_tr[j]]], "-")
    Dmat[, j] <- rowSums((diff %*% Sigma_inv) * diff)
  }
  pred <- factor(cls_tr[max.col(-Dmat)], levels = levels(y_te))
  mean(pred == y_te)
}

folds <- make_strat_folds(labels, K_FOLD, seed = 1001)
records <- list()

for (ii in seq_along(t_end_grid)) {
  t_end <- t_end_grid[ii]
  idx_use <- which(thr_full >= t_start & thr_full <= t_end)
  cat(sprintf("Exploratory EK | t_end %.1f (%d/%d)\n", t_end, ii, length(t_end_grid)))

  smoothed <- smooth_all_channels(ilr_arr, cov_arr, idx_use)

  for (fold in seq_len(K_FOLD)) {
    test_idx <- folds[[fold]]
    train_idx <- setdiff(seq_along(labels), test_idx)
    scores <- build_mfpca_scores(smoothed, train_idx, test_idx)
    y_tr <- droplevels(labels[train_idx])
    y_te <- factor(labels[test_idx], levels = levels(y_tr))

    records[[length(records) + 1L]] <- data.frame(
      t_end = t_end,
      fold = fold,
      classifier = "RF",
      accuracy = rf_accuracy(scores$Ztr, y_tr, scores$Zte, y_te, seed = 3000 + 100 * ii + fold),
      k_use = scores$k_use,
      stringsAsFactors = FALSE
    )
    records[[length(records) + 1L]] <- data.frame(
      t_end = t_end,
      fold = fold,
      classifier = "MultinomialLasso",
      accuracy = lasso_accuracy(scores$Ztr, y_tr, scores$Zte, y_te),
      k_use = scores$k_use,
      stringsAsFactors = FALSE
    )
    records[[length(records) + 1L]] <- data.frame(
      t_end = t_end,
      fold = fold,
      classifier = "Mahalanobis",
      accuracy = mah_accuracy(scores$Ztr, y_tr, scores$Zte, y_te),
      k_use = scores$k_use,
      stringsAsFactors = FALSE
    )
  }
}

fold_df <- do.call(rbind, records)
summary_df <- aggregate(
  cbind(mean_accuracy = accuracy, mean_k_use = k_use) ~ t_end + classifier,
  data = fold_df,
  FUN = mean
)
sd_df <- aggregate(
  cbind(sd_accuracy = accuracy, sd_k_use = k_use) ~ t_end + classifier,
  data = fold_df,
  FUN = sd
)
summary_df <- merge(summary_df, sd_df, by = c("t_end", "classifier"))
summary_df <- summary_df[order(summary_df$classifier, summary_df$t_end), ]

fold_file <- file.path(out_dir, "EK_exploratory_accuracy_vs_tend_5fold_by0p2_fold_results.csv")
summary_file <- file.path(out_dir, "EK_exploratory_accuracy_vs_tend_5fold_by0p2_summary.csv")
plot_file <- file.path(fig_dir, "EK_exploratory_accuracy_vs_tend_5fold_by0p2.png")

write.csv(fold_df, fold_file, row.names = FALSE)
write.csv(summary_df, summary_file, row.names = FALSE)

plot_wide <- reshape(
  summary_df[, c("t_end", "classifier", "mean_accuracy")],
  idvar = "t_end",
  timevar = "classifier",
  direction = "wide"
)
plot_wide <- plot_wide[order(plot_wide$t_end), ]
method_names <- c("RF", "MultinomialLasso", "Mahalanobis")
plot_mat <- as.matrix(plot_wide[paste0("mean_accuracy.", method_names)])
colnames(plot_mat) <- c("RF", "MLR-l1", "NMM")

png(plot_file, width = 7.5, height = 5.0, units = "in", res = 600)
op <- par(mar = c(5, 5, 1.5, 1.5), cex.lab = 1.25, cex.axis = 1.05)
plot_cols <- c("#000000", "#1b9e77", "#d95f02")
plot_pch <- rep(16, 3)
plot_lty <- rep(1, 3)
matplot(
  x = plot_wide$t_end,
  y = plot_mat,
  type = "b",
  pch = plot_pch,
  lty = plot_lty,
  lwd = 2.2,
  col = plot_cols,
  xaxt = "n",
  xlab = expression(t[end]),
  ylab = "Mean CV accuracy"
)
axis(1, at = seq(6, 20, by = 2), labels = seq(6, 20, by = 2))
grid(col = "gray85", lty = 1)
legend(
  "right",
  legend = expression(RF, MLR-l[1], NMM),
  col = plot_cols,
  pch = plot_pch,
  lty = plot_lty,
  lwd = 2.2,
  bty = "n",
  cex = 1.0,
  inset = c(0.02, 0)
)
par(op)
dev.off()

cat("\n===== EK exploratory accuracy vs t_end summary =====\n")
print(summary_df)
cat("\nSaved outputs:\n")
cat(" ", fold_file, "\n")
cat(" ", summary_file, "\n")
cat(" ", plot_file, "\n")
