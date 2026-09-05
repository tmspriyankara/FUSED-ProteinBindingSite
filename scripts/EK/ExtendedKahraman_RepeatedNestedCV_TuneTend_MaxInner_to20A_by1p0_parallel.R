# =====================================================================
# Extended Kahraman | Repeated nested CV for threshold and model tuning
#   - t_start is fixed at 4.8 A.
#   - t_end and model hyperparameters are tuned only inside inner CV.
#   - Candidate t_end values use a 1.0 A grid from 6.0 to 20.0 A.
#   - Selection rule: maximum mean inner-CV accuracy.
#   - Outer CV is used only for final performance estimation.
# =====================================================================

rm(list = ls()); invisible(gc())
set.seed(1)
options(repos = c(CRAN = "https://cran.rstudio.com/"))
script_start_time <- Sys.time()
script_start_proc <- proc.time()

suppressPackageStartupMessages({
  library(fda)
  library(ranger)
  library(glmnet)
  library(parallel)
})

# ------------------- paths -------------------
root_dir <- file.path("data", "EK", "Extended Kahraman Proteins Sorted")
cache_dir <- file.path("data", "cache")
out_dir   <- file.path("results", "EK")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
timing_file <- file.path(out_dir, "EK_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_timing.csv")

# ------------------- analysis settings -------------------
thr_full   <- seq(4.8, 20.0, by = 0.1)
t_start    <- 4.8
t_end_grid <- seq(6.0, 20.0, by = 1.0)

K_OUTER <- 5
K_INNER <- 5
N_REPEATS <- 5

MAX_PC     <- 20
VAR_TARGET <- 0.95

SPLINE_DEGREE <- 3
SMOOTH_LAMBDA <- 1e-3
BASIS_DIVS    <- c(2, 3, 4, 5)
NORDER        <- SPLINE_DEGREE + 1

RF_NUM_TREES <- 400
RF_MIN_NODE_SIZE_GRID <- c(2, 5, 7)
N_CORES <- as.integer(Sys.getenv("RR_N_CORES", unset = "4"))
available_cores <- parallel::detectCores(logical = FALSE)
if (is.na(available_cores)) available_cores <- parallel::detectCores(logical = TRUE)
if (is.na(available_cores)) available_cores <- 1L
N_CORES <- max(1L, min(N_CORES, available_cores))

LASSO_LAMBDA_GRID <- 10^seq(-2, -4, length.out = 10)
MAH_EPS <- 1e-6

# ------------------- ligand aliases -------------------
ligand_aliases <- list(
  AMP = c("AMP"),
  ATP = c("ATP", "AGS"),
  FAD = c("FAD"),
  FMN = c("FMN"),
  GLC = c("GLC", "IMD", "LAT", "NAG", "BGC"),
  HEM = c("HEM"),
  NAD = c("NAD", "NAI"),
  PO4 = c("PO4")
)
target_classes <- c("AMP", "ATP", "FAD", "FMN", "GLC", "HEM", "NAD", "PO4")

# =====================================================================
# helpers: PDB parsing and descriptors
# =====================================================================
read_pdb_minimal <- function(pdb_path) {
  ln <- readLines(pdb_path, warn = FALSE)
  ln <- ln[grepl("^(ATOM  |HETATM)", ln)]
  if (length(ln) == 0) return(NULL)

  rec   <- substr(ln, 1, 6)
  resn  <- trimws(substr(ln, 18, 20))
  chain <- substr(ln, 22, 22)
  x     <- as.numeric(substr(ln, 31, 38))
  y     <- as.numeric(substr(ln, 39, 46))
  z     <- as.numeric(substr(ln, 47, 54))
  elem_raw <- trimws(substr(ln, 77, 78))
  aname    <- trimws(substr(ln, 13, 16))
  elem <- ifelse(elem_raw == "", toupper(substr(aname, 1, 1)), toupper(elem_raw))

  keepH <- elem != "H"
  data.frame(
    rec = rec[keepH], resn = resn[keepH], chain = chain[keepH],
    x = x[keepH], y = y[keepH], z = z[keepH], elem = elem[keepH],
    stringsAsFactors = FALSE
  )
}

get_site_atoms <- function(atom_df, ligand_vec, dthr, min_atoms = 10) {
  lig <- atom_df[atom_df$resn %in% ligand_vec, , drop = FALSE]
  if (nrow(lig) == 0) return(NULL)

  lig_chain <- lig$chain[1]
  prot <- atom_df[
    atom_df$rec == "ATOM  " &
      atom_df$chain == lig_chain &
      !(atom_df$resn %in% ligand_vec),
    , drop = FALSE
  ]
  if (nrow(prot) == 0) return(NULL)

  pick <- logical(nrow(prot))
  for (i in seq_len(nrow(lig))) {
    dx <- prot$x - lig$x[i]
    dy <- prot$y - lig$y[i]
    dz <- prot$z - lig$z[i]
    pick <- pick | (sqrt(dx*dx + dy*dy + dz*dz) <= dthr)
  }

  site <- prot[pick, , drop = FALSE]
  if (nrow(site) < min_atoms) return(NULL)
  site
}

compute_ilr_noeps <- function(elems) {
  cntC <- sum(elems == "C")
  cntO <- sum(elems == "O")
  cntN <- sum(elems == "N")
  tot  <- cntC + cntO + cntN
  if (tot == 0) return(c(NA_real_, NA_real_))

  xC <- cntC / tot
  xO <- cntO / tot
  xN <- cntN / tot
  if (xC == 0 || xO == 0 || xN == 0) return(c(NA_real_, NA_real_))

  c(
    sqrt(1/2) * log(xC / xO),
    sqrt(2/3) * log(sqrt(xC * xO) / xN)
  )
}

distances_to_local_axes <- function(coords_mat) {
  pr <- prcomp(coords_mat, center = TRUE, scale. = FALSE)
  center <- colMeans(coords_mat)
  axes <- pr$rotation[, 1:3, drop = FALSE]
  dmat <- matrix(NA_real_, nrow(coords_mat), 3)

  for (i in seq_len(nrow(coords_mat))) {
    v <- coords_mat[i, ] - center
    for (k in 1:3) {
      a <- axes[, k]
      resid <- v - sum(v * a) * a
      dmat[i, k] <- sqrt(sum(resid^2))
    }
  }
  dmat
}

cov_from_axis_dist <- function(coords_mat) {
  cv <- cov(distances_to_local_axes(coords_mat))
  c(
    cv[1, 2] * sqrt(2),
    cv[1, 3] * sqrt(2),
    cv[2, 3] * sqrt(2),
    cv[1, 1],
    cv[2, 2],
    cv[3, 3]
  )
}

# =====================================================================
# helpers: folds, smoothing, MFPCA scores, classifiers
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

rf_predict_accuracy <- function(Ztr, y_tr, Zte, y_te, num_trees, min_node_size, seed = 1) {
  p <- ncol(Ztr)
  mtry_use <- max(1, floor(sqrt(p)))

  df_tr <- data.frame(y = y_tr, Ztr, check.names = FALSE)
  fit <- ranger(
    dependent.variable.name = "y",
    data = df_tr,
    num.trees = num_trees,
    mtry = mtry_use,
    min.node.size = min_node_size,
    probability = TRUE,
    oob.error = FALSE,
    num.threads = 1,
    seed = seed
  )

  df_te <- data.frame(Zte, check.names = FALSE)
  pp <- predict(fit, data = df_te)$predictions
  pred <- factor(colnames(pp)[max.col(pp)], levels = levels(y_te))
  mean(pred == y_te)
}

lasso_predict_accuracies <- function(Ztr, y_tr, Zte, y_te, lambda_values) {
  lambda_values <- sort(unique(lambda_values), decreasing = TRUE)
  fit_lambda_values <- sort(unique(c(1, 1e-1, lambda_values)), decreasing = TRUE)
  fit <- glmnet(
    x = as.matrix(Ztr),
    y = y_tr,
    family = "multinomial",
    alpha = 1,
    lambda = fit_lambda_values,
    standardize = TRUE,
    control = list(maxit = 1000000)
  )

  pred <- predict(
    fit,
    newx = as.matrix(Zte),
    s = lambda_values,
    type = "class"
  )
  pred_mat <- as.matrix(pred)
  if (ncol(pred_mat) != length(lambda_values)) {
    pred_mat <- matrix(as.vector(pred), ncol = length(lambda_values))
  }

  vapply(seq_along(lambda_values), function(j) {
    pred_j <- factor(pred_mat[, j], levels = levels(y_te))
    mean(pred_j == y_te)
  }, numeric(1))
}

lasso_predict_accuracy <- function(Ztr, y_tr, Zte, y_te, lambda_value) {
  lambda_path <- sort(unique(c(LASSO_LAMBDA_GRID, lambda_value)), decreasing = TRUE)
  acc <- lasso_predict_accuracies(Ztr, y_tr, Zte, y_te, lambda_path)
  acc[which.min(abs(lambda_path - lambda_value))]
}

mah_predict_accuracy <- function(Ztr, y_tr, Zte, y_te, eps_value) {
  Ztr_mat <- as.matrix(Ztr)
  Zte_mat <- as.matrix(Zte)
  cls_tr <- levels(droplevels(y_tr))

  mu_list <- lapply(cls_tr, function(cl) {
    colMeans(Ztr_mat[y_tr == cl, , drop = FALSE])
  })
  names(mu_list) <- cls_tr

  Sigma <- cov(Ztr_mat)
  Sigma_reg <- Sigma + diag(eps_value, ncol(Sigma))
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

summarize_by_grid <- function(df, group_cols) {
  key <- do.call(interaction, c(df[group_cols], list(drop = TRUE, sep = "|")))
  mean_acc <- tapply(df$accuracy, key, mean)
  best_key <- names(which.max(mean_acc))
  parts <- strsplit(best_key, "\\|", fixed = FALSE)[[1]]
  out <- as.list(parts)
  names(out) <- group_cols
  out$mean_inner_accuracy <- unname(max(mean_acc))
  out
}

# =====================================================================
# build/load raw descriptor cache
# =====================================================================
cache_file <- file.path(cache_dir, "EK_raw_ILR_CDPA_thr4p8_20p0_by0p1.rds")

descriptor_start_time <- Sys.time()
descriptor_start_proc <- proc.time()
if (file.exists(cache_file)) {
  cat("Loading raw descriptor cache:", cache_file, "\n")
  raw <- readRDS(cache_file)
  ilr_arr <- raw$ilr_arr
  cov_arr <- raw$cov_arr
  labels  <- raw$labels
  files_df <- raw$files_df
} else {
  cat("Building raw descriptor cache...\n")

  all_files <- list()
  for (cls in target_classes) {
    fdir <- file.path(root_dir, cls)
    pdbs <- list.files(fdir, pattern = "\\.pdb$", full.names = TRUE)
    if (length(pdbs) > 0) {
      all_files[[cls]] <- data.frame(path = pdbs, label = cls, stringsAsFactors = FALSE)
    }
  }

  files_df <- do.call(rbind, all_files)
  row.names(files_df) <- NULL
  labels <- factor(files_df$label, levels = target_classes)

  n0 <- nrow(files_df)
  n_thr_full <- length(thr_full)
  ilr_arr <- array(NA_real_, dim = c(n0, n_thr_full, 2))
  cov_arr <- array(NA_real_, dim = c(n0, n_thr_full, 6))

  cat("Total PDBs found:", n0, "\n")
  for (i in seq_len(n0)) {
    fpath <- files_df$path[i]
    lbl <- files_df$label[i]
    aliases <- ligand_aliases[[lbl]]
    adf <- tryCatch(read_pdb_minimal(fpath), error = function(e) NULL)
    if (is.null(adf)) next

    ok_all <- TRUE
    for (tt in seq_along(thr_full)) {
      site <- get_site_atoms(adf, aliases, thr_full[tt], min_atoms = 10)
      if (is.null(site)) { ok_all <- FALSE; break }

      ilr12 <- compute_ilr_noeps(site$elem)
      if (any(is.na(ilr12))) { ok_all <- FALSE; break }

      cov6 <- cov_from_axis_dist(as.matrix(site[, c("x", "y", "z")]))
      ilr_arr[i, tt, ] <- ilr12
      cov_arr[i, tt, ] <- cov6
    }

    if (!ok_all) {
      ilr_arr[i, , ] <- NA_real_
      cov_arr[i, , ] <- NA_real_
    }
    if (i %% 50 == 0) cat("processed", i, "of", n0, "\n")
  }

  ok_row <- apply(ilr_arr, 1, function(v) all(!is.na(v)))
  cat("Dropping", sum(!ok_row), "structures with bad site(s)\n")

  ilr_arr <- ilr_arr[ok_row, , , drop = FALSE]
  cov_arr <- cov_arr[ok_row, , , drop = FALSE]
  labels <- droplevels(labels[ok_row])
  files_df <- files_df[ok_row, , drop = FALSE]

  saveRDS(
    list(ilr_arr = ilr_arr, cov_arr = cov_arr, labels = labels, files_df = files_df),
    cache_file
  )
  cat("Saved raw descriptor cache:", cache_file, "\n")
}
descriptor_end_time <- Sys.time()
descriptor_end_proc <- proc.time()
descriptor_elapsed <- descriptor_end_proc - descriptor_start_proc

n <- length(labels)
cat("Used structures:", n, "\n")
print(table(labels))

# =====================================================================
# repeated nested CV
# =====================================================================
modeling_start_time <- Sys.time()
modeling_start_proc <- proc.time()

rf_outer_results <- list()
lasso_outer_results <- list()
mah_outer_results <- list()
rf_inner_records <- list()
lasso_inner_records <- list()
mah_inner_records <- list()

for (rep_id in seq_len(N_REPEATS)) {
  cat("\n############################################################\n")
  cat("Repeat", rep_id, "of", N_REPEATS, "\n")
  cat("############################################################\n")

  outer_folds <- make_strat_folds(labels, K_OUTER, seed = 1000 + rep_id)

for (outer in seq_len(K_OUTER)) {
  cat("\n============================================================\n")
  cat("Repeat", rep_id, "| Outer fold", outer, "of", K_OUTER, "\n")
  cat("============================================================\n")

  outer_test <- outer_folds[[outer]]
  outer_train <- setdiff(seq_len(n), outer_test)
  y_outer_train <- droplevels(labels[outer_train])
  inner_folds_rel <- make_strat_folds(y_outer_train, K_INNER, seed = 2000 + 100 * rep_id + outer)

  cat(sprintf(
    "Repeat %d | Outer %d | evaluating %d candidate t_end values on %d cores\n",
    rep_id, outer, length(t_end_grid), N_CORES
  ))

  evaluate_t_end_candidate <- function(ii) {
    t_end <- t_end_grid[ii]
    idx_use <- which(thr_full >= t_start & thr_full <= t_end)
    if (length(idx_use) < 5) return(NULL)

    smoothed <- smooth_all_channels(ilr_arr, cov_arr, idx_use)
    rf_inner <- list()
    lasso_inner <- list()
    mah_inner <- list()

    for (inner in seq_len(K_INNER)) {
      inner_valid <- outer_train[inner_folds_rel[[inner]]]
      inner_train <- setdiff(outer_train, inner_valid)

      scores <- build_mfpca_scores(smoothed, inner_train, inner_valid)
      y_tr <- droplevels(labels[inner_train])
      y_va <- factor(labels[inner_valid], levels = levels(y_tr))

      for (min_node_size in RF_MIN_NODE_SIZE_GRID) {
        acc <- rf_predict_accuracy(
          scores$Ztr, y_tr, scores$Zte, y_va,
          num_trees = RF_NUM_TREES,
          min_node_size = min_node_size,
          seed = 3000 + 1000 * rep_id + 100 * outer + 10 * inner + min_node_size
        )
        rf_inner[[length(rf_inner) + 1L]] <- data.frame(
          repeat_id = rep_id,
          outer_fold = outer,
          inner_fold = inner,
          t_end = t_end,
          num.trees = RF_NUM_TREES,
          mtry_rule = "floor_sqrt_p",
          min.node.size = min_node_size,
          k_use = scores$k_use,
          mtry = max(1, floor(sqrt(scores$k_use))),
          accuracy = acc,
          stringsAsFactors = FALSE
        )
      }

      lasso_acc <- lasso_predict_accuracies(
        scores$Ztr, y_tr, scores$Zte, y_va,
        lambda_values = LASSO_LAMBDA_GRID
      )
      for (ll in seq_along(LASSO_LAMBDA_GRID)) {
        lambda_value <- LASSO_LAMBDA_GRID[ll]
        lasso_inner[[length(lasso_inner) + 1L]] <- data.frame(
          repeat_id = rep_id,
          outer_fold = outer,
          inner_fold = inner,
          t_end = t_end,
          lambda = lambda_value,
          log10_lambda = log10(lambda_value),
          k_use = scores$k_use,
          accuracy = lasso_acc[ll],
          stringsAsFactors = FALSE
        )
      }

      acc <- mah_predict_accuracy(
        scores$Ztr, y_tr, scores$Zte, y_va,
        eps_value = MAH_EPS
      )
      mah_inner[[length(mah_inner) + 1L]] <- data.frame(
        repeat_id = rep_id,
        outer_fold = outer,
        inner_fold = inner,
        t_end = t_end,
        eps = MAH_EPS,
        log10_eps = log10(MAH_EPS),
        k_use = scores$k_use,
        accuracy = acc,
        stringsAsFactors = FALSE
      )
    }

    list(
      rf = do.call(rbind, rf_inner),
      lasso = do.call(rbind, lasso_inner),
      mah = do.call(rbind, mah_inner)
    )
  }

  candidate_results <- parallel::mclapply(
    seq_along(t_end_grid),
    evaluate_t_end_candidate,
    mc.cores = N_CORES,
    mc.preschedule = FALSE
  )
  candidate_results <- Filter(Negate(is.null), candidate_results)

  rf_inner <- lapply(candidate_results, `[[`, "rf")
  lasso_inner <- lapply(candidate_results, `[[`, "lasso")
  mah_inner <- lapply(candidate_results, `[[`, "mah")

  rf_inner_df <- do.call(rbind, rf_inner)
  lasso_inner_df <- do.call(rbind, lasso_inner)
  mah_inner_df <- do.call(rbind, mah_inner)
  rf_inner_records[[length(rf_inner_records) + 1L]] <- rf_inner_df
  lasso_inner_records[[length(lasso_inner_records) + 1L]] <- lasso_inner_df
  mah_inner_records[[length(mah_inner_records) + 1L]] <- mah_inner_df

  rf_best <- summarize_by_grid(
    rf_inner_df,
    c("t_end", "min.node.size")
  )
  lasso_best <- summarize_by_grid(
    lasso_inner_df,
    c("t_end", "lambda")
  )
  mah_best <- summarize_by_grid(
    mah_inner_df,
    c("t_end")
  )

  rf_best$t_end <- as.numeric(rf_best$t_end)
  rf_best$min.node.size <- as.integer(rf_best$min.node.size)
  lasso_best$t_end <- as.numeric(lasso_best$t_end)
  lasso_best$lambda <- as.numeric(lasso_best$lambda)
  mah_best$t_end <- as.numeric(mah_best$t_end)

  cat("\nSelected RF hyperparameters:\n")
  print(rf_best)
  cat("\nSelected lasso hyperparameters:\n")
  print(lasso_best)
  cat("\nSelected Mahalanobis hyperparameters:\n")
  print(mah_best)

  # Final outer evaluation: refit once on full outer-training data.
  rf_idx_use <- which(thr_full >= t_start & thr_full <= rf_best$t_end)
  rf_smoothed <- smooth_all_channels(ilr_arr, cov_arr, rf_idx_use)
  rf_scores <- build_mfpca_scores(rf_smoothed, outer_train, outer_test)
  y_tr_outer <- droplevels(labels[outer_train])
  y_te_outer <- factor(labels[outer_test], levels = levels(y_tr_outer))

  rf_final_acc <- rf_predict_accuracy(
    rf_scores$Ztr, y_tr_outer, rf_scores$Zte, y_te_outer,
    num_trees = RF_NUM_TREES,
    min_node_size = rf_best$min.node.size,
    seed = 9000 + 100 * rep_id + outer
  )

  lasso_idx_use <- which(thr_full >= t_start & thr_full <= lasso_best$t_end)
  lasso_smoothed <- smooth_all_channels(ilr_arr, cov_arr, lasso_idx_use)
  lasso_scores <- build_mfpca_scores(lasso_smoothed, outer_train, outer_test)
  lasso_final_acc <- lasso_predict_accuracy(
    lasso_scores$Ztr, y_tr_outer, lasso_scores$Zte, y_te_outer,
    lambda_value = lasso_best$lambda
  )

  mah_idx_use <- which(thr_full >= t_start & thr_full <= mah_best$t_end)
  mah_smoothed <- smooth_all_channels(ilr_arr, cov_arr, mah_idx_use)
  mah_scores <- build_mfpca_scores(mah_smoothed, outer_train, outer_test)
  mah_final_acc <- mah_predict_accuracy(
    mah_scores$Ztr, y_tr_outer, mah_scores$Zte, y_te_outer,
    eps_value = MAH_EPS
  )

  rf_outer_results[[length(rf_outer_results) + 1L]] <- data.frame(
    model = "RF",
    repeat_id = rep_id,
    outer_fold = outer,
    t_end = rf_best$t_end,
    num.trees = RF_NUM_TREES,
    mtry_rule = "floor_sqrt_p",
    min.node.size = rf_best$min.node.size,
    mean_inner_accuracy = rf_best$mean_inner_accuracy,
    outer_accuracy = rf_final_acc,
    k_use_final = rf_scores$k_use,
    mtry_final = max(1, floor(sqrt(rf_scores$k_use))),
    stringsAsFactors = FALSE
  )

  lasso_outer_results[[length(lasso_outer_results) + 1L]] <- data.frame(
    model = "MultinomialLasso",
    repeat_id = rep_id,
    outer_fold = outer,
    t_end = lasso_best$t_end,
    lambda = lasso_best$lambda,
    log10_lambda = log10(lasso_best$lambda),
    mean_inner_accuracy = lasso_best$mean_inner_accuracy,
    outer_accuracy = lasso_final_acc,
    k_use_final = lasso_scores$k_use,
    stringsAsFactors = FALSE
  )

  mah_outer_results[[length(mah_outer_results) + 1L]] <- data.frame(
    model = "Mahalanobis",
    repeat_id = rep_id,
    outer_fold = outer,
    t_end = mah_best$t_end,
    eps = MAH_EPS,
    log10_eps = log10(MAH_EPS),
    mean_inner_accuracy = mah_best$mean_inner_accuracy,
    outer_accuracy = mah_final_acc,
    k_use_final = mah_scores$k_use,
    stringsAsFactors = FALSE
  )

  cat(sprintf("\nRepeat %d | Outer %d final accuracy | RF=%.3f | Lasso=%.3f | Mahalanobis=%.3f\n",
              rep_id, outer, rf_final_acc, lasso_final_acc, mah_final_acc))

  rm(rf_smoothed, lasso_smoothed, mah_smoothed, rf_scores, lasso_scores, mah_scores)
  invisible(gc())
}
} # end repeat loop

modeling_end_time <- Sys.time()
modeling_end_proc <- proc.time()
modeling_elapsed <- modeling_end_proc - modeling_start_proc

rf_outer_df <- do.call(rbind, rf_outer_results)
lasso_outer_df <- do.call(rbind, lasso_outer_results)
mah_outer_df <- do.call(rbind, mah_outer_results)

rf_inner_df_all <- do.call(rbind, rf_inner_records)
lasso_inner_df_all <- do.call(rbind, lasso_inner_records)
mah_inner_df_all <- do.call(rbind, mah_inner_records)

write.csv(rf_outer_df, file.path(out_dir, "EK_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_RF_outer_results.csv"), row.names = FALSE)
write.csv(lasso_outer_df, file.path(out_dir, "EK_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_MultinomialLasso_outer_results.csv"), row.names = FALSE)
write.csv(mah_outer_df, file.path(out_dir, "EK_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_Mahalanobis_outer_results.csv"), row.names = FALSE)
write.csv(rf_inner_df_all, file.path(out_dir, "EK_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_RF_inner_results.csv"), row.names = FALSE)
write.csv(lasso_inner_df_all, file.path(out_dir, "EK_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_MultinomialLasso_inner_results.csv"), row.names = FALSE)
write.csv(mah_inner_df_all, file.path(out_dir, "EK_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_Mahalanobis_inner_results.csv"), row.names = FALSE)

summary_df <- rbind(
  data.frame(
    model = "RF",
    selection_rule = "max_inner_accuracy",
    n_repeats = N_REPEATS,
    n_outer_evals = nrow(rf_outer_df),
    mean_outer_accuracy = mean(rf_outer_df$outer_accuracy),
    sd_outer_accuracy = sd(rf_outer_df$outer_accuracy),
    mean_selected_t_end = mean(rf_outer_df$t_end),
    sd_selected_t_end = sd(rf_outer_df$t_end),
    mean_selected_log10_lambda = NA_real_,
    sd_selected_log10_lambda = NA_real_,
    mean_selected_log10_eps = NA_real_,
    sd_selected_log10_eps = NA_real_,
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "MultinomialLasso",
    selection_rule = "max_inner_accuracy",
    n_repeats = N_REPEATS,
    n_outer_evals = nrow(lasso_outer_df),
    mean_outer_accuracy = mean(lasso_outer_df$outer_accuracy),
    sd_outer_accuracy = sd(lasso_outer_df$outer_accuracy),
    mean_selected_t_end = mean(lasso_outer_df$t_end),
    sd_selected_t_end = sd(lasso_outer_df$t_end),
    mean_selected_log10_lambda = mean(lasso_outer_df$log10_lambda),
    sd_selected_log10_lambda = sd(lasso_outer_df$log10_lambda),
    mean_selected_log10_eps = NA_real_,
    sd_selected_log10_eps = NA_real_,
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "Mahalanobis",
    selection_rule = "max_inner_accuracy",
    n_repeats = N_REPEATS,
    n_outer_evals = nrow(mah_outer_df),
    mean_outer_accuracy = mean(mah_outer_df$outer_accuracy),
    sd_outer_accuracy = sd(mah_outer_df$outer_accuracy),
    mean_selected_t_end = mean(mah_outer_df$t_end),
    sd_selected_t_end = sd(mah_outer_df$t_end),
    mean_selected_log10_lambda = NA_real_,
    sd_selected_log10_lambda = NA_real_,
    mean_selected_log10_eps = mean(mah_outer_df$log10_eps),
    sd_selected_log10_eps = sd(mah_outer_df$log10_eps),
    stringsAsFactors = FALSE
  )
)

write.csv(summary_df, file.path(out_dir, "EK_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_summary.csv"), row.names = FALSE)

script_end_time <- Sys.time()
script_end_proc <- proc.time()
script_elapsed <- script_end_proc - script_start_proc
timing_df <- data.frame(
  step = c("DESCRIPTOR_LOAD_OR_BUILD", "NESTED_CV_MODELING", "ALL"),
  start_time = c(
    format(descriptor_start_time, "%Y-%m-%d %H:%M:%S %Z"),
    format(modeling_start_time, "%Y-%m-%d %H:%M:%S %Z"),
    format(script_start_time, "%Y-%m-%d %H:%M:%S %Z")
  ),
  end_time = c(
    format(descriptor_end_time, "%Y-%m-%d %H:%M:%S %Z"),
    format(modeling_end_time, "%Y-%m-%d %H:%M:%S %Z"),
    format(script_end_time, "%Y-%m-%d %H:%M:%S %Z")
  ),
  wall_clock_seconds = c(
    as.numeric(difftime(descriptor_end_time, descriptor_start_time, units = "secs")),
    as.numeric(difftime(modeling_end_time, modeling_start_time, units = "secs")),
    as.numeric(difftime(script_end_time, script_start_time, units = "secs"))
  ),
  cpu_user_seconds = c(
    unname(descriptor_elapsed["user.self"]),
    unname(modeling_elapsed["user.self"]),
    unname(script_elapsed["user.self"])
  ),
  cpu_system_seconds = c(
    unname(descriptor_elapsed["sys.self"]),
    unname(modeling_elapsed["sys.self"]),
    unname(script_elapsed["sys.self"])
  ),
  cpu_elapsed_seconds = c(
    unname(descriptor_elapsed["elapsed"]),
    unname(modeling_elapsed["elapsed"]),
    unname(script_elapsed["elapsed"])
  ),
  n_repeats = c(NA_integer_, N_REPEATS, N_REPEATS),
  k_outer = c(NA_integer_, K_OUTER, K_OUTER),
  k_inner = c(NA_integer_, K_INNER, K_INNER),
  n_cores = c(N_CORES, N_CORES, N_CORES),
  n_sites = c(n, n, n),
  stringsAsFactors = FALSE
)
write.csv(timing_df, timing_file, row.names = FALSE)

cat("\n===== Repeated nested CV summary: max-inner-accuracy to 20A, 1.0A grid, parallel =====\n")
print(summary_df)

cat("\n===== EK repeated nested CV max-inner-accuracy timing summary =====\n")
print(timing_df)

cat("\nSaved outputs under:", out_dir, "\n")
