# =====================================================================
# TOUGH-C1 | Repeated nested CV with max-inner-accuracy threshold selection
#   - Binary tasks:
#       1) NUC vs CONTROL
#       2) HEME vs CONTROL
#   - Repeated stratified nested CV:
#       outer 5-fold CV repeated 5 times for performance estimation
#       inner 5-fold CV for t_end and model-parameter selection
#   - Candidate t_end values use a 1.0 A grid from 6.0 to 20.0 A.
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
rr_rds_path <- file.path("data", "cache", "TOUGH_C1_LBS_ligand_ILR_CDPA_to20A.rds")
pairs_rds <- file.path("data", "TOUGH-C1", "TOUGH_C1_protein_ligand_pairs.rds")
out_dir <- file.path("results", "TOUGH-C1")
cache_dir <- file.path("data", "cache")
checkpoint_dir <- file.path(out_dir, "checkpoints")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

# Set RR_START_FROM_SCRATCH=1 before running if you want to recompute all
# checkpoints instead of resuming completed outer folds.
START_FROM_SCRATCH <- identical(Sys.getenv("RR_START_FROM_SCRATCH"), "1")

N_CORES <- as.integer(Sys.getenv("RR_N_CORES", unset = "4"))
available_cores <- parallel::detectCores(logical = FALSE)
if (is.na(available_cores)) available_cores <- parallel::detectCores(logical = TRUE)
if (is.na(available_cores)) available_cores <- 1L
N_CORES <- max(1L, min(N_CORES, available_cores))

# ------------------- descriptor-builder helpers -------------------
read_pdb_like_minimal <- function(path) {
  ln <- tryCatch(readLines(path, warn = FALSE), error = function(e) NULL)
  if (is.null(ln)) return(NULL)
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

min_dist_to_ligand <- function(prot_coords, lig_coords) {
  nP <- nrow(prot_coords)
  nL <- nrow(lig_coords)
  dmin <- numeric(nP)

  for (i in seq_len(nP)) {
    v <- prot_coords[i, ]
    dd <- sqrt(rowSums((lig_coords - matrix(v, nL, 3, byrow = TRUE))^2))
    dmin[i] <- min(dd)
  }
  dmin
}

resolve_existing_path <- function(path_value) {
  candidates <- unique(c(path_value, file.path("data", path_value)))
  hit <- candidates[file.exists(candidates)][1]
  ifelse(is.na(hit), candidates[1], hit)
}

build_toughc1_descriptor_object <- function(save_path, thr_max0 = 20.0) {
  cat("Building TOUGH-C1 descriptor object to", thr_max0, "Angstrom...\n")

  pairs_df <- readRDS(pairs_rds)
  pairs_df$p_file <- vapply(pairs_df$p_file, resolve_existing_path, character(1))
  pairs_df$l_file <- vapply(pairs_df$l_file, resolve_existing_path, character(1))
  stopifnot(all(c("class", "p_file", "l_file") %in% names(pairs_df)))

  if (!all(file.exists(pairs_df$p_file)) || !all(file.exists(pairs_df$l_file))) {
    stop(
      "Cannot rebuild TOUGH-C1 descriptors to the requested upper threshold because some protein/ligand files are missing. ",
      "Expected folders such as TOUGH-C1_Proteins_Sorted and TOUGH-C1_Ligands_Sorted ",
      "at the project root or under data/."
    )
  }

  thr_min0 <- 3.0
  thr_step <- 0.1
  thr_full0 <- seq(thr_min0, thr_max0, by = thr_step)
  n_thr_full0 <- length(thr_full0)
  thr_max_complete <- thr_max0

  classes <- sort(unique(pairs_df$class))
  n0 <- nrow(pairs_df)
  min_atoms <- 10
  comp_names <- c("c12", "c13", "c23", "v1", "v2", "v3")

  ilr_arr0 <- array(NA_real_, dim = c(n0, n_thr_full0, 2))
  cov_arr0 <- array(NA_real_, dim = c(n0, n_thr_full0, 6))
  labels0 <- factor(pairs_df$class, levels = classes)

  for (i in seq_len(n0)) {
    p_df <- read_pdb_like_minimal(pairs_df$p_file[i])
    l_df <- read_pdb_like_minimal(pairs_df$l_file[i])
    if (is.null(p_df) || is.null(l_df)) next

    p_coords <- as.matrix(p_df[, c("x", "y", "z")])
    l_coords <- as.matrix(l_df[, c("x", "y", "z")])
    if (nrow(p_coords) == 0 || nrow(l_coords) == 0) next

    dmin_all <- min_dist_to_ligand(p_coords, l_coords)

    for (tt in seq_len(n_thr_full0)) {
      pick <- dmin_all <= thr_full0[tt]
      if (!any(pick)) next
      site <- p_df[pick, , drop = FALSE]
      if (nrow(site) < min_atoms) next

      ilr12 <- compute_ilr_noeps(site$elem)
      if (any(is.na(ilr12))) next

      cov6 <- cov_from_axis_dist(as.matrix(site[, c("x", "y", "z")]))
      ilr_arr0[i, tt, ] <- ilr12
      cov_arr0[i, tt, ] <- cov6
    }

    if (i %% 50 == 0) cat("processed", i, "of", n0, "pairs\n")
  }

  valid_radius <- array(FALSE, dim = c(n0, n_thr_full0))
  for (i in seq_len(n0)) {
    for (tt in seq_len(n_thr_full0)) {
      valid_radius[i, tt] <- all(is.finite(c(ilr_arr0[i, tt, ], cov_arr0[i, tt, ])))
    }
  }

  thr_min_grid <- thr_full0[thr_full0 <= thr_max_complete]
  frac_all <- numeric(length(thr_min_grid))
  frac_cls <- matrix(NA_real_, nrow = length(classes), ncol = length(thr_min_grid),
                     dimnames = list(classes, NULL))

  for (k in seq_along(thr_min_grid)) {
    dmin <- thr_min_grid[k]
    idx_t <- which(thr_full0 >= dmin & thr_full0 <= thr_max_complete)
    good_all <- apply(valid_radius[, idx_t, drop = FALSE], 1, all)
    frac_all[k] <- mean(good_all)
    for (j in seq_along(classes)) {
      frac_cls[j, k] <- mean(good_all[labels0 == classes[j]])
    }
  }

  thr_95_all <- thr_min_grid[which(frac_all >= 0.95)[1]]
  thr_95_each <- apply(frac_cls, 1, function(v) {
    idx <- which(v >= 0.95)[1]
    if (is.na(idx)) NA_real_ else thr_min_grid[idx]
  })

  thr_min_sel <- max(c(thr_95_all, thr_95_each), na.rm = TRUE)
  cat(sprintf("Chosen global d_min for [d_min, %.1f]: %.2f Angstrom\n", thr_max0, thr_min_sel))

  idx_use <- which(thr_full0 >= thr_min_sel & thr_full0 <= thr_max0)
  thr_full <- thr_full0[idx_use]
  ilr_arr <- ilr_arr0[, idx_use, , drop = FALSE]
  cov_arr <- cov_arr0[, idx_use, , drop = FALSE]
  labels <- labels0

  ok_row <- logical(n0)
  for (i in seq_len(n0)) {
    ok_row[i] <- all(is.finite(c(ilr_arr[i, , ], cov_arr[i, , ])))
  }

  cat("Dropping", sum(!ok_row), "pairs with missing values on selected interval\n")
  out <- list(
    thr_full = thr_full,
    ilr_arr = ilr_arr[ok_row, , , drop = FALSE],
    cov_arr = cov_arr[ok_row, , , drop = FALSE],
    labels = droplevels(labels[ok_row]),
    comp_names = comp_names,
    pairs_df = pairs_df[ok_row, , drop = FALSE],
    d_min_sel = thr_min_sel
  )

  saveRDS(out, save_path)
  cat("Saved TOUGH-C1 descriptor object:", save_path, "\n")
  out
}

# ------------------- load/build TOUGH-C1 descriptor object -------------------
target_t_end_max <- 20.0

descriptor_start_time <- Sys.time()
descriptor_start_proc <- proc.time()
if (file.exists(rr_rds_path)) {
  obj <- readRDS(rr_rds_path)
  cat("Loaded TOUGH-C1 descriptor object:", rr_rds_path, "\n")
} else {
  obj <- build_toughc1_descriptor_object(rr_rds_path, thr_max0 = target_t_end_max)
}
descriptor_end_time <- Sys.time()
descriptor_end_proc <- proc.time()
descriptor_elapsed <- descriptor_end_proc - descriptor_start_proc

thr_full   <- as.numeric(obj$thr_full)
ilr_arr    <- obj$ilr_arr
cov_arr    <- obj$cov_arr
labels     <- droplevels(obj$labels)
d_min_sel  <- obj$d_min_sel

descriptor_timing <- data.frame(
  task = "DESCRIPTOR_OBJECT_LOAD_OR_BUILD",
  start_time = format(descriptor_start_time, "%Y-%m-%d %H:%M:%S %Z"),
  end_time = format(descriptor_end_time, "%Y-%m-%d %H:%M:%S %Z"),
  wall_clock_seconds = as.numeric(difftime(descriptor_end_time, descriptor_start_time, units = "secs")),
  cpu_user_seconds = unname(descriptor_elapsed["user.self"]),
  cpu_system_seconds = unname(descriptor_elapsed["sys.self"]),
  cpu_elapsed_seconds = unname(descriptor_elapsed["elapsed"]),
  n_repeats = NA_integer_,
  k_outer = NA_integer_,
  k_inner = NA_integer_,
  n_cores = N_CORES,
  n_sites = length(labels),
  stringsAsFactors = FALSE
)

cat("Loaded/built:", rr_rds_path, "\n")
cat("Radius grid range:", range(thr_full), "Angstrom\n")
cat("d_min_sel:", d_min_sel, "Angstrom\n")
cat("Number of pockets:", length(labels), "\n")
print(table(labels))

# ------------------- analysis settings -------------------
t_start <- d_min_sel
t_end_grid <- seq(6.0, target_t_end_max, by = 1.0)

N_REPEATS <- 5
K_OUTER <- 5
K_INNER <- 5

MAX_PC     <- 20
VAR_TARGET <- 0.95

SPLINE_DEGREE <- 3
SMOOTH_LAMBDA <- 1e-3
BASIS_DIVS    <- c(2, 3, 4, 5)
NORDER        <- SPLINE_DEGREE + 1

RF_NUM_TREES <- 400
RF_MIN_NODE_SIZE_GRID <- c(2, 5, 7)

LASSO_LAMBDA_GRID <- 10^seq(-2, -4, length.out = 10)
MAH_EPS <- 1e-6

TIMING_FILE <- file.path(out_dir, "TOUGHC1_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_Resumable_timing.csv")

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

smooth_all_channels <- function(ilr_arr_task, cov_arr_task, idx_use) {
  t_vec <- thr_full[idx_use]
  nb_use <- choose_nbasis(t_vec)

  list(
    t_vec = t_vec,
    nbasis = nb_use,
    channels = list(
      ilr1 = smooth_channel(ilr_arr_task[, idx_use, 1], t_vec, nb_use, NORDER, SMOOTH_LAMBDA),
      ilr2 = smooth_channel(ilr_arr_task[, idx_use, 2], t_vec, nb_use, NORDER, SMOOTH_LAMBDA),
      c12  = smooth_channel(cov_arr_task[, idx_use, 1], t_vec, nb_use, NORDER, SMOOTH_LAMBDA),
      c13  = smooth_channel(cov_arr_task[, idx_use, 2], t_vec, nb_use, NORDER, SMOOTH_LAMBDA),
      c23  = smooth_channel(cov_arr_task[, idx_use, 3], t_vec, nb_use, NORDER, SMOOTH_LAMBDA),
      v1   = smooth_channel(cov_arr_task[, idx_use, 4], t_vec, nb_use, NORDER, SMOOTH_LAMBDA),
      v2   = smooth_channel(cov_arr_task[, idx_use, 5], t_vec, nb_use, NORDER, SMOOTH_LAMBDA),
      v3   = smooth_channel(cov_arr_task[, idx_use, 6], t_vec, nb_use, NORDER, SMOOTH_LAMBDA)
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

auc_rank <- function(y_true, score, positive_label) {
  y01 <- as.integer(y_true == positive_label)
  pos <- score[y01 == 1]
  neg <- score[y01 == 0]
  n_pos <- length(pos)
  n_neg <- length(neg)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)

  ranks <- rank(c(pos, neg), ties.method = "average")
  rank_pos <- ranks[seq_len(n_pos)]
  (sum(rank_pos) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

binary_metrics <- function(y_true, pred, score, positive_label) {
  y_true <- as.character(y_true)
  pred <- as.character(pred)
  pos <- positive_label
  neg <- setdiff(unique(y_true), pos)
  if (length(neg) != 1) stop("Expected exactly one negative class.")
  neg <- neg[1]

  tp <- sum(y_true == pos & pred == pos)
  fp <- sum(y_true == neg & pred == pos)
  tn <- sum(y_true == neg & pred == neg)
  fn <- sum(y_true == pos & pred == neg)
  safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)

  data.frame(
    accuracy = safe_div(tp + tn, tp + tn + fp + fn),
    ppv = safe_div(tp, tp + fp),
    tpr = safe_div(tp, tp + fn),
    tnr = safe_div(tn, tn + fp),
    auc = auc_rank(y_true, score, positive_label = pos),
    tp = tp,
    fp = fp,
    tn = tn,
    fn = fn,
    stringsAsFactors = FALSE
  )
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

rf_predict_metrics <- function(Ztr, y_tr, Zte, y_te, positive_label,
                               num_trees, min_node_size, seed = 1) {
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

  pp <- predict(fit, data = data.frame(Zte, check.names = FALSE))$predictions
  pred <- factor(colnames(pp)[max.col(pp)], levels = levels(y_te))
  score <- pp[, positive_label]
  binary_metrics(y_te, pred, score, positive_label = positive_label)
}

lasso_predict_accuracies <- function(Ztr, y_tr, Zte, y_te, lambda_values) {
  lambda_values <- sort(unique(lambda_values), decreasing = TRUE)
  fit_lambda_values <- sort(unique(c(1, 1e-1, lambda_values)), decreasing = TRUE)
  fit <- glmnet(
    x = as.matrix(Ztr),
    y = y_tr,
    family = "binomial",
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

lasso_predict_metrics <- function(Ztr, y_tr, Zte, y_te, positive_label, lambda_value) {
  lambda_path <- sort(unique(c(LASSO_LAMBDA_GRID, lambda_value)), decreasing = TRUE)
  fit <- glmnet(
    x = as.matrix(Ztr),
    y = y_tr,
    family = "binomial",
    alpha = 1,
    lambda = lambda_path,
    standardize = TRUE,
    control = list(maxit = 1000000)
  )

  pred <- predict(
    fit,
    newx = as.matrix(Zte),
    s = lambda_value,
    type = "class"
  )
  pred <- factor(as.vector(pred), levels = levels(y_te))

  prob <- as.vector(predict(
    fit,
    newx = as.matrix(Zte),
    s = lambda_value,
    type = "response"
  ))
  score <- if (identical(levels(y_tr)[2], positive_label)) prob else 1 - prob
  binary_metrics(y_te, pred, score, positive_label = positive_label)
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

mah_predict_metrics <- function(Ztr, y_tr, Zte, y_te, positive_label, eps_value) {
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
  score <- -Dmat[, positive_label]
  binary_metrics(y_te, pred, score, positive_label = positive_label)
}

summarize_by_grid_max_inner <- function(df, group_cols) {
  grid_df <- aggregate(
    df$accuracy,
    by = df[group_cols],
    FUN = mean
  )
  names(grid_df)[names(grid_df) == "x"] <- "mean_inner_accuracy"

  grid_df$t_end <- as.numeric(as.character(grid_df$t_end))
  if ("min.node.size" %in% names(grid_df)) {
    grid_df$min.node.size <- as.integer(as.character(grid_df$min.node.size))
  }
  if ("lambda" %in% names(grid_df)) {
    grid_df$lambda <- as.numeric(as.character(grid_df$lambda))
  }

  best_i <- which.max(grid_df$mean_inner_accuracy)
  selected <- grid_df[best_i, , drop = FALSE]
  out <- as.list(selected[group_cols])
  out$mean_inner_accuracy <- selected$mean_inner_accuracy
  out
}

checkpoint_path <- function(suffix, rep_id, outer) {
  file.path(
    checkpoint_dir,
    sprintf("TOUGHC1_MaxInner_to20A_by1p0_parallel_%s_repeat%02d_outer%02d.rds", suffix, rep_id, outer)
  )
}

run_repeated_binary_task <- function(cls_pos, cls_neg, suffix) {
  task_start_time <- Sys.time()
  task_start_proc <- proc.time()

  cat("\n============================================================\n")
  cat("Repeated nested CV task:", cls_pos, "vs", cls_neg, "\n")
  cat("============================================================\n")

  idx_task <- which(labels %in% c(cls_pos, cls_neg))
  y_task <- droplevels(labels[idx_task])
  ilr_arr_task <- ilr_arr[idx_task, , , drop = FALSE]
  cov_arr_task <- cov_arr[idx_task, , , drop = FALSE]
  n_task <- length(y_task)

  cat("Number of samples:", n_task, "\n")
  print(table(y_task))

  rf_outer_results <- list()
  lasso_outer_results <- list()
  mah_outer_results <- list()
  rf_inner_records <- list()
  lasso_inner_records <- list()
  mah_inner_records <- list()

  for (rep_id in seq_len(N_REPEATS)) {
    cat("\n############################################################\n")
    cat("Task", suffix, "| Repeat", rep_id, "of", N_REPEATS, "\n")
    cat("############################################################\n")

    outer_folds <- make_strat_folds(y_task, K_OUTER, seed = 1000 + rep_id)

    for (outer in seq_len(K_OUTER)) {
      cp_file <- checkpoint_path(suffix, rep_id, outer)
      if (!START_FROM_SCRATCH && file.exists(cp_file)) {
        cp <- readRDS(cp_file)
        required_cols <- c("outer_accuracy", "outer_ppv", "outer_tpr", "outer_tnr", "outer_auc")
        if (all(required_cols %in% names(cp$rf_outer)) &&
            all(required_cols %in% names(cp$lasso_outer)) &&
            all(required_cols %in% names(cp$mah_outer))) {
          cat(sprintf(
            "\nCheckpoint found; skipping Task %s | Repeat %d | Outer %d\n",
            suffix, rep_id, outer
          ))
          rf_outer_results[[length(rf_outer_results) + 1L]] <- cp$rf_outer
          lasso_outer_results[[length(lasso_outer_results) + 1L]] <- cp$lasso_outer
          mah_outer_results[[length(mah_outer_results) + 1L]] <- cp$mah_outer
          rf_inner_records[[length(rf_inner_records) + 1L]] <- cp$rf_inner
          lasso_inner_records[[length(lasso_inner_records) + 1L]] <- cp$lasso_inner
          mah_inner_records[[length(mah_inner_records) + 1L]] <- cp$mah_inner
          next
        }
        cat(sprintf(
          "\nOld checkpoint lacks full metrics; recomputing Task %s | Repeat %d | Outer %d\n",
          suffix, rep_id, outer
        ))
      }

      cat("\n============================================================\n")
      cat("Task", suffix, "| Repeat", rep_id, "| Outer fold", outer, "of", K_OUTER, "\n")
      cat("============================================================\n")

      outer_test <- outer_folds[[outer]]
      outer_train <- setdiff(seq_len(n_task), outer_test)
      y_outer_train <- droplevels(y_task[outer_train])
      inner_folds_rel <- make_strat_folds(y_outer_train, K_INNER, seed = 2000 + 100 * rep_id + outer)

      cat(sprintf(
        "Task %s | Repeat %d | Outer %d | evaluating %d candidate t_end values on %d cores\n",
        suffix, rep_id, outer, length(t_end_grid), N_CORES
      ))

      evaluate_t_end_candidate <- function(ii) {
        t_end <- t_end_grid[ii]
        idx_use <- which(thr_full >= t_start & thr_full <= t_end)
        if (length(idx_use) < 5) return(NULL)

        smoothed <- smooth_all_channels(ilr_arr_task, cov_arr_task, idx_use)
        rf_inner <- list()
        lasso_inner <- list()
        mah_inner <- list()

        for (inner in seq_len(K_INNER)) {
          inner_valid <- outer_train[inner_folds_rel[[inner]]]
          inner_train <- setdiff(outer_train, inner_valid)

          scores <- build_mfpca_scores(smoothed, inner_train, inner_valid)
          y_tr <- droplevels(y_task[inner_train])
          y_va <- factor(y_task[inner_valid], levels = levels(y_tr))

          for (min_node_size in RF_MIN_NODE_SIZE_GRID) {
            acc <- rf_predict_accuracy(
              scores$Ztr, y_tr, scores$Zte, y_va,
              num_trees = RF_NUM_TREES,
              min_node_size = min_node_size,
              seed = 3000 + 1000 * rep_id + 100 * outer + 10 * inner + min_node_size
            )
            rf_inner[[length(rf_inner) + 1L]] <- data.frame(
              task = suffix,
              repeat_id = rep_id,
              outer_fold = outer,
              inner_fold = inner,
              t_end = t_end,
              num.trees = RF_NUM_TREES,
              mtry_rule = "floor_sqrt_p",
              min.node.size = min_node_size,
              k_use = scores$k_use,
              mtry = max(1, floor(sqrt(ncol(scores$Ztr)))),
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
              task = suffix,
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
            task = suffix,
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

      rf_best <- summarize_by_grid_max_inner(rf_inner_df, c("t_end", "min.node.size"))
      lasso_best <- summarize_by_grid_max_inner(lasso_inner_df, c("t_end", "lambda"))
      mah_best <- summarize_by_grid_max_inner(mah_inner_df, c("t_end"))

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

      y_tr_outer <- droplevels(y_task[outer_train])
      y_te_outer <- factor(y_task[outer_test], levels = levels(y_tr_outer))

      rf_idx_use <- which(thr_full >= t_start & thr_full <= rf_best$t_end)
      rf_smoothed <- smooth_all_channels(ilr_arr_task, cov_arr_task, rf_idx_use)
      rf_scores <- build_mfpca_scores(rf_smoothed, outer_train, outer_test)
      rf_final_metrics <- rf_predict_metrics(
        rf_scores$Ztr, y_tr_outer, rf_scores$Zte, y_te_outer,
        positive_label = cls_pos,
        num_trees = RF_NUM_TREES,
        min_node_size = rf_best$min.node.size,
        seed = 9000 + 100 * rep_id + outer
      )

      lasso_idx_use <- which(thr_full >= t_start & thr_full <= lasso_best$t_end)
      lasso_smoothed <- smooth_all_channels(ilr_arr_task, cov_arr_task, lasso_idx_use)
      lasso_scores <- build_mfpca_scores(lasso_smoothed, outer_train, outer_test)
      lasso_final_metrics <- lasso_predict_metrics(
        lasso_scores$Ztr, y_tr_outer, lasso_scores$Zte, y_te_outer,
        positive_label = cls_pos,
        lambda_value = lasso_best$lambda
      )

      mah_idx_use <- which(thr_full >= t_start & thr_full <= mah_best$t_end)
      mah_smoothed <- smooth_all_channels(ilr_arr_task, cov_arr_task, mah_idx_use)
      mah_scores <- build_mfpca_scores(mah_smoothed, outer_train, outer_test)
      mah_final_metrics <- mah_predict_metrics(
        mah_scores$Ztr, y_tr_outer, mah_scores$Zte, y_te_outer,
        positive_label = cls_pos,
        eps_value = MAH_EPS
      )

      rf_outer <- data.frame(
        task = suffix,
        model = "RF",
        repeat_id = rep_id,
        outer_fold = outer,
        t_end = rf_best$t_end,
        num.trees = RF_NUM_TREES,
        mtry_rule = "floor_sqrt_p",
        min.node.size = rf_best$min.node.size,
        mean_inner_accuracy = rf_best$mean_inner_accuracy,
        outer_accuracy = rf_final_metrics$accuracy,
        outer_ppv = rf_final_metrics$ppv,
        outer_tpr = rf_final_metrics$tpr,
        outer_tnr = rf_final_metrics$tnr,
        outer_auc = rf_final_metrics$auc,
        tp = rf_final_metrics$tp,
        fp = rf_final_metrics$fp,
        tn = rf_final_metrics$tn,
        fn = rf_final_metrics$fn,
        k_use_final = rf_scores$k_use,
        mtry_final = max(1, floor(sqrt(ncol(rf_scores$Ztr)))),
        stringsAsFactors = FALSE
      )
      rf_outer_results[[length(rf_outer_results) + 1L]] <- rf_outer

      lasso_outer <- data.frame(
        task = suffix,
        model = "BinomialLasso",
        repeat_id = rep_id,
        outer_fold = outer,
        t_end = lasso_best$t_end,
        lambda = lasso_best$lambda,
        log10_lambda = log10(lasso_best$lambda),
        mean_inner_accuracy = lasso_best$mean_inner_accuracy,
        outer_accuracy = lasso_final_metrics$accuracy,
        outer_ppv = lasso_final_metrics$ppv,
        outer_tpr = lasso_final_metrics$tpr,
        outer_tnr = lasso_final_metrics$tnr,
        outer_auc = lasso_final_metrics$auc,
        tp = lasso_final_metrics$tp,
        fp = lasso_final_metrics$fp,
        tn = lasso_final_metrics$tn,
        fn = lasso_final_metrics$fn,
        k_use_final = lasso_scores$k_use,
        stringsAsFactors = FALSE
      )
      lasso_outer_results[[length(lasso_outer_results) + 1L]] <- lasso_outer

      mah_outer <- data.frame(
        task = suffix,
        model = "Mahalanobis",
        repeat_id = rep_id,
        outer_fold = outer,
        t_end = mah_best$t_end,
        eps = MAH_EPS,
        log10_eps = log10(MAH_EPS),
        mean_inner_accuracy = mah_best$mean_inner_accuracy,
        outer_accuracy = mah_final_metrics$accuracy,
        outer_ppv = mah_final_metrics$ppv,
        outer_tpr = mah_final_metrics$tpr,
        outer_tnr = mah_final_metrics$tnr,
        outer_auc = mah_final_metrics$auc,
        tp = mah_final_metrics$tp,
        fp = mah_final_metrics$fp,
        tn = mah_final_metrics$tn,
        fn = mah_final_metrics$fn,
        k_use_final = mah_scores$k_use,
        stringsAsFactors = FALSE
      )
      mah_outer_results[[length(mah_outer_results) + 1L]] <- mah_outer

      cat(sprintf(
        "\nTask %s | Repeat %d | Outer %d final metrics | RF Acc=%.3f AUC=%.3f | Lasso Acc=%.3f AUC=%.3f | Mahalanobis Acc=%.3f AUC=%.3f\n",
        suffix, rep_id, outer,
        rf_final_metrics$accuracy, rf_final_metrics$auc,
        lasso_final_metrics$accuracy, lasso_final_metrics$auc,
        mah_final_metrics$accuracy, mah_final_metrics$auc
      ))

      saveRDS(
        list(
          task = suffix,
          repeat_id = rep_id,
          outer_fold = outer,
          selection_rule = "max_inner_accuracy",
          rf_outer = rf_outer,
          lasso_outer = lasso_outer,
          mah_outer = mah_outer,
          rf_inner = rf_inner_df,
          lasso_inner = lasso_inner_df,
          mah_inner = mah_inner_df
        ),
        cp_file
      )
      cat("Saved checkpoint:", cp_file, "\n")

      rm(rf_smoothed, lasso_smoothed, mah_smoothed, rf_scores, lasso_scores, mah_scores,
         rf_final_metrics, lasso_final_metrics, mah_final_metrics)
      invisible(gc())
    }
  }

  rf_outer_df <- do.call(rbind, rf_outer_results)
  lasso_outer_df <- do.call(rbind, lasso_outer_results)
  mah_outer_df <- do.call(rbind, mah_outer_results)
  rf_inner_df_all <- do.call(rbind, rf_inner_records)
  lasso_inner_df_all <- do.call(rbind, lasso_inner_records)
  mah_inner_df_all <- do.call(rbind, mah_inner_records)

  write.csv(rf_outer_df, file.path(out_dir, sprintf("TOUGHC1_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_Resumable_%s_RF_outer_results.csv", suffix)), row.names = FALSE)
  write.csv(lasso_outer_df, file.path(out_dir, sprintf("TOUGHC1_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_Resumable_%s_BinomialLasso_outer_results.csv", suffix)), row.names = FALSE)
  write.csv(mah_outer_df, file.path(out_dir, sprintf("TOUGHC1_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_Resumable_%s_Mahalanobis_outer_results.csv", suffix)), row.names = FALSE)
  write.csv(rf_inner_df_all, file.path(out_dir, sprintf("TOUGHC1_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_Resumable_%s_RF_inner_results.csv", suffix)), row.names = FALSE)
  write.csv(lasso_inner_df_all, file.path(out_dir, sprintf("TOUGHC1_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_Resumable_%s_BinomialLasso_inner_results.csv", suffix)), row.names = FALSE)
  write.csv(mah_inner_df_all, file.path(out_dir, sprintf("TOUGHC1_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_Resumable_%s_Mahalanobis_inner_results.csv", suffix)), row.names = FALSE)

  summary_df <- rbind(
    data.frame(
      task = suffix,
      model = "RF",
      selection_rule = "max_inner_accuracy",
      n_repeats = N_REPEATS,
      n_outer_evals = nrow(rf_outer_df),
      mean_outer_accuracy = mean(rf_outer_df$outer_accuracy),
      sd_outer_accuracy = sd(rf_outer_df$outer_accuracy),
      mean_outer_ppv = mean(rf_outer_df$outer_ppv),
      sd_outer_ppv = sd(rf_outer_df$outer_ppv),
      mean_outer_tpr = mean(rf_outer_df$outer_tpr),
      sd_outer_tpr = sd(rf_outer_df$outer_tpr),
      mean_outer_tnr = mean(rf_outer_df$outer_tnr),
      sd_outer_tnr = sd(rf_outer_df$outer_tnr),
      mean_outer_auc = mean(rf_outer_df$outer_auc),
      sd_outer_auc = sd(rf_outer_df$outer_auc),
      mean_selected_t_end = mean(rf_outer_df$t_end),
      sd_selected_t_end = sd(rf_outer_df$t_end),
      mean_selected_log10_lambda = NA_real_,
      sd_selected_log10_lambda = NA_real_,
      mean_selected_log10_eps = NA_real_,
      sd_selected_log10_eps = NA_real_,
      stringsAsFactors = FALSE
    ),
    data.frame(
      task = suffix,
      model = "BinomialLasso",
      selection_rule = "max_inner_accuracy",
      n_repeats = N_REPEATS,
      n_outer_evals = nrow(lasso_outer_df),
      mean_outer_accuracy = mean(lasso_outer_df$outer_accuracy),
      sd_outer_accuracy = sd(lasso_outer_df$outer_accuracy),
      mean_outer_ppv = mean(lasso_outer_df$outer_ppv),
      sd_outer_ppv = sd(lasso_outer_df$outer_ppv),
      mean_outer_tpr = mean(lasso_outer_df$outer_tpr),
      sd_outer_tpr = sd(lasso_outer_df$outer_tpr),
      mean_outer_tnr = mean(lasso_outer_df$outer_tnr),
      sd_outer_tnr = sd(lasso_outer_df$outer_tnr),
      mean_outer_auc = mean(lasso_outer_df$outer_auc),
      sd_outer_auc = sd(lasso_outer_df$outer_auc),
      mean_selected_t_end = mean(lasso_outer_df$t_end),
      sd_selected_t_end = sd(lasso_outer_df$t_end),
      mean_selected_log10_lambda = mean(lasso_outer_df$log10_lambda),
      sd_selected_log10_lambda = sd(lasso_outer_df$log10_lambda),
      mean_selected_log10_eps = NA_real_,
      sd_selected_log10_eps = NA_real_,
      stringsAsFactors = FALSE
    ),
    data.frame(
      task = suffix,
      model = "Mahalanobis",
      selection_rule = "max_inner_accuracy",
      n_repeats = N_REPEATS,
      n_outer_evals = nrow(mah_outer_df),
      mean_outer_accuracy = mean(mah_outer_df$outer_accuracy),
      sd_outer_accuracy = sd(mah_outer_df$outer_accuracy),
      mean_outer_ppv = mean(mah_outer_df$outer_ppv),
      sd_outer_ppv = sd(mah_outer_df$outer_ppv),
      mean_outer_tpr = mean(mah_outer_df$outer_tpr),
      sd_outer_tpr = sd(mah_outer_df$outer_tpr),
      mean_outer_tnr = mean(mah_outer_df$outer_tnr),
      sd_outer_tnr = sd(mah_outer_df$outer_tnr),
      mean_outer_auc = mean(mah_outer_df$outer_auc),
      sd_outer_auc = sd(mah_outer_df$outer_auc),
      mean_selected_t_end = mean(mah_outer_df$t_end),
      sd_selected_t_end = sd(mah_outer_df$t_end),
      mean_selected_log10_lambda = NA_real_,
      sd_selected_log10_lambda = NA_real_,
      mean_selected_log10_eps = mean(mah_outer_df$log10_eps),
      sd_selected_log10_eps = sd(mah_outer_df$log10_eps),
      stringsAsFactors = FALSE
    )
  )

  write.csv(summary_df, file.path(out_dir, sprintf("TOUGHC1_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_Resumable_%s_summary.csv", suffix)), row.names = FALSE)

  cat("\n===== Repeated nested CV max-inner-accuracy resumable summary, 1.0A grid, parallel:", suffix, "=====\n")
  print(summary_df)

  task_end_time <- Sys.time()
  task_end_proc <- proc.time()
  task_elapsed <- task_end_proc - task_start_proc
  timing_df <- data.frame(
    task = suffix,
    start_time = format(task_start_time, "%Y-%m-%d %H:%M:%S %Z"),
    end_time = format(task_end_time, "%Y-%m-%d %H:%M:%S %Z"),
    wall_clock_seconds = as.numeric(difftime(task_end_time, task_start_time, units = "secs")),
    cpu_user_seconds = unname(task_elapsed["user.self"]),
    cpu_system_seconds = unname(task_elapsed["sys.self"]),
    cpu_elapsed_seconds = unname(task_elapsed["elapsed"]),
    n_repeats = N_REPEATS,
    k_outer = K_OUTER,
    k_inner = K_INNER,
    n_cores = N_CORES,
    n_sites = n_task,
    stringsAsFactors = FALSE
  )

  invisible(list(summary = summary_df, timing = timing_df))
}

summary_NUC <- run_repeated_binary_task("NUC", "CONTROL", "NUC_vs_CONTROL")
summary_HEME <- run_repeated_binary_task("HEME", "CONTROL", "HEME_vs_CONTROL")

summary_all <- rbind(summary_NUC$summary, summary_HEME$summary)
write.csv(summary_all, file.path(out_dir, "TOUGHC1_repeatedNestedCV_MaxInner_to20A_by1p0_parallel_Resumable_all_tasks_summary.csv"), row.names = FALSE)

script_end_time <- Sys.time()
script_end_proc <- proc.time()
script_elapsed <- script_end_proc - script_start_proc
timing_all <- rbind(
  descriptor_timing,
  summary_NUC$timing,
  summary_HEME$timing,
  data.frame(
    task = "ALL",
    start_time = format(script_start_time, "%Y-%m-%d %H:%M:%S %Z"),
    end_time = format(script_end_time, "%Y-%m-%d %H:%M:%S %Z"),
    wall_clock_seconds = as.numeric(difftime(script_end_time, script_start_time, units = "secs")),
    cpu_user_seconds = unname(script_elapsed["user.self"]),
    cpu_system_seconds = unname(script_elapsed["sys.self"]),
    cpu_elapsed_seconds = unname(script_elapsed["elapsed"]),
    n_repeats = N_REPEATS,
    k_outer = K_OUTER,
    k_inner = K_INNER,
    n_cores = N_CORES,
    n_sites = length(labels),
    stringsAsFactors = FALSE
  )
)
write.csv(timing_all, TIMING_FILE, row.names = FALSE)

cat("\n===== TOUGH-C1 repeated nested CV max-inner-accuracy resumable all-task summary, 1.0A grid, parallel =====\n")
print(summary_all)

cat("\n===== TOUGH-C1 repeated nested CV max-inner-accuracy timing summary =====\n")
print(timing_all)

cat("\nSaved TOUGH-C1 repeated nested CV max-inner-accuracy resumable outputs under:", out_dir, "\n")
