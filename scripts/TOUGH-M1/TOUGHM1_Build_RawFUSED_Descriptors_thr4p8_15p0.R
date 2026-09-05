# =====================================================================
# TOUGH-M1 | FUSED/MFPC pairwise pocket-matching AUC
#   - Uses TOUGH-M1 official positive/negative pair lists.
#   - Fixed interval: t_start = 4.8 A, t_end = 15.0 A.
# =====================================================================

rm(list = ls()); invisible(gc())
set.seed(1)
options(repos = c(CRAN = "https://cran.rstudio.com/"))

suppressPackageStartupMessages({
  library(fda)
})

# ------------------- paths -------------------
inspect_dir <- file.path("data", "TOUGH-M1")
archive_file <- file.path(inspect_dir, "TOUGH-M1_dataset.tar.gz")
dataset_dir <- file.path(inspect_dir, "TOUGH-M1_dataset")

positive_file <- file.path(inspect_dir, "TOUGH-M1_positive.list")
negative_file <- file.path(inspect_dir, "TOUGH-M1_negative.list")
pocket_file <- file.path(inspect_dir, "TOUGH-M1_pocket.list")

cache_dir <- file.path("data", "cache")
out_dir <- file.path("results", "TOUGH-M1")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw_cache_file <- file.path(cache_dir, "TOUGHM1_raw_ILR_CDPA_thr4p8_15p0_by0p1.rds")
raw_progress_file <- file.path(cache_dir, "TOUGHM1_raw_ILR_CDPA_thr4p8_15p0_by0p1_PROGRESS.rds")
score_cache_file <- file.path(cache_dir, "TOUGHM1_MFPC_scores_thr4p8_15p0_var95.rds")

pair_scores_file <- file.path(out_dir, "TOUGHM1_FUSED_pair_scores_thr4p8_15p0.csv")
summary_file <- file.path(out_dir, "TOUGHM1_FUSED_pairwise_AUC_thr4p8_15p0_summary.csv")
score_summary_file <- file.path(out_dir, "TOUGHM1_FUSED_pair_score_distribution_thr4p8_15p0.csv")
timing_file <- file.path(out_dir, "TOUGHM1_FUSED_pairwise_AUC_thr4p8_15p0_timing.csv")

# ------------------- analysis settings -------------------
t_start <- 4.8
t_end <- 15.0
t_step <- 0.1
thr_full <- seq(t_start, t_end, by = t_step)

MIN_SITE_ATOMS <- 10
MAX_PC <- 20
VAR_TARGET <- 0.95

SPLINE_DEGREE <- 3
SMOOTH_LAMBDA <- 1e-3
BASIS_DIVS <- c(2, 3, 4, 5)
NORDER <- SPLINE_DEGREE + 1

SAVE_EVERY <- 100
PAIR_CHUNK_SIZE <- 100000

# =====================================================================
# helpers
# =====================================================================
ensure_toughm1_extracted <- function() {
  if (dir.exists(dataset_dir)) {
    n_dirs <- length(list.dirs(dataset_dir, recursive = FALSE, full.names = FALSE))
    if (n_dirs >= 7000) {
      cat("Using extracted TOUGH-M1 dataset:", dataset_dir, "\n")
      return(invisible(TRUE))
    }
  }

  if (!file.exists(archive_file)) {
    stop("Missing TOUGH-M1 archive: ", archive_file)
  }

  cat("Extracting TOUGH-M1 archive to", inspect_dir, "\n")
  utils::untar(archive_file, exdir = inspect_dir)
  invisible(TRUE)
}

read_pdb_like_minimal <- function(path) {
  ln <- tryCatch(readLines(path, warn = FALSE), error = function(e) NULL)
  if (is.null(ln)) return(NULL)
  ln <- ln[grepl("^(ATOM  |HETATM)", ln)]
  if (length(ln) == 0) return(NULL)

  x <- suppressWarnings(as.numeric(substr(ln, 31, 38)))
  y <- suppressWarnings(as.numeric(substr(ln, 39, 46)))
  z <- suppressWarnings(as.numeric(substr(ln, 47, 54)))

  elem_raw <- trimws(substr(ln, 77, 78))
  aname <- trimws(substr(ln, 13, 16))
  elem <- ifelse(elem_raw == "", toupper(substr(aname, 1, 1)), toupper(elem_raw))

  keep <- is.finite(x) & is.finite(y) & is.finite(z) & elem != "H"
  data.frame(
    x = x[keep], y = y[keep], z = z[keep], elem = elem[keep],
    stringsAsFactors = FALSE
  )
}

compute_ilr_noeps <- function(elems) {
  cntC <- sum(elems == "C")
  cntO <- sum(elems == "O")
  cntN <- sum(elems == "N")
  tot <- cntC + cntO + cntN
  if (tot == 0) return(c(NA_real_, NA_real_))

  xC <- cntC / tot
  xO <- cntO / tot
  xN <- cntN / tot
  if (xC == 0 || xO == 0 || xN == 0) return(c(NA_real_, NA_real_))

  c(
    sqrt(1 / 2) * log(xC / xO),
    sqrt(2 / 3) * log(sqrt(xC * xO) / xN)
  )
}

distances_to_local_axes <- function(coords_mat) {
  if (nrow(coords_mat) < 3) {
    return(matrix(NA_real_, nrow(coords_mat), 3))
  }

  center <- colMeans(coords_mat)
  X <- sweep(coords_mat, 2, center, "-")
  cv <- cov(X)
  eg <- tryCatch(eigen(cv, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eg)) return(matrix(NA_real_, nrow(coords_mat), 3))

  axes <- eg$vectors[, 1:3, drop = FALSE]
  dmat <- matrix(NA_real_, nrow(coords_mat), 3)

  for (k in 1:3) {
    a <- axes[, k]
    proj <- as.numeric(X %*% a)
    resid <- X - tcrossprod(proj, a)
    dmat[, k] <- sqrt(rowSums(resid * resid))
  }

  dmat
}

cov_from_axis_dist <- function(coords_mat) {
  if (nrow(coords_mat) < 3) return(rep(NA_real_, 6))

  dmat <- distances_to_local_axes(coords_mat)
  if (!all(is.finite(dmat))) return(rep(NA_real_, 6))

  cv <- cov(dmat)
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
  d2min <- rep(Inf, nrow(prot_coords))
  for (j in seq_len(nrow(lig_coords))) {
    dx <- prot_coords[, 1] - lig_coords[j, 1]
    dy <- prot_coords[, 2] - lig_coords[j, 2]
    dz <- prot_coords[, 3] - lig_coords[j, 3]
    d2min <- pmin(d2min, dx * dx + dy * dy + dz * dz)
  }
  sqrt(d2min)
}

load_pocket_metadata <- function() {
  if (!file.exists(pocket_file)) stop("Missing pocket list: ", pocket_file)

  tab <- read.table(
    pocket_file, header = FALSE, stringsAsFactors = FALSE,
    col.names = c("code5", "selected_fpocket_number", "overlap_score")
  )

  tab$protein_file <- file.path(dataset_dir, tab$code5, paste0(tab$code5, ".pdb"))
  tab$ligand_file <- file.path(dataset_dir, tab$code5, paste0(tab$code5, "00.pdb"))
  tab
}

build_or_load_raw_descriptors <- function() {
  if (file.exists(raw_cache_file)) {
    cat("Loading raw TOUGH-M1 descriptor cache:", raw_cache_file, "\n")
    return(readRDS(raw_cache_file))
  }

  ensure_toughm1_extracted()
  meta <- load_pocket_metadata()
  n <- nrow(meta)
  n_thr <- length(thr_full)

  if (file.exists(raw_progress_file)) {
    cat("Resuming raw descriptor progress:", raw_progress_file, "\n")
    prog <- readRDS(raw_progress_file)
    ilr_arr <- prog$ilr_arr
    cov_arr <- prog$cov_arr
    processed <- prog$processed
    meta <- prog$metadata
  } else {
    cat("Starting raw descriptor cache for", n, "TOUGH-M1 entries.\n")
    ilr_arr <- array(NA_real_, dim = c(n, n_thr, 2))
    cov_arr <- array(NA_real_, dim = c(n, n_thr, 6))
    processed <- rep(FALSE, n)
  }

  for (i in seq_len(n)) {
    if (processed[i]) next

    p_df <- read_pdb_like_minimal(meta$protein_file[i])
    l_df <- read_pdb_like_minimal(meta$ligand_file[i])

    if (!is.null(p_df) && !is.null(l_df) && nrow(p_df) > 0 && nrow(l_df) > 0) {
      p_coords <- as.matrix(p_df[, c("x", "y", "z")])
      l_coords <- as.matrix(l_df[, c("x", "y", "z")])
      dmin_all <- min_dist_to_ligand(p_coords, l_coords)

      for (tt in seq_len(n_thr)) {
        pick <- dmin_all <= thr_full[tt]
        if (sum(pick) < MIN_SITE_ATOMS) next

        site <- p_df[pick, , drop = FALSE]
        ilr12 <- compute_ilr_noeps(site$elem)
        if (any(!is.finite(ilr12))) next

        cov6 <- cov_from_axis_dist(as.matrix(site[, c("x", "y", "z")]))
        if (any(!is.finite(cov6))) next

        ilr_arr[i, tt, ] <- ilr12
        cov_arr[i, tt, ] <- cov6
      }
    }

    processed[i] <- TRUE

    if (i %% SAVE_EVERY == 0 || i == n) {
      saveRDS(
        list(
          thr_full = thr_full,
          ilr_arr = ilr_arr,
          cov_arr = cov_arr,
          processed = processed,
          metadata = meta
        ),
        raw_progress_file
      )
      cat("processed", i, "of", n, "\n")
    }
  }

  ok_row <- logical(n)
  for (i in seq_len(n)) {
    ok_row[i] <- all(is.finite(c(ilr_arr[i, , ], cov_arr[i, , ])))
  }

  out <- list(
    thr_full = thr_full,
    ilr_arr = ilr_arr[ok_row, , , drop = FALSE],
    cov_arr = cov_arr[ok_row, , , drop = FALSE],
    metadata = meta[ok_row, , drop = FALSE],
    dropped_metadata = meta[!ok_row, , drop = FALSE],
    comp_names = c("c12", "c13", "c23", "v1", "v2", "v3"),
    t_start = t_start,
    t_end = t_end,
    min_site_atoms = MIN_SITE_ATOMS
  )

  saveRDS(out, raw_cache_file)
  cat("Saved raw descriptor cache:", raw_cache_file, "\n")
  cat("Used entries:", nrow(out$metadata), "| Dropped:", nrow(out$dropped_metadata), "\n")
  out
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

smooth_all_channels <- function(ilr_arr, cov_arr, tvec) {
  nb_use <- choose_nbasis(tvec)
  cat("Smoothing with nbasis =", nb_use, "and lambda =", SMOOTH_LAMBDA, "\n")

  list(
    t_vec = tvec,
    nbasis = nb_use,
    channels = list(
      ilr1 = smooth_channel(ilr_arr[, , 1], tvec, nb_use, NORDER, SMOOTH_LAMBDA),
      ilr2 = smooth_channel(ilr_arr[, , 2], tvec, nb_use, NORDER, SMOOTH_LAMBDA),
      c12 = smooth_channel(cov_arr[, , 1], tvec, nb_use, NORDER, SMOOTH_LAMBDA),
      c13 = smooth_channel(cov_arr[, , 2], tvec, nb_use, NORDER, SMOOTH_LAMBDA),
      c23 = smooth_channel(cov_arr[, , 3], tvec, nb_use, NORDER, SMOOTH_LAMBDA),
      v1 = smooth_channel(cov_arr[, , 4], tvec, nb_use, NORDER, SMOOTH_LAMBDA),
      v2 = smooth_channel(cov_arr[, , 5], tvec, nb_use, NORDER, SMOOTH_LAMBDA),
      v3 = smooth_channel(cov_arr[, , 6], tvec, nb_use, NORDER, SMOOTH_LAMBDA)
    )
  )
}

build_or_load_mfpc_scores <- function(raw) {
  if (file.exists(score_cache_file)) {
    cat("Loading MFPC score cache:", score_cache_file, "\n")
    return(readRDS(score_cache_file))
  }

  cat("Building global MFPC scores for TOUGH-M1.\n")
  sm <- smooth_all_channels(raw$ilr_arr, raw$cov_arr, raw$thr_full)
  ch <- sm$channels

  ilr_energy <- mean(abs(c(ch$ilr1, ch$ilr2)))
  cdpa_energy <- mean(abs(c(ch$c12, ch$c13, ch$c23, ch$v1, ch$v2, ch$v3)))
  scale_ilr <- if (ilr_energy == 0) 1 else cdpa_energy / ilr_energy

  X <- cbind(
    ch$ilr1 * scale_ilr,
    ch$ilr2 * scale_ilr,
    ch$c12, ch$c13, ch$c23,
    ch$v1, ch$v2, ch$v3
  )

  if (!all(is.finite(X))) stop("Non-finite values before PCA.")

  center <- colMeans(X)
  Xc <- sweep(X, 2, center, "-")
  eg <- eigen(cov(Xc), symmetric = TRUE)
  keep <- which(eg$values > .Machine$double.eps)
  values <- pmax(eg$values[keep], 0)
  rotation <- eg$vectors[, keep, drop = FALSE]
  prop_var <- values / sum(values)
  cum_var <- cumsum(prop_var)
  k_use <- which(cum_var >= VAR_TARGET)[1]
  if (is.na(k_use)) k_use <- length(cum_var)
  k_use <- min(k_use, MAX_PC)

  Z <- Xc %*% rotation[, seq_len(k_use), drop = FALSE]
  colnames(Z) <- paste0("PC", seq_len(k_use))
  rownames(Z) <- raw$metadata$code5

  out <- list(
    scores = Z,
    metadata = raw$metadata,
    k_use = k_use,
    var_target = VAR_TARGET,
    explained_variance = sum(prop_var[seq_len(k_use)]),
    scale_ilr = scale_ilr,
    center = center,
    rotation = rotation[, seq_len(k_use), drop = FALSE],
    eigenvalues = values
  )

  saveRDS(out, score_cache_file)
  cat("Saved MFPC score cache:", score_cache_file, "\n")
  cat("Retained MFPC scores:", k_use, "| explained variance:", out$explained_variance, "\n")
  out
}

read_pair_file <- function(path, label_value) {
  tab <- read.table(path, header = FALSE, stringsAsFactors = FALSE)
  data.frame(
    id1 = tab[[1]],
    id2 = tab[[2]],
    benchmark_score = tab[[5]],
    label = label_value,
    stringsAsFactors = FALSE
  )
}

auc_rank <- function(labels, scores) {
  labels <- as.integer(labels)
  ok <- is.finite(scores) & labels %in% c(0L, 1L)
  labels <- labels[ok]
  scores <- scores[ok]

  n_pos <- as.numeric(sum(labels == 1L))
  n_neg <- as.numeric(sum(labels == 0L))
  if (n_pos == 0 || n_neg == 0) return(NA_real_)

  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[labels == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

score_official_pairs <- function(score_obj) {
  pos <- read_pair_file(positive_file, 1L)
  neg <- read_pair_file(negative_file, 0L)
  pairs <- rbind(pos, neg)

  id_to_row <- setNames(seq_len(nrow(score_obj$scores)), rownames(score_obj$scores))
  idx1 <- unname(id_to_row[pairs$id1])
  idx2 <- unname(id_to_row[pairs$id2])
  ok <- is.finite(idx1) & is.finite(idx2)

  if (sum(!ok) > 0) {
    cat("Dropping", sum(!ok), "official pairs with missing FUSED scores.\n")
  }

  pairs <- pairs[ok, , drop = FALSE]
  idx1 <- idx1[ok]
  idx2 <- idx2[ok]

  n_pairs <- nrow(pairs)
  fused_distance <- numeric(n_pairs)

  for (start in seq(1, n_pairs, by = PAIR_CHUNK_SIZE)) {
    end <- min(start + PAIR_CHUNK_SIZE - 1, n_pairs)
    rows <- start:end
    diff <- score_obj$scores[idx1[rows], , drop = FALSE] -
      score_obj$scores[idx2[rows], , drop = FALSE]
    fused_distance[rows] <- sqrt(rowSums(diff * diff))
    cat("scored pairs", start, "to", end, "of", n_pairs, "\n")
  }

  pairs$fused_distance <- fused_distance
  pairs$fused_similarity <- -fused_distance
  pairs
}

summarize_pair_scores <- function(pair_scores, score_obj, raw) {
  auc <- auc_rank(pair_scores$label, pair_scores$fused_similarity)

  summary <- data.frame(
    dataset = "TOUGH-M1",
    method = "FUSED_MFPC_Euclidean",
    t_start = t_start,
    t_end = t_end,
    n_structures_in_pocket_list = length(readLines(pocket_file, warn = FALSE)),
    n_structures_with_fused_scores = nrow(score_obj$scores),
    n_structures_dropped = nrow(raw$dropped_metadata),
    n_positive_pairs_used = sum(pair_scores$label == 1L),
    n_negative_pairs_used = sum(pair_scores$label == 0L),
    retained_mfpc_scores = score_obj$k_use,
    explained_variance = score_obj$explained_variance,
    roc_auc = auc,
    stringsAsFactors = FALSE
  )

  score_summary <- do.call(
    rbind,
    lapply(split(pair_scores, pair_scores$label), function(df) {
      data.frame(
        label = unique(df$label),
        pair_type = ifelse(unique(df$label) == 1L, "positive", "negative"),
        n_pairs = nrow(df),
        mean_distance = mean(df$fused_distance),
        sd_distance = sd(df$fused_distance),
        median_distance = median(df$fused_distance),
        q25_distance = unname(quantile(df$fused_distance, 0.25)),
        q75_distance = unname(quantile(df$fused_distance, 0.75)),
        mean_similarity = mean(df$fused_similarity),
        sd_similarity = sd(df$fused_similarity),
        stringsAsFactors = FALSE
      )
    })
  )

  list(summary = summary, score_summary = score_summary)
}

# =====================================================================
# run analysis
# =====================================================================
stopifnot(file.exists(positive_file), file.exists(negative_file), file.exists(pocket_file))

script_start_wall <- Sys.time()
script_start_cpu <- proc.time()

descriptor_start_wall <- Sys.time()
descriptor_start_cpu <- proc.time()
raw <- build_or_load_raw_descriptors()
descriptor_elapsed_cpu <- proc.time() - descriptor_start_cpu
descriptor_end_wall <- Sys.time()

mfpc_start_wall <- Sys.time()
mfpc_start_cpu <- proc.time()
score_obj <- build_or_load_mfpc_scores(raw)
mfpc_elapsed_cpu <- proc.time() - mfpc_start_cpu
mfpc_end_wall <- Sys.time()

pair_start_wall <- Sys.time()
pair_start_cpu <- proc.time()
pair_scores <- score_official_pairs(score_obj)
pair_elapsed_cpu <- proc.time() - pair_start_cpu
pair_end_wall <- Sys.time()

summ <- summarize_pair_scores(pair_scores, score_obj, raw)

write.csv(pair_scores, pair_scores_file, row.names = FALSE)
write.csv(summ$summary, summary_file, row.names = FALSE)
write.csv(summ$score_summary, score_summary_file, row.names = FALSE)

script_elapsed_cpu <- proc.time() - script_start_cpu
script_end_wall <- Sys.time()
timing_summary <- data.frame(
  step = c("DESCRIPTOR_CONSTRUCTION", "MFPC_CONSTRUCTION", "PAIR_SCORING", "ALL"),
  start_time = c(descriptor_start_wall, mfpc_start_wall, pair_start_wall, script_start_wall),
  end_time = c(descriptor_end_wall, mfpc_end_wall, pair_end_wall, script_end_wall),
  wall_clock_seconds = as.numeric(difftime(
    c(descriptor_end_wall, mfpc_end_wall, pair_end_wall, script_end_wall),
    c(descriptor_start_wall, mfpc_start_wall, pair_start_wall, script_start_wall),
    units = "secs"
  )),
  cpu_user_seconds = c(
    descriptor_elapsed_cpu[["user.self"]],
    mfpc_elapsed_cpu[["user.self"]],
    pair_elapsed_cpu[["user.self"]],
    script_elapsed_cpu[["user.self"]]
  ),
  cpu_system_seconds = c(
    descriptor_elapsed_cpu[["sys.self"]],
    mfpc_elapsed_cpu[["sys.self"]],
    pair_elapsed_cpu[["sys.self"]],
    script_elapsed_cpu[["sys.self"]]
  ),
  n_pockets_used = c(nrow(raw$metadata), nrow(raw$metadata), nrow(raw$metadata), nrow(raw$metadata)),
  n_pockets_dropped = c(nrow(raw$dropped_metadata), nrow(raw$dropped_metadata), nrow(raw$dropped_metadata), nrow(raw$dropped_metadata)),
  stringsAsFactors = FALSE
)
write.csv(timing_summary, timing_file, row.names = FALSE)

cat("\n===== TOUGH-M1 FUSED pairwise AUC summary =====\n")
print(summ$summary)

cat("\n===== TOUGH-M1 FUSED pair score distribution =====\n")
print(summ$score_summary)

cat("\n===== TOUGH-M1 FUSED 4.8-15.0 timing summary =====\n")
print(timing_summary)

cat("\nSaved outputs:\n")
cat("  ", pair_scores_file, "\n")
cat("  ", summary_file, "\n")
cat("  ", score_summary_file, "\n")
cat("  ", timing_file, "\n")
