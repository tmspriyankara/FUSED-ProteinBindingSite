# =====================================================================
# Extended Kahraman | Group mean functional curves with raw ILR values
#   - Smooths ILR and CDPA channels over 4.8--20.0 Angstrom.
#   - Plots group mean curves.
# =====================================================================

rm(list = ls()); invisible(gc())
set.seed(1)
options(repos = c(CRAN = "https://cran.rstudio.com/"))

suppressPackageStartupMessages({
  library(fda)
})

cache_file <- file.path(
  "data", "cache",
  "EK_raw_ILR_CDPA_thr4p8_20p0_by0p1.rds"
)
out_dir <- file.path("results", "EK")
fig_dir <- file.path("figures", "EK")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cache_file)) {
  stop("Missing EK raw descriptor cache: ", cache_file)
}

T_START <- 4.8
T_END <- 20.0
T_STEP <- 0.1
SPLINE_DEGREE <- 3
SMOOTH_LAMBDA <- 1e-3
BASIS_DIVS <- c(2, 3, 4, 5)
NORDER <- SPLINE_DEGREE + 1

target_classes <- c("AMP", "ATP", "FAD", "FMN", "GLC", "HEM", "NAD", "PO4")
chan_names <- c("ILR1", "ILR2", "c12", "c13", "c23", "v1", "v2", "v3")

choose_nbasis <- function(tvec) {
  n_t <- length(tvec)
  cand_nb <- sort(unique(pmax(NORDER, round(n_t / BASIS_DIVS))))
  cand_nb[1]
}

smooth_channel <- function(mat, tvec, nbasis, norder, lambda) {
  basis <- create.bspline.basis(range(tvec), nbasis, norder)
  fd_param <- fdPar(basis, int2Lfd(2), lambda)
  fit <- smooth.basis(tvec, t(mat), fd_param)
  t(eval.fd(tvec, fit$fd))
}

lab_fun <- function(nm) {
  switch(
    nm,
    "ILR1" = expression(ILR[1]),
    "ILR2" = expression(ILR[2]),
    "c12"  = expression(c[12]),
    "c13"  = expression(c[13]),
    "c23"  = expression(c[23]),
    "v1"   = expression(v[1]),
    "v2"   = expression(v[2]),
    "v3"   = expression(v[3]),
    nm
  )
}

obj <- readRDS(cache_file)
thr_full <- seq(T_START, by = T_STEP, length.out = dim(obj$ilr_arr)[2])
idx_use <- which(thr_full >= T_START & thr_full <= T_END)
thr_use <- thr_full[idx_use]
nbasis <- choose_nbasis(thr_use)

labels <- factor(as.character(obj$labels), levels = target_classes)
classes <- levels(droplevels(labels))
n <- length(labels)

cat("\n===== EK group mean curves using raw ILR + CDPA =====\n")
cat("Sites:", n, "\n")
print(table(labels))
cat("Threshold interval:", min(thr_use), "to", max(thr_use), "\n")
cat("Threshold grid points:", length(thr_use), "\n")
cat("B-spline nbasis:", nbasis, "\n")
cat("Smoothing lambda:", SMOOTH_LAMBDA, "\n")
cat("Important: ILR group mean curves are plotted without ILR/CDPA scaling.\n\n")

ilr1_s <- smooth_channel(obj$ilr_arr[, idx_use, 1], thr_use, nbasis, NORDER, SMOOTH_LAMBDA)
ilr2_s <- smooth_channel(obj$ilr_arr[, idx_use, 2], thr_use, nbasis, NORDER, SMOOTH_LAMBDA)

cov_s <- vector("list", 6)
for (j in seq_len(6)) {
  cov_s[[j]] <- smooth_channel(obj$cov_arr[, idx_use, j], thr_use, nbasis, NORDER, SMOOTH_LAMBDA)
}

chan_list <- list(
  ILR1 = ilr1_s,
  ILR2 = ilr2_s,
  c12  = cov_s[[1]],
  c13  = cov_s[[2]],
  c23  = cov_s[[3]],
  v1   = cov_s[[4]],
  v2   = cov_s[[5]],
  v3   = cov_s[[6]]
)

col_vec <- c(
  "#1b9e77",  # AMP
  "#d95f02",  # ATP
  "#7570b3",  # FAD
  "#e7298a",  # FMN
  "#66a61e",  # GLC
  "#e6ab02",  # HEM
  "#a6761d",  # NAD
  "#666666"   # PO4
)
names(col_vec) <- target_classes

lty_vec <- c(1, 1, 1, 1, 2, 2, 2, 2)
names(lty_vec) <- target_classes
x_ticks <- seq(from = 5, to = 19, by = 2)

group_mean_rows <- list()

for (nm in chan_names) {
  mat <- chan_list[[nm]]

  means_by_cls <- lapply(classes, function(cl) {
    idx <- which(labels == cl)
    colMeans(mat[idx, , drop = FALSE])
  })
  means_mat <- do.call(rbind, means_by_cls)
  rownames(means_mat) <- classes

  group_mean_rows[[nm]] <- data.frame(
    channel = nm,
    label = rep(classes, each = length(thr_use)),
    threshold = rep(thr_use, times = length(classes)),
    mean_value = as.vector(t(means_mat)),
    stringsAsFactors = FALSE
  )

  ymin <- min(means_mat, na.rm = TRUE)
  ymax <- max(means_mat, na.rm = TRUE)
  fname <- file.path(fig_dir, paste0("mean_", nm, "_groupcurves.png"))

  png(
    filename = fname,
    width = 7,
    height = 7,
    units = "in",
    res = 600,
    bg = "white"
  )

  par(mar = c(4.5, 4.8, 0.5, 0.5), cex.lab = 1.6, cex.axis = 1.3)

  plot(
    thr_use, means_mat[1, ],
    type = "n",
    ylim = c(ymin, ymax),
    xlab = "Distance threshold (\u00C5)",
    ylab = lab_fun(nm),
    xaxt = "n"
  )

  axis(1, at = x_ticks, labels = x_ticks)

  for (j in seq_along(classes)) {
    cl <- classes[j]
    lines(
      thr_use,
      means_mat[j, ],
      lwd = 4.5,
      col = col_vec[cl],
      lty = lty_vec[cl]
    )
  }

  dev.off()
  cat("Saved:", fname, "\n")
}

group_means_df <- do.call(rbind, group_mean_rows)
means_csv <- file.path(out_dir, "EK_group_mean_functions_rawILR_CDPA_4p8to20A.csv")
write.csv(group_means_df, means_csv, row.names = FALSE)

metadata <- data.frame(
  dataset = "Extended Kahraman",
  t_start = T_START,
  t_end = T_END,
  t_step = T_STEP,
  n_thresholds = length(thr_use),
  n_sites = n,
  spline_degree = SPLINE_DEGREE,
  nbasis = nbasis,
  smoothing_lambda = SMOOTH_LAMBDA,
  ilr_scaling_applied = FALSE,
  output_directory = out_dir,
  stringsAsFactors = FALSE
)
metadata_csv <- file.path(out_dir, "EK_group_mean_functions_rawILR_CDPA_4p8to20A_metadata.csv")
write.csv(metadata, metadata_csv, row.names = FALSE)

cat("\nSaved group mean values:\n  ", means_csv, "\n", sep = "")
cat("Saved metadata:\n  ", metadata_csv, "\n", sep = "")
