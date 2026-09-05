# =====================================================================
# Extended Kahraman | AMP within-group ILR1(t) and c12(t) curves
#   - Smooths raw ILR1 and CDPA c12 over 4.8--20.0 Angstrom.
#   - Plots AMP individual curves, mean, and mean and SD for both channels.

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
    "c12" = expression(c[12]),
    nm
  )
}

ac <- function(col, alpha = 0.25) {
  grDevices::adjustcolor(col, alpha.f = alpha)
}

obj <- readRDS(cache_file)
thr_full <- seq(T_START, by = T_STEP, length.out = dim(obj$ilr_arr)[2])
idx_use <- which(thr_full >= T_START & thr_full <= T_END)
thr_use <- thr_full[idx_use]
nbasis <- choose_nbasis(thr_use)

labels <- factor(as.character(obj$labels))
amp_idx <- which(labels == "AMP")
if (length(amp_idx) == 0L) {
  stop("No AMP binding sites found in EK cache.")
}

cat("\n===== EK AMP ILR1 and c12 within-group curves =====\n")
cat("AMP sites:", length(amp_idx), "\n")
cat("Threshold interval:", min(thr_use), "to", max(thr_use), "\n")
cat("Threshold grid points:", length(thr_use), "\n")
cat("B-spline nbasis:", nbasis, "\n")
cat("Smoothing lambda:", SMOOTH_LAMBDA, "\n")
cat("Important: ILR1 is plotted without ILR/CDPA scaling.\n\n")

x_ticks <- seq(from = 5, to = 19, by = 2)

plot_within_group <- function(mat, channel_name, fig_file) {
  mfun <- colMeans(mat, na.rm = TRUE)
  sfun <- apply(mat, 2, sd, na.rm = TRUE)

  ymin <- min(mfun - 2 * sfun, mat, na.rm = TRUE)
  ymax <- max(mfun + 2 * sfun, mat, na.rm = TRUE)

  png(
    filename = fig_file,
    width = 7,
    height = 7,
    units = "in",
    res = 600,
    bg = "white"
  )

  par(
    mar = c(4.5, 4.8, 0.5, 0.5),
    cex.lab = 1.6,
    cex.axis = 1.3
  )

  plot(
    thr_use, mfun,
    type = "n",
    ylim = c(ymin, ymax),
    xlab = "Distance threshold (\u00C5)",
    ylab = lab_fun(channel_name),
    xaxt = "n"
  )
  axis(1, at = x_ticks, labels = x_ticks)

  apply(mat, 1, function(x) {
    lines(thr_use, x, col = ac("grey", 0.7))
  })

  lines(thr_use, mfun, lwd = 3, col = "black")
  lines(thr_use, mfun + sfun, lwd = 2, lty = 2)
  lines(thr_use, mfun - sfun, lwd = 2, lty = 2)

  dev.off()

  list(mean = mfun, sd = sfun)
}

ilr1_all <- obj$ilr_arr[, idx_use, 1]
c12_all <- obj$cov_arr[, idx_use, 1]

ilr1_s <- smooth_channel(ilr1_all, thr_use, nbasis, NORDER, SMOOTH_LAMBDA)
c12_s <- smooth_channel(c12_all, thr_use, nbasis, NORDER, SMOOTH_LAMBDA)

ilr1_mat <- ilr1_s[amp_idx, , drop = FALSE]
c12_mat <- c12_s[amp_idx, , drop = FALSE]

ilr1_fig_file <- file.path(fig_dir, "AMP_ILR1_curves.png")
c12_fig_file <- file.path(fig_dir, "AMP_c12_curves.png")
csv_file <- file.path(out_dir, "EK_AMP_ILR1_c12_within_group_curves.csv")
metadata_file <- file.path(out_dir, "EK_AMP_ILR1_c12_within_group_curves_metadata.csv")

ilr1_summary <- plot_within_group(ilr1_mat, "ILR1", ilr1_fig_file)
c12_summary <- plot_within_group(c12_mat, "c12", c12_fig_file)

curve_df <- rbind(
  data.frame(
    channel = "ILR1",
    binding_site_index = rep(amp_idx, each = length(thr_use)),
    ligand_group = "AMP",
    threshold = rep(thr_use, times = length(amp_idx)),
    smooth_value = as.vector(t(ilr1_mat)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    channel = "c12",
    binding_site_index = rep(amp_idx, each = length(thr_use)),
    ligand_group = "AMP",
    threshold = rep(thr_use, times = length(amp_idx)),
    smooth_value = as.vector(t(c12_mat)),
    stringsAsFactors = FALSE
  )
)

summary_df <- rbind(
  data.frame(
    channel = "ILR1",
    threshold = thr_use,
    mean_value = ilr1_summary$mean,
    sd_value = ilr1_summary$sd
  ),
  data.frame(
    channel = "c12",
    threshold = thr_use,
    mean_value = c12_summary$mean,
    sd_value = c12_summary$sd
  )
)

out_df <- merge(curve_df, summary_df, by = c("channel", "threshold"), all.x = TRUE)
write.csv(out_df, csv_file, row.names = FALSE)

metadata <- data.frame(
  dataset = "Extended Kahraman",
  ligand_group = "AMP",
  n_amp_sites = length(amp_idx),
  t_start = T_START,
  t_end = T_END,
  t_step = T_STEP,
  n_thresholds = length(thr_use),
  spline_degree = SPLINE_DEGREE,
  nbasis = nbasis,
  smoothing_lambda = SMOOTH_LAMBDA,
  ilr_scaling_applied = FALSE,
  ilr1_figure_file = ilr1_fig_file,
  c12_figure_file = c12_fig_file,
  stringsAsFactors = FALSE
)
write.csv(metadata, metadata_file, row.names = FALSE)

cat("Saved figures:\n  ", ilr1_fig_file, "\n  ", c12_fig_file, "\n", sep = "")
cat("Saved plotted data:\n  ", csv_file, "\n", sep = "")
cat("Saved metadata:\n  ", metadata_file, "\n", sep = "")
