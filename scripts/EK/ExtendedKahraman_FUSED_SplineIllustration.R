# =====================================================================
# FUSED | Raw threshold trajectory and cubic B-spline illustration
#   - Dataset: Extended Kahraman
#   - Creates an illustrative figure.
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
SMOOTH_LAMBDA <- 1e-3
SPLINE_DEGREE <- 3
NORDER <- SPLINE_DEGREE + 1
BASIS_DIVS <- c(2, 3, 4, 5)

# Optional controls:
#   SITE_INDEX=120 Rscript scripts/EK/ExtendedKahraman_FUSED_SplineIllustration.R
#   SITE_LABEL=HEM  Rscript scripts/EK/ExtendedKahraman_FUSED_SplineIllustration.R

env_site_index <- suppressWarnings(as.integer(Sys.getenv("SITE_INDEX", "")))
env_site_label <- toupper(trimws(Sys.getenv("SITE_LABEL", "")))
if (is.na(env_site_index)) env_site_index <- NA_integer_

choose_nbasis <- function(tvec) {
  n_t <- length(tvec)
  cand_nb <- sort(unique(pmax(NORDER, round(n_t / BASIS_DIVS))))
  cand_nb[1]
}

fit_spline_values <- function(y, tvec, nbasis, norder, lambda) {
  basis <- create.bspline.basis(range(tvec), nbasis, norder)
  fd_param <- fdPar(basis, int2Lfd(2), lambda)
  fit <- smooth.basis(tvec, y, fd_param)
  as.numeric(eval.fd(tvec, fit$fd))
}

select_representative_site <- function(ilr1, v1, candidate_idx = seq_len(nrow(ilr1))) {
  variability <- apply(cbind(ilr1, v1), 1, function(z) {
    z1 <- z[seq_along(ilr1[1, ])]
    z2 <- z[-seq_along(ilr1[1, ])]
    sd(z1, na.rm = TRUE) + sd(z2, na.rm = TRUE)
  })
  finite_idx <- intersect(candidate_idx, which(is.finite(variability)))
  finite_idx[which.min(abs(variability[finite_idx] - median(variability[finite_idx])))]
}

obj <- readRDS(cache_file)
thr <- seq(T_START, by = T_STEP, length.out = dim(obj$ilr_arr)[2])
idx_use <- which(thr >= T_START & thr <= T_END)
tvec <- thr[idx_use]
nbasis <- choose_nbasis(tvec)

ilr1_mat <- obj$ilr_arr[, idx_use, 1]
v1_mat <- obj$cov_arr[, idx_use, 4]

if (!is.na(env_site_index)) {
  if (env_site_index < 1L || env_site_index > length(obj$labels)) {
    stop("SITE_INDEX must be between 1 and ", length(obj$labels), ".")
  }
  site_idx <- env_site_index
  selection_method <- "user_supplied_SITE_INDEX"
} else if (nzchar(env_site_label)) {
  label_idx <- which(toupper(as.character(obj$labels)) == env_site_label)
  if (length(label_idx) == 0L) {
    stop("No EK binding sites found for SITE_LABEL=", env_site_label, ".")
  }
  site_idx <- select_representative_site(ilr1_mat, v1_mat, candidate_idx = label_idx)
  selection_method <- paste0("representative_within_SITE_LABEL_", env_site_label)
} else {
  site_idx <- select_representative_site(ilr1_mat, v1_mat)
  selection_method <- "representative_all_sites_median_variability"
}

site_label <- as.character(obj$labels[site_idx])
site_file <- if (!is.null(obj$files_df$path)) obj$files_df$path[site_idx] else NA_character_

plot_df <- data.frame(
  threshold = tvec,
  ILR1_raw = as.numeric(ilr1_mat[site_idx, ]),
  ILR1_spline = fit_spline_values(
    as.numeric(ilr1_mat[site_idx, ]), tvec, nbasis, NORDER, SMOOTH_LAMBDA
  ),
  v1_raw = as.numeric(v1_mat[site_idx, ]),
  v1_spline = fit_spline_values(
    as.numeric(v1_mat[site_idx, ]), tvec, nbasis, NORDER, SMOOTH_LAMBDA
  )
)

metadata <- data.frame(
  dataset = "Extended Kahraman",
  selection_method = selection_method,
  selected_site_index = site_idx,
  selected_ligand_group = site_label,
  selected_pdb_file = site_file,
  t_start = T_START,
  t_end = T_END,
  t_step = T_STEP,
  n_thresholds = length(tvec),
  spline_degree = SPLINE_DEGREE,
  nbasis = nbasis,
  smoothing_lambda = SMOOTH_LAMBDA,
  stringsAsFactors = FALSE
)

csv_file <- file.path(out_dir, "EK_FUSED_spline_illustration_data.csv")
metadata_file <- file.path(out_dir, "EK_FUSED_spline_illustration_metadata.csv")
fig_file <- file.path(fig_dir, "EK_FUSED_spline_illustration.png")

write.csv(plot_df, csv_file, row.names = FALSE)
write.csv(metadata, metadata_file, row.names = FALSE)

png(fig_file, width = 7.2, height = 4.0, units = "in", res = 600, bg = "white")
op <- par(
  mfrow = c(1, 2),
  mar = c(4.6, 4.8, 1.0, 0.6),
  oma = c(0, 0, 0, 0),
  cex.lab = 0.70,
  cex.axis = 0.60,
  family = "Helvetica"
)

plot_component <- function(x, raw, smooth, ylab) {
  yr <- range(c(raw, smooth), finite = TRUE)
  ypad <- 0.08 * diff(yr)
  if (!is.finite(ypad) || ypad == 0) ypad <- 1

  plot(
    x, raw,
    type = "s",
    lwd = 1.4,
    col = "grey45",
    ylim = c(yr[1] - ypad, yr[2] + ypad),
    xlab = expression(paste("Distance threshold (", ring(A), ")")),
    ylab = ylab,
    las = 1,
    bty = "l"
  )
  points(x, raw, pch = 16, cex = 0.33, col = adjustcolor("grey35", alpha.f = 0.70))
  lines(x, smooth, lwd = 3.0, col = "#1b9e77")
  grid(col = adjustcolor("grey80", alpha.f = 0.65), lty = "dotted")
  box(bty = "l", lwd = 1.1)
}

plot_component(plot_df$threshold, plot_df$ILR1_raw, plot_df$ILR1_spline, expression(ILR[1](t)))
legend(
  "topright",
  legend = c("Raw values", "Cubic B-spline"),
  col = c("grey45", "#1b9e77"),
  lwd = c(1.4, 3.0),
  pch = c(16, NA),
  pt.cex = c(0.7, NA),
  bty = "n",
  cex = 0.9
)

plot_component(plot_df$threshold, plot_df$v1_raw, plot_df$v1_spline, expression(v[1](t)))

par(op)
dev.off()

cat("\n===== FUSED spline illustration =====\n")
print(metadata)
cat("\nSaved figure:\n  ", fig_file, "\n", sep = "")
cat("Saved plotted data:\n  ", csv_file, "\n", sep = "")
cat("Saved metadata:\n  ", metadata_file, "\n", sep = "")
