# =====================================================================
# TOUGH-M1 | Export raw FUSED curve matrix through 15 A
#   - Uses the raw TOUGH-M1 ILR/CDPA cache over [4.8, 15.0] A.
# =====================================================================

rm(list = ls()); invisible(gc())

cache_file <- file.path("data", "cache", "TOUGHM1_raw_ILR_CDPA_thr4p8_15p0_by0p1.rds")
out_file <- file.path("data", "cache", "TOUGHM1_rawFUSED_features_thr4p8_15p0.csv")

if (!file.exists(cache_file)) {
  stop("Missing raw TOUGH-M1 cache: ", cache_file,
       "\nRun scripts/TOUGH-M1/TOUGHM1_Build_RawFUSED_Descriptors_thr4p8_15p0.R first.")
}

raw <- readRDS(cache_file)
tvec <- raw$thr_full
ilr_arr <- raw$ilr_arr
cov_arr <- raw$cov_arr
ids <- raw$metadata$code5

channel_names <- c("ILR1", "ILR2", "c12", "c13", "c23", "v1", "v2", "v3")

X <- cbind(
  ilr_arr[, , 1],
  ilr_arr[, , 2],
  cov_arr[, , 1],
  cov_arr[, , 2],
  cov_arr[, , 3],
  cov_arr[, , 4],
  cov_arr[, , 5],
  cov_arr[, , 6]
)

make_names <- function(channel) {
  paste0(channel, "_t", gsub("\\.", "p", sprintf("%.1f", tvec)))
}

colnames(X) <- unlist(lapply(channel_names, make_names), use.names = FALSE)

out <- data.frame(code5 = ids, X, check.names = FALSE)
write.csv(out, out_file, row.names = FALSE)

cat("Exported raw FUSED feature matrix:\n")
cat("  rows:", nrow(out), "\n")
cat("  feature columns:", ncol(out) - 1, "\n")
cat("  file:", out_file, "\n")
