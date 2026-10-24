# Lightweight checks; no descriptor construction or model training.
script_dir <- file.path("scripts", "TOUGH-C1")
load_functions <- function(filename) {
  env <- new.env(parent = globalenv())
  for (expr in parse(file.path(script_dir, filename))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
      eval(expr, env)
    }
  }
  env
}
rf <- load_functions("TOUGHC1_RF_ThreeClass_RepeatedNestedCV.R")
mlr <- load_functions("TOUGHC1_MLR_NMM_ThreeClass_RepeatedNestedCV.R")
for (env in list(rf, mlr)) {
  stopifnot(env$auc_rank(c(TRUE, FALSE), c(0.5, 0.5)) == 0.5)
  stopifnot(env$auc_rank(c(TRUE, FALSE), c(1, 0)) == 1)
  y <- c("NUC", "CONTROL", "HEME", "NUC")
  predicted <- c("NUC", "NUC", "CONTROL", "CONTROL")
  scores <- cbind(CONTROL = c(.1, .2, .6, .65),
                  HEME = c(.2, .2, .3, .1), NUC = c(.7, .6, .1, .25))
  input_scores <- if (identical(env, rf)) scores else scores[, "NUC"]
  metrics <- env$one_vs_rest_metrics(y, predicted, input_scores, "NUC")
  stopifnot(metrics$accuracy == .5, metrics$ppv == .5,
            metrics$tpr == .5, metrics$tnr == .5, metrics$auc == .75)
}
full <- read.csv(file.path(script_dir, "folds",
  "TOUGHC1_threeclass_5x5_outer_fold_assignments.csv"))
subset <- read.csv(file.path(script_dir, "folds",
  "TOUGHC1_threeclass_5x5_outer_fold_assignments_common4000.csv"))
for (folds in list(full, subset)) {
  stopifnot(!anyDuplicated(folds[c("repetition", "descriptor_row_index")]),
            setequal(unique(folds$repetition), 1:5),
            setequal(unique(folds$outer_test_fold), 1:5))
  for (r in 1:5) {
    rows <- folds[folds$repetition == r, ]
    stopifnot(nrow(rows) == nrow(folds) / 5)
    for (f in 1:5) {
      train <- rows$descriptor_row_index[rows$outer_test_fold != f]
      test <- rows$descriptor_row_index[rows$outer_test_fold == f]
      stopifnot(length(intersect(train, test)) == 0,
                setequal(c(train, test), rows$descriptor_row_index))
    }
  }
}
stopifnot(nrow(full) == 4002 * 5, nrow(subset) == 4000 * 5,
          setequal(setdiff(full$p_base, subset$p_base), c("1ukmB", "2o8cD")))
expected <- full[!full$p_base %in% c("1ukmB", "2o8cD"), ]
rownames(expected) <- NULL
stopifnot(identical(expected, subset))
cat("PASS: three-class metrics, AUC ties, outer folds, and retained subset.\n")
