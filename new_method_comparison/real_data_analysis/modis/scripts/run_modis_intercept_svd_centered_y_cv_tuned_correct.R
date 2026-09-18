#!/usr/bin/env Rscript

options(digits = 16, warn = 1)

suppressPackageStartupMessages({
  library(mgcv)
})

root_dir <- normalizePath(Sys.getenv("GDSSGL_ROOT", getwd()), winslash = "/", mustWork = TRUE)
comparison_root <- file.path(root_dir, "new_method_comparison")
raw_run_root <- file.path(
  comparison_root, "real_data_analysis",
  "modis_n10000_intercept_svd_current_methods_run"
)
out_root <- file.path(
  comparison_root, "real_data_analysis",
  "modis_n10000_intercept_svd_centered_y_cv_tuned_correct"
)
modis_path <- Sys.getenv("MODIS_RDATA", file.path(root_dir, "data", "data_cleaned_small_expanded.RData"))
dirs <- c("scripts", "logs", "data", "results", "results/cv_fits", "reports")
invisible(lapply(file.path(out_root, dirs), dir.create, recursive = TRUE, showWarnings = FALSE))

log_file <- file.path(out_root, "logs", "run_modis_intercept_svd_centered_y_cv_tuned_correct.log")
zz <- file(log_file, open = "at")
sink(zz, split = TRUE)
sink(zz, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(zz), silent = TRUE)
}, add = TRUE)

cat("\nMODIS correct intercept/SVD centered-y CV-tuned run started:", format(Sys.time()), "\n")

`%||%` <- function(x, y) if (is.null(x)) y else x

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Could not write ", path)
  invisible(path)
}

append_row_csv <- function(row, path, key_cols) {
  old <- if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else row[0, ]
  if (nrow(old)) {
    key_old <- do.call(paste, c(old[, key_cols, drop = FALSE], sep = "\r"))
    key_new <- do.call(paste, c(row[, key_cols, drop = FALSE], sep = "\r"))
    old <- old[!(key_old %in% key_new), , drop = FALSE]
  }
  write_csv(rbind(old, row), path)
}

cache <- function(path, expr) {
  if (file.exists(path)) return(readRDS(path))
  value <- force(expr)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(value, tmp)
  if (!file.rename(tmp, path)) stop("Could not write ", path)
  value
}

metric_row <- function(method, pred, y, runtime_sec, status = "completed",
                       notes = "", beta0_centered = NA_real_) {
  resid <- y - pred
  data.frame(
    method = method,
    MSPE = mean(resid^2),
    RMSE = sqrt(mean(resid^2)),
    MAE = mean(abs(resid)),
    bias = mean(resid),
    test_R2 = 1 - sum(resid^2) / sum((y - mean(y))^2),
    correlation = suppressWarnings(stats::cor(y, pred)),
    beta0_centered = beta0_centered,
    runtime_sec = runtime_sec,
    status = status,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

scale_train_test <- function(train_x, test_x) {
  mu <- colMeans(train_x, na.rm = TRUE)
  sig <- apply(train_x, 2, stats::sd, na.rm = TRUE)
  sig[!is.finite(sig) | sig == 0] <- 1
  list(
    train = sweep(sweep(train_x, 2, mu, "-"), 2, sig, "/"),
    test = sweep(sweep(test_x, 2, mu, "-"), 2, sig, "/"),
    mean = mu,
    sd = sig
  )
}

make_lc_dummy_design <- function(train_lc, test_lc, levels_full = NULL, reference = NULL) {
  train_levels <- if (is.null(levels_full)) sort(unique(train_lc)) else levels_full
  ref_level <- reference %||% train_levels[1]
  dummy_levels <- train_levels[train_levels != ref_level]
  train_dummy <- sapply(dummy_levels, function(level) as.numeric(train_lc == level))
  test_dummy <- sapply(dummy_levels, function(level) as.numeric(test_lc == level))
  if (length(dummy_levels) == 1) {
    train_dummy <- matrix(train_dummy, ncol = 1)
    test_dummy <- matrix(test_dummy, ncol = 1)
  }
  colnames(train_dummy) <- paste0("LC_Type4_", dummy_levels)
  colnames(test_dummy) <- paste0("LC_Type4_", dummy_levels)
  list(train = train_dummy, test = test_dummy, reference = ref_level,
       levels = train_levels, dummy_levels = dummy_levels)
}

make_grid_data <- function(train_coords, train_X, grid_n = 40L) {
  grid_coords <- expand.grid(
    x = seq(min(train_coords[, 1]), max(train_coords[, 1]), length.out = grid_n),
    y = seq(min(train_coords[, 2]), max(train_coords[, 2]), length.out = grid_n)
  )
  grid_X <- matrix(rep(apply(train_X, 2, median, na.rm = TRUE), each = nrow(grid_coords)),
                   nrow = nrow(grid_coords), ncol = ncol(train_X))
  colnames(grid_X) <- colnames(train_X)
  list(coords = as.matrix(grid_coords), X = grid_X)
}

cv_metric <- function(method, pred_centered, fold_obj, runtime, params = list(),
                      status = "completed", error = "") {
  pred_raw <- pred_centered + fold_obj$y_mean_fold
  data.frame(
    method = method,
    fold = params$fold %||% NA_integer_,
    K_B = params$K_B %||% NA_real_,
    H = params$H %||% NA_real_,
    effective_H = params$effective_H %||% NA_real_,
    lambda0 = params$lambda0 %||% NA_real_,
    lambda1 = params$lambda1 %||% NA_real_,
    kappa2 = params$kappa2 %||% NA_real_,
    gam_k = params$gam_k %||% NA_real_,
    validation_mse_centered = mean((fold_obj$val$y - pred_centered)^2),
    validation_mse_raw = mean((fold_obj$y_val_raw - pred_raw)^2),
    fit_runtime_sec = runtime,
    status = status,
    error_message = error,
    n_train = fold_obj$n_train,
    n_val = fold_obj$n_val,
    stringsAsFactors = FALSE
  )
}

select_config <- function(df, param_cols, tie_tol = 0.01) {
  ok <- df[df$status == "completed", , drop = FALSE]
  ag <- aggregate(ok$validation_mse_raw, ok[, param_cols, drop = FALSE],
                  function(z) c(mean = mean(z), sd = sd(z), n = length(z)))
  metrics <- do.call(data.frame, ag)
  names(metrics)[seq_along(param_cols)] <- param_cols
  names(metrics)[(length(param_cols) + 1):ncol(metrics)] <-
    c("mean_validation_mse", "sd_validation_mse", "successful_folds")
  metrics$se_validation_mse <- metrics$sd_validation_mse / sqrt(metrics$successful_folds)
  metrics <- metrics[order(metrics$mean_validation_mse), ]
  best_mse <- metrics$mean_validation_mse[1]
  cand <- metrics[metrics$mean_validation_mse <= best_mse * (1 + tie_tol), , drop = FALSE]
  if ("K_B" %in% names(cand)) cand <- cand[order(cand$K_B), , drop = FALSE]
  if ("gam_k" %in% names(cand)) cand <- cand[order(cand$gam_k), , drop = FALSE]
  if ("lambda0" %in% names(cand)) cand <- cand[order(-cand$lambda0), , drop = FALSE]
  if ("lambda1" %in% names(cand)) cand <- cand[order(-cand$lambda1), , drop = FALSE]
  if ("kappa2" %in% names(cand)) cand <- cand[order(cand$kappa2), , drop = FALSE]
  list(summary = metrics, selected = cand[1, , drop = FALSE],
       exact_min = metrics[1, , drop = FALSE], tied = cand)
}

continuous_vars <- c(
  "red_reflectance", "NIR_reflectance", "blue_reflectance",
  "MIR_reflectance", "GPP", "LE", "view_zenith_angle",
  "sun_zenith_angle", "relative_azimuth_angle"
)
display_vars <- c(continuous_vars, "LC_Type4")
required_vars <- c("EVI", "scaled_x", "scaled_y", "LC_Type4", continuous_vars)

split_df <- utils::read.csv(file.path(raw_run_root, "data", "modis_n10000_intercept_svd_train_test_split.csv"))
load(modis_path)
complete <- stats::complete.cases(data_cleaned_small[, required_vars]) &
  is.finite(data_cleaned_small$EVI) & data_cleaned_small$EVI > -1
real_data_clean <- data_cleaned_small[complete, ]
sample_data <- real_data_clean[split_df$original_clean_row, ]
train_indices <- which(split_df$split == "train")
test_indices <- which(split_df$split == "test")

train_y_raw <- log(sample_data$EVI[train_indices] + 1)
test_y_raw <- log(sample_data$EVI[test_indices] + 1)
y_train_mean <- mean(train_y_raw)
train_y <- train_y_raw - y_train_mean
test_y <- test_y_raw - y_train_mean

lc <- make_lc_dummy_design(sample_data$LC_Type4[train_indices], sample_data$LC_Type4[test_indices])
train_X_raw <- cbind(as.matrix(sample_data[train_indices, continuous_vars]), lc$train)
test_X_raw <- cbind(as.matrix(sample_data[test_indices, continuous_vars]), lc$test)
var_names <- colnames(train_X_raw)
scaled_X <- scale_train_test(train_X_raw, test_X_raw)
train_coords <- as.matrix(sample_data[train_indices, c("scaled_x", "scaled_y")])
test_coords <- as.matrix(sample_data[test_indices, c("scaled_x", "scaled_y")])
colnames(train_coords) <- c("s1", "s2")
colnames(test_coords) <- c("s1", "s2")
grid_data <- make_grid_data(train_coords, scaled_X$train, 40L)
train_data <- list(y = train_y, X = scaled_X$train, coords = train_coords)
test_data <- list(y = test_y, X = scaled_X$test, coords = test_coords)

source(file.path(comparison_root, "R", "methods", "accelerated", "newssgl_intercept.R"))
source(file.path(comparison_root, "R", "methods", "wrappers.R"))
Rcpp::sourceCpp(file.path(comparison_root, "R", "methods", "accelerated", "original_ssgl_intercept.cpp"))
Rcpp::sourceCpp(file.path(comparison_root, "R", "methods", "accelerated", "full_svc_intercept.cpp"))

cv_fold_seed <- 123L
cv_folds_n <- 5L
set.seed(cv_fold_seed)
cv_folds <- sample(rep(seq_len(cv_folds_n), length.out = length(train_indices)))
fold_df <- data.frame(
  train_position = seq_along(train_indices),
  row_in_sample = train_indices,
  original_clean_row = split_df$original_clean_row[train_indices],
  cv_fold = cv_folds
)
write_csv(fold_df, file.path(out_root, "data", "modis_correct_cv_folds_seed123.csv"))

prepare_cv_fold <- function(fold) {
  val_pos <- which(cv_folds == fold)
  tr_pos <- setdiff(seq_along(train_indices), val_pos)
  y_tr_raw <- train_y_raw[tr_pos]
  y_val_raw <- train_y_raw[val_pos]
  y_mean_fold <- mean(y_tr_raw)
  X_scaled <- scale_train_test(train_X_raw[tr_pos, , drop = FALSE],
                               train_X_raw[val_pos, , drop = FALSE])
  tr_coords <- train_coords[tr_pos, , drop = FALSE]
  val_coords <- train_coords[val_pos, , drop = FALSE]
  list(
    train = list(y = y_tr_raw - y_mean_fold, X = X_scaled$train, coords = tr_coords),
    val = list(y = y_val_raw - y_mean_fold, X = X_scaled$test, coords = val_coords),
    grid = make_grid_data(tr_coords, X_scaled$train, 25L),
    y_val_raw = y_val_raw,
    y_mean_fold = y_mean_fold,
    n_train = length(tr_pos),
    n_val = length(val_pos)
  )
}

validation <- data.frame(
  check = c("same_split_file_exists", "train_y_centered", "x_train_means_zero",
            "x_train_sds_one", "gd_k6_svd_fullrank", "svc_k7_svd_fullrank",
            "cv_folds_n"),
  pass = c(
    file.exists(file.path(raw_run_root, "data", "modis_n10000_intercept_svd_train_test_split.csv")),
    abs(mean(train_y)) < 1e-14,
    max(abs(colMeans(scaled_X$train))) < 1e-10,
    max(abs(apply(scaled_X$train, 2, sd) - 1)) < 1e-10,
    {
      meta <- fit_centered_basis_metadata_fullrank(train_coords, 6L, method = "svd")
      phi <- apply_centered_basis_metadata_fullrank(train_coords, meta)
      max(abs(colMeans(phi))) < 1e-8 && qr(phi)$rank == ncol(phi) && ncol(phi) == 35L
    },
    {
      meta <- fit_centered_basis_metadata_fullrank(train_coords, 7L, method = "svd")
      phi <- apply_centered_basis_metadata_fullrank(train_coords, meta)
      max(abs(colMeans(phi))) < 1e-8 && qr(phi)$rank == ncol(phi) && ncol(phi) == 48L
    },
    length(unique(cv_folds)) == cv_folds_n
  )
)
write_csv(validation, file.path(out_root, "data", "modis_correct_cv_validation.csv"))
if (!all(validation$pass)) stop("Validation failed.")

write_csv(data.frame(
  item = c("data_source", "source_split", "train_n", "test_n", "response",
           "response_centering", "x_standardization", "cv_folds", "cv_fold_seed",
           "cv_mcmc", "final_mcmc", "gd_basis", "svc_basis", "ws_basis"),
  value = c(modis_path, raw_run_root, length(train_indices), length(test_indices),
            "log(EVI + 1)", "subtract fold/full training response mean; add back for prediction",
            "fold-training/full-training mean and SD", cv_folds_n, cv_fold_seed,
            "n_iter=1500 burn_in=500", "n_iter=3000 burn_in=1000",
            "SVD full-rank centered basis", "SVD full-rank centered basis",
            "original uncentered whole-surface basis")
), file.path(out_root, "data", "modis_correct_cv_manifest.csv"))

cv_dir <- file.path(out_root, "results", "cv_fits")
cv_mcmc <- list(n_iter = 1500L, burn_in = 500L)
final_mcmc <- list(n_iter = 3000L, burn_in = 1000L)

run_gd_cv <- function() {
  path <- file.path(out_root, "data", "gdssgl_correct_cv_grid_results.csv")
  grid <- expand.grid(K_B = c(5L, 6L, 7L), lambda0 = c(15, 20, 30),
                      lambda1 = 2, fold = seq_len(cv_folds_n))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    existing <- if (file.exists(path)) read.csv(path) else NULL
    if (!is.null(existing) && any(existing$K_B == g$K_B & existing$lambda0 == g$lambda0 &
                                  existing$lambda1 == g$lambda1 & existing$fold == g$fold)) next
    cat("CV GD-SSGL", i, "/", nrow(grid), "K", g$K_B, "lambda0", g$lambda0, "fold", g$fold, "\n")
    fobj <- prepare_cv_fold(g$fold)
    fit_path <- file.path(cv_dir, sprintf("gdssgl_K%s_l0_%s_l1_%s_fold%s.rds",
                                          g$K_B, g$lambda0, g$lambda1, g$fold))
    row <- tryCatch({
      fit <- cache(fit_path, {
        fit_newssgl_intercept_fast(
          fobj$train, fobj$val, fobj$grid,
          basis_config = list(n_basis = as.integer(g$K_B), full_rank_method = "svd"),
          model_config = list(lambda0 = g$lambda0, lambda1 = g$lambda1,
                              b_sigma = var(fobj$train$y) / 2),
          mcmc_config = cv_mcmc,
          seed = 20260914L + 1000L * g$fold + 10L * g$K_B + g$lambda0
        )
      })
      cv_metric("GD-SSGL", fit$pred_test, fobj, fit$runtime,
                list(fold = g$fold, K_B = g$K_B, H = g$K_B^2,
                     effective_H = fit$diagnostics$effective_basis_dimension,
                     lambda0 = g$lambda0, lambda1 = g$lambda1))
    }, error = function(e) {
      cv_metric("GD-SSGL", rep(NA_real_, fobj$n_val), fobj, NA_real_,
                list(fold = g$fold, K_B = g$K_B, H = g$K_B^2,
                     effective_H = g$K_B^2 - 1, lambda0 = g$lambda0,
                     lambda1 = g$lambda1),
                status = "failed", error = conditionMessage(e))
    })
    append_row_csv(row, path, c("method", "fold", "K_B", "lambda0", "lambda1"))
  }
}

run_ws_cv <- function() {
  path <- file.path(out_root, "data", "wsssgl_correct_cv_grid_results.csv")
  grid <- expand.grid(K_B = c(5L, 6L, 7L), lambda0 = c(15, 20, 30),
                      lambda1 = c(2, 3), fold = seq_len(cv_folds_n))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    existing <- if (file.exists(path)) read.csv(path) else NULL
    if (!is.null(existing) && any(existing$K_B == g$K_B & existing$lambda0 == g$lambda0 &
                                  existing$lambda1 == g$lambda1 & existing$fold == g$fold)) next
    cat("CV WS-SSGL", i, "/", nrow(grid), "K", g$K_B, "lambda0", g$lambda0, "lambda1", g$lambda1, "fold", g$fold, "\n")
    fobj <- prepare_cv_fold(g$fold)
    fit_path <- file.path(cv_dir, sprintf("wsssgl_K%s_l0_%s_l1_%s_fold%s.rds",
                                          g$K_B, g$lambda0, g$lambda1, g$fold))
    row <- tryCatch({
      fit <- cache(fit_path, {
        fit_original_ssgl_intercept(
          fobj$train, fobj$val, fobj$grid,
          basis_config = list(n_basis = as.integer(g$K_B)),
          model_config = list(lambda0 = g$lambda0, lambda1 = g$lambda1,
                              b_sigma = var(fobj$train$y) / 2),
          mcmc_config = cv_mcmc,
          seed = 20261914L + 1000L * g$fold + 10L * g$K_B + g$lambda0 + g$lambda1
        )
      })
      cv_metric("WS-SSGL", fit$pred_test, fobj, fit$runtime,
                list(fold = g$fold, K_B = g$K_B, H = g$K_B^2,
                     effective_H = g$K_B^2, lambda0 = g$lambda0,
                     lambda1 = g$lambda1))
    }, error = function(e) {
      cv_metric("WS-SSGL", rep(NA_real_, fobj$n_val), fobj, NA_real_,
                list(fold = g$fold, K_B = g$K_B, H = g$K_B^2,
                     effective_H = g$K_B^2, lambda0 = g$lambda0,
                     lambda1 = g$lambda1),
                status = "failed", error = conditionMessage(e))
    })
    append_row_csv(row, path, c("method", "fold", "K_B", "lambda0", "lambda1"))
  }
}

run_svc_cv <- function() {
  path <- file.path(out_root, "data", "gaussian_svc_correct_cv_grid_results.csv")
  grid <- expand.grid(K_B = c(5L, 6L, 7L), kappa2 = c(4, 8, 12),
                      fold = seq_len(cv_folds_n))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    existing <- if (file.exists(path)) read.csv(path) else NULL
    if (!is.null(existing) && any(existing$K_B == g$K_B & existing$kappa2 == g$kappa2 & existing$fold == g$fold)) next
    cat("CV Gaussian SVC", i, "/", nrow(grid), "K", g$K_B, "kappa2", g$kappa2, "fold", g$fold, "\n")
    fobj <- prepare_cv_fold(g$fold)
    fit_path <- file.path(cv_dir, sprintf("gaussian_svc_K%s_kappa_%s_fold%s.rds",
                                          g$K_B, g$kappa2, g$fold))
    row <- tryCatch({
      fit <- cache(fit_path, {
        fit_full_svc_no_selection_intercept(
          fobj$train, fobj$val, fobj$grid,
          basis_config = list(n_basis = as.integer(g$K_B),
                              full_rank_centered = TRUE,
                              full_rank_method = "svd"),
          model_config = list(kappa2_alpha = g$kappa2,
                              b_sigma = var(fobj$train$y) / 2),
          mcmc_config = cv_mcmc,
          seed = 20262914L + 1000L * g$fold + 10L * g$K_B + g$kappa2
        )
      })
      cv_metric("Gaussian SVC", fit$pred_test, fobj, fit$runtime,
                list(fold = g$fold, K_B = g$K_B, H = g$K_B^2,
                     effective_H = fit$diagnostics$effective_basis_dimension,
                     kappa2 = g$kappa2))
    }, error = function(e) {
      cv_metric("Gaussian SVC", rep(NA_real_, fobj$n_val), fobj, NA_real_,
                list(fold = g$fold, K_B = g$K_B, H = g$K_B^2,
                     effective_H = g$K_B^2 - 1, kappa2 = g$kappa2),
                status = "failed", error = conditionMessage(e))
    })
    append_row_csv(row, path, c("method", "fold", "K_B", "kappa2"))
  }
}

run_gam_cv <- function() {
  path <- file.path(out_root, "data", "gam_correct_cv_grid_results.csv")
  grid <- expand.grid(gam_k = c(8L, 10L, 12L, 15L), fold = seq_len(cv_folds_n))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    existing <- if (file.exists(path)) read.csv(path) else NULL
    if (!is.null(existing) && any(existing$gam_k == g$gam_k & existing$fold == g$fold)) next
    cat("CV GAM", i, "/", nrow(grid), "k", g$gam_k, "fold", g$fold, "\n")
    fobj <- prepare_cv_fold(g$fold)
    row <- tryCatch({
      train_df <- as.data.frame(fobj$train$X); names(train_df) <- paste0("X", seq_along(var_names))
      train_df$y <- fobj$train$y; train_df$s1 <- fobj$train$coords[, 1]; train_df$s2 <- fobj$train$coords[, 2]
      val_df <- as.data.frame(fobj$val$X); names(val_df) <- paste0("X", seq_along(var_names))
      val_df$s1 <- fobj$val$coords[, 1]; val_df$s2 <- fobj$val$coords[, 2]
      linear_terms <- paste0("X", seq_along(var_names), collapse = " + ")
      vc_terms <- paste0("s(s1, s2, by = X", seq_along(continuous_vars),
                         ", bs = 'tp', k = ", g$gam_k, ")",
                         collapse = " + ")
      form <- as.formula(paste0("y ~ ", linear_terms, " + ", vc_terms))
      st <- proc.time()[3]
      fit <- tryCatch(mgcv::gam(form, data = train_df, method = "REML"),
                      error = function(e) mgcv::bam(form, data = train_df,
                                                    method = "fREML", discrete = TRUE))
      pred <- as.numeric(predict(fit, newdata = val_df))
      cv_metric("Practical GAM", pred, fobj, proc.time()[3] - st,
                list(fold = g$fold, gam_k = g$gam_k))
    }, error = function(e) {
      cv_metric("Practical GAM", rep(NA_real_, fobj$n_val), fobj, NA_real_,
                list(fold = g$fold, gam_k = g$gam_k),
                status = "failed", error = conditionMessage(e))
    })
    append_row_csv(row, path, c("method", "fold", "gam_k"))
  }
}

run_gd_cv()
run_ws_cv()
run_svc_cv()
run_gam_cv()

gd_cv <- read.csv(file.path(out_root, "data", "gdssgl_correct_cv_grid_results.csv"), stringsAsFactors = FALSE)
ws_cv <- read.csv(file.path(out_root, "data", "wsssgl_correct_cv_grid_results.csv"), stringsAsFactors = FALSE)
svc_cv <- read.csv(file.path(out_root, "data", "gaussian_svc_correct_cv_grid_results.csv"), stringsAsFactors = FALSE)
gam_cv <- read.csv(file.path(out_root, "data", "gam_correct_cv_grid_results.csv"), stringsAsFactors = FALSE)

sel_gd <- select_config(gd_cv, c("K_B", "H", "effective_H", "lambda0", "lambda1"))
sel_ws <- select_config(ws_cv, c("K_B", "H", "effective_H", "lambda0", "lambda1"))
sel_svc <- select_config(svc_cv, c("K_B", "H", "effective_H", "kappa2"))
sel_gam <- select_config(gam_cv, c("gam_k"))
write_csv(sel_gd$summary, file.path(out_root, "data", "gdssgl_correct_cv_configuration_summary.csv"))
write_csv(sel_ws$summary, file.path(out_root, "data", "wsssgl_correct_cv_configuration_summary.csv"))
write_csv(sel_svc$summary, file.path(out_root, "data", "gaussian_svc_correct_cv_configuration_summary.csv"))
write_csv(sel_gam$summary, file.path(out_root, "data", "gam_correct_cv_configuration_summary.csv"))
write_csv(sel_gd$selected, file.path(out_root, "data", "gdssgl_correct_cv_selected_hyperparameters.csv"))
write_csv(sel_ws$selected, file.path(out_root, "data", "wsssgl_correct_cv_selected_hyperparameters.csv"))
write_csv(sel_svc$selected, file.path(out_root, "data", "gaussian_svc_correct_cv_selected_hyperparameters.csv"))
write_csv(sel_gam$selected, file.path(out_root, "data", "gam_correct_cv_selected_hyperparameters.csv"))

selected_table <- data.frame(
  method = c("GD-SSGL", "WS-SSGL", "Gaussian SVC", "Practical GAM", "Bayesian Lasso"),
  selected_hyperparameters = c(
    sprintf("K_B=%s, raw H=%s, effective H=%s, lambda0=%s, lambda1=%s",
            sel_gd$selected$K_B[1], sel_gd$selected$H[1],
            sel_gd$selected$effective_H[1], sel_gd$selected$lambda0[1],
            sel_gd$selected$lambda1[1]),
    sprintf("K_B=%s, H=%s, lambda0=%s, lambda1=%s",
            sel_ws$selected$K_B[1], sel_ws$selected$H[1],
            sel_ws$selected$lambda0[1], sel_ws$selected$lambda1[1]),
    sprintf("K_B=%s, raw H=%s, effective H=%s, kappa2=%s",
            sel_svc$selected$K_B[1], sel_svc$selected$H[1],
            sel_svc$selected$effective_H[1], sel_svc$selected$kappa2[1]),
    sprintf("k=%s", sel_gam$selected$gam_k[1]),
    "default wrapper priors; no external CV grid"
  ),
  mean_cv_mse = c(
    sel_gd$selected$mean_validation_mse[1],
    sel_ws$selected$mean_validation_mse[1],
    sel_svc$selected$mean_validation_mse[1],
    sel_gam$selected$mean_validation_mse[1],
    NA_real_
  ),
  stringsAsFactors = FALSE
)
write_csv(selected_table, file.path(out_root, "data", "modis_correct_cv_selected_hyperparameters_all_methods.csv"))

cat("Final refits started:", format(Sys.time()), "\n")
fits <- list()
fits[["GD-SSGL"]] <- cache(file.path(out_root, "results", "gdssgl_correct_cv_tuned_final_fit.rds"), {
  g <- sel_gd$selected[1, ]
  fit_newssgl_intercept_fast(
    train_data, test_data, grid_data,
    basis_config = list(n_basis = as.integer(g$K_B), full_rank_method = "svd"),
    model_config = list(lambda0 = g$lambda0, lambda1 = g$lambda1,
                        b_sigma = var(train_y) / 2),
    mcmc_config = final_mcmc,
    seed = 2026091401L
  )
})
fits[["WS-SSGL"]] <- cache(file.path(out_root, "results", "wsssgl_correct_cv_tuned_final_fit.rds"), {
  g <- sel_ws$selected[1, ]
  fit_original_ssgl_intercept(
    train_data, test_data, grid_data,
    basis_config = list(n_basis = as.integer(g$K_B)),
    model_config = list(lambda0 = g$lambda0, lambda1 = g$lambda1,
                        b_sigma = var(train_y) / 2),
    mcmc_config = final_mcmc,
    seed = 2026091402L
  )
})
fits[["Gaussian SVC"]] <- cache(file.path(out_root, "results", "gaussian_svc_correct_cv_tuned_final_fit.rds"), {
  g <- sel_svc$selected[1, ]
  fit_full_svc_no_selection_intercept(
    train_data, test_data, grid_data,
    basis_config = list(n_basis = as.integer(g$K_B),
                        full_rank_centered = TRUE,
                        full_rank_method = "svd"),
    model_config = list(kappa2_alpha = g$kappa2, b_sigma = var(train_y) / 2),
    mcmc_config = final_mcmc,
    seed = 2026091403L
  )
})
fits[["Bayesian Lasso"]] <- cache(file.path(out_root, "results", "bayesian_lasso_correct_final_fit.rds"), {
  fit_global_only_blasso_intercept(
    train_data, test_data, grid_data,
    model_config = list(b_sigma = var(train_y) / 2),
    mcmc_config = final_mcmc,
    seed = 2026091404L
  )
})
fits[["Practical GAM"]] <- cache(file.path(out_root, "results", "gam_correct_cv_tuned_final_fit.rds"), {
  g <- sel_gam$selected[1, ]
  train_df <- as.data.frame(train_data$X); names(train_df) <- paste0("X", seq_along(var_names))
  train_df$y <- train_data$y; train_df$s1 <- train_data$coords[, 1]; train_df$s2 <- train_data$coords[, 2]
  test_df <- as.data.frame(test_data$X); names(test_df) <- paste0("X", seq_along(var_names))
  test_df$s1 <- test_data$coords[, 1]; test_df$s2 <- test_data$coords[, 2]
  linear_terms <- paste0("X", seq_along(var_names), collapse = " + ")
  vc_terms <- paste0("s(s1, s2, by = X", seq_along(continuous_vars),
                     ", bs = 'tp', k = ", g$gam_k, ")",
                     collapse = " + ")
  form <- as.formula(paste0("y ~ ", linear_terms, " + ", vc_terms))
  st <- proc.time()[3]
  fit <- tryCatch(mgcv::gam(form, data = train_df, method = "REML"),
                  error = function(e) mgcv::bam(form, data = train_df,
                                                method = "fREML", discrete = TRUE))
  list(fit = fit, pred_test = as.numeric(predict(fit, newdata = test_df)),
       runtime = proc.time()[3] - st, formula = paste(deparse(form), collapse = " "))
})

predictions <- list()
metrics <- list()
for (method in names(fits)) {
  fit <- fits[[method]]
  pred_centered <- if (method == "Practical GAM") fit$pred_test else fit$pred_test
  pred_raw <- pred_centered + y_train_mean
  predictions[[method]] <- pred_raw
  beta0_centered <- if (method == "Practical GAM") {
    unname(stats::coef(fit$fit)["(Intercept)"] %||% NA_real_)
  } else {
    fit$beta0_mean %||% fit$diagnostics$beta0_mean %||% NA_real_
  }
  metrics[[method]] <- metric_row(
    method, pred_raw, test_y_raw,
    fit$runtime %||% NA_real_,
    beta0_centered = beta0_centered,
    notes = selected_table$selected_hyperparameters[match(method, selected_table$method)]
  )
}
metrics_df <- do.call(rbind, metrics)
metrics_df <- metrics_df[order(metrics_df$MSPE), ]
write_csv(metrics_df, file.path(out_root, "data", "modis_correct_cv_tuned_prediction_metrics_logscale.csv"))

pred_long <- do.call(rbind, lapply(names(predictions), function(method) {
  data.frame(method = method, observed_log = test_y_raw, predicted_log = predictions[[method]],
             residual_log = test_y_raw - predictions[[method]],
             observed_evi = sample_data$EVI[test_indices],
             predicted_evi = exp(predictions[[method]]) - 1,
             stringsAsFactors = FALSE)
}))
write_csv(pred_long, file.path(out_root, "data", "modis_correct_cv_tuned_predictions_long.csv"))
evi_metrics <- do.call(rbind, lapply(split(pred_long, pred_long$method), function(d) {
  resid <- d$observed_evi - d$predicted_evi
  data.frame(
    method = d$method[1],
    MSPE_EVI = mean(resid^2),
    RMSE_EVI = sqrt(mean(resid^2)),
    MAE_EVI = mean(abs(resid)),
    bias_EVI = mean(resid),
    test_R2_EVI = 1 - sum(resid^2) / sum((d$observed_evi - mean(d$observed_evi))^2),
    correlation_EVI = suppressWarnings(stats::cor(d$observed_evi, d$predicted_evi)),
    stringsAsFactors = FALSE
  )
}))
evi_metrics <- evi_metrics[order(evi_metrics$MSPE_EVI), ]
write_csv(evi_metrics, file.path(out_root, "data", "modis_correct_cv_tuned_prediction_metrics_eviscale.csv"))

gd <- fits[["GD-SSGL"]]
gd_pip <- data.frame(
  method = "GD-SSGL",
  variable = var_names,
  sampled_pip = gd$pip,
  rb_pip = gd$diagnostics$rb_pip,
  selected = gd$diagnostics$rb_pip >= 0.5,
  selection_target = "spatial deviation",
  stringsAsFactors = FALSE
)
write_csv(gd_pip, file.path(out_root, "data", "modis_correct_cv_tuned_gdssgl_pip.csv"))

ws <- fits[["WS-SSGL"]]
ws_pip <- data.frame(
  method = "WS-SSGL",
  variable = var_names,
  sampled_pip = ws$pip,
  rb_pip = ws$diagnostics$rb_pip,
  selected = ws$diagnostics$rb_pip >= 0.5,
  selection_target = "whole coefficient surface",
  stringsAsFactors = FALSE
)
write_csv(ws_pip, file.path(out_root, "data", "modis_correct_cv_tuned_wsssgl_pip.csv"))

group_gd <- do.call(rbind, lapply(display_vars, function(v) {
  if (v == "LC_Type4") {
    rows <- gd_pip[grepl("^LC_Type4_", gd_pip$variable), ]
    data.frame(variable = v, rb_pip = max(rows$rb_pip),
               selected = max(rows$rb_pip) >= 0.5,
               component = rows$variable[which.max(rows$rb_pip)],
               selection_target = "spatial deviation")
  } else {
    rows <- gd_pip[gd_pip$variable == v, ]
    data.frame(variable = v, rb_pip = rows$rb_pip,
               selected = rows$selected,
               component = v,
               selection_target = "spatial deviation")
  }
}))
group_ws <- do.call(rbind, lapply(display_vars, function(v) {
  if (v == "LC_Type4") {
    rows <- ws_pip[grepl("^LC_Type4_", ws_pip$variable), ]
    data.frame(variable = v, rb_pip = max(rows$rb_pip),
               selected = max(rows$rb_pip) >= 0.5,
               component = rows$variable[which.max(rows$rb_pip)],
               selection_target = "whole coefficient surface")
  } else {
    rows <- ws_pip[ws_pip$variable == v, ]
    data.frame(variable = v, rb_pip = rows$rb_pip,
               selected = rows$selected,
               component = v,
               selection_target = "whole coefficient surface")
  }
}))
write_csv(rbind(group_gd, group_ws),
          file.path(out_root, "data", "modis_correct_cv_tuned_grouped_pip_summary.csv"))

report <- c(
  "# MODIS n=10000 Correct Intercept/SVD Centered-y CV-tuned Results",
  "",
  paste0("- Output root: `", out_root, "`."),
  "- This rerun uses the current explicit-intercept Bayesian implementations.",
  "- GD-SSGL and Gaussian SVC use SVD full-rank centered bases.",
  "- WS-SSGL uses the original uncentered whole-surface basis.",
  "- CV is 5-fold training-only with fold-specific response centering and predictor standardization.",
  "- CV MCMC: n_iter=1500, burn_in=500. Final MCMC: n_iter=3000, burn_in=1000.",
  "",
  "## Selected configurations",
  "",
  paste(capture.output(print(selected_table)), collapse = "\n"),
  "",
  "## Test performance on log(EVI+1) scale",
  "",
  paste(capture.output(print(metrics_df)), collapse = "\n"),
  "",
  "## Test performance on EVI scale",
  "",
  paste(capture.output(print(evi_metrics)), collapse = "\n"),
  "",
  "## GD-SSGL grouped spatial-deviation PIPs",
  "",
  paste(capture.output(print(group_gd)), collapse = "\n"),
  "",
  "## WS-SSGL grouped whole-surface PIPs",
  "",
  paste(capture.output(print(group_ws)), collapse = "\n")
)
writeLines(report, file.path(out_root, "reports", "modis_correct_cv_tuned_report.md"))

cat("MODIS correct intercept/SVD centered-y CV-tuned run finished:", format(Sys.time()), "\n")
cat("Output root:", out_root, "\n")
cat("Selected configurations:\n")
print(selected_table)
cat("EVI metrics:\n")
print(evi_metrics)
cat("GD grouped PIPs:\n")
print(group_gd)
cat("WS grouped PIPs:\n")
print(group_ws)
