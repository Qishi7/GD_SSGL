#!/usr/bin/env Rscript

options(digits = 16, warn = 1)

suppressPackageStartupMessages({
  library(GWmodel)
  library(sp)
})

root_dir <- normalizePath(Sys.getenv("GDSSGL_ROOT", getwd()),
                          winslash = "/", mustWork = TRUE)
comparison_root <- file.path(root_dir, "new_method_comparison")
out_root <- file.path(
  comparison_root, "real_data_analysis",
  "modis_n10000_intercept_svd_centered_y_check"
)
raw_run_root <- file.path(
  comparison_root, "real_data_analysis",
  "modis_n10000_intercept_svd_current_methods_run"
)
modis_path <- Sys.getenv("MODIS_RDATA",
                         file.path(root_dir, "data",
                                   "data_cleaned_small_expanded.RData"))
dirs <- c("logs", "data", "results", "reports")
invisible(lapply(file.path(out_root, dirs), dir.create, recursive = TRUE, showWarnings = FALSE))

mgwr_mode <- Sys.getenv("MODIS_MGWR_MODE", unset = "auto")
if (!mgwr_mode %in% c("auto", "practical_fixedbw", "practical_localwls")) {
  stop("Unsupported MODIS_MGWR_MODE: ", mgwr_mode)
}
suffix <- if (mgwr_mode == "auto") "" else paste0("_", mgwr_mode)
localwls_k <- as.integer(Sys.getenv("MODIS_MGWR_LOCALWLS_K", unset = "1000"))
mgwr_nlower <- as.integer(Sys.getenv("MODIS_MGWR_NLOWER", unset = "10"))
mgwr_bws0 <- Sys.getenv("MODIS_MGWR_BWS0", unset = "")
mgwr_bws0 <- if (nzchar(mgwr_bws0)) as.integer(mgwr_bws0) else NA_integer_
mgwr_hatmatrix <- tolower(Sys.getenv("MODIS_MGWR_HATMATRIX", unset = "false")) %in%
  c("true", "t", "1", "yes", "y")
mgwr_force_armadillo <- tolower(Sys.getenv("MODIS_MGWR_FORCE_ARMADILLO", unset = "true")) %in%
  c("true", "t", "1", "yes", "y")

log_file <- file.path(out_root, "logs", paste0("run_modis_centered_y_mgwr_only_evi_metrics", suffix, ".log"))
zz <- file(log_file, open = "at")
sink(zz, split = TRUE)
sink(zz, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(zz), silent = TRUE)
}, add = TRUE)

cat("\nMODIS centered-y MGWR-only run started:", format(Sys.time()), "\n")
cat("MGWR mode:", mgwr_mode, "\n")
cat("MGWR nlower:", mgwr_nlower, "\n")
cat("MGWR warm-start bws0:", ifelse(is.na(mgwr_bws0), "NULL", mgwr_bws0), "\n")
cat("MGWR hatmatrix:", mgwr_hatmatrix, "\n")
cat("MGWR force.armadillo:", mgwr_force_armadillo, "\n")

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Could not write ", path)
  invisible(path)
}

save_rds_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Could not write ", path)
  invisible(path)
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

make_lc_dummy_design <- function(train_lc, test_lc) {
  train_levels <- sort(unique(train_lc))
  ref_level <- train_levels[1]
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

make_spdf <- function(X, y, coords) {
  data_df <- data.frame(y = as.numeric(y))
  for (j in seq_len(ncol(X))) data_df[[paste0("X", j)]] <- X[, j]
  sp::coordinates(data_df) <- coords
  data_df
}

extract_mgwr_coef <- function(fit, p) {
  d <- fit$SDF@data
  intercept_col <- intersect(c("Intercept", "(Intercept)"), names(d))[1]
  x_cols <- paste0("X", seq_len(p))
  missing <- setdiff(c(intercept_col, x_cols), names(d))
  if (length(missing)) stop("Missing MGWR coefficient columns: ", paste(missing, collapse = ", "))
  out <- as.matrix(d[, c(intercept_col, x_cols), drop = FALSE])
  colnames(out) <- c("Intercept", x_cols)
  out
}

idw_coefficients <- function(train_coords, coef_mat, new_coords, power = 2) {
  train_coords <- as.matrix(train_coords)
  new_coords <- as.matrix(new_coords)
  coef_mat <- as.matrix(coef_mat)
  out <- matrix(NA_real_, nrow(new_coords), ncol(coef_mat))
  colnames(out) <- colnames(coef_mat)
  for (i in seq_len(nrow(new_coords))) {
    d <- sqrt(rowSums((train_coords - matrix(new_coords[i, ],
                                             nrow(train_coords),
                                             ncol(train_coords),
                                             byrow = TRUE))^2))
    exact <- which(d < 1e-12)
    if (length(exact)) {
      out[i, ] <- coef_mat[exact[1], ]
    } else {
      w <- 1 / d^power
      w <- w / sum(w)
      out[i, ] <- colSums(coef_mat * w)
    }
  }
  out
}

predict_local_wls <- function(train_X, train_y, train_coords, test_X, test_coords,
                              k = 1000L, ridge = 1e-8) {
  train_X <- as.matrix(train_X)
  test_X <- as.matrix(test_X)
  train_coords <- as.matrix(train_coords)
  test_coords <- as.matrix(test_coords)
  n_train <- nrow(train_X)
  n_test <- nrow(test_X)
  p <- ncol(train_X)
  k <- max(p + 2L, min(as.integer(k), n_train - 1L))
  X_design <- cbind(Intercept = 1, train_X)
  pred <- numeric(n_test)
  coef_test <- matrix(NA_real_, nrow = n_test, ncol = p + 1L)
  colnames(coef_test) <- colnames(X_design)
  for (i in seq_len(n_test)) {
    d <- sqrt(rowSums((train_coords - matrix(test_coords[i, ],
                                             nrow(train_coords),
                                             ncol(train_coords),
                                             byrow = TRUE))^2))
    kth <- sort(d, partial = k)[k]
    if (!is.finite(kth) || kth <= 0) kth <- max(d)
    idx <- which(d <= kth)
    if (length(idx) < p + 2L) idx <- order(d)[seq_len(min(n_train, p + 2L))]
    di <- d[idx]
    bw <- max(di)
    if (!is.finite(bw) || bw <= 0) bw <- max(d)
    wi <- (1 - (di / bw)^2)^2
    wi[di > bw] <- 0
    if (!any(wi > 0)) wi[] <- 1
    Xi <- X_design[idx, , drop = FALSE]
    yi <- train_y[idx]
    sw <- sqrt(wi)
    Xw <- Xi * sw
    yw <- yi * sw
    XtX <- crossprod(Xw)
    diag(XtX) <- diag(XtX) + ridge
    beta <- tryCatch(
      solve(XtX, crossprod(Xw, yw)),
      error = function(e) qr.solve(XtX, crossprod(Xw, yw))
    )
    beta <- as.numeric(beta)
    coef_test[i, ] <- beta
    pred[i] <- beta[1] + sum(test_X[i, ] * beta[-1])
  }
  list(pred = pred, coef_test = coef_test, k = k)
}

metric_row <- function(method, observed, predicted, runtime_sec, scale_name, beta0_hat = NA_real_) {
  resid <- observed - predicted
  data.frame(
    method = method,
    scale = scale_name,
    MSPE = mean(resid^2),
    RMSE = sqrt(mean(resid^2)),
    MAE = mean(abs(resid)),
    bias = mean(resid),
    test_R2 = 1 - sum(resid^2) / sum((observed - mean(observed))^2),
    correlation = suppressWarnings(stats::cor(observed, predicted)),
    beta0_hat = beta0_hat,
    runtime_sec = runtime_sec,
    stringsAsFactors = FALSE
  )
}

continuous_vars <- c(
  "red_reflectance", "NIR_reflectance", "blue_reflectance",
  "MIR_reflectance", "GPP", "LE", "view_zenith_angle",
  "sun_zenith_angle", "relative_azimuth_angle"
)
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

lc <- make_lc_dummy_design(sample_data$LC_Type4[train_indices], sample_data$LC_Type4[test_indices])
train_X_raw <- cbind(as.matrix(sample_data[train_indices, continuous_vars]), lc$train)
test_X_raw <- cbind(as.matrix(sample_data[test_indices, continuous_vars]), lc$test)
var_names <- colnames(train_X_raw)
scaled_X <- scale_train_test(train_X_raw, test_X_raw)
train_coords <- as.matrix(sample_data[train_indices, c("scaled_x", "scaled_y")])
test_coords <- as.matrix(sample_data[test_indices, c("scaled_x", "scaled_y")])
colnames(train_coords) <- c("s1", "s2")
colnames(test_coords) <- c("s1", "s2")

validation <- data.frame(
  check = c("same_split_file_exists", "train_y_centered", "x_train_means_zero",
            "x_train_sds_one", "mgwr_formula_has_default_intercept"),
  pass = c(
    file.exists(file.path(raw_run_root, "data", "modis_n10000_intercept_svd_train_test_split.csv")),
    abs(mean(train_y)) < 1e-14,
    max(abs(colMeans(scaled_X$train))) < 1e-10,
    max(abs(apply(scaled_X$train, 2, sd) - 1)) < 1e-10,
    attr(stats::terms(as.formula(paste("y ~", paste(paste0("X", seq_len(ncol(scaled_X$train))),
                                                    collapse = " + ")))), "intercept") == 1
  ),
  stringsAsFactors = FALSE
)
write_csv(validation, file.path(out_root, "data", "modis_centered_y_mgwr_validation.csv"))
if (!all(validation$pass)) stop("MGWR validation failed.")

fit_path <- file.path(out_root, "results", paste0("mgwr_intercept_centered_y", suffix, "_fit.rds"))
pred_path <- file.path(out_root, "data", paste0("modis_centered_y_mgwr", suffix, "_predictions.csv"))
metrics_path <- file.path(out_root, "data", paste0("modis_centered_y_mgwr", suffix, "_prediction_metrics_log_and_evi_scale.csv"))

formula_obj <- as.formula(
  paste("y ~", paste(paste0("X", seq_len(ncol(scaled_X$train))), collapse = " + "))
)
cat("MGWR formula:", paste(deparse(formula_obj), collapse = " "), "\n")
cat("Formula intercept attr:", attr(stats::terms(formula_obj), "intercept"), "\n")

if (file.exists(fit_path)) {
  obj <- readRDS(fit_path)
  cat("Loaded existing MGWR fit cache:", fit_path, "\n")
} else {
  st <- proc.time()[3]
  if (mgwr_mode == "practical_localwls") {
    lw <- predict_local_wls(
      scaled_X$train, train_y, train_coords,
      scaled_X$test, test_coords,
      k = localwls_k
    )
    fit <- list(
      mode = mgwr_mode,
      formula = paste(deparse(formula_obj), collapse = " "),
      adaptive = TRUE,
      kernel = "bisquare",
      localwls_k = lw$k,
      note = "Practical MGWR/GWR-family local weighted least squares approximation; no multiscale bandwidth search."
    )
    runtime <- proc.time()[3] - st
    coef_test <- lw$coef_test
    pred_log_centered <- lw$pred
  } else {
    train_spdf <- make_spdf(scaled_X$train, train_y, train_coords)
    mgwr_var_n <- ncol(scaled_X$train) + 1L
    fit_args <- list(
      formula = formula_obj, data = train_spdf, adaptive = TRUE,
      kernel = "bisquare", max.iterations = 1000,
      criterion = "CVR", verbose = TRUE, nlower = mgwr_nlower,
      hatmatrix = mgwr_hatmatrix,
      force.armadillo = mgwr_force_armadillo
    )
    if (!is.na(mgwr_bws0)) {
      fit_args$bws0 <- rep(mgwr_bws0, mgwr_var_n)
      fit_args$bw.seled <- rep(FALSE, mgwr_var_n)
    }
    if (mgwr_mode == "practical_fixedbw") {
      fit_args$bws0 <- rep(length(train_y) - 1L, mgwr_var_n)
      fit_args$bw.seled <- rep(TRUE, mgwr_var_n)
      fit_args$hatmatrix <- FALSE
      fit_args$verbose <- FALSE
      fit_args$force.armadillo <- TRUE
    }
    fit <- do.call(GWmodel::gwr.multiscale, fit_args)
    runtime <- proc.time()[3] - st
    coef_train <- extract_mgwr_coef(fit, ncol(scaled_X$train))
    coef_test <- idw_coefficients(sp::coordinates(fit$SDF), coef_train, test_coords)
    pred_log_centered <- as.numeric(coef_test[, "Intercept"] +
                                      rowSums(scaled_X$test * coef_test[, -1, drop = FALSE]))
  }
  obj <- list(
    fit = fit,
    formula = paste(deparse(formula_obj), collapse = " "),
    intercept_attr = attr(stats::terms(formula_obj), "intercept"),
    y_train_mean = y_train_mean,
    pred_test_log_centered = pred_log_centered,
    pred_test_log = pred_log_centered + y_train_mean,
    pred_test_evi = exp(pred_log_centered + y_train_mean) - 1,
    beta0_centered_mean = mean(coef_test[, "Intercept"]),
    beta0_raw_log_scale_mean = mean(coef_test[, "Intercept"]) + y_train_mean,
    runtime_sec = runtime,
    var_names = var_names,
    mgwr_mode = mgwr_mode
  )
  save_rds_atomic(obj, fit_path)
}

observed_evi <- sample_data$EVI[test_indices]
method_label <- switch(
  mgwr_mode,
  auto = "MGWR",
  practical_fixedbw = "MGWR practical fixed-bw",
  practical_localwls = paste0("MGWR practical local-WLS k=", obj$fit$localwls_k)
)
pred_df <- data.frame(
  method = method_label,
  observed_log = test_y_raw,
  predicted_log = obj$pred_test_log,
  residual_log = test_y_raw - obj$pred_test_log,
  observed_evi = observed_evi,
  predicted_evi = obj$pred_test_evi,
  residual_evi = observed_evi - obj$pred_test_evi,
  s1 = test_coords[, 1],
  s2 = test_coords[, 2]
)
write_csv(pred_df, pred_path)

metrics <- rbind(
  metric_row(method_label, test_y_raw, obj$pred_test_log, obj$runtime_sec,
             "log(EVI+1)", beta0_hat = obj$beta0_raw_log_scale_mean),
  metric_row(method_label, observed_evi, obj$pred_test_evi, obj$runtime_sec,
             "EVI", beta0_hat = obj$beta0_raw_log_scale_mean)
)
write_csv(metrics, metrics_path)

all_evi_path <- file.path(out_root, "data", "modis_n10000_centered_y_check_prediction_metrics_evi_scale.csv")
if (file.exists(all_evi_path)) {
  all_evi <- utils::read.csv(all_evi_path)
  mgwr_evi <- metrics[metrics$scale == "EVI", ]
  mgwr_evi_out <- data.frame(
    method = method_label,
    MSPE_EVI = mgwr_evi$MSPE,
    RMSE_EVI = mgwr_evi$RMSE,
    MAE_EVI = mgwr_evi$MAE,
    bias_EVI = mgwr_evi$bias,
    test_R2_EVI = mgwr_evi$test_R2,
    correlation_EVI = mgwr_evi$correlation,
    min_pred_EVI = min(obj$pred_test_evi),
    max_pred_EVI = max(obj$pred_test_evi)
  )
  combined <- rbind(all_evi, mgwr_evi_out[, names(all_evi)])
  combined <- combined[order(combined$MSPE_EVI), ]
  write_csv(combined, file.path(out_root, "data", paste0("modis_n10000_centered_y_check_prediction_metrics_evi_scale_with_mgwr", suffix, ".csv")))
}

writeLines(c(
  "# MODIS centered-y MGWR EVI-scale metrics",
  "",
  paste0("- Output root: `", out_root, "`."),
  paste0("- MGWR mode: `", mgwr_mode, "`."),
  paste0("- MGWR adaptive bandwidth lower bound `nlower`: ", mgwr_nlower, "."),
  paste0("- MGWR warm-start `bws0`: ", ifelse(is.na(mgwr_bws0), "NULL", mgwr_bws0), "."),
  paste0("- MGWR `hatmatrix`: ", mgwr_hatmatrix, "."),
  paste0("- MGWR `force.armadillo`: ", mgwr_force_armadillo, "."),
  "- MGWR was fit on centered `log(EVI+1)` using an explicit-intercept formula.",
  paste0("- Formula: `", obj$formula, "`."),
  paste0("- Formula intercept attribute: ", obj$intercept_attr, "."),
  paste0("- Training response mean added back before EVI transform: ", sprintf("%.12f", obj$y_train_mean), "."),
  "- EVI-scale predictions use `exp(predicted_log) - 1`.",
  "",
  "## Metrics",
  "",
  paste(capture.output(print(metrics, row.names = FALSE)), collapse = "\n")
), file.path(out_root, "reports", paste0("modis_centered_y_mgwr", suffix, "_evi_metrics_report.md")))

cat("\nMGWR metrics:\n")
print(metrics, row.names = FALSE)
cat("\nSaved fit:", fit_path, "\n")
cat("Saved metrics:", metrics_path, "\n")
cat("Report:", file.path(out_root, "reports", paste0("modis_centered_y_mgwr", suffix, "_evi_metrics_report.md")), "\n")
cat("MODIS centered-y MGWR-only run finished:", format(Sys.time()), "\n")
