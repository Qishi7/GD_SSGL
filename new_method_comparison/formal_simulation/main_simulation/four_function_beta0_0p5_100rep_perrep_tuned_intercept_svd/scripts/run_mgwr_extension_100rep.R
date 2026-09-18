#!/usr/bin/env Rscript

options(digits = 16, warn = 1)

suppressPackageStartupMessages({
  library(GWmodel)
  library(sp)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = "") {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
workers <- max(1L, as.integer(arg_value("workers", "4")))

project_root <- normalizePath(Sys.getenv("GDSSGL_ROOT", getwd()), winslash = "/", mustWork = TRUE)
root <- file.path(project_root, "new_method_comparison", "formal_simulation", "main_simulation", "four_function_beta0_0p5_100rep_perrep_tuned_intercept_svd")
out_root <- file.path(root, "mgwr_extension_100rep")
invisible(lapply(file.path(out_root, c("fits", "results", "reports", "logs")),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

predictor_names <- paste0("X", 1:10)
theta_true <- c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0)
active_u <- 3:6

atomic_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Atomic rename failed: ", path)
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Atomic rename failed: ", path)
}

make_spdf <- function(X, y, coords) {
  data_df <- data.frame(y = as.numeric(y))
  for (j in seq_len(ncol(X))) data_df[[paste0("X", j)]] <- X[, j]
  sp::coordinates(data_df) <- coords
  data_df
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

extract_coef_matrix <- function(fit, p) {
  d <- fit$SDF@data
  intercept_col <- intersect(c("Intercept", "(Intercept)"), names(d))[1]
  x_cols <- paste0("X", seq_len(p))
  missing <- setdiff(c(intercept_col, x_cols), names(d))
  if (length(missing)) stop("Missing coefficient columns: ",
                            paste(missing, collapse = ", "))
  as.matrix(d[, c(intercept_col, x_cols), drop = FALSE])
}

predict_from_local_coefficients <- function(fit, X_new, coords_new, p) {
  coef_train <- extract_coef_matrix(fit, p)
  train_coords <- sp::coordinates(fit$SDF)
  coef_new <- idw_coefficients(train_coords, coef_train, coords_new)
  as.numeric(coef_new[, 1] + rowSums(X_new * coef_new[, -1, drop = FALSE]))
}

grid_coefficients <- function(fit, grid_coords, p) {
  coef_train <- extract_coef_matrix(fit, p)
  train_coords <- sp::coordinates(fit$SDF)
  coef_grid <- idw_coefficients(train_coords, coef_train, grid_coords)
  beta <- coef_grid[, -1, drop = FALSE]
  colnames(beta) <- paste0("X", seq_len(p))
  list(intercept = coef_grid[, 1], beta = beta)
}

run_one_rep <- function(rep_id) {
  rep_tag <- sprintf("rep_%03d", rep_id)
  rep_dir <- file.path(out_root, "fits", rep_tag)
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  complete_path <- file.path(rep_dir, "complete.rds")
  if (file.exists(complete_path)) return(readRDS(complete_path))

  log_file <- file.path(out_root, "logs", paste0(rep_tag, ".log"))
  log_con <- file(log_file, open = "at")
  sink(log_con, type = "output")
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(log_con)
  }, add = TRUE)

  result <- tryCatch({
    dat <- readRDS(file.path(root, "fits", rep_tag,
                             "shared_dataset_beta0_0p5.rds"))
    p <- ncol(dat$train$X)
    train_spdf <- make_spdf(dat$train$X, dat$train$y, dat$train$coords)
    formula_obj <- as.formula(
      paste("y ~", paste(paste0("X", seq_len(p)), collapse = " + "))
    )
    has_formula_intercept <- attr(terms(formula_obj), "intercept")

    fit_path <- file.path(rep_dir, "mgwr_fit.rds")
    started <- proc.time()[3]
    if (file.exists(fit_path)) {
      fit <- readRDS(fit_path)
      fit_runtime <- NA_real_
    } else {
      fit <- GWmodel::gwr.multiscale(
        formula_obj, data = train_spdf, adaptive = TRUE,
        kernel = "bisquare", max.iterations = 1000,
        criterion = "CVR", verbose = FALSE
      )
      fit_runtime <- proc.time()[3] - started
      atomic_save_rds(fit, fit_path)
    }

    post_started <- proc.time()[3]
    pred <- predict_from_local_coefficients(fit, dat$test$X, dat$test$coords, p)
    grid_coef <- grid_coefficients(fit, dat$grid$coords, p)
    beta_grid <- grid_coef$beta
    post_runtime <- proc.time()[3] - post_started

    beta_mise_j <- colMeans((beta_grid - dat$grid$true_beta)^2)
    theta_hat <- colMeans(beta_grid)
    u_hat <- sweep(beta_grid, 2, theta_hat, "-")
    theta_mse <- mean((theta_hat - dat$truth$theta)^2)
    u_mise_j <- colMeans((u_hat - dat$grid$true_u)^2)

    d <- fit$SDF@data
    intercept_col <- intersect(c("Intercept", "(Intercept)"), names(d))[1]
    bw <- if (!is.null(fit$GW.arguments$bws)) {
      fit$GW.arguments$bws
    } else if (!is.null(fit$GW.arguments$bw)) {
      fit$GW.arguments$bw
    } else {
      NA_real_
    }

    metrics <- data.frame(
      replicate = rep_id,
      method = "MGWR",
      status = "success",
      true_beta0 = dat$truth$beta0,
      beta0_hat = mean(grid_coef$intercept),
      beta0_error = mean(grid_coef$intercept) - dat$truth$beta0,
      mspe = mean((dat$test$y - pred)^2),
      rmse = sqrt(mean((dat$test$y - pred)^2)),
      mae = mean(abs(dat$test$y - pred)),
      bias = mean(dat$test$y - pred),
      beta_mise = mean(beta_mise_j),
      mise_global_only = mean(beta_mise_j[1:2]),
      mise_spatial_only = mean(beta_mise_j[3:4]),
      mise_global_plus_spatial = mean(beta_mise_j[5:6]),
      mise_null = mean(beta_mise_j[7:10]),
      theta_mse = theta_mse,
      u_mise_x3_x6 = mean(u_mise_j[active_u]),
      runtime_sec = fit_runtime,
      postprocess_runtime_sec = post_runtime,
      total_runtime_sec = fit_runtime + post_runtime,
      formula = deparse(formula_obj),
      terms_intercept = has_formula_intercept,
      sdf_has_intercept = !is.na(intercept_col),
      intercept_column = intercept_col,
      bandwidth_summary = paste(signif(as.numeric(bw), 6), collapse = ";"),
      stringsAsFactors = FALSE
    )
    by_predictor <- data.frame(
      replicate = rep_id,
      method = "MGWR",
      predictor = paste0("X", seq_len(p)),
      theta_hat = theta_hat,
      theta_true = dat$truth$theta,
      theta_error = theta_hat - dat$truth$theta,
      beta_mise = beta_mise_j,
      u_mise = u_mise_j,
      stringsAsFactors = FALSE
    )
    atomic_write_csv(metrics, file.path(rep_dir, "mgwr_metrics.csv"))
    atomic_write_csv(by_predictor, file.path(rep_dir, "mgwr_by_predictor.csv"))
    out <- list(metrics = metrics, by_predictor = by_predictor,
                status = data.frame(replicate = rep_id, status = "success",
                                    error_message = "", stringsAsFactors = FALSE))
    atomic_save_rds(out, complete_path)
    out
  }, error = function(e) {
    out <- list(
      metrics = data.frame(replicate = rep_id, method = "MGWR",
                           status = "failed", error_message = conditionMessage(e)),
      by_predictor = data.frame(),
      status = data.frame(replicate = rep_id, status = "failed",
                          error_message = conditionMessage(e),
                          stringsAsFactors = FALSE)
    )
    atomic_write_csv(out$status, file.path(rep_dir, "mgwr_status.csv"))
    atomic_save_rds(out, complete_path)
    out
  })
  atomic_write_csv(result$status, file.path(rep_dir, "mgwr_status.csv"))
  result
}

replicates <- 1:100
tasks <- as.list(replicates)
results <- if (workers > 1L) {
  parallel::mclapply(tasks, run_one_rep, mc.cores = workers)
} else {
  lapply(tasks, run_one_rep)
}

complete_files <- file.path(out_root, "fits", sprintf("rep_%03d", replicates),
                            "complete.rds")
objects <- lapply(complete_files[file.exists(complete_files)], readRDS)
nonempty_bind <- function(objects, field) {
  pieces <- lapply(objects, `[[`, field)
  pieces <- pieces[vapply(pieces, nrow, integer(1)) > 0L]
  if (!length(pieces)) data.frame() else do.call(rbind, pieces)
}
metrics <- nonempty_bind(objects, "metrics")
by_predictor <- nonempty_bind(objects, "by_predictor")
status <- nonempty_bind(objects, "status")

summarize_vec <- function(x) {
  x <- x[is.finite(x)]
  c(mean = mean(x), sd = sd(x), mcse = sd(x) / sqrt(length(x)),
    n = length(x))
}

metric_names <- c("mspe", "rmse", "mae", "bias", "beta_mise",
                  "mise_global_only", "mise_spatial_only",
                  "mise_global_plus_spatial", "mise_null",
                  "theta_mse", "u_mise_x3_x6", "beta0_hat",
                  "runtime_sec", "postprocess_runtime_sec",
                  "total_runtime_sec")
summary <- do.call(rbind, lapply(metric_names, function(nm) {
  vals <- summarize_vec(metrics[[nm]])
  data.frame(method = "MGWR", metric = nm,
             mean = vals["mean"], sd = vals["sd"], mcse = vals["mcse"],
             n = vals["n"], stringsAsFactors = FALSE)
}))

predictor_summary <- do.call(rbind, lapply(split(by_predictor,
                                                 by_predictor$predictor),
                                           function(x) {
  data.frame(
    method = "MGWR",
    predictor = x$predictor[1],
    theta_hat_mean = mean(x$theta_hat),
    theta_hat_sd = sd(x$theta_hat),
    theta_hat_mcse = sd(x$theta_hat) / sqrt(nrow(x)),
    theta_true = x$theta_true[1],
    beta_mise_mean = mean(x$beta_mise),
    beta_mise_sd = sd(x$beta_mise),
    beta_mise_mcse = sd(x$beta_mise) / sqrt(nrow(x)),
    u_mise_mean = mean(x$u_mise),
    u_mise_sd = sd(x$u_mise),
    u_mise_mcse = sd(x$u_mise) / sqrt(nrow(x)),
    n = nrow(x),
    stringsAsFactors = FALSE
  )
}))
predictor_summary <- predictor_summary[order(as.integer(sub("X", "",
                                                            predictor_summary$predictor))), ]

existing_path <- file.path(root, "results",
                           "beta0_0p5_100rep_metric_summary_with_wsssgl_posthoc_theta_u.csv")
if (!file.exists(existing_path)) {
  existing_path <- file.path(root, "results", "beta0_0p5_100rep_metric_summary.csv")
}
existing_summary <- if (file.exists(existing_path)) {
  read.csv(existing_path, stringsAsFactors = FALSE)
} else {
  data.frame()
}

atomic_write_csv(metrics, file.path(out_root, "results",
                                    "mgwr_100rep_metrics_by_replicate.csv"))
atomic_write_csv(by_predictor, file.path(out_root, "results",
                                         "mgwr_100rep_by_predictor.csv"))
atomic_write_csv(status, file.path(out_root, "results",
                                   "mgwr_100rep_status.csv"))
atomic_write_csv(summary, file.path(out_root, "results",
                                    "mgwr_100rep_metric_summary.csv"))
atomic_write_csv(predictor_summary, file.path(out_root, "results",
                                              "mgwr_100rep_predictor_summary.csv"))

fmt <- function(mean, sd, mcse) sprintf("%.4f (%.4f; %.4f)", mean, sd, mcse)
metric_wide <- data.frame(metric = metric_names)
metric_wide$MGWR <- vapply(metric_names, function(nm) {
  x <- summary[summary$metric == nm, ]
  fmt(x$mean, x$sd, x$mcse)
}, character(1))

completed <- sum(status$status == "success")
failed <- sum(status$status == "failed")
report <- c(
  "# MGWR extension for beta0=0.5 four-function 100-replicate simulation",
  "",
  paste0("- Source root: `", root, "`"),
  paste0("- Output root: `", out_root, "`"),
  paste0("- Completed MGWR fits: ", completed, " / 100"),
  paste0("- Failed MGWR fits: ", failed),
  "- Formula: `y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10`.",
  "- R formula intercept flag is 1; fitted `SDF` contains an `Intercept` coefficient column.",
  "- MGWR fit runtime includes `GWmodel::gwr.multiscale()` bandwidth selection/calibration.",
  "- Prediction uses interpolated local intercept plus interpolated local coefficients.",
  "- Theta/u summaries are post-hoc decompositions from the reconstructed MGWR beta surfaces.",
  "",
  "## MGWR Metric Summary: mean (SD; MCSE)",
  "",
  paste(capture.output(print(metric_wide, row.names = FALSE)), collapse = "\n"),
  "",
  "## MGWR Predictor Summary",
  "",
  paste(capture.output(print(predictor_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Output files",
  "",
  paste0("- `", file.path(out_root, "results",
                         "mgwr_100rep_metrics_by_replicate.csv"), "`"),
  paste0("- `", file.path(out_root, "results",
                         "mgwr_100rep_metric_summary.csv"), "`"),
  paste0("- `", file.path(out_root, "results",
                         "mgwr_100rep_predictor_summary.csv"), "`")
)
writeLines(report, file.path(out_root, "reports",
                             "mgwr_100rep_extension_report.md"))
cat(paste(report, collapse = "\n"))
