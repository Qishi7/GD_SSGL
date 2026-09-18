project_root <- normalizePath(Sys.getenv("GDSSGL_ROOT", getwd()), winslash = "/", mustWork = TRUE)
root <- file.path(project_root, "new_method_comparison", "formal_simulation", "main_simulation", "four_function_beta0_0p5_100rep_perrep_tuned_intercept_svd")

results_dir <- file.path(root, "results")
reports_dir <- file.path(root, "reports")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)

theta_true <- c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0)
groups <- list(
  global_only = 1:2,
  spatial_only = 3:4,
  global_plus_spatial = 5:6,
  null = 7:10
)

mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

summarize_vec <- function(x) {
  x <- x[is.finite(x)]
  data.frame(
    n = length(x),
    mean = if (length(x)) mean(x) else NA_real_,
    sd = if (length(x) > 1) stats::sd(x) else NA_real_,
    mcse = mcse(x),
    median = if (length(x)) stats::median(x) else NA_real_,
    q025 = if (length(x)) unname(stats::quantile(x, 0.025)) else NA_real_,
    q975 = if (length(x)) unname(stats::quantile(x, 0.975)) else NA_real_
  )
}

rep_rows <- list()
pred_rows <- list()

for (rep_id in seq_len(100)) {
  rep_dir <- file.path(root, "fits", sprintf("rep_%03d", rep_id))
  fit_path <- file.path(rep_dir, "wsssgl_final_fit.rds")
  data_path <- file.path(rep_dir, "shared_dataset_beta0_0p5.rds")
  if (!file.exists(fit_path) || !file.exists(data_path)) {
    stop("Missing WS-SSGL fit or data for replicate ", rep_id)
  }

  fit <- readRDS(fit_path)
  dat <- readRDS(data_path)

  beta_hat <- fit$beta_grid
  if (is.null(beta_hat)) beta_hat <- fit$beta_hat_grid
  if (is.null(beta_hat)) stop("No beta grid found for replicate ", rep_id)
  if (!all(dim(beta_hat) == dim(dat$grid$true_beta))) {
    stop("Beta grid dimension mismatch for replicate ", rep_id)
  }

  theta_hat <- colMeans(beta_hat)
  u_hat <- sweep(beta_hat, 2, theta_hat, "-")

  theta_err2 <- (theta_hat - theta_true)^2
  u_ise <- colMeans((u_hat - dat$grid$true_u)^2)
  beta_ise <- colMeans((beta_hat - dat$grid$true_beta)^2)

  rep_rows[[rep_id]] <- data.frame(
    replicate = rep_id,
    method = "WS-SSGL",
    decomposition = "posthoc_grid_average",
    theta_mse = mean(theta_err2),
    theta_mse_global_only = mean(theta_err2[groups$global_only]),
    theta_mse_spatial_only = mean(theta_err2[groups$spatial_only]),
    theta_mse_global_plus_spatial = mean(theta_err2[groups$global_plus_spatial]),
    theta_mse_null = mean(theta_err2[groups$null]),
    u_mise_x3_x6 = mean(u_ise[3:6]),
    u_mise_spatial_only = mean(u_ise[groups$spatial_only]),
    u_mise_global_plus_spatial = mean(u_ise[groups$global_plus_spatial]),
    beta_mise_check = mean(beta_ise),
    stringsAsFactors = FALSE
  )

  pred_rows[[rep_id]] <- data.frame(
    replicate = rep_id,
    predictor = paste0("X", seq_len(10)),
    target_theta = theta_true,
    theta_hat = theta_hat,
    theta_error = theta_hat - theta_true,
    theta_sq_error = theta_err2,
    u_mise = u_ise,
    beta_ise = beta_ise,
    stringsAsFactors = FALSE
  )
}

rep_df <- do.call(rbind, rep_rows)
pred_df <- do.call(rbind, pred_rows)

metric_cols <- setdiff(names(rep_df), c("replicate", "method", "decomposition"))
summary_df <- do.call(rbind, lapply(metric_cols, function(nm) {
  cbind(method = "WS-SSGL", metric = nm, summarize_vec(rep_df[[nm]]))
}))
row.names(summary_df) <- NULL

predictor_summary <- do.call(rbind, lapply(split(pred_df, pred_df$predictor), function(z) {
  data.frame(
    predictor = z$predictor[1],
    target_theta = z$target_theta[1],
    theta_hat_mean = mean(z$theta_hat),
    theta_hat_sd = stats::sd(z$theta_hat),
    theta_mse_mean = mean(z$theta_sq_error),
    theta_mse_sd = stats::sd(z$theta_sq_error),
    theta_mse_mcse = mcse(z$theta_sq_error),
    u_mise_mean = mean(z$u_mise),
    u_mise_sd = stats::sd(z$u_mise),
    u_mise_mcse = mcse(z$u_mise),
    beta_ise_mean = mean(z$beta_ise),
    beta_ise_sd = stats::sd(z$beta_ise),
    beta_ise_mcse = mcse(z$beta_ise),
    stringsAsFactors = FALSE
  )
}))
row.names(predictor_summary) <- NULL
predictor_summary <- predictor_summary[order(as.integer(sub("X", "", predictor_summary$predictor))), ]

write.csv(rep_df,
          file.path(results_dir, "beta0_0p5_100rep_wsssgl_posthoc_theta_u_by_replicate.csv"),
          row.names = FALSE)
write.csv(pred_df,
          file.path(results_dir, "beta0_0p5_100rep_wsssgl_posthoc_theta_u_by_predictor_replicate.csv"),
          row.names = FALSE)
write.csv(summary_df,
          file.path(results_dir, "beta0_0p5_100rep_wsssgl_posthoc_theta_u_summary.csv"),
          row.names = FALSE)
write.csv(predictor_summary,
          file.path(results_dir, "beta0_0p5_100rep_wsssgl_posthoc_theta_u_predictor_summary.csv"),
          row.names = FALSE)

metric_path <- file.path(results_dir, "beta0_0p5_100rep_metric_summary.csv")
if (file.exists(metric_path)) {
  metric_summary <- read.csv(metric_path)
  fill_metric <- function(metric_name, value_metric) {
    idx <- metric_summary$method == "WS-SSGL" & metric_summary$metric == metric_name
    src <- summary_df[summary_df$metric == value_metric, ]
    if (sum(idx) == 1L && nrow(src) == 1L) {
      metric_summary$n[idx] <<- src$n
      metric_summary$mean[idx] <<- src$mean
      metric_summary$sd[idx] <<- src$sd
      metric_summary$mcse[idx] <<- src$mcse
    }
  }
  fill_metric("theta_mse", "theta_mse")
  fill_metric("u_mise_x3_x6", "u_mise_x3_x6")
  write.csv(metric_summary,
            file.path(results_dir, "beta0_0p5_100rep_metric_summary_with_wsssgl_posthoc_theta_u.csv"),
            row.names = FALSE)
}

fmt <- function(x) sprintf("%.6f", x)
main_summary <- summary_df[summary_df$metric %in% c("theta_mse", "u_mise_x3_x6",
                                                    "u_mise_spatial_only",
                                                    "u_mise_global_plus_spatial"), ]

report <- c(
  "# WS-SSGL post-hoc theta/u decomposition for beta0=0.5 100-rep run",
  "",
  "This analysis does not refit any model. For each replicate, it loads the saved WS-SSGL whole-surface estimate `beta_grid`, then computes:",
  "",
  "- `theta_hat_j = mean_grid beta_hat_j(s)`",
  "- `u_hat_j(s) = beta_hat_j(s) - theta_hat_j`",
  "",
  "This is a post-hoc grid-average decomposition, analogous to the GAM theta/u post-processing. WS-SSGL itself remains a whole coefficient-surface selection model.",
  "",
  "## Main summary",
  "",
  paste(capture.output(print(main_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Predictor-level summary",
  "",
  paste(capture.output(print(predictor_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Output files",
  "",
  paste0("- ", file.path(results_dir, "beta0_0p5_100rep_wsssgl_posthoc_theta_u_by_replicate.csv")),
  paste0("- ", file.path(results_dir, "beta0_0p5_100rep_wsssgl_posthoc_theta_u_by_predictor_replicate.csv")),
  paste0("- ", file.path(results_dir, "beta0_0p5_100rep_wsssgl_posthoc_theta_u_summary.csv")),
  paste0("- ", file.path(results_dir, "beta0_0p5_100rep_wsssgl_posthoc_theta_u_predictor_summary.csv")),
  paste0("- ", file.path(results_dir, "beta0_0p5_100rep_metric_summary_with_wsssgl_posthoc_theta_u.csv"))
)
writeLines(report, file.path(reports_dir, "beta0_0p5_100rep_wsssgl_posthoc_theta_u_report.md"))

cat("\nWS-SSGL post-hoc theta/u summary:\n")
print(main_summary, row.names = FALSE)
cat("\nPredictor-level summary:\n")
print(predictor_summary, row.names = FALSE)
cat("\nSaved outputs under:\n", results_dir, "\n", sep = "")
