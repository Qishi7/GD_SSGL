comparison_root <- normalizePath(
  file.path(getwd(), "new_method_comparison"),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(comparison_root, "R", "methods", "accelerated",
                 "newssgl_intercept.R"))

out_root <- file.path(comparison_root, "results_debug",
                      "intercept_comparators_minimal")
dir.create(file.path(out_root, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "reports"), recursive = TRUE, showWarnings = FALSE)

set.seed(20260817)
true_beta0 <- 1.5
dat <- simulate_p10_pilot_split(
  n_train = 350,
  n_test = 120,
  grid_size = 20,
  sigma = 0.25,
  u_target_norm = 1.5,
  seed = 20260817
)
dat$train$y <- dat$train$y + true_beta0
dat$test$y <- dat$test$y + true_beta0
dat$truth$beta0 <- true_beta0

basis_config <- list(n_basis = 6, full_rank_centered = TRUE,
                     full_rank_method = "svd")
mcmc_config <- list(n_iter = 2500, burn_in = 700)
blasso_config <- list(a_sigma = 0.5, b_sigma = var(dat$train$y) / 2,
                      a_theta = 1, b_theta = 1)
svc_config <- list(a_sigma = 0.5, b_sigma = var(dat$train$y) / 2,
                   a_theta = 1, b_theta = 1, kappa2_alpha = 8)

center_data <- dat
y_mean <- mean(center_data$train$y)
center_data$train$y <- center_data$train$y - y_mean
center_data$test$y <- center_data$test$y - y_mean

fits <- list(
  blasso_raw = fit_global_only_blasso(
    dat$train, dat$test, dat$grid, list(), blasso_config, mcmc_config,
    seed = 101
  ),
  blasso_center_add_mean = fit_global_only_blasso(
    center_data$train, center_data$test, center_data$grid, list(),
    blasso_config, mcmc_config, seed = 102
  ),
  blasso_intercept = fit_global_only_blasso_intercept(
    dat$train, dat$test, dat$grid, list(), blasso_config, mcmc_config,
    seed = 103
  ),
  svc_raw = fit_full_svc_no_selection(
    dat$train, dat$test, dat$grid, list(n_basis = 6),
    svc_config, mcmc_config, seed = 201
  ),
  svc_center_add_mean = fit_full_svc_no_selection(
    center_data$train, center_data$test, center_data$grid, list(n_basis = 6),
    svc_config, mcmc_config, seed = 202
  ),
  svc_intercept = fit_full_svc_no_selection_intercept(
    dat$train, dat$test, dat$grid, basis_config,
    svc_config, mcmc_config, seed = 203
  )
)
fits$blasso_center_add_mean$pred_test <- fits$blasso_center_add_mean$pred_test +
  y_mean
fits$blasso_center_add_mean$diagnostics$y_mean_added_back <- y_mean
fits$svc_center_add_mean$pred_test <- fits$svc_center_add_mean$pred_test +
  y_mean
fits$svc_center_add_mean$diagnostics$y_mean_added_back <- y_mean

method_kind <- c(
  blasso_raw = "Bayesian Lasso old raw-y",
  blasso_center_add_mean = "Bayesian Lasso old center-y + add mean",
  blasso_intercept = "Bayesian Lasso explicit beta0",
  svc_raw = "Gaussian SVC old raw-y",
  svc_center_add_mean = "Gaussian SVC old center-y + add mean",
  svc_intercept = "Gaussian SVC explicit beta0"
)

prediction_identity <- function(name, fit) {
  if (name == "blasso_intercept") {
    return(max(abs(fit$pred_test -
      (fit$beta0_mean + as.vector(dat$test$X %*% fit$theta_mean)))))
  }
  if (name == "svc_intercept") {
    phi_test <- apply_centered_basis_metadata_fullrank(
      dat$test$coords, fit$config$basis
    )
    beta_test <- surface_from_components(
      fit$theta_mean, fit$alpha_mean, phi_test, ncol(dat$train$X)
    )
    return(max(abs(fit$pred_test -
      (fit$beta0_mean + rowSums(dat$test$X * beta_test)))))
  }
  NA_real_
}

metrics <- do.call(rbind, lapply(names(fits), function(nm) {
  fit <- fits[[nm]]
  beta_mse <- colMeans((fit$beta_hat_grid - dat$grid$true_beta)^2)
  u_mse <- if (is.null(fit$u_hat_grid)) rep(NA_real_, ncol(dat$train$X)) else {
    colMeans((fit$u_hat_grid - dat$grid$true_u)^2)
  }
  data.frame(
    method = nm,
    label = method_kind[[nm]],
    beta0_type = if (!is.null(fit$beta0_mean)) {
      "explicit p(beta0) proportional to 1"
    } else if (!is.null(fit$diagnostics$y_mean_added_back)) {
      "train y mean added after centered fit"
    } else {
      "none"
    },
    true_beta0 = true_beta0,
    beta0_estimate = fit$beta0_mean %||%
      fit$diagnostics$y_mean_added_back %||% NA_real_,
    beta0_error = (fit$beta0_mean %||%
      fit$diagnostics$y_mean_added_back %||% NA_real_) - true_beta0,
    mspe = mean((dat$test$y - fit$pred_test)^2),
    theta_mse = mean((fit$theta_mean - dat$truth$theta)^2),
    beta_mise = mean(beta_mse),
    beta_mise_global_only = mean(beta_mse[1:2]),
    beta_mise_spatial_only = mean(beta_mse[3:4]),
    beta_mise_global_spatial = mean(beta_mse[5:6]),
    beta_mise_null = mean(beta_mse[7:10]),
    u_mise_x3_x6 = mean(u_mse[3:6], na.rm = TRUE),
    sigma2_mean = mean(fit$diagnostics$sigma2_draws %||% NA_real_,
                       na.rm = TRUE),
    runtime_sec = fit$runtime,
    prediction_identity_max_abs_error = prediction_identity(nm, fit),
    raw_basis_dimension = fit$diagnostics$raw_basis_dimension %||% NA_integer_,
    effective_basis_dimension =
      fit$diagnostics$effective_basis_dimension %||% NA_integer_
  )
}))

theta_summary <- do.call(rbind, lapply(names(fits), function(nm) {
  fit <- fits[[nm]]
  data.frame(
    method = nm,
    predictor = paste0("X", seq_along(fit$theta_mean)),
    true_theta = dat$truth$theta,
    theta_mean = fit$theta_mean,
    theta_error = fit$theta_mean - dat$truth$theta
  )
}))

checks <- data.frame(
  check = c(
    "blasso_explicit_prediction_identity",
    "svc_explicit_prediction_identity",
    "blasso_beta0_finite",
    "svc_beta0_finite",
    "svc_full_rank_effective_h_35"
  ),
  value = c(
    metrics$prediction_identity_max_abs_error[
      metrics$method == "blasso_intercept"] < 1e-10,
    metrics$prediction_identity_max_abs_error[
      metrics$method == "svc_intercept"] < 1e-10,
    is.finite(fits$blasso_intercept$beta0_mean),
    is.finite(fits$svc_intercept$beta0_mean),
    fits$svc_intercept$diagnostics$effective_basis_dimension == 35
  )
)

write.csv(metrics, file.path(out_root, "data",
                             "intercept_comparator_metrics.csv"),
          row.names = FALSE)
write.csv(theta_summary, file.path(out_root, "data",
                                   "intercept_comparator_theta_summary.csv"),
          row.names = FALSE)
write.csv(checks, file.path(out_root, "data",
                            "intercept_comparator_checks.csv"),
          row.names = FALSE)
saveRDS(list(data = dat, fits = fits, metrics = metrics,
             theta_summary = theta_summary, checks = checks),
        file.path(out_root, "data", "intercept_comparator_test_objects.rds"))

fmt <- function(x) paste(capture.output(print(x, row.names = FALSE)),
                         collapse = "\n")
report <- c(
  "# Bayesian Lasso and Gaussian SVC Explicit-Intercept Comparator Test",
  "",
  paste0("- Output root: `", out_root, "`"),
  "- Scope: one small nonzero-intercept p=10 simulation.",
  "- New Bayesian Lasso comparator uses an unpenalized scalar beta0 with p(beta0) proportional to 1.",
  "- New Gaussian SVC comparator uses an unpenalized scalar beta0 with p(beta0) proportional to 1; other priors unchanged.",
  "- Gaussian SVC uses the same SVD full-rank centered basis representation for centered spatial basis.",
  "",
  "## Checks",
  fmt(checks),
  "",
  "## Metrics",
  fmt(metrics),
  "",
  "## Theta Summary",
  fmt(theta_summary)
)
writeLines(report, file.path(out_root, "reports",
                             "intercept_comparator_test_report.md"))

cat("Output root:", out_root, "\n")
cat("All checks passed:", all(checks$value), "\n")
cat("\nChecks:\n")
print(checks, row.names = FALSE)
cat("\nMetrics:\n")
print(metrics, row.names = FALSE)
cat("\nReport:", file.path(out_root, "reports",
                          "intercept_comparator_test_report.md"), "\n")
