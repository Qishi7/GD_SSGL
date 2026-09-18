comparison_root <- normalizePath(
  file.path(getwd(), "new_method_comparison"),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(comparison_root, "R", "methods", "accelerated",
                 "newssgl_intercept.R"))

out_root <- file.path(comparison_root, "results_debug",
                      "newssgl_intercept_minimal_test")
dir.create(file.path(out_root, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "reports"), recursive = TRUE, showWarnings = FALSE)

empty_selection_label <- function(pip, threshold = 0.5) {
  selected <- which(pip >= threshold)
  if (!length(selected)) return("")
  paste0("X", selected, collapse = ",")
}

run_intercept_only_case <- function() {
  set.seed(20260816)
  n_train <- 180
  n_test <- 80
  p <- 10
  true_beta0 <- -1.75
  sigma <- 0.15
  coords <- cbind(runif(n_train + n_test), runif(n_train + n_test))
  X_raw <- matrix(rnorm((n_train + n_test) * p), n_train + n_test, p)
  x_mean <- colMeans(X_raw[seq_len(n_train), , drop = FALSE])
  x_sd <- apply(X_raw[seq_len(n_train), , drop = FALSE], 2, sd)
  X <- sweep(sweep(X_raw, 2, x_mean, "-"), 2, x_sd, "/")
  y <- true_beta0 + rnorm(n_train + n_test, sd = sigma)
  grid_coords <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 20),
    y = seq(0, 1, length.out = 20)
  ))
  dat0 <- list(
    train = list(y = y[seq_len(n_train)],
                 X = X[seq_len(n_train), , drop = FALSE],
                 coords = coords[seq_len(n_train), , drop = FALSE]),
    test = list(y = y[n_train + seq_len(n_test)],
                X = X[n_train + seq_len(n_test), , drop = FALSE],
                coords = coords[n_train + seq_len(n_test), , drop = FALSE]),
    grid = list(coords = grid_coords,
                true_beta = matrix(0, nrow(grid_coords), p),
                true_u = matrix(0, nrow(grid_coords), p)),
    truth = list(beta0 = true_beta0, theta = rep(0, p),
                 spatial_deviation = rep(FALSE, p))
  )
  fit0 <- fit_newssgl_intercept_fast(
    dat0$train, dat0$test, dat0$grid,
    basis_config = list(n_basis = 4),
    model_config = list(lambda0 = 20, lambda1 = 2,
                        a_sigma = 0.5, b_sigma = var(dat0$train$y) / 2,
                        a_theta = 1, b_theta = 1, a_gamma = 1, b_gamma = 10),
    mcmc_config = list(n_iter = 700, burn_in = 250),
    seed = 260816
  )
  data.frame(
    case = "intercept_only_null",
    true_beta0 = true_beta0,
    beta0_mean = fit0$beta0_mean,
    beta0_abs_error = abs(fit0$beta0_mean - true_beta0),
    mspe = mean((dat0$test$y - fit0$pred_test)^2),
    noise_variance = sigma^2,
    theta_mse = mean(fit0$theta_mean^2),
    beta_mise = mean((fit0$beta_hat_grid - dat0$grid$true_beta)^2),
    max_rb_pip = max(fit0$diagnostics$rb_pip),
    selected = empty_selection_label(fit0$diagnostics$rb_pip),
    runtime_sec = fit0$runtime
  )
}

intercept_only_metrics <- run_intercept_only_case()

set.seed(20260815)
true_beta0 <- 2.25
dat <- simulate_p10_pilot_split(
  n_train = 220,
  n_test = 100,
  grid_size = 30,
  sigma = 0.20,
  u_target_norm = 0.8,
  seed = 20260815
)

mu_train_no_intercept <- dat$train$y -
  as.numeric(rnorm(length(dat$train$y), sd = 0))
dat$train$y <- dat$train$y + true_beta0
dat$test$y <- dat$test$y + true_beta0
dat$truth$beta0 <- true_beta0

prepared <- prepare_newssgl_intercept(
  dat$train,
  basis_config = list(n_basis = 4),
  model_config = list(lambda0 = 20, lambda1 = 2)
)

fit <- fit_newssgl_intercept_fast(
  dat$train, dat$test, dat$grid,
  basis_config = list(n_basis = 4),
  model_config = list(lambda0 = 20, lambda1 = 2,
                      a_sigma = 0.5, b_sigma = var(dat$train$y) / 2,
                      a_theta = 1, b_theta = 1, a_gamma = 1, b_gamma = 10),
  mcmc_config = list(n_iter = 900, burn_in = 300),
  seed = 260815
)

phi_test <- apply_basis_metadata(dat$test$coords, fit$config$basis)
beta_test <- surface_from_components(
  fit$theta_mean, fit$alpha_mean, phi_test, ncol(dat$train$X)
)
pred_identity <- fit$beta0_mean + rowSums(dat$test$X * beta_test)

beta_mse <- colMeans((fit$beta_hat_grid - dat$grid$true_beta)^2)
u_mse <- colMeans((fit$u_hat_grid - dat$grid$true_u)^2)
theta_error <- fit$theta_mean - dat$truth$theta

checks <- data.frame(
  check = c(
    "design_has_intercept_plus_X_plus_Z",
    "centered_basis_training_col_means",
    "prediction_identity_uses_beta0",
    "beta0_absolute_error_below_0p35",
    "theta_mse_below_0p10",
    "no_missing_outputs"
  ),
  value = c(
    prepared$q == 1 + prepared$p + prepared$p * prepared$h,
    prepared$center_error < 1e-8,
    max(abs(fit$pred_test - pred_identity)) < 1e-10,
    abs(fit$beta0_mean - true_beta0) < 0.35,
    mean(theta_error^2) < 0.10,
    all(is.finite(c(fit$beta0_mean, fit$theta_mean, fit$pred_test,
                    fit$beta_hat_grid, fit$u_hat_grid)))
  ),
  numeric_value = c(
    prepared$q,
    prepared$center_error,
    max(abs(fit$pred_test - pred_identity)),
    abs(fit$beta0_mean - true_beta0),
    mean(theta_error^2),
    NA_real_
  )
)

theta_table <- data.frame(
  predictor = paste0("X", seq_along(dat$truth$theta)),
  true_theta = dat$truth$theta,
  theta_mean = fit$theta_mean,
  theta_error = theta_error,
  rb_pip = fit$diagnostics$rb_pip,
  sampled_pip = fit$pip,
  beta_mse = beta_mse,
  u_mse = u_mse
)

metrics <- data.frame(
  method = fit$method_name,
  true_beta0 = true_beta0,
  beta0_mean = fit$beta0_mean,
  beta0_error = fit$beta0_mean - true_beta0,
  mspe = mean((dat$test$y - fit$pred_test)^2),
  theta_mse = mean(theta_error^2),
  beta_mise = mean(beta_mse),
  u_mise_x3_x6 = mean(u_mse[3:6]),
  runtime_sec = fit$runtime,
  selected = empty_selection_label(fit$diagnostics$rb_pip)
)

write.csv(checks, file.path(out_root, "data", "intercept_minimal_checks.csv"),
          row.names = FALSE)
write.csv(theta_table, file.path(out_root, "data", "intercept_theta_pip_table.csv"),
          row.names = FALSE)
write.csv(metrics, file.path(out_root, "data", "intercept_minimal_metrics.csv"),
          row.names = FALSE)
write.csv(intercept_only_metrics,
          file.path(out_root, "data", "intercept_only_null_metrics.csv"),
          row.names = FALSE)
saveRDS(list(data = dat, prepared = prepared, fit = fit, checks = checks,
             theta_table = theta_table, metrics = metrics),
        file.path(out_root, "data", "intercept_minimal_test_fit.rds"))

report <- c(
  "# New SSGL Explicit-Intercept Minimal Test",
  "",
  paste0("- Output root: `", out_root, "`"),
  "- Scope: one small p=10 diagnostic dataset, one GD-SSGL intercept chain.",
  "- Model prediction checked as beta0 + rowSums(X * beta(s)).",
  "- Beta surfaces are reconstructed from theta + Bc alpha only; beta0 is not included in beta_j(s).",
  "",
  "## Intercept-Only Null Check",
  paste(capture.output(print(intercept_only_metrics, row.names = FALSE)),
        collapse = "\n"),
  "",
  "## Metrics",
  paste(capture.output(print(metrics, row.names = FALSE)), collapse = "\n"),
  "",
  "## Checks",
  paste(capture.output(print(checks, row.names = FALSE)), collapse = "\n"),
  "",
  "## Theta and PIP",
  paste(capture.output(print(theta_table, row.names = FALSE)), collapse = "\n")
)
writeLines(report, file.path(out_root, "reports",
                             "newssgl_intercept_minimal_test_report.md"))

cat("Output root:", out_root, "\n")
cat("All checks passed:", all(checks$value), "\n")
cat("Intercept-only null check:\n")
print(intercept_only_metrics, row.names = FALSE)
cat("Metrics:\n")
print(metrics, row.names = FALSE)
cat("Checks:\n")
print(checks, row.names = FALSE)
cat("Theta/PIP:\n")
print(theta_table, row.names = FALSE)
