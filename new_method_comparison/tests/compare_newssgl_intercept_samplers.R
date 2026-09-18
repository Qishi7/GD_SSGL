comparison_root <- normalizePath(
  file.path(getwd(), "new_method_comparison"),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(comparison_root, "R", "methods", "accelerated",
                 "newssgl_intercept.R"))

out_root <- file.path(comparison_root, "results_debug",
                      "newssgl_intercept_sampler_comparison")
dir.create(file.path(out_root, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "reports"), recursive = TRUE, showWarnings = FALSE)

set.seed(20260816)
true_beta0 <- 1.75
dat <- simulate_p10_pilot_split(
  n_train = 320,
  n_test = 140,
  grid_size = 30,
  sigma = 0.25,
  u_target_norm = 1.2,
  seed = 20260816
)
dat$train$y <- dat$train$y + true_beta0
dat$test$y <- dat$test$y + true_beta0
dat$truth$beta0 <- true_beta0

basis_config <- list(n_basis = 4)
model_config <- list(
  lambda0 = 20,
  lambda1 = 2,
  a_sigma = 0.5,
  b_sigma = var(dat$train$y) / 2,
  a_theta = 1,
  b_theta = 1,
  a_gamma = 1,
  b_gamma = 10,
  initialization = list(
    beta0 = mean(dat$train$y),
    theta = rep(0, ncol(dat$train$X)),
    alpha = rep(0, ncol(dat$train$X) * basis_config$n_basis^2),
    gamma = rep(0, ncol(dat$train$X)),
    sigma2 = var(dat$train$y),
    lambda_theta2 = 1,
    pi_gamma = 0.5
  )
)
mcmc_config <- list(n_iter = 3500, burn_in = 1000)

prepared <- prepare_newssgl_intercept(dat$train, basis_config, model_config)
W_nonintercept <- prepared$W[, -1, drop = FALSE]
y_centered <- dat$train$y - mean(dat$train$y)
W_centered <- sweep(W_nonintercept, 2, colMeans(W_nonintercept), "-")
projection_checks <- data.frame(
  check = c(
    "direct_W_has_intercept",
    "nonintercept_centered_design_col_means",
    "centered_y_mean",
    "basis_center_error"
  ),
  value = c(
    ncol(prepared$W) == 1 + prepared$p + prepared$p * prepared$h,
    max(abs(colMeans(W_centered))) < 1e-12,
    abs(mean(y_centered)) < 1e-12,
    prepared$center_error < 1e-8
  ),
  numeric_value = c(
    ncol(prepared$W),
    max(abs(colMeans(W_centered))),
    abs(mean(y_centered)),
    prepared$center_error
  )
)

direct <- fit_newssgl_intercept_fast(
  dat$train, dat$test, dat$grid,
  basis_config = basis_config,
  model_config = model_config,
  mcmc_config = mcmc_config,
  seed = 2026081601
)

collapsed <- fit_newssgl_intercept_collapsed_fast(
  dat$train, dat$test, dat$grid,
  basis_config = basis_config,
  model_config = model_config,
  mcmc_config = mcmc_config,
  seed = 2026081602
)

summarize_fit <- function(fit, label) {
  beta_mse <- colMeans((fit$beta_hat_grid - dat$grid$true_beta)^2)
  data.frame(
    sampler = label,
    beta0_mean = fit$beta0_mean,
    beta0_error = fit$beta0_mean - dat$truth$beta0,
    theta_mse = mean((fit$theta_mean - dat$truth$theta)^2),
    beta_mise = mean(beta_mse),
    u_mise_x3_x6 = mean(colMeans((fit$u_hat_grid - dat$grid$true_u)^2)[3:6]),
    mspe = mean((dat$test$y - fit$pred_test)^2),
    rb_pip_x3 = fit$diagnostics$rb_pip[3],
    rb_pip_x4 = fit$diagnostics$rb_pip[4],
    rb_pip_x5 = fit$diagnostics$rb_pip[5],
    rb_pip_x6 = fit$diagnostics$rb_pip[6],
    max_rb_pip_null = max(fit$diagnostics$rb_pip[c(1, 2, 7:10)]),
    sigma2_mean = mean(fit$diagnostics$sigma2_draws),
    pi_gamma_mean = mean(fit$diagnostics$pi_gamma_draws),
    runtime_sec = fit$runtime
  )
}

fit_summary <- rbind(
  summarize_fit(direct, "direct_beta0_sampled"),
  summarize_fit(collapsed, "collapsed_beta0_integrated")
)

theta_compare <- data.frame(
  predictor = paste0("X", seq_along(direct$theta_mean)),
  true_theta = dat$truth$theta,
  direct_theta = direct$theta_mean,
  collapsed_theta = collapsed$theta_mean,
  difference_collapsed_minus_direct = collapsed$theta_mean - direct$theta_mean,
  direct_rb_pip = direct$diagnostics$rb_pip,
  collapsed_rb_pip = collapsed$diagnostics$rb_pip,
  rb_pip_difference = collapsed$diagnostics$rb_pip - direct$diagnostics$rb_pip
)

surface_diff <- collapsed$beta_hat_grid - direct$beta_hat_grid
u_diff <- collapsed$u_hat_grid - direct$u_hat_grid
pred_diff <- collapsed$pred_test - direct$pred_test
alpha_diff <- collapsed$alpha_mean - direct$alpha_mean
beta0_draw_sd <- sd(direct$beta0_draws)

difference_summary <- data.frame(
  quantity = c(
    "beta0_mean",
    "theta_mean_max_abs_diff",
    "theta_mean_rmse_diff",
    "alpha_mean_max_abs_diff",
    "alpha_mean_rmse_diff",
    "rb_pip_max_abs_diff",
    "rb_pip_rmse_diff",
    "beta_grid_max_abs_diff",
    "beta_grid_rmse_diff",
    "u_grid_max_abs_diff",
    "u_grid_rmse_diff",
    "pred_test_max_abs_diff",
    "pred_test_rmse_diff",
    "sigma2_mean",
    "pi_gamma_mean",
    "runtime_ratio_collapsed_over_direct"
  ),
  direct = c(
    direct$beta0_mean,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    mean(direct$diagnostics$sigma2_draws),
    mean(direct$diagnostics$pi_gamma_draws),
    NA
  ),
  collapsed = c(
    collapsed$beta0_mean,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    mean(collapsed$diagnostics$sigma2_draws),
    mean(collapsed$diagnostics$pi_gamma_draws),
    NA
  ),
  absolute_difference = c(
    abs(collapsed$beta0_mean - direct$beta0_mean),
    max(abs(collapsed$theta_mean - direct$theta_mean)),
    sqrt(mean((collapsed$theta_mean - direct$theta_mean)^2)),
    max(abs(alpha_diff)),
    sqrt(mean(alpha_diff^2)),
    max(abs(collapsed$diagnostics$rb_pip - direct$diagnostics$rb_pip)),
    sqrt(mean((collapsed$diagnostics$rb_pip - direct$diagnostics$rb_pip)^2)),
    max(abs(surface_diff)),
    sqrt(mean(surface_diff^2)),
    max(abs(u_diff)),
    sqrt(mean(u_diff^2)),
    max(abs(pred_diff)),
    sqrt(mean(pred_diff^2)),
    abs(mean(collapsed$diagnostics$sigma2_draws) -
          mean(direct$diagnostics$sigma2_draws)),
    abs(mean(collapsed$diagnostics$pi_gamma_draws) -
          mean(direct$diagnostics$pi_gamma_draws)),
    collapsed$runtime / direct$runtime
  )
)

mc_reference <- data.frame(
  quantity = c(
    "direct_beta0_draw_sd",
    "direct_beta0_mcse",
    "collapsed_beta0_cond_mean_sd",
    "collapsed_beta0_cond_mean_mcse",
    "n_kept_draws"
  ),
  value = c(
    beta0_draw_sd,
    beta0_draw_sd / sqrt(length(direct$beta0_draws)),
    sd(collapsed$beta0_cond_mean_draws),
    sd(collapsed$beta0_cond_mean_draws) /
      sqrt(length(collapsed$beta0_cond_mean_draws)),
    length(direct$beta0_draws)
  )
)

equivalence_checks <- data.frame(
  check = c(
    "projection_checks_pass",
    "beta0_means_within_3_direct_mcse_plus_0p02",
    "theta_max_abs_diff_below_0p10",
    "rb_pip_max_abs_diff_below_0p10",
    "prediction_rmse_diff_below_0p10",
    "beta_grid_rmse_diff_below_0p10",
    "same_selected_set_at_0p5"
  ),
  value = c(
    all(projection_checks$value),
    abs(collapsed$beta0_mean - direct$beta0_mean) <
      3 * beta0_draw_sd / sqrt(length(direct$beta0_draws)) + 0.02,
    max(abs(collapsed$theta_mean - direct$theta_mean)) < 0.10,
    max(abs(collapsed$diagnostics$rb_pip - direct$diagnostics$rb_pip)) < 0.10,
    sqrt(mean(pred_diff^2)) < 0.10,
    sqrt(mean(surface_diff^2)) < 0.10,
    identical(which(direct$diagnostics$rb_pip >= 0.5),
              which(collapsed$diagnostics$rb_pip >= 0.5))
  )
)

write.csv(projection_checks,
          file.path(out_root, "data", "sampler_projection_checks.csv"),
          row.names = FALSE)
write.csv(fit_summary,
          file.path(out_root, "data", "sampler_fit_summary.csv"),
          row.names = FALSE)
write.csv(theta_compare,
          file.path(out_root, "data", "sampler_theta_pip_comparison.csv"),
          row.names = FALSE)
write.csv(difference_summary,
          file.path(out_root, "data", "sampler_difference_summary.csv"),
          row.names = FALSE)
write.csv(mc_reference,
          file.path(out_root, "data", "sampler_mc_reference.csv"),
          row.names = FALSE)
write.csv(equivalence_checks,
          file.path(out_root, "data", "sampler_equivalence_checks.csv"),
          row.names = FALSE)
saveRDS(list(data = dat, direct = direct, collapsed = collapsed,
             projection_checks = projection_checks,
             fit_summary = fit_summary,
             theta_compare = theta_compare,
             difference_summary = difference_summary,
             mc_reference = mc_reference,
             equivalence_checks = equivalence_checks),
        file.path(out_root, "data",
                  "direct_vs_collapsed_intercept_sampler_comparison.rds"))

report <- c(
  "# Direct vs Collapsed Explicit-Intercept GD-SSGL Sampler Check",
  "",
  paste0("- Output root: `", out_root, "`"),
  "- Scope: one small p=10 diagnostic dataset.",
  "- Direct sampler: samples beta0 jointly with theta and alpha.",
  "- Collapsed sampler: integrates beta0 out by centering y and all non-intercept design columns; beta0 is not sampled.",
  "- Prediction for collapsed sampler uses the posterior mean of the conditional beta0 mean.",
  "",
  "## Projection Checks",
  paste(capture.output(print(projection_checks, row.names = FALSE)), collapse = "\n"),
  "",
  "## Fit Summary",
  paste(capture.output(print(fit_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Difference Summary",
  paste(capture.output(print(difference_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## MC Reference",
  paste(capture.output(print(mc_reference, row.names = FALSE)), collapse = "\n"),
  "",
  "## Equivalence Checks",
  paste(capture.output(print(equivalence_checks, row.names = FALSE)), collapse = "\n"),
  "",
  "## Theta/PIP Comparison",
  paste(capture.output(print(theta_compare, row.names = FALSE)), collapse = "\n")
)
writeLines(report, file.path(out_root, "reports",
                             "direct_vs_collapsed_sampler_comparison_report.md"))

cat("Output root:", out_root, "\n")
cat("Projection checks pass:", all(projection_checks$value), "\n")
cat("Equivalence checks pass:", all(equivalence_checks$value), "\n")
cat("\nFit summary:\n")
print(fit_summary, row.names = FALSE)
cat("\nDifference summary:\n")
print(difference_summary, row.names = FALSE)
cat("\nEquivalence checks:\n")
print(equivalence_checks, row.names = FALSE)
cat("\nTheta/PIP comparison:\n")
print(theta_compare, row.names = FALSE)
