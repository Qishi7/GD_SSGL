comparison_root <- normalizePath(
  file.path(getwd(), "new_method_comparison"),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(comparison_root, "R", "methods", "accelerated",
                 "newssgl_intercept.R"))

out_root <- file.path(comparison_root, "results_debug",
                      "intercept_direct_collapsed_mcmc_svd")
dir.create(file.path(out_root, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "reports"), recursive = TRUE, showWarnings = FALSE)

ess <- function(x) {
  x <- as.numeric(x)
  if (length(unique(x[is.finite(x)])) < 2) return(NA_real_)
  as.numeric(coda::effectiveSize(x))
}

summ <- function(x) {
  c(mean = mean(x), sd = sd(x),
    q025 = unname(quantile(x, 0.025)),
    q975 = unname(quantile(x, 0.975)))
}

surface_draws_for_predictor <- function(theta_draws, alpha_draws, phi, j, p, h) {
  idx <- ((j - 1L) * h + 1L):(j * h)
  sweep(phi %*% alpha_draws[idx, , drop = FALSE], 2,
        theta_draws[j, ], "+")
}

u_draws_for_predictor <- function(alpha_draws, phi, j, h) {
  idx <- ((j - 1L) * h + 1L):(j * h)
  phi %*% alpha_draws[idx, , drop = FALSE]
}

surface_summary <- function(fit, phi_grid, predictors, locations, sampler) {
  p <- length(fit$theta_mean)
  h <- ncol(phi_grid)
  do.call(rbind, lapply(predictors, function(j) {
    beta_draws <- surface_draws_for_predictor(
      fit$theta_draws, fit$diagnostics$alpha_draws, phi_grid, j, p, h
    )
    u_draws <- u_draws_for_predictor(
      fit$diagnostics$alpha_draws, phi_grid, j, h
    )
    do.call(rbind, lapply(seq_along(locations), function(a) {
      loc <- locations[[a]]
      data.frame(
        sampler = sampler,
        predictor = paste0("X", j),
        location = names(locations)[a],
        grid_index = loc,
        component = c("beta", "u"),
        rbind(summ(beta_draws[loc, ]), summ(u_draws[loc, ])),
        row.names = NULL
      )
    }))
  }))
}

surface_aggregate <- function(fit, phi_grid, true_beta, true_u, predictors,
                              sampler) {
  p <- length(fit$theta_mean)
  h <- ncol(phi_grid)
  do.call(rbind, lapply(predictors, function(j) {
    beta_draws <- surface_draws_for_predictor(
      fit$theta_draws, fit$diagnostics$alpha_draws, phi_grid, j, p, h
    )
    u_draws <- u_draws_for_predictor(
      fit$diagnostics$alpha_draws, phi_grid, j, h
    )
    beta_mean <- rowMeans(beta_draws)
    u_mean <- rowMeans(u_draws)
    beta_sd <- apply(beta_draws, 1, sd)
    u_sd <- apply(u_draws, 1, sd)
    beta_ci <- t(apply(beta_draws, 1, quantile, probs = c(0.025, 0.975)))
    u_ci <- t(apply(u_draws, 1, quantile, probs = c(0.025, 0.975)))
    data.frame(
      sampler = sampler,
      predictor = paste0("X", j),
      beta_mise = mean((beta_mean - true_beta[, j])^2),
      u_mise = mean((u_mean - true_u[, j])^2),
      beta_mean_posterior_sd = mean(beta_sd),
      u_mean_posterior_sd = mean(u_sd),
      beta_mean_interval_width = mean(beta_ci[, 2] - beta_ci[, 1]),
      u_mean_interval_width = mean(u_ci[, 2] - u_ci[, 1]),
      beta_pointwise_coverage = mean(beta_ci[, 1] <= true_beta[, j] &
                                       true_beta[, j] <= beta_ci[, 2]),
      u_pointwise_coverage = mean(u_ci[, 1] <= true_u[, j] &
                                    true_u[, j] <= u_ci[, 2])
    )
  }))
}

set.seed(20260817)
true_beta0 <- 1.5
dat <- simulate_p10_pilot_split(
  n_train = 500,
  n_test = 150,
  grid_size = 25,
  sigma = 0.25,
  u_target_norm = 2.0,
  seed = 20260817
)
dat$train$y <- dat$train$y + true_beta0
dat$test$y <- dat$test$y + true_beta0
dat$truth$beta0 <- true_beta0

basis_config <- list(n_basis = 6, full_rank_centered = TRUE,
                     full_rank_method = "svd")
tmp_prepare <- prepare_newssgl_intercept(
  dat$train, basis_config, list(lambda0 = 20, lambda1 = 2)
)
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
    theta = rep(0, 10),
    alpha = rep(0, 10 * tmp_prepare$h),
    gamma = rep(0, 10),
    sigma2 = var(dat$train$y),
    lambda_theta2 = 1,
    pi_gamma = 0.5
  )
)
mcmc_config <- list(n_iter = 10000, burn_in = 3000)

direct <- fit_newssgl_intercept_fast(
  dat$train, dat$test, dat$grid,
  basis_config, model_config, mcmc_config, seed = 17001
)
collapsed <- fit_newssgl_intercept_collapsed_fast(
  dat$train, dat$test, dat$grid,
  basis_config, model_config, mcmc_config, seed = 17002
)

set.seed(17003)
n_train <- length(dat$train$y)
collapsed_beta0_recovered <- collapsed$beta0_cond_mean_draws +
  rnorm(length(collapsed$beta0_cond_mean_draws),
        sd = sqrt(collapsed$diagnostics$sigma2_draws / n_train))

phi_grid <- apply_centered_basis_metadata_fullrank(
  dat$grid$coords, direct$config$basis
)
selected_predictors <- c(1, 3, 5, 6)
surface_predictors <- c(3, 4, 5, 6)
loc_targets <- rbind(c(0.25, 0.25), c(0.5, 0.5), c(0.75, 0.75))
grid_coords <- as.matrix(dat$grid$coords)
loc_ids <- apply(loc_targets, 1, function(z) {
  which.min((grid_coords[, 1] - z[1])^2 + (grid_coords[, 2] - z[2])^2)
})
locations <- as.list(loc_ids)
names(locations) <- c("near_0.25_0.25", "near_0.50_0.50", "near_0.75_0.75")

theta_summary <- do.call(rbind, lapply(seq_len(10), function(j) {
  d <- summ(direct$theta_draws[j, ])
  c <- summ(collapsed$theta_draws[j, ])
  data.frame(
    predictor = paste0("X", j),
    true_theta = dat$truth$theta[j],
    direct_mean = d["mean"],
    collapsed_mean = c["mean"],
    mean_diff = c["mean"] - d["mean"],
    direct_sd = d["sd"],
    collapsed_sd = c["sd"],
    sd_diff = c["sd"] - d["sd"],
    direct_q025 = d["q025"],
    direct_q975 = d["q975"],
    collapsed_q025 = c["q025"],
    collapsed_q975 = c["q975"],
    direct_ess = ess(direct$theta_draws[j, ]),
    collapsed_ess = ess(collapsed$theta_draws[j, ]),
    direct_ess_per_sec = ess(direct$theta_draws[j, ]) / direct$runtime,
    collapsed_ess_per_sec = ess(collapsed$theta_draws[j, ]) /
      collapsed$runtime
  )
}))

beta0_summary <- rbind(
  data.frame(sampler = "direct",
             t(summ(direct$beta0_draws)),
             ess = ess(direct$beta0_draws),
             ess_per_sec = ess(direct$beta0_draws) / direct$runtime),
  data.frame(sampler = "collapsed_recovered",
             t(summ(collapsed_beta0_recovered)),
             ess = ess(collapsed_beta0_recovered),
             ess_per_sec = ess(collapsed_beta0_recovered) /
               collapsed$runtime),
  data.frame(sampler = "collapsed_conditional_mean_only",
             t(summ(collapsed$beta0_cond_mean_draws)),
             ess = ess(collapsed$beta0_cond_mean_draws),
             ess_per_sec = ess(collapsed$beta0_cond_mean_draws) /
               collapsed$runtime)
)
beta0_summary$true_beta0 <- true_beta0

sigma_pip_summary <- rbind(
  data.frame(
    sampler = "direct",
    sigma2_mean = mean(direct$diagnostics$sigma2_draws),
    sigma2_sd = sd(direct$diagnostics$sigma2_draws),
    sigma2_q025 = quantile(direct$diagnostics$sigma2_draws, 0.025),
    sigma2_q975 = quantile(direct$diagnostics$sigma2_draws, 0.975),
    sigma2_ess = ess(direct$diagnostics$sigma2_draws),
    sigma2_ess_per_sec = ess(direct$diagnostics$sigma2_draws) /
      direct$runtime,
    runtime_sec = direct$runtime
  ),
  data.frame(
    sampler = "collapsed",
    sigma2_mean = mean(collapsed$diagnostics$sigma2_draws),
    sigma2_sd = sd(collapsed$diagnostics$sigma2_draws),
    sigma2_q025 = quantile(collapsed$diagnostics$sigma2_draws, 0.025),
    sigma2_q975 = quantile(collapsed$diagnostics$sigma2_draws, 0.975),
    sigma2_ess = ess(collapsed$diagnostics$sigma2_draws),
    sigma2_ess_per_sec = ess(collapsed$diagnostics$sigma2_draws) /
      collapsed$runtime,
    runtime_sec = collapsed$runtime
  )
)

pip_summary <- data.frame(
  predictor = paste0("X", 1:10),
  target_spatial = 1:10 %in% 3:6,
  direct_sampled_pip = direct$pip,
  collapsed_sampled_pip = collapsed$pip,
  sampled_pip_diff = collapsed$pip - direct$pip,
  direct_rb_pip = direct$diagnostics$rb_pip,
  collapsed_rb_pip = collapsed$diagnostics$rb_pip,
  rb_pip_diff = collapsed$diagnostics$rb_pip - direct$diagnostics$rb_pip,
  direct_gamma_ess = vapply(seq_len(10), function(j) {
    ess(direct$diagnostics$gamma_draws[j, ])
  }, numeric(1)),
  collapsed_gamma_ess = vapply(seq_len(10), function(j) {
    ess(collapsed$diagnostics$gamma_draws[j, ])
  }, numeric(1))
)
pip_summary$direct_gamma_ess_per_sec <- pip_summary$direct_gamma_ess /
  direct$runtime
pip_summary$collapsed_gamma_ess_per_sec <- pip_summary$collapsed_gamma_ess /
  collapsed$runtime

point_surface_summary <- rbind(
  surface_summary(direct, phi_grid, surface_predictors, locations, "direct"),
  surface_summary(collapsed, phi_grid, surface_predictors, locations,
                  "collapsed")
)

surface_agg <- rbind(
  surface_aggregate(direct, phi_grid, dat$grid$true_beta, dat$grid$true_u,
                    surface_predictors, "direct"),
  surface_aggregate(collapsed, phi_grid, dat$grid$true_beta, dat$grid$true_u,
                    surface_predictors, "collapsed")
)

performance <- rbind(
  data.frame(
    sampler = "direct",
    mspe = mean((dat$test$y - direct$pred_test)^2),
    beta_mise = mean((direct$beta_hat_grid - dat$grid$true_beta)^2),
    u_mise_x3_x6 = mean((direct$u_hat_grid[, 3:6] -
                           dat$grid$true_u[, 3:6])^2),
    theta_mse = mean((direct$theta_mean - dat$truth$theta)^2),
    runtime_sec = direct$runtime
  ),
  data.frame(
    sampler = "collapsed",
    mspe = mean((dat$test$y - collapsed$pred_test)^2),
    beta_mise = mean((collapsed$beta_hat_grid - dat$grid$true_beta)^2),
    u_mise_x3_x6 = mean((collapsed$u_hat_grid[, 3:6] -
                           dat$grid$true_u[, 3:6])^2),
    theta_mse = mean((collapsed$theta_mean - dat$truth$theta)^2),
    runtime_sec = collapsed$runtime
  )
)

conditional_audit <- {
  prep <- tmp_prepare
  p <- prep$p
  h <- prep$h
  omega <- rep(1, p)
  tau <- rep((h + 1) / 20^2, p)
  sigma2 <- var(dat$train$y)
  prior_non <- c(1 / pmax(omega, 1e-10), rep(1 / pmax(tau, 1e-10), each = h))
  Wfull <- prep$W
  y <- dat$train$y
  prec_full <- crossprod(Wfull) + diag(c(0, prior_non), ncol(Wfull))
  rhs_full <- crossprod(Wfull, y)
  mean_full <- solve(prec_full, rhs_full)
  cov_full <- solve(prec_full) * sigma2
  W <- Wfull[, -1, drop = FALSE]
  yc <- y - mean(y)
  Wc <- sweep(W, 2, colMeans(W), "-")
  prec_coll <- crossprod(Wc) + diag(prior_non, ncol(Wc))
  rhs_coll <- crossprod(Wc, yc)
  mean_coll <- solve(prec_coll, rhs_coll)
  cov_coll <- solve(prec_coll) * sigma2
  data.frame(
    compared_quantity = c(
      "marginal p(delta | others) mean",
      "marginal p(delta | others) covariance",
      "conditional beta0 mean"
    ),
    max_abs_difference = c(
      max(abs(mean_full[-1] - mean_coll)),
      max(abs(cov_full[-1, -1] - cov_coll)),
      abs(as.numeric(mean_full[1] - mean(y - W %*% mean_coll)))
    ),
    note = c(
      "This is direct joint conditional with beta0 integrated out, not p(delta | beta0, others).",
      "Schur-complement marginal covariance equals centered-design covariance.",
      "Collapsed sampler recovers beta0 from mean residual conditional on delta."
    )
  )
}

basis_audit <- {
  raw_meta <- fit_basis_metadata(dat$train$coords, basis_config$n_basis,
                                 centered = TRUE)
  phi_raw <- apply_basis_metadata(dat$train$coords, raw_meta)
  P_raw <- phi_raw %*% MASS::ginv(phi_raw)
  P_svd <- tmp_prepare$phi %*% MASS::ginv(tmp_prepare$phi)
  cp <- crossprod(tmp_prepare$phi)
  offdiag <- cp
  diag(offdiag) <- 0
  data.frame(
    K = basis_config$n_basis,
    raw_H = tmp_prepare$meta$raw_h,
    effective_h = tmp_prepare$h,
    centered_raw_rank = qr(phi_raw)$rank,
    svd_basis_rank = qr(tmp_prepare$phi)$rank,
    svd_full_rank = qr(tmp_prepare$phi)$rank == ncol(tmp_prepare$phi),
    max_training_column_mean = max(abs(colMeans(tmp_prepare$phi))),
    max_crossprod_offdiag = max(abs(offdiag)),
    crossprod_diag_min = min(diag(cp)),
    crossprod_diag_max = max(diag(cp)),
    column_space_projector_max_diff = max(abs(P_raw - P_svd))
  )
}

write.csv(basis_audit, file.path(out_root, "data", "svd_basis_audit.csv"),
          row.names = FALSE)
write.csv(conditional_audit,
          file.path(out_root, "data", "direct_collapsed_conditional_audit.csv"),
          row.names = FALSE)
write.csv(performance, file.path(out_root, "data", "sampler_performance.csv"),
          row.names = FALSE)
write.csv(theta_summary, file.path(out_root, "data", "theta_summary.csv"),
          row.names = FALSE)
write.csv(beta0_summary, file.path(out_root, "data", "beta0_summary.csv"),
          row.names = FALSE)
write.csv(sigma_pip_summary,
          file.path(out_root, "data", "sigma2_summary.csv"),
          row.names = FALSE)
write.csv(pip_summary, file.path(out_root, "data", "pip_summary.csv"),
          row.names = FALSE)
write.csv(point_surface_summary,
          file.path(out_root, "data", "point_surface_posterior_summary.csv"),
          row.names = FALSE)
write.csv(surface_agg,
          file.path(out_root, "data", "surface_aggregate_summary.csv"),
          row.names = FALSE)
saveRDS(list(data = dat, direct = direct, collapsed = collapsed,
             collapsed_beta0_recovered = collapsed_beta0_recovered,
             basis_audit = basis_audit,
             conditional_audit = conditional_audit,
             performance = performance,
             theta_summary = theta_summary,
             beta0_summary = beta0_summary,
             sigma2_summary = sigma_pip_summary,
             pip_summary = pip_summary,
             point_surface_summary = point_surface_summary,
             surface_aggregate = surface_agg),
        file.path(out_root, "data",
                  "direct_collapsed_mcmc_svd_comparison.rds"))

fmt <- function(x) paste(capture.output(print(x, row.names = FALSE)),
                         collapse = "\n")
report <- c(
  "# Direct vs Collapsed Intercept Sampler MCMC Comparison with SVD Basis",
  "",
  paste0("- Output root: `", out_root, "`"),
  "- Scope: one simulated p=10 data set, K=6, SVD full-rank centered basis.",
  "- No full simulation or real-data analysis was rerun.",
  "",
  "## Basis Audit",
  fmt(basis_audit),
  "",
  "## Conditional Audit Clarification",
  "The comparison is between the direct joint Gaussian conditional p(beta0, delta | others), marginalized over beta0, and the collapsed conditional p(delta | others). It is not comparing p(delta | beta0, others).",
  "",
  fmt(conditional_audit),
  "",
  "## Performance",
  fmt(performance),
  "",
  "## beta0",
  fmt(beta0_summary),
  "",
  "## sigma2",
  fmt(sigma_pip_summary),
  "",
  "## PIPs",
  fmt(pip_summary),
  "",
  "## theta",
  fmt(theta_summary),
  "",
  "## Surface Aggregates",
  fmt(surface_agg),
  "",
  "## Pointwise u/beta Posterior Summaries",
  fmt(point_surface_summary)
)
writeLines(report, file.path(out_root, "reports",
                             "direct_collapsed_mcmc_svd_report.md"))

cat("Output root:", out_root, "\n")
cat("\nBasis audit:\n")
print(basis_audit, row.names = FALSE)
cat("\nConditional audit:\n")
print(conditional_audit, row.names = FALSE)
cat("\nPerformance:\n")
print(performance, row.names = FALSE)
cat("\nBeta0:\n")
print(beta0_summary, row.names = FALSE)
cat("\nSigma2:\n")
print(sigma_pip_summary, row.names = FALSE)
cat("\nPIPs:\n")
print(pip_summary, row.names = FALSE)
cat("\nTheta selected predictors:\n")
print(theta_summary[theta_summary$predictor %in% paste0("X", selected_predictors), ],
      row.names = FALSE)
cat("\nSurface aggregates:\n")
print(surface_agg, row.names = FALSE)
cat("\nReport:", file.path(out_root, "reports",
                          "direct_collapsed_mcmc_svd_report.md"), "\n")
