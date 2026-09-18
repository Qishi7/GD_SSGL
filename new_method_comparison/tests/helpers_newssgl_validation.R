max_abs_difference <- function(x, y) {
  max(abs(as.numeric(x) - as.numeric(y)))
}

relative_difference <- function(x, y) {
  abs(x - y) / max(abs(x), abs(y), 1e-12)
}

surface_rmse <- function(x, y) {
  sqrt(mean((x - y)^2))
}

posterior_summary_newssgl <- function(fit) {
  list(
    theta = fit$theta_mean,
    alpha = fit$alpha_mean,
    sigma2 = mean(fit$diagnostics$sigma2_draws),
    omega = rowMeans(fit$diagnostics$omega_draws),
    tau = rowMeans(fit$diagnostics$tau_draws),
    lambda_theta2 = mean(fit$diagnostics$lambda_theta2_draws),
    pi_gamma = mean(fit$diagnostics$pi_gamma_draws),
    sampled_pip = fit$pip,
    rb_pip = fit$diagnostics$rb_pip,
    beta_grid = fit$beta_hat_grid,
    pred_test = fit$pred_test
  )
}

aggregate_posterior_summaries <- function(summaries) {
  fields <- names(summaries[[1]])
  out <- list()
  for (field in fields) {
    values <- lapply(summaries, `[[`, field)
    out[[field]] <- Reduce("+", values) / length(values)
  }
  out
}

compute_newssgl_metrics <- function(fit, data) {
  beta_error <- colMeans((fit$beta_hat_grid - data$grid$true_beta)^2)
  u_error <- colMeans((fit$u_hat_grid - data$grid$true_u)^2)
  c(
    mspe = mean((data$test$y - fit$pred_test)^2),
    beta_mise = mean(beta_error),
    mise_global_only = mean(beta_error[1:2]),
    mise_spatial_only = mean(beta_error[3:4]),
    mise_global_plus_spatial = mean(beta_error[5:6]),
    mise_null = mean(beta_error[7:10]),
    theta_mse = mean((fit$theta_mean - data$truth$theta)^2),
    u_mise_x3_x6 = mean(u_error[3:6])
  )
}

make_test_result <- function(group, test, value, tolerance, passed,
                             detail = "") {
  data.frame(
    group = group,
    test = test,
    value = as.numeric(value),
    tolerance = as.numeric(tolerance),
    passed = isTRUE(passed),
    detail = detail,
    stringsAsFactors = FALSE
  )
}
