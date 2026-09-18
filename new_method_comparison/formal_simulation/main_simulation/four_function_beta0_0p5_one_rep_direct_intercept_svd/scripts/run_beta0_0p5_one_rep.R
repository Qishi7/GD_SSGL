#!/usr/bin/env Rscript

options(digits = 17)
suppressPackageStartupMessages({
  library(Rcpp)
  library(GIGrvg)
  library(mgcv)
})

comparison_root <- normalizePath("new_method_comparison", winslash = "/",
                                 mustWork = TRUE)
formal_root <- file.path(comparison_root, "formal_simulation")
out_root <- file.path(formal_root, "main_simulation",
                      "four_function_beta0_0p5_one_rep_direct_intercept_svd")
dirs <- c("config", "data", "fits", "results", "reports", "logs", "scripts")
invisible(lapply(file.path(out_root, dirs), dir.create, recursive = TRUE,
                 showWarnings = FALSE))

Sys.setenv(SSGL_SOURCE_FUNCTIONS_ONLY = "1")
source(file.path(formal_root, "scripts",
                 "run_four_function_10rep_uncertainty_calibration_Kmax6.R"))
source(file.path(comparison_root, "R", "methods", "accelerated",
                 "newssgl_intercept.R"))
source(file.path(comparison_root, "R", "methods", "wrappers.R"))
load_newssgl_intercept_fast()
Rcpp::sourceCpp(file.path(comparison_root, "R", "methods", "accelerated",
                          "original_ssgl_intercept.cpp"))
Rcpp::sourceCpp(file.path(comparison_root, "R", "methods", "accelerated",
                          "full_svc_intercept.cpp"))

true_beta0 <- 0.5
replicate <- 1L
p <- 10L
K <- 6L
theta_true <- c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0)
predictor_names <- paste0("X", seq_len(p))
active_deviation <- c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE,
                      FALSE, FALSE, FALSE, FALSE)
active_surface <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                    FALSE, FALSE, FALSE, FALSE)
seeds <- data.frame(
  replicate = 1L,
  data_seed = 2026100002L,
  proposed_seed = 2026410001L,
  original_seed_base = 2026420100L,
  full_svc_seed_base = 2026430100L,
  blasso_seed_base = 2026440100L
)

write.csv(data.frame(
  item = c("true_beta0", "replicate", "data_seed", "K_B_raw", "H_raw",
           "effective_centered_H", "intercept_rule"),
  value = c(true_beta0, replicate, seeds$data_seed, K, K^2L, K^2L - 1L,
            "direct scalar beta0 with p(beta0) proportional to 1")
), file.path(out_root, "config", "manifest.csv"), row.names = FALSE)

apply_phi <- function(coords, meta) {
  if (isTRUE(meta$full_rank_centered)) {
    apply_centered_basis_metadata_fullrank(coords, meta)
  } else {
    apply_basis_metadata(coords, meta)
  }
}

surface_theta_alpha <- function(theta, alpha, phi, p) {
  h <- ncol(phi)
  out <- matrix(0, nrow(phi), p)
  for (j in seq_len(p)) {
    idx <- ((j - 1L) * h + 1L):(j * h)
    out[, j] <- theta[j] + as.vector(phi %*% alpha[idx])
  }
  out
}

surface_eta <- function(eta, phi, p) {
  h <- ncol(phi)
  out <- matrix(0, nrow(phi), p)
  for (j in seq_len(p)) {
    idx <- ((j - 1L) * h + 1L):(j * h)
    out[, j] <- as.vector(phi %*% eta[idx])
  }
  out
}

metric_row <- function(method, pred_test, beta_grid, theta_hat, u_grid, data,
                       runtime_sec, beta0_hat, sampler, basis_note) {
  beta_error <- colMeans((beta_grid - data$grid$true_beta)^2)
  data.frame(
    replicate = replicate,
    method = method,
    true_beta0 = true_beta0,
    beta0_hat = beta0_hat,
    beta0_error = beta0_hat - true_beta0,
    mspe = mean((data$test$y - pred_test)^2),
    rmse = sqrt(mean((data$test$y - pred_test)^2)),
    mae = mean(abs(data$test$y - pred_test)),
    bias = mean(data$test$y - pred_test),
    beta_mise = mean(beta_error),
    theta_mse = if (is.null(theta_hat)) NA_real_ else
      mean((theta_hat - theta_true)^2),
    u_mise_x3_x6 = if (is.null(u_grid)) NA_real_ else
      mean(colMeans((u_grid[, 3:6, drop = FALSE] -
                       data$grid$true_u[, 3:6, drop = FALSE])^2)),
    runtime_sec = runtime_sec,
    sampler = sampler,
    basis_note = basis_note,
    stringsAsFactors = FALSE
  )
}

fit_original_three_chain <- function(data) {
  started <- proc.time()[3]
  chains <- lapply(1:3, function(ch) {
    fit_original_ssgl_intercept(
      data$train, data$test, data$grid,
      list(n_basis = K),
      list(lambda0 = 15, lambda1 = 3, a_gamma = 1, b_gamma = 10,
           a_sigma = 0.5, b_sigma = var(data$train$y) / 2,
           zeta0 = 0.1, zeta1 = 50),
      list(n_iter = 5000L, burn_in = 2000L),
      seeds$original_seed_base + ch
    )
  })
  meta <- chains[[1]]$meta
  eta <- do.call(cbind, lapply(chains, `[[`, "eta"))
  gamma <- do.call(cbind, lapply(chains, `[[`, "gamma"))
  gamma_prob <- do.call(cbind, lapply(chains, `[[`, "gamma_prob"))
  beta0 <- unlist(lapply(chains, `[[`, "beta0_draws"))
  eta_mean <- rowMeans(eta)
  beta0_mean <- mean(beta0)
  phi_test <- apply_basis_metadata(data$test$coords, meta)
  phi_grid <- apply_basis_metadata(data$grid$coords, meta)
  beta_test <- surface_eta(eta_mean, phi_test, p)
  beta_grid <- surface_eta(eta_mean, phi_grid, p)
  list(meta = meta, eta = eta, gamma = gamma, gamma_prob = gamma_prob,
       beta0_draws = beta0, beta0_mean = beta0_mean,
       beta_grid = beta_grid, pred_test = beta0_mean + rowSums(data$test$X * beta_test),
       runtime = proc.time()[3] - started)
}

fit_full_svc_three_chain <- function(data) {
  started <- proc.time()[3]
  chains <- lapply(1:3, function(ch) {
    fit_full_svc_no_selection_intercept(
      data$train, data$test, data$grid,
      list(n_basis = K, full_rank_centered = TRUE, full_rank_method = "svd"),
      list(kappa2_alpha = 8, a_sigma = 0.5,
           b_sigma = var(data$train$y) / 2, a_theta = 1, b_theta = 1),
      list(n_iter = 5000L, burn_in = 2000L),
      seeds$full_svc_seed_base + ch
    )
  })
  meta <- chains[[1]]$config$basis
  theta <- do.call(cbind, lapply(chains, `[[`, "theta_draws"))
  alpha <- do.call(cbind, lapply(chains, function(x) x$diagnostics$alpha_draws))
  beta0 <- unlist(lapply(chains, `[[`, "beta0_draws"))
  theta_mean <- rowMeans(theta)
  alpha_mean <- rowMeans(alpha)
  beta0_mean <- mean(beta0)
  phi_test <- apply_phi(data$test$coords, meta)
  phi_grid <- apply_phi(data$grid$coords, meta)
  beta_test <- surface_theta_alpha(theta_mean, alpha_mean, phi_test, p)
  beta_grid <- surface_theta_alpha(theta_mean, alpha_mean, phi_grid, p)
  list(meta = meta, theta = theta, alpha = alpha, beta0_draws = beta0,
       beta0_mean = beta0_mean, theta_mean = theta_mean,
       beta_grid = beta_grid, u_grid = sweep(beta_grid, 2, theta_mean, "-"),
       pred_test = beta0_mean + rowSums(data$test$X * beta_test),
       runtime = proc.time()[3] - started)
}

fit_blasso_three_chain <- function(data) {
  started <- proc.time()[3]
  chains <- lapply(1:3, function(ch) {
    fit_global_only_blasso_intercept(
      data$train, data$test, data$grid,
      list(),
      list(a_sigma = 0.5, b_sigma = var(data$train$y) / 2,
           a_theta = 1, b_theta = 1),
      list(n_iter = 5000L, burn_in = 2000L),
      seeds$blasso_seed_base + ch
    )
  })
  theta <- do.call(cbind, lapply(chains, `[[`, "theta_draws"))
  beta0 <- unlist(lapply(chains, `[[`, "beta0_draws"))
  theta_mean <- rowMeans(theta)
  beta0_mean <- mean(beta0)
  beta_grid <- matrix(rep(theta_mean, each = nrow(data$grid$coords)),
                      nrow(data$grid$coords), p)
  list(theta = theta, beta0_draws = beta0, beta0_mean = beta0_mean,
       theta_mean = theta_mean, beta_grid = beta_grid,
       u_grid = matrix(0, nrow(data$grid$coords), p),
       pred_test = beta0_mean + as.vector(data$test$X %*% theta_mean),
       runtime = proc.time()[3] - started)
}

fit_gam_intercept <- function(data) {
  started <- proc.time()[3]
  train_df <- as.data.frame(data$train$X)
  names(train_df) <- predictor_names
  train_df$y <- data$train$y
  train_df$s1 <- data$train$coords[, 1]
  train_df$s2 <- data$train$coords[, 2]
  form <- as.formula(paste0(
    "y ~ ", paste(predictor_names, collapse = " + "), " + ",
    paste0("s(s1, s2, bs = 'tp', by = ", predictor_names, ")",
           collapse = " + ")
  ))
  fit <- mgcv::gam(form, data = train_df, method = "REML")
  test_df <- as.data.frame(data$test$X)
  names(test_df) <- predictor_names
  test_df$s1 <- data$test$coords[, 1]
  test_df$s2 <- data$test$coords[, 2]
  pred_test <- as.numeric(predict(fit, newdata = test_df))
  base_grid <- data.frame(s1 = data$grid$coords[, 1],
                          s2 = data$grid$coords[, 2])
  for (jj in seq_len(p)) base_grid[[predictor_names[jj]]] <- 0
  baseline_grid <- as.numeric(predict(fit, newdata = base_grid))
  beta_grid <- matrix(0, nrow(data$grid$coords), p)
  for (j in seq_len(p)) {
    gd <- base_grid
    gd[[predictor_names[j]]] <- 1
    beta_grid[, j] <- as.numeric(predict(fit, newdata = gd)) - baseline_grid
  }
  theta_hat <- colMeans(beta_grid)
  list(fit = fit, beta0_mean = unname(coef(fit)["(Intercept)"]),
       theta_mean = theta_hat, beta_grid = beta_grid,
       u_grid = sweep(beta_grid, 2, theta_hat, "-"),
       pred_test = pred_test, runtime = proc.time()[3] - started)
}

log_file <- file.path(out_root, "logs", "beta0_0p5_one_rep.log")
sink(log_file, append = TRUE, split = TRUE)
on.exit(sink(), add = TRUE)

data <- generate_four_function_data(seeds$data_seed)
data$train$y <- data$train$y + true_beta0
data$test$y <- data$test$y + true_beta0
data$truth$beta0 <- true_beta0
saveRDS(data, file.path(out_root, "data", "beta0_0p5_rep001_data.rds"))

cat("GD-SSGL...\n")
gd <- fit_newssgl_intercept_fast(
  data$train, data$test, data$grid,
  list(n_basis = K, full_rank_centered = TRUE, full_rank_method = "svd"),
  list(lambda0 = 20, lambda1 = 2, a_sigma = 0.5,
       b_sigma = var(data$train$y) / 2,
       a_theta = 1, b_theta = 1, a_gamma = 1, b_gamma = 10),
  list(n_iter = 10000L, burn_in = 3000L),
  seeds$proposed_seed
)
saveRDS(gd, file.path(out_root, "fits", "gdssgl_fit.rds"))

cat("WS-SSGL...\n")
ws <- fit_original_three_chain(data)
saveRDS(ws, file.path(out_root, "fits", "wsssgl_fit.rds"))

cat("Gaussian SVC...\n")
svc <- fit_full_svc_three_chain(data)
saveRDS(svc, file.path(out_root, "fits", "gaussian_svc_fit.rds"))

cat("Bayesian Lasso...\n")
bl <- fit_blasso_three_chain(data)
saveRDS(bl, file.path(out_root, "fits", "bayesian_lasso_fit.rds"))

cat("GAM...\n")
gam <- fit_gam_intercept(data)
saveRDS(gam, file.path(out_root, "fits", "gam_intercept_fit.rds"))

metrics <- rbind(
  metric_row("GD-SSGL", gd$pred_test, gd$beta_hat_grid, gd$theta_mean,
             gd$u_hat_grid, data, gd$runtime, gd$beta0_mean,
             "direct beta0", "SVD full-rank centered H=35"),
  metric_row("WS-SSGL", ws$pred_test, ws$beta_grid, NULL, NULL, data,
             ws$runtime, ws$beta0_mean, "direct beta0",
             "original uncentered tensor basis H=36"),
  metric_row("Gaussian SVC", svc$pred_test, svc$beta_grid, svc$theta_mean,
             svc$u_grid, data, svc$runtime, svc$beta0_mean,
             "direct beta0", "SVD full-rank centered H=35"),
  metric_row("Bayesian Lasso", bl$pred_test, bl$beta_grid, bl$theta_mean,
             bl$u_grid, data, bl$runtime, bl$beta0_mean,
             "direct beta0", "no spatial basis"),
  metric_row("Standard GAM + intercept", gam$pred_test, gam$beta_grid,
             gam$theta_mean, gam$u_grid, data, gam$runtime, gam$beta0_mean,
             "mgcv REML with formula intercept", "mgcv thin-plate default")
)

pips <- rbind(
  data.frame(method = "GD-SSGL", predictor = predictor_names,
             target = active_deviation, sampled_pip = gd$pip,
             rb_pip = gd$diagnostics$rb_pip,
             selection_target = "spatial deviation"),
  data.frame(method = "WS-SSGL", predictor = predictor_names,
             target = active_surface, sampled_pip = rowMeans(ws$gamma),
             rb_pip = rowMeans(ws$gamma_prob),
             selection_target = "whole coefficient surface")
)

write.csv(metrics, file.path(out_root, "results",
                             "beta0_0p5_one_rep_metrics.csv"),
          row.names = FALSE)
write.csv(pips, file.path(out_root, "results",
                          "beta0_0p5_one_rep_pips.csv"),
          row.names = FALSE)

sink(file.path(out_root, "reports", "beta0_0p5_one_rep_report.md"))
cat("# beta0 = 0.5 one-replicate direct-intercept check\n\n")
cat("True beta0:", true_beta0, "\n\n")
print(metrics)
cat("\n## PIPs\n\n")
print(pips)
sink()

cat("Saved to: ", out_root, "\n", sep = "")
print(metrics)
