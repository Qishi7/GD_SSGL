#!/usr/bin/env Rscript

options(digits = 17, warn = 1)

suppressPackageStartupMessages({
  library(Rcpp)
})

project_root <- normalizePath(Sys.getenv("GDSSGL_ROOT", getwd()), winslash = "/", mustWork = TRUE)
comparison_root <- file.path(project_root, "new_method_comparison")
formal_root <- file.path(comparison_root, "formal_simulation")
out_root <- file.path(
  formal_root, "main_simulation",
  "w6_onbasis_K6_beta0_0p5_intercept_svd_one_rep_pilot"
)
invisible(lapply(file.path(out_root, c("config", "metadata", "truth", "data",
                                       "results", "reports", "logs",
                                       "replicate_fits/replicate_001")),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

source(file.path(comparison_root, "R", "methods", "accelerated",
                 "newssgl_intercept.R"))
load_newssgl_intercept_fast()

p <- 10L
K_B <- 6L
true_beta0 <- 0.5
n_train <- 700L
n_test <- 250L
sigma_epsilon <- 0.35
grid_size_eval <- 50L
grid_size_reference <- 101L
theta_true <- c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0)
active_u <- c(3L, 4L, 5L, 6L)
predictor_names <- paste0("X", seq_len(p))
lambda0 <- 20
lambda1 <- 2
n_iter <- 10000L
burn_in <- 3000L
data_seed <- 2026100002L
mcmc_seed <- 2026410001L

writeLines(c(
  "scenario: W6_onbasis_K6_beta0_0p5_intercept_svd_one_rep_pilot",
  "replicates: 1",
  "method: GD-SSGL only",
  sprintf("true_beta0: %.1f", true_beta0),
  sprintf("n_train: %d", n_train),
  sprintf("n_test: %d", n_test),
  sprintf("p: %d", p),
  sprintf("sigma_epsilon: %.2f", sigma_epsilon),
  sprintf("K_B: %d", K_B),
  "raw_H: 36",
  "effective_H: 35",
  sprintf("lambda0: %.0f", lambda0),
  sprintf("lambda1: %.0f", lambda1),
  sprintf("n_iter: %d", n_iter),
  sprintf("burn_in: %d", burn_in),
  "thin: 1",
  "sampler: fit_newssgl_intercept_fast direct scalar beta0",
  "basis: SVD full-rank centered basis, Phi_full = Phi_c %*% V_r",
  paste0("theta_true: [", paste(theta_true, collapse = ", "), "]"),
  "active_spatial_deviation_set: [3, 4, 5, 6]",
  sprintf("data_seed: %d", data_seed),
  sprintf("mcmc_seed: %d", mcmc_seed)
), file.path(out_root, "config", "w6_beta0_0p5_manifest.yml"))

g1 <- function(s1, s2) sin(2 * pi * s1) * cos(2 * pi * s2)
g2 <- function(s1, s2) cos(2 * pi * s1) * sin(2 * pi * s2)
g3 <- function(s1, s2) exp(-((s1 - 0.35)^2 + (s2 - 0.65)^2) / 0.04)
g4 <- function(s1, s2) {
  exp(-((s1 - 0.70)^2 + (s2 - 0.30)^2) / 0.03) -
    exp(-((s1 - 0.25)^2 + (s2 - 0.25)^2) / 0.03)
}
analytic_funs <- list(`3` = g1, `4` = g2, `5` = g3, `6` = g4)

weighted_rms <- function(x, w) sqrt(sum(w * x^2))
weighted_cor <- function(x, y, w) {
  sum(w * x * y) / (weighted_rms(x, w) * weighted_rms(y, w))
}
pinv_solve <- function(A, y) {
  s <- svd(A)
  tol <- max(dim(A)) * max(s$d) * .Machine$double.eps
  rank <- sum(s$d > tol)
  coef <- s$v[, seq_len(rank), drop = FALSE] %*%
    ((t(s$u[, seq_len(rank), drop = FALSE]) %*% y) / s$d[seq_len(rank)])
  list(coef = as.numeric(coef), rank = rank, tol = tol)
}
surface_from_alpha <- function(theta, alpha, phi, p) {
  h <- ncol(phi)
  out <- matrix(0, nrow(phi), p)
  for (j in seq_len(p)) {
    idx <- ((j - 1L) * h + 1L):(j * h)
    out[, j] <- theta[j] + as.vector(phi %*% alpha[idx])
  }
  out
}
band_summary_for <- function(draw_mat, truth_vec) {
  mean_s <- rowMeans(draw_mat)
  sd_s <- apply(draw_mat, 1, sd)
  sd_s[sd_s < 1e-12] <- 1e-12
  max_std <- apply(abs(sweep(draw_mat, 1, mean_s, "-") / sd_s), 2, max)
  crit <- quantile(max_std, 0.95, names = FALSE)
  lower <- mean_s - crit * sd_s
  upper <- mean_s + crit * sd_s
  data.frame(
    simultaneous_coverage = all(truth_vec >= lower & truth_vec <= upper),
    simultaneous_mean_width = mean(upper - lower),
    simultaneous_crit = crit
  )
}

set.seed(data_seed)
n <- n_train + n_test
coords <- cbind(runif(n), runif(n))
X_raw <- matrix(rnorm(n * p), n, p)
epsilon <- rnorm(n, sd = sigma_epsilon)
train_id <- seq_len(n_train)
test_id <- n_train + seq_len(n_test)
x_mean <- colMeans(X_raw[train_id, , drop = FALSE])
x_sd <- apply(X_raw[train_id, , drop = FALSE], 2, sd)
X <- sweep(sweep(X_raw, 2, x_mean, "-"), 2, x_sd, "/")
train_coords <- coords[train_id, , drop = FALSE]
test_coords <- coords[test_id, , drop = FALSE]
X_train <- X[train_id, , drop = FALSE]
X_test <- X[test_id, , drop = FALSE]

eval_grid <- as.matrix(expand.grid(
  s1 = seq(0, 1, length.out = grid_size_eval),
  s2 = seq(0, 1, length.out = grid_size_eval)
))
ref_grid <- as.matrix(expand.grid(
  s1 = seq(0, 1, length.out = grid_size_reference),
  s2 = seq(0, 1, length.out = grid_size_reference)
))
ref_w <- rep(1 / nrow(ref_grid), nrow(ref_grid))

basis_meta <- fit_centered_basis_metadata_fullrank(
  train_coords, K_B, method = "svd"
)
B_ref <- apply_centered_basis_metadata_fullrank(ref_grid, basis_meta)
B_train <- apply_centered_basis_metadata_fullrank(train_coords, basis_meta)
B_test <- apply_centered_basis_metadata_fullrank(test_coords, basis_meta)
B_eval <- apply_centered_basis_metadata_fullrank(eval_grid, basis_meta)
h_eff <- ncol(B_train)

alpha_star <- matrix(0, h_eff, p)
truth_rows <- vector("list", length(active_u))
for (ii in seq_along(active_u)) {
  j <- active_u[ii]
  g <- analytic_funs[[as.character(j)]](ref_grid[, 1], ref_grid[, 2])
  g_c <- g - sum(ref_w * g)
  fit0 <- pinv_solve(sqrt(ref_w) * B_ref, sqrt(ref_w) * g_c)
  h0 <- as.vector(B_ref %*% fit0$coef)
  scale_j <- weighted_rms(g_c, ref_w) / weighted_rms(h0, ref_w)
  a_star <- scale_j * fit0$coef
  h_star <- as.vector(B_ref %*% a_star)
  alpha_star[, j] <- a_star
  truth_rows[[ii]] <- data.frame(
    predictor = paste0("X", j),
    original_RMS = weighted_rms(g_c, ref_w),
    W6_RMS = weighted_rms(h_star, ref_w),
    RMS_ratio = weighted_rms(h_star, ref_w) / weighted_rms(g_c, ref_w),
    pre_rescaling_projection_RMSE = weighted_rms(g_c - h0, ref_w),
    post_rescaling_RMSE = weighted_rms(g_c - h_star, ref_w),
    shape_correlation = weighted_cor(g_c, h_star, ref_w),
    effective_basis_rank = h_eff,
    raw_H = basis_meta$raw_h,
    centered_rank = basis_meta$centered_rank,
    stringsAsFactors = FALSE
  )
}
truth_validation <- do.call(rbind, truth_rows)
write.csv(truth_validation, file.path(out_root, "truth",
                                      "w6_beta0_0p5_truth_validation_one_rep.csv"),
          row.names = FALSE)

u_train <- B_train %*% alpha_star
u_test <- B_test %*% alpha_star
u_eval <- B_eval %*% alpha_star
beta_train <- sweep(u_train, 2, theta_true, "+")
beta_test <- sweep(u_test, 2, theta_true, "+")
beta_eval <- sweep(u_eval, 2, theta_true, "+")
y_train <- true_beta0 + rowSums(X_train * beta_train) + epsilon[train_id]
y_test <- true_beta0 + rowSums(X_test * beta_test) + epsilon[test_id]

data_obj <- list(
  train = list(y = y_train, X = X_train, coords = train_coords),
  test = list(y = y_test, X = X_test, coords = test_coords),
  grid = list(coords = eval_grid, true_beta = beta_eval,
              true_u = u_eval, grid_size = grid_size_eval),
  truth = list(theta = theta_true, beta0 = true_beta0,
               spatial_deviation = seq_len(p) %in% active_u,
               active_surface = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                                  FALSE, FALSE, FALSE, FALSE),
               alpha_star_svd = alpha_star,
               basis = basis_meta),
  preprocessing = list(x_mean = x_mean, x_sd = x_sd),
  config = list(n_train = n_train, n_test = n_test, p = p,
                sigma_epsilon = sigma_epsilon, K_B = K_B,
                raw_H = basis_meta$raw_h, effective_H = h_eff)
)
saveRDS(data_obj, file.path(out_root, "data", "w6_beta0_0p5_rep_001.rds"))

response_identity_error <- max(abs(y_train - (
  true_beta0 + rowSums(X_train * beta_train) + epsilon[train_id]
)))
validation_checks <- data.frame(
  check = c("true_beta0_0p5", "effective_H_35", "training_basis_centering",
            "training_u_centering", "training_beta_mean_equals_theta",
            "active_set_identity", "response_identity", "finite_values"),
  value = c(
    true_beta0,
    h_eff,
    max(abs(colMeans(B_train))),
    max(abs(colMeans(u_train[, active_u, drop = FALSE]))),
    max(abs(colMeans(beta_train[, active_u, drop = FALSE]) - theta_true[active_u])),
    identical(which(colSums(abs(u_eval)) > 1e-12), active_u),
    response_identity_error,
    all(is.finite(unlist(list(y_train, y_test, X_train, X_test, beta_eval, u_eval))))
  ),
  pass = c(
    identical(true_beta0, 0.5),
    h_eff == 35L,
    max(abs(colMeans(B_train))) < 1e-10,
    max(abs(colMeans(u_train[, active_u, drop = FALSE]))) < 1e-10,
    max(abs(colMeans(beta_train[, active_u, drop = FALSE]) - theta_true[active_u])) < 1e-10,
    identical(which(colSums(abs(u_eval)) > 1e-12), active_u),
    response_identity_error < 1e-12,
    all(is.finite(unlist(list(y_train, y_test, X_train, X_test, beta_eval, u_eval))))
  ),
  stringsAsFactors = FALSE
)
write.csv(validation_checks, file.path(out_root, "results",
                                       "w6_beta0_0p5_validation_checks.csv"),
          row.names = FALSE)
if (!all(validation_checks$pass)) stop("Validation failed; not fitting.")

started <- proc.time()[3]
fit <- fit_newssgl_intercept_fast(
  data_obj$train, data_obj$test, data_obj$grid,
  list(n_basis = K_B, full_rank_centered = TRUE, full_rank_method = "svd"),
  list(lambda0 = lambda0, lambda1 = lambda1, a_sigma = 0.5,
       b_sigma = var(y_train) / 2, a_theta = 1, b_theta = 1,
       a_gamma = 1, b_gamma = 10),
  list(n_iter = n_iter, burn_in = burn_in),
  mcmc_seed
)
runtime_sec <- proc.time()[3] - started
saveRDS(fit, file.path(out_root, "replicate_fits", "replicate_001",
                       "gdssgl_intercept_svd_fit.rds"))

phi_grid <- apply_centered_basis_metadata_fullrank(eval_grid, fit$config$basis)
alpha_draws <- fit$diagnostics$alpha_draws
n_draw <- ncol(fit$theta_draws)
coverage_rows <- list()
beta_mean_all <- matrix(0, nrow(eval_grid), p)
u_mean_all <- matrix(0, nrow(eval_grid), p)
for (j in seq_len(p)) {
  idx <- ((j - 1L) * h_eff + 1L):(j * h_eff)
  u_draws <- phi_grid %*% alpha_draws[idx, , drop = FALSE]
  beta_draws <- sweep(u_draws, 2, fit$theta_draws[j, ], "+")
  beta_mean_all[, j] <- rowMeans(beta_draws)
  u_mean_all[, j] <- rowMeans(u_draws)
  if (j %in% active_u) {
    uq <- rbind(lower = apply(u_draws, 1, quantile, 0.025, names = FALSE),
                upper = apply(u_draws, 1, quantile, 0.975, names = FALSE))
    bq <- rbind(lower = apply(beta_draws, 1, quantile, 0.025, names = FALSE),
                upper = apply(beta_draws, 1, quantile, 0.975, names = FALSE))
    ub <- band_summary_for(u_draws, u_eval[, j])
    bb <- band_summary_for(beta_draws, beta_eval[, j])
    coverage_rows[[length(coverage_rows) + 1L]] <- data.frame(
      predictor = predictor_names[j],
      u_pointwise_coverage = mean(u_eval[, j] >= uq["lower", ] & u_eval[, j] <= uq["upper", ]),
      beta_pointwise_coverage = mean(beta_eval[, j] >= bq["lower", ] & beta_eval[, j] <= bq["upper", ]),
      u_mean_width = mean(uq["upper", ] - uq["lower", ]),
      beta_mean_width = mean(bq["upper", ] - bq["lower", ]),
      u_simultaneous_coverage = ub$simultaneous_coverage,
      beta_simultaneous_coverage = bb$simultaneous_coverage,
      u_simultaneous_mean_width = ub$simultaneous_mean_width,
      beta_simultaneous_mean_width = bb$simultaneous_mean_width,
      stringsAsFactors = FALSE
    )
  }
}
coverage <- do.call(rbind, coverage_rows)
beta_mise_j <- colMeans((beta_mean_all - beta_eval)^2)
u_mise_j <- colMeans((u_mean_all - u_eval)^2)
rb_pip <- fit$diagnostics$rb_pip
selected <- rb_pip >= 0.5
active <- seq_len(p) %in% active_u
pip_table <- data.frame(predictor = predictor_names, target = as.integer(active),
                        rb_pip = rb_pip, sampled_pip = fit$pip,
                        selected = selected)
rep_level <- data.frame(
  replicate = 1L,
  method = "GD-SSGL intercept SVD W6",
  true_beta0 = true_beta0,
  beta0_hat = fit$beta0_mean,
  mspe = mean((y_test - fit$pred_test)^2),
  beta_mise = mean(beta_mise_j),
  u_mise_x3_x6 = mean(u_mise_j[active_u]),
  theta_mse = mean((fit$theta_mean - theta_true)^2),
  runtime_sec = runtime_sec,
  selected_predictors = paste(predictor_names[selected], collapse = ", "),
  TPR = sum(selected & active) / sum(active),
  FPR = sum(selected & !active) / sum(!active),
  FDR = if (sum(selected) > 0) sum(selected & !active) / sum(selected) else 0,
  exact_recovery = identical(selected, active),
  selected_size = sum(selected)
)
write.csv(coverage, file.path(out_root, "results",
                              "w6_beta0_0p5_coverage_one_rep.csv"),
          row.names = FALSE)
write.csv(pip_table, file.path(out_root, "results",
                               "w6_beta0_0p5_pips_one_rep.csv"),
          row.names = FALSE)
write.csv(rep_level, file.path(out_root, "results",
                               "w6_beta0_0p5_replicate_level_one_rep.csv"),
          row.names = FALSE)

report <- c(
  "# W6 beta0=0.5 intercept/SVD one-replicate coverage pilot",
  "",
  paste0("- Output root: `", out_root, "`"),
  "- R = 1 only; GD-SSGL only.",
  "- This tests the previous W6 on-basis coverage method under the current direct-intercept SVD implementation.",
  "",
  "## Truth validation",
  "",
  paste(capture.output(print(truth_validation, row.names = FALSE)), collapse = "\n"),
  "",
  "## Checks",
  "",
  paste(capture.output(print(validation_checks, row.names = FALSE)), collapse = "\n"),
  "",
  "## Performance",
  "",
  paste(capture.output(print(rep_level, row.names = FALSE)), collapse = "\n"),
  "",
  "## PIPs",
  "",
  paste(capture.output(print(pip_table, row.names = FALSE)), collapse = "\n"),
  "",
  "## Coverage",
  "",
  paste(capture.output(print(coverage, row.names = FALSE)), collapse = "\n")
)
writeLines(report, file.path(out_root, "reports",
                             "w6_beta0_0p5_intercept_svd_one_rep_report.md"))
cat(paste(report, collapse = "\n"))
