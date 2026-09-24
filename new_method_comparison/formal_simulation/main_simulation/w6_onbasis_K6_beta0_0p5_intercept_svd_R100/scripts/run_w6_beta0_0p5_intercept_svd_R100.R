#!/usr/bin/env Rscript

options(digits = 17, warn = 1)

suppressPackageStartupMessages({
  library(Rcpp)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

project_root <- normalizePath(Sys.getenv("GDSSGL_ROOT", getwd()), winslash = "/", mustWork = TRUE)
comparison_root <- file.path(project_root, "new_method_comparison")
formal_root <- file.path(comparison_root, "formal_simulation")
out_root <- file.path(
  formal_root, "main_simulation",
  "w6_onbasis_K6_beta0_0p5_intercept_svd_R100"
)

dir_names <- c("config", "metadata", "truth", "data", "results", "reports",
               "logs", "fits/production", "scripts")
invisible(lapply(file.path(out_root, dir_names), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

source(file.path(comparison_root, "R", "methods", "accelerated",
                 "newssgl_intercept.R"))
load_newssgl_intercept_fast()

p <- 10L
K_B <- 6L
raw_H <- K_B^2L
R <- 100L
true_beta0 <- 0.5
n_train <- 700L
n_test <- 250L
sigma_epsilon <- 0.35
grid_size_eval <- 50L
grid_size_reference <- 101L
theta_true <- c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0)
active_u <- c(3L, 4L, 5L, 6L)
predictor_names <- paste0("X", seq_len(p))
predictor_group <- c("global_only", "global_only", "spatial_only",
                     "spatial_only", "global_plus_spatial",
                     "global_plus_spatial", rep("null", 4))
lambda0 <- 20
lambda1 <- 2
n_iter <- 5000L
burn_in <- 2000L
n_chains <- 3L

seed_ledger <- data.frame(
  replicate = seq_len(R),
  data_seed = 2026100001L + seq_len(R),
  mcmc_seed = 2026410000L + seq_len(R),
  stringsAsFactors = FALSE
)

write.csv(seed_ledger, file.path(out_root, "metadata",
                                 "w6_beta0_0p5_seed_ledger.csv"),
          row.names = FALSE)

writeLines(c(
  "scenario: W6_onbasis_K6_beta0_0p5_intercept_svd_R100",
  "replicates: 100",
  "method: GD-SSGL only",
  sprintf("true_beta0: %.1f", true_beta0),
  sprintf("n_train: %d", n_train),
  sprintf("n_test: %d", n_test),
  sprintf("p: %d", p),
  sprintf("sigma_epsilon: %.2f", sigma_epsilon),
  sprintf("K_B: %d", K_B),
  sprintf("raw_H: %d", raw_H),
  "effective_H: 35",
  sprintf("lambda0: %.0f", lambda0),
  sprintf("lambda1: %.0f", lambda1),
  sprintf("n_iter: %d", n_iter),
  sprintf("burn_in: %d", burn_in),
  sprintf("chains: %d", n_chains),
  "thin: 1",
  "sampler: fit_newssgl_intercept_fast direct scalar beta0",
  "basis: SVD full-rank centered basis, Phi_full = Phi_c %*% V_r",
  "truth: W6 on-basis truth projected into each replicate training-basis SVD space",
  paste0("theta_true: [", paste(theta_true, collapse = ", "), "]"),
  "active_spatial_deviation_set: [3, 4, 5, 6]"
), file.path(out_root, "config", "w6_beta0_0p5_intercept_svd_manifest.yml"))

atomic_write_csv <- function(x, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Atomic rename failed: ", path)
}

atomic_save_rds <- function(x, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Atomic rename failed: ", path)
}

g1 <- function(s1, s2) sin(2 * pi * s1) * cos(2 * pi * s2)
g2 <- function(s1, s2) cos(2 * pi * s1) * sin(2 * pi * s2)
g3 <- function(s1, s2) exp(-((s1 - 0.35)^2 + (s2 - 0.65)^2) / 0.04)
g4 <- function(s1, s2) {
  exp(-((s1 - 0.70)^2 + (s2 - 0.30)^2) / 0.03) -
    exp(-((s1 - 0.25)^2 + (s2 - 0.25)^2) / 0.03)
}
analytic_funs <- list(`3` = g1, `4` = g2, `5` = g3, `6` = g4)

eval_grid <- as.matrix(expand.grid(
  s1 = seq(0, 1, length.out = grid_size_eval),
  s2 = seq(0, 1, length.out = grid_size_eval)
))
ref_grid <- as.matrix(expand.grid(
  s1 = seq(0, 1, length.out = grid_size_reference),
  s2 = seq(0, 1, length.out = grid_size_reference)
))
ref_w <- rep(1 / nrow(ref_grid), nrow(ref_grid))

weighted_rms <- function(x, w) sqrt(sum(w * x^2))

weighted_cor <- function(x, y, w) {
  denom <- weighted_rms(x, w) * weighted_rms(y, w)
  if (denom < 1e-14) return(NA_real_)
  sum(w * x * y) / denom
}

pinv_solve <- function(A, y) {
  s <- svd(A)
  tol <- max(dim(A)) * max(s$d) * .Machine$double.eps
  rank <- sum(s$d > tol)
  coef <- s$v[, seq_len(rank), drop = FALSE] %*%
    ((t(s$u[, seq_len(rank), drop = FALSE]) %*% y) / s$d[seq_len(rank)])
  list(coef = as.numeric(coef), rank = rank, tol = tol)
}

make_random_components <- function(seed) {
  set.seed(seed)
  n <- n_train + n_test
  coords <- cbind(runif(n), runif(n))
  X_raw <- matrix(rnorm(n * p), n, p)
  epsilon <- rnorm(n, sd = sigma_epsilon)
  train_id <- seq_len(n_train)
  test_id <- n_train + seq_len(n_test)
  x_mean <- colMeans(X_raw[train_id, , drop = FALSE])
  x_sd <- apply(X_raw[train_id, , drop = FALSE], 2, sd)
  X <- sweep(sweep(X_raw, 2, x_mean, "-"), 2, x_sd, "/")
  list(coords = coords, X = X, X_raw = X_raw, epsilon = epsilon,
       train_id = train_id, test_id = test_id, x_mean = x_mean, x_sd = x_sd)
}

construct_w6_truth <- function(train_coords, basis_meta) {
  B_ref <- apply_centered_basis_metadata_fullrank(ref_grid, basis_meta)
  h_eff <- ncol(B_ref)
  alpha_star <- matrix(0, h_eff, p)
  truth_rows <- vector("list", length(active_u))
  for (ii in seq_along(active_u)) {
    j <- active_u[ii]
    g <- analytic_funs[[as.character(j)]](ref_grid[, 1], ref_grid[, 2])
    g_c <- g - sum(ref_w * g)
    fit0 <- pinv_solve(sqrt(ref_w) * B_ref, sqrt(ref_w) * g_c)
    h0 <- as.vector(B_ref %*% fit0$coef)
    if (weighted_rms(h0, ref_w) < 1e-12) stop("Projection RMS is zero for X", j)
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
  list(alpha_star = alpha_star, validation = do.call(rbind, truth_rows))
}

make_dataset <- function(r) {
  row <- seed_ledger[seed_ledger$replicate == r, ]
  comp <- make_random_components(row$data_seed)
  train_id <- comp$train_id
  test_id <- comp$test_id
  train_coords <- comp$coords[train_id, , drop = FALSE]
  test_coords <- comp$coords[test_id, , drop = FALSE]
  X_train <- comp$X[train_id, , drop = FALSE]
  X_test <- comp$X[test_id, , drop = FALSE]
  basis_meta <- fit_centered_basis_metadata_fullrank(
    train_coords, K_B, tolerance = 1e-10, method = "svd"
  )
  B_train <- apply_centered_basis_metadata_fullrank(train_coords, basis_meta)
  B_test <- apply_centered_basis_metadata_fullrank(test_coords, basis_meta)
  B_eval <- apply_centered_basis_metadata_fullrank(eval_grid, basis_meta)
  h_eff <- ncol(B_train)
  truth <- construct_w6_truth(train_coords, basis_meta)
  u_train <- B_train %*% truth$alpha_star
  u_test <- B_test %*% truth$alpha_star
  u_eval <- B_eval %*% truth$alpha_star
  beta_train <- sweep(u_train, 2, theta_true, "+")
  beta_test <- sweep(u_test, 2, theta_true, "+")
  beta_eval <- sweep(u_eval, 2, theta_true, "+")
  eps_train <- comp$epsilon[train_id]
  eps_test <- comp$epsilon[test_id]
  y_train <- true_beta0 + rowSums(X_train * beta_train) + eps_train
  y_test <- true_beta0 + rowSums(X_test * beta_test) + eps_test
  response_identity <- max(abs(c(
    y_train - (true_beta0 + rowSums(X_train * beta_train) + eps_train),
    y_test - (true_beta0 + rowSums(X_test * beta_test) + eps_test)
  )))
  onbasis_error <- max(abs(rbind(
    B_train %*% truth$alpha_star - u_train,
    B_test %*% truth$alpha_star - u_test,
    B_eval %*% truth$alpha_star - u_eval
  )[, active_u, drop = FALSE]))
  validation <- data.frame(
    replicate = r,
    true_beta0 = true_beta0,
    effective_H = h_eff,
    raw_H = basis_meta$raw_h,
    centered_rank = basis_meta$centered_rank,
    rank_loss_from_centering = basis_meta$rank_loss_from_centering,
    training_basis_centering_max_error = max(abs(colMeans(B_train))),
    training_u_centering_max_error =
      max(abs(colMeans(u_train[, active_u, drop = FALSE]))),
    training_beta_mean_equals_theta_max_error =
      max(abs(colMeans(beta_train[, active_u, drop = FALSE]) - theta_true[active_u])),
    active_set_identity =
      identical(which(colSums(abs(u_eval)) > 1e-12), active_u),
    onbasis_max_error = onbasis_error,
    response_identity_max_error = response_identity,
    finite_values = all(is.finite(unlist(list(y_train, y_test, X_train, X_test,
                                             beta_eval, u_eval)))),
    dimension_ok = nrow(X_train) == n_train && nrow(X_test) == n_test &&
      ncol(X_train) == p && h_eff == raw_H - 1L,
    stringsAsFactors = FALSE
  )
  validation$pass <- with(validation,
    effective_H == 35L &&
      training_basis_centering_max_error < 1e-10 &&
      training_u_centering_max_error < 1e-10 &&
      training_beta_mean_equals_theta_max_error < 1e-10 &&
      active_set_identity &&
      onbasis_max_error < 1e-12 &&
      response_identity_max_error < 1e-12 &&
      finite_values &&
      dimension_ok
  )
  data_obj <- list(
    train = list(y = y_train, X = X_train, coords = train_coords),
    test = list(y = y_test, X = X_test, coords = test_coords),
    grid = list(coords = eval_grid, true_beta = beta_eval,
                true_u = u_eval, grid_size = grid_size_eval),
    truth = list(theta = theta_true, beta0 = true_beta0,
                 spatial_deviation = seq_len(p) %in% active_u,
                 active_surface = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                                    FALSE, FALSE, FALSE, FALSE),
                 alpha_star_svd = truth$alpha_star,
                 basis = basis_meta,
                 truth_validation = truth$validation),
    preprocessing = list(x_mean = comp$x_mean, x_sd = comp$x_sd),
    config = list(n_train = n_train, n_test = n_test, p = p,
                  sigma_epsilon = sigma_epsilon, K_B = K_B,
                  raw_H = raw_H, effective_H = h_eff,
                  data_seed = row$data_seed, mcmc_seed = row$mcmc_seed)
  )
  list(data = data_obj, validation = validation, truth_validation = truth$validation)
}

qint <- function(mat) {
  rbind(lower = apply(mat, 1, quantile, 0.025, names = FALSE),
        upper = apply(mat, 1, quantile, 0.975, names = FALSE))
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

fit_one <- function(data_obj, seed) {
  started <- proc.time()[3]
  chains <- lapply(seq_len(n_chains), function(ch) {
    fit_newssgl_intercept_fast(
      data_obj$train, data_obj$test, data_obj$grid,
      list(n_basis = K_B, full_rank_centered = TRUE, full_rank_method = "svd"),
      list(lambda0 = lambda0, lambda1 = lambda1, a_sigma = 0.5,
           b_sigma = var(data_obj$train$y) / 2,
           a_theta = 1, b_theta = 1, a_gamma = 1, b_gamma = 10),
      list(n_iter = n_iter, burn_in = burn_in),
      seed + ch
    )
  })
  meta <- chains[[1]]$config$basis
  theta <- do.call(cbind, lapply(chains, `[[`, "theta_draws"))
  alpha <- do.call(cbind, lapply(chains, function(x) x$diagnostics$alpha_draws))
  gamma <- do.call(cbind, lapply(chains, function(x) x$diagnostics$gamma_draws))
  gamma_prob <- do.call(cbind, lapply(chains, function(x) x$diagnostics$gamma_prob_draws))
  beta0 <- unlist(lapply(chains, `[[`, "beta0_draws"))
  theta_mean <- rowMeans(theta)
  alpha_mean <- rowMeans(alpha)
  beta0_mean <- mean(beta0)
  phi_test <- apply_centered_basis_metadata_fullrank(data_obj$test$coords, meta)
  phi_grid <- apply_centered_basis_metadata_fullrank(data_obj$grid$coords, meta)
  h <- ncol(phi_grid)
  beta_test <- matrix(0, nrow(data_obj$test$X), p)
  beta_grid <- matrix(0, nrow(data_obj$grid$coords), p)
  for (j in seq_len(p)) {
    idx <- ((j - 1L) * h + 1L):(j * h)
    beta_test[, j] <- theta_mean[j] + as.vector(phi_test %*% alpha_mean[idx])
    beta_grid[, j] <- theta_mean[j] + as.vector(phi_grid %*% alpha_mean[idx])
  }
  fit <- list(
    method_name = "GD-SSGL intercept/SVD W6",
    theta_mean = theta_mean,
    theta_draws = theta,
    pip = rowMeans(gamma),
    beta_hat_grid = beta_grid,
    u_hat_grid = sweep(beta_grid, 2, theta_mean, "-"),
    pred_test = beta0_mean + rowSums(data_obj$test$X * beta_test),
    runtime = proc.time()[3] - started,
    diagnostics = list(
      beta0_mean = beta0_mean,
      beta0_draws = beta0,
      rb_pip = rowMeans(gamma_prob),
      gamma_draws = gamma,
      gamma_prob_draws = gamma_prob,
      alpha_draws = alpha,
      chain_count = n_chains
    ),
    config = list(basis = meta, mcmc = list(n_iter = n_iter, burn_in = burn_in,
                                            chains = n_chains)),
    beta0_mean = beta0_mean,
    beta0_draws = beta0,
    alpha_mean = alpha_mean
  )
  fit$runtime_sec_external <- fit$runtime
  fit
}

postprocess_fit <- function(r, data_obj, fit, runtime_sec) {
  phi_grid <- apply_centered_basis_metadata_fullrank(data_obj$grid$coords,
                                                     fit$config$basis)
  h_eff <- ncol(phi_grid)
  alpha_draws <- fit$diagnostics$alpha_draws
  beta_mean_all <- matrix(0, nrow(data_obj$grid$coords), p)
  u_mean_all <- matrix(0, nrow(data_obj$grid$coords), p)
  beta_rows <- vector("list", p)
  u_rows <- vector("list", length(active_u))
  band_rows <- list()
  for (j in seq_len(p)) {
    idx <- ((j - 1L) * h_eff + 1L):(j * h_eff)
    u_draws <- phi_grid %*% alpha_draws[idx, , drop = FALSE]
    beta_draws <- sweep(u_draws, 2, fit$theta_draws[j, ], "+")
    u_mean <- rowMeans(u_draws)
    beta_mean <- rowMeans(beta_draws)
    u_mean_all[, j] <- u_mean
    beta_mean_all[, j] <- beta_mean
    bq <- qint(beta_draws)
    beta_rows[[j]] <- data.frame(
      replicate = r,
      predictor = predictor_names[j],
      predictor_group = predictor_group[j],
      beta_pointwise_coverage_95 =
        mean(data_obj$grid$true_beta[, j] >= bq["lower", ] &
               data_obj$grid$true_beta[, j] <= bq["upper", ]),
      beta_mean_interval_width = mean(bq["upper", ] - bq["lower", ]),
      beta_mean_posterior_sd = mean(apply(beta_draws, 1, sd)),
      beta_mise = mean((beta_mean - data_obj$grid$true_beta[, j])^2),
      stringsAsFactors = FALSE
    )
    if (j %in% active_u) {
      uq <- qint(u_draws)
      u_rows[[match(j, active_u)]] <- data.frame(
        replicate = r,
        predictor = predictor_names[j],
        predictor_group = predictor_group[j],
        u_pointwise_coverage_95 =
          mean(data_obj$grid$true_u[, j] >= uq["lower", ] &
                 data_obj$grid$true_u[, j] <= uq["upper", ]),
        u_mean_interval_width = mean(uq["upper", ] - uq["lower", ]),
        u_mean_posterior_sd = mean(apply(u_draws, 1, sd)),
        u_mise = mean((u_mean - data_obj$grid$true_u[, j])^2),
        stringsAsFactors = FALSE
      )
      ub <- band_summary_for(u_draws, data_obj$grid$true_u[, j])
      bb <- band_summary_for(beta_draws, data_obj$grid$true_beta[, j])
      band_rows[[length(band_rows) + 1L]] <- data.frame(
        replicate = r, predictor = predictor_names[j],
        predictor_group = predictor_group[j],
        estimand = "u",
        simultaneous_coverage = ub$simultaneous_coverage,
        simultaneous_mean_width = ub$simultaneous_mean_width,
        simultaneous_crit = ub$simultaneous_crit,
        stringsAsFactors = FALSE
      )
      band_rows[[length(band_rows) + 1L]] <- data.frame(
        replicate = r, predictor = predictor_names[j],
        predictor_group = predictor_group[j],
        estimand = "beta",
        simultaneous_coverage = bb$simultaneous_coverage,
        simultaneous_mean_width = bb$simultaneous_mean_width,
        simultaneous_crit = bb$simultaneous_crit,
        stringsAsFactors = FALSE
      )
    }
  }
  theta_q <- rbind(
    lower = apply(fit$theta_draws, 1, quantile, 0.025, names = FALSE),
    upper = apply(fit$theta_draws, 1, quantile, 0.975, names = FALSE)
  )
  theta_rows <- data.frame(
    replicate = r,
    predictor = predictor_names,
    predictor_group = predictor_group,
    theta_covered_95 = theta_true >= theta_q["lower", ] &
      theta_true <= theta_q["upper", ],
    theta_interval_width = theta_q["upper", ] - theta_q["lower", ],
    theta_posterior_sd = apply(fit$theta_draws, 1, sd),
    stringsAsFactors = FALSE
  )
  beta_mise_j <- colMeans((beta_mean_all - data_obj$grid$true_beta)^2)
  u_mise_j <- colMeans((u_mean_all - data_obj$grid$true_u)^2)
  rb_pip <- fit$diagnostics$rb_pip
  selected <- rb_pip >= 0.5
  active <- seq_len(p) %in% active_u
  pip_rows <- data.frame(
    replicate = r,
    predictor = predictor_names,
    target = as.integer(active),
    rb_pip = rb_pip,
    sampled_pip = fit$pip,
    selected = selected,
    stringsAsFactors = FALSE
  )
  selection_row <- data.frame(
    replicate = r,
    selected_predictors = paste(predictor_names[selected], collapse = ", "),
    TPR = sum(selected & active) / sum(active),
    FPR = sum(selected & !active) / sum(!active),
    FDR = if (sum(selected) > 0) sum(selected & !active) / sum(selected) else 0,
    exact_recovery = identical(selected, active),
    selected_size = sum(selected),
    false_positives = sum(selected & !active),
    false_negatives = sum(!selected & active),
    stringsAsFactors = FALSE
  )
  replicate_row <- data.frame(
    replicate = r,
    method = "GD-SSGL intercept/SVD W6",
    true_beta0 = true_beta0,
    beta0_hat = fit$beta0_mean,
    mspe = mean((data_obj$test$y - fit$pred_test)^2),
    rmse = sqrt(mean((data_obj$test$y - fit$pred_test)^2)),
    beta_mise = mean(beta_mise_j),
    u_mise_x3_x6 = mean(u_mise_j[active_u]),
    theta_mse = mean((fit$theta_mean - theta_true)^2),
    runtime_sec = runtime_sec,
    retained_draws = ncol(fit$theta_draws),
    selected_predictors = selection_row$selected_predictors,
    TPR = selection_row$TPR,
    FPR = selection_row$FPR,
    FDR = selection_row$FDR,
    exact_recovery = selection_row$exact_recovery,
    selected_size = selection_row$selected_size,
    stringsAsFactors = FALSE
  )
  list(replicate = replicate_row,
       beta = do.call(rbind, beta_rows),
       u = do.call(rbind, u_rows),
       theta = theta_rows,
       bands = do.call(rbind, band_rows),
       pip = pip_rows,
       selection = selection_row)
}

append_registry <- function(row) {
  path <- file.path(out_root, "results", "w6_beta0_0p5_fit_registry.csv")
  old <- if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else NULL
  new <- if (is.null(old)) row else rbind(old[old$replicate != row$replicate, ], row)
  atomic_write_csv(new[order(new$replicate), ], path)
}

for (r in seq_len(R)) {
  data_path <- file.path(out_root, "data", sprintf("w6_beta0_0p5_rep_%03d.rds", r))
  validation_path <- file.path(out_root, "metadata",
                               sprintf("w6_beta0_0p5_rep_%03d_validation.csv", r))
  truth_path <- file.path(out_root, "truth",
                          sprintf("w6_beta0_0p5_truth_rep_%03d.csv", r))
  fit_path <- file.path(out_root, "fits", "production",
                        sprintf("w6_beta0_0p5_rep_%03d_gdssgl_fit.rds", r))
  if (!file.exists(data_path) || !file.exists(validation_path) ||
      !file.exists(truth_path)) {
    made <- make_dataset(r)
    if (!all(made$validation$pass)) {
      atomic_write_csv(made$validation, validation_path)
      stop("Validation failed for replicate ", r)
    }
    tv <- made$truth_validation
    tv$replicate <- r
    atomic_save_rds(made$data, data_path)
    atomic_write_csv(made$validation, validation_path)
    atomic_write_csv(tv, truth_path)
  }
  data_obj <- readRDS(data_path)
  status <- "success"
  err <- ""
  runtime <- NA_real_
  if (!file.exists(fit_path)) {
    message(sprintf("[%s] fitting W6 beta0=0.5 intercept/SVD replicate %03d / %03d",
                    format(Sys.time(), "%F %T"), r, R))
    fit <- tryCatch(fit_one(data_obj, seed_ledger$mcmc_seed[r]), error = identity)
    if (inherits(fit, "error")) {
      status <- "failed"
      err <- conditionMessage(fit)
    } else {
      runtime <- fit$runtime_sec_external
      atomic_save_rds(fit, fit_path)
    }
  } else {
    fit <- readRDS(fit_path)
    runtime <- fit$runtime_sec_external %||% fit$runtime %||% NA_real_
  }
  append_registry(data.frame(
    replicate = r,
    data_seed = seed_ledger$data_seed[r],
    mcmc_seed = seed_ledger$mcmc_seed[r],
    status = status,
    fit_path = fit_path,
    data_path = data_path,
    runtime_sec = runtime,
    retained_draws = if (exists("fit") && !inherits(fit, "error")) ncol(fit$theta_draws) else NA_integer_,
    error_message = err,
    stringsAsFactors = FALSE
  ))
  if (status != "success") warning("Fit failed for replicate ", r, ": ", err)
}

validation_all <- do.call(rbind, lapply(seq_len(R), function(r) {
  read.csv(file.path(out_root, "metadata",
                     sprintf("w6_beta0_0p5_rep_%03d_validation.csv", r)),
           stringsAsFactors = FALSE)
}))
truth_all <- do.call(rbind, lapply(seq_len(R), function(r) {
  read.csv(file.path(out_root, "truth",
                     sprintf("w6_beta0_0p5_truth_rep_%03d.csv", r)),
           stringsAsFactors = FALSE)
}))
atomic_write_csv(validation_all, file.path(out_root, "metadata",
                                           "w6_beta0_0p5_data_validation.csv"))
atomic_write_csv(truth_all, file.path(out_root, "truth",
                                      "w6_beta0_0p5_truth_validation_100rep.csv"))

registry <- read.csv(file.path(out_root, "results",
                               "w6_beta0_0p5_fit_registry.csv"),
                     stringsAsFactors = FALSE)
success_reps <- registry$replicate[registry$status == "success"]

objects <- list()
for (r in success_reps) {
  message(sprintf("[%s] postprocessing W6 replicate %03d / %03d",
                  format(Sys.time(), "%F %T"), r, R))
  data_obj <- readRDS(file.path(out_root, "data",
                                sprintf("w6_beta0_0p5_rep_%03d.rds", r)))
  fit <- readRDS(file.path(out_root, "fits", "production",
                           sprintf("w6_beta0_0p5_rep_%03d_gdssgl_fit.rds", r)))
  runtime <- registry$runtime_sec[registry$replicate == r][1]
  objects[[as.character(r)]] <- postprocess_fit(r, data_obj, fit, runtime)
}

bind_field <- function(field) do.call(rbind, lapply(objects, `[[`, field))
replicate_level <- bind_field("replicate")
beta <- bind_field("beta")
u <- bind_field("u")
theta <- bind_field("theta")
bands <- bind_field("bands")
pip <- bind_field("pip")
selection <- bind_field("selection")

summarize_cols <- function(dat, by, cols) {
  if (length(by) == 0L) {
    out <- data.frame(scope = "overall", stringsAsFactors = FALSE)
    for (cc in cols) {
      z <- dat[[cc]]
      out[[paste0(cc, "_mean")]] <- mean(z)
      out[[paste0(cc, "_sd")]] <- sd(z)
      out[[paste0(cc, "_mcse")]] <- sd(z) / sqrt(length(z))
      out[[paste0(cc, "_n")]] <- length(z)
    }
    return(out)
  }
  split_key <- interaction(dat[by], drop = TRUE, sep = " | ")
  pieces <- lapply(split(dat, split_key), function(x) {
    vals <- unlist(strsplit(as.character(interaction(x[1, by, drop = FALSE],
                                                     drop = TRUE, sep = " | ")),
                            " \\| "))
    out <- as.data.frame(as.list(vals), stringsAsFactors = FALSE)
    names(out) <- by
    for (cc in cols) {
      z <- x[[cc]]
      out[[paste0(cc, "_mean")]] <- mean(z)
      out[[paste0(cc, "_sd")]] <- sd(z)
      out[[paste0(cc, "_mcse")]] <- sd(z) / sqrt(length(z))
      out[[paste0(cc, "_n")]] <- length(z)
    }
    out
  })
  do.call(rbind, pieces)
}

beta_predictor_summary <- summarize_cols(
  beta, c("predictor", "predictor_group"),
  c("beta_pointwise_coverage_95", "beta_mean_interval_width",
    "beta_mean_posterior_sd", "beta_mise")
)
beta_group_summary <- summarize_cols(
  beta, "predictor_group",
  c("beta_pointwise_coverage_95", "beta_mean_interval_width",
    "beta_mean_posterior_sd", "beta_mise")
)
u_predictor_summary <- summarize_cols(
  u, c("predictor", "predictor_group"),
  c("u_pointwise_coverage_95", "u_mean_interval_width",
    "u_mean_posterior_sd", "u_mise")
)
theta_predictor_summary <- summarize_cols(
  theta, c("predictor", "predictor_group"),
  c("theta_covered_95", "theta_interval_width", "theta_posterior_sd")
)
band_summary <- summarize_cols(
  bands, c("estimand", "predictor", "predictor_group"),
  c("simultaneous_coverage", "simultaneous_mean_width")
)
pip_summary <- summarize_cols(
  pip, c("predictor", "target"),
  c("rb_pip", "sampled_pip", "selected")
)
selection_summary <- summarize_cols(
  selection, character(0),
  c("TPR", "FPR", "FDR", "exact_recovery", "selected_size",
    "false_positives", "false_negatives")
)
overall_summary <- summarize_cols(
  replicate_level, character(0),
  c("mspe", "rmse", "beta_mise", "u_mise_x3_x6", "theta_mse",
    "beta0_hat", "runtime_sec", "TPR", "FPR", "FDR", "exact_recovery",
    "selected_size")
)

atomic_write_csv(replicate_level, file.path(out_root, "results",
                                            "w6_beta0_0p5_replicate_level_100rep.csv"))
atomic_write_csv(beta, file.path(out_root, "results",
                                 "w6_beta0_0p5_beta_coverage_by_replicate_100rep.csv"))
atomic_write_csv(u, file.path(out_root, "results",
                              "w6_beta0_0p5_u_coverage_by_replicate_100rep.csv"))
atomic_write_csv(theta, file.path(out_root, "results",
                                  "w6_beta0_0p5_theta_coverage_by_replicate_100rep.csv"))
atomic_write_csv(bands, file.path(out_root, "results",
                                  "w6_beta0_0p5_simultaneous_band_coverage_100rep.csv"))
atomic_write_csv(pip, file.path(out_root, "results",
                                "w6_beta0_0p5_predictor_pips_by_replicate_100rep.csv"))
atomic_write_csv(selection, file.path(out_root, "results",
                                      "w6_beta0_0p5_selection_by_replicate_100rep.csv"))
atomic_write_csv(beta_predictor_summary, file.path(out_root, "results",
                                                   "w6_beta0_0p5_beta_coverage_predictor_summary_100rep.csv"))
atomic_write_csv(beta_group_summary, file.path(out_root, "results",
                                               "w6_beta0_0p5_beta_coverage_group_summary_100rep.csv"))
atomic_write_csv(u_predictor_summary, file.path(out_root, "results",
                                                "w6_beta0_0p5_u_coverage_predictor_summary_100rep.csv"))
atomic_write_csv(theta_predictor_summary, file.path(out_root, "results",
                                                    "w6_beta0_0p5_theta_coverage_predictor_summary_100rep.csv"))
atomic_write_csv(band_summary, file.path(out_root, "results",
                                         "w6_beta0_0p5_simultaneous_band_summary_100rep.csv"))
atomic_write_csv(pip_summary, file.path(out_root, "results",
                                        "w6_beta0_0p5_predictor_pip_summary_100rep.csv"))
atomic_write_csv(selection_summary, file.path(out_root, "results",
                                             "w6_beta0_0p5_selection_summary_100rep.csv"))
atomic_write_csv(overall_summary, file.path(out_root, "results",
                                           "w6_beta0_0p5_overall_summary_100rep.csv"))

truth_summary <- summarize_cols(
  truth_all, "predictor",
  c("original_RMS", "W6_RMS", "RMS_ratio",
    "pre_rescaling_projection_RMSE", "post_rescaling_RMSE",
    "shape_correlation", "effective_basis_rank")
)
atomic_write_csv(truth_summary, file.path(out_root, "results",
                                          "w6_beta0_0p5_truth_validation_summary_100rep.csv"))

active_beta <- beta_predictor_summary[beta_predictor_summary$predictor %in%
                                        c("X3", "X4", "X5", "X6"), ]
active_u_summary <- u_predictor_summary[u_predictor_summary$predictor %in%
                                          c("X3", "X4", "X5", "X6"), ]
x5x6_bands <- band_summary[band_summary$predictor %in% c("X5", "X6"), ]

report <- c(
  "# W6 on-basis beta0=0.5 intercept/SVD GD-SSGL R=100",
  "",
  paste0("- Output root: `", out_root, "`"),
  "- Method: GD-SSGL only.",
  "- No CV, no retuning, no other methods.",
  "- Truth: W6 on-basis, projected into the same SVD full-rank centered K=6 basis space used by the sampler.",
  "- DGP includes true constant intercept beta0 = 0.5.",
  sprintf("- Completed fits: %d / 100", length(success_reps)),
  sprintf("- Failed fits: %d", sum(registry$status != "success")),
  sprintf("- Validation passed: %d / 100", sum(validation_all$pass)),
  "",
  "## Overall summary",
  "",
  paste(capture.output(print(overall_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Active beta pointwise coverage",
  "",
  paste(capture.output(print(active_beta, row.names = FALSE)), collapse = "\n"),
  "",
  "## Active u pointwise coverage",
  "",
  paste(capture.output(print(active_u_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## X5/X6 simultaneous bands",
  "",
  paste(capture.output(print(x5x6_bands, row.names = FALSE)), collapse = "\n"),
  "",
  "## PIP summary",
  "",
  paste(capture.output(print(pip_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Selection summary",
  "",
  paste(capture.output(print(selection_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Truth validation summary",
  "",
  paste(capture.output(print(truth_summary, row.names = FALSE)), collapse = "\n")
)
writeLines(report, file.path(out_root, "reports",
                             "w6_beta0_0p5_intercept_svd_R100_report.md"))
cat(paste(report, collapse = "\n"))
