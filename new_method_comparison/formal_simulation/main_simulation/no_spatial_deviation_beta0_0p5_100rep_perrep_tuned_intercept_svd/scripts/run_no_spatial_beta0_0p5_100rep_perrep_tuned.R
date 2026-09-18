#!/usr/bin/env Rscript

options(digits = 17)
suppressPackageStartupMessages({
  library(Rcpp)
  library(GIGrvg)
  library(mgcv)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = "") {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
only_rep <- arg_value("replicate", "")
only_rep <- if (nzchar(only_rep)) as.integer(only_rep) else NA_integer_
aggregate_only <- any(args == "--aggregate-only")
workers <- max(1L, as.integer(arg_value("workers", "1")))

comparison_root <- normalizePath("new_method_comparison", winslash = "/",
                                 mustWork = TRUE)
formal_root <- file.path(comparison_root, "formal_simulation")
out_root <- file.path(formal_root, "main_simulation",
                      "no_spatial_deviation_beta0_0p5_100rep_perrep_tuned_intercept_svd")
dirs <- c("config", "data", "cv", "fits", "results", "reports", "logs", "scripts")
invisible(lapply(file.path(out_root, dirs), dir.create, recursive = TRUE,
                 showWarnings = FALSE))

Sys.setenv(SSGL_SOURCE_FUNCTIONS_ONLY = "1")
source(file.path(formal_root, "scripts",
                 "run_four_function_10rep_uncertainty_calibration_Kmax6.R"))
source(file.path(formal_root, "R", "tuning.R"))
source(file.path(formal_root, "configs", "frozen_tuning_protocol.R"))
source(file.path(comparison_root, "R", "methods", "accelerated",
                 "newssgl_intercept.R"))
source(file.path(comparison_root, "R", "methods", "wrappers.R"))
load_newssgl_intercept_fast()
Rcpp::sourceCpp(file.path(comparison_root, "R", "methods", "accelerated",
                          "original_ssgl_intercept.cpp"))
Rcpp::sourceCpp(file.path(comparison_root, "R", "methods", "accelerated",
                          "full_svc_intercept.cpp"))

true_beta0 <- 0.5
n_rep <- 100L
p <- 10L
K_values <- c(4L, 5L, 6L)
theta_true <- c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0)
predictor_names <- paste0("X", seq_len(p))
active_deviation <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                      FALSE, FALSE, FALSE, FALSE)
active_surface <- c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE,
                    FALSE, FALSE, FALSE, FALSE)
fold_seed <- 123L
cv_mcmc <- list(n_iter = 500L, burn_in = 100L)
final_mcmc <- list(
  proposed = list(n_iter = 10000L, burn_in = 3000L, chains = 1L),
  original = list(n_iter = 5000L, burn_in = 2000L, chains = 3L),
  full_svc = list(n_iter = 5000L, burn_in = 2000L, chains = 3L),
  blasso = list(n_iter = 5000L, burn_in = 2000L, chains = 3L)
)
protocol <- frozen_tuning_protocol

replicate_seeds <- data.frame(
  replicate = seq_len(n_rep),
  data_seed = 2026100001L + seq_len(n_rep),
  proposed_seed = 2026410000L + seq_len(n_rep),
  original_seed_base = 2026420000L + 100L * seq_len(n_rep),
  full_svc_seed_base = 2026430000L + 100L * seq_len(n_rep),
  blasso_seed_base = 2026440000L + 100L * seq_len(n_rep)
)
write.csv(replicate_seeds, file.path(out_root, "config", "replicate_seeds.csv"),
          row.names = FALSE)
saveRDS(list(
  true_beta0 = true_beta0,
  n_rep = n_rep,
  K_values = K_values,
  H_raw_values = K_values^2L,
  H_effective_centered = K_values^2L - 1L,
  fold_seed = fold_seed,
  cv_mcmc = cv_mcmc,
  final_mcmc = final_mcmc,
  grids = list(
    proposed = expand_basis_hyperparameter_grid(K_values, protocol$proposed_ssgl$candidates),
    original = expand_basis_hyperparameter_grid(K_values, protocol$original_ssgl$candidates),
    full_svc = expand_basis_hyperparameter_grid(
      K_values, data.frame(kappa2_alpha = protocol$full_svc_no_selection$kappa2_alpha)
    )
  ),
  intercept_rule = "direct scalar beta0 with p(beta0) proportional to 1",
  basis_rule = "SVD full-rank centered basis for GD-SSGL and Gaussian SVC"
), file.path(out_root, "config", "manifest.rds"))

atomic_save_rds <- function(x, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Atomic rename failed: ", path)
}
atomic_write_csv <- function(x, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Atomic rename failed: ", path)
}

make_data <- function(seed) {
  dat <- generate_four_function_data(seed)
  old_u_train <- make_u_from_training(dat$train$coords, dat$train$coords)$u
  old_u_test <- make_u_from_training(dat$train$coords, dat$test$coords)$u
  old_beta_train <- sweep(old_u_train, 2, theta_true, "+")
  old_beta_test <- sweep(old_u_test, 2, theta_true, "+")
  train_mu <- as.vector(dat$train$X %*% theta_true)
  test_mu <- as.vector(dat$test$X %*% theta_true)
  train_eps <- dat$train$y - rowSums(dat$train$X * old_beta_train)
  test_eps <- dat$test$y - rowSums(dat$test$X * old_beta_test)
  dat$train$y <- true_beta0 + train_mu + train_eps
  dat$test$y <- true_beta0 + test_mu + test_eps
  dat$train$true_beta <- matrix(rep(theta_true, each = nrow(dat$train$X)),
                                nrow(dat$train$X), p)
  dat$test$true_beta <- matrix(rep(theta_true, each = nrow(dat$test$X)),
                               nrow(dat$test$X), p)
  dat$grid$true_u <- matrix(0, nrow(dat$grid$coords), p)
  dat$grid$true_beta <- matrix(rep(theta_true, each = nrow(dat$grid$coords)),
                               nrow(dat$grid$coords), p)
  dat$truth$theta <- theta_true
  dat$truth$spatial_deviation <- active_deviation
  dat$truth$active_surface <- active_surface
  dat$truth$beta0 <- true_beta0
  dat$truth$epsilon_train <- train_eps
  dat$truth$epsilon_test <- test_eps
  dat
}

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

cv_data <- function(dat, folds, fold) cv_data_objects_leakage_free(dat, folds, fold)

run_cv_fit <- function(dat, folds, method, candidate_id, params, fold, rep_id) {
  fd <- cv_data(dat, folds, fold)
  K <- as.integer(params$K)
  H_raw <- K^2L
  seed <- switch(
    method,
    proposed_ssgl = 8810000L,
    original_ssgl = 8820000L,
    full_svc_no_selection = 8830000L
  ) + 100000L * rep_id + 1000L * candidate_id + fold
  model_config <- list(a_sigma = 0.5, b_sigma = var(fd$train$y) / 2,
                       a_theta = 1, b_theta = 1,
                       a_gamma = 1, b_gamma = 10)
  fit_started <- proc.time()[3]
  fit <- tryCatch({
    if (method == "proposed_ssgl") {
      model_config$lambda0 <- params$lambda0
      model_config$lambda1 <- params$lambda1
      fit_newssgl_intercept_fast(
        fd$train, fd$validation, list(coords = fd$validation$coords),
        list(n_basis = K, full_rank_centered = TRUE, full_rank_method = "svd"),
        model_config, cv_mcmc, seed
      )
    } else if (method == "original_ssgl") {
      model_config$lambda0 <- params$lambda0
      model_config$lambda1 <- params$lambda1
      model_config$zeta0 <- protocol$original_ssgl$zeta0
      model_config$zeta1 <- protocol$original_ssgl$zeta1
      fit_original_ssgl_intercept(
        fd$train, fd$validation, list(coords = fd$validation$coords),
        list(n_basis = K),
        model_config, cv_mcmc, seed
      )
    } else {
      model_config$kappa2_alpha <- params$kappa2_alpha
      fit_full_svc_no_selection_intercept(
        fd$train, fd$validation, list(coords = fd$validation$coords),
        list(n_basis = K, full_rank_centered = TRUE, full_rank_method = "svd"),
        model_config, cv_mcmc, seed
      )
    }
  }, error = identity)
  success <- !inherits(fit, "error")
  data.frame(
    replicate = rep_id,
    method = method,
    candidate_id = candidate_id,
    K = K,
    H_raw = H_raw,
    H_effective = if (method %in% c("proposed_ssgl", "full_svc_no_selection"))
      H_raw - 1L else H_raw,
    lambda0 = if ("lambda0" %in% names(params)) params$lambda0 else NA_real_,
    lambda1 = if ("lambda1" %in% names(params)) params$lambda1 else NA_real_,
    kappa2_alpha = if ("kappa2_alpha" %in% names(params))
      params$kappa2_alpha else NA_real_,
    fold = fold,
    validation_mse = if (success) mean((fd$validation$y - fit$pred_test)^2) else Inf,
    fit_runtime_sec = proc.time()[3] - fit_started,
    convergence_status = if (success) "success" else "failed",
    error_message = if (success) "" else conditionMessage(fit),
    n_iter = cv_mcmc$n_iter,
    burn_in = cv_mcmc$burn_in,
    stringsAsFactors = FALSE
  )
}

summarize_cv <- function(fit_results) {
  ok <- fit_results[fit_results$convergence_status == "success", , drop = FALSE]
  split_key <- paste(ok$method, ok$candidate_id, sep = "__")
  out <- do.call(rbind, lapply(split(ok, split_key), function(x) {
    data.frame(
      replicate = x$replicate[1], method = x$method[1],
      candidate_id = x$candidate_id[1],
      K = x$K[1], H_raw = x$H_raw[1], H_effective = x$H_effective[1],
      lambda0 = x$lambda0[1], lambda1 = x$lambda1[1],
      kappa2_alpha = x$kappa2_alpha[1],
      mean_validation_mse = mean(x$validation_mse),
      sd_validation_mse = sd(x$validation_mse),
      se_validation_mse = sd(x$validation_mse) / sqrt(nrow(x)),
      total_cv_runtime_sec = sum(x$fit_runtime_sec),
      successful_folds = nrow(x),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$method, out$mean_validation_mse, out$K, out$candidate_id), ]
}

pick_best <- function(x) {
  min_mse <- min(x$mean_validation_mse)
  tied <- x[abs(x$mean_validation_mse - min_mse) <= 1e-10, , drop = FALSE]
  if (x$method[1] %in% c("proposed_ssgl", "original_ssgl")) {
    tied <- tied[order(tied$K, -tied$lambda0, -tied$lambda1,
                       tied$candidate_id), , drop = FALSE]
  } else {
    tied <- tied[order(tied$K, tied$kappa2_alpha, tied$candidate_id), , drop = FALSE]
  }
  tied[1, , drop = FALSE]
}

tune_one_replicate <- function(rep_id, dat, rep_out) {
  cv_done <- file.path(rep_out, "cv_selected_configurations.csv")
  if (file.exists(cv_done)) {
    return(list(
      fit_results = read.csv(file.path(rep_out, "cv_fit_results.csv"), stringsAsFactors = FALSE),
      summary = read.csv(file.path(rep_out, "cv_configuration_summary.csv"), stringsAsFactors = FALSE),
      selected = read.csv(cv_done, stringsAsFactors = FALSE)
    ))
  }
  folds <- make_fixed_cv_folds(nrow(dat$train$X), 5L, fold_seed)
  atomic_write_csv(data.frame(row_id = seq_along(folds), fold = folds),
                   file.path(rep_out, "cv_fold_assignment_seed123.csv"))
  grids <- list(
    proposed_ssgl = expand_basis_hyperparameter_grid(
      K_values, protocol$proposed_ssgl$candidates
    ),
    original_ssgl = expand_basis_hyperparameter_grid(
      K_values, protocol$original_ssgl$candidates
    ),
    full_svc_no_selection = expand_basis_hyperparameter_grid(
      K_values, data.frame(kappa2_alpha = protocol$full_svc_no_selection$kappa2_alpha)
    )
  )
  rows <- list()
  row_id <- 1L
  fit_results_path <- file.path(rep_out, "cv_fit_results.csv")
  existing <- if (file.exists(fit_results_path)) {
    read.csv(fit_results_path, stringsAsFactors = FALSE)
  } else data.frame()
  done <- if (nrow(existing)) {
    paste(existing$method, existing$candidate_id, existing$fold, sep = "__")
  } else character()
  if (nrow(existing)) {
    rows <- split(existing, seq_len(nrow(existing)))
    row_id <- length(rows) + 1L
  }
  for (method in names(grids)) {
    grid <- grids[[method]]
    for (candidate_id in seq_len(nrow(grid))) {
      params <- as.list(grid[candidate_id, , drop = FALSE])
      for (fold in sort(unique(folds))) {
        key <- paste(method, candidate_id, fold, sep = "__")
        if (key %in% done) next
        row <- run_cv_fit(dat, folds, method, candidate_id, params, fold, rep_id)
        rows[[row_id]] <- row
        row_id <- row_id + 1L
        fit_results <- do.call(rbind, rows)
        atomic_write_csv(fit_results, fit_results_path)
        cat(sprintf(
          "rep %03d %s cand %02d fold %d: %s mse=%.6f %.2fs\n",
          rep_id, method, candidate_id, fold, row$convergence_status,
          row$validation_mse, row$fit_runtime_sec
        ))
      }
    }
  }
  fit_results <- read.csv(fit_results_path, stringsAsFactors = FALSE)
  config_summary <- summarize_cv(fit_results)
  selected <- do.call(rbind, lapply(split(config_summary, config_summary$method),
                                    pick_best))
  atomic_write_csv(config_summary, file.path(rep_out, "cv_configuration_summary.csv"))
  atomic_write_csv(selected, cv_done)
  list(fit_results = fit_results, summary = config_summary, selected = selected)
}

fit_gd_final <- function(dat, cfg, seed) {
  fit_newssgl_intercept_fast(
    dat$train, dat$test, dat$grid,
    list(n_basis = cfg$K, full_rank_centered = TRUE, full_rank_method = "svd"),
    list(lambda0 = cfg$lambda0, lambda1 = cfg$lambda1,
         a_sigma = 0.5, b_sigma = var(dat$train$y) / 2,
         a_theta = 1, b_theta = 1, a_gamma = 1, b_gamma = 10),
    list(n_iter = final_mcmc$proposed$n_iter,
         burn_in = final_mcmc$proposed$burn_in),
    seed
  )
}

fit_ws_final <- function(dat, cfg, seed_base) {
  started <- proc.time()[3]
  chains <- lapply(1:3, function(ch) {
    fit_original_ssgl_intercept(
      dat$train, dat$test, dat$grid,
      list(n_basis = cfg$K),
      list(lambda0 = cfg$lambda0, lambda1 = cfg$lambda1,
           a_gamma = 1, b_gamma = 10, a_sigma = 0.5,
           b_sigma = var(dat$train$y) / 2,
           zeta0 = protocol$original_ssgl$zeta0,
           zeta1 = protocol$original_ssgl$zeta1),
      list(n_iter = final_mcmc$original$n_iter,
           burn_in = final_mcmc$original$burn_in),
      seed_base + ch
    )
  })
  meta <- chains[[1]]$meta
  eta <- do.call(cbind, lapply(chains, `[[`, "eta"))
  gamma <- do.call(cbind, lapply(chains, `[[`, "gamma"))
  gamma_prob <- do.call(cbind, lapply(chains, `[[`, "gamma_prob"))
  beta0 <- unlist(lapply(chains, `[[`, "beta0_draws"))
  eta_mean <- rowMeans(eta)
  beta0_mean <- mean(beta0)
  phi_test <- apply_basis_metadata(dat$test$coords, meta)
  phi_grid <- apply_basis_metadata(dat$grid$coords, meta)
  beta_test <- surface_eta(eta_mean, phi_test, p)
  beta_grid <- surface_eta(eta_mean, phi_grid, p)
  list(meta = meta, eta = eta, gamma = gamma, gamma_prob = gamma_prob,
       beta0_draws = beta0, beta0_mean = beta0_mean,
       beta_grid = beta_grid,
       pred_test = beta0_mean + rowSums(dat$test$X * beta_test),
       runtime = proc.time()[3] - started)
}

fit_svc_final <- function(dat, cfg, seed_base) {
  started <- proc.time()[3]
  chains <- lapply(1:3, function(ch) {
    fit_full_svc_no_selection_intercept(
      dat$train, dat$test, dat$grid,
      list(n_basis = cfg$K, full_rank_centered = TRUE, full_rank_method = "svd"),
      list(kappa2_alpha = cfg$kappa2_alpha,
           a_sigma = 0.5, b_sigma = var(dat$train$y) / 2,
           a_theta = 1, b_theta = 1),
      list(n_iter = final_mcmc$full_svc$n_iter,
           burn_in = final_mcmc$full_svc$burn_in),
      seed_base + ch
    )
  })
  meta <- chains[[1]]$config$basis
  theta <- do.call(cbind, lapply(chains, `[[`, "theta_draws"))
  alpha <- do.call(cbind, lapply(chains, function(x) x$diagnostics$alpha_draws))
  beta0 <- unlist(lapply(chains, `[[`, "beta0_draws"))
  theta_mean <- rowMeans(theta)
  alpha_mean <- rowMeans(alpha)
  beta0_mean <- mean(beta0)
  phi_test <- apply_phi(dat$test$coords, meta)
  phi_grid <- apply_phi(dat$grid$coords, meta)
  beta_test <- surface_theta_alpha(theta_mean, alpha_mean, phi_test, p)
  beta_grid <- surface_theta_alpha(theta_mean, alpha_mean, phi_grid, p)
  list(meta = meta, theta = theta, alpha = alpha, beta0_draws = beta0,
       beta0_mean = beta0_mean, theta_mean = theta_mean,
       beta_grid = beta_grid, u_grid = sweep(beta_grid, 2, theta_mean, "-"),
       pred_test = beta0_mean + rowSums(dat$test$X * beta_test),
       runtime = proc.time()[3] - started)
}

fit_blasso_final <- function(dat, seed_base) {
  started <- proc.time()[3]
  chains <- lapply(1:3, function(ch) {
    fit_global_only_blasso_intercept(
      dat$train, dat$test, dat$grid,
      list(),
      list(a_sigma = 0.5, b_sigma = var(dat$train$y) / 2,
           a_theta = 1, b_theta = 1),
      list(n_iter = final_mcmc$blasso$n_iter,
           burn_in = final_mcmc$blasso$burn_in),
      seed_base + ch
    )
  })
  theta <- do.call(cbind, lapply(chains, `[[`, "theta_draws"))
  beta0 <- unlist(lapply(chains, `[[`, "beta0_draws"))
  theta_mean <- rowMeans(theta)
  beta0_mean <- mean(beta0)
  beta_grid <- matrix(rep(theta_mean, each = nrow(dat$grid$coords)),
                      nrow(dat$grid$coords), p)
  list(theta = theta, beta0_draws = beta0, beta0_mean = beta0_mean,
       theta_mean = theta_mean, beta_grid = beta_grid,
       u_grid = matrix(0, nrow(dat$grid$coords), p),
       pred_test = beta0_mean + as.vector(dat$test$X %*% theta_mean),
       runtime = proc.time()[3] - started)
}

fit_gam_final <- function(dat) {
  started <- proc.time()[3]
  train_df <- as.data.frame(dat$train$X)
  names(train_df) <- predictor_names
  train_df$y <- dat$train$y
  train_df$s1 <- dat$train$coords[, 1]
  train_df$s2 <- dat$train$coords[, 2]
  form <- as.formula(paste0(
    "y ~ ", paste(predictor_names, collapse = " + "), " + ",
    paste0("s(s1, s2, bs = 'tp', by = ", predictor_names, ")",
           collapse = " + ")
  ))
  fit <- mgcv::gam(form, data = train_df, method = "REML")
  test_df <- as.data.frame(dat$test$X)
  names(test_df) <- predictor_names
  test_df$s1 <- dat$test$coords[, 1]
  test_df$s2 <- dat$test$coords[, 2]
  pred_test <- as.numeric(predict(fit, newdata = test_df))
  base_grid <- data.frame(s1 = dat$grid$coords[, 1], s2 = dat$grid$coords[, 2])
  for (jj in seq_len(p)) base_grid[[predictor_names[jj]]] <- 0
  baseline_grid <- as.numeric(predict(fit, newdata = base_grid))
  beta_grid <- matrix(0, nrow(dat$grid$coords), p)
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

metric_row <- function(rep_id, method, pred_test, beta_grid, theta_hat,
                       u_grid, dat, runtime_sec, beta0_hat,
                       config_text, selected_by_cv) {
  beta_error <- colMeans((beta_grid - dat$grid$true_beta)^2)
  data.frame(
    replicate = rep_id,
    method = method,
    true_beta0 = true_beta0,
    beta0_hat = beta0_hat,
    beta0_error = beta0_hat - true_beta0,
    mspe = mean((dat$test$y - pred_test)^2),
    rmse = sqrt(mean((dat$test$y - pred_test)^2)),
    mae = mean(abs(dat$test$y - pred_test)),
    bias = mean(dat$test$y - pred_test),
    beta_mise = mean(beta_error),
    mise_global_only = mean(beta_error[c(1:2, 5:6)]),
    mise_spatial_only = mean(beta_error[3:4]),
    mise_global_plus_spatial = mean(beta_error[5:6]),
    mise_null = mean(beta_error[c(3:4, 7:10)]),
    mise_true_global_only = mean(beta_error[c(1:2, 5:6)]),
    mise_true_null = mean(beta_error[c(3:4, 7:10)]),
    theta_mse = if (is.null(theta_hat)) NA_real_ else
      mean((theta_hat - theta_true)^2),
    u_mise_x3_x6 = if (is.null(u_grid)) NA_real_ else
      mean(colMeans((u_grid[, 3:6, drop = FALSE] -
                       dat$grid$true_u[, 3:6, drop = FALSE])^2)),
    u_mise_all = if (is.null(u_grid)) NA_real_ else
      mean(colMeans((u_grid - dat$grid$true_u)^2)),
    runtime_sec = runtime_sec,
    final_config = config_text,
    selected_by_cv = selected_by_cv,
    stringsAsFactors = FALSE
  )
}

run_one_replicate <- function(rep_id) {
  rep_out <- file.path(out_root, "fits", sprintf("rep_%03d", rep_id))
  dir.create(rep_out, recursive = TRUE, showWarnings = FALSE)
  done_file <- file.path(rep_out, "replicate_complete.rds")
  if (file.exists(done_file)) {
    cat("Skipping complete replicate ", rep_id, "\n", sep = "")
    return(readRDS(done_file))
  }
  log_file <- file.path(out_root, "logs", sprintf("rep_%03d.log", rep_id))
  sink(log_file, append = TRUE, split = TRUE)
  on.exit(sink(), add = TRUE)
  seeds <- replicate_seeds[replicate_seeds$replicate == rep_id, ]
  dat <- make_data(seeds$data_seed)
  atomic_save_rds(dat, file.path(rep_out, "shared_dataset_beta0_0p5.rds"))
  validation <- data.frame(
    check = c("true_beta0_0p5", "n_train_700", "n_test_250",
              "p_10", "active_deviation_set", "active_surface_set",
              "true_u_zero", "response_identity"),
    passed = c(
      identical(dat$truth$beta0, true_beta0),
      nrow(dat$train$X) == 700L,
      nrow(dat$test$X) == 250L,
      ncol(dat$train$X) == 10L,
      identical(dat$truth$spatial_deviation, active_deviation),
      identical(dat$truth$active_surface, active_surface),
      max(abs(dat$grid$true_u)) < 1e-12,
      max(abs(dat$train$y - (true_beta0 + as.vector(dat$train$X %*% theta_true) +
                              dat$truth$epsilon_train))) < 1e-12
    )
  )
  atomic_write_csv(validation, file.path(rep_out, "validation_checks.csv"))
  if (!all(validation$passed)) stop("Validation failed for replicate ", rep_id)

  cat(sprintf("[%s] rep %03d tuning\n", format(Sys.time(), "%F %T"), rep_id))
  tuning <- tune_one_replicate(rep_id, dat, rep_out)
  selected <- tuning$selected
  get_cfg <- function(method) selected[selected$method == method, , drop = FALSE][1, ]

  metrics <- list()
  pips <- list()

  gd_path <- file.path(rep_out, "gdssgl_final_fit.rds")
  cfg <- get_cfg("proposed_ssgl")
  if (file.exists(gd_path)) gd <- readRDS(gd_path) else {
    gd <- fit_gd_final(dat, cfg, seeds$proposed_seed)
    atomic_save_rds(gd, gd_path)
  }
  metrics[[length(metrics) + 1L]] <- metric_row(
    rep_id, "GD-SSGL", gd$pred_test, gd$beta_hat_grid, gd$theta_mean,
    gd$u_hat_grid, dat, gd$runtime, gd$beta0_mean,
    sprintf("K=%d,H_eff=%d,lambda0=%g,lambda1=%g",
            cfg$K, cfg$H_effective, cfg$lambda0, cfg$lambda1),
    TRUE
  )
  pips[[length(pips) + 1L]] <- data.frame(
    replicate = rep_id, method = "GD-SSGL", predictor = predictor_names,
    target = active_deviation, sampled_pip = gd$pip,
    rb_pip = gd$diagnostics$rb_pip,
    selection_target = "spatial deviation"
  )

  ws_path <- file.path(rep_out, "wsssgl_final_fit.rds")
  cfg <- get_cfg("original_ssgl")
  if (file.exists(ws_path)) ws <- readRDS(ws_path) else {
    ws <- fit_ws_final(dat, cfg, seeds$original_seed_base)
    atomic_save_rds(ws, ws_path)
  }
  metrics[[length(metrics) + 1L]] <- metric_row(
    rep_id, "WS-SSGL", ws$pred_test, ws$beta_grid, colMeans(ws$beta_grid),
    sweep(ws$beta_grid, 2, colMeans(ws$beta_grid), "-"), dat,
    ws$runtime, ws$beta0_mean,
    sprintf("K=%d,H=%d,lambda0=%g,lambda1=%g",
            cfg$K, cfg$H_effective, cfg$lambda0, cfg$lambda1),
    TRUE
  )
  pips[[length(pips) + 1L]] <- data.frame(
    replicate = rep_id, method = "WS-SSGL", predictor = predictor_names,
    target = active_surface, sampled_pip = rowMeans(ws$gamma),
    rb_pip = rowMeans(ws$gamma_prob),
    selection_target = "whole coefficient surface"
  )

  svc_path <- file.path(rep_out, "gaussian_svc_final_fit.rds")
  cfg <- get_cfg("full_svc_no_selection")
  if (file.exists(svc_path)) svc <- readRDS(svc_path) else {
    svc <- fit_svc_final(dat, cfg, seeds$full_svc_seed_base)
    atomic_save_rds(svc, svc_path)
  }
  metrics[[length(metrics) + 1L]] <- metric_row(
    rep_id, "Gaussian SVC", svc$pred_test, svc$beta_grid, svc$theta_mean,
    svc$u_grid, dat, svc$runtime, svc$beta0_mean,
    sprintf("K=%d,H_eff=%d,kappa2_alpha=%g",
            cfg$K, cfg$H_effective, cfg$kappa2_alpha),
    TRUE
  )

  bl_path <- file.path(rep_out, "bayesian_lasso_final_fit.rds")
  if (file.exists(bl_path)) bl <- readRDS(bl_path) else {
    bl <- fit_blasso_final(dat, seeds$blasso_seed_base)
    atomic_save_rds(bl, bl_path)
  }
  metrics[[length(metrics) + 1L]] <- metric_row(
    rep_id, "Bayesian Lasso", bl$pred_test, bl$beta_grid, bl$theta_mean,
    bl$u_grid, dat, bl$runtime, bl$beta0_mean,
    "internal shrinkage,no spatial basis", FALSE
  )

  gam_path <- file.path(rep_out, "gam_final_fit.rds")
  if (file.exists(gam_path)) gam <- readRDS(gam_path) else {
    gam <- fit_gam_final(dat)
    atomic_save_rds(gam, gam_path)
  }
  metrics[[length(metrics) + 1L]] <- metric_row(
    rep_id, "Standard GAM + intercept", gam$pred_test, gam$beta_grid,
    gam$theta_mean, gam$u_grid, dat, gam$runtime, gam$beta0_mean,
    "mgcv REML default k with formula intercept", FALSE
  )

  metrics <- do.call(rbind, metrics)
  pips <- do.call(rbind, pips)
  atomic_write_csv(metrics, file.path(rep_out, "final_metrics.csv"))
  atomic_write_csv(pips, file.path(rep_out, "final_pips.csv"))
  out <- list(replicate = rep_id, validation = validation,
              selected = selected, metrics = metrics, pips = pips)
  atomic_save_rds(out, done_file)
  cat(sprintf("[%s] rep %03d complete\n", format(Sys.time(), "%F %T"), rep_id))
  out
}

summ_stats <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(n = 0, mean = NA, sd = NA, mcse = NA))
  c(n = length(x), mean = mean(x), sd = sd(x), mcse = sd(x) / sqrt(length(x)))
}

aggregate_outputs <- function() {
  files <- file.path(out_root, "fits", sprintf("rep_%03d", seq_len(n_rep)),
                     "replicate_complete.rds")
  if (!all(file.exists(files))) {
    stop("Missing complete reps: ", paste(which(!file.exists(files)), collapse = ", "))
  }
  reps <- lapply(files, readRDS)
  metrics <- do.call(rbind, lapply(reps, `[[`, "metrics"))
  pips <- do.call(rbind, lapply(reps, `[[`, "pips"))
  selected <- do.call(rbind, lapply(reps, `[[`, "selected"))
  atomic_write_csv(metrics, file.path(out_root, "results",
                                      "no_spatial_deviation_beta0_0p5_100rep_metrics_by_replicate.csv"))
  atomic_write_csv(pips, file.path(out_root, "results",
                                   "no_spatial_deviation_beta0_0p5_100rep_pips_by_replicate.csv"))
  atomic_write_csv(selected, file.path(out_root, "results",
                                       "no_spatial_deviation_beta0_0p5_100rep_selected_configs.csv"))

  metric_cols <- c("beta0_hat", "beta0_error", "mspe", "rmse", "mae", "bias",
                   "beta_mise", "mise_global_only", "mise_spatial_only",
                   "mise_global_plus_spatial", "mise_null",
                   "mise_true_global_only", "mise_true_null", "theta_mse",
                   "u_mise_x3_x6", "u_mise_all", "runtime_sec")
  metric_summary <- do.call(rbind, lapply(split(metrics, metrics$method), function(x) {
    do.call(rbind, lapply(metric_cols, function(v) {
      data.frame(method = x$method[1], metric = v,
                 as.list(summ_stats(x[[v]])), stringsAsFactors = FALSE)
    }))
  }))
  atomic_write_csv(metric_summary, file.path(out_root, "results",
                                             "no_spatial_deviation_beta0_0p5_100rep_metric_summary.csv"))

  pips$selected <- as.integer(pips$rb_pip >= 0.5)
  pip_summary <- aggregate(cbind(sampled_pip, rb_pip, selected) ~
                             method + predictor + target + selection_target,
                           pips, function(x) c(mean = mean(x), sd = sd(x)))
  pip_summary <- do.call(data.frame, pip_summary)
  names(pip_summary) <- sub("\\.", "_", names(pip_summary))
  atomic_write_csv(pip_summary, file.path(out_root, "results",
                                          "no_spatial_deviation_beta0_0p5_100rep_pip_summary.csv"))

  sel_summary <- do.call(rbind, lapply(split(selected, selected$method), function(x) {
    data.frame(
      method = x$method[1],
      selected_K_table = paste(capture.output(print(table(x$K))), collapse = " "),
      selected_lambda0_table = if (all(is.na(x$lambda0))) "" else
        paste(capture.output(print(table(x$lambda0))), collapse = " "),
      selected_lambda1_table = if (all(is.na(x$lambda1))) "" else
        paste(capture.output(print(table(x$lambda1))), collapse = " "),
      selected_kappa2_table = if (all(is.na(x$kappa2_alpha))) "" else
        paste(capture.output(print(table(x$kappa2_alpha))), collapse = " "),
      stringsAsFactors = FALSE
    )
  }))
  atomic_write_csv(sel_summary, file.path(out_root, "results",
                                          "no_spatial_deviation_beta0_0p5_100rep_config_frequency_summary.csv"))

  sink(file.path(out_root, "reports", "no_spatial_deviation_beta0_0p5_100rep_perrep_tuned_report.md"))
  cat("# no-spatial-deviation beta0=0.5 100-rep per-replicate tuning diagnostic\n\n")
  cat("Every replicate independently tunes GD-SSGL, WS-SSGL, and Gaussian SVC using 5-fold CV.\n\n")
  cat("## Selected configurations\n\n")
  print(selected)
  cat("\n## Configuration frequencies\n\n")
  print(sel_summary)
  cat("\n## Metric summary\n\n")
  print(metric_summary)
  cat("\n## PIP summary\n\n")
  print(pip_summary)
  sink()
  invisible(list(metrics = metrics, pips = pips, selected = selected,
                 metric_summary = metric_summary, pip_summary = pip_summary))
}

if (!aggregate_only) {
  targets <- if (is.na(only_rep)) seq_len(n_rep) else only_rep
  if (workers > 1L && length(targets) > 1L &&
      .Platform$OS.type != "windows") {
    invisible(parallel::mclapply(targets, run_one_replicate, mc.cores = workers,
                                 mc.preschedule = FALSE, mc.set.seed = FALSE))
  } else {
    for (r in targets) run_one_replicate(r)
  }
}
if (is.na(only_rep)) {
  aggregate_outputs()
  cat("Saved beta0=0.5 100-rep per-rep tuned simulation to:\n",
      out_root, "\n", sep = "")
}
