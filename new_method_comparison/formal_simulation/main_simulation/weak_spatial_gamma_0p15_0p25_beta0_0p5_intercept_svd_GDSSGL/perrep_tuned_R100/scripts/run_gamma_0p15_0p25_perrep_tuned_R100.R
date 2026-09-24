#!/usr/bin/env Rscript

options(digits = 17, warn = 1)

suppressPackageStartupMessages({
  library(Rcpp)
  library(GIGrvg)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = "") {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
workers <- max(1L, as.integer(arg_value("workers", "4")))

project_root <- normalizePath(Sys.getenv("GDSSGL_ROOT", getwd()), winslash = "/", mustWork = TRUE)
comparison_root <- file.path(project_root, "new_method_comparison")
formal_root <- file.path(comparison_root, "formal_simulation")
fixed_root <- file.path(
  formal_root, "main_simulation",
  "weak_spatial_gamma_0p15_0p25_beta0_0p5_intercept_svd_GDSSGL",
  "smoke_R2"
)
out_root <- file.path(
  formal_root, "main_simulation",
  "weak_spatial_gamma_0p15_0p25_beta0_0p5_intercept_svd_GDSSGL",
  "perrep_tuned_R100"
)
reuse_r10_root <- file.path(
  formal_root, "main_simulation",
  "weak_spatial_gamma_0p15_0p25_beta0_0p5_intercept_svd_GDSSGL",
  "perrep_tuned_R10"
)

dirs <- file.path(out_root, c("config", "cv", "data", "fits", "results",
                              "reports", "logs", "scripts"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

source(file.path(comparison_root, "R", "methods", "wrappers.R"))
source(file.path(comparison_root, "R", "methods", "accelerated",
                 "newssgl_intercept.R"))
source(file.path(formal_root, "R", "tuning.R"))
source(file.path(formal_root, "configs", "frozen_tuning_protocol.R"))
load_newssgl_intercept_fast()

`%||%` <- function(x, y) if (is.null(x)) y else x

gamma_grid <- seq(0.15, 0.25, by = 0.01)
replicates <- 1:100
p <- 10L
n_train <- 700L
n_test <- 250L
sigma_epsilon <- 0.35
true_beta0 <- 0.5
theta_true <- c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0)
active_u <- c(3L, 4L, 5L, 6L)
active_deviation <- seq_len(p) %in% active_u
predictor_names <- paste0("X", seq_len(p))
K_values <- c(4L, 5L, 6L)
fold_seed <- 123L
cv_mcmc <- list(n_iter = 500L, burn_in = 100L)
final_mcmc <- list(n_iter = 5000L, burn_in = 2000L, chains = 3L)
grid_size <- 50L
grid_coords <- as.matrix(expand.grid(
  s1 = seq(0, 1, length.out = grid_size),
  s2 = seq(0, 1, length.out = grid_size)
))
grid <- expand_basis_hyperparameter_grid(
  K_values, frozen_tuning_protocol$proposed_ssgl$candidates
)

writeLines(c(
  "scenario: weak_gamma_beta0_0p5_intercept_svd_GDSSGL_perrep_tuned_R100",
  "method: GD-SSGL only",
  "data_source: rep1-10 reused from perrep_tuned_R10 when available; rep11-100 generated with the same seed rule",
  paste0("source_root: ", fixed_root),
  paste0("reuse_R10_root: ", file.path(formal_root, "main_simulation", "weak_spatial_gamma_0p15_0p25_beta0_0p5_intercept_svd_GDSSGL", "perrep_tuned_R10")),
  paste0("gamma_grid: [", paste(sprintf("%.2f", gamma_grid), collapse = ", "), "]"),
  "replicates: [1, ..., 100]",
  "per_replicate_tuning: true",
  "K_values: [4, 5, 6]",
  "lambda_grid: frozen_tuning_protocol$proposed_ssgl$candidates",
  "cv_folds: 5",
  paste0("cv_seed: ", fold_seed),
  paste0("cv_n_iter: ", cv_mcmc$n_iter),
  paste0("cv_burn_in: ", cv_mcmc$burn_in),
  paste0("final_n_iter: ", final_mcmc$n_iter),
  paste0("final_burn_in: ", final_mcmc$burn_in),
  "basis: SVD full-rank centered basis",
  paste0("true_beta0: ", true_beta0),
  paste0("workers: ", workers)
), file.path(out_root, "config", "manifest.yml"))

gamma_tag <- function(gamma) sub("[.]", "p", sprintf("gamma_%0.3f", gamma))

atomic_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Atomic rename failed: ", path)
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Atomic rename failed: ", path)
}

g1 <- function(s1, s2) sin(2 * pi * s1) * cos(2 * pi * s2)
g2 <- function(s1, s2) cos(2 * pi * s1) * sin(2 * pi * s2)
g3 <- function(s1, s2) exp(-((s1 - 0.35)^2 + (s2 - 0.65)^2) / 0.04)
g4 <- function(s1, s2) {
  exp(-((s1 - 0.70)^2 + (s2 - 0.30)^2) / 0.03) -
    exp(-((s1 - 0.25)^2 + (s2 - 0.25)^2) / 0.03)
}
four_functions <- list(`3` = g1, `4` = g2, `5` = g3, `6` = g4)

make_base_components <- function(rep_id) {
  data_seed <- 2026100001L + rep_id
  epsilon_seed <- 2026200001L + rep_id
  set.seed(data_seed)
  n <- n_train + n_test
  coords <- cbind(runif(n), runif(n))
  X_raw <- matrix(rnorm(n * p), n, p)
  train_id <- seq_len(n_train)
  test_id <- n_train + seq_len(n_test)
  x_mean <- colMeans(X_raw[train_id, , drop = FALSE])
  x_sd <- apply(X_raw[train_id, , drop = FALSE], 2, sd)
  X <- sweep(sweep(X_raw, 2, x_mean, "-"), 2, x_sd, "/")
  set.seed(epsilon_seed)
  epsilon <- rnorm(n, sd = sigma_epsilon)
  train_coords <- coords[train_id, , drop = FALSE]
  make_u_unit <- function(target_coords) {
    u_unit <- matrix(0, nrow(target_coords), p)
    for (j in active_u) {
      f <- four_functions[[as.character(j)]]
      raw_train <- f(train_coords[, 1], train_coords[, 2])
      center <- mean(raw_train)
      norm <- sqrt(mean((raw_train - center)^2))
      raw_target <- f(target_coords[, 1], target_coords[, 2])
      u_unit[, j] <- (raw_target - center) / norm
    }
    u_unit
  }
  list(
    coords = coords,
    X = X,
    epsilon = epsilon,
    train_id = train_id,
    test_id = test_id,
    preprocessing = list(x_mean = x_mean, x_sd = x_sd),
    u_unit_train = make_u_unit(coords[train_id, , drop = FALSE]),
    u_unit_test = make_u_unit(coords[test_id, , drop = FALSE]),
    u_unit_grid = make_u_unit(grid_coords),
    seeds = list(data_seed = data_seed, epsilon_seed = epsilon_seed)
  )
}

make_gamma_dataset <- function(base, gamma) {
  u_train <- gamma * base$u_unit_train
  u_test <- gamma * base$u_unit_test
  u_grid <- gamma * base$u_unit_grid
  beta_train <- sweep(u_train, 2, theta_true, "+")
  beta_test <- sweep(u_test, 2, theta_true, "+")
  beta_grid <- sweep(u_grid, 2, theta_true, "+")
  X_train <- base$X[base$train_id, , drop = FALSE]
  X_test <- base$X[base$test_id, , drop = FALSE]
  eps_train <- base$epsilon[base$train_id]
  eps_test <- base$epsilon[base$test_id]
  y_train <- true_beta0 + rowSums(X_train * beta_train) + eps_train
  y_test <- true_beta0 + rowSums(X_test * beta_test) + eps_test
  list(
    train = list(y = y_train, X = X_train,
                 coords = base$coords[base$train_id, , drop = FALSE]),
    test = list(y = y_test, X = X_test,
                coords = base$coords[base$test_id, , drop = FALSE]),
    grid = list(coords = grid_coords, true_beta = beta_grid,
                true_u = u_grid, grid_size = grid_size),
    truth = list(beta0 = true_beta0, theta = theta_true, gamma = gamma,
                 spatial_deviation = active_deviation,
                 active_surface = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                                    FALSE, FALSE, FALSE, FALSE)),
    preprocessing = base$preprocessing,
    config = list(n_train = n_train, n_test = n_test, p = p,
                  sigma_epsilon = sigma_epsilon, true_beta0 = true_beta0,
                  K_B = 6L, H_raw = 36L, H_effective = 35L,
                  gamma = gamma, data_seed = base$seeds$data_seed,
                  epsilon_seed = base$seeds$epsilon_seed)
  )
}

load_dataset <- function(gamma, rep_id) {
  source_path <- file.path(
    fixed_root, "data", gamma_tag(gamma),
    sprintf("weak_gamma_beta0_0p5_rep_%03d.rds", rep_id)
  )
  out_path <- file.path(
    out_root, "data", gamma_tag(gamma),
    sprintf("weak_gamma_beta0_0p5_rep_%03d.rds", rep_id)
  )
  if (file.exists(out_path)) return(readRDS(out_path))
  if (file.exists(source_path)) {
    dat <- readRDS(source_path)
  } else {
    dat <- make_gamma_dataset(make_base_components(rep_id), gamma)
  }
  atomic_save_rds(dat, out_path)
  dat
}

run_cv_fit <- function(dat, folds, candidate_id, params, fold, gamma, rep_id) {
  fd <- cv_data_objects_leakage_free(dat, folds, fold)
  K <- as.integer(params$K)
  seed <- 2030200000L + as.integer(round(gamma * 1000)) * 100000L +
    10000L * rep_id + 100L * candidate_id + fold
  model_config <- list(
    lambda0 = params$lambda0,
    lambda1 = params$lambda1,
    a_sigma = 0.5,
    b_sigma = var(fd$train$y) / 2,
    a_theta = 1,
    b_theta = 1,
    a_gamma = 1,
    b_gamma = 10
  )
  fit_started <- proc.time()[3]
  fit <- tryCatch(
    fit_newssgl_intercept_fast(
      fd$train, fd$validation, list(coords = fd$validation$coords),
      list(n_basis = K, full_rank_centered = TRUE, full_rank_method = "svd"),
      model_config,
      cv_mcmc,
      seed
    ),
    error = identity
  )
  success <- !inherits(fit, "error")
  data.frame(
    gamma = gamma,
    replicate = rep_id,
    candidate_id = candidate_id,
    K = K,
    H_raw = K^2L,
    H_effective = K^2L - 1L,
    lambda0 = params$lambda0,
    lambda1 = params$lambda1,
    fold = fold,
    validation_mse = if (success) mean((fd$validation$y - fit$pred_test)^2) else Inf,
    fit_runtime_sec = proc.time()[3] - fit_started,
    status = if (success) "success" else "failed",
    error_message = if (success) "" else conditionMessage(fit),
    stringsAsFactors = FALSE
  )
}

summarize_cv <- function(cv_rows) {
  ok <- cv_rows[cv_rows$status == "success", , drop = FALSE]
  out <- do.call(rbind, lapply(split(ok, ok$candidate_id), function(x) {
    data.frame(
      gamma = x$gamma[1],
      replicate = x$replicate[1],
      candidate_id = x$candidate_id[1],
      K = x$K[1],
      H_raw = x$H_raw[1],
      H_effective = x$H_effective[1],
      lambda0 = x$lambda0[1],
      lambda1 = x$lambda1[1],
      mean_validation_mse = mean(x$validation_mse),
      sd_validation_mse = sd(x$validation_mse),
      se_validation_mse = sd(x$validation_mse) / sqrt(nrow(x)),
      total_cv_runtime_sec = sum(x$fit_runtime_sec),
      successful_folds = nrow(x),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$mean_validation_mse, out$K, -out$lambda0, -out$lambda1,
            out$candidate_id), , drop = FALSE]
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

fit_final <- function(dat, cfg, gamma, rep_id) {
  seed_base <- 2040200000L + as.integer(round(gamma * 1000)) * 10000L + rep_id
  started <- proc.time()[3]
  chains <- lapply(seq_len(final_mcmc$chains), function(ch) {
    fit_newssgl_intercept_fast(
      dat$train, dat$test, dat$grid,
      list(n_basis = cfg$K, full_rank_centered = TRUE,
           full_rank_method = "svd"),
      list(lambda0 = cfg$lambda0, lambda1 = cfg$lambda1,
           a_sigma = 0.5, b_sigma = var(dat$train$y) / 2,
           a_theta = 1, b_theta = 1, a_gamma = 1, b_gamma = 10),
      list(n_iter = final_mcmc$n_iter, burn_in = final_mcmc$burn_in),
      seed_base + ch
    )
  })
  meta <- chains[[1]]$config$basis
  theta <- do.call(cbind, lapply(chains, `[[`, "theta_draws"))
  alpha <- do.call(cbind, lapply(chains, function(x) x$diagnostics$alpha_draws))
  gamma_draws <- do.call(cbind, lapply(chains, function(x) x$diagnostics$gamma_draws))
  gamma_prob <- do.call(cbind, lapply(chains, function(x) x$diagnostics$gamma_prob_draws))
  beta0 <- unlist(lapply(chains, `[[`, "beta0_draws"))
  theta_mean <- rowMeans(theta)
  alpha_mean <- rowMeans(alpha)
  beta0_mean <- mean(beta0)
  phi_test <- apply_phi(dat$test$coords, meta)
  phi_grid <- apply_phi(dat$grid$coords, meta)
  beta_test <- surface_theta_alpha(theta_mean, alpha_mean, phi_test, p)
  beta_grid <- surface_theta_alpha(theta_mean, alpha_mean, phi_grid, p)
  list(
    method_name = "GD-SSGL",
    theta_mean = theta_mean,
    theta_draws = theta,
    pip = rowMeans(gamma_draws),
    beta_hat_grid = beta_grid,
    u_hat_grid = sweep(beta_grid, 2, theta_mean, "-"),
    pred_test = beta0_mean + rowSums(dat$test$X * beta_test),
    runtime = proc.time()[3] - started,
    diagnostics = list(
      beta0_mean = beta0_mean,
      beta0_draws = beta0,
      rb_pip = rowMeans(gamma_prob),
      gamma_draws = gamma_draws,
      gamma_prob_draws = gamma_prob,
      alpha_draws = alpha,
      chain_count = final_mcmc$chains
    ),
    config = list(basis = meta, mcmc = final_mcmc),
    beta0_mean = beta0_mean,
    beta0_draws = beta0,
    alpha_mean = alpha_mean
  )
}

postprocess <- function(dat, fit, cfg, gamma, rep_id, final_runtime) {
  rb_pip <- fit$diagnostics$rb_pip
  selected <- rb_pip >= 0.5
  beta_mise_j <- colMeans((fit$beta_hat_grid - dat$grid$true_beta)^2)
  u_mise_j <- colMeans((fit$u_hat_grid - dat$grid$true_u)^2)
  list(
    metrics = data.frame(
      gamma = gamma,
      replicate = rep_id,
      K = cfg$K,
      H_raw = cfg$H_raw,
      H_effective = cfg$H_effective,
      lambda0 = cfg$lambda0,
      lambda1 = cfg$lambda1,
      cv_mse = cfg$mean_validation_mse,
      true_beta0 = true_beta0,
      beta0_hat = fit$beta0_mean,
      mspe = mean((dat$test$y - fit$pred_test)^2),
      beta_mise = mean(beta_mise_j),
      theta_mse = mean((fit$theta_mean - theta_true)^2),
      u_mise_x3_x6 = mean(u_mise_j[active_u]),
      selected_predictors = paste(predictor_names[selected], collapse = ", "),
      TPR = sum(selected & active_deviation) / sum(active_deviation),
      FPR = sum(selected & !active_deviation) / sum(!active_deviation),
      exact_recovery = identical(selected, active_deviation),
      selected_size = sum(selected),
      runtime_sec = final_runtime,
      stringsAsFactors = FALSE
    ),
    pips = data.frame(
      gamma = gamma,
      replicate = rep_id,
      predictor = predictor_names,
      target = active_deviation,
      rb_pip = rb_pip,
      sampled_pip = fit$pip,
      selected = selected,
      stringsAsFactors = FALSE
    )
  )
}

validate_dataset <- function(dat, gamma, rep_id) {
  ok <- TRUE
  messages <- character()
  if (length(dat$train$y) != n_train || length(dat$test$y) != n_test) {
    ok <- FALSE
    messages <- c(messages, "unexpected train/test sample size")
  }
  if (!isTRUE(all.equal(dat$truth$beta0, true_beta0))) {
    ok <- FALSE
    messages <- c(messages, "beta0 truth mismatch")
  }
  if (!isTRUE(all.equal(dat$truth$gamma, gamma))) {
    ok <- FALSE
    messages <- c(messages, "gamma truth mismatch")
  }
  if (!identical(as.logical(dat$truth$spatial_deviation), active_deviation)) {
    ok <- FALSE
    messages <- c(messages, "active deviation set mismatch")
  }
  if (any(!is.finite(dat$train$y)) || any(!is.finite(dat$test$y)) ||
      any(!is.finite(dat$train$X)) || any(!is.finite(dat$test$X))) {
    ok <- FALSE
    messages <- c(messages, "non-finite values")
  }
  u_means <- colMeans(dat$grid$true_u[, active_u, drop = FALSE])
  data.frame(
    gamma = gamma,
    replicate = rep_id,
    passed = ok,
    message = paste(messages, collapse = "; "),
    max_abs_grid_u_mean_x3_x6 = max(abs(u_means)),
    stringsAsFactors = FALSE
  )
}

run_one_task <- function(task) {
  gamma <- as.numeric(task$gamma)
  rep_id <- as.integer(task$replicate)
  tag <- gamma_tag(gamma)
  rep_dir <- file.path(out_root, "fits", tag, sprintf("rep_%03d", rep_id))
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  complete_path <- file.path(rep_dir, "complete.rds")
  if (file.exists(complete_path)) {
    return(readRDS(complete_path))
  }
  source_complete <- file.path(reuse_r10_root, "fits", tag,
                               sprintf("rep_%03d", rep_id), "complete.rds")
  if (file.exists(source_complete)) {
    reused <- readRDS(source_complete)
    atomic_save_rds(reused, complete_path)
    return(reused)
  }

  log_file <- file.path(out_root, "logs",
                        sprintf("%s_rep_%03d.log", tag, rep_id))
  log_con <- file(log_file, open = "at")
  sink(log_con, type = "output")
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(log_con)
  }, add = TRUE)

  status <- data.frame(
    gamma = gamma, replicate = rep_id, status = "started",
    error_message = "", stringsAsFactors = FALSE
  )
  status_path <- file.path(rep_dir, "task_status.csv")
  atomic_write_csv(status, status_path)

  result <- tryCatch({
    dat <- load_dataset(gamma, rep_id)
    validation <- validate_dataset(dat, gamma, rep_id)
    atomic_write_csv(validation, file.path(rep_dir, "data_validation.csv"))
    if (!isTRUE(validation$passed)) {
      stop("Dataset validation failed: ", validation$message)
    }
    folds <- make_fixed_cv_folds(nrow(dat$train$X), 5L, fold_seed)
    atomic_write_csv(data.frame(row_id = seq_along(folds), fold = folds),
                     file.path(rep_dir, "cv_fold_assignment_seed123.csv"))

    cv_path <- file.path(rep_dir, "cv_fit_results.csv")
    if (file.exists(cv_path)) {
      rows_done <- read.csv(cv_path, stringsAsFactors = FALSE)
    } else {
      rows_done <- data.frame()
    }
    rows <- if (nrow(rows_done)) split(rows_done, seq_len(nrow(rows_done))) else list()
    done_keys <- if (nrow(rows_done)) {
      paste(rows_done$candidate_id, rows_done$fold, sep = ":")
    } else character()
    row_id <- length(rows) + 1L
    for (candidate_id in seq_len(nrow(grid))) {
      params <- as.list(grid[candidate_id, , drop = FALSE])
      for (fold in sort(unique(folds))) {
        key <- paste(candidate_id, fold, sep = ":")
        if (key %in% done_keys) next
        row <- run_cv_fit(dat, folds, candidate_id, params, fold, gamma, rep_id)
        rows[[row_id]] <- row
        row_id <- row_id + 1L
        cv_current <- do.call(rbind, rows)
        atomic_write_csv(cv_current, cv_path)
        cat(sprintf("gamma %.2f rep %03d cand %02d fold %d: %s mse=%.6f %.2fs\n",
                    gamma, rep_id, candidate_id, fold, row$status,
                    row$validation_mse, row$fit_runtime_sec))
      }
    }
    cv_rows <- do.call(rbind, rows)
    cfg_summary <- summarize_cv(cv_rows)
    cfg <- cfg_summary[1, , drop = FALSE]
    atomic_write_csv(cfg_summary, file.path(rep_dir, "cv_configuration_summary.csv"))
    atomic_write_csv(cfg, file.path(rep_dir, "selected_configuration.csv"))

    final_path <- file.path(rep_dir, "gdssgl_tuned_final_fit.rds")
    started <- proc.time()[3]
    if (file.exists(final_path)) {
      fit <- readRDS(final_path)
      final_runtime <- fit$runtime %||% NA_real_
    } else {
      fit <- fit_final(dat, cfg, gamma, rep_id)
      final_runtime <- fit$runtime %||% (proc.time()[3] - started)
      atomic_save_rds(fit, final_path)
    }
    pp <- postprocess(dat, fit, cfg, gamma, rep_id, final_runtime)
    atomic_write_csv(pp$metrics, file.path(rep_dir, "tuned_final_metrics.csv"))
    atomic_write_csv(pp$pips, file.path(rep_dir, "tuned_final_pips.csv"))

    out <- list(
      validation = validation,
      cv = cv_rows,
      cfg = cfg,
      metrics = pp$metrics,
      pips = pp$pips,
      status = data.frame(
        gamma = gamma, replicate = rep_id, status = "success",
        error_message = "", stringsAsFactors = FALSE
      )
    )
    atomic_save_rds(out, complete_path)
    out
  }, error = function(e) {
    out <- list(
      validation = data.frame(gamma = gamma, replicate = rep_id,
                              passed = FALSE,
                              message = conditionMessage(e),
                              max_abs_grid_u_mean_x3_x6 = NA_real_),
      cv = data.frame(),
      cfg = data.frame(),
      metrics = data.frame(),
      pips = data.frame(),
      status = data.frame(
        gamma = gamma, replicate = rep_id, status = "failed",
        error_message = conditionMessage(e), stringsAsFactors = FALSE
      )
    )
    atomic_write_csv(out$status, status_path)
    atomic_save_rds(out, complete_path)
    out
  })
  atomic_write_csv(result$status, status_path)
  result
}

tasks <- expand.grid(gamma = gamma_grid, replicate = replicates)
task_list <- split(tasks, seq_len(nrow(tasks)))
if (workers > 1L) {
  results <- parallel::mclapply(task_list, run_one_task, mc.cores = workers)
} else {
  results <- lapply(task_list, run_one_task)
}

complete_files <- unlist(lapply(gamma_grid, function(gamma) {
  file.path(out_root, "fits", gamma_tag(gamma),
            sprintf("rep_%03d", replicates), "complete.rds")
}))
objects <- lapply(complete_files[file.exists(complete_files)], readRDS)

nonempty_bind <- function(objects, field) {
  pieces <- lapply(objects, `[[`, field)
  pieces <- pieces[vapply(pieces, nrow, integer(1)) > 0L]
  if (!length(pieces)) data.frame() else do.call(rbind, pieces)
}

validation_all <- nonempty_bind(objects, "validation")
status_all <- nonempty_bind(objects, "status")
cv_all <- nonempty_bind(objects, "cv")
cfg_all <- nonempty_bind(objects, "cfg")
metrics_all <- nonempty_bind(objects, "metrics")
pips_all <- nonempty_bind(objects, "pips")

summary_stat <- function(x) {
  data.frame(mean = mean(x), sd = sd(x), mcse = sd(x) / sqrt(length(x)))
}

pip_summary <- do.call(rbind, lapply(split(
  pips_all, list(pips_all$gamma, pips_all$predictor), drop = TRUE
), function(x) {
  data.frame(
    gamma = x$gamma[1],
    predictor = x$predictor[1],
    target = x$target[1],
    rb_pip_mean = mean(x$rb_pip),
    rb_pip_sd = sd(x$rb_pip),
    rb_pip_mcse = sd(x$rb_pip) / sqrt(nrow(x)),
    sampled_pip_mean = mean(x$sampled_pip),
    selected_frequency = mean(x$selected),
    n_replicates = nrow(x),
    stringsAsFactors = FALSE
  )
}))
pip_summary <- pip_summary[order(pip_summary$gamma,
                                 as.integer(sub("X", "", pip_summary$predictor))), ]

selection_summary <- do.call(rbind, lapply(split(metrics_all, metrics_all$gamma),
                                           function(x) {
  data.frame(
    gamma = x$gamma[1],
    TPR_mean = mean(x$TPR),
    FPR_mean = mean(x$FPR),
    FDR_mean = mean(ifelse(x$selected_size > 0,
                           (x$selected_size - x$TPR * sum(active_deviation)) /
                             x$selected_size, 0)),
    exact_recovery_frequency = mean(x$exact_recovery),
    selected_size_mean = mean(x$selected_size),
    selected_size_sd = sd(x$selected_size),
    MSPE_mean = mean(x$mspe),
    beta_MISE_mean = mean(x$beta_mise),
    u_MISE_x3_x6_mean = mean(x$u_mise_x3_x6),
    n_replicates = nrow(x),
    stringsAsFactors = FALSE
  )
}))
selection_summary <- selection_summary[order(selection_summary$gamma), ]

cfg_freq <- as.data.frame(table(cfg_all$gamma, cfg_all$K, cfg_all$lambda0,
                                cfg_all$lambda1), stringsAsFactors = FALSE)
names(cfg_freq) <- c("gamma", "K", "lambda0", "lambda1", "n")
cfg_freq <- cfg_freq[cfg_freq$n > 0, ]
cfg_freq$gamma <- as.numeric(as.character(cfg_freq$gamma))
cfg_freq$K <- as.integer(as.character(cfg_freq$K))
cfg_freq$lambda0 <- as.numeric(as.character(cfg_freq$lambda0))
cfg_freq$lambda1 <- as.numeric(as.character(cfg_freq$lambda1))
cfg_freq$frequency <- cfg_freq$n / length(replicates)
cfg_freq <- cfg_freq[order(cfg_freq$gamma, -cfg_freq$n, cfg_freq$K), ]

atomic_write_csv(validation_all, file.path(out_root, "results", "gdssgl_tuned_data_validation.csv"))
atomic_write_csv(status_all, file.path(out_root, "results", "gdssgl_tuned_status.csv"))
atomic_write_csv(cv_all, file.path(out_root, "results", "gdssgl_tuned_cv_fit_results.csv"))
atomic_write_csv(cfg_all, file.path(out_root, "results", "gdssgl_tuned_selected_configs.csv"))
atomic_write_csv(metrics_all, file.path(out_root, "results", "gdssgl_tuned_metrics_by_replicate.csv"))
atomic_write_csv(pips_all, file.path(out_root, "results", "gdssgl_tuned_pips_by_replicate.csv"))
atomic_write_csv(pip_summary, file.path(out_root, "results", "gdssgl_tuned_pip_summary.csv"))
atomic_write_csv(selection_summary, file.path(out_root, "results", "gdssgl_tuned_selection_summary.csv"))
atomic_write_csv(cfg_freq, file.path(out_root, "results", "gdssgl_tuned_config_frequency_summary.csv"))

active_pip <- pip_summary[pip_summary$predictor %in% paste0("X", active_u), ]
active_wide <- reshape(active_pip[, c("gamma", "predictor", "rb_pip_mean",
                                      "rb_pip_sd", "selected_frequency")],
                       idvar = "gamma", timevar = "predictor",
                       direction = "wide")
active_wide <- active_wide[order(active_wide$gamma), ]

completed <- sum(status_all$status == "success")
failed <- sum(status_all$status == "failed")
report <- c(
  "# GD-SSGL gamma 0.15-0.25 R100 per-replicate tuning",
  "",
  paste0("- Output root: `", out_root, "`"),
  paste0("- Completed tasks: ", completed, " / ", nrow(tasks)),
  paste0("- Failed tasks: ", failed),
  paste0("- Data validation passed: ", all(validation_all$passed)),
  "- Model: GD-SSGL only, direct intercept sampler, true beta0 = 0.5.",
  "- Basis: SVD full-rank centered basis.",
  "- Tuning: independent 5-fold CV for each gamma/replicate over K={4,5,6} and the 8 proposed SSGL lambda pairs.",
  "- Final fit: selected per-replicate configuration.",
  "",
  "## Active RB-PIP Means",
  "",
  paste(capture.output(print(active_wide, row.names = FALSE)), collapse = "\n"),
  "",
  "## Selection Summary",
  "",
  paste(capture.output(print(selection_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Configuration Frequencies",
  "",
  paste(capture.output(print(cfg_freq, row.names = FALSE)), collapse = "\n"),
  "",
  "## Metrics By Replicate",
  "",
  paste(capture.output(print(metrics_all, row.names = FALSE)), collapse = "\n")
)
writeLines(report, file.path(out_root, "reports",
                             "gdssgl_gamma_0p15_0p25_perrep_tuned_R100_report.md"))
cat(paste(report, collapse = "\n"))
