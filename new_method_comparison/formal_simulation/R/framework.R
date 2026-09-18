find_formal_root <- function() {
  candidates <- c(
    file.path(getwd(), "new_method_comparison", "formal_simulation"),
    file.path(getwd(), "formal_simulation"),
    getwd()
  )
  hit <- candidates[file.exists(file.path(candidates, "R", "framework.R"))]
  if (!length(hit)) stop("Cannot locate new_method_comparison/formal_simulation.")
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

.formal_root <- find_formal_root()
.comparison_root <- normalizePath(file.path(.formal_root, ".."),
                                  winslash = "/", mustWork = TRUE)
source(file.path(.comparison_root, "R", "methods", "wrappers.R"))

formal_default_config <- function() {
  list(
    n_train = 700L,
    n_test = 250L,
    p = 10L,
    theta_strength = 1,
    u_target_norm = 0.5,
    rho_x = 0,
    noise_sd = 0.35,
    basis_dimension = 5L,
    n_iter = 1200L,
    burn_in = 400L,
    n_replicates = 5L,
    seed_start = 2026061900L,
    grid_size = 50L,
    lambda0 = 10,
    lambda1 = 1,
    kappa2_alpha = 1,
    zeta0 = 0.1,
    zeta1 = 50
  )
}

validate_formal_config <- function(config) {
  required <- c(
    "n_train", "n_test", "p", "theta_strength", "u_target_norm",
    "rho_x", "noise_sd", "basis_dimension", "n_iter", "burn_in",
    "n_replicates", "seed_start"
  )
  missing <- setdiff(required, names(config))
  if (length(missing)) stop("Missing config fields: ", paste(missing, collapse = ", "))
  if (config$p < 10L) stop("The manuscript design requires p >= 10.")
  if (abs(config$rho_x) >= 1) stop("rho_x must lie strictly between -1 and 1.")
  if (config$burn_in >= config$n_iter) stop("burn_in must be smaller than n_iter.")
  invisible(config)
}

merge_formal_config <- function(overrides = list()) {
  config <- modifyList(formal_default_config(), overrides)
  validate_formal_config(config)
  config
}

make_ar1_predictors <- function(n, p, rho_x) {
  sigma_x <- outer(seq_len(p), seq_len(p),
                   function(i, j) rho_x^abs(i - j))
  z <- matrix(rnorm(n * p), n, p)
  z %*% chol(sigma_x)
}

generate_formal_dataset <- function(config, replicate_id = 1L) {
  validate_formal_config(config)
  data_seed <- as.integer(config$seed_start + replicate_id)
  set.seed(data_seed)
  n <- config$n_train + config$n_test
  train_id <- seq_len(config$n_train)
  test_id <- config$n_train + seq_len(config$n_test)
  coords <- cbind(runif(n), runif(n))
  X_raw <- make_ar1_predictors(n, config$p, config$rho_x)
  x_mean <- colMeans(X_raw[train_id, , drop = FALSE])
  x_sd <- apply(X_raw[train_id, , drop = FALSE], 2, stats::sd)
  X <- sweep(sweep(X_raw, 2, x_mean, "-"), 2, x_sd, "/")

  bump_raw <- function(s) {
    exp(-((s[, 1] - 0.35)^2 + (s[, 2] - 0.65)^2) / 0.04)
  }
  raw_train <- bump_raw(coords[train_id, , drop = FALSE])
  bump_train_mean <- mean(raw_train)
  bump_train_norm <- sqrt(mean((raw_train - bump_train_mean)^2))
  bump <- function(s) {
    config$u_target_norm *
      (bump_raw(s) - bump_train_mean) / bump_train_norm
  }

  theta <- numeric(config$p)
  theta[c(1, 2, 5, 6)] <- config$theta_strength * c(1, -1, 1, -1)
  spatial_sign <- numeric(config$p)
  spatial_sign[3:6] <- c(1, -1, 1, -1)
  truth_at <- function(s) {
    u <- outer(bump(s), spatial_sign)
    list(u = u, beta = sweep(u, 2, theta, "+"))
  }
  truth_all <- truth_at(coords)
  mu <- rowSums(X * truth_all$beta)
  noise_z <- rnorm(n)
  y <- mu + config$noise_sd * noise_z

  grid_coords <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = config$grid_size),
    y = seq(0, 1, length.out = config$grid_size)
  ))
  truth_grid <- truth_at(grid_coords)
  list(
    train = list(y = y[train_id], X = X[train_id, , drop = FALSE],
                 coords = coords[train_id, , drop = FALSE]),
    test = list(y = y[test_id], X = X[test_id, , drop = FALSE],
                coords = coords[test_id, , drop = FALSE]),
    grid = list(coords = grid_coords, true_beta = truth_grid$beta,
                true_u = truth_grid$u, grid_size = config$grid_size),
    truth = list(
      theta = theta,
      spatial_sign = spatial_sign,
      spatial_deviation = spatial_sign != 0,
      active_surface = seq_len(config$p) <= 6,
      bump_train_mean = bump_train_mean,
      bump_train_norm = bump_train_norm,
      u_target_norm = config$u_target_norm
    ),
    preprocessing = list(x_mean = x_mean, x_sd = x_sd),
    random_numbers = list(data_seed = data_seed, noise_z = noise_z),
    config = config
  )
}

formal_method_registry <- function() {
  list(
    proposed_ssgl = fit_proposed_ssgl,
    original_ssgl = fit_original_ssgl,
    global_only_blasso = fit_global_only_blasso,
    full_svc_no_selection = fit_full_svc_no_selection,
    standard_thin_plate_gam = fit_original_gam
  )
}

formal_fit_configs <- function(data, config) {
  list(
    basis = list(n_basis = config$basis_dimension),
    model = list(
      lambda0 = config$lambda0, lambda1 = config$lambda1,
      kappa2_alpha = config$kappa2_alpha,
      a_sigma = 0.5, b_sigma = stats::var(data$train$y) / 2,
      a_theta = 1, b_theta = 1, a_gamma = 1, b_gamma = config$p,
      zeta0 = config$zeta0, zeta1 = config$zeta1
    ),
    mcmc = list(n_iter = config$n_iter, burn_in = config$burn_in)
  )
}

gamma_switches <- function(gamma_draws) {
  if (is.null(gamma_draws)) return(NULL)
  apply(gamma_draws, 1, function(x) sum(diff(x) != 0))
}

summarize_formal_fit <- function(fit, data, replicate_id, method_name,
                                 theta_strength, u_target_norm) {
  beta_mse_j <- colMeans((fit$beta_hat_grid - data$grid$true_beta)^2)
  u_mse_j <- if (is.null(fit$u_hat_grid)) rep(NA_real_, ncol(data$train$X)) else {
    colMeans((fit$u_hat_grid - data$grid$true_u)^2)
  }
  theta_bias_j <- if (is.null(fit$theta_mean)) rep(NA_real_, ncol(data$train$X)) else {
    fit$theta_mean - data$truth$theta
  }
  rb_pip <- fit$diagnostics$rb_pip %||% rep(NA_real_, ncol(data$train$X))
  sampled_pip <- fit$pip %||% rep(NA_real_, ncol(data$train$X))
  switches <- gamma_switches(fit$diagnostics$gamma_draws)
  if (is.null(switches)) switches <- rep(NA_real_, ncol(data$train$X))
  predictor <- data.frame(
    method = method_name, replicate = replicate_id,
    theta_strength = theta_strength, u_target_norm = u_target_norm,
    predictor_index = seq_len(ncol(data$train$X)),
    predictor = paste0("X", seq_len(ncol(data$train$X))),
    group = c(rep("global_only", 2), rep("spatial_only", 2),
              rep("global_plus_spatial", 2),
              rep("null", ncol(data$train$X) - 6)),
    true_theta = data$truth$theta,
    theta_estimate = fit$theta_mean %||% rep(NA_real_, ncol(data$train$X)),
    theta_bias = theta_bias_j,
    theta_sq_error = theta_bias_j^2,
    sampled_pip = sampled_pip,
    rb_pip = rb_pip,
    gamma_switches = switches,
    beta_mse = beta_mse_j,
    u_mse = u_mse_j,
    stringsAsFactors = FALSE
  )
  run <- data.frame(
    method = method_name, replicate = replicate_id,
    theta_strength = theta_strength, u_target_norm = u_target_norm,
    mspe = mean((data$test$y - fit$pred_test)^2),
    theta_bias = mean(theta_bias_j, na.rm = TRUE),
    theta_mse = mean(theta_bias_j^2, na.rm = TRUE),
    beta_mse = mean(beta_mse_j),
    beta_mse_active = mean(beta_mse_j[1:6]),
    beta_mse_null = mean(beta_mse_j[7:10]),
    u_mse = mean(u_mse_j, na.rm = TRUE),
    runtime_seconds = fit$runtime,
    stringsAsFactors = FALSE
  )
  list(run = run, predictor = predictor)
}

safe_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  saveRDS(object, tmp)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not atomically save ", path)
  invisible(path)
}

append_failure <- function(path, row) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  old <- if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else NULL
  write.csv(rbind(old, row), path, row.names = FALSE)
}

run_formal_replicate <- function(config, replicate_id, methods = NULL,
                                 output_dir, resume = TRUE) {
  if (is.null(methods)) methods <- names(formal_method_registry())
  registry <- formal_method_registry()
  unknown <- setdiff(methods, names(registry))
  if (length(unknown)) stop("Unknown methods: ", paste(unknown, collapse = ", "))
  data <- generate_formal_dataset(config, replicate_id)
  fit_configs <- formal_fit_configs(data, config)
  replicate_dir <- file.path(output_dir, sprintf("replicate_%03d", replicate_id))
  dir.create(replicate_dir, recursive = TRUE, showWarnings = FALSE)
  safe_save_rds(data, file.path(replicate_dir, "shared_dataset.rds"))
  rows <- list()
  for (i in seq_along(methods)) {
    method_name <- methods[i]
    result_path <- file.path(replicate_dir, paste0(method_name, ".rds"))
    if (resume && file.exists(result_path)) {
      saved <- readRDS(result_path)
      rows[[method_name]] <- saved$metrics
      next
    }
    sampler_seed <- as.integer(config$seed_start + 100000L + replicate_id)
    started <- Sys.time()
    result <- tryCatch({
      fit <- registry[[method_name]](
        data$train, data$test, data$grid,
        fit_configs$basis, fit_configs$model, fit_configs$mcmc,
        seed = sampler_seed
      )
      metrics <- summarize_formal_fit(
        fit, data, replicate_id, method_name,
        config$theta_strength, config$u_target_norm
      )
      object <- list(fit = fit, metrics = metrics, config = config,
                     data_seed = data$random_numbers$data_seed,
                     sampler_seed = sampler_seed)
      safe_save_rds(object, result_path)
      rows[[method_name]] <- metrics
      data.frame(method = method_name, replicate = replicate_id,
                 status = "success", started = as.character(started),
                 finished = as.character(Sys.time()),
                 runtime_seconds = fit$runtime)
    }, error = function(e) {
      failure <- data.frame(
        timestamp = as.character(Sys.time()), method = method_name,
        replicate = replicate_id, theta_strength = config$theta_strength,
        u_target_norm = config$u_target_norm,
        message = conditionMessage(e), stringsAsFactors = FALSE
      )
      append_failure(file.path(output_dir, "failure_log.csv"), failure)
      data.frame(method = method_name, replicate = replicate_id,
                 status = "failure", started = as.character(started),
                 finished = as.character(Sys.time()),
                 runtime_seconds = NA_real_)
    })
    write.table(result, file.path(output_dir, "runtime_log.csv"),
                sep = ",", row.names = FALSE,
                col.names = !file.exists(file.path(output_dir, "runtime_log.csv")),
                append = file.exists(file.path(output_dir, "runtime_log.csv")))
  }
  invisible(rows)
}

aggregate_saved_results <- function(output_dir) {
  files <- list.files(output_dir, pattern = "\\.rds$", recursive = TRUE,
                      full.names = TRUE)
  files <- files[basename(files) != "shared_dataset.rds"]
  objects <- lapply(files, readRDS)
  objects <- objects[vapply(objects, function(x) !is.null(x$metrics), logical(1))]
  if (!length(objects)) stop("No completed method result files found.")
  run <- do.call(rbind, lapply(objects, function(x) x$metrics$run))
  predictor <- do.call(rbind, lapply(objects, function(x) x$metrics$predictor))
  write.csv(run, file.path(output_dir, "all_run_metrics.csv"), row.names = FALSE)
  write.csv(predictor, file.path(output_dir, "all_predictor_metrics.csv"),
            row.names = FALSE)
  list(run = run, predictor = predictor)
}

plot_surface_panel <- function(beta, grid_size, path, title) {
  png(path, width = 2200, height = 900, res = 160)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(2, 5), mar = c(2.2, 2.2, 2.8, 2.2), oma = c(1, 1, 3, 1))
  axis_values <- seq(0, 1, length.out = grid_size)
  for (j in 1:10) {
    z <- matrix(beta[, j], nrow = grid_size, ncol = grid_size)
    image(axis_values, axis_values, z, col = hcl.colors(80, "Blue-Red 3"),
          axes = FALSE, xlab = "", ylab = "", main = paste0("X", j))
    box()
  }
  mtext(title, outer = TRUE, cex = 1.3, font = 2)
}

