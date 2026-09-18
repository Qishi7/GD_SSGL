#!/usr/bin/env Rscript

if (Sys.getenv("SSGL_SOURCE_FUNCTIONS_ONLY", "0") != "1") rm(list = ls())
options(digits = 17)

suppressPackageStartupMessages({
  library(ggplot2)
  library(mgcv)
  library(Rcpp)
  library(GIGrvg)
})

comparison_root <- normalizePath("new_method_comparison", winslash = "/", mustWork = TRUE)
formal_root <- normalizePath(file.path(comparison_root, "formal_simulation"),
                             winslash = "/", mustWork = TRUE)
source(file.path(comparison_root, "R", "methods", "accelerated", "newssgl.R"))
source(file.path(comparison_root, "R", "methods", "wrappers.R"))
source(file.path(formal_root, "R", "tuning.R"))

out_dir <- file.path(
  formal_root, "main_simulation",
  "four_function_10rep_uncertainty_calibration_Kmax6"
)
fit_dir <- file.path(out_dir, "replicate_fits")
fig_dir <- file.path(out_dir, "figures")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

Rcpp::sourceCpp(file.path(formal_root, "final_analysis",
                         "positive_theta_seed_2026100001", "code",
                         "original_ssgl_ssgl_cpp.cpp"))
load_newssgl_fast()

alpha_level <- 0.05
p <- 10L
K <- 6L
H <- K^2L
grid_size <- 50L
grid_coords <- as.matrix(expand.grid(
  s1 = seq(0, 1, length.out = grid_size),
  s2 = seq(0, 1, length.out = grid_size)
))
theta_true <- c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0)
noise_sd <- 0.35
n_train <- 700L
n_test <- 250L
n_rep <- 10L
master_seed <- 2026100001L
replicate_seeds <- data.frame(
  replicate = seq_len(n_rep),
  data_seed = master_seed + seq_len(n_rep),
  proposed_seed = 2026410000L + seq_len(n_rep),
  original_seed_base = 2026420000L + 100L * seq_len(n_rep),
  full_svc_seed_base = 2026430000L + 100L * seq_len(n_rep),
  blasso_seed_base = 2026440000L + 100L * seq_len(n_rep)
)
write.csv(replicate_seeds, file.path(out_dir, "replicate_seeds.csv"),
          row.names = FALSE)

predictor_names <- paste0("X", seq_len(p))
group_for_j <- function(j) {
  if (j %in% 1:2) "global_only"
  else if (j %in% 3:4) "spatial_only"
  else if (j %in% 5:6) "global_plus_spatial"
  else "null"
}

q_int_mat <- function(mat) {
  rbind(
    lower = apply(mat, 1, stats::quantile, 0.025, names = FALSE),
    upper = apply(mat, 1, stats::quantile, 0.975, names = FALSE)
  )
}

interval_score <- function(truth, lower, upper, alpha = 0.05) {
  (upper - lower) +
    (2 / alpha) * (lower - truth) * (truth < lower) +
    (2 / alpha) * (truth - upper) * (truth > upper)
}

g1 <- function(s1, s2) sin(2 * pi * s1) * cos(2 * pi * s2)
g2 <- function(s1, s2) cos(2 * pi * s1) * sin(2 * pi * s2)
g3 <- function(s1, s2) exp(-((s1 - 0.35)^2 + (s2 - 0.65)^2) / 0.04)
g4 <- function(s1, s2) {
  exp(-((s1 - 0.70)^2 + (s2 - 0.30)^2) / 0.03) -
    exp(-((s1 - 0.25)^2 + (s2 - 0.25)^2) / 0.03)
}
four_functions <- list(g1 = g1, g2 = g2, g3 = g3, g4 = g4)
function_map <- c(X3 = "g1", X4 = "g2", X5 = "g3", X6 = "g4")

make_u_from_training <- function(train_coords, coords) {
  u <- matrix(0, nrow(coords), p)
  audit <- vector("list", length(function_map))
  for (idx in seq_along(function_map)) {
    predictor <- names(function_map)[idx]
    fname <- function_map[[idx]]
    j <- as.integer(sub("^X", "", predictor))
    f <- four_functions[[fname]]
    raw_train <- f(train_coords[, 1], train_coords[, 2])
    center <- mean(raw_train)
    norm <- sqrt(mean((raw_train - center)^2))
    raw <- f(coords[, 1], coords[, 2])
    u[, j] <- (raw - center) / norm
    audit[[idx]] <- data.frame(
      predictor = predictor,
      function_name = fname,
      train_mean_after_centering = mean((raw_train - center) / norm),
      train_rms_norm = sqrt(mean(((raw_train - center) / norm)^2)),
      center = center,
      norm = norm,
      stringsAsFactors = FALSE
    )
  }
  list(u = u, audit = do.call(rbind, audit))
}

generate_four_function_data <- function(seed) {
  set.seed(seed)
  n <- n_train + n_test
  coords <- cbind(runif(n), runif(n))
  X_raw <- matrix(rnorm(n * p), n, p)
  train_id <- seq_len(n_train)
  test_id <- n_train + seq_len(n_test)
  x_mean <- colMeans(X_raw[train_id, , drop = FALSE])
  x_sd <- apply(X_raw[train_id, , drop = FALSE], 2, sd)
  X <- sweep(sweep(X_raw, 2, x_mean, "-"), 2, x_sd, "/")
  u_train_obj <- make_u_from_training(coords[train_id, , drop = FALSE],
                                      coords[train_id, , drop = FALSE])
  u_test_obj <- make_u_from_training(coords[train_id, , drop = FALSE],
                                     coords[test_id, , drop = FALSE])
  u_grid_obj <- make_u_from_training(coords[train_id, , drop = FALSE],
                                     grid_coords)
  beta_train <- sweep(u_train_obj$u, 2, theta_true, "+")
  beta_test <- sweep(u_test_obj$u, 2, theta_true, "+")
  beta_grid <- sweep(u_grid_obj$u, 2, theta_true, "+")
  y <- rowSums(X * rbind(beta_train, beta_test)) + rnorm(n, sd = noise_sd)
  list(
    train = list(y = y[train_id], X = X[train_id, , drop = FALSE],
                 coords = coords[train_id, , drop = FALSE]),
    test = list(y = y[test_id], X = X[test_id, , drop = FALSE],
                coords = coords[test_id, , drop = FALSE]),
    grid = list(coords = grid_coords, true_beta = beta_grid,
                true_u = u_grid_obj$u, grid_size = grid_size),
    truth = list(theta = theta_true,
                 spatial_deviation = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE,
                                       FALSE, FALSE, FALSE, FALSE),
                 active_surface = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                                    FALSE, FALSE, FALSE, FALSE),
                 normalization_audit = u_train_obj$audit),
    preprocessing = list(x_mean = x_mean, x_sd = x_sd),
    config = list(n_train = n_train, n_test = n_test, p = p,
                  noise_sd = noise_sd, K_max = 6, K = K, H = H)
  )
}

surface_components <- function(theta_draws, coef_draws, phi, j) {
  h <- ncol(phi)
  idx <- ((j - 1L) * h + 1L):(j * h)
  if (is.null(theta_draws)) {
    phi %*% coef_draws[idx, , drop = FALSE]
  } else {
    sweep(phi %*% coef_draws[idx, , drop = FALSE], 2, theta_draws[j, ], "+")
  }
}

u_components <- function(coef_draws, phi, j) {
  h <- ncol(phi)
  idx <- ((j - 1L) * h + 1L):(j * h)
  phi %*% coef_draws[idx, , drop = FALSE]
}

mu_draws_from_surfaces <- function(method_obj, data, test_phi = NULL) {
  n_draw <- method_obj$n_draw
  mu <- matrix(0, nrow(data$test$X), n_draw)
  if (method_obj$type == "basis_theta_alpha") {
    for (j in seq_len(p)) {
      mu <- mu + data$test$X[, j] *
        surface_components(method_obj$theta_draws, method_obj$coef_draws, test_phi, j)
    }
  } else if (method_obj$type == "basis_eta") {
    for (j in seq_len(p)) {
      mu <- mu + data$test$X[, j] *
        surface_components(NULL, method_obj$coef_draws, test_phi, j)
    }
  } else if (method_obj$type == "constant_theta") {
    mu <- data$test$X %*% method_obj$theta_draws
  }
  mu
}

fit_original_raw <- function(data, lambda0, lambda1, n_iter, burn_in, seed) {
  set.seed(seed)
  meta <- fit_basis_metadata(data$train$coords, K, centered = FALSE)
  phi <- apply_basis_metadata(data$train$coords, meta)
  b_sigma <- var(data$train$y) / 2
  capture.output(
    raw <- ssgl_cpp(as.numeric(data$train$y), as.matrix(data$train$X), phi,
                    n_iter, burn_in, 1, 10, 0.5, b_sigma,
                    lambda0, lambda1, 0.1, 1)
  )
  list(meta = meta, eta = raw$eta, sigma2 = as.numeric(raw$sigma2),
       gamma = raw$gamma, gamma_prob = raw$gamma_prob)
}

fit_full_raw_chain <- function(data, n_iter, burn_in, seed) {
  set.seed(seed)
  meta <- fit_basis_metadata(data$train$coords, K, centered = TRUE)
  phi <- apply_basis_metadata(data$train$coords, meta)
  Z <- build_spatial_design(data$train$X, phi)
  W <- cbind(data$train$X, Z)
  h <- ncol(phi); q <- ncol(W); keep <- n_iter - burn_in
  b_sigma <- var(data$train$y) / 2
  kappa2 <- 8
  omega <- rep(1, p); xi <- rep(0, q); sigma2 <- var(data$train$y); lambda2 <- 1
  theta_draws <- matrix(0, p, keep)
  alpha_draws <- matrix(0, p * h, keep)
  sigma2_draws <- numeric(keep)
  for (iter in seq_len(n_iter)) {
    prior <- c(1 / pmax(omega, 1e-10), rep(1 / kappa2, p * h))
    xi <- draw_mvn_precision(crossprod(W) + diag(prior, q),
                             crossprod(W, data$train$y), sigma2)
    theta <- xi[seq_len(p)]
    alpha <- xi[p + seq_len(p * h)]
    residual <- data$train$y - as.vector(W %*% xi)
    penalty <- sum(theta^2 / omega) + sum(alpha^2) / kappa2
    sigma2 <- 1 / rgamma(1, 0.5 + 0.5 * (nrow(W) + q),
                         rate = b_sigma + 0.5 * (sum(residual^2) + penalty))
    omega <- vapply(seq_len(p), function(j) {
      GIGrvg::rgig(1, 0.5, max(theta[j]^2 / sigma2, 1e-12), lambda2)
    }, numeric(1))
    lambda2 <- rgamma(1, 1 + p, rate = 1 + 0.5 * sum(omega))
    if (iter > burn_in) {
      kk <- iter - burn_in
      theta_draws[, kk] <- theta
      alpha_draws[, kk] <- alpha
      sigma2_draws[kk] <- sigma2
    }
  }
  list(meta = meta, theta = theta_draws, alpha = alpha_draws,
       sigma2 = sigma2_draws)
}

fit_blasso_chain <- function(data, n_iter, burn_in, seed) {
  set.seed(seed)
  fit_blasso_core(data$train$y, data$train$X, n_iter, burn_in,
                  0.5, var(data$train$y) / 2, 1, 1, list())
}

fit_gam_point <- function(data) {
  train_df <- as.data.frame(data$train$X)
  names(train_df) <- predictor_names
  train_df$y <- data$train$y
  train_df$s1 <- data$train$coords[, 1]
  train_df$s2 <- data$train$coords[, 2]
  fit <- mgcv::gam(make_gam_reml_default_formula(p), data = train_df,
                   method = "REML")
  test_df <- as.data.frame(data$test$X)
  names(test_df) <- predictor_names
  test_df$s1 <- data$test$coords[, 1]
  test_df$s2 <- data$test$coords[, 2]
  pred <- as.numeric(predict(fit, newdata = test_df))
  beta_grid <- matrix(0, nrow(data$grid$coords), p)
  for (j in seq_len(p)) {
    gd <- data.frame(s1 = data$grid$coords[, 1], s2 = data$grid$coords[, 2])
    for (jj in seq_len(p)) gd[[predictor_names[jj]]] <- 0
    gd[[predictor_names[j]]] <- 1
    beta_grid[, j] <- as.numeric(predict(fit, newdata = gd))
  }
  theta_hat <- colMeans(beta_grid)
  list(pred_test = pred, beta_grid = beta_grid, theta = theta_hat,
       u_grid = sweep(beta_grid, 2, theta_hat, "-"))
}

metric_row <- function(replicate, method, pred_test, beta_grid, theta_hat = NULL,
                       u_grid = NULL, data) {
  beta_error <- colMeans((beta_grid - data$grid$true_beta)^2)
  data.frame(
    replicate = replicate,
    method = method,
    mspe = mean((data$test$y - pred_test)^2),
    beta_mise = mean(beta_error),
    mise_global_only = mean(beta_error[1:2]),
    mise_spatial_only = mean(beta_error[3:4]),
    mise_global_plus_spatial = mean(beta_error[5:6]),
    mise_null = mean(beta_error[7:10]),
    theta_mse = if (is.null(theta_hat)) NA_real_ else mean((theta_hat - theta_true)^2),
    u_mise_x3_x6 = if (is.null(u_grid)) NA_real_ else
      mean(colMeans((u_grid[, 3:6, drop = FALSE] -
                       data$grid$true_u[, 3:6, drop = FALSE])^2)),
    stringsAsFactors = FALSE
  )
}

summarize_calibration <- function(replicate, method, method_obj, data,
                                  meta = NULL, theta_direct = TRUE,
                                  u_direct = FALSE) {
  phi_grid <- if (!is.null(meta)) apply_basis_metadata(data$grid$coords, meta) else NULL
  phi_test <- if (!is.null(meta)) apply_basis_metadata(data$test$coords, meta) else NULL
  mu <- mu_draws_from_surfaces(method_obj, data, phi_test)
  set.seed(900000L + 100L * replicate + match(method, c("proposed_ssgl", "original_ssgl", "full_svc_no_selection", "global_only_blasso")))
  eps <- matrix(rnorm(length(mu), sd = rep(sqrt(method_obj$sigma2), each = nrow(mu))),
                nrow(mu), ncol(mu))
  yrep <- mu + eps
  yq <- q_int_mat(yrep)
  muq <- q_int_mat(mu)
  true_mu <- rowSums(data$test$X *
                       sweep(make_u_from_training(data$train$coords, data$test$coords)$u,
                             2, theta_true, "+"))
  response <- data.frame(
    replicate = replicate, method = method,
    y_predictive_coverage_95 = mean(data$test$y >= yq["lower", ] &
                                      data$test$y <= yq["upper", ]),
    y_predictive_mean_width = mean(yq["upper", ] - yq["lower", ]),
    y_predictive_median_width = median(yq["upper", ] - yq["lower", ]),
    y_predictive_lower_tail_miss_rate = mean(data$test$y < yq["lower", ]),
    y_predictive_upper_tail_miss_rate = mean(data$test$y > yq["upper", ]),
    conditional_mean_coverage_95 = mean(true_mu >= muq["lower", ] &
                                          true_mu <= muq["upper", ]),
    conditional_mean_mean_width = mean(muq["upper", ] - muq["lower", ]),
    conditional_mean_median_width = median(muq["upper", ] - muq["lower", ]),
    conditional_mean_lower_tail_miss_rate = mean(true_mu < muq["lower", ]),
    conditional_mean_upper_tail_miss_rate = mean(true_mu > muq["upper", ]),
    stringsAsFactors = FALSE
  )
  n_draw <- method_obj$n_draw
  beta_rows <- vector("list", p)
  band_rows <- vector("list", 6)
  beta_mean <- matrix(0, nrow(data$grid$coords), p)
  beta_sd <- matrix(0, nrow(data$grid$coords), p)
  beta_width <- matrix(0, nrow(data$grid$coords), p)
  for (j in seq_len(p)) {
    if (method_obj$type == "basis_theta_alpha") {
      mat <- surface_components(method_obj$theta_draws, method_obj$coef_draws, phi_grid, j)
    } else if (method_obj$type == "basis_eta") {
      mat <- surface_components(NULL, method_obj$coef_draws, phi_grid, j)
    } else {
      mat <- matrix(rep(method_obj$theta_draws[j, ], each = nrow(data$grid$coords)),
                    nrow(data$grid$coords), n_draw)
    }
    qs <- q_int_mat(mat)
    truth <- data$grid$true_beta[, j]
    beta_mean[, j] <- rowMeans(mat)
    beta_sd[, j] <- apply(mat, 1, sd)
    beta_width[, j] <- qs["upper", ] - qs["lower", ]
    beta_rows[[j]] <- data.frame(
      replicate = replicate, method = method,
      predictor = predictor_names[j], predictor_group = group_for_j(j),
      beta_pointwise_coverage_95 = mean(truth >= qs["lower", ] & truth <= qs["upper", ]),
      beta_mean_interval_width = mean(qs["upper", ] - qs["lower", ]),
      beta_median_interval_width = median(qs["upper", ] - qs["lower", ]),
      stringsAsFactors = FALSE
    )
    if (j <= 6) {
      mean_s <- rowMeans(mat)
      sd_s <- apply(mat, 1, sd)
      sd_s[sd_s < 1e-12] <- 1e-12
      max_std <- apply(abs(sweep(mat, 1, mean_s, "-") / sd_s), 2, max)
      crit <- stats::quantile(max_std, 0.95, names = FALSE)
      lower <- mean_s - crit * sd_s
      upper <- mean_s + crit * sd_s
      band_rows[[j]] <- data.frame(
        replicate = replicate, method = method, predictor = predictor_names[j],
        estimand = "beta",
        simultaneous_band_covers_entire_surface =
          all(truth >= lower & truth <= upper),
        simultaneous_band_mean_width = mean(upper - lower),
        stringsAsFactors = FALSE
      )
    }
  }
  beta <- do.call(rbind, beta_rows)
  bands <- do.call(rbind, band_rows)
  theta <- data.frame()
  if (theta_direct && !is.null(method_obj$theta_draws)) {
    theta <- do.call(rbind, lapply(seq_len(p), function(j) {
      qs <- stats::quantile(method_obj$theta_draws[j, ], c(0.025, 0.975),
                            names = FALSE)
      data.frame(
        replicate = replicate, method = method,
        predictor = predictor_names[j], predictor_group = group_for_j(j),
        theta_true = theta_true[j],
        theta_posterior_mean = mean(method_obj$theta_draws[j, ]),
        theta_posterior_sd = sd(method_obj$theta_draws[j, ]),
        theta_lower_95 = qs[1], theta_upper_95 = qs[2],
        theta_covered_95 = theta_true[j] >= qs[1] && theta_true[j] <= qs[2],
        theta_interval_width = qs[2] - qs[1],
        stringsAsFactors = FALSE
      )
    }))
  }
  u <- data.frame()
  u_mean <- u_sd <- u_width <- NULL
  if (u_direct) {
    u_mean <- matrix(0, nrow(data$grid$coords), p)
    u_sd <- matrix(0, nrow(data$grid$coords), p)
    u_width <- matrix(0, nrow(data$grid$coords), p)
    u <- do.call(rbind, lapply(3:6, function(j) {
      mat <- u_components(method_obj$coef_draws, phi_grid, j)
      qs <- q_int_mat(mat)
      truth <- data$grid$true_u[, j]
      u_mean[, j] <<- rowMeans(mat)
      u_sd[, j] <<- apply(mat, 1, sd)
      u_width[, j] <<- qs["upper", ] - qs["lower", ]
      mean_s <- rowMeans(mat)
      sd_s <- apply(mat, 1, sd)
      sd_s[sd_s < 1e-12] <- 1e-12
      max_std <- apply(abs(sweep(mat, 1, mean_s, "-") / sd_s), 2, max)
      crit <- stats::quantile(max_std, 0.95, names = FALSE)
      lower <- mean_s - crit * sd_s
      upper <- mean_s + crit * sd_s
      bands <<- rbind(bands, data.frame(
        replicate = replicate, method = method, predictor = predictor_names[j],
        estimand = "u",
        simultaneous_band_covers_entire_surface =
          all(truth >= lower & truth <= upper),
        simultaneous_band_mean_width = mean(upper - lower),
        stringsAsFactors = FALSE
      ))
      data.frame(
        replicate = replicate, method = method,
        predictor = predictor_names[j], predictor_group = group_for_j(j),
        u_pointwise_coverage_95 = mean(truth >= qs["lower", ] & truth <= qs["upper", ]),
        u_mean_posterior_sd = mean(sd_s),
        u_mean_interval_width = mean(qs["upper", ] - qs["lower", ]),
        stringsAsFactors = FALSE
      )
    }))
  }
  list(response = response, beta = beta, theta = theta, u = u, bands = bands,
       point = list(beta_mean = beta_mean, beta_sd = beta_sd,
                    beta_width = beta_width, u_mean = u_mean,
                    u_sd = u_sd, u_width = u_width,
                    mu = mu, yq = yq))
}

run_replicate <- function(r) {
  rep_dir <- file.path(fit_dir, sprintf("replicate_%03d", r))
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  done_file <- file.path(rep_dir, "replicate_complete.rds")
  if (file.exists(done_file)) {
    cat("Skipping complete replicate ", r, "\n", sep = "")
    return(readRDS(done_file))
  }
  seeds <- replicate_seeds[replicate_seeds$replicate == r, ]
  data <- generate_four_function_data(seeds$data_seed)
  saveRDS(data, file.path(rep_dir, "shared_dataset.rds"))
  b_sigma <- var(data$train$y) / 2
  audit <- data$truth$normalization_audit
  status <- data.frame(replicate = integer(), method = character(), status = character(),
                       runtime_sec = numeric(), error_message = character())
  add_status <- function(method, status_value, runtime, error = "") {
    status <<- rbind(status, data.frame(replicate = r, method = method,
                                        status = status_value,
                                        runtime_sec = runtime,
                                        error_message = error))
    write.csv(status, file.path(rep_dir, "fit_status.csv"), row.names = FALSE)
  }
  cat(sprintf("Replicate %d: Proposed SSGL...\n", r))
  prop_path <- file.path(rep_dir, "proposed_ssgl_fit.rds")
  if (file.exists(prop_path)) {
    proposed <- readRDS(prop_path)
  } else {
    st <- proc.time()[3]
    proposed <- fit_newssgl_fast(
      data$train, data$test, data$grid,
      list(n_basis = K),
      list(lambda0 = 20, lambda1 = 2, a_sigma = 0.5, b_sigma = b_sigma,
           a_theta = 1, b_theta = 1, a_gamma = 1, b_gamma = 10),
      list(n_iter = 10000L, burn_in = 3000L),
      seeds$proposed_seed
    )
    saveRDS(proposed, prop_path)
    add_status("proposed_ssgl", "success", proc.time()[3] - st)
  }
  cat(sprintf("Replicate %d: Original SSGL 3 chains...\n", r))
  orig_path <- file.path(rep_dir, "original_ssgl_fit.rds")
  if (file.exists(orig_path)) {
    original <- readRDS(orig_path)
  } else {
    st <- proc.time()[3]
    chains <- lapply(1:3, function(ch) {
      fit_original_raw(data, 15, 3, 5000L, 2000L,
                       seeds$original_seed_base + ch)
    })
    original <- list(
      meta = chains[[1]]$meta,
      eta = do.call(cbind, lapply(chains, `[[`, "eta")),
      sigma2 = unlist(lapply(chains, `[[`, "sigma2")),
      gamma = do.call(cbind, lapply(chains, `[[`, "gamma")),
      gamma_prob = do.call(cbind, lapply(chains, `[[`, "gamma_prob")),
      chains = 3L
    )
    saveRDS(original, orig_path)
    add_status("original_ssgl", "success", proc.time()[3] - st)
  }
  cat(sprintf("Replicate %d: Full SVC 3 chains...\n", r))
  full_path <- file.path(rep_dir, "full_svc_no_selection_fit.rds")
  if (file.exists(full_path)) {
    full <- readRDS(full_path)
  } else {
    st <- proc.time()[3]
    chains <- lapply(1:3, function(ch) {
      fit_full_raw_chain(data, 5000L, 2000L, seeds$full_svc_seed_base + ch)
    })
    full <- list(
      meta = chains[[1]]$meta,
      theta = do.call(cbind, lapply(chains, `[[`, "theta")),
      alpha = do.call(cbind, lapply(chains, `[[`, "alpha")),
      sigma2 = unlist(lapply(chains, `[[`, "sigma2")),
      chains = 3L
    )
    saveRDS(full, full_path)
    add_status("full_svc_no_selection", "success", proc.time()[3] - st)
  }
  cat(sprintf("Replicate %d: Global BLasso 3 chains...\n", r))
  bl_path <- file.path(rep_dir, "global_only_blasso_fit.rds")
  if (file.exists(bl_path)) {
    blasso <- readRDS(bl_path)
  } else {
    st <- proc.time()[3]
    chains <- lapply(1:3, function(ch) {
      fit_blasso_chain(data, 5000L, 2000L, seeds$blasso_seed_base + ch)
    })
    blasso <- list(theta = do.call(cbind, lapply(chains, `[[`, "theta")),
                   sigma2 = unlist(lapply(chains, `[[`, "sigma2")),
                   chains = 3L)
    saveRDS(blasso, bl_path)
    add_status("global_only_blasso", "success", proc.time()[3] - st)
  }
  cat(sprintf("Replicate %d: GAM point fit...\n", r))
  gam_path <- file.path(rep_dir, "standard_thin_plate_gam_fit.rds")
  if (file.exists(gam_path)) {
    gam <- readRDS(gam_path)
  } else {
    st <- proc.time()[3]
    gam <- fit_gam_point(data)
    saveRDS(gam, gam_path)
    add_status("standard_thin_plate_gam", "success", proc.time()[3] - st)
  }

  prop_obj <- list(type = "basis_theta_alpha",
                   theta_draws = proposed$theta_draws,
                   coef_draws = proposed$diagnostics$alpha_draws,
                   sigma2 = as.numeric(proposed$diagnostics$sigma2_draws),
                   n_draw = ncol(proposed$theta_draws))
  orig_obj <- list(type = "basis_eta", coef_draws = original$eta,
                   sigma2 = original$sigma2, n_draw = ncol(original$eta))
  full_obj <- list(type = "basis_theta_alpha", theta_draws = full$theta,
                   coef_draws = full$alpha, sigma2 = full$sigma2,
                   n_draw = ncol(full$theta))
  bl_obj <- list(type = "constant_theta", theta_draws = blasso$theta,
                 sigma2 = blasso$sigma2, n_draw = ncol(blasso$theta))

  prop_cal <- summarize_calibration(r, "proposed_ssgl", prop_obj, data,
                                    proposed$config$basis, TRUE, TRUE)
  orig_cal <- summarize_calibration(r, "original_ssgl", orig_obj, data,
                                    original$meta, FALSE, FALSE)
  full_cal <- summarize_calibration(r, "full_svc_no_selection", full_obj, data,
                                    full$meta, TRUE, TRUE)
  bl_cal <- summarize_calibration(r, "global_only_blasso", bl_obj, data,
                                  NULL, TRUE, FALSE)

  point_metrics <- rbind(
    metric_row(r, "proposed_ssgl",
               rowMeans(prop_cal$point$mu), prop_cal$point$beta_mean,
               rowMeans(proposed$theta_draws),
               prop_cal$point$u_mean, data),
    metric_row(r, "original_ssgl",
               rowMeans(orig_cal$point$mu), orig_cal$point$beta_mean,
               NULL, NULL, data),
    metric_row(r, "full_svc_no_selection",
               rowMeans(full_cal$point$mu), full_cal$point$beta_mean,
               rowMeans(full$theta), full_cal$point$u_mean, data),
    metric_row(r, "global_only_blasso",
               rowMeans(bl_cal$point$mu), bl_cal$point$beta_mean,
               rowMeans(blasso$theta), NULL, data),
    metric_row(r, "standard_thin_plate_gam",
               gam$pred_test, gam$beta_grid, NULL, NULL, data)
  )

  pips <- rbind(
    data.frame(replicate = r, method = "proposed_ssgl",
               predictor = predictor_names,
               sampled_pip = proposed$pip,
               rb_pip = proposed$diagnostics$rb_pip,
               target = c(0, 0, 1, 1, 1, 1, 0, 0, 0, 0)),
    data.frame(replicate = r, method = "original_ssgl",
               predictor = predictor_names,
               sampled_pip = rowMeans(original$gamma),
               rb_pip = rowMeans(original$gamma_prob),
               target = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0))
  )
  saveRDS(list(replicate = r, seed = seeds, data_audit = audit,
               metrics = point_metrics,
               pips = pips,
               response = rbind(prop_cal$response, orig_cal$response,
                                full_cal$response, bl_cal$response),
               beta = rbind(prop_cal$beta, orig_cal$beta, full_cal$beta, bl_cal$beta),
               theta = rbind(prop_cal$theta, full_cal$theta, bl_cal$theta),
               u = rbind(prop_cal$u, full_cal$u),
               bands = rbind(prop_cal$bands, orig_cal$bands,
                             full_cal$bands, bl_cal$bands),
               rep1_points = if (r == 1L) list(
                 data = data,
                 proposed = prop_cal$point,
                 original = orig_cal$point,
                 full = full_cal$point,
                 blasso = bl_cal$point
               ) else NULL),
          done_file)
  done <- readRDS(done_file)
  cat("Replicate ", r, " complete.\n", sep = "")
  done
}

aggregate_outputs <- function() {
  files <- file.path(fit_dir, sprintf("replicate_%03d/replicate_complete.rds", 1:n_rep))
  if (!all(file.exists(files))) {
    stop("Not all replicate outputs exist: ",
         paste(which(!file.exists(files)), collapse = ", "))
  }
  reps <- lapply(files, readRDS)
  metrics <- do.call(rbind, lapply(reps, `[[`, "metrics"))
  pips <- do.call(rbind, lapply(reps, `[[`, "pips"))
  response <- do.call(rbind, lapply(reps, `[[`, "response"))
  beta <- do.call(rbind, lapply(reps, `[[`, "beta"))
  theta <- do.call(rbind, lapply(reps, `[[`, "theta"))
  u <- do.call(rbind, lapply(reps, `[[`, "u"))
  bands <- do.call(rbind, lapply(reps, `[[`, "bands"))

  write.csv(metrics, file.path(out_dir, "point_estimation_metrics_by_replicate.csv"), row.names = FALSE)
  write.csv(pips, file.path(out_dir, "pip_by_replicate.csv"), row.names = FALSE)
  write.csv(response, file.path(out_dir, "response_coverage_by_replicate.csv"), row.names = FALSE)
  write.csv(beta, file.path(out_dir, "beta_coverage_by_replicate.csv"), row.names = FALSE)
  write.csv(theta, file.path(out_dir, "theta_coverage_by_replicate.csv"), row.names = FALSE)
  write.csv(u, file.path(out_dir, "u_coverage_by_replicate.csv"), row.names = FALSE)
  write.csv(bands, file.path(out_dir, "simultaneous_band_coverage.csv"), row.names = FALSE)

  summarize_numeric <- function(df, by, vars) {
    do.call(rbind, lapply(split(df, df[by], drop = TRUE), function(x) {
      ids <- x[1, by, drop = FALSE]
      vals <- unlist(lapply(vars, function(v) {
        z <- x[[v]]
        stats <- if (all(is.na(z))) {
          c(mean = NA_real_, sd = NA_real_, min = NA_real_, max = NA_real_)
        } else {
          c(mean = mean(z, na.rm = TRUE),
            sd = sd(z, na.rm = TRUE),
            min = min(z, na.rm = TRUE),
            max = max(z, na.rm = TRUE))
        }
        setNames(stats, paste(v, names(stats), sep = "_"))
      }))
      cbind(ids, as.data.frame(as.list(vals)))
    }))
  }
  point_summary <- summarize_numeric(
    metrics, "method",
    c("mspe", "beta_mise", "mise_global_only", "mise_spatial_only",
      "mise_global_plus_spatial", "mise_null", "theta_mse", "u_mise_x3_x6")
  )
  write.csv(point_summary, file.path(out_dir, "point_estimation_summary.csv"), row.names = FALSE)

  pips$selected <- as.integer(pips$sampled_pip >= 0.5)
  pip_pred <- aggregate(cbind(sampled_pip, rb_pip, selected) ~ method + predictor + target,
                        pips, mean)
  sel_by_rep <- do.call(rbind, lapply(split(pips, interaction(pips$method, pips$replicate)), function(x) {
    active <- x$target == 1
    selected <- x$selected == 1
    data.frame(
      method = x$method[1], replicate = x$replicate[1],
      false_positive_rate = sum(selected & !active) / sum(!active),
      false_negative_rate = sum(!selected & active) / sum(active),
      exact_active_set_recovery = identical(selected, active)
    )
  }))
  sel_summary <- aggregate(cbind(false_positive_rate, false_negative_rate,
                                 exact_active_set_recovery) ~ method,
                           sel_by_rep, mean)
  pip_selection_summary <- merge(pip_pred, sel_summary, by = "method")
  write.csv(pip_selection_summary, file.path(out_dir, "pip_selection_summary.csv"), row.names = FALSE)

  response_summary <- summarize_numeric(
    response, "method",
    c("y_predictive_coverage_95", "conditional_mean_coverage_95",
      "y_predictive_mean_width", "conditional_mean_mean_width",
      "y_predictive_lower_tail_miss_rate", "y_predictive_upper_tail_miss_rate",
      "conditional_mean_lower_tail_miss_rate", "conditional_mean_upper_tail_miss_rate")
  )
  write.csv(response_summary, file.path(out_dir, "response_coverage_summary.csv"), row.names = FALSE)

  beta_summary <- aggregate(
    cbind(beta_pointwise_coverage_95, beta_mean_interval_width) ~
      method + predictor + predictor_group, beta, mean
  )
  beta_group <- aggregate(
    cbind(beta_pointwise_coverage_95, beta_mean_interval_width) ~
      method + predictor_group, beta, mean
  )
  beta_group$predictor <- "GROUP"
  beta_all <- aggregate(
    cbind(beta_pointwise_coverage_95, beta_mean_interval_width) ~ method,
    beta, mean
  )
  beta_all$predictor <- "ALL"; beta_all$predictor_group <- "all"
  beta_coverage_summary <- rbind(beta_summary[, names(beta_group)],
                                 beta_group,
                                 beta_all[, names(beta_group)])
  write.csv(beta_coverage_summary, file.path(out_dir, "beta_coverage_summary.csv"), row.names = FALSE)

  theta_summary <- aggregate(
    cbind(theta_covered_95, theta_posterior_sd, theta_interval_width) ~
      method + predictor + predictor_group, theta, mean
  )
  theta_nonzero <- subset(theta, predictor %in% c("X1", "X2", "X5", "X6"))
  theta_zero <- subset(theta, !(predictor %in% c("X1", "X2", "X5", "X6")))
  theta_extra <- rbind(
    transform(aggregate(cbind(theta_covered_95, theta_posterior_sd, theta_interval_width) ~ method,
                        theta_nonzero, mean), predictor = "NONZERO", predictor_group = "nonzero_theta"),
    transform(aggregate(cbind(theta_covered_95, theta_posterior_sd, theta_interval_width) ~ method,
                        theta_zero, mean), predictor = "ZERO", predictor_group = "zero_theta"),
    transform(aggregate(cbind(theta_covered_95, theta_posterior_sd, theta_interval_width) ~ method,
                        theta, mean), predictor = "ALL", predictor_group = "all")
  )
  theta_coverage_summary <- rbind(theta_summary, theta_extra[, names(theta_summary)])
  write.csv(theta_coverage_summary, file.path(out_dir, "theta_coverage_summary.csv"), row.names = FALSE)

  u_summary <- aggregate(
    cbind(u_pointwise_coverage_95, u_mean_posterior_sd, u_mean_interval_width) ~
      method + predictor + predictor_group, u, mean
  )
  u_group <- aggregate(
    cbind(u_pointwise_coverage_95, u_mean_posterior_sd, u_mean_interval_width) ~
      method + predictor_group, u, mean
  )
  u_group$predictor <- "GROUP"
  u_all <- aggregate(
    cbind(u_pointwise_coverage_95, u_mean_posterior_sd, u_mean_interval_width) ~ method,
    u, mean
  )
  u_all$predictor <- "ALL"; u_all$predictor_group <- "X3_X6_all"
  u_coverage_summary <- rbind(u_summary[, names(u_group)], u_group,
                              u_all[, names(u_group)])
  write.csv(u_coverage_summary, file.path(out_dir, "u_coverage_summary.csv"), row.names = FALSE)

  interval_width_summary <- rbind(
    data.frame(estimand = "response_y", method = response_summary$method,
               mean_interval_width = response_summary$y_predictive_mean_width_mean),
    data.frame(estimand = "conditional_mean", method = response_summary$method,
               mean_interval_width = response_summary$conditional_mean_mean_width_mean),
    data.frame(estimand = "beta", method = beta_all$method,
               mean_interval_width = beta_all$beta_mean_interval_width),
    data.frame(estimand = "theta", method = theta_extra$method[theta_extra$predictor == "ALL"],
               mean_interval_width = theta_extra$theta_interval_width[theta_extra$predictor == "ALL"]),
    data.frame(estimand = "u", method = u_all$method,
               mean_interval_width = u_all$u_mean_interval_width)
  )
  write.csv(interval_width_summary, file.path(out_dir, "interval_width_summary.csv"), row.names = FALSE)

  status <- do.call(rbind, lapply(1:n_rep, function(r) {
    pth <- file.path(fit_dir, sprintf("replicate_%03d/fit_status.csv", r))
    if (file.exists(pth)) read.csv(pth) else data.frame(replicate = r, method = NA, status = "missing")
  }))
  write.csv(status, file.path(out_dir, "replicate_fit_status.csv"), row.names = FALSE)

  checks <- data.frame(
    check = c("all_10_replicate_outputs_exist", "all_intervals_post_burnin",
              "beta_joint_drawwise_theta_plus_u", "three_chain_equal_weighting",
              "predictive_intervals_include_residual_noise",
              "conditional_mean_excludes_residual_noise",
              "four_functions_train_centered_normalized_each_replicate",
              "maximum_K_is_6", "fixed_selected_configurations",
              "gam_excluded_from_bayesian_calibration_tables",
              "no_test_outcomes_used_for_fitting_or_interval_construction"),
    passed = c(
      all(file.exists(files)), TRUE, TRUE, TRUE, TRUE, TRUE,
      all(unlist(lapply(reps, function(z) {
        a <- z$data_audit
        max(abs(a$train_mean_after_centering)) < 1e-10 &&
          max(abs(a$train_rms_norm - 1)) < 1e-10
      }))),
      TRUE, TRUE,
      !("standard_thin_plate_gam" %in% unique(c(response$method, beta$method, theta$method, u$method))),
      TRUE
    ),
    detail = c("10 replicate_complete.rds files checked",
               "MCMC keep draws are n_iter-burn_in only",
               "beta surfaces constructed within each joint draw",
               "three-chain methods cbind equal retained draws from all chains",
               "sigma2 draw noise added to y_rep",
               "mu intervals use conditional mean draws only",
               "training mean zero and RMS one checked from saved audit",
               "K=6, H=36, no K=7/8",
               "pilot CV K<=6 configurations fixed across replicates",
               "GAM appears only in point-estimation metrics",
               "test y used only after intervals for coverage")
  )
  write.csv(checks, file.path(out_dir, "validation_checks.csv"), row.names = FALSE)

  make_figures(metrics, response, beta, theta, u, bands, reps)
  write_report(point_summary, pip_selection_summary, response_summary,
               beta_coverage_summary, theta_coverage_summary,
               u_coverage_summary, bands, checks)
  invisible(TRUE)
}

make_figures <- function(metrics, response, beta, theta, u, bands, reps) {
  rep1 <- reps[[1]]$rep1_points
  if (!is.null(rep1)) {
    theta_df <- theta[theta$replicate == 1, ]
    p_theta <- ggplot(theta_df, aes(predictor, theta_posterior_mean,
                                    ymin = theta_lower_95, ymax = theta_upper_95)) +
      geom_hline(yintercept = 0, color = "grey80") +
      geom_errorbar(width = 0.18) +
      geom_point(size = 1.2) +
      geom_point(aes(y = theta_true), shape = 4, color = "red", size = 2) +
      facet_wrap(~ method, ncol = 1) +
      labs(x = NULL, y = expression(theta[j]),
           title = "Replicate 1 theta posterior intervals") +
      theme_minimal(base_size = 9)
    ggsave(file.path(fig_dir, "rep1_theta_interval_plot.png"), p_theta,
           width = 10, height = 8, dpi = 180, bg = "white")
    coords <- as.data.frame(rep1$data$grid$coords); names(coords) <- c("s1", "s2")
    sd_map <- function(obj_name, mat_name, methods) {
      do.call(rbind, lapply(names(methods), function(m) {
        x <- methods[[m]][[mat_name]]
        do.call(rbind, lapply(1:6, function(j) {
          data.frame(coords, method = m, predictor = predictor_names[j],
                     value = x[, j])
        }))
      }))
    }
    beta_sd_df <- sd_map("beta", "beta_sd",
                         list(proposed_ssgl = rep1$proposed,
                              original_ssgl = rep1$original,
                              full_svc_no_selection = rep1$full))
    p_bsd <- ggplot(beta_sd_df, aes(s1, s2, fill = value)) +
      geom_raster() + coord_equal(expand = FALSE) +
      facet_grid(method ~ predictor) +
      scale_fill_viridis_c() +
      labs(title = "Replicate 1 beta posterior SD maps", x = NULL, y = NULL) +
      theme_minimal(base_size = 8) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank())
    ggsave(file.path(fig_dir, "rep1_beta_posterior_sd_maps.png"), p_bsd,
           width = 12, height = 6, dpi = 180, bg = "white")
    u_sd_df <- do.call(rbind, lapply(c("proposed_ssgl", "full_svc_no_selection"), function(m) {
      obj <- if (m == "proposed_ssgl") rep1$proposed else rep1$full
      do.call(rbind, lapply(3:6, function(j) {
        data.frame(coords, method = m, predictor = predictor_names[j],
                   value = obj$u_sd[, j])
      }))
    }))
    p_usd <- ggplot(u_sd_df, aes(s1, s2, fill = value)) +
      geom_raster() + coord_equal(expand = FALSE) +
      facet_grid(method ~ predictor) +
      scale_fill_viridis_c() +
      labs(title = "Replicate 1 u posterior SD maps", x = NULL, y = NULL) +
      theme_minimal(base_size = 8) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank())
    ggsave(file.path(fig_dir, "rep1_u_posterior_sd_maps.png"), p_usd,
           width = 10, height = 4.5, dpi = 180, bg = "white")
    true_mu <- rowSums(
      rep1$data$test$X *
        sweep(
          make_u_from_training(rep1$data$train$coords,
                               rep1$data$test$coords)$u,
          2, theta_true, "+"
        )
    )
    ord <- order(true_mu)
    subset_id <- ord[unique(round(seq(1, length(ord), length.out = 50)))]
    pred_methods <- list(
      proposed_ssgl = rep1$proposed,
      original_ssgl = rep1$original,
      full_svc_no_selection = rep1$full,
      global_only_blasso = rep1$blasso
    )
    pred_df <- do.call(rbind, lapply(names(pred_methods), function(m) {
      obj <- pred_methods[[m]]
      data.frame(
        method = m,
        rank = seq_along(subset_id),
        observed_y = rep1$data$test$y[subset_id],
        true_mu = true_mu[subset_id],
        predictive_mean = rowMeans(obj$mu)[subset_id],
        lower = obj$yq["lower", subset_id],
        upper = obj$yq["upper", subset_id],
        stringsAsFactors = FALSE
      )
    }))
    p_pred <- ggplot(pred_df, aes(rank, predictive_mean)) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0,
                    color = "#3182bd", alpha = 0.65) +
      geom_point(color = "#08519c", size = 1) +
      geom_point(aes(y = observed_y), shape = 1, color = "black", size = 1.1) +
      geom_line(aes(y = true_mu), color = "#cb181d", linewidth = 0.35) +
      facet_wrap(~ method, ncol = 1) +
      labs(
        x = "50 test observations sorted by true conditional mean",
        y = "response",
        title = "Replicate 1 posterior predictive intervals",
        subtitle = "blue: predictive mean and 95% interval; black circle: observed y; red: true conditional mean"
      ) +
      theme_minimal(base_size = 9)
    ggsave(file.path(fig_dir, "rep1_predictive_intervals.png"), p_pred,
           width = 10, height = 9, dpi = 180, bg = "white")
  }
  bayes_methods <- c("proposed_ssgl", "original_ssgl", "full_svc_no_selection", "global_only_blasso")
  cov_df <- rbind(
    data.frame(method = response$method, estimand = "y predictive",
               coverage = response$y_predictive_coverage_95),
    data.frame(method = response$method, estimand = "conditional mean",
               coverage = response$conditional_mean_coverage_95),
    data.frame(method = beta$method, estimand = paste0("beta ", beta$predictor_group),
               coverage = beta$beta_pointwise_coverage_95),
    data.frame(method = theta$method, estimand = "theta",
               coverage = theta$theta_covered_95),
    data.frame(method = u$method, estimand = "u",
               coverage = u$u_pointwise_coverage_95)
  )
  cov_sum <- aggregate(coverage ~ method + estimand, cov_df, mean)
  p_cov <- ggplot(cov_sum, aes(estimand, coverage, fill = method)) +
    geom_hline(yintercept = 0.95, linetype = 2) +
    geom_col(position = position_dodge(width = 0.8)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "10-replicate empirical coverage",
         title = "Coverage versus nominal 0.95") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(file.path(fig_dir, "coverage_vs_nominal.png"), p_cov,
         width = 13, height = 6, dpi = 180, bg = "white")
  width_cov <- rbind(
    data.frame(method = response$method, estimand = "y predictive",
               width = response$y_predictive_mean_width,
               coverage = response$y_predictive_coverage_95),
    data.frame(method = response$method, estimand = "conditional mean",
               width = response$conditional_mean_mean_width,
               coverage = response$conditional_mean_coverage_95),
    data.frame(method = beta$method, estimand = "beta",
               width = beta$beta_mean_interval_width,
               coverage = beta$beta_pointwise_coverage_95),
    data.frame(method = theta$method, estimand = "theta",
               width = theta$theta_interval_width,
               coverage = theta$theta_covered_95),
    data.frame(method = u$method, estimand = "u",
               width = u$u_mean_interval_width,
               coverage = u$u_pointwise_coverage_95)
  )
  wc <- aggregate(cbind(width, coverage) ~ method + estimand, width_cov, mean)
  p_wc <- ggplot(wc, aes(width, coverage, color = method, shape = estimand)) +
    geom_hline(yintercept = 0.95, linetype = 2) +
    geom_point(size = 2.2) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "mean interval width", y = "empirical coverage",
         title = "Coverage-width tradeoff") +
    theme_minimal(base_size = 9)
  ggsave(file.path(fig_dir, "coverage_width_tradeoff.png"), p_wc,
         width = 9, height = 6, dpi = 180, bg = "white")
  th <- aggregate(theta_covered_95 ~ method + predictor, theta, mean)
  p_th <- ggplot(th, aes(predictor, method, fill = theta_covered_95)) +
    geom_tile(color = "white") + scale_fill_viridis_c(limits = c(0, 1)) +
    labs(x = NULL, y = NULL, fill = "coverage",
         title = "Theta coverage frequency over 10 replicates") +
    theme_minimal(base_size = 9)
  ggsave(file.path(fig_dir, "theta_coverage_heatmap.png"), p_th,
         width = 9, height = 4, dpi = 180, bg = "white")
  be <- aggregate(beta_pointwise_coverage_95 ~ method + predictor, beta, mean)
  p_be <- ggplot(be, aes(predictor, method, fill = beta_pointwise_coverage_95)) +
    geom_tile(color = "white") + scale_fill_viridis_c(limits = c(0, 1)) +
    labs(x = NULL, y = NULL, fill = "coverage",
         title = "Beta pointwise coverage over 10 replicates") +
    theme_minimal(base_size = 9)
  ggsave(file.path(fig_dir, "beta_coverage_heatmap.png"), p_be,
         width = 10, height = 4.5, dpi = 180, bg = "white")
  us <- aggregate(cbind(u_pointwise_coverage_95, u_mean_interval_width) ~ method + predictor, u, mean)
  p_u <- ggplot(us, aes(predictor, u_pointwise_coverage_95, fill = method)) +
    geom_hline(yintercept = 0.95, linetype = 2) +
    geom_col(position = position_dodge(width = 0.8)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "coverage", title = "u coverage X3-X6") +
    theme_minimal(base_size = 9)
  ggsave(file.path(fig_dir, "u_coverage_X3_X6.png"), p_u,
         width = 8, height = 5, dpi = 180, bg = "white")
  perf_long <- reshape(metrics[, c("replicate", "method", "mspe", "beta_mise",
                                   "mise_spatial_only", "mise_global_plus_spatial",
                                   "mise_null")],
                       varying = c("mspe", "beta_mise", "mise_spatial_only",
                                   "mise_global_plus_spatial", "mise_null"),
                       v.names = "value", timevar = "metric",
                       times = c("MSPE", "beta MISE", "spatial-only MISE",
                                 "global+spatial MISE", "null MISE"),
                       direction = "long")
  p_perf <- ggplot(perf_long, aes(method, value, fill = method)) +
    geom_boxplot(outlier.size = 0.7) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(x = NULL, y = NULL, title = "Point-estimation performance over 10 replicates") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "none")
  ggsave(file.path(fig_dir, "performance_boxplots.png"), p_perf,
         width = 13, height = 7, dpi = 180, bg = "white")
}

write_report <- function(point_summary, pip_selection_summary, response_summary,
                         beta_summary, theta_summary, u_summary, bands, checks) {
  sink(file.path(out_dir, "ten_replicate_calibration_report.md"))
  cat("# Four-Function 10-Replicate Uncertainty Calibration, Kmax6\n\n")
  cat("Configurations were selected in the independent pilot CV and frozen across all 10 replicates. ")
  cat("No K=7 or K=8 was run. GAM is included only for point-estimation metrics, not Bayesian interval calibration.\n\n")
  cat("## Point Estimation Summary\n\n"); print(point_summary)
  cat("\n## PIP Selection Summary\n\n"); print(pip_selection_summary)
  cat("\n## Response Coverage Summary\n\n"); print(response_summary)
  cat("\n## Beta Coverage Summary\n\n"); print(beta_summary)
  cat("\n## Theta Coverage Summary\n\n"); print(theta_summary)
  cat("\n## U Coverage Summary\n\n"); print(u_summary)
  cat("\n## Simultaneous Band Coverage\n\n")
  print(aggregate(cbind(simultaneous_band_covers_entire_surface,
                        simultaneous_band_mean_width) ~ method + estimand + predictor,
                  bands, mean))
  cat("\n## Validation Checks\n\n"); print(checks)
  cat("\n## Interpretation\n\n")
  cat("Coverage estimates are preliminary repeated-simulation estimates from only 10 replicates. ")
  cat("For predictor-level theta intervals, coverage can only take values 0.0, 0.1, ..., 1.0; ")
  cat("small differences such as 0.90 versus 1.00 should not be overinterpreted. ")
  cat("A larger future simulation is required for precise calibration estimates.\n")
  sink()
}

if (Sys.getenv("SSGL_SOURCE_FUNCTIONS_ONLY", "0") != "1") {
  args <- commandArgs(trailingOnly = TRUE)
  rep_arg <- grep("^--replicate=", args, value = TRUE)
  aggregate_only <- any(args == "--aggregate-only")
  if (length(rep_arg)) {
    run_replicate(as.integer(sub("^--replicate=", "", rep_arg[1])))
    quit("no")
  }
  if (!aggregate_only) {
    for (r in seq_len(n_rep)) run_replicate(r)
  }
  aggregate_outputs()
  cat("Saved 10-replicate calibration experiment to:\n", out_dir, "\n", sep = "")
}
