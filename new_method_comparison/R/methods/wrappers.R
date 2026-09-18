find_comparison_root <- function() {
  candidates <- c(getwd(), file.path(getwd(), "new_method_comparison"))
  hit <- candidates[file.exists(file.path(candidates, "R", "common", "common.R"))]
  if (!length(hit)) stop("Run from the repository root or new_method_comparison.")
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

.comparison_root <- find_comparison_root()
.repo_root <- normalizePath(file.path(.comparison_root, ".."), winslash = "/")
source(file.path(.comparison_root, "R", "common", "common.R"))

draw_mvn_precision <- function(precision, rhs, sigma2 = 1) {
  R <- chol(precision)
  mean <- backsolve(R, forwardsolve(t(R), rhs))
  mean + backsolve(R, rnorm(length(rhs)) * sqrt(sigma2))
}

fit_blasso_core <- function(y, X, n_iter, burn_in, a_sigma, b_sigma,
                            a_theta, b_theta, initialization = list()) {
  n <- nrow(X); p <- ncol(X); keep <- n_iter - burn_in
  theta <- initialization$theta %||% rep(0, p)
  omega <- initialization$omega %||% rep(1, p)
  sigma2 <- initialization$sigma2 %||% var(y)
  lambda2 <- initialization$lambda2 %||% 1
  draws <- matrix(0, p, keep); sigma_draws <- numeric(keep)
  for (iter in seq_len(n_iter)) {
    precision <- crossprod(X) + diag(1 / pmax(omega, 1e-10), p)
    theta <- draw_mvn_precision(precision, crossprod(X, y), sigma2)
    residual <- y - as.vector(X %*% theta)
    rate <- b_sigma + 0.5 * (sum(residual^2) + sum(theta^2 / omega))
    sigma2 <- 1 / rgamma(1, a_sigma + 0.5 * (n + p), rate = rate)
    omega <- vapply(seq_len(p), function(j) {
      GIGrvg::rgig(1, lambda = 0.5,
                   chi = max(theta[j]^2 / sigma2, 1e-12),
                   psi = lambda2)
    }, numeric(1))
    lambda2 <- rgamma(1, a_theta + p,
                      rate = b_theta + 0.5 * sum(omega))
    if (iter > burn_in) {
      k <- iter - burn_in
      draws[, k] <- theta
      sigma_draws[k] <- sigma2
    }
  }
  list(theta = draws, sigma2 = sigma_draws)
}

fit_blasso_core_intercept <- function(y, X, n_iter, burn_in, a_sigma, b_sigma,
                                      a_theta, b_theta,
                                      initialization = list()) {
  n <- nrow(X); p <- ncol(X); keep <- n_iter - burn_in
  W <- cbind(intercept = 1, X)
  beta0 <- initialization$beta0 %||% mean(y)
  theta <- initialization$theta %||% rep(0, p)
  omega <- initialization$omega %||% rep(1, p)
  sigma2 <- initialization$sigma2 %||% var(y)
  lambda2 <- initialization$lambda2 %||% 1
  beta0_draws <- numeric(keep)
  theta_draws <- matrix(0, p, keep)
  sigma_draws <- numeric(keep)
  for (iter in seq_len(n_iter)) {
    prior <- c(0, 1 / pmax(omega, 1e-10))
    xi <- draw_mvn_precision(
      crossprod(W) + diag(prior, p + 1L),
      crossprod(W, y),
      sigma2
    )
    beta0 <- xi[1]
    theta <- xi[-1]
    residual <- y - as.vector(W %*% xi)
    rate <- b_sigma + 0.5 * (sum(residual^2) + sum(theta^2 / omega))
    sigma2 <- 1 / rgamma(1, a_sigma + 0.5 * (n + p), rate = rate)
    omega <- vapply(seq_len(p), function(j) {
      GIGrvg::rgig(1, lambda = 0.5,
                   chi = max(theta[j]^2 / sigma2, 1e-12),
                   psi = lambda2)
    }, numeric(1))
    lambda2 <- rgamma(1, a_theta + p,
                      rate = b_theta + 0.5 * sum(omega))
    if (iter > burn_in) {
      k <- iter - burn_in
      beta0_draws[k] <- beta0
      theta_draws[, k] <- theta
      sigma_draws[k] <- sigma2
    }
  }
  list(beta0 = beta0_draws, theta = theta_draws, sigma2 = sigma_draws)
}

fit_proposed_ssgl <- function(train_data, test_data, grid_data,
                              basis_config = list(), model_config = list(),
                              mcmc_config = list(), seed = 1L) {
  set.seed(seed)
  start <- proc.time()[3]
  p <- ncol(train_data$X)
  n_iter <- mcmc_config$n_iter %||% 400L
  burn_in <- mcmc_config$burn_in %||% 150L
  meta <- fit_basis_metadata(train_data$coords,
                             basis_config$n_basis %||% 4L, centered = TRUE)
  phi <- apply_basis_metadata(train_data$coords, meta)
  center_error <- max(abs(colMeans(phi)))
  if (center_error >= 1e-8) stop("Centered training basis check failed.")
  Z <- build_spatial_design(train_data$X, phi)
  W <- cbind(train_data$X, Z)
  h <- ncol(phi); q <- p + p * h; keep <- n_iter - burn_in
  lambda0 <- model_config$lambda0 %||% 10
  lambda1 <- model_config$lambda1 %||% 1
  a_sigma <- model_config$a_sigma %||% 0.5
  b_sigma <- model_config$b_sigma %||% (var(train_data$y) / 2)
  a_theta <- model_config$a_theta %||% 1
  b_theta <- model_config$b_theta %||% 1
  a_gamma <- model_config$a_gamma %||% 1
  b_gamma <- model_config$b_gamma %||% p
  initialization <- model_config$initialization %||% list()
  xi <- initialization$xi %||% rep(0, q)
  omega <- initialization$omega %||% rep(1, p)
  gamma <- initialization$gamma %||% rep(0, p)
  tau <- initialization$tau %||% vapply(gamma, function(g) {
    lambda <- if (g > 0.5) lambda1 else lambda0
    (h + 1) / lambda^2
  }, numeric(1))
  sigma2 <- initialization$sigma2 %||% var(train_data$y)
  lambda_theta2 <- initialization$lambda_theta2 %||% 1
  pi_gamma <- initialization$pi_gamma %||% 0.5
  theta_draws <- matrix(0, p, keep)
  alpha_draws <- matrix(0, p * h, keep)
  gamma_draws <- matrix(0, p, keep)
  gamma_prob_draws <- matrix(0, p, keep)
  for (iter in seq_len(n_iter)) {
    prior <- c(1 / pmax(omega, 1e-10),
               rep(1 / pmax(tau, 1e-10), each = h))
    xi <- draw_mvn_precision(crossprod(W) + diag(prior, q),
                             crossprod(W, train_data$y), sigma2)
    theta <- xi[seq_len(p)]
    alpha <- xi[p + seq_len(p * h)]
    residual <- train_data$y - as.vector(W %*% xi)
    penalty <- sum(theta^2 / omega)
    for (j in seq_len(p)) {
      idx <- ((j - 1L) * h + 1L):(j * h)
      penalty <- penalty + sum(alpha[idx]^2) / tau[j]
    }
    sigma2 <- 1 / rgamma(
      1, a_sigma + 0.5 * (nrow(W) + q),
      rate = b_sigma + 0.5 * (sum(residual^2) + penalty)
    )
    omega <- vapply(seq_len(p), function(j) {
      GIGrvg::rgig(1, 0.5, max(theta[j]^2 / sigma2, 1e-12),
                   lambda_theta2)
    }, numeric(1))
    lambda_theta2 <- rgamma(1, a_theta + p,
                            rate = b_theta + 0.5 * sum(omega))
    pi_gamma <- rbeta(1, a_gamma + sum(gamma),
                      b_gamma + p - sum(gamma))
    for (j in seq_len(p)) {
      idx <- ((j - 1L) * h + 1L):(j * h)
      norm_scaled <- sqrt(sum(alpha[idx]^2) / sigma2)
      log_slab <- log(max(pi_gamma, 1e-12)) +
        h * log(lambda1) - lambda1 * norm_scaled
      log_spike <- log(max(1 - pi_gamma, 1e-12)) +
        h * log(lambda0) - lambda0 * norm_scaled
      prob <- plogis(log_slab - log_spike)
      gamma[j] <- rbinom(1, 1, prob)
      lambda <- if (gamma[j] == 1) lambda1 else lambda0
      tau[j] <- GIGrvg::rgig(1, 0.5,
                             max(sum(alpha[idx]^2) / sigma2, 1e-12),
                             lambda^2)
      if (iter > burn_in) gamma_prob_draws[j, iter - burn_in] <- prob
    }
    if (iter > burn_in) {
      k <- iter - burn_in
      theta_draws[, k] <- theta
      alpha_draws[, k] <- alpha
      gamma_draws[, k] <- gamma
    }
  }
  theta_mean <- rowMeans(theta_draws)
  alpha_mean <- rowMeans(alpha_draws)
  phi_test <- apply_basis_metadata(test_data$coords, meta)
  phi_grid <- apply_basis_metadata(grid_data$coords, meta)
  beta_test <- surface_from_components(theta_mean, alpha_mean, phi_test, p)
  beta_grid <- surface_from_components(theta_mean, alpha_mean, phi_grid, p)
  standard_result(
    "proposed_ssgl", theta_mean, theta_draws, rowMeans(gamma_draws),
    beta_grid, sweep(beta_grid, 2, theta_mean, "-"),
    rowSums(test_data$X * beta_test), proc.time()[3] - start,
    diagnostics = list(
      center_error = center_error,
      rb_pip = rowMeans(gamma_prob_draws),
      gamma_draws = gamma_draws,
      gamma_prob_draws = gamma_prob_draws,
      alpha_draws = alpha_draws
    ),
    config = list(basis = meta, model = model_config, mcmc = mcmc_config)
  )
}

fit_original_ssgl <- function(train_data, test_data, grid_data,
                              basis_config = list(), model_config = list(),
                              mcmc_config = list(), seed = 1L) {
  set.seed(seed)
  start <- proc.time()[3]
  suppressPackageStartupMessages(library(GIGrvg))
  if (!exists("ssgl_cpp", mode = "function")) {
    rtools_bins <- c(
      "C:/RBuildTools/4.4/usr/bin",
      "C:/RBuildTools/4.4/x86_64-w64-mingw32.static.posix/bin"
    )
    existing_bins <- rtools_bins[file.exists(rtools_bins)]
    if (length(existing_bins)) {
      Sys.setenv(PATH = paste(c(existing_bins, Sys.getenv("PATH")),
                              collapse = .Platform$path.sep))
    }
    makevars <- file.path(.comparison_root, "configs", "Makevars.win")
    if (.Platform$OS.type == "windows" && file.exists(makevars)) {
      Sys.setenv(R_MAKEVARS_USER = makevars)
      toolchain_bin <- "C:/RBuildTools/4.4/x86_64-w64-mingw32.static.posix/bin"
      Sys.setenv(COMPILER_PATH = toolchain_bin)
    }
    Rcpp::sourceCpp(file.path(.repo_root, "ssgl", "spatial_ssgl",
                              "ssgl_cpp.cpp"))
  }
  meta <- fit_basis_metadata(train_data$coords,
                             basis_config$n_basis %||% 4L, centered = FALSE)
  phi <- apply_basis_metadata(train_data$coords, meta)
  fit <- ssgl_cpp(
    as.numeric(train_data$y), as.matrix(train_data$X), phi,
    mcmc_config$n_iter %||% 400L, mcmc_config$burn_in %||% 150L,
    model_config$a_gamma %||% 1, model_config$b_gamma %||% ncol(train_data$X),
    model_config$a_sigma %||% 0.5,
    model_config$b_sigma %||% (var(train_data$y) / 2),
    model_config$lambda0 %||% 10, model_config$lambda1 %||% 1,
    model_config$zeta0 %||% 0.1, model_config$zeta1 %||% 50
  )
  eta <- rowMeans(fit$eta)
  p <- ncol(train_data$X); h <- ncol(phi)
  make_beta <- function(coords) {
    ph <- apply_basis_metadata(coords, meta)
    out <- matrix(0, nrow(ph), p)
    for (j in seq_len(p)) {
      idx <- ((j - 1L) * h + 1L):(j * h)
      out[, j] <- as.vector(ph %*% eta[idx])
    }
    out
  }
  beta_test <- make_beta(test_data$coords)
  beta_grid <- make_beta(grid_data$coords)
  standard_result(
    "original_ssgl", NULL, NULL, rowMeans(fit$gamma),
    beta_grid, NULL, rowSums(test_data$X * beta_test),
    proc.time()[3] - start,
    diagnostics = list(
      pip_target = "whole coefficient surface",
      rb_pip = rowMeans(fit$gamma_prob),
      gamma_draws = fit$gamma,
      gamma_prob_draws = fit$gamma_prob
    ),
    config = list(basis = meta, source = "ssgl/spatial_ssgl/ssgl_cpp.cpp")
  )
}

fit_original_ssgl_intercept <- function(train_data, test_data, grid_data,
                                        basis_config = list(),
                                        model_config = list(),
                                        mcmc_config = list(),
                                        seed = 1L) {
  set.seed(seed)
  start <- proc.time()[3]
  suppressPackageStartupMessages(library(GIGrvg))
  p <- ncol(train_data$X)
  n_iter <- mcmc_config$n_iter %||% 400L
  burn_in <- mcmc_config$burn_in %||% 150L
  keep <- n_iter - burn_in
  meta <- fit_basis_metadata(train_data$coords,
                             basis_config$n_basis %||% 4L,
                             centered = FALSE)
  phi <- apply_basis_metadata(train_data$coords, meta)
  h <- ncol(phi)
  initialization <- model_config$initialization %||% list()
  if (!exists("ssgl_intercept_cpp", mode = "function")) {
    Rcpp::sourceCpp(file.path(.comparison_root, "R", "methods",
                              "accelerated",
                              "original_ssgl_intercept.cpp"))
  }
  fit <- ssgl_intercept_cpp(
    as.numeric(train_data$y), as.matrix(train_data$X), phi,
    n_iter, burn_in,
    model_config$a_gamma %||% 1, model_config$b_gamma %||% p,
    model_config$a_sigma %||% 0.5,
    model_config$b_sigma %||% (var(train_data$y) / 2),
    model_config$lambda0 %||% 10, model_config$lambda1 %||% 1,
    initialization$beta0 %||% mean(train_data$y),
    model_config$zeta0 %||% 0.1, model_config$zeta1 %||% 50
  )
  eta_mean <- rowMeans(fit$eta)
  make_beta <- function(coords) {
    ph <- apply_basis_metadata(coords, meta)
    out <- matrix(0, nrow(ph), p)
    for (j in seq_len(p)) {
      idx <- ((j - 1L) * h + 1L):(j * h)
      out[, j] <- as.vector(ph %*% eta_mean[idx])
    }
    out
  }
  beta_test <- make_beta(test_data$coords)
  beta_grid <- make_beta(grid_data$coords)
  beta0_mean <- mean(fit$beta0)
  result <- standard_result(
    "original_ssgl_intercept",
    NULL, NULL, rowMeans(fit$gamma), beta_grid, NULL,
    beta0_mean + rowSums(test_data$X * beta_test),
    proc.time()[3] - start,
    diagnostics = list(
      beta0_mean = beta0_mean,
      beta0_draws = fit$beta0,
      eta_draws = fit$eta,
      sigma2_draws = fit$sigma2,
      gamma_draws = fit$gamma,
      gamma_prob_draws = fit$gamma_prob,
      rb_pip = rowMeans(fit$gamma_prob),
      theta_inclusion_draws = fit$theta,
      zeta0_draws = fit$zeta0,
      zeta1_draws = fit$zeta1,
      intercept_prior = "p(beta0) proportional to 1",
      beta0_penalized = FALSE,
      selection_target = "whole coefficient surface",
      sampler = "ssgl_intercept_cpp"
    ),
    config = list(basis = meta, model = model_config, mcmc = mcmc_config)
  )
  result$beta0_mean <- beta0_mean
  result$beta0_draws <- fit$beta0
  result$eta <- fit$eta
  result$gamma <- fit$gamma
  result$gamma_prob <- fit$gamma_prob
  result$sigma2 <- fit$sigma2
  result$meta <- meta
  return(result)

  Z <- build_spatial_design(train_data$X, phi)
  W <- cbind(intercept = 1, Z)
  q <- ncol(W)

  a_gamma <- model_config$a_gamma %||% 1
  b_gamma <- model_config$b_gamma %||% p
  a_sigma <- model_config$a_sigma %||% 0.5
  b_sigma <- model_config$b_sigma %||% (var(train_data$y) / 2)
  lambda0 <- model_config$lambda0 %||% 10
  lambda1 <- model_config$lambda1 %||% 1
  initialization <- model_config$initialization %||% list()

  beta0 <- initialization$beta0 %||% mean(train_data$y)
  eta <- initialization$eta %||% rep(0, p * h)
  xi <- initialization$xi %||% c(beta0, eta)
  if (length(xi) != q) {
    stop("Initial xi must have length 1 + p*h for intercept Original SSGL.")
  }
  gamma <- initialization$gamma %||% rbinom(p, 1, 0.5)
  zeta0 <- initialization$zeta0 %||% (model_config$zeta0 %||% 0.1)
  zeta1 <- initialization$zeta1 %||% (model_config$zeta1 %||% 50)
  sigma2 <- initialization$sigma2 %||% var(train_data$y)
  theta_inclusion <- initialization$theta %||% 0.5

  beta0_draws <- numeric(keep)
  eta_draws <- matrix(0, p * h, keep)
  gamma_draws <- matrix(0, p, keep)
  gamma_prob_draws <- matrix(0, p, keep)
  sigma2_draws <- numeric(keep)
  theta_draws <- numeric(keep)
  zeta0_draws <- numeric(keep)
  zeta1_draws <- numeric(keep)

  for (iter in seq_len(n_iter)) {
    zeta_by_group <- ifelse(gamma > 0.5, zeta1, zeta0)
    prior <- c(0, rep(1 / pmax(zeta_by_group, 1e-12), each = h))
    xi <- draw_mvn_precision(crossprod(W) / sigma2 + diag(prior, q),
                             crossprod(W, train_data$y) / sigma2,
                             sigma2 = 1)
    beta0 <- xi[1]
    eta <- xi[-1]

    residual <- train_data$y - as.vector(W %*% xi)
    sigma2 <- 1 / rgamma(
      1,
      a_sigma + 0.5 * nrow(W),
      rate = b_sigma + 0.5 * sum(residual^2)
    )

    theta_inclusion <- rbeta(1, a_gamma + sum(gamma),
                             b_gamma + p - sum(gamma))
    eta_squares <- vapply(seq_len(p), function(j) {
      idx <- ((j - 1L) * h + 1L):(j * h)
      sum(eta[idx]^2)
    }, numeric(1))

    if (any(gamma == 1)) {
      sum_eta_squared <- sum(eta_squares[gamma == 1])
      n_selected <- sum(gamma)
      zeta_lambda <- (h + 1.0 - n_selected * h) / 2.0
      zeta1 <- GIGrvg::rgig(
        1, lambda = zeta_lambda,
        chi = max(sum_eta_squared, 1e-12),
        psi = lambda1 * lambda1
      )
    }
    if (any(gamma == 0)) {
      sum_eta_squared <- sum(eta_squares[gamma == 0])
      n_unselected <- p - sum(gamma)
      zeta_lambda <- (h + 1.0 - n_unselected * h) / 2.0
      zeta0 <- GIGrvg::rgig(
        1, lambda = zeta_lambda,
        chi = max(sum_eta_squared, 1e-12),
        psi = lambda0 * lambda0
      )
    }

    for (j in seq_len(p)) {
      log_slab <- log(max(theta_inclusion, 1e-12)) -
        0.5 * h * log(max(zeta1, 1e-12)) -
        0.5 * eta_squares[j] / max(zeta1, 1e-12)
      log_spike <- log(max(1 - theta_inclusion, 1e-12)) -
        0.5 * h * log(max(zeta0, 1e-12)) -
        0.5 * eta_squares[j] / max(zeta0, 1e-12)
      prob <- plogis(log_slab - log_spike)
      gamma[j] <- rbinom(1, 1, prob)
      if (iter > burn_in) gamma_prob_draws[j, iter - burn_in] <- prob
    }

    if (iter > burn_in) {
      k <- iter - burn_in
      beta0_draws[k] <- beta0
      eta_draws[, k] <- eta
      gamma_draws[, k] <- gamma
      sigma2_draws[k] <- sigma2
      theta_draws[k] <- theta_inclusion
      zeta0_draws[k] <- zeta0
      zeta1_draws[k] <- zeta1
    }
  }

  eta_mean <- rowMeans(eta_draws)
  make_beta <- function(coords) {
    ph <- apply_basis_metadata(coords, meta)
    out <- matrix(0, nrow(ph), p)
    for (j in seq_len(p)) {
      idx <- ((j - 1L) * h + 1L):(j * h)
      out[, j] <- as.vector(ph %*% eta_mean[idx])
    }
    out
  }
  beta_test <- make_beta(test_data$coords)
  beta_grid <- make_beta(grid_data$coords)
  beta0_mean <- mean(beta0_draws)
  result <- standard_result(
    "original_ssgl_intercept",
    NULL, NULL, rowMeans(gamma_draws), beta_grid, NULL,
    beta0_mean + rowSums(test_data$X * beta_test),
    proc.time()[3] - start,
    diagnostics = list(
      beta0_mean = beta0_mean,
      beta0_draws = beta0_draws,
      eta_draws = eta_draws,
      sigma2_draws = sigma2_draws,
      gamma_draws = gamma_draws,
      gamma_prob_draws = gamma_prob_draws,
      rb_pip = rowMeans(gamma_prob_draws),
      theta_inclusion_draws = theta_draws,
      zeta0_draws = zeta0_draws,
      zeta1_draws = zeta1_draws,
      intercept_prior = "p(beta0) proportional to 1",
      beta0_penalized = FALSE,
      selection_target = "whole coefficient surface"
    ),
    config = list(basis = meta, model = model_config, mcmc = mcmc_config)
  )
  result$beta0_mean <- beta0_mean
  result$beta0_draws <- beta0_draws
  result$eta <- eta_draws
  result$gamma <- gamma_draws
  result$gamma_prob <- gamma_prob_draws
  result$sigma2 <- sigma2_draws
  result$meta <- meta
  result
}

fit_global_only_blasso <- function(train_data, test_data, grid_data,
                                   basis_config = list(), model_config = list(),
                                   mcmc_config = list(), seed = 1L) {
  set.seed(seed)
  start <- proc.time()[3]
  fit <- fit_blasso_core(
    train_data$y, train_data$X,
    mcmc_config$n_iter %||% 400L, mcmc_config$burn_in %||% 150L,
    model_config$a_sigma %||% 0.5,
    model_config$b_sigma %||% (var(train_data$y) / 2),
    model_config$a_theta %||% 1, model_config$b_theta %||% 1,
    model_config$initialization %||% list()
  )
  theta <- rowMeans(fit$theta)
  beta_grid <- matrix(rep(theta, each = nrow(grid_data$coords)),
                      nrow(grid_data$coords), length(theta))
  standard_result(
    "global_only_blasso", theta, fit$theta, NULL, beta_grid,
    matrix(0, nrow(beta_grid), ncol(beta_grid)),
    as.vector(test_data$X %*% theta), proc.time()[3] - start,
    config = list(model = model_config, mcmc = mcmc_config)
  )
}

fit_global_only_blasso_intercept <- function(train_data, test_data, grid_data,
                                             basis_config = list(),
                                             model_config = list(),
                                             mcmc_config = list(),
                                             seed = 1L) {
  set.seed(seed)
  start <- proc.time()[3]
  fit <- fit_blasso_core_intercept(
    train_data$y, train_data$X,
    mcmc_config$n_iter %||% 400L, mcmc_config$burn_in %||% 150L,
    model_config$a_sigma %||% 0.5,
    model_config$b_sigma %||% (var(train_data$y) / 2),
    model_config$a_theta %||% 1, model_config$b_theta %||% 1,
    model_config$initialization %||% list()
  )
  beta0 <- mean(fit$beta0)
  theta <- rowMeans(fit$theta)
  beta_grid <- matrix(rep(theta, each = nrow(grid_data$coords)),
                      nrow(grid_data$coords), length(theta))
  pred_test <- beta0 + as.vector(test_data$X %*% theta)
  result <- standard_result(
    "global_only_blasso_intercept", theta, fit$theta, NULL, beta_grid,
    matrix(0, nrow(beta_grid), ncol(beta_grid)),
    pred_test, proc.time()[3] - start,
    diagnostics = list(
      beta0_mean = beta0,
      beta0_draws = fit$beta0,
      sigma2_draws = fit$sigma2,
      intercept_prior = "p(beta0) proportional to 1",
      beta0_penalized = FALSE
    ),
    config = list(model = model_config, mcmc = mcmc_config)
  )
  result$beta0_mean <- beta0
  result$beta0_draws <- fit$beta0
  result
}

fit_full_svc_no_selection <- function(train_data, test_data, grid_data,
                                      basis_config = list(),
                                      model_config = list(),
                                      mcmc_config = list(), seed = 1L) {
  set.seed(seed)
  start <- proc.time()[3]
  p <- ncol(train_data$X)
  n_iter <- mcmc_config$n_iter %||% 400L
  burn_in <- mcmc_config$burn_in %||% 150L
  meta <- fit_basis_metadata(train_data$coords,
                             basis_config$n_basis %||% 4L, centered = TRUE)
  phi <- apply_basis_metadata(train_data$coords, meta)
  center_error <- max(abs(colMeans(phi)))
  if (center_error >= 1e-8) stop("Centered training basis check failed.")
  Z <- build_spatial_design(train_data$X, phi)
  W <- cbind(train_data$X, Z)
  h <- ncol(phi); q <- ncol(W); keep <- n_iter - burn_in
  kappa2 <- model_config$kappa2_alpha %||% 1
  initialization <- model_config$initialization %||% list()
  omega <- initialization$omega %||% rep(1, p)
  xi <- initialization$xi %||% rep(0, q)
  sigma2 <- initialization$sigma2 %||% var(train_data$y)
  lambda2 <- initialization$lambda2 %||% 1
  theta_draws <- matrix(0, p, keep)
  alpha_draws <- matrix(0, p * h, keep)
  for (iter in seq_len(n_iter)) {
    prior <- c(1 / pmax(omega, 1e-10), rep(1 / kappa2, p * h))
    xi <- draw_mvn_precision(crossprod(W) + diag(prior, q),
                             crossprod(W, train_data$y), sigma2)
    theta <- xi[seq_len(p)]
    alpha <- xi[p + seq_len(p * h)]
    residual <- train_data$y - as.vector(W %*% xi)
    penalty <- sum(theta^2 / omega) + sum(alpha^2) / kappa2
    sigma2 <- 1 / rgamma(
      1, (model_config$a_sigma %||% 0.5) + 0.5 * (nrow(W) + q),
      rate = (model_config$b_sigma %||% (var(train_data$y) / 2)) +
        0.5 * (sum(residual^2) + penalty)
    )
    omega <- vapply(seq_len(p), function(j) {
      GIGrvg::rgig(1, 0.5, max(theta[j]^2 / sigma2, 1e-12), lambda2)
    }, numeric(1))
    lambda2 <- rgamma(
      1, (model_config$a_theta %||% 1) + p,
      rate = (model_config$b_theta %||% 1) + 0.5 * sum(omega)
    )
    if (iter > burn_in) {
      k <- iter - burn_in
      theta_draws[, k] <- theta
      alpha_draws[, k] <- alpha
    }
  }
  theta <- rowMeans(theta_draws); alpha <- rowMeans(alpha_draws)
  phi_test <- apply_basis_metadata(test_data$coords, meta)
  phi_grid <- apply_basis_metadata(grid_data$coords, meta)
  beta_test <- surface_from_components(theta, alpha, phi_test, p)
  beta_grid <- surface_from_components(theta, alpha, phi_grid, p)
  standard_result(
    "full_svc_no_selection", theta, theta_draws, NULL, beta_grid,
    sweep(beta_grid, 2, theta, "-"), rowSums(test_data$X * beta_test),
    proc.time()[3] - start,
    diagnostics = list(center_error = center_error),
    config = list(basis = meta, kappa2_alpha = kappa2)
  )
}

fit_full_svc_no_selection_intercept <- function(train_data, test_data, grid_data,
                                                basis_config = list(),
                                                model_config = list(),
                                                mcmc_config = list(),
                                                seed = 1L) {
  set.seed(seed)
  start <- proc.time()[3]
  p <- ncol(train_data$X)
  n_iter <- mcmc_config$n_iter %||% 400L
  burn_in <- mcmc_config$burn_in %||% 150L
  full_rank_centered <- basis_config$full_rank_centered %||% TRUE
  if (isTRUE(full_rank_centered) &&
      exists("fit_centered_basis_metadata_fullrank", mode = "function") &&
      exists("apply_centered_basis_metadata_fullrank", mode = "function")) {
    meta <- fit_centered_basis_metadata_fullrank(
      train_data$coords,
      basis_config$n_basis %||% 4L,
      basis_config$rank_tolerance %||% 1e-10,
      basis_config$full_rank_method %||% "svd"
    )
    phi <- apply_centered_basis_metadata_fullrank(train_data$coords, meta)
  } else {
    meta <- fit_basis_metadata(train_data$coords,
                               basis_config$n_basis %||% 4L,
                               centered = TRUE)
    phi <- apply_basis_metadata(train_data$coords, meta)
  }
  center_error <- max(abs(colMeans(phi)))
  if (center_error >= 1e-8) stop("Centered training basis check failed.")
  Z <- build_spatial_design(train_data$X, phi)
  W <- cbind(intercept = 1, train_data$X, Z)
  h <- ncol(phi); q <- ncol(W); keep <- n_iter - burn_in
  q_penalized <- q - 1L
  kappa2 <- model_config$kappa2_alpha %||% 1
  initialization <- model_config$initialization %||% list()
  omega <- initialization$omega %||% rep(1, p)
  init_alpha <- initialization$alpha
  if (is.null(init_alpha)) init_alpha <- rep(0, p * h)
  xi <- initialization$xi %||% c(
    initialization$beta0 %||% mean(train_data$y),
    initialization$theta %||% rep(0, p),
    init_alpha
  )
  if (length(xi) != q) {
    stop("Initial xi must have length 1 + p + p*h for intercept Full SVC.")
  }
  sigma2 <- initialization$sigma2 %||% var(train_data$y)
  lambda2 <- initialization$lambda2 %||% 1
  if (!exists("full_svc_intercept_cpp", mode = "function")) {
    Rcpp::sourceCpp(file.path(.comparison_root, "R", "methods",
                              "accelerated", "full_svc_intercept.cpp"))
  }
  draws <- full_svc_intercept_cpp(
    as.numeric(train_data$y), W, p, h, n_iter, burn_in,
    model_config$a_sigma %||% 0.5,
    model_config$b_sigma %||% (var(train_data$y) / 2),
    model_config$a_theta %||% 1,
    model_config$b_theta %||% 1,
    kappa2, xi, omega, sigma2, lambda2
  )
  beta0 <- mean(draws$beta0)
  theta <- rowMeans(draws$theta)
  alpha <- rowMeans(draws$alpha)
  if (isTRUE(full_rank_centered) &&
      exists("apply_centered_basis_metadata_fullrank", mode = "function") &&
      isTRUE(meta$full_rank_centered)) {
    phi_test <- apply_centered_basis_metadata_fullrank(test_data$coords, meta)
    phi_grid <- apply_centered_basis_metadata_fullrank(grid_data$coords, meta)
  } else {
    phi_test <- apply_basis_metadata(test_data$coords, meta)
    phi_grid <- apply_basis_metadata(grid_data$coords, meta)
  }
  beta_test <- surface_from_components(theta, alpha, phi_test, p)
  beta_grid <- surface_from_components(theta, alpha, phi_grid, p)
  result <- standard_result(
    "full_svc_no_selection_intercept", theta, draws$theta, NULL, beta_grid,
    sweep(beta_grid, 2, theta, "-"),
    beta0 + rowSums(test_data$X * beta_test),
    proc.time()[3] - start,
    diagnostics = list(
      beta0_mean = beta0,
      beta0_draws = draws$beta0,
      alpha_draws = draws$alpha,
      sigma2_draws = draws$sigma2,
      lambda2_draws = draws$lambda2,
      omega_draws = draws$omega,
      center_error = center_error,
      intercept_prior = "p(beta0) proportional to 1",
      beta0_penalized = FALSE,
      sampler = "full_svc_intercept_cpp",
      raw_basis_dimension = meta$raw_h %||% h,
      effective_basis_dimension = h,
      centered_basis_rank = meta$centered_rank %||% qr(phi)$rank,
      rank_loss_from_centering = meta$rank_loss_from_centering %||% NA_integer_,
      full_rank_method = meta$full_rank_method %||% NA_character_
    ),
    config = list(basis = meta, kappa2_alpha = kappa2,
                  model = model_config, mcmc = mcmc_config)
  )
  result$beta0_mean <- beta0
  result$beta0_draws <- draws$beta0
  result$alpha_mean <- alpha
  return(result)

  beta0_draws <- numeric(keep)
  theta_draws <- matrix(0, p, keep)
  alpha_draws <- matrix(0, p * h, keep)
  sigma2_draws <- numeric(keep)
  for (iter in seq_len(n_iter)) {
    prior <- c(0, 1 / pmax(omega, 1e-10), rep(1 / kappa2, p * h))
    xi <- draw_mvn_precision(crossprod(W) + diag(prior, q),
                             crossprod(W, train_data$y), sigma2)
    beta0 <- xi[1]
    theta <- xi[1L + seq_len(p)]
    alpha <- xi[1L + p + seq_len(p * h)]
    residual <- train_data$y - as.vector(W %*% xi)
    penalty <- sum(theta^2 / omega) + sum(alpha^2) / kappa2
    sigma2 <- 1 / rgamma(
      1, (model_config$a_sigma %||% 0.5) +
        0.5 * (nrow(W) + q_penalized),
      rate = (model_config$b_sigma %||% (var(train_data$y) / 2)) +
        0.5 * (sum(residual^2) + penalty)
    )
    omega <- vapply(seq_len(p), function(j) {
      GIGrvg::rgig(1, 0.5, max(theta[j]^2 / sigma2, 1e-12), lambda2)
    }, numeric(1))
    lambda2 <- rgamma(
      1, (model_config$a_theta %||% 1) + p,
      rate = (model_config$b_theta %||% 1) + 0.5 * sum(omega)
    )
    if (iter > burn_in) {
      k <- iter - burn_in
      beta0_draws[k] <- beta0
      theta_draws[, k] <- theta
      alpha_draws[, k] <- alpha
      sigma2_draws[k] <- sigma2
    }
  }
  beta0 <- mean(beta0_draws)
  theta <- rowMeans(theta_draws)
  alpha <- rowMeans(alpha_draws)
  if (isTRUE(full_rank_centered) &&
      exists("apply_centered_basis_metadata_fullrank", mode = "function") &&
      isTRUE(meta$full_rank_centered)) {
    phi_test <- apply_centered_basis_metadata_fullrank(test_data$coords, meta)
    phi_grid <- apply_centered_basis_metadata_fullrank(grid_data$coords, meta)
  } else {
    phi_test <- apply_basis_metadata(test_data$coords, meta)
    phi_grid <- apply_basis_metadata(grid_data$coords, meta)
  }
  beta_test <- surface_from_components(theta, alpha, phi_test, p)
  beta_grid <- surface_from_components(theta, alpha, phi_grid, p)
  result <- standard_result(
    "full_svc_no_selection_intercept", theta, theta_draws, NULL, beta_grid,
    sweep(beta_grid, 2, theta, "-"),
    beta0 + rowSums(test_data$X * beta_test),
    proc.time()[3] - start,
    diagnostics = list(
      beta0_mean = beta0,
      beta0_draws = beta0_draws,
      alpha_draws = alpha_draws,
      sigma2_draws = sigma2_draws,
      center_error = center_error,
      intercept_prior = "p(beta0) proportional to 1",
      beta0_penalized = FALSE,
      raw_basis_dimension = meta$raw_h %||% h,
      effective_basis_dimension = h,
      centered_basis_rank = meta$centered_rank %||% qr(phi)$rank,
      rank_loss_from_centering = meta$rank_loss_from_centering %||% NA_integer_,
      full_rank_method = meta$full_rank_method %||% NA_character_
    ),
    config = list(basis = meta, kappa2_alpha = kappa2,
                  model = model_config, mcmc = mcmc_config)
  )
  result$beta0_mean <- beta0
  result$beta0_draws <- beta0_draws
  result$alpha_mean <- alpha
  result
}

fit_original_gam <- function(train_data, test_data, grid_data,
                             basis_config = list(), model_config = list(),
                             mcmc_config = list(), seed = 1L) {
  set.seed(seed)
  start <- proc.time()[3]
  p <- ncol(train_data$X)
  data <- as.data.frame(train_data$X)
  names(data) <- paste0("X", seq_len(p))
  data$y <- train_data$y
  data$x <- train_data$coords[, 1]
  data$y_coord <- train_data$coords[, 2]
  linear_terms <- paste0("X", seq_len(p), collapse = " + ")
  smooth_terms <- paste0(
    "s(x, y_coord, bs='tp', by=X", seq_len(p), ")",
    collapse = " + "
  )
  formula_str <- paste0("y ~ 0 + ", linear_terms, " + ", smooth_terms)
  fit <- mgcv::gam(as.formula(formula_str), data = data)
  test_df <- as.data.frame(test_data$X)
  names(test_df) <- paste0("X", seq_len(p))
  test_df$x <- test_data$coords[, 1]
  test_df$y_coord <- test_data$coords[, 2]
  pred <- predict(fit, newdata = test_df)

  beta_grid <- matrix(0, nrow(grid_data$coords), p)
  for (j in seq_len(p)) {
    grid_df <- data.frame(
      x = grid_data$coords[, 1],
      y_coord = grid_data$coords[, 2]
    )
    for (jj in seq_len(p)) grid_df[[paste0("X", jj)]] <- 0
    grid_df[[paste0("X", j)]] <- 1
    beta_grid[, j] <- predict(fit, newdata = grid_df)
  }
  coef_fit <- stats::coef(fit)
  theta_mean <- as.numeric(coef_fit[paste0("X", seq_len(p))])
  standard_result(
    "original_gam", theta_mean, NULL, NULL, beta_grid,
    sweep(beta_grid, 2, theta_mean, "-"), pred,
    proc.time()[3] - start,
    diagnostics = list(converged = fit$converged %||% TRUE),
    config = list(
      source = paste0(
        "standard thin-plate counterpart of ",
        "ssgl/spatial_ssgl/helpers_gam.R"
      ),
      formula = formula_str,
      basis = "thin-plate varying-coefficient smooth",
      k = "mgcv default",
      method = "GCV.Cp (mgcv::gam default)",
      coefficient_extraction = "Xj=1; all other predictors=0"
    )
  )
}
