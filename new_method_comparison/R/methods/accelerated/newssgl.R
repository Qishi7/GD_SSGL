find_accelerated_root <- function() {
  candidates <- c(
    file.path(getwd(), "new_method_comparison"),
    getwd()
  )
  hit <- candidates[file.exists(file.path(
    candidates, "R", "methods", "accelerated",
    "newssgl.cpp"
  ))]
  if (!length(hit)) stop("Cannot locate new_method_comparison.")
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

.accelerated_root <- find_accelerated_root()
source(file.path(.accelerated_root, "R", "methods", "ssgl.R"))

load_newssgl_fast <- function(rebuild = FALSE) {
  if (exists("newssgl_gibbs_cpp", mode = "function") && !rebuild) {
    return(invisible(TRUE))
  }
  makevars <- file.path(.accelerated_root, "configs", "Makevars.win")
  if (.Platform$OS.type == "windows" && file.exists(makevars)) {
    Sys.setenv(R_MAKEVARS_USER = makevars)
    rtools_bins <- c(
      "C:/RBuildTools/4.4/usr/bin",
      "C:/RBuildTools/4.4/x86_64-w64-mingw32.static.posix/bin"
    )
    rtools_bins <- rtools_bins[file.exists(rtools_bins)]
    if (length(rtools_bins)) {
      Sys.setenv(PATH = paste(c(rtools_bins, Sys.getenv("PATH")),
                              collapse = .Platform$path.sep))
    }
  }
  Rcpp::sourceCpp(file.path(
    .accelerated_root, "R", "methods", "accelerated",
    "newssgl.cpp"
  ), rebuild = rebuild)
  invisible(TRUE)
}

prepare_newssgl <- function(train_data,
                            basis_config = list(),
                            model_config = list()) {
  p <- ncol(train_data$X)
  meta <- fit_basis_metadata(
    train_data$coords,
    basis_config$n_basis %||% 4L,
    centered = TRUE
  )
  phi <- apply_basis_metadata(train_data$coords, meta)
  center_error <- max(abs(colMeans(phi)))
  if (center_error >= 1e-8) stop("Centered training basis check failed.")
  Z <- build_spatial_design(train_data$X, phi)
  W <- cbind(train_data$X, Z)
  h <- ncol(phi)
  q <- ncol(W)
  lambda0 <- model_config$lambda0 %||% 10
  lambda1 <- model_config$lambda1 %||% 1
  initialization <- model_config$initialization %||% list()
  gamma <- initialization$gamma %||% rep(0, p)
  tau <- initialization$tau %||% vapply(gamma, function(g) {
    lambda <- if (g > 0.5) lambda1 else lambda0
    (h + 1) / lambda^2
  }, numeric(1))

  list(
    p = p, h = h, q = q, meta = meta, phi = phi, W = W,
    center_error = center_error,
    parameters = list(
      lambda0 = lambda0,
      lambda1 = lambda1,
      a_sigma = model_config$a_sigma %||% 0.5,
      b_sigma = model_config$b_sigma %||% (var(train_data$y) / 2),
      a_theta = model_config$a_theta %||% 1,
      b_theta = model_config$b_theta %||% 1,
      a_gamma = model_config$a_gamma %||% 1,
      b_gamma = model_config$b_gamma %||% p
    ),
    initialization = list(
      xi = initialization$xi %||% rep(0, q),
      omega = initialization$omega %||% rep(1, p),
      gamma = gamma,
      tau = tau,
      sigma2 = initialization$sigma2 %||% var(train_data$y),
      lambda_theta2 = initialization$lambda_theta2 %||% 1,
      pi_gamma = initialization$pi_gamma %||% 0.5
    )
  )
}

fit_newssgl_fast <- function(
    train_data, test_data, grid_data,
    basis_config = list(), model_config = list(),
    mcmc_config = list(), seed = 1L) {
  load_newssgl_fast()
  set.seed(seed)
  started <- proc.time()[3]
  prepared <- prepare_newssgl(
    train_data, basis_config, model_config
  )
  pars <- prepared$parameters
  init <- prepared$initialization
  n_iter <- mcmc_config$n_iter %||% 400L
  burn_in <- mcmc_config$burn_in %||% 150L

  draws <- newssgl_gibbs_cpp(
    y = as.numeric(train_data$y),
    W = prepared$W,
    p = prepared$p,
    h = prepared$h,
    n_iter = n_iter,
    burn_in = burn_in,
    lambda0 = pars$lambda0,
    lambda1 = pars$lambda1,
    a_sigma = pars$a_sigma,
    b_sigma = pars$b_sigma,
    a_theta = pars$a_theta,
    b_theta = pars$b_theta,
    a_gamma = pars$a_gamma,
    b_gamma = pars$b_gamma,
    xi = init$xi,
    omega = init$omega,
    gamma = init$gamma,
    tau = init$tau,
    sigma2 = init$sigma2,
    lambda_theta2 = init$lambda_theta2,
    pi_gamma = init$pi_gamma
  )

  theta_mean <- rowMeans(draws$theta)
  alpha_mean <- rowMeans(draws$alpha)
  phi_test <- apply_basis_metadata(test_data$coords, prepared$meta)
  phi_grid <- apply_basis_metadata(grid_data$coords, prepared$meta)
  beta_test <- surface_from_components(
    theta_mean, alpha_mean, phi_test, prepared$p
  )
  beta_grid <- surface_from_components(
    theta_mean, alpha_mean, phi_grid, prepared$p
  )

  result <- standard_result(
    "newssgl",
    theta_mean,
    draws$theta,
    rowMeans(draws$gamma),
    beta_grid,
    sweep(beta_grid, 2, theta_mean, "-"),
    rowSums(test_data$X * beta_test),
    proc.time()[3] - started,
    diagnostics = list(
      center_error = prepared$center_error,
      rb_pip = rowMeans(draws$gamma_prob),
      gamma_draws = draws$gamma,
      gamma_prob_draws = draws$gamma_prob,
      alpha_draws = draws$alpha,
      sigma2_draws = draws$sigma2,
      omega_draws = draws$omega,
      tau_draws = draws$tau,
      lambda_theta2_draws = draws$lambda_theta2,
      pi_gamma_draws = draws$pi_gamma
    ),
    config = list(
      basis = prepared$meta,
      model = model_config,
      mcmc = mcmc_config,
      implementation = "RcppArmadillo accelerated"
    )
  )
  result$alpha_mean <- alpha_mean
  result
}

# Reference diagnostic sampler: same R update order and conditionals as the
# validated wrapper, with additional draws retained only for equivalence tests.
fit_newssgl_reference_test <- function(
    train_data, test_data, grid_data,
    basis_config = list(), model_config = list(),
    mcmc_config = list(), seed = 1L) {
  set.seed(seed)
  started <- proc.time()[3]
  prepared <- prepare_newssgl(
    train_data, basis_config, model_config
  )
  p <- prepared$p
  h <- prepared$h
  q <- prepared$q
  W <- prepared$W
  pars <- prepared$parameters
  init <- prepared$initialization
  n_iter <- mcmc_config$n_iter %||% 400L
  burn_in <- mcmc_config$burn_in %||% 150L
  keep <- n_iter - burn_in

  xi <- init$xi
  omega <- init$omega
  gamma <- init$gamma
  tau <- init$tau
  sigma2 <- init$sigma2
  lambda_theta2 <- init$lambda_theta2
  pi_gamma <- init$pi_gamma
  WtW <- crossprod(W)
  Wty <- crossprod(W, train_data$y)

  theta_draws <- matrix(0, p, keep)
  alpha_draws <- matrix(0, p * h, keep)
  gamma_draws <- matrix(0, p, keep)
  gamma_prob_draws <- matrix(0, p, keep)
  omega_draws <- matrix(0, p, keep)
  tau_draws <- matrix(0, p, keep)
  sigma2_draws <- numeric(keep)
  lambda_theta2_draws <- numeric(keep)
  pi_gamma_draws <- numeric(keep)

  for (iter in seq_len(n_iter)) {
    prior <- c(
      1 / pmax(omega, 1e-10),
      rep(1 / pmax(tau, 1e-10), each = h)
    )
    xi <- draw_mvn_precision(
      WtW + diag(prior, q), Wty, sigma2
    )
    theta <- xi[seq_len(p)]
    alpha <- xi[p + seq_len(p * h)]
    residual <- train_data$y - as.vector(W %*% xi)
    penalty <- sum(theta^2 / omega)
    for (j in seq_len(p)) {
      idx <- ((j - 1L) * h + 1L):(j * h)
      penalty <- penalty + sum(alpha[idx]^2) / tau[j]
    }
    sigma2 <- 1 / rgamma(
      1, pars$a_sigma + 0.5 * (nrow(W) + q),
      rate = pars$b_sigma +
        0.5 * (sum(residual^2) + penalty)
    )
    omega <- vapply(seq_len(p), function(j) {
      GIGrvg::rgig(
        1, 0.5, max(theta[j]^2 / sigma2, 1e-12),
        lambda_theta2
      )
    }, numeric(1))
    lambda_theta2 <- rgamma(
      1, pars$a_theta + p,
      rate = pars$b_theta + 0.5 * sum(omega)
    )
    pi_gamma <- rbeta(
      1, pars$a_gamma + sum(gamma),
      pars$b_gamma + p - sum(gamma)
    )
    for (j in seq_len(p)) {
      idx <- ((j - 1L) * h + 1L):(j * h)
      norm_scaled <- sqrt(sum(alpha[idx]^2) / sigma2)
      log_slab <- log(max(pi_gamma, 1e-12)) +
        h * log(pars$lambda1) - pars$lambda1 * norm_scaled
      log_spike <- log(max(1 - pi_gamma, 1e-12)) +
        h * log(pars$lambda0) - pars$lambda0 * norm_scaled
      probability <- plogis(log_slab - log_spike)
      gamma[j] <- rbinom(1, 1, probability)
      lambda <- if (gamma[j] == 1) pars$lambda1 else pars$lambda0
      tau[j] <- GIGrvg::rgig(
        1, 0.5, max(sum(alpha[idx]^2) / sigma2, 1e-12),
        lambda^2
      )
      if (iter > burn_in) {
        gamma_prob_draws[j, iter - burn_in] <- probability
      }
    }
    if (iter > burn_in) {
      k <- iter - burn_in
      theta_draws[, k] <- theta
      alpha_draws[, k] <- alpha
      gamma_draws[, k] <- gamma
      omega_draws[, k] <- omega
      tau_draws[, k] <- tau
      sigma2_draws[k] <- sigma2
      lambda_theta2_draws[k] <- lambda_theta2
      pi_gamma_draws[k] <- pi_gamma
    }
  }

  theta_mean <- rowMeans(theta_draws)
  alpha_mean <- rowMeans(alpha_draws)
  phi_test <- apply_basis_metadata(test_data$coords, prepared$meta)
  phi_grid <- apply_basis_metadata(grid_data$coords, prepared$meta)
  beta_test <- surface_from_components(
    theta_mean, alpha_mean, phi_test, p
  )
  beta_grid <- surface_from_components(
    theta_mean, alpha_mean, phi_grid, p
  )
  result <- standard_result(
    "newssgl_reference",
    theta_mean, theta_draws, rowMeans(gamma_draws),
    beta_grid, sweep(beta_grid, 2, theta_mean, "-"),
    rowSums(test_data$X * beta_test),
    proc.time()[3] - started,
    diagnostics = list(
      center_error = prepared$center_error,
      rb_pip = rowMeans(gamma_prob_draws),
      gamma_draws = gamma_draws,
      gamma_prob_draws = gamma_prob_draws,
      alpha_draws = alpha_draws,
      sigma2_draws = sigma2_draws,
      omega_draws = omega_draws,
      tau_draws = tau_draws,
      lambda_theta2_draws = lambda_theta2_draws,
      pi_gamma_draws = pi_gamma_draws
    ),
    config = list(
      basis = prepared$meta,
      model = model_config,
      mcmc = mcmc_config,
      implementation = "R reference diagnostic"
    )
  )
  result$alpha_mean <- alpha_mean
  result
}

# Current project default after validation:
# keep the original R implementation available, but route future newssgl calls
# through the validated RcppArmadillo implementation.
fit_newssgl_reference <- fit_proposed_ssgl
fit_newssgl <- fit_newssgl_fast
