find_matched_wsssgl_root <- function() {
  candidates <- c(
    file.path(getwd(), "new_method_comparison"),
    getwd()
  )
  hit <- candidates[file.exists(file.path(
    candidates, "R", "methods", "accelerated",
    "matched_wsssgl_intercept.cpp"
  ))]
  if (!length(hit)) stop("Cannot locate new_method_comparison.")
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

.matched_wsssgl_root <- find_matched_wsssgl_root()
source(file.path(.matched_wsssgl_root, "R", "methods", "accelerated",
                 "newssgl_intercept.R"))

load_matched_wsssgl_intercept_fast <- function(rebuild = FALSE) {
  if (exists("matched_wsssgl_intercept_gibbs_cpp", mode = "function") &&
      !rebuild) {
    return(invisible(TRUE))
  }
  makevars <- file.path(.matched_wsssgl_root, "configs", "Makevars.win")
  if (.Platform$OS.type == "windows" && file.exists(makevars)) {
    Sys.setenv(R_MAKEVARS_USER = makevars)
  }
  Rcpp::sourceCpp(file.path(
    .matched_wsssgl_root, "R", "methods", "accelerated",
    "matched_wsssgl_intercept.cpp"
  ), rebuild = rebuild)
  invisible(TRUE)
}

build_whole_surface_basis <- function(coords, metadata) {
  phi <- apply_centered_basis_metadata_fullrank(coords, metadata)
  cbind(constant = 1, phi)
}

build_whole_surface_design <- function(X, q_basis) {
  build_spatial_design(X, q_basis)
}

surface_from_eta <- function(eta, q_basis, p) {
  H <- ncol(q_basis)
  beta <- matrix(0, nrow(q_basis), p)
  for (j in seq_len(p)) {
    idx <- ((j - 1L) * H + 1L):(j * H)
    beta[, j] <- as.vector(q_basis %*% eta[idx])
  }
  beta
}

split_eta_whole_surface <- function(eta, p, H) {
  theta <- numeric(p)
  alpha <- numeric(p * (H - 1L))
  for (j in seq_len(p)) {
    eta_idx <- ((j - 1L) * H + 1L):(j * H)
    alpha_idx <- ((j - 1L) * (H - 1L) + 1L):(j * (H - 1L))
    theta[j] <- eta[eta_idx[1]]
    alpha[alpha_idx] <- eta[eta_idx[-1]]
  }
  list(theta = theta, alpha = alpha)
}

prepare_matched_wsssgl_intercept <- function(train_data,
                                             basis_config = list(),
                                             model_config = list()) {
  p <- ncol(train_data$X)
  meta <- fit_centered_basis_metadata_fullrank(
    train_data$coords,
    basis_config$n_basis %||% 6L,
    basis_config$rank_tolerance %||% 1e-10,
    basis_config$full_rank_method %||% "svd"
  )
  phi <- apply_centered_basis_metadata_fullrank(train_data$coords, meta)
  center_error <- max(abs(colMeans(phi)))
  if (center_error >= 1e-8) stop("Centered training basis check failed.")
  q_basis <- cbind(constant = 1, phi)
  H <- ncol(q_basis)
  Z <- build_whole_surface_design(train_data$X, q_basis)
  W <- cbind(intercept = 1, Z)
  lambda0 <- model_config$lambda0 %||% 10
  lambda1 <- model_config$lambda1 %||% 1
  initialization <- model_config$initialization %||% list()
  gamma <- initialization$gamma %||% rep(0, p)
  tau <- initialization$tau %||% vapply(gamma, function(g) {
    lambda <- if (g > 0.5) lambda1 else lambda0
    (H + 1) / lambda^2
  }, numeric(1))
  xi <- initialization$xi %||% c(
    initialization$beta0 %||% mean(train_data$y),
    initialization$eta %||% rep(0, p * H)
  )
  if (length(xi) != ncol(W)) {
    stop("Initial xi must have length 1 + p*H for matched WS-SSGL.")
  }
  list(
    p = p,
    r = ncol(phi),
    H = H,
    q = ncol(W),
    meta = meta,
    phi = phi,
    q_basis = q_basis,
    W = W,
    center_error = center_error,
    parameters = list(
      lambda0 = lambda0,
      lambda1 = lambda1,
      a_sigma = model_config$a_sigma %||% 0.5,
      b_sigma = model_config$b_sigma %||% (var(train_data$y) / 2),
      a_gamma = model_config$a_gamma %||% 1,
      b_gamma = model_config$b_gamma %||% p
    ),
    initialization = list(
      xi = xi,
      gamma = gamma,
      tau = tau,
      sigma2 = initialization$sigma2 %||% var(train_data$y),
      pi_gamma = initialization$pi_gamma %||% 0.5
    )
  )
}

fit_matched_wsssgl_intercept_fast <- function(
    train_data, test_data, grid_data,
    basis_config = list(n_basis = 6L, full_rank_centered = TRUE,
                        full_rank_method = "svd", rank_tolerance = 1e-10),
    model_config = list(), mcmc_config = list(), seed = 1L) {
  load_matched_wsssgl_intercept_fast()
  set.seed(seed)
  started <- proc.time()[3]
  prepared <- prepare_matched_wsssgl_intercept(
    train_data, basis_config, model_config
  )
  pars <- prepared$parameters
  init <- prepared$initialization
  n_iter <- mcmc_config$n_iter %||% 400L
  burn_in <- mcmc_config$burn_in %||% 150L

  draws <- matched_wsssgl_intercept_gibbs_cpp(
    y = as.numeric(train_data$y),
    W = prepared$W,
    p = prepared$p,
    H = prepared$H,
    n_iter = n_iter,
    burn_in = burn_in,
    lambda0 = pars$lambda0,
    lambda1 = pars$lambda1,
    a_sigma = pars$a_sigma,
    b_sigma = pars$b_sigma,
    a_gamma = pars$a_gamma,
    b_gamma = pars$b_gamma,
    xi = init$xi,
    gamma = init$gamma,
    tau = init$tau,
    sigma2 = init$sigma2,
    pi_gamma = init$pi_gamma
  )

  eta_mean <- rowMeans(draws$eta)
  split <- split_eta_whole_surface(eta_mean, prepared$p, prepared$H)
  theta_draws <- draws$eta[seq(1L, by = prepared$H, length.out = prepared$p),
                           , drop = FALSE]
  q_test <- build_whole_surface_basis(test_data$coords, prepared$meta)
  q_grid <- build_whole_surface_basis(grid_data$coords, prepared$meta)
  phi_grid <- q_grid[, -1, drop = FALSE]
  beta_test <- surface_from_eta(eta_mean, q_test, prepared$p)
  beta_grid <- surface_from_eta(eta_mean, q_grid, prepared$p)
  u_grid <- matrix(0, nrow(beta_grid), prepared$p)
  for (j in seq_len(prepared$p)) {
    alpha_idx <- ((j - 1L) * prepared$r + 1L):(j * prepared$r)
    u_grid[, j] <- as.vector(phi_grid %*% split$alpha[alpha_idx])
  }
  beta0_mean <- mean(draws$beta0)
  pred_test <- beta0_mean + rowSums(test_data$X * beta_test)

  result <- standard_result(
    "matched_wsssgl_intercept",
    split$theta,
    theta_draws,
    rowMeans(draws$gamma),
    beta_grid,
    u_grid,
    pred_test,
    proc.time()[3] - started,
    diagnostics = list(
      beta0_mean = beta0_mean,
      beta0_draws = draws$beta0,
      center_error = prepared$center_error,
      rb_pip = rowMeans(draws$gamma_prob),
      gamma_draws = draws$gamma,
      gamma_prob_draws = draws$gamma_prob,
      eta_draws = draws$eta,
      eta_mean = eta_mean,
      tau_draws = draws$tau,
      sigma2_draws = draws$sigma2,
      pi_gamma_draws = draws$pi_gamma,
      pip_target = "whole coefficient surface",
      theta_u_decomposition = "post-hoc from eta=(theta, alpha)",
      design_column_order = c("beta0", "eta_by_predictor_blocks"),
      raw_basis_dimension = prepared$meta$raw_h %||% NA_integer_,
      centered_svd_rank = prepared$r,
      whole_surface_block_dimension = prepared$H,
      rank_loss_from_centering =
        prepared$meta$rank_loss_from_centering %||% NA_integer_,
      full_rank_method = prepared$meta$full_rank_method %||% NA_character_
    ),
    config = list(
      basis = prepared$meta,
      model = model_config,
      mcmc = mcmc_config,
      implementation = "matched whole-surface SSGL with explicit intercept and GD centered-SVD geometry"
    )
  )
  result$beta0_mean <- beta0_mean
  result$beta0_draws <- draws$beta0
  result$eta_mean <- eta_mean
  result$eta_draws <- draws$eta
  result$alpha_mean <- split$alpha
  result
}
