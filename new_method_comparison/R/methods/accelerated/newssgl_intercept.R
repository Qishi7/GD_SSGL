find_accelerated_intercept_root <- function() {
  candidates <- c(
    file.path(getwd(), "new_method_comparison"),
    getwd()
  )
  hit <- candidates[file.exists(file.path(
    candidates, "R", "methods", "accelerated",
    "newssgl_intercept.cpp"
  ))]
  if (!length(hit)) stop("Cannot locate new_method_comparison.")
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

.accelerated_intercept_root <- find_accelerated_intercept_root()
source(file.path(.accelerated_intercept_root, "R", "methods", "ssgl.R"))

load_newssgl_intercept_fast <- function(rebuild = FALSE) {
  if (exists("newssgl_intercept_gibbs_cpp", mode = "function") && !rebuild) {
    return(invisible(TRUE))
  }
  makevars <- file.path(.accelerated_intercept_root, "configs", "Makevars.win")
  if (.Platform$OS.type == "windows" && file.exists(makevars)) {
    Sys.setenv(R_MAKEVARS_USER = makevars)
  }
  Rcpp::sourceCpp(file.path(
    .accelerated_intercept_root, "R", "methods", "accelerated",
    "newssgl_intercept.cpp"
  ), rebuild = rebuild)
  invisible(TRUE)
}

fit_centered_basis_metadata_fullrank <- function(coords, n_basis = 4L,
                                                 tolerance = 1e-10,
                                                 method = c("svd", "qr")) {
  method <- match.arg(method)
  meta <- fit_basis_metadata(coords, n_basis, centered = TRUE)
  phi <- apply_basis_metadata(coords, meta)
  if (method == "svd") {
    sv <- svd(phi)
    rank <- sum(sv$d > tolerance * max(sv$d))
    selected <- seq_len(rank)
    meta$svd_v <- sv$v[, seq_len(rank), drop = FALSE]
    meta$svd_d <- sv$d[seq_len(rank)]
  } else {
    qr_phi <- qr(phi, tol = tolerance)
    rank <- qr_phi$rank
    selected <- sort(qr_phi$pivot[seq_len(rank)])
  }
  meta$raw_h <- ncol(phi)
  meta$raw_rank <- qr(apply_basis_metadata(coords, fit_basis_metadata(
    coords, n_basis, centered = FALSE
  )), tol = tolerance)$rank
  meta$centered_rank <- rank
  meta$basis_rank_tolerance <- tolerance
  meta$full_rank_centered <- TRUE
  meta$full_rank_method <- method
  meta$selected_columns <- selected
  meta$rank_loss_from_centering <- ncol(phi) - rank
  meta
}

apply_centered_basis_metadata_fullrank <- function(coords, metadata) {
  phi <- apply_basis_metadata(coords, metadata)
  if (isTRUE(metadata$full_rank_centered)) {
    if (identical(metadata$full_rank_method, "svd")) {
      phi <- phi %*% metadata$svd_v
    } else {
      phi <- phi[, metadata$selected_columns, drop = FALSE]
    }
  }
  phi
}

coerce_intercept_alpha_init <- function(alpha_init, p, h, metadata) {
  if (is.null(alpha_init)) return(rep(0, p * h))
  raw_h <- metadata$raw_h %||% h
  if (length(alpha_init) == p * h) return(alpha_init)
  if (isTRUE(metadata$full_rank_centered) && length(alpha_init) == p * raw_h) {
    out <- numeric(p * h)
    for (j in seq_len(p)) {
      raw_idx <- ((j - 1L) * raw_h + 1L):(j * raw_h)
      new_idx <- ((j - 1L) * h + 1L):(j * h)
      if (identical(metadata$full_rank_method, "svd")) {
        out[new_idx] <- as.numeric(t(metadata$svd_v) %*% alpha_init[raw_idx])
      } else {
        out[new_idx] <- alpha_init[raw_idx][metadata$selected_columns]
      }
    }
    return(out)
  }
  stop("Initial alpha must have length p*h, or p*raw_h for full-rank basis.")
}

prepare_newssgl_intercept <- function(train_data,
                                      basis_config = list(),
                                      model_config = list()) {
  p <- ncol(train_data$X)
  full_rank_centered <- basis_config$full_rank_centered %||% TRUE
  meta <- if (isTRUE(full_rank_centered)) {
    fit_centered_basis_metadata_fullrank(
      train_data$coords,
      basis_config$n_basis %||% 4L,
      basis_config$rank_tolerance %||% 1e-10,
      basis_config$full_rank_method %||% "svd"
    )
  } else {
    fit_basis_metadata(
      train_data$coords,
      basis_config$n_basis %||% 4L,
      centered = TRUE
    )
  }
  phi <- if (isTRUE(full_rank_centered)) {
    apply_centered_basis_metadata_fullrank(train_data$coords, meta)
  } else {
    apply_basis_metadata(train_data$coords, meta)
  }
  center_error <- max(abs(colMeans(phi)))
  if (center_error >= 1e-8) stop("Centered training basis check failed.")
  Z <- build_spatial_design(train_data$X, phi)
  W <- cbind(intercept = 1, train_data$X, Z)
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

  init_xi <- initialization$xi
  if (is.null(init_xi)) {
    init_alpha <- coerce_intercept_alpha_init(
      initialization$alpha, p, h, meta
    )
    init_xi <- c(
      initialization$beta0 %||% mean(train_data$y),
      initialization$theta %||% rep(0, p),
      init_alpha
    )
  }
  if (length(init_xi) != q) {
    stop("Initial xi must have length 1 + p + p*h for intercept model.")
  }

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
      xi = init_xi,
      omega = initialization$omega %||% rep(1, p),
      gamma = gamma,
      tau = tau,
      sigma2 = initialization$sigma2 %||% var(train_data$y),
      lambda_theta2 = initialization$lambda_theta2 %||% 1,
      pi_gamma = initialization$pi_gamma %||% 0.5
    )
  )
}

fit_newssgl_intercept_fast <- function(
    train_data, test_data, grid_data,
    basis_config = list(), model_config = list(),
    mcmc_config = list(), seed = 1L) {
  load_newssgl_intercept_fast()
  set.seed(seed)
  started <- proc.time()[3]
  prepared <- prepare_newssgl_intercept(
    train_data, basis_config, model_config
  )
  pars <- prepared$parameters
  init <- prepared$initialization
  n_iter <- mcmc_config$n_iter %||% 400L
  burn_in <- mcmc_config$burn_in %||% 150L

  draws <- newssgl_intercept_gibbs_cpp(
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

  beta0_mean <- mean(draws$beta0)
  theta_mean <- rowMeans(draws$theta)
  alpha_mean <- rowMeans(draws$alpha)
  phi_test <- apply_centered_basis_metadata_fullrank(
    test_data$coords, prepared$meta
  )
  phi_grid <- apply_centered_basis_metadata_fullrank(
    grid_data$coords, prepared$meta
  )
  beta_test <- surface_from_components(
    theta_mean, alpha_mean, phi_test, prepared$p
  )
  beta_grid <- surface_from_components(
    theta_mean, alpha_mean, phi_grid, prepared$p
  )
  pred_test <- beta0_mean + rowSums(test_data$X * beta_test)

  result <- standard_result(
    "newssgl_intercept",
    theta_mean,
    draws$theta,
    rowMeans(draws$gamma),
    beta_grid,
    sweep(beta_grid, 2, theta_mean, "-"),
    pred_test,
    proc.time()[3] - started,
    diagnostics = list(
      beta0_mean = beta0_mean,
      beta0_draws = draws$beta0,
      center_error = prepared$center_error,
      rb_pip = rowMeans(draws$gamma_prob),
      gamma_draws = draws$gamma,
      gamma_prob_draws = draws$gamma_prob,
      alpha_draws = draws$alpha,
      sigma2_draws = draws$sigma2,
      omega_draws = draws$omega,
      tau_draws = draws$tau,
      lambda_theta2_draws = draws$lambda_theta2,
      pi_gamma_draws = draws$pi_gamma,
      design_column_order = c("beta0", paste0("theta", seq_len(prepared$p)),
                              "alpha_by_predictor_blocks"),
      raw_basis_dimension = prepared$meta$raw_h %||% prepared$h,
      effective_basis_dimension = prepared$h,
      centered_basis_rank = prepared$meta$centered_rank %||% NA_integer_,
      rank_loss_from_centering =
        prepared$meta$rank_loss_from_centering %||% NA_integer_,
      full_rank_method = prepared$meta$full_rank_method %||% NA_character_,
      selected_basis_columns = prepared$meta$selected_columns %||% seq_len(prepared$h)
    ),
    config = list(
      basis = prepared$meta,
      model = model_config,
      mcmc = mcmc_config,
      implementation = "RcppArmadillo accelerated with explicit intercept"
    )
  )
  result$beta0_mean <- beta0_mean
  result$beta0_draws <- draws$beta0
  result$alpha_mean <- alpha_mean
  result
}

fit_newssgl_intercept_collapsed_fast <- function(
    train_data, test_data, grid_data,
    basis_config = list(), model_config = list(),
    mcmc_config = list(), seed = 1L) {
  load_newssgl_intercept_fast()
  set.seed(seed)
  started <- proc.time()[3]
  prepared <- prepare_newssgl_intercept(
    train_data, basis_config, model_config
  )
  pars <- prepared$parameters
  init <- prepared$initialization
  n_iter <- mcmc_config$n_iter %||% 400L
  burn_in <- mcmc_config$burn_in %||% 150L
  W_nonintercept <- prepared$W[, -1, drop = FALSE]
  xi_nonintercept <- init$xi[-1]

  draws <- newssgl_intercept_collapsed_gibbs_cpp(
    y = as.numeric(train_data$y),
    W = W_nonintercept,
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
    xi = xi_nonintercept,
    omega = init$omega,
    gamma = init$gamma,
    tau = init$tau,
    sigma2 = init$sigma2,
    lambda_theta2 = init$lambda_theta2,
    pi_gamma = init$pi_gamma
  )

  beta0_mean <- mean(draws$beta0_cond_mean)
  theta_mean <- rowMeans(draws$theta)
  alpha_mean <- rowMeans(draws$alpha)
  phi_test <- apply_centered_basis_metadata_fullrank(
    test_data$coords, prepared$meta
  )
  phi_grid <- apply_centered_basis_metadata_fullrank(
    grid_data$coords, prepared$meta
  )
  beta_test <- surface_from_components(
    theta_mean, alpha_mean, phi_test, prepared$p
  )
  beta_grid <- surface_from_components(
    theta_mean, alpha_mean, phi_grid, prepared$p
  )
  pred_test <- beta0_mean + rowSums(test_data$X * beta_test)

  result <- standard_result(
    "newssgl_intercept_collapsed",
    theta_mean,
    draws$theta,
    rowMeans(draws$gamma),
    beta_grid,
    sweep(beta_grid, 2, theta_mean, "-"),
    pred_test,
    proc.time()[3] - started,
    diagnostics = list(
      beta0_mean = beta0_mean,
      beta0_cond_mean_draws = draws$beta0_cond_mean,
      beta0_sampled = FALSE,
      center_error = prepared$center_error,
      rb_pip = rowMeans(draws$gamma_prob),
      gamma_draws = draws$gamma,
      gamma_prob_draws = draws$gamma_prob,
      alpha_draws = draws$alpha,
      sigma2_draws = draws$sigma2,
      omega_draws = draws$omega,
      tau_draws = draws$tau,
      lambda_theta2_draws = draws$lambda_theta2,
      pi_gamma_draws = draws$pi_gamma,
      design_column_order = c(paste0("theta", seq_len(prepared$p)),
                              "alpha_by_predictor_blocks"),
      intercept_handling = "beta0 integrated out; prediction uses posterior mean conditional beta0"
      ,
      raw_basis_dimension = prepared$meta$raw_h %||% prepared$h,
      effective_basis_dimension = prepared$h,
      centered_basis_rank = prepared$meta$centered_rank %||% NA_integer_,
      rank_loss_from_centering =
        prepared$meta$rank_loss_from_centering %||% NA_integer_,
      full_rank_method = prepared$meta$full_rank_method %||% NA_character_,
      selected_basis_columns = prepared$meta$selected_columns %||% seq_len(prepared$h)
    ),
    config = list(
      basis = prepared$meta,
      model = model_config,
      mcmc = mcmc_config,
      implementation = "RcppArmadillo accelerated with collapsed explicit intercept"
    )
  )
  result$beta0_mean <- beta0_mean
  result$beta0_cond_mean_draws <- draws$beta0_cond_mean
  result$alpha_mean <- alpha_mean
  result
}
