`%||%` <- function(x, y) if (is.null(x)) y else x

tensor_basis <- function(bx, by) {
  out <- matrix(0, nrow(bx), ncol(bx) * ncol(by))
  for (i in seq_len(ncol(bx))) {
    for (j in seq_len(ncol(by))) {
      out[, (i - 1L) * ncol(by) + j] <- bx[, i] * by[, j]
    }
  }
  out
}

fit_basis_metadata <- function(coords, n_basis = 4L, centered = TRUE) {
  coords <- as.matrix(coords)
  boundary_x <- if (min(coords[, 1]) >= 0 && max(coords[, 1]) <= 1) {
    c(0, 1)
  } else {
    range(coords[, 1])
  }
  boundary_y <- if (min(coords[, 2]) >= 0 && max(coords[, 2]) <= 1) {
    c(0, 1)
  } else {
    range(coords[, 2])
  }
  bx <- splines::bs(coords[, 1], df = n_basis, intercept = TRUE,
                    Boundary.knots = boundary_x)
  by <- splines::bs(coords[, 2], df = n_basis, intercept = TRUE,
                    Boundary.knots = boundary_y)
  phi <- tensor_basis(bx, by)
  list(
    n_basis = n_basis,
    knots_x = attr(bx, "knots"),
    knots_y = attr(by, "knots"),
    boundary_x = boundary_x,
    boundary_y = boundary_y,
    train_mean = if (centered) colMeans(phi) else rep(0, ncol(phi)),
    centered = centered
  )
}

apply_basis_metadata <- function(coords, metadata) {
  coords <- as.matrix(coords)
  bx <- splines::bs(coords[, 1], knots = metadata$knots_x,
                    Boundary.knots = metadata$boundary_x, intercept = TRUE)
  by <- splines::bs(coords[, 2], knots = metadata$knots_y,
                    Boundary.knots = metadata$boundary_y, intercept = TRUE)
  phi <- tensor_basis(bx, by)
  sweep(phi, 2, metadata$train_mean, "-")
}

build_spatial_design <- function(X, phi) {
  X <- as.matrix(X)
  z <- matrix(0, nrow(X), ncol(X) * ncol(phi))
  for (j in seq_len(ncol(X))) {
    idx <- ((j - 1L) * ncol(phi) + 1L):(j * ncol(phi))
    z[, idx] <- phi * X[, j]
  }
  z
}

surface_from_components <- function(theta, alpha, phi, p) {
  h <- ncol(phi)
  beta <- matrix(0, nrow(phi), p)
  for (j in seq_len(p)) {
    idx <- ((j - 1L) * h + 1L):(j * h)
    beta[, j] <- theta[j] + as.vector(phi %*% alpha[idx])
  }
  beta
}

standard_result <- function(method_name, theta_mean = NULL, theta_draws = NULL,
                            pip = NULL, beta_hat_grid, u_hat_grid = NULL,
                            pred_test, runtime, diagnostics = list(),
                            config = list()) {
  list(
    method_name = method_name,
    theta_mean = theta_mean,
    theta_draws = theta_draws,
    pip = pip,
    beta_hat_grid = beta_hat_grid,
    u_hat_grid = u_hat_grid,
    pred_test = as.numeric(pred_test),
    runtime = as.numeric(runtime),
    diagnostics = diagnostics,
    config = config
  )
}

simulate_four_predictor_split <- function(n_train = 200L, n_test = 100L,
                                          grid_size = 20L, sigma = 0.5,
                                          theta_value = 1,
                                          u_target_norm = 1,
                                          seed = 20260618L) {
  set.seed(seed)
  n <- n_train + n_test
  coords <- cbind(runif(n), runif(n))
  X_raw <- matrix(rnorm(n * 4L), n, 4L)
  train_id <- seq_len(n_train)
  test_id <- n_train + seq_len(n_test)
  x_mean <- colMeans(X_raw[train_id, , drop = FALSE])
  x_sd <- apply(X_raw[train_id, , drop = FALSE], 2, sd)
  X <- sweep(sweep(X_raw, 2, x_mean, "-"), 2, x_sd, "/")

  bump_raw <- function(s) {
    exp(-((s[, 1] - 0.35)^2 + (s[, 2] - 0.65)^2) / 0.04)
  }
  bump_train <- bump_raw(coords[train_id, , drop = FALSE])
  bump_mean <- mean(bump_train)
  bump_norm <- sqrt(mean((bump_train - bump_mean)^2))
  bump <- function(s) u_target_norm * (bump_raw(s) - bump_mean) / bump_norm

  theta <- c(theta_value, 0, theta_value, 0)
  make_truth <- function(s) {
    u <- bump(s)
    u_mat <- cbind(0, u, u, 0)
    list(
      theta = theta,
      u = u_mat,
      beta = sweep(u_mat, 2, theta, "+")
    )
  }
  truth_all <- make_truth(coords)
  mu <- rowSums(X * truth_all$beta)
  y <- mu + rnorm(n, sd = sigma)
  grid_coords <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = grid_size),
    y = seq(0, 1, length.out = grid_size)
  ))
  truth_grid <- make_truth(grid_coords)

  list(
    train = list(y = y[train_id], X = X[train_id, , drop = FALSE],
                 coords = coords[train_id, , drop = FALSE]),
    test = list(y = y[test_id], X = X[test_id, , drop = FALSE],
                coords = coords[test_id, , drop = FALSE]),
    grid = list(coords = grid_coords, true_beta = truth_grid$beta,
                true_u = truth_grid$u),
    truth = list(theta = theta,
                 spatial_deviation = c(FALSE, TRUE, TRUE, FALSE),
                 bump_train_mean = bump_mean,
                 bump_train_norm = bump_norm),
    preprocessing = list(x_mean = x_mean, x_sd = x_sd),
    seed = seed
  )
}

simulate_p10_pilot_split <- function(n_train = 700L, n_test = 250L,
                                     grid_size = 50L, sigma = 0.35,
                                     u_target_norm = 1,
                                     seed = 20260619L) {
  set.seed(seed)
  n <- n_train + n_test
  train_id <- seq_len(n_train)
  test_id <- n_train + seq_len(n_test)
  coords <- cbind(runif(n), runif(n))
  X_raw <- matrix(rnorm(n * 10L), n, 10L)
  x_mean <- colMeans(X_raw[train_id, , drop = FALSE])
  x_sd <- apply(X_raw[train_id, , drop = FALSE], 2, sd)
  X <- sweep(sweep(X_raw, 2, x_mean, "-"), 2, x_sd, "/")

  bump_raw <- function(s) {
    exp(-((s[, 1] - 0.35)^2 + (s[, 2] - 0.65)^2) / 0.04)
  }
  raw_train <- bump_raw(coords[train_id, , drop = FALSE])
  bump_mean <- mean(raw_train)
  bump_norm <- sqrt(mean((raw_train - bump_mean)^2))
  bump <- function(s) {
    u_target_norm * (bump_raw(s) - bump_mean) / bump_norm
  }

  theta <- c(1, -1, 0, 0, 1, -1, rep(0, 4))
  signs <- c(0, 0, 1, -1, 1, -1, rep(0, 4))
  make_truth <- function(s) {
    g <- bump(s)
    u <- outer(g, signs)
    list(
      theta = theta,
      u = u,
      beta = sweep(u, 2, theta, "+")
    )
  }
  truth_all <- make_truth(coords)
  mu <- rowSums(X * truth_all$beta)
  y <- mu + rnorm(n, sd = sigma)
  grid_coords <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = grid_size),
    y = seq(0, 1, length.out = grid_size)
  ))
  truth_grid <- make_truth(grid_coords)

  list(
    train = list(
      y = y[train_id], X = X[train_id, , drop = FALSE],
      coords = coords[train_id, , drop = FALSE]
    ),
    test = list(
      y = y[test_id], X = X[test_id, , drop = FALSE],
      coords = coords[test_id, , drop = FALSE]
    ),
    grid = list(
      coords = grid_coords, true_beta = truth_grid$beta,
      true_u = truth_grid$u, grid_size = grid_size
    ),
    truth = list(
      theta = theta,
      spatial_deviation = signs != 0,
      bump_train_mean = bump_mean,
      bump_train_norm = bump_norm,
      u_target_norm = u_target_norm,
      active_surface = seq_len(10) <= 6
    ),
    preprocessing = list(x_mean = x_mean, x_sd = x_sd),
    seed = seed
  )
}

compute_debug_metrics <- function(fit, data) {
  beta_mse <- colMeans((fit$beta_hat_grid - data$grid$true_beta)^2)
  out <- data.frame(
    method = fit$method_name,
    finished = TRUE,
    mspe = mean((data$test$y - fit$pred_test)^2),
    theta_mse = if (is.null(fit$theta_mean)) NA_real_
      else mean((fit$theta_mean - data$truth$theta)^2),
    beta_mse_x1 = beta_mse[1],
    beta_mse_x2 = beta_mse[2],
    beta_mse_x3 = beta_mse[3],
    beta_mse_x4 = beta_mse[4],
    runtime_seconds = fit$runtime,
    pip = if (is.null(fit$pip)) NA_character_
      else paste(round(fit$pip, 3), collapse = ","),
    stringsAsFactors = FALSE
  )
  out
}
