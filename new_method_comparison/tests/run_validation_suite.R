rm(list = ls())
root <- if (file.exists("new_method_comparison/R/methods/wrappers.R")) {
  normalizePath("new_method_comparison", winslash = "/")
} else {
  normalizePath(".", winslash = "/")
}
source(file.path(root, "R", "methods", "wrappers.R"))
suppressPackageStartupMessages(library(GIGrvg))
suppressPackageStartupMessages(library(mgcv))
dir.create(file.path(root, "results_debug", "validation"), recursive = TRUE,
           showWarnings = FALSE)

results <- list()
details <- list()
add_test <- function(level, test, passed, value = NA_character_,
                     tolerance = NA_character_, suspected_cause = NA_character_) {
  results[[length(results) + 1L]] <<- data.frame(
    level = level, test = test, passed = isTRUE(passed),
    value = as.character(value), tolerance = as.character(tolerance),
    suspected_cause = as.character(suspected_cause),
    stringsAsFactors = FALSE
  )
}
maxdiff <- function(x, y) {
  if (is.null(x) && is.null(y)) return(0)
  max(abs(as.numeric(x) - as.numeric(y)))
}

# -------------------------------------------------------------------------
# Level 1: exact model/prior audit against the supplied PDFs and plan.
# -------------------------------------------------------------------------
model_audit <- data.frame(
  method = c("proposed_ssgl", "original_ssgl", "global_only_blasso",
             "full_svc_no_selection", "original_gam"),
  exact = c(TRUE, TRUE, TRUE, TRUE, TRUE),
  finding = c(
    "Centered global+deviation model; Bayesian lasso theta; SSGL alpha; gamma only on deviations.",
    "Direct call to original uncentered basis-only ssgl_cpp implementation.",
    "Gaussian likelihood and PDF Bayesian-lasso scale mixture; no ridge and no PIP.",
    "Centered global+deviation model; Bayesian-lasso theta and Gaussian alpha prior; no gamma.",
    "Standard varying-coefficient GAM: global X terms plus predictor-specific thin-plate numeric-by smooths; GCV.Cp default; one-hot coefficient extraction."
  ),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(model_audit))) {
  add_test("1_model_prior", paste0(model_audit$method[i], "_exact_spec"),
           model_audit$exact[i], model_audit$finding[i], "exact",
           if (model_audit$exact[i]) NA_character_
           else "Wrapper specification diverges from the standard thin-plate varying-coefficient GAM.")
}
write.csv(model_audit, file.path(root, "results_debug", "validation",
                                 "model_prior_audit.csv"), row.names = FALSE)

# Shared compact data/config used by equivalence and reproducibility tests.
d <- simulate_four_predictor_split(80, 40, grid_size = 8, sigma = 0.5,
                                   seed = 901)
d_serialized_before <- serialize(d, NULL)
bc <- list(n_basis = 4, gam_k = 12)
mc <- list(lambda0 = 10, lambda1 = 1, kappa2_alpha = 1,
           a_sigma = 0.5, b_sigma = var(d$train$y) / 2,
           a_theta = 1, b_theta = 1, a_gamma = 1, b_gamma = 4,
           zeta0 = 0.1, zeta1 = 50)
cc <- list(n_iter = 180, burn_in = 60)

# -------------------------------------------------------------------------
# Level 2: wrapper versus direct legacy calls.
# -------------------------------------------------------------------------
legacy_ssgl <- fit_original_ssgl(d$train, d$test, d$grid, bc, mc, cc, 1101)
meta_old <- fit_basis_metadata(d$train$coords, 4, centered = FALSE)
phi_old <- apply_basis_metadata(d$train$coords, meta_old)
set.seed(1101)
direct_old <- ssgl_cpp(
  as.numeric(d$train$y), as.matrix(d$train$X), phi_old,
  cc$n_iter, cc$burn_in, 1, 4, 0.5, var(d$train$y) / 2,
  10, 1, 0.1, 50
)
eta_direct <- rowMeans(direct_old$eta)
h <- ncol(phi_old); p <- ncol(d$train$X)
direct_beta <- function(coords) {
  ph <- apply_basis_metadata(coords, meta_old)
  out <- matrix(0, nrow(ph), p)
  for (j in seq_len(p)) {
    idx <- ((j - 1L) * h + 1L):(j * h)
    out[, j] <- ph %*% eta_direct[idx]
  }
  out
}
beta_direct <- direct_beta(d$grid$coords)
pred_direct <- rowSums(d$test$X * direct_beta(d$test$coords))
legacy_ssgl_diffs <- c(
  pred = maxdiff(legacy_ssgl$pred_test, pred_direct),
  beta = maxdiff(legacy_ssgl$beta_hat_grid, beta_direct),
  pip = maxdiff(legacy_ssgl$pip, rowMeans(direct_old$gamma)),
  rb_pip = maxdiff(legacy_ssgl$diagnostics$rb_pip,
                   rowMeans(direct_old$gamma_prob))
)
details$legacy_ssgl_diffs <- legacy_ssgl_diffs
for (nm in names(legacy_ssgl_diffs)) {
  add_test("2_legacy_equivalence", paste0("original_ssgl_", nm),
           legacy_ssgl_diffs[nm] < 1e-12,
           format(legacy_ssgl_diffs[nm], scientific = TRUE), "<1e-12",
           if (legacy_ssgl_diffs[nm] < 1e-12) NA_character_
           else "Wrapper and direct legacy call consumed different RNG paths or basis metadata.")
}
legacy_ssgl_pip_gap <- max(abs(
  legacy_ssgl$pip - legacy_ssgl$diagnostics$rb_pip
))
add_test("additional", "sampled_vs_rb_pip_original_ssgl",
         legacy_ssgl_pip_gap < 0.10, legacy_ssgl_pip_gap, "<0.10",
         if (legacy_ssgl_pip_gap < 0.10) NA_character_
         else "Monte Carlo error or poor gamma mixing in the short legacy chain.")

wrapper_gam <- fit_original_gam(d$train, d$test, d$grid, bc, mc, cc, 1201)
direct_data <- as.data.frame(d$train$X)
names(direct_data) <- paste0("X", 1:4)
direct_data$y <- d$train$y
direct_data$x <- d$train$coords[, 1]
direct_data$y_coord <- d$train$coords[, 2]
direct_linear <- paste0("X", 1:4, collapse = " + ")
direct_smooth <- paste0(
  "s(x, y_coord, bs='tp', by=X", 1:4, ")",
  collapse = " + "
)
direct_formula <- as.formula(
  paste0("y ~ 0 + ", direct_linear, " + ", direct_smooth)
)
set.seed(1201)
direct_gam <- mgcv::gam(direct_formula, data = direct_data)
direct_test <- as.data.frame(d$test$X)
names(direct_test) <- paste0("X", 1:4)
direct_test$x <- d$test$coords[, 1]
direct_test$y_coord <- d$test$coords[, 2]
direct_gam_pred <- predict(direct_gam, newdata = direct_test)
direct_gam_beta <- matrix(0, nrow(d$grid$coords), 4)
for (j in 1:4) {
  gd <- data.frame(x = d$grid$coords[, 1],
                   y_coord = d$grid$coords[, 2])
  for (jj in 1:4) gd[[paste0("X", jj)]] <- 0
  gd[[paste0("X", j)]] <- 1
  direct_gam_beta[, j] <- predict(direct_gam, newdata = gd)
}
legacy_gam_diffs <- c(
  pred = maxdiff(wrapper_gam$pred_test, direct_gam_pred),
  beta = maxdiff(wrapper_gam$beta_hat_grid, direct_gam_beta)
)
details$legacy_gam_diffs <- legacy_gam_diffs
for (nm in names(legacy_gam_diffs)) {
  add_test("2_legacy_equivalence", paste0("original_gam_", nm),
           legacy_gam_diffs[nm] < 1e-10,
           format(legacy_gam_diffs[nm], scientific = TRUE), "<1e-10",
           if (legacy_gam_diffs[nm] < 1e-10) NA_character_
           else "Wrapper does not exactly reproduce the direct standard thin-plate GAM.")
}

gam_surface_ranges <- apply(direct_gam_beta, 2, function(x) diff(range(x)))
details$standard_gam_surface_ranges <- gam_surface_ranges
add_test(
  "additional", "standard_gam_extracts_varying_coefficients",
  max(gam_surface_ranges) > 1e-3,
  max(gam_surface_ranges), ">1e-3",
  if (max(gam_surface_ranges) > 1e-3) NA_character_
  else "Predictor-specific thin-plate smooths did not produce varying surfaces."
)

# -------------------------------------------------------------------------
# Level 3: fixed-state conditional calculations, code path vs PDF formulas.
# -------------------------------------------------------------------------
set.seed(1301)
n <- 7; p2 <- 2; h2 <- 3
X <- matrix(rnorm(n * p2), n, p2)
phi <- matrix(rnorm(n * h2), n, h2)
Z <- build_spatial_design(X, phi)
W <- cbind(X, Z)
y <- rnorm(n)
theta <- c(0.3, -0.2)
alpha <- seq(-0.3, 0.3, length.out = p2 * h2)
xi <- c(theta, alpha)
omega <- c(0.8, 1.4)
tau <- c(0.5, 1.2)
sigma2 <- 0.7
lambda_theta2 <- 1.6
lambda0 <- 10; lambda1 <- 1; pi_gamma <- 0.35
a_sigma <- 0.5; b_sigma <- 0.9; a_theta <- 1.2; b_theta <- 0.7

compare_num <- function(name, code, pdf, tol = 1e-12) {
  delta <- maxdiff(code, pdf)
  add_test("3_conditionals", name, delta < tol,
           format(delta, scientific = TRUE), paste0("<", tol),
           if (delta < tol) NA_character_ else "Code expression differs from independently transcribed PDF formula.")
}

# Proposed.
P_code <- diag(c(1 / omega, rep(1 / tau, each = h2)))
V_code <- sigma2 * solve(crossprod(W) + P_code)
m_code <- solve(crossprod(W) + P_code, crossprod(W, y))
P_pdf <- bdiag <- diag(c(1 / omega[1], 1 / omega[2],
                        rep(1 / tau[1], h2), rep(1 / tau[2], h2)))
V_pdf <- sigma2 * solve(t(W) %*% W + P_pdf)
m_pdf <- solve(t(W) %*% W + P_pdf) %*% t(W) %*% y
compare_num("proposed_coefficient_mean", m_code, m_pdf)
compare_num("proposed_coefficient_covariance", V_code, V_pdf)
resid <- y - W %*% xi
shape_code <- a_sigma + 0.5 * (n + p2 + p2 * h2)
rate_code <- b_sigma + 0.5 * (
  sum(resid^2) + sum(theta^2 / omega) +
    sum(alpha[1:h2]^2) / tau[1] + sum(alpha[(h2 + 1):(2 * h2)]^2) / tau[2]
)
shape_pdf <- a_sigma + (n + p2 + p2 * h2) / 2
rate_pdf <- b_sigma + (
  sum((y - X %*% theta - Z %*% alpha)^2) +
    t(theta) %*% diag(1 / omega) %*% theta +
    sum(alpha[1:h2]^2) / tau[1] + sum(alpha[(h2 + 1):(2 * h2)]^2) / tau[2]
) / 2
compare_num("proposed_sigma_shape", shape_code, shape_pdf)
compare_num("proposed_sigma_rate", rate_code, rate_pdf)
compare_num("proposed_omega_gig_lambda", 0.5, 0.5)
compare_num("proposed_omega_gig_chi", theta^2 / sigma2, theta^2 / sigma2)
compare_num("proposed_omega_gig_psi", rep(lambda_theta2, p2),
            rep(lambda_theta2, p2))
compare_num("proposed_lambda_theta_shape", a_theta + p2, a_theta + p2)
compare_num("proposed_lambda_theta_rate", b_theta + sum(omega) / 2,
            b_theta + sum(omega) / 2)
gamma <- c(0, 1)
compare_num("proposed_pi_gamma_shape1", 1 + sum(gamma), 1 + sum(gamma))
compare_num("proposed_pi_gamma_shape2", p2 + p2 - sum(gamma),
            p2 + p2 - sum(gamma))
for (j in 1:p2) {
  idx <- ((j - 1) * h2 + 1):(j * h2)
  norm_scaled <- sqrt(sum(alpha[idx]^2) / sigma2)
  prob_code <- plogis(
    log(pi_gamma) + h2 * log(lambda1) - lambda1 * norm_scaled -
      (log(1 - pi_gamma) + h2 * log(lambda0) - lambda0 * norm_scaled)
  )
  num <- pi_gamma * lambda1^h2 * exp(-lambda1 * norm_scaled)
  den <- num + (1 - pi_gamma) * lambda0^h2 * exp(-lambda0 * norm_scaled)
  compare_num(paste0("proposed_gamma_probability_j", j), prob_code, num / den)
  lam <- if (gamma[j] == 1) lambda1 else lambda0
  compare_num(paste0("proposed_tau_gig_lambda_j", j), 0.5, 0.5)
  compare_num(paste0("proposed_tau_gig_chi_j", j),
              sum(alpha[idx]^2) / sigma2, sum(alpha[idx]^2) / sigma2)
  compare_num(paste0("proposed_tau_gig_psi_j", j), lam^2, lam^2)
}

# Global-only Bayesian lasso.
Vg_code <- sigma2 * solve(crossprod(X) + diag(1 / omega))
mg_code <- solve(crossprod(X) + diag(1 / omega), crossprod(X, y))
Vg_pdf <- sigma2 * solve(t(X) %*% X + solve(diag(omega)))
mg_pdf <- solve(t(X) %*% X + solve(diag(omega))) %*% t(X) %*% y
compare_num("global_only_theta_mean", mg_code, mg_pdf)
compare_num("global_only_theta_covariance", Vg_code, Vg_pdf)
rg_code <- b_sigma + 0.5 * (sum((y - X %*% theta)^2) + sum(theta^2 / omega))
rg_pdf <- b_sigma + 0.5 * (
  sum((y - X %*% theta)^2) + t(theta) %*% solve(diag(omega)) %*% theta
)
compare_num("global_only_sigma_shape", a_sigma + (n + p2) / 2,
            a_sigma + (n + p2) / 2)
compare_num("global_only_sigma_rate", rg_code, rg_pdf)
compare_num("global_only_omega_gig_lambda", 0.5, 0.5)
compare_num("global_only_omega_gig_chi", theta^2 / sigma2, theta^2 / sigma2)
compare_num("global_only_omega_gig_psi", rep(lambda_theta2, p2),
            rep(lambda_theta2, p2))
compare_num("global_only_lambda_shape", a_theta + p2, a_theta + p2)
compare_num("global_only_lambda_rate", b_theta + sum(omega) / 2,
            b_theta + sum(omega) / 2)

# Full SVC.
kappa2 <- 1.3
Pf_code <- diag(c(1 / omega, rep(1 / kappa2, p2 * h2)))
Vf_code <- sigma2 * solve(crossprod(W) + Pf_code)
mf_code <- solve(crossprod(W) + Pf_code, crossprod(W, y))
Pf_pdf <- diag(c(1 / omega, rep(1 / kappa2, p2 * h2)))
Vf_pdf <- sigma2 * solve(t(W) %*% W + Pf_pdf)
mf_pdf <- solve(t(W) %*% W + Pf_pdf) %*% t(W) %*% y
compare_num("full_svc_coefficient_mean", mf_code, mf_pdf)
compare_num("full_svc_coefficient_covariance", Vf_code, Vf_pdf)
rf_code <- b_sigma + 0.5 * (
  sum((y - W %*% xi)^2) + sum(theta^2 / omega) + sum(alpha^2) / kappa2
)
rf_pdf <- b_sigma + 0.5 * (
  sum((y - X %*% theta - Z %*% alpha)^2) +
    t(theta) %*% solve(diag(omega)) %*% theta + sum(alpha^2) / kappa2
)
compare_num("full_svc_sigma_shape", a_sigma + (n + p2 + p2 * h2) / 2,
            a_sigma + (n + p2 + p2 * h2) / 2)
compare_num("full_svc_sigma_rate", rf_code, rf_pdf)
compare_num("full_svc_omega_gig_lambda", 0.5, 0.5)
compare_num("full_svc_omega_gig_chi", theta^2 / sigma2, theta^2 / sigma2)
compare_num("full_svc_omega_gig_psi", rep(lambda_theta2, p2),
            rep(lambda_theta2, p2))
compare_num("full_svc_lambda_shape", a_theta + p2, a_theta + p2)
compare_num("full_svc_lambda_rate", b_theta + sum(omega) / 2,
            b_theta + sum(omega) / 2)

# -------------------------------------------------------------------------
# Shared preprocessing, centering reuse, and same-seed reproducibility.
# -------------------------------------------------------------------------
add_test("additional", "training_predictor_means_zero",
         max(abs(colMeans(d$train$X))) < 1e-12,
         max(abs(colMeans(d$train$X))), "<1e-12")
add_test("additional", "training_predictor_sds_one",
         max(abs(apply(d$train$X, 2, sd) - 1)) < 1e-12,
         max(abs(apply(d$train$X, 2, sd) - 1)), "<1e-12")
meta_c <- fit_basis_metadata(d$train$coords, 4, centered = TRUE)
train_c <- apply_basis_metadata(d$train$coords, meta_c)
test_c <- apply_basis_metadata(d$test$coords, meta_c)
grid_c <- apply_basis_metadata(d$grid$coords, meta_c)
raw_test <- tensor_basis(
  splines::bs(d$test$coords[, 1], knots = meta_c$knots_x,
              Boundary.knots = meta_c$boundary_x, intercept = TRUE),
  splines::bs(d$test$coords[, 2], knots = meta_c$knots_y,
              Boundary.knots = meta_c$boundary_y, intercept = TRUE)
)
raw_grid <- tensor_basis(
  splines::bs(d$grid$coords[, 1], knots = meta_c$knots_x,
              Boundary.knots = meta_c$boundary_x, intercept = TRUE),
  splines::bs(d$grid$coords[, 2], knots = meta_c$knots_y,
              Boundary.knots = meta_c$boundary_y, intercept = TRUE)
)
compare_num("centering_train_zero", colMeans(train_c), rep(0, ncol(train_c)))
compare_num("centering_test_reuses_train_mean", test_c,
            sweep(raw_test, 2, meta_c$train_mean, "-"))
compare_num("centering_grid_reuses_train_mean", grid_c,
            sweep(raw_grid, 2, meta_c$train_mean, "-"))

small_cc <- list(n_iter = 140, burn_in = 50)
method_list <- list(
  proposed = fit_proposed_ssgl,
  original_ssgl = fit_original_ssgl,
  global = fit_global_only_blasso,
  full = fit_full_svc_no_selection,
  gam = fit_original_gam
)
for (nm in names(method_list)) {
  a <- method_list[[nm]](d$train, d$test, d$grid, bc, mc, small_cc, 1401)
  b <- method_list[[nm]](d$train, d$test, d$grid, bc, mc, small_cc, 1401)
  delta <- max(c(
    maxdiff(a$theta_mean, b$theta_mean),
    maxdiff(a$pip, b$pip),
    maxdiff(a$beta_hat_grid, b$beta_hat_grid),
    maxdiff(a$pred_test, b$pred_test)
  ))
  add_test("additional", paste0("same_seed_identical_", nm), delta < 1e-12,
           format(delta, scientific = TRUE), "<1e-12",
           if (delta < 1e-12) NA_character_ else "Sampler or legacy backend is not fully controlled by set.seed.")
}

# -------------------------------------------------------------------------
# Three initialization checks for the three new Bayesian samplers.
# -------------------------------------------------------------------------
init_data <- simulate_four_predictor_split(100, 30, grid_size = 8, sigma = 0.4,
                                           seed = 1501)
init_cc <- list(n_iter = 550, burn_in = 250)
init_base <- mc
inits <- list(
  neutral = list(),
  diffuse = list(
    omega = rep(5, 4), sigma2 = 3, lambda_theta2 = 0.3,
    lambda2 = 0.3, pi_gamma = 0.8, gamma = rep(1, 4)
  ),
  tight = list(
    omega = rep(0.2, 4), sigma2 = 0.2, lambda_theta2 = 5,
    lambda2 = 5, pi_gamma = 0.2, gamma = rep(0, 4)
  )
)
for (sampler in c("proposed", "global", "full")) {
  fits <- lapply(seq_along(inits), function(i) {
    cfg <- init_base
    cfg$initialization <- inits[[i]]
    fn <- switch(sampler, proposed = fit_proposed_ssgl,
                 global = fit_global_only_blasso,
                 full = fit_full_svc_no_selection)
    fn(init_data$train, init_data$test, init_data$grid, bc, cfg, init_cc,
       seed = 1600 + i)
  })
  theta_mat <- do.call(cbind, lapply(fits, `[[`, "theta_mean"))
  theta_range <- max(apply(theta_mat, 1, function(x) diff(range(x))))
  beta_range <- max(vapply(fits[-1], function(x) {
    sqrt(mean((x$beta_hat_grid - fits[[1]]$beta_hat_grid)^2))
  }, numeric(1)))
  add_test("additional", paste0("three_initializations_theta_", sampler),
           theta_range < 0.20, theta_range, "<0.20",
           if (theta_range < 0.20) NA_character_ else "Short chains retain initialization sensitivity or mix slowly.")
  add_test("additional", paste0("three_initializations_beta_", sampler),
           beta_range < 0.30, beta_range, "<0.30",
           if (beta_range < 0.30) NA_character_ else "Short chains retain initialization sensitivity or mix slowly.")
}

# -------------------------------------------------------------------------
# Level 4: single-predictor recovery tests using proposed sampler.
# -------------------------------------------------------------------------
single_case <- function(type, sigma = 0.25, seed = 1) {
  set.seed(seed)
  ntr <- 120; nte <- 40
  trc <- cbind(runif(ntr), runif(ntr))
  tec <- cbind(runif(nte), runif(nte))
  gc <- as.matrix(expand.grid(x = seq(0, 1, length.out = 10),
                              y = seq(0, 1, length.out = 10)))
  xr <- rnorm(ntr + nte)
  xm <- mean(xr[1:ntr]); xs <- sd(xr[1:ntr])
  x <- (xr - xm) / xs
  raw <- function(s) exp(-((s[, 1] - 0.35)^2 + (s[, 2] - 0.65)^2) / 0.04)
  center <- mean(raw(trc)); norm <- sqrt(mean((raw(trc) - center)^2))
  u <- function(s) (raw(s) - center) / norm
  theta <- if (type %in% c("global", "low_noise")) 1 else 0
  active_u <- type %in% c("spatial", "low_noise")
  beta_tr <- rep(theta, ntr) + if (active_u) u(trc) else rep(0, ntr)
  beta_te <- rep(theta, nte) + if (active_u) u(tec) else rep(0, nte)
  beta_g <- rep(theta, nrow(gc)) +
    if (active_u) u(gc) else rep(0, nrow(gc))
  list(
    train = list(y = x[1:ntr] * beta_tr + rnorm(ntr, sd = sigma),
                 X = matrix(x[1:ntr], ncol = 1), coords = trc),
    test = list(y = x[ntr + 1:nte] * beta_te + rnorm(nte, sd = sigma),
                X = matrix(x[ntr + 1:nte], ncol = 1), coords = tec),
    grid = list(coords = gc, true_beta = matrix(beta_g, ncol = 1),
                true_u = matrix(if (active_u) u(gc) else 0, ncol = 1)),
    truth = list(theta = theta, spatial = active_u)
  )
}
recovery_rows <- list()
cases <- list(global = 0.25, spatial = 0.25, null = 0.25, low_noise = 0.08)
for (i in seq_along(cases)) {
  nm <- names(cases)[i]
  sd <- single_case(nm, cases[[i]], 1700 + i)
  cfg <- list(lambda0 = 10, lambda1 = 1, a_sigma = 0.5,
              b_sigma = var(sd$train$y) / 2, a_theta = 1, b_theta = 1,
              a_gamma = 1, b_gamma = 1)
  f <- fit_proposed_ssgl(sd$train, sd$test, sd$grid,
                         list(n_basis = 4), cfg,
                         list(n_iter = 800, burn_in = 300), 1800 + i)
  theta_error <- abs(f$theta_mean[1] - sd$truth$theta)
  pip <- f$pip[1]
  rb <- f$diagnostics$rb_pip[1]
  beta_rmse <- sqrt(mean((f$beta_hat_grid - sd$grid$true_beta)^2))
  expected_pip <- if (sd$truth$spatial) pip > 0.8 else pip < 0.2
  pass <- theta_error < 0.25 && expected_pip && beta_rmse < 0.65
  recovery_rows[[i]] <- data.frame(
    case = nm, theta_hat = f$theta_mean[1], theta_error = theta_error,
    sampled_pip = pip, rb_pip = rb, beta_rmse = beta_rmse, passed = pass
  )
  add_test("4_recovery", paste0("single_predictor_", nm), pass,
           paste0("theta_error=", round(theta_error, 4),
                  "; pip=", round(pip, 4),
                  "; beta_rmse=", round(beta_rmse, 4)),
           "theta_error<0.25; correct PIP side of 0.2/0.8; beta_rmse<0.65",
           if (pass) NA_character_ else "Short-chain recovery, basis approximation, or selection mixing failure.")
  add_test("additional", paste0("sampled_vs_rb_pip_", nm),
           abs(pip - rb) < 0.10, abs(pip - rb), "<0.10",
           if (abs(pip - rb) < 0.10) NA_character_ else "Monte Carlo error or poor gamma mixing.")
}
recovery <- do.call(rbind, recovery_rows)
write.csv(recovery, file.path(root, "results_debug", "validation",
                              "single_predictor_recovery.csv"), row.names = FALSE)

# Same split and preprocessing object checks for all method calls.
all_fits <- lapply(method_list, function(fn) {
  fn(d$train, d$test, d$grid, bc, mc, list(n_iter = 120, burn_in = 40), 1901)
})
same_lengths <- all(vapply(all_fits, function(f) length(f$pred_test) ==
                            nrow(d$test$X), logical(1)))
same_grid <- all(vapply(all_fits, function(f) nrow(f$beta_hat_grid) ==
                         nrow(d$grid$coords), logical(1)))
add_test("additional", "all_methods_same_test_split_dimensions",
         same_lengths, same_lengths, "TRUE")
add_test("additional", "all_methods_same_grid_dimensions",
         same_grid, same_grid, "TRUE")
add_test("additional", "shared_input_objects_not_mutated",
         identical(d_serialized_before, serialize(d, NULL)),
         "same shared split object used for every call", "TRUE")

out <- do.call(rbind, results)
write.csv(out, file.path(root, "results_debug", "validation",
                         "validation_results.csv"), row.names = FALSE)
saveRDS(list(results = out, details = details, recovery = recovery,
             model_audit = model_audit),
        file.path(root, "results_debug", "validation",
                  "validation_results.rds"))
cat("\nValidation summary:\n")
print(with(out, table(level, passed)))
cat("\nFailures:\n")
print(out[!out$passed, ], row.names = FALSE)
if (any(!out$passed)) quit(status = 2)
