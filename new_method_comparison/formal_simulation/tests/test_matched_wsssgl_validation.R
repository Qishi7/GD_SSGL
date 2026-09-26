#!/usr/bin/env Rscript

options(digits = 16, warn = 1)

project_root <- "/Users/bayeslab/Desktop/ssgl"
comparison_root <- file.path(project_root, "new_method_comparison")

source(file.path(comparison_root, "R", "methods", "accelerated",
                 "matched_wsssgl_intercept.R"))
load_matched_wsssgl_intercept_fast(rebuild = TRUE)

out_root <- file.path(comparison_root, "formal_simulation", "tests",
                      "matched_wsssgl_validation")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

check <- function(name, value, tolerance = NA_real_, passed = isTRUE(value)) {
  value_text <- if (is.numeric(value)) {
    if (length(value) == 1) as.character(signif(value, 8)) else
      paste(signif(value, 8), collapse = ",")
  } else {
    paste(as.character(value), collapse = ",")
  }
  data.frame(
    check = name,
    value = value_text,
    tolerance = if (is.na(tolerance)) "" else as.character(tolerance),
    passed = passed,
    stringsAsFactors = FALSE
  )
}

set.seed(20260925)
dat <- simulate_p10_pilot_split(n_train = 120, n_test = 50, grid_size = 12,
                                sigma = 0.35, seed = 20260925)
dat$truth$beta0 <- 0.5
dat$train$y <- dat$train$y + 0.5
dat$test$y <- dat$test$y + 0.5

prep <- prepare_matched_wsssgl_intercept(
  dat$train,
  list(n_basis = 6, full_rank_method = "svd", rank_tolerance = 1e-10),
  list(lambda0 = 20, lambda1 = 2)
)
p <- ncol(dat$train$X)
r <- prep$r
H <- prep$H

theta <- rnorm(p)
alpha <- rnorm(p * r)
eta <- numeric(p * H)
for (j in seq_len(p)) {
  eta_idx <- ((j - 1L) * H + 1L):(j * H)
  alpha_idx <- ((j - 1L) * r + 1L):(j * r)
  eta[eta_idx] <- c(theta[j], alpha[alpha_idx])
}
beta_gd_train <- surface_from_components(theta, alpha, prep$phi, p)
beta_ws_train <- surface_from_eta(eta, prep$q_basis, p)
q_test <- build_whole_surface_basis(dat$test$coords, prep$meta)
phi_test <- apply_centered_basis_metadata_fullrank(dat$test$coords, prep$meta)
beta_gd_test <- surface_from_components(theta, alpha, phi_test, p)
beta_ws_test <- surface_from_eta(eta, q_test, p)

lambda0 <- 20
lambda1 <- 2
pi_gamma <- 0.3
eta_norm_scaled <- 1.7
logit_q <- log(pi_gamma / (1 - pi_gamma)) +
  H * log(lambda1 / lambda0) - (lambda1 - lambda0) * eta_norm_scaled
prob_q <- plogis(logit_q)
prob_equal_lambda <- plogis(log(pi_gamma / (1 - pi_gamma)) +
                            H * log(lambda1 / lambda1))
prob_larger_signal <- plogis(log(pi_gamma / (1 - pi_gamma)) +
                             H * log(lambda1 / lambda0) -
                             (lambda1 - lambda0) * (eta_norm_scaled + 0.5))

fit <- fit_matched_wsssgl_intercept_fast(
  dat$train, dat$test, dat$grid,
  list(n_basis = 6, full_rank_centered = TRUE, full_rank_method = "svd",
       rank_tolerance = 1e-10),
  list(lambda0 = 20, lambda1 = 2, a_sigma = 0.5,
       b_sigma = var(dat$train$y) / 2, a_gamma = 1, b_gamma = p),
  list(n_iter = 150, burn_in = 50),
  seed = 20260925
)

results <- rbind(
  check("centered_svd_rank_expected_35", r, 35, r == 35),
  check("whole_surface_block_dimension_expected_36", H, 36, H == 36),
  check("training_center_error", prep$center_error, 1e-8,
        prep$center_error < 1e-8),
  check("gd_ws_basis_equivalence_train",
        max(abs(beta_gd_train - beta_ws_train)), 1e-10,
        max(abs(beta_gd_train - beta_ws_train)) < 1e-10),
  check("gd_ws_basis_equivalence_test",
        max(abs(beta_gd_test - beta_ws_test)), 1e-10,
        max(abs(beta_gd_test - beta_ws_test)) < 1e-10),
  check("lambda_equal_probability_reduces_to_pi", abs(prob_equal_lambda - pi_gamma),
        1e-12, abs(prob_equal_lambda - pi_gamma) < 1e-12),
  check("q_increases_with_eta_norm_when_lambda0_gt_lambda1",
        prob_larger_signal - prob_q, 0, prob_larger_signal > prob_q),
  check("sampler_returns_finite_prediction",
        all(is.finite(fit$pred_test)), NA_real_, all(is.finite(fit$pred_test))),
  check("sampler_eta_draw_dimension",
        paste(dim(fit$eta_draws), collapse = "x"), NA_real_,
        identical(dim(fit$eta_draws), c(p * H, 100L))),
  check("sampler_rb_pip_dimension",
        length(fit$diagnostics$rb_pip), p,
        length(fit$diagnostics$rb_pip) == p)
)

write.csv(results, file.path(out_root, "matched_wsssgl_validation_checks.csv"),
          row.names = FALSE)
saveRDS(fit, file.path(out_root, "matched_wsssgl_short_fit.rds"))

cat("# Matched WS-SSGL validation checks\n\n")
print(results, row.names = FALSE)
cat("\nAll checks passed:", all(results$passed), "\n")
