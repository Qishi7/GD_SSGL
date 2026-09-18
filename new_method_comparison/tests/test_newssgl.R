rm(list = ls())
options(digits = 17)

comparison_root <- normalizePath(
  "new_method_comparison", winslash = "/", mustWork = TRUE
)
source(file.path(
  comparison_root, "R", "methods", "accelerated",
  "newssgl.R"
))
source(file.path(
  comparison_root, "tests", "helpers_newssgl_validation.R"
))
suppressPackageStartupMessages(library(GIGrvg))

out_dir <- file.path(
  comparison_root, "results_debug", "newssgl_validation"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
load_newssgl_fast()

frozen_root <- file.path(
  comparison_root, "formal_simulation", "final_analysis",
  "positive_theta_seed_2026100001"
)
data <- readRDS(file.path(frozen_root, "inputs", "dataset.rds"))
frozen <- readRDS(file.path(frozen_root, "results", "complete_results.rds"))
selected <- read.csv(
  file.path(frozen_root, "results", "selected_settings.csv"),
  stringsAsFactors = FALSE
)

model_config <- list(
  lambda0 = 20, lambda1 = 1,
  a_sigma = 0.5, b_sigma = var(data$train$y) / 2,
  a_theta = 1, b_theta = 1,
  a_gamma = 1, b_gamma = 10
)
basis_config <- list(n_basis = 6)
results <- list()
add <- function(...) {
  results[[length(results) + 1L]] <<- make_test_result(...)
}

# -------------------------------------------------------------------------
# 1. Fixed-state conditional parameters.
# -------------------------------------------------------------------------
prepared <- prepare_newssgl(
  data$train, basis_config, model_config
)
set.seed(81001)
p <- prepared$p
h <- prepared$h
q <- prepared$q
W <- prepared$W
y <- data$train$y
state <- list(
  xi = rnorm(q, sd = 0.08),
  omega = exp(rnorm(p, sd = 0.2)),
  tau = exp(rnorm(p, mean = -1, sd = 0.2)),
  gamma = rep(c(0, 1), length.out = p),
  sigma2 = 0.31,
  lambda_theta2 = 1.7,
  pi_gamma = 0.37
)

cpp_cond <- newssgl_conditionals_cpp(
  y, W, p, h, state$xi, state$omega, state$tau, state$gamma,
  state$sigma2, state$lambda_theta2, state$pi_gamma,
  model_config$a_sigma, model_config$b_sigma,
  model_config$a_theta, model_config$b_theta,
  model_config$a_gamma, model_config$b_gamma,
  model_config$lambda0, model_config$lambda1
)

prior_r <- c(
  1 / state$omega,
  rep(1 / state$tau, each = h)
)
precision_r <- crossprod(W) + diag(prior_r, q)
mean_r <- solve(precision_r, crossprod(W, y))
covariance_r <- state$sigma2 * solve(precision_r)
theta_r <- state$xi[seq_len(p)]
alpha_r <- state$xi[p + seq_len(p * h)]
residual_r <- y - as.vector(W %*% state$xi)
penalty_r <- sum(theta_r^2 / state$omega)
alpha_sq_r <- numeric(p)
for (j in seq_len(p)) {
  idx <- ((j - 1L) * h + 1L):(j * h)
  alpha_sq_r[j] <- sum(alpha_r[idx]^2)
  penalty_r <- penalty_r + alpha_sq_r[j] / state$tau[j]
}
sigma_shape_r <- model_config$a_sigma + 0.5 * (nrow(W) + q)
sigma_rate_r <- model_config$b_sigma +
  0.5 * (sum(residual_r^2) + penalty_r)
norm_scaled_r <- sqrt(alpha_sq_r / state$sigma2)
log_slab_r <- log(state$pi_gamma) +
  h * log(model_config$lambda1) -
  model_config$lambda1 * norm_scaled_r
log_spike_r <- log(1 - state$pi_gamma) +
  h * log(model_config$lambda0) -
  model_config$lambda0 * norm_scaled_r
gamma_prob_r <- plogis(log_slab_r - log_spike_r)

conditional_checks <- list(
  prior_diagonal = c(cpp = cpp_cond$prior_diagonal, r = prior_r),
  precision = c(cpp = cpp_cond$precision, r = precision_r),
  xi_mean = c(cpp = cpp_cond$xi_mean, r = mean_r),
  xi_covariance = c(cpp = cpp_cond$xi_covariance, r = covariance_r),
  sigma_shape = c(cpp = cpp_cond$sigma_shape, r = sigma_shape_r),
  sigma_rate = c(cpp = cpp_cond$sigma_rate, r = sigma_rate_r),
  omega_chi = c(
    cpp = cpp_cond$omega_chi,
    r = theta_r^2 / state$sigma2
  ),
  omega_psi = c(
    cpp = cpp_cond$omega_psi,
    r = rep(state$lambda_theta2, p)
  ),
  lambda_theta_shape = c(
    cpp = cpp_cond$lambda_theta_shape,
    r = model_config$a_theta + p
  ),
  lambda_theta_rate = c(
    cpp = cpp_cond$lambda_theta_rate,
    r = model_config$b_theta + 0.5 * sum(state$omega)
  ),
  pi_shape1 = c(
    cpp = cpp_cond$pi_shape1,
    r = model_config$a_gamma + sum(state$gamma)
  ),
  pi_shape2 = c(
    cpp = cpp_cond$pi_shape2,
    r = model_config$b_gamma + p - sum(state$gamma)
  ),
  gamma_probability = c(
    cpp = cpp_cond$gamma_probability,
    r = gamma_prob_r
  ),
  tau_chi = c(
    cpp = cpp_cond$tau_chi,
    r = alpha_sq_r / state$sigma2
  )
)
for (name in names(conditional_checks)) {
  z <- conditional_checks[[name]]
  half <- length(z) / 2
  difference <- max_abs_difference(z[seq_len(half)], z[half + seq_len(half)])
  add(
    "fixed_conditionals", name, difference, 1e-9,
    difference < 1e-9
  )
}

# -------------------------------------------------------------------------
# 2. Exact reproducibility of newssgl_fast.
# -------------------------------------------------------------------------
repro_mcmc <- list(n_iter = 500L, burn_in = 200L)
repro1 <- fit_newssgl_fast(
  data$train, data$test, data$grid,
  basis_config, model_config, repro_mcmc, 82001L
)
repro2 <- fit_newssgl_fast(
  data$train, data$test, data$grid,
  basis_config, model_config, repro_mcmc, 82001L
)
repro_fields <- list(
  theta = max_abs_difference(repro1$theta_draws, repro2$theta_draws),
  alpha = max_abs_difference(
    repro1$diagnostics$alpha_draws,
    repro2$diagnostics$alpha_draws
  ),
  sigma2 = max_abs_difference(
    repro1$diagnostics$sigma2_draws,
    repro2$diagnostics$sigma2_draws
  ),
  gamma = max_abs_difference(
    repro1$diagnostics$gamma_draws,
    repro2$diagnostics$gamma_draws
  ),
  gamma_probability = max_abs_difference(
    repro1$diagnostics$gamma_prob_draws,
    repro2$diagnostics$gamma_prob_draws
  )
)
for (name in names(repro_fields)) {
  add(
    "reproducibility", name, repro_fields[[name]], 0,
    identical(repro_fields[[name]], 0)
  )
}

# -------------------------------------------------------------------------
# 3. Posterior summary equivalence, including shrinkage draws.
# Two independent chains per implementation are averaged.
# -------------------------------------------------------------------------
validation_mcmc <- list(n_iter = 5000L, burn_in = 2000L)
validation_seeds <- c(83001L, 83002L)
reference_fits <- vector("list", length(validation_seeds))
accelerated_fits <- vector("list", length(validation_seeds))
for (i in seq_along(validation_seeds)) {
  reference_fits[[i]] <- fit_newssgl_reference_test(
    data$train, data$test, data$grid,
    basis_config, model_config, validation_mcmc,
    validation_seeds[i]
  )
  accelerated_fits[[i]] <- fit_newssgl_fast(
    data$train, data$test, data$grid,
    basis_config, model_config, validation_mcmc,
    validation_seeds[i]
  )
  saveRDS(
    list(reference = reference_fits[[i]],
         accelerated = accelerated_fits[[i]]),
    file.path(out_dir, sprintf("validation_chain_%02d.rds", i))
  )
}
reference_summary <- aggregate_posterior_summaries(lapply(
  reference_fits, posterior_summary_newssgl
))
accelerated_summary <- aggregate_posterior_summaries(lapply(
  accelerated_fits, posterior_summary_newssgl
))

posterior_tests <- data.frame(
  test = c(
    "theta_max_abs", "alpha_rmse", "sigma2_relative",
    "omega_relative_rmse", "tau_relative_rmse",
    "lambda_theta2_relative", "pi_gamma_absolute",
    "sampled_pip_max_abs", "rb_pip_max_abs",
    "surface_rmse", "prediction_rmse"
  ),
  value = c(
    max_abs_difference(
      reference_summary$theta, accelerated_summary$theta
    ),
    surface_rmse(
      reference_summary$alpha, accelerated_summary$alpha
    ),
    relative_difference(
      reference_summary$sigma2, accelerated_summary$sigma2
    ),
    surface_rmse(
      reference_summary$omega, accelerated_summary$omega
    ) / max(sqrt(mean(reference_summary$omega^2)), 1e-12),
    surface_rmse(
      reference_summary$tau, accelerated_summary$tau
    ) / max(sqrt(mean(reference_summary$tau^2)), 1e-12),
    relative_difference(
      reference_summary$lambda_theta2,
      accelerated_summary$lambda_theta2
    ),
    abs(
      reference_summary$pi_gamma - accelerated_summary$pi_gamma
    ),
    max_abs_difference(
      reference_summary$sampled_pip,
      accelerated_summary$sampled_pip
    ),
    max_abs_difference(
      reference_summary$rb_pip,
      accelerated_summary$rb_pip
    ),
    surface_rmse(
      reference_summary$beta_grid,
      accelerated_summary$beta_grid
    ),
    surface_rmse(
      reference_summary$pred_test,
      accelerated_summary$pred_test
    )
  ),
  tolerance = c(
    0.05, 0.10, 0.15, 0.30, 0.30, 0.20,
    0.10, 0.10, 0.10, 0.05, 0.08
  ),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(posterior_tests))) {
  add(
    "posterior_equivalence",
    posterior_tests$test[i],
    posterior_tests$value[i],
    posterior_tests$tolerance[i],
    posterior_tests$value[i] <= posterior_tests$tolerance[i]
  )
}

# -------------------------------------------------------------------------
# 4. Agreement with the frozen final primary fit.
# -------------------------------------------------------------------------
frozen_primary <- frozen$fits$proposed_ssgl
frozen_config <- list(n_iter = 20000L, burn_in = 5000L)
accelerated_primary <- fit_newssgl_fast(
  data$train, data$test, data$grid,
  basis_config, model_config, frozen_config, 2026210101L
)
saveRDS(
  accelerated_primary,
  file.path(out_dir, "newssgl_fit.rds")
)
frozen_metrics <- compute_newssgl_metrics(frozen_primary, data)
accelerated_metrics <- compute_newssgl_metrics(
  accelerated_primary, data
)
metric_difference <- abs(frozen_metrics - accelerated_metrics)
metric_tolerance <- c(
  mspe = 0.03, beta_mise = 0.01,
  mise_global_only = 0.01, mise_spatial_only = 0.02,
  mise_global_plus_spatial = 0.02, mise_null = 0.01,
  theta_mse = 0.005, u_mise_x3_x6 = 0.02
)
for (name in names(metric_difference)) {
  add(
    "frozen_primary_metrics", name,
    metric_difference[name], metric_tolerance[name],
    metric_difference[name] <= metric_tolerance[name]
  )
}
frozen_comparison <- c(
  theta_max_abs = max_abs_difference(
    frozen_primary$theta_mean, accelerated_primary$theta_mean
  ),
  alpha_rmse = surface_rmse(
    rowMeans(frozen_primary$diagnostics$alpha_draws),
    accelerated_primary$alpha_mean
  ),
  sampled_pip_max_abs = max_abs_difference(
    frozen_primary$pip, accelerated_primary$pip
  ),
  rb_pip_max_abs = max_abs_difference(
    frozen_primary$diagnostics$rb_pip,
    accelerated_primary$diagnostics$rb_pip
  ),
  surface_rmse = surface_rmse(
    frozen_primary$beta_hat_grid,
    accelerated_primary$beta_hat_grid
  )
)
frozen_tolerance <- c(
  theta_max_abs = 0.05, alpha_rmse = 0.10,
  sampled_pip_max_abs = 0.05, rb_pip_max_abs = 0.05,
  surface_rmse = 0.05
)
for (name in names(frozen_comparison)) {
  add(
    "frozen_primary_summary", name,
    frozen_comparison[name], frozen_tolerance[name],
    frozen_comparison[name] <= frozen_tolerance[name]
  )
}

# -------------------------------------------------------------------------
# 5. Runtime comparison against the untouched validated wrapper.
# Compilation is already complete and excluded.
# -------------------------------------------------------------------------
runtime_mcmc <- list(n_iter = 1200L, burn_in = 400L)
runtime_seed <- 84001L
reference_runtime <- system.time({
  runtime_reference_fit <- fit_newssgl_reference(
    data$train, data$test, data$grid,
    basis_config, model_config, runtime_mcmc, runtime_seed
  )
})[["elapsed"]]
accelerated_runtime <- system.time({
  runtime_accelerated_fit <- fit_newssgl_fast(
    data$train, data$test, data$grid,
    basis_config, model_config, runtime_mcmc, runtime_seed
  )
})[["elapsed"]]
speedup <- reference_runtime / accelerated_runtime
add(
  "runtime", "speedup", speedup, 1,
  speedup > 1,
  sprintf(
    "newssgl %.3fs; newssgl_fast %.3fs",
    reference_runtime, accelerated_runtime
  )
)

runtime_table <- data.frame(
  implementation = c("newssgl", "newssgl_fast"),
  n_iter = runtime_mcmc$n_iter,
  burn_in = runtime_mcmc$burn_in,
  elapsed_seconds = c(reference_runtime, accelerated_runtime),
  speedup_vs_reference = c(1, speedup)
)
write.csv(
  runtime_table, file.path(out_dir, "runtime_comparison.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    metric = names(frozen_metrics),
    frozen_reference = as.numeric(frozen_metrics),
    newssgl_fast = as.numeric(accelerated_metrics),
    absolute_difference = as.numeric(metric_difference),
    tolerance = as.numeric(metric_tolerance)
  ),
  file.path(out_dir, "metric_comparison.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    predictor = paste0("X", 1:10),
    frozen_sampled_pip = frozen_primary$pip,
    newssgl_fast_sampled_pip = accelerated_primary$pip,
    frozen_rb_pip = frozen_primary$diagnostics$rb_pip,
    newssgl_fast_rb_pip =
      accelerated_primary$diagnostics$rb_pip
  ),
  file.path(out_dir, "pip_comparison.csv"),
  row.names = FALSE
)
write.csv(
  posterior_tests,
  file.path(out_dir, "posterior_summary_comparison.csv"),
  row.names = FALSE
)

test_results <- do.call(rbind, results)
write.csv(
  test_results, file.path(out_dir, "validation_tests.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    test_results = test_results,
    posterior_tests = posterior_tests,
    frozen_metrics = frozen_metrics,
    accelerated_metrics = accelerated_metrics,
    frozen_comparison = frozen_comparison,
    runtime = runtime_table,
    selected_settings = selected
  ),
  file.path(out_dir, "validation_results.rds")
)

cat("\nValidation summary\n")
print(aggregate(passed ~ group, test_results, function(z) {
  sprintf("%d/%d", sum(z), length(z))
}))
cat("\nRuntime\n")
print(runtime_table)
if (!all(test_results$passed)) {
  print(test_results[!test_results$passed, ])
  stop("newssgl validation failed.")
}
cat("\nALL_NEWSSGL_TESTS_PASSED\n")
