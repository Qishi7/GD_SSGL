# Prepared configuration only. Do not execute until the theta diagnostic has
# been reviewed and the manuscript-scale experiment is explicitly authorized.
main_experiment_config <- list(
  n_train = 700L,
  n_test = 250L,
  p = 10L,
  theta_strength = NA_real_, # fill from theta diagnostic recommendation
  u_target_norm = 0.5,
  rho_x = 0,
  noise_sd = 0.35,
  basis_dimension = 5L,
  n_iter = 1200L,
  burn_in = 400L,
  n_replicates = 100L,
  seed_start = 2026070000L,
  grid_size = 50L,
  lambda0 = 10,
  lambda1 = 1
)

