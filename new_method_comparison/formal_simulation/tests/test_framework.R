rm(list = ls())
source("new_method_comparison/formal_simulation/R/framework.R")

config <- merge_formal_config(list(
  n_train = 40L, n_test = 20L, grid_size = 8L,
  n_iter = 30L, burn_in = 10L, basis_dimension = 4L,
  theta_strength = 1, u_target_norm = 0.5, rho_x = 0.3
))
data1 <- generate_formal_dataset(config, 1L)
config2 <- config
config2$theta_strength <- 2
data2 <- generate_formal_dataset(config2, 1L)

stopifnot(
  ncol(data1$train$X) == 10L,
  max(abs(colMeans(data1$train$X))) < 1e-12,
  identical(data1$train$X, data2$train$X),
  identical(data1$train$coords, data2$train$coords),
  identical(data1$random_numbers$noise_z, data2$random_numbers$noise_z)
)

fit_configs <- formal_fit_configs(data1, config)
fit <- fit_proposed_ssgl(
  data1$train, data1$test, data1$grid,
  fit_configs$basis, fit_configs$model, fit_configs$mcmc,
  seed = 77L
)
metrics <- summarize_formal_fit(
  fit, data1, 1L, "proposed_ssgl", 1, 0.5
)
stopifnot(nrow(metrics$predictor) == 10L, nrow(metrics$run) == 1L)
source("new_method_comparison/formal_simulation/R/reporting.R")
source("new_method_comparison/formal_simulation/configs/main_experiment_config.R")
stopifnot(is.na(main_experiment_config$theta_strength))
cat("formal framework reduced-fit test: PASS\n")
