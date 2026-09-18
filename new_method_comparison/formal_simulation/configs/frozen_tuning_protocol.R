# Frozen after the training-only three-replicate calibration.
# H = 36 is common to all tensor-basis methods.

frozen_tuning_protocol <- list(
  folds = list(
    k = 5L,
    seed = 123L,
    construction = "sample(rep(1:k, length.out=n_train))",
    loss = "mean validation response squared error",
    tie_rule = "first minimum in listed candidate order"
  ),
  cv_mcmc = list(n_iter = 500L, burn_in = 100L),
  proposed_ssgl = list(
    H = 36L,
    candidates = data.frame(
      lambda0 = c(5, 10, 10, 20, 20, 20, 30, 30),
      lambda1 = c(0.5, 1, 2, 1, 2, 3, 2, 4)
    )
  ),
  original_ssgl = list(
    H = 36L,
    zeta0 = 0.1,
    zeta1 = 1,
    candidates = data.frame(
      lambda0 = c(5, 10, 15, 20, 15, 20, 20, 30),
      lambda1 = c(0.5, 1, 1.5, 2, 3, 3, 4, 5)
    )
  ),
  full_svc_no_selection = list(
    H = 36L,
    kappa2_alpha = c(0.5, 1, 2, 4, 8, 16, 32)
  ),
  global_only_blasso = list(
    external_cv = FALSE,
    lambda_theta2 = "sampled internally",
    a_theta = 1,
    b_theta = 1
  ),
  standard_thin_plate_gam = list(
    method = "GCV.Cp",
    k = c(20L, 30L, 40L, 50L)
  )
)

