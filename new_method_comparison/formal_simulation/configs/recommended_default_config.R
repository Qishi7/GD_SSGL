# Recommendation from the completed 4 x 3 x 5 proposed-model theta diagnostic.
# This file records the selected default but does not launch any experiment.
recommended_default_config <- list(
  theta_strength = 1.0,
  u_target_norm = 0.5,
  rationale = paste(
    "Moderate global signal; global-only spatial PIPs remain near zero;",
    "u=0.5 spatial PIPs remain near one; theta MSE is stable;",
    "the u=0.25 transition is less compressed than at theta=2."
  )
)

