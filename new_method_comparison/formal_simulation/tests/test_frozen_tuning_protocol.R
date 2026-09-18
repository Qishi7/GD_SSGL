rm(list = ls())
source("new_method_comparison/formal_simulation/configs/frozen_tuning_protocol.R")
x <- read.csv(
  "new_method_comparison/formal_simulation/results_tuning_calibration_first3/all_cv_candidates_frozen.csv"
)

selected_within_grid <- function(method, grid) {
  z <- x[x$method == method, , drop = FALSE]
  keep <- rep(FALSE, nrow(z))
  for (i in seq_len(nrow(grid))) {
    candidate_match <- rep(TRUE, nrow(z))
    for (name in names(grid)) {
      candidate_match <- candidate_match & z[[name]] == grid[[name]][i]
    }
    keep <- keep | candidate_match
  }
  z <- z[keep, , drop = FALSE]
  do.call(rbind, lapply(split(z, z$replicate), function(q) {
    q[which.min(q$cv_mse),
      c("replicate", names(grid), "cv_mse"), drop = FALSE]
  }))
}

proposed <- selected_within_grid(
  "proposed_ssgl",
  frozen_tuning_protocol$proposed_ssgl$candidates
)
original <- selected_within_grid(
  "original_ssgl",
  frozen_tuning_protocol$original_ssgl$candidates
)
full_svc <- selected_within_grid(
  "full_svc_no_selection",
  data.frame(
    kappa2_alpha =
      frozen_tuning_protocol$full_svc_no_selection$kappa2_alpha
  )
)

stopifnot(
  identical(as.numeric(proposed$lambda0), c(20, 20, 10)),
  identical(as.numeric(proposed$lambda1), c(2, 2, 1)),
  identical(as.numeric(original$lambda0), c(15, 20, 20)),
  identical(as.numeric(original$lambda1), c(1.5, 3, 4)),
  identical(as.numeric(full_svc$kappa2_alpha), c(16, 8, 16))
)
cat("frozen tuning protocol: PASS\n")
print(proposed)
print(original)
print(full_svc)

