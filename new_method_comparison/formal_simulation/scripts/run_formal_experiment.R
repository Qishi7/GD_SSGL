rm(list = ls())
formal_root <- if (file.exists("new_method_comparison/formal_simulation/R/framework.R")) {
  normalizePath("new_method_comparison/formal_simulation", winslash = "/")
} else {
  normalizePath("formal_simulation", winslash = "/")
}
source(file.path(formal_root, "R", "framework.R"))
source(file.path(formal_root, "R", "reporting.R"))
source(file.path(formal_root, "configs", "main_experiment_config.R"))

if (is.na(main_experiment_config$theta_strength)) {
  stop("Set theta_strength after reviewing the theta diagnostic.")
}
config <- merge_formal_config(main_experiment_config)
output_dir <- file.path(formal_root, "results", "main_experiment")
for (replicate_id in seq_len(config$n_replicates)) {
  message("Formal replicate ", replicate_id, "/", config$n_replicates)
  run_formal_replicate(config, replicate_id, output_dir = output_dir,
                       resume = TRUE)
}
build_formal_reports(output_dir)
