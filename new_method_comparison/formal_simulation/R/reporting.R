if (!exists(".formal_root")) {
  source(file.path("new_method_comparison", "formal_simulation",
                   "R", "framework.R"))
}

formal_method_summary <- function(aggregated, output_dir) {
  keys <- split(aggregated$run, aggregated$run$method)
  summary <- do.call(rbind, lapply(keys, function(x) {
    data.frame(
      method = x$method[1],
      n_completed = nrow(x),
      mspe_mean = mean(x$mspe), mspe_sd = sd(x$mspe),
      beta_mse_active_mean = mean(x$beta_mse_active),
      beta_mse_active_sd = sd(x$beta_mse_active),
      beta_mse_null_mean = mean(x$beta_mse_null),
      beta_mse_null_sd = sd(x$beta_mse_null),
      theta_mse_mean = mean(x$theta_mse, na.rm = TRUE),
      theta_mse_sd = sd(x$theta_mse, na.rm = TRUE),
      runtime_seconds_mean = mean(x$runtime_seconds),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(summary, file.path(output_dir, "method_comparison_summary.csv"),
            row.names = FALSE)
  summary
}

formal_pip_summary <- function(aggregated, output_dir) {
  x <- aggregated$predictor[
    is.finite(aggregated$predictor$rb_pip), ]
  if (!nrow(x)) return(NULL)
  keys <- interaction(x$method, x$predictor, drop = TRUE)
  summary <- do.call(rbind, lapply(split(x, keys), function(z) {
    data.frame(
      method = z$method[1], predictor = z$predictor[1],
      predictor_index = z$predictor_index[1], group = z$group[1],
      rb_pip_mean = mean(z$rb_pip), rb_pip_sd = sd(z$rb_pip),
      sampled_pip_mean = mean(z$sampled_pip),
      sampled_pip_sd = sd(z$sampled_pip),
      stringsAsFactors = FALSE
    )
  }))
  summary <- summary[order(summary$method, summary$predictor_index), ]
  write.csv(summary, file.path(output_dir, "pip_summary.csv"),
            row.names = FALSE)
  png(file.path(output_dir, "pip_summary.png"),
      width = 1500, height = 900, res = 150)
  methods <- unique(summary$method)
  par(mfrow = c(length(methods), 1), mar = c(4, 4, 3, 1))
  for (method in methods) {
    z <- summary[summary$method == method, ]
    plot(z$predictor_index, z$rb_pip_mean, ylim = c(0, 1),
         type = "b", pch = 19, xlab = "Predictor",
         ylab = "Mean RB PIP", main = method, xaxt = "n")
    axis(1, at = z$predictor_index, labels = z$predictor)
    arrows(z$predictor_index, pmax(0, z$rb_pip_mean - z$rb_pip_sd),
           z$predictor_index, pmin(1, z$rb_pip_mean + z$rb_pip_sd),
           angle = 90, code = 3, length = 0.04)
  }
  dev.off()
  summary
}

formal_reconstruction_figures <- function(output_dir, replicate_id = 1L) {
  replicate_dir <- file.path(
    output_dir, sprintf("replicate_%03d", replicate_id)
  )
  data <- readRDS(file.path(replicate_dir, "shared_dataset.rds"))
  figure_dir <- file.path(output_dir, "reconstruction_figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  plot_surface_panel(
    data$grid$true_beta, data$grid$grid_size,
    file.path(figure_dir, sprintf("replicate_%03d_truth.png", replicate_id)),
    sprintf("Replicate %03d truth", replicate_id)
  )
  method_files <- setdiff(
    list.files(replicate_dir, pattern = "\\.rds$", full.names = TRUE),
    file.path(replicate_dir, "shared_dataset.rds")
  )
  for (path in method_files) {
    object <- readRDS(path)
    method <- object$fit$method_name
    plot_surface_panel(
      object$fit$beta_hat_grid, data$grid$grid_size,
      file.path(
        figure_dir,
        sprintf("replicate_%03d_%s_reconstruction.png",
                replicate_id, method)
      ),
      sprintf("Replicate %03d: %s reconstruction", replicate_id, method)
    )
    plot_surface_panel(
      object$fit$beta_hat_grid - data$grid$true_beta,
      data$grid$grid_size,
      file.path(
        figure_dir,
        sprintf("replicate_%03d_%s_error.png", replicate_id, method)
      ),
      sprintf("Replicate %03d: %s error", replicate_id, method)
    )
  }
  invisible(figure_dir)
}

build_formal_reports <- function(output_dir, reconstruction_replicate = 1L) {
  aggregated <- aggregate_saved_results(output_dir)
  method_table <- formal_method_summary(aggregated, output_dir)
  pip_table <- formal_pip_summary(aggregated, output_dir)
  formal_reconstruction_figures(output_dir, reconstruction_replicate)
  invisible(list(method_table = method_table, pip_table = pip_table))
}

