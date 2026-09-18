#!/usr/bin/env Rscript

options(warn = 1)

suppressPackageStartupMessages({
  library(sp)
})

main_root <- file.path(normalizePath(Sys.getenv("GDSSGL_ROOT", getwd()), winslash = "/", mustWork = TRUE), "new_method_comparison", "formal_simulation", "main_simulation")
out_root <- file.path(main_root, "beta_reconstruction_single_beta_figure_pack")
script_out_dir <- file.path(out_root, "scripts")
dir.create(script_out_dir, recursive = TRUE, showWarnings = FALSE)

vars <- paste0("X", 1:10)
method_order <- c("True", "GD-SSGL", "WS-SSGL", "Gaussian SVC", "MGWR", "Bayesian Lasso")
pal <- grDevices::colorRampPalette(c("#27316f", "#2b8cbe", "#1fa187", "#fde725"))(256)

scenarios <- list(
  list(
    label = "Main four-function",
    slug = "main_four_function",
    root = file.path(main_root, "four_function_beta0_0p5_100rep_perrep_tuned_intercept_svd"),
    data_file = file.path("fits", "rep_001", "shared_dataset_beta0_0p5.rds"),
    fit_dir = file.path("fits", "rep_001"),
    mgwr_dir = file.path("mgwr_extension_100rep", "fits", "rep_001"),
    methods = method_order
  ),
  list(
    label = "No spatial deviation",
    slug = "no_spatial_deviation",
    root = file.path(main_root, "no_spatial_deviation_beta0_0p5_100rep_perrep_tuned_intercept_svd"),
    data_file = file.path("fits", "rep_001", "shared_dataset_beta0_0p5.rds"),
    fit_dir = file.path("fits", "rep_001"),
    mgwr_dir = file.path("mgwr_extension_100rep", "fits", "rep_001"),
    methods = method_order
  ),
  list(
    label = "Weak spatial signal gamma=0.20",
    slug = "weak_spatial_gamma0p20",
    root = file.path(main_root, "weak_spatial_gamma_0p15_0p25_beta0_0p5_intercept_svd_GDSSGL", "perrep_tuned_R100"),
    data_file = file.path("data", "gamma_0p200", "weak_gamma_beta0_0p5_rep_011.rds"),
    fit_dir = file.path("fits", "gamma_0p200", "rep_011"),
    mgwr_dir = NA_character_,
    methods = c("True", "GD-SSGL")
  ),
  list(
    label = "Correlated predictors rho=0.7",
    slug = "rho07_correlated_predictors",
    root = file.path(main_root, "correlated_predictor_rho07_beta0_0p5_50rep_perrep_tuned_intercept_svd"),
    data_file = file.path("fits", "rep_001", "shared_dataset_beta0_0p5.rds"),
    fit_dir = file.path("fits", "rep_001"),
    mgwr_dir = file.path("mgwr_extension_R50", "fits", "rep_001"),
    methods = method_order
  )
)

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Could not write ", path)
  invisible(path)
}

load_beta_field <- function(path) {
  if (!file.exists(path)) stop("Missing beta fit file: ", path)
  obj <- readRDS(path)
  if (!is.null(obj$beta_grid)) return(obj$beta_grid)
  if (!is.null(obj$beta_hat_grid)) return(obj$beta_hat_grid)
  stop("No beta_grid or beta_hat_grid found in ", path)
}

idw_coefficients <- function(train_coords, coef_mat, new_coords, power = 2) {
  train_coords <- as.matrix(train_coords)
  new_coords <- as.matrix(new_coords)
  coef_mat <- as.matrix(coef_mat)
  out <- matrix(NA_real_, nrow(new_coords), ncol(coef_mat))
  colnames(out) <- colnames(coef_mat)
  for (i in seq_len(nrow(new_coords))) {
    d <- sqrt(rowSums((train_coords - matrix(new_coords[i, ],
                                             nrow(train_coords),
                                             ncol(train_coords),
                                             byrow = TRUE))^2))
    exact <- which(d < 1e-12)
    if (length(exact)) {
      out[i, ] <- coef_mat[exact[1], ]
    } else {
      w <- 1 / d^power
      w <- w / sum(w)
      out[i, ] <- colSums(coef_mat * w)
    }
  }
  out
}

mgwr_beta_grid <- function(path, grid_coords, p) {
  if (!file.exists(path)) stop("Missing MGWR fit file: ", path)
  fit <- readRDS(path)
  d <- fit$SDF@data
  intercept_col <- intersect(c("Intercept", "(Intercept)"), names(d))[1]
  x_cols <- paste0("X", seq_len(p))
  missing <- setdiff(c(intercept_col, x_cols), names(d))
  if (length(missing)) stop("Missing MGWR coefficient columns: ", paste(missing, collapse = ", "))
  coef_train <- as.matrix(d[, c(intercept_col, x_cols), drop = FALSE])
  train_coords <- sp::coordinates(fit$SDF)
  coef_grid <- idw_coefficients(train_coords, coef_train, grid_coords)
  beta <- coef_grid[, -1, drop = FALSE]
  colnames(beta) <- x_cols
  beta
}

load_scenario_betas <- function(sc, dat) {
  fit_dir <- file.path(sc$root, sc$fit_dir)
  out <- list(True = dat$grid$true_beta)

  paths <- list(
    "GD-SSGL" = c(file.path(fit_dir, "gdssgl_final_fit.rds"),
                  file.path(fit_dir, "gdssgl_tuned_final_fit.rds")),
    "WS-SSGL" = file.path(fit_dir, "wsssgl_final_fit.rds"),
    "Gaussian SVC" = file.path(fit_dir, "gaussian_svc_final_fit.rds"),
    "Bayesian Lasso" = file.path(fit_dir, "bayesian_lasso_final_fit.rds")
  )

  for (nm in names(paths)) {
    if (nm %in% sc$methods) {
      candidates <- paths[[nm]]
      path <- candidates[file.exists(candidates)][1]
      if (!is.na(path)) out[[nm]] <- load_beta_field(path)
    }
  }

  if ("MGWR" %in% sc$methods && !is.na(sc$mgwr_dir)) {
    path <- file.path(sc$root, sc$mgwr_dir, "mgwr_fit.rds")
    if (file.exists(path)) out[["MGWR"]] <- mgwr_beta_grid(path, dat$grid$coords, length(vars))
  }

  out[sc$methods[sc$methods %in% names(out)]]
}

draw_vertical_colorbar <- function(lim, pal) {
  format_tick <- function(x) {
    ax <- max(abs(x), na.rm = TRUE)
    if (ax < 0.01) return(sprintf("%.3f", x))
    if (ax < 0.1) return(sprintf("%.2f", x))
    if (ax < 2) return(sprintf("%.1f", x))
    sprintf("%.0f", x)
  }
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(-lim, lim), xaxs = "i", yaxs = "i")
  bar_low <- -0.83 * lim
  bar_high <- 0.83 * lim
  rasterImage(as.raster(matrix(rev(pal), ncol = 1)), 0.12, bar_low, 0.42, bar_high)
  ticks <- pretty(c(-lim, lim), n = 5)
  ticks <- ticks[ticks >= bar_low & ticks <= bar_high]
  axis(4, at = ticks, labels = format_tick(ticks), cex.axis = 2.45,
       lwd = 0.75, lwd.ticks = 0.75, las = 1, line = -2.25)
}

plot_one_beta <- function(beta_list, dat, sc, j, out_dir) {
  coords <- as.data.frame(dat$grid$coords)
  names(coords) <- c("s1", "s2")
  grid_size <- dat$grid$grid_size
  methods <- names(beta_list)
  vals <- unlist(lapply(beta_list, function(mat) mat[, j]))
  lim <- max(abs(vals), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1

  base <- sprintf("%s_beta_%02d", sc$slug, j)
  path <- file.path(out_dir, paste0(base, ".png"))
  png(path, width = max(2700, 760 * length(methods) + 470), height = 1010, res = 180)
  on.exit(dev.off(), add = TRUE)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  layout(matrix(seq_len(length(methods) + 2), nrow = 1),
         widths = c(0.18, rep(1, length(methods)), 0.34))
  par(oma = c(0.08, 0.08, 0.35, 0.08), mar = c(0.08, 0.02, 2.65, 0.02))

  plot.new()
  text(0.62, 0.50, sprintf("Beta %d", j), srt = 90, cex = 4.85, font = 2)

  for (idx_m in seq_along(methods)) {
    m <- methods[idx_m]
    mat <- beta_list[[m]]
    z <- matrix(mat[, j], nrow = grid_size, ncol = grid_size)
    image(sort(unique(coords$s1)), sort(unique(coords$s2)), z,
          col = pal, zlim = c(-lim, lim), axes = FALSE,
          xlab = "", ylab = "", useRaster = TRUE, asp = 1)
    title(main = m, line = -0.92, cex.main = 3.95, font.main = 2)
  }

  par(mar = c(0.08, 0.02, 2.65, 2.50))
  draw_vertical_colorbar(lim, pal)
  invisible(data.frame(
    scenario = sc$slug,
    scenario_label = sc$label,
    beta = j,
    variable = vars[j],
    file = path,
    methods = paste(methods, collapse = "; "),
    min_beta_across_methods = min(vals, na.rm = TRUE),
    max_beta_across_methods = max(vals, na.rm = TRUE),
    symmetric_abs_limit = lim,
    stringsAsFactors = FALSE
  ))
}

plot_scenario <- function(sc) {
  dat_path <- file.path(sc$root, sc$data_file)
  if (!file.exists(dat_path)) stop("Missing scenario data file: ", dat_path)
  dat <- readRDS(dat_path)
  beta_list <- load_scenario_betas(sc, dat)
  out_dir <- file.path(out_root, sc$slug)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  inv <- list()
  k <- 1L
  for (j in seq_along(vars)) {
    inv[[k]] <- plot_one_beta(beta_list, dat, sc, j, out_dir)
    k <- k + 1L
  }
  do.call(rbind, inv)
}

inventory <- do.call(rbind, lapply(scenarios, plot_scenario))
inventory_path <- file.path(out_root, "single_beta_figure_inventory.csv")
atomic_write_csv(inventory, inventory_path)

report_path <- file.path(out_root, "single_beta_figure_pack_report.md")
report_lines <- c(
  "# Single-beta reconstruction figure pack",
  "",
  paste0("- Created: ", format(Sys.time())),
  paste0("- Output root: `", out_root, "`."),
  "- Each beta is plotted in a separate file.",
  "- Each beta uses its own symmetric color range across the available methods for that scenario.",
  "- MGWR is included where the scenario has completed MGWR fit objects.",
  "",
  "## Scenario folders",
  "",
  paste0("- `", unique(inventory$scenario), "`", collapse = "\n"),
  "",
  "## Inventory",
  "",
  paste(capture.output(print(inventory[, c("scenario", "beta", "file", "symmetric_abs_limit")],
                             row.names = FALSE)), collapse = "\n")
)
writeLines(report_lines, report_path)

cat("Generated single-beta figure pack\n")
cat("Output root:", out_root, "\n")
cat("Inventory:", inventory_path, "\n")
cat("Report:", report_path, "\n")
print(inventory[, c("scenario", "beta", "file", "symmetric_abs_limit")], row.names = FALSE)
