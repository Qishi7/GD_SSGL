suppressPackageStartupMessages({
  library(ggplot2)
  library(sp)
})

project_root <- normalizePath(Sys.getenv("GDSSGL_ROOT", getwd()), winslash = "/", mustWork = TRUE)
root <- file.path(project_root, "new_method_comparison", "formal_simulation", "main_simulation", "four_function_beta0_0p5_100rep_perrep_tuned_intercept_svd")
rep_id <- "rep_001"
rep_dir <- file.path(root, "fits", rep_id)
out_dir <- file.path(root, "figures_per_beta_range")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(file.path(rep_dir, "shared_dataset_beta0_0p5.rds"))
coords <- as.data.frame(dat$grid$coords)
names(coords) <- c("s1", "s2")
vars <- paste0("X", 1:10)

load_beta <- function(file, field = "beta_grid") {
  obj <- readRDS(file.path(rep_dir, file))
  if (!is.null(obj[[field]])) return(obj[[field]])
  if (!is.null(obj$beta_hat_grid)) return(obj$beta_hat_grid)
  stop("Cannot find beta grid in ", file)
}

beta_list <- list(
  True = dat$grid$true_beta,
  `GD-SSGL` = load_beta("gdssgl_final_fit.rds"),
  `WS-SSGL` = load_beta("wsssgl_final_fit.rds"),
  `Gaussian SVC` = load_beta("gaussian_svc_final_fit.rds"),
  GAM = load_beta("gam_final_fit.rds"),
  MGWR = NULL,
  BLasso = load_beta("bayesian_lasso_final_fit.rds")
)

mgwr_fit_path <- file.path(root, "mgwr_extension_100rep", "fits", rep_id, "mgwr_fit.rds")
if (file.exists(mgwr_fit_path)) {
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
  extract_coef_matrix <- function(fit, p) {
    d <- fit$SDF@data
    intercept_col <- intersect(c("Intercept", "(Intercept)"), names(d))[1]
    x_cols <- paste0("X", seq_len(p))
    missing <- setdiff(c(intercept_col, x_cols), names(d))
    if (length(missing)) stop("Missing MGWR coefficient columns: ",
                              paste(missing, collapse = ", "))
    as.matrix(d[, c(intercept_col, x_cols), drop = FALSE])
  }
  mgwr_fit <- readRDS(mgwr_fit_path)
  mgwr_coef_train <- extract_coef_matrix(mgwr_fit, length(vars))
  mgwr_train_coords <- sp::coordinates(mgwr_fit$SDF)
  mgwr_coef_grid <- idw_coefficients(mgwr_train_coords, mgwr_coef_train, dat$grid$coords)
  beta_list$MGWR <- mgwr_coef_grid[, -1, drop = FALSE]
  colnames(beta_list$MGWR) <- vars
} else {
  warning("MGWR fit not found; MGWR row will be omitted: ", mgwr_fit_path)
  beta_list$MGWR <- NULL
}

beta_list <- beta_list[!vapply(beta_list, is.null, logical(1))]

method_levels <- names(beta_list)

surface_df <- do.call(rbind, lapply(method_levels, function(m) {
  mat <- beta_list[[m]]
  do.call(rbind, lapply(seq_along(vars), function(j) {
    data.frame(
      method = m,
      variable = vars[j],
      s1 = coords$s1,
      s2 = coords$s2,
      beta = mat[, j],
      stringsAsFactors = FALSE
    )
  }))
}))

surface_df$method <- factor(surface_df$method, levels = rev(method_levels))
surface_df$variable <- factor(surface_df$variable, levels = vars)

png_file <- file.path(out_dir, "main_four_function_rep001_beta_reconstruction_per_beta_range.png")
pdf_file <- file.path(out_dir, "main_four_function_rep001_beta_reconstruction_per_beta_range.pdf")

range_table <- do.call(rbind, lapply(vars, function(v) {
  sub <- subset(surface_df, variable == v)
  data.frame(
    variable = v,
    min_beta_across_methods = min(sub$beta, na.rm = TRUE),
    max_beta_across_methods = max(sub$beta, na.rm = TRUE),
    symmetric_abs_limit = max(abs(sub$beta), na.rm = TRUE)
  )
}))
write.csv(range_table, file.path(out_dir, "main_four_function_rep001_per_beta_color_ranges.csv"), row.names = FALSE)

pal_fun <- grDevices::colorRampPalette(c("#27316f", "#2b8cbe", "#1fa187", "#fde725"))
pal <- pal_fun(256)
grid_size <- dat$grid$grid_size

draw_horizontal_colorbar <- function(lim, pal) {
  plot.new()
  plot.window(xlim = c(-lim, lim), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  rasterImage(as.raster(matrix(pal, nrow = 1)), -lim, 0.15, lim, 0.85)
  ticks <- pretty(c(-lim, lim), n = 3)
  axis(1, at = ticks, labels = sprintf("%.2g", ticks), cex.axis = 0.45,
       lwd = 0.35, lwd.ticks = 0.35, padj = -0.55)
  box(lwd = 0.35, col = "grey35")
}

draw_panel <- function(device_file = NULL, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(device_file, width = 7200, height = 3900, res = 300)
  } else {
    pdf(device_file, width = 24, height = 13)
  }
  on.exit(dev.off(), add = TRUE)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  layout_mat <- matrix(seq_len((length(method_levels) + 1) * length(vars)),
                       nrow = length(method_levels) + 1,
                       byrow = TRUE)
  layout(layout_mat, heights = c(rep(1, length(method_levels)), 0.18))
  par(oma = c(1.2, 5.8, 0.6, 0.5), mar = c(0.15, 0.15, 0.85, 0.15))

  for (m in method_levels) {
    mat <- beta_list[[m]]
    for (j in seq_along(vars)) {
      v <- vars[j]
      lim <- range_table$symmetric_abs_limit[range_table$variable == v]
      z <- matrix(mat[, j], nrow = grid_size, ncol = grid_size)
      image(
        x = sort(unique(coords$s1)),
        y = sort(unique(coords$s2)),
        z = z,
        col = pal,
        zlim = c(-lim, lim),
        axes = FALSE,
        xlab = "",
        ylab = "",
        useRaster = TRUE
      )
      box(lwd = 0.7, col = "grey25")
      if (m == method_levels[1]) {
        title(main = sprintf("Beta %d", j), line = 0.05, cex.main = 0.95, font.main = 2)
      }
      if (j == 1) {
        mtext(m, side = 2, line = 0.55, las = 1, cex = 0.85, font = 2)
      }
    }
  }
  par(mar = c(1.05, 0.15, 0.15, 0.15))
  for (j in seq_along(vars)) {
    lim <- range_table$symmetric_abs_limit[range_table$variable == vars[j]]
    draw_horizontal_colorbar(lim, pal)
  }
}

draw_panel(png_file, "png")
draw_panel(pdf_file, "pdf")

cat("Generated:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
cat(file.path(out_dir, "main_four_function_rep001_per_beta_color_ranges.csv"), "\n")
