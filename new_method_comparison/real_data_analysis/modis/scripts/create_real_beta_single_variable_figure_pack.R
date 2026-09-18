#!/usr/bin/env Rscript

options(digits = 16, warn = 1)

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
})

root_dir <- normalizePath(Sys.getenv("GDSSGL_ROOT", getwd()), winslash = "/", mustWork = TRUE)
comparison_root <- file.path(root_dir, "new_method_comparison")
run_root <- file.path(
  comparison_root, "real_data_analysis",
  "modis_n10000_intercept_svd_centered_y_cv_tuned_correct"
)
raw_root <- file.path(
  comparison_root, "real_data_analysis",
  "modis_n10000_intercept_svd_current_methods_run"
)
out_root <- file.path(
  comparison_root, "real_data_analysis",
  "modis_n10000_intercept_svd_centered_y_check",
  "real_beta_single_variable_figure_pack"
)
fig_dir <- file.path(out_root, "figures")
data_dir <- file.path(out_root, "data")
report_dir <- file.path(out_root, "reports")
script_dir <- file.path(out_root, "scripts")
invisible(lapply(c(fig_dir, data_dir, report_dir, script_dir),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

continuous_vars <- c(
  "red_reflectance", "NIR_reflectance", "blue_reflectance",
  "MIR_reflectance", "GPP", "LE", "view_zenith_angle",
  "sun_zenith_angle", "relative_azimuth_angle"
)
internal_vars <- c(continuous_vars, paste0("LC_Type4_", c(2, 4, 5, 6, 7, 8)))
group_vars <- c(continuous_vars, "LC_Type4")
display_labels <- c(
  red_reflectance = "Red reflectance",
  NIR_reflectance = "NIR reflectance",
  blue_reflectance = "Blue reflectance",
  MIR_reflectance = "MIR reflectance",
  GPP = "GPP",
  LE = "LE",
  view_zenith_angle = "View zenith angle",
  sun_zenith_angle = "Sun zenith angle",
  relative_azimuth_angle = "Relative azimuth angle",
  LC_Type4 = "LC Type4"
)
method_order <- c("GD-SSGL", "WS-SSGL", "Gaussian SVC", "MGWR", "Bayesian Lasso")
blue_green_cols <- grDevices::colorRampPalette(c("#253494", "#2C7FB8", "#41B6C4", "#2CA25F", "#FDE725"))(256)

atomic_write_csv <- function(x, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Could not write ", path)
  invisible(path)
}

scale_train_test <- function(train_x, test_x) {
  mu <- colMeans(train_x, na.rm = TRUE)
  sig <- apply(train_x, 2, stats::sd, na.rm = TRUE)
  sig[!is.finite(sig) | sig == 0] <- 1
  list(
    train = sweep(sweep(train_x, 2, mu, "-"), 2, sig, "/"),
    test = sweep(sweep(test_x, 2, mu, "-"), 2, sig, "/"),
    mean = mu,
    sd = sig
  )
}

make_lc_dummy_design <- function(train_lc, test_lc) {
  train_levels <- sort(unique(train_lc))
  ref_level <- train_levels[1]
  dummy_levels <- train_levels[train_levels != ref_level]
  train_dummy <- sapply(dummy_levels, function(level) as.numeric(train_lc == level))
  test_dummy <- sapply(dummy_levels, function(level) as.numeric(test_lc == level))
  if (length(dummy_levels) == 1) {
    train_dummy <- matrix(train_dummy, ncol = 1)
    test_dummy <- matrix(test_dummy, ncol = 1)
  }
  colnames(train_dummy) <- paste0("LC_Type4_", dummy_levels)
  colnames(test_dummy) <- paste0("LC_Type4_", dummy_levels)
  list(train = train_dummy, test = test_dummy, reference = ref_level,
       levels = train_levels, dummy_levels = dummy_levels)
}

predict_local_wls_coef <- function(train_X, train_y, train_coords, new_X, new_coords,
                                   k = 1000L, ridge = 1e-8) {
  train_X <- as.matrix(train_X)
  new_X <- as.matrix(new_X)
  train_coords <- as.matrix(train_coords)
  new_coords <- as.matrix(new_coords)
  n_train <- nrow(train_X)
  p <- ncol(train_X)
  k <- max(p + 2L, min(as.integer(k), n_train - 1L))
  X_design <- cbind(Intercept = 1, train_X)
  coef_new <- matrix(NA_real_, nrow = nrow(new_coords), ncol = p + 1L)
  colnames(coef_new) <- colnames(X_design)
  for (i in seq_len(nrow(new_coords))) {
    d <- sqrt(rowSums((train_coords - matrix(new_coords[i, ],
                                             nrow(train_coords),
                                             ncol(train_coords),
                                             byrow = TRUE))^2))
    kth <- sort(d, partial = k)[k]
    if (!is.finite(kth) || kth <= 0) kth <- max(d)
    idx <- which(d <= kth)
    if (length(idx) < p + 2L) idx <- order(d)[seq_len(min(n_train, p + 2L))]
    di <- d[idx]
    bw <- max(di)
    if (!is.finite(bw) || bw <= 0) bw <- max(d)
    wi <- (1 - (di / bw)^2)^2
    wi[di > bw] <- 0
    if (!any(wi > 0)) wi[] <- 1
    Xi <- X_design[idx, , drop = FALSE]
    yi <- train_y[idx]
    sw <- sqrt(wi)
    Xw <- Xi * sw
    yw <- yi * sw
    XtX <- crossprod(Xw)
    diag(XtX) <- diag(XtX) + ridge
    beta <- tryCatch(
      solve(XtX, crossprod(Xw, yw)),
      error = function(e) qr.solve(XtX, crossprod(Xw, yw))
    )
    coef_new[i, ] <- as.numeric(beta)
  }
  coef_new
}

signed_max_abs_group <- function(mat, cols) {
  vals <- as.matrix(mat[, cols, drop = FALSE])
  vals[cbind(seq_len(nrow(vals)), max.col(abs(vals), ties.method = "first"))]
}

surface_values <- function(mat, variable) {
  if (variable == "LC_Type4") {
    if ("LC_Type4" %in% colnames(mat)) return(mat[, "LC_Type4"])
    signed_max_abs_group(mat, grep("^LC_Type4_", colnames(mat), value = TRUE))
  } else {
    mat[, variable]
  }
}

idw_predict_grid <- function(train_coords, train_values, grid_coords, k = 160L,
                             power = 2, eps = 1e-10) {
  train_coords <- as.matrix(train_coords)
  grid_coords <- as.matrix(grid_coords)
  train_values <- as.matrix(train_values)
  k <- max(5L, min(as.integer(k), nrow(train_coords)))
  pred <- matrix(NA_real_, nrow = nrow(grid_coords), ncol = ncol(train_values))
  for (i in seq_len(nrow(grid_coords))) {
    d <- sqrt(rowSums((train_coords - matrix(grid_coords[i, ],
                                             nrow(train_coords),
                                             ncol(train_coords),
                                             byrow = TRUE))^2))
    idx <- order(d)[seq_len(k)]
    if (d[idx[1]] < eps) {
      pred[i, ] <- train_values[idx[1], ]
    } else {
      w <- 1 / pmax(d[idx], eps)^power
      w <- w / sum(w)
      pred[i, ] <- as.numeric(crossprod(w, train_values[idx, , drop = FALSE]))
    }
  }
  colnames(pred) <- colnames(train_values)
  pred
}

read_ne_zip <- function(file_name) {
  zip_path <- normalizePath(file.path(Sys.getenv("NATURALEARTH_DIR", file.path(root_dir, "data", "geo_naturalearth_10m")), file_name))
  shp_name <- sub("\\.zip$", ".shp", file_name)
  st_read(file.path("/vsizip", zip_path, shp_name), quiet = TRUE)
}

scaled_to_lonlat <- function(df, modis_all_data) {
  x_original <- df$s1 * 1e6 + min(modis_all_data$x)
  y_original <- df$s2 * 1e6 + min(modis_all_data$y)
  modis_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
  pts <- st_as_sf(data.frame(x = x_original, y = y_original), coords = c("x", "y"), crs = modis_crs)
  xy <- st_coordinates(st_transform(pts, 4326))
  df$lon <- xy[, 1]
  df$lat <- xy[, 2]
  df
}

upsample_surface_lonlat <- function(df, modis_all_data, n_fine = 180L) {
  xs <- sort(unique(df$s1))
  ys <- sort(unique(df$s2))
  z <- matrix(NA_real_, nrow = length(ys), ncol = length(xs))
  z[cbind(match(df$s2, ys), match(df$s1, xs))] <- df$value
  fine_x <- seq(min(xs), max(xs), length.out = n_fine)
  fine_y <- seq(min(ys), max(ys), length.out = n_fine)
  z_x <- t(apply(z, 1, function(row) approx(xs, row, xout = fine_x, rule = 2)$y))
  z_xy <- apply(z_x, 2, function(col) approx(ys, col, xout = fine_y, rule = 2)$y)
  out <- expand.grid(s1 = fine_x, s2 = fine_y)
  out$value <- as.vector(t(z_xy))
  out <- scaled_to_lonlat(out, modis_all_data)
  row_widths <- unlist(tapply(out$lon, round(out$s2, 8), function(v) diff(sort(v))))
  out$tile_width <- median(row_widths, na.rm = TRUE) * 1.10
  out$tile_height <- median(diff(sort(unique(round(out$lat, 6)))), na.rm = TRUE) * 1.10
  out
}

load(Sys.getenv("MODIS_RDATA", file.path(root_dir, "data", "data_cleaned_small_expanded.RData")))
required_vars <- c("EVI", "scaled_x", "scaled_y", "LC_Type4", continuous_vars)
complete <- stats::complete.cases(data_cleaned_small[, required_vars]) &
  is.finite(data_cleaned_small$EVI) & data_cleaned_small$EVI > -1
real_data_clean <- data_cleaned_small[complete, ]

split_df <- read.csv(file.path(raw_root, "data/modis_n10000_intercept_svd_train_test_split.csv"))
sample_data <- real_data_clean[split_df$original_clean_row, ]
train_indices <- which(split_df$split == "train")
test_indices <- which(split_df$split == "test")

train_coords <- as.matrix(sample_data[train_indices, c("scaled_x", "scaled_y")])
colnames(train_coords) <- c("s1", "s2")
base_grid <- expand.grid(
  s1 = seq(min(train_coords[, 1]), max(train_coords[, 1]), length.out = 40L),
  s2 = seq(min(train_coords[, 2]), max(train_coords[, 2]), length.out = 40L)
)

lc <- make_lc_dummy_design(sample_data$LC_Type4[train_indices], sample_data$LC_Type4[test_indices])
train_X_raw <- cbind(as.matrix(sample_data[train_indices, continuous_vars]), lc$train)
test_X_raw <- cbind(as.matrix(sample_data[test_indices, continuous_vars]), lc$test)
scaled_X <- scale_train_test(train_X_raw, test_X_raw)
train_y_raw <- log(sample_data$EVI[train_indices] + 1)
train_y_centered <- train_y_raw - mean(train_y_raw)

gd <- readRDS(file.path(run_root, "results/gdssgl_correct_cv_tuned_final_fit.rds"))
ws <- readRDS(file.path(run_root, "results/wsssgl_correct_cv_tuned_final_fit.rds"))
svc <- readRDS(file.path(run_root, "results/gaussian_svc_correct_cv_tuned_final_fit.rds"))
bl <- readRDS(file.path(run_root, "results/bayesian_lasso_correct_final_fit.rds"))

beta_mats <- list(
  "GD-SSGL" = gd$beta_hat_grid,
  "WS-SSGL" = ws$beta_hat_grid,
  "Gaussian SVC" = svc$beta_hat_grid,
  "Bayesian Lasso" = bl$beta_hat_grid
)
for (nm in names(beta_mats)) colnames(beta_mats[[nm]]) <- internal_vars

mgwr_param_path <- file.path(
  run_root, "data",
  "modis_correct_cv_tuned_mgwr_grouped_lc_train_params.csv"
)
mgwr_params <- read.csv(mgwr_param_path, check.names = FALSE)
mgwr_beta_grid <- idw_predict_grid(
  train_coords,
  as.matrix(mgwr_params[, setdiff(colnames(mgwr_params), "Intercept"), drop = FALSE]),
  as.matrix(base_grid)
)
beta_mats[["MGWR"]] <- mgwr_beta_grid
beta_mats <- beta_mats[method_order]

sf_use_s2(FALSE)
study_area <- st_as_sfc(st_bbox(c(xmin = -125, xmax = -103, ymin = 29, ymax = 41), crs = st_crs(4326)))
geo_layers <- list(
  ocean = st_intersection(read_ne_zip("ne_10m_ocean.zip"), study_area),
  lakes = st_intersection(read_ne_zip("ne_10m_lakes.zip"), study_area),
  rivers = st_intersection(read_ne_zip("ne_10m_rivers_lake_centerlines.zip"), study_area),
  coast = st_intersection(read_ne_zip("ne_10m_coastline.zip"), study_area)
)
sf_use_s2(TRUE)

make_beta_df <- function(variable) {
  out <- do.call(rbind, lapply(names(beta_mats), function(m) {
    df <- base_grid
    df$value <- surface_values(beta_mats[[m]], variable)
    ff <- upsample_surface_lonlat(df, data_cleaned_small)
    ff$method <- m
    ff
  }))
  out$method <- factor(out$method, levels = method_order)
  out$variable <- variable
  out$label <- unname(display_labels[[variable]])
  out
}

theme_real_beta <- function(base_size = 20) {
  theme_void(base_size = base_size) +
    theme(
      strip.text = element_text(face = "bold", size = base_size + 7, margin = margin(t = 0, b = -2)),
      legend.text = element_text(size = base_size + 1),
      legend.title = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(0, 1, 0, 1),
      panel.spacing = unit(0.025, "in")
    )
}

plot_one_variable <- function(variable) {
  df <- make_beta_df(variable)
  lim <- max(abs(df$value), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  label <- unname(display_labels[[variable]])
  p <- ggplot() +
    geom_sf(data = geo_layers$ocean, fill = "#C6DBEF", color = NA, alpha = 0.55) +
    geom_tile(data = df, aes(lon, lat, fill = value),
              width = unique(df$tile_width)[1],
              height = unique(df$tile_height)[1]) +
    geom_sf(data = geo_layers$lakes, fill = "#4292C6", color = "#2171B5", linewidth = 0.10) +
    geom_sf(data = geo_layers$rivers, color = "#2171B5", linewidth = 0.16, alpha = 0.55) +
    geom_sf(data = geo_layers$coast, color = "grey28", linewidth = 0.26, fill = NA) +
    coord_sf(xlim = c(-125, -103), ylim = c(29, 41), expand = FALSE) +
    facet_wrap(~ method, nrow = 1) +
    scale_fill_gradientn(
      colors = blue_green_cols,
      limits = c(-lim, lim),
      breaks = c(-lim, 0, lim),
      labels = function(x) {
        ax <- max(abs(x), na.rm = TRUE)
        if (ax < 0.01) sprintf("%.3f", x)
        else if (ax < 0.1) sprintf("%.2f", x)
        else if (ax < 1) sprintf("%.2f", x)
        else if (ax < 2) sprintf("%.1f", x)
        else sprintf("%.0f", x)
      },
      guide = guide_colorbar(
        barheight = unit(0.92, "in"),
        barwidth = unit(0.12, "in"),
        ticks = TRUE
      )
    ) +
    labs(title = NULL) +
    theme_real_beta(20) +
    theme(
      legend.position = "right",
      legend.margin = margin(0, 0, 0, -6),
      legend.box.margin = margin(0, 0, 0, -8)
    )

  png_path <- file.path(fig_dir, paste0("modis_centered_y_beta_", variable, "_no_gam.png"))
  pdf_path <- file.path(fig_dir, paste0("modis_centered_y_beta_", variable, "_no_gam.pdf"))
  ggsave(png_path, p, width = 18.2, height = 2.85, dpi = 260, bg = "white")
  if (file.exists(pdf_path)) unlink(pdf_path)
  data.frame(
    variable = variable,
    label = label,
    png = png_path,
    methods = paste(method_order, collapse = "; "),
    symmetric_abs_limit = lim,
    stringsAsFactors = FALSE
  )
}

inventory <- do.call(rbind, lapply(group_vars, plot_one_variable))
atomic_write_csv(inventory, file.path(data_dir, "real_beta_single_variable_figure_inventory.csv"))

writeLines(c(
  "# Real MODIS beta single-variable figure pack",
  "",
  paste0("- Output root: `", out_root, "`."),
  "- Protocol: centered-y, explicit-intercept, CV-tuned current intercept/SVD fits.",
  "- GAM is intentionally excluded.",
  "- LC Type4 is plotted as one grouped factor using the signed maximum-absolute dummy coefficient at each grid location.",
  "- MGWR uses the completed Python mgwr grouped-LC coefficients; train-location coefficients are interpolated to the common plotting grid.",
  "",
  "## Methods",
  paste0("- ", method_order),
  "",
  "## Figures",
  paste0("- `", inventory$png, "`")
), file.path(report_dir, "real_beta_single_variable_figure_pack_report.md"))

cat("Generated real-data beta figure pack:\n")
cat("Output root:", out_root, "\n")
print(inventory, row.names = FALSE)
