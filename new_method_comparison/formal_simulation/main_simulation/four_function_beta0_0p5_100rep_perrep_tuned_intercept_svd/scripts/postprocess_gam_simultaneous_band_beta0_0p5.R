#!/usr/bin/env Rscript

options(digits = 17, warn = 1)

suppressPackageStartupMessages({
  library(mgcv)
})

project_root <- normalizePath(Sys.getenv("GDSSGL_ROOT", getwd()), winslash = "/", mustWork = TRUE)
root <- file.path(project_root, "new_method_comparison", "formal_simulation", "main_simulation", "four_function_beta0_0p5_100rep_perrep_tuned_intercept_svd")
out_root <- file.path(root, "gam_simultaneous_band_beta0_0p5_intercept")
dir.create(file.path(out_root, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "reports"), recursive = TRUE, showWarnings = FALSE)

p <- 10L
R <- 100L
B <- 1000L
seed_base <- 2026083101L
predictor_names <- paste0("X", seq_len(p))
predictor_group <- c("global_only", "global_only", "spatial_only",
                     "spatial_only", "global_plus_spatial",
                     "global_plus_spatial", rep("null", 4))

summary_vec <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  s <- if (n > 1L) stats::sd(x) else NA_real_
  data.frame(n = n, mean = mean(x), sd = s,
             mcse = if (n > 1L) s / sqrt(n) else NA_real_,
             median = stats::median(x), iqr = stats::IQR(x),
             min = min(x), max = max(x))
}

summarise_split <- function(df, groups, value_cols) {
  key <- if (length(groups)) interaction(df[groups], drop = TRUE, lex.order = TRUE) else
    factor(rep("overall", nrow(df)))
  rows <- list()
  for (idx in split(seq_len(nrow(df)), key)) {
    for (v in value_cols) {
      head <- if (length(groups)) df[idx[1], groups, drop = FALSE] else
        data.frame(level = "overall", subgroup = "overall")
      rows[[length(rows) + 1L]] <- cbind(head, metric = v, summary_vec(df[[v]][idx]))
    }
  }
  do.call(rbind, rows)
}

quad_se <- function(L, Vp) sqrt(pmax(0, rowSums((L %*% Vp) * L)))

coef_noise <- function(Vp, B, seed) {
  set.seed(seed)
  q <- ncol(Vp)
  z <- matrix(stats::rnorm(B * q), nrow = B, ncol = q)
  Rchol <- tryCatch(chol(Vp), error = function(e) NULL)
  if (!is.null(Rchol)) return(t(z %*% Rchol))
  ee <- eigen((Vp + t(Vp)) / 2, symmetric = TRUE)
  vals <- pmax(ee$values, 0)
  t(z %*% (ee$vectors %*% diag(sqrt(vals), length(vals), length(vals))))
}

make_grid <- function(coords) {
  gd <- data.frame(s1 = coords[, 1], s2 = coords[, 2])
  for (j in seq_len(p)) gd[[predictor_names[j]]] <- 0
  gd
}

process_rep <- function(r) {
  rep_dir <- file.path(root, "fits", sprintf("rep_%03d", r))
  dat <- readRDS(file.path(rep_dir, "shared_dataset_beta0_0p5.rds"))
  gf <- readRDS(file.path(rep_dir, "gam_final_fit.rds"))
  fit <- if (inherits(gf, "gam")) gf else gf$fit
  if (is.null(fit$Vp)) stop("No Vp in GAM fit for replicate ", r)
  cf <- stats::coef(fit)
  Vp <- fit$Vp
  base_grid <- make_grid(dat$grid$coords)
  L0 <- predict(fit, newdata = base_grid, type = "lpmatrix")
  noise <- coef_noise(Vp, B, seed_base + r)
  beta_rows <- vector("list", p)
  sim_rows <- vector("list", p)
  for (j in seq_len(p)) {
    gd <- base_grid
    gd[[predictor_names[j]]] <- 1
    L1 <- predict(fit, newdata = gd, type = "lpmatrix")
    Lb <- L1 - L0
    bhat <- as.numeric(Lb %*% cf)
    bse <- quad_se(Lb, Vp)
    blo <- bhat - 1.96 * bse
    bhi <- bhat + 1.96 * bse
    truth_beta <- dat$grid$true_beta[, j]
    beta_draws_centered <- Lb %*% noise
    max_std <- apply(abs(beta_draws_centered / pmax(bse, 1e-12)), 2, max)
    crit <- stats::quantile(max_std, 0.95, names = FALSE)
    slo <- bhat - crit * bse
    shi <- bhat + crit * bse
    beta_rows[[j]] <- data.frame(
      replicate = r,
      predictor = predictor_names[j],
      predictor_group = predictor_group[j],
      pointwise_grid_coverage = mean(truth_beta >= blo & truth_beta <= bhi),
      pointwise_full_surface_coverage = all(truth_beta >= blo & truth_beta <= bhi),
      mean_pointwise_band_width = mean(bhi - blo),
      stringsAsFactors = FALSE
    )
    sim_rows[[j]] <- data.frame(
      replicate = r,
      predictor = predictor_names[j],
      predictor_group = predictor_group[j],
      simultaneous_coverage = all(truth_beta >= slo & truth_beta <= shi),
      mean_simultaneous_band_width = mean(shi - slo),
      c_0.95 = crit,
      stringsAsFactors = FALSE
    )
  }
  beta <- do.call(rbind, beta_rows)
  sim <- do.call(rbind, sim_rows)
  overall <- data.frame(
    replicate = r,
    pointwise_grid_coverage = mean(beta$pointwise_grid_coverage),
    pointwise_full_surface_coverage = mean(beta$pointwise_full_surface_coverage),
    simultaneous_coverage = mean(sim$simultaneous_coverage),
    mean_pointwise_band_width = mean(beta$mean_pointwise_band_width),
    mean_simultaneous_band_width = mean(sim$mean_simultaneous_band_width),
    mean_c_0.95 = mean(sim$c_0.95),
    stringsAsFactors = FALSE
  )
  group_rows <- merge(beta, sim, by = c("replicate", "predictor", "predictor_group"))
  group <- do.call(rbind, lapply(split(group_rows, group_rows$predictor_group), function(x) {
    data.frame(
      replicate = r,
      predictor_group = x$predictor_group[1],
      pointwise_grid_coverage = mean(x$pointwise_grid_coverage),
      pointwise_full_surface_coverage = mean(x$pointwise_full_surface_coverage),
      simultaneous_coverage = mean(x$simultaneous_coverage),
      mean_pointwise_band_width = mean(x$mean_pointwise_band_width),
      mean_simultaneous_band_width = mean(x$mean_simultaneous_band_width),
      mean_c_0.95 = mean(x$c_0.95),
      stringsAsFactors = FALSE
    )
  }))
  list(beta = beta, sim = sim, overall = overall, group = group,
       audit = data.frame(
         replicate = r,
         formula = paste(deparse(formula(fit)), collapse = " "),
         method = fit$method,
         has_intercept = "(Intercept)" %in% names(cf),
         intercept_hat = unname(cf["(Intercept)"]),
         coef_n = length(cf),
         warnings = length(fit$warnings %||% character()),
         stringsAsFactors = FALSE
       ))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

objects <- vector("list", R)
for (r in seq_len(R)) {
  message(sprintf("GAM beta0=0.5 simultaneous band postprocess: rep %03d/%03d", r, R))
  objects[[r]] <- process_rep(r)
}

bind_field <- function(field) do.call(rbind, lapply(objects, `[[`, field))
beta <- bind_field("beta")
sim <- bind_field("sim")
overall_by_rep <- bind_field("overall")
group_by_rep <- bind_field("group")
audit <- bind_field("audit")
replicate_level <- merge(beta, sim, by = c("replicate", "predictor", "predictor_group"))

overall_summary <- summarise_split(
  overall_by_rep, character(0),
  c("simultaneous_coverage", "pointwise_grid_coverage",
    "pointwise_full_surface_coverage", "mean_simultaneous_band_width",
    "mean_pointwise_band_width", "mean_c_0.95")
)
group_summary <- summarise_split(
  group_by_rep, "predictor_group",
  c("simultaneous_coverage", "pointwise_grid_coverage",
    "pointwise_full_surface_coverage", "mean_simultaneous_band_width",
    "mean_pointwise_band_width", "mean_c_0.95")
)
predictor_summary <- summarise_split(
  replicate_level, c("predictor", "predictor_group"),
  c("simultaneous_coverage", "pointwise_grid_coverage",
    "pointwise_full_surface_coverage", "mean_simultaneous_band_width",
    "mean_pointwise_band_width", "c_0.95")
)

write.csv(beta, file.path(out_root, "data", "gam_beta0_0p5_pointwise_by_predictor_replicate.csv"),
          row.names = FALSE)
write.csv(sim, file.path(out_root, "data", "gam_beta0_0p5_simband_by_predictor_replicate.csv"),
          row.names = FALSE)
write.csv(overall_by_rep, file.path(out_root, "data", "gam_beta0_0p5_simband_overall_by_replicate.csv"),
          row.names = FALSE)
write.csv(group_by_rep, file.path(out_root, "data", "gam_beta0_0p5_simband_group_by_replicate.csv"),
          row.names = FALSE)
write.csv(replicate_level, file.path(out_root, "data", "gam_beta0_0p5_simband_replicate_level.csv"),
          row.names = FALSE)
write.csv(overall_summary, file.path(out_root, "data", "gam_beta0_0p5_simband_overall_summary.csv"),
          row.names = FALSE)
write.csv(group_summary, file.path(out_root, "data", "gam_beta0_0p5_simband_group_summary.csv"),
          row.names = FALSE)
write.csv(predictor_summary, file.path(out_root, "data", "gam_beta0_0p5_simband_predictor_summary.csv"),
          row.names = FALSE)
write.csv(audit, file.path(out_root, "data", "gam_beta0_0p5_simband_audit.csv"),
          row.names = FALSE)

report <- c(
  "# GAM simultaneous-band diagnostic for beta0=0.5 intercept-version run",
  "",
  paste0("- Source root: `", root, "`"),
  paste0("- Output root: `", out_root, "`"),
  "- Used saved GAM fit objects; no models were refit.",
  "- Each saved GAM fit contains an mgcv covariance matrix `Vp`.",
  sprintf("- Replicates processed: %d / 100", length(unique(audit$replicate))),
  sprintf("- All fits have intercept: %s", all(audit$has_intercept)),
  "- Coefficient surfaces are reconstructed by prediction contrast: X_j=1 minus all X=0, so the scalar intercept cancels out.",
  "- Simultaneous bands use coefficient simulation from mgcv Vp with B=1000.",
  "",
  "## Overall summary",
  "",
  paste(capture.output(print(overall_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Group summary",
  "",
  paste(capture.output(print(group_summary, row.names = FALSE)), collapse = "\n")
)
writeLines(report, file.path(out_root, "reports",
                             "gam_beta0_0p5_simultaneous_band_report.md"))
cat(paste(report, collapse = "\n"))
