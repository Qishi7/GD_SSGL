#!/usr/bin/env Rscript

options(digits = 16, warn = 1)

root_dir <- normalizePath(Sys.getenv("GDSSGL_ROOT", getwd()),
                          winslash = "/", mustWork = TRUE)
comparison_root <- file.path(root_dir, "new_method_comparison")
out_root <- file.path(
  comparison_root, "real_data_analysis",
  "modis_n10000_intercept_svd_centered_y_check"
)
raw_run_root <- Sys.getenv(
  "MODIS_SPLIT_RUN_ROOT",
  file.path(comparison_root, "real_data_analysis",
            "modis_n10000_intercept_svd_current_methods_run")
)
modis_path <- Sys.getenv("MODIS_RDATA",
                         file.path(root_dir, "data",
                                   "data_cleaned_small_expanded.RData"))
input_tag <- Sys.getenv("PY_MGWR_INPUT_TAG", "grouped_lc")
lc_encoding <- Sys.getenv("MODIS_PY_MGWR_LC_ENCODING", "grouped_numeric")
if (!lc_encoding %in% c("grouped_numeric", "dummy")) {
  stop("Unsupported MODIS_PY_MGWR_LC_ENCODING: ", lc_encoding)
}

dir.create(file.path(out_root, "data"), recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
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

continuous_vars <- c(
  "red_reflectance", "NIR_reflectance", "blue_reflectance",
  "MIR_reflectance", "GPP", "LE", "view_zenith_angle",
  "sun_zenith_angle", "relative_azimuth_angle"
)
required_vars <- c("EVI", "scaled_x", "scaled_y", "LC_Type4", continuous_vars)

split_path <- file.path(raw_run_root, "data",
                        "modis_n10000_intercept_svd_train_test_split.csv")
if (!file.exists(split_path)) stop("Missing frozen split file: ", split_path)
if (!file.exists(modis_path)) stop("Missing MODIS RData. Set MODIS_RDATA: ", modis_path)

split_df <- utils::read.csv(split_path)
load(modis_path)
if (!exists("data_cleaned_small")) {
  stop("MODIS RData must contain object `data_cleaned_small`.")
}

complete <- stats::complete.cases(data_cleaned_small[, required_vars]) &
  is.finite(data_cleaned_small$EVI) & data_cleaned_small$EVI > -1
real_data_clean <- data_cleaned_small[complete, ]
sample_data <- real_data_clean[split_df$original_clean_row, ]
train_indices <- which(split_df$split == "train")
test_indices <- which(split_df$split == "test")

train_y_log <- log(sample_data$EVI[train_indices] + 1)
test_y_log <- log(sample_data$EVI[test_indices] + 1)
y_train_mean <- mean(train_y_log)
train_y_centered <- train_y_log - y_train_mean
test_y_centered <- test_y_log - y_train_mean

if (lc_encoding == "grouped_numeric") {
  train_X_raw <- cbind(
    as.matrix(sample_data[train_indices, continuous_vars]),
    LC_Type4 = as.numeric(factor(sample_data$LC_Type4[train_indices]))
  )
  test_X_raw <- cbind(
    as.matrix(sample_data[test_indices, continuous_vars]),
    LC_Type4 = as.numeric(factor(sample_data$LC_Type4[test_indices],
                                 levels = levels(factor(sample_data$LC_Type4[train_indices]))))
  )
} else {
  lc <- make_lc_dummy_design(sample_data$LC_Type4[train_indices],
                             sample_data$LC_Type4[test_indices])
  train_X_raw <- cbind(as.matrix(sample_data[train_indices, continuous_vars]), lc$train)
  test_X_raw <- cbind(as.matrix(sample_data[test_indices, continuous_vars]), lc$test)
}

scaled_X <- scale_train_test(train_X_raw, test_X_raw)
train_df <- data.frame(
  y_centered = train_y_centered,
  y_log = train_y_log,
  EVI = sample_data$EVI[train_indices],
  s1 = sample_data$scaled_x[train_indices],
  s2 = sample_data$scaled_y[train_indices],
  scaled_X$train,
  check.names = FALSE
)
test_df <- data.frame(
  y_centered = test_y_centered,
  y_log = test_y_log,
  EVI = sample_data$EVI[test_indices],
  s1 = sample_data$scaled_x[test_indices],
  s2 = sample_data$scaled_y[test_indices],
  scaled_X$test,
  check.names = FALSE
)
metadata <- data.frame(
  key = c("y_train_mean", "n_train", "n_test", "p", "lc_encoding"),
  value = c(y_train_mean, nrow(train_df), nrow(test_df), ncol(scaled_X$train),
            if (lc_encoding == "grouped_numeric") {
              "LC_Type4 numeric factor code standardized by train mean/sd"
            } else {
              "LC_Type4 dummy variables standardized by train mean/sd"
            }),
  stringsAsFactors = FALSE
)

mid <- if (nzchar(input_tag)) paste0("_", input_tag) else ""
train_path <- file.path(out_root, "data",
                        paste0("modis_centered_y_mgwr_python", mid, "_train.csv"))
test_path <- file.path(out_root, "data",
                       paste0("modis_centered_y_mgwr_python", mid, "_test.csv"))
metadata_path <- file.path(out_root, "data",
                           paste0("modis_centered_y_mgwr_python", mid, "_metadata.csv"))
write_csv(train_df, train_path)
write_csv(test_df, test_path)
write_csv(metadata, metadata_path)

cat("Saved Python MGWR inputs:\n")
cat("- train: ", train_path, "\n", sep = "")
cat("- test: ", test_path, "\n", sep = "")
cat("- metadata: ", metadata_path, "\n", sep = "")
