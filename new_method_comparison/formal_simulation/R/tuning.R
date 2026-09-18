if (!exists(".formal_root")) {
  source(file.path("new_method_comparison", "formal_simulation",
                   "R", "framework.R"))
}

make_historical_cv_folds <- function(n, k_folds = 5L, seed = 123L) {
  set.seed(seed)
  sample(rep(seq_len(k_folds), length.out = n))
}

cv_data_objects <- function(data, train_idx, validation_idx) {
  list(
    train = list(
      y = data$train$y[train_idx],
      X = data$train$X[train_idx, , drop = FALSE],
      coords = data$train$coords[train_idx, , drop = FALSE]
    ),
    validation = list(
      y = data$train$y[validation_idx],
      X = data$train$X[validation_idx, , drop = FALSE],
      coords = data$train$coords[validation_idx, , drop = FALSE]
    )
  )
}

cv_wrapper_grid <- function(data, method_function, param_grid, folds,
                            basis_dimension = 6L,
                            cv_n_iter = 500L, cv_burn_in = 100L,
                            fixed_model_config = list(),
                            seed_start = 740000L) {
  results <- list()
  row_id <- 1L
  started <- proc.time()[3]
  for (candidate_id in seq_len(nrow(param_grid))) {
    params <- as.list(param_grid[candidate_id, , drop = FALSE])
    for (fold in sort(unique(folds))) {
      validation_idx <- which(folds == fold)
      train_idx <- which(folds != fold)
      fold_data <- cv_data_objects(data, train_idx, validation_idx)
      model_config <- modifyList(fixed_model_config, params)
      model_config$b_sigma <- stats::var(fold_data$train$y) / 2
      fit_started <- proc.time()[3]
      fit <- tryCatch(
        method_function(
          fold_data$train,
          fold_data$validation,
          list(coords = fold_data$validation$coords),
          list(n_basis = basis_dimension),
          model_config,
          list(n_iter = cv_n_iter, burn_in = cv_burn_in),
          seed = as.integer(seed_start + fold)
        ),
        error = identity
      )
      success <- !inherits(fit, "error")
      results[[row_id]] <- data.frame(
        candidate_id = candidate_id,
        fold = fold,
        validation_mse = if (success) {
          mean((fold_data$validation$y - fit$pred_test)^2)
        } else {
          Inf
        },
        success = success,
        error = if (success) "" else conditionMessage(fit),
        runtime_seconds = proc.time()[3] - fit_started,
        stringsAsFactors = FALSE
      )
      for (name in names(params)) results[[row_id]][[name]] <- params[[name]]
      row_id <- row_id + 1L
    }
  }
  fold_results <- do.call(rbind, results)
  candidate_results <- do.call(rbind, lapply(
    split(fold_results, fold_results$candidate_id),
    function(x) {
      parameter_columns <- setdiff(
        names(x),
        c("candidate_id", "fold", "validation_mse", "success",
          "error", "runtime_seconds")
      )
      out <- data.frame(
        candidate_id = x$candidate_id[1],
        cv_mse = if (all(x$success)) mean(x$validation_mse) else Inf,
        cv_mse_sd = if (all(x$success)) sd(x$validation_mse) else NA_real_,
        successful_folds = sum(x$success),
        total_runtime_seconds = sum(x$runtime_seconds),
        stringsAsFactors = FALSE
      )
      for (name in parameter_columns) out[[name]] <- x[[name]][1]
      out
    }
  ))
  candidate_results <- candidate_results[
    order(candidate_results$candidate_id), ]
  best_index <- which.min(candidate_results$cv_mse)
  list(
    folds = folds,
    fold_results = fold_results,
    candidate_results = candidate_results,
    best = candidate_results[best_index, , drop = FALSE],
    runtime_seconds = proc.time()[3] - started
  )
}

make_gam_formula <- function(p, k) {
  linear_terms <- paste0("X", seq_len(p), collapse = " + ")
  smooth_terms <- paste0(
    "s(x, y_coord, bs='tp', k=", k, ", by=X",
    seq_len(p), ")",
    collapse = " + "
  )
  as.formula(paste0(
    "y ~ 0 + ", linear_terms, " + ", smooth_terms
  ))
}

fit_gam_with_k <- function(train, validation, k) {
  p <- ncol(train$X)
  train_df <- as.data.frame(train$X)
  names(train_df) <- paste0("X", seq_len(p))
  train_df$y <- train$y
  train_df$x <- train$coords[, 1]
  train_df$y_coord <- train$coords[, 2]
  fit <- mgcv::gam(
    make_gam_formula(p, k), data = train_df,
    method = "GCV.Cp"
  )
  validation_df <- as.data.frame(validation$X)
  names(validation_df) <- paste0("X", seq_len(p))
  validation_df$x <- validation$coords[, 1]
  validation_df$y_coord <- validation$coords[, 2]
  list(
    fit = fit,
    pred = as.numeric(stats::predict(fit, newdata = validation_df))
  )
}

cv_gam_k <- function(data, k_values, folds) {
  rows <- list()
  row_id <- 1L
  started <- proc.time()[3]
  for (candidate_id in seq_along(k_values)) {
    k <- k_values[candidate_id]
    for (fold in sort(unique(folds))) {
      validation_idx <- which(folds == fold)
      train_idx <- which(folds != fold)
      fold_data <- cv_data_objects(data, train_idx, validation_idx)
      fit_started <- proc.time()[3]
      result <- tryCatch(
        fit_gam_with_k(fold_data$train, fold_data$validation, k),
        error = identity
      )
      success <- !inherits(result, "error")
      rows[[row_id]] <- data.frame(
        candidate_id = candidate_id, k = k, fold = fold,
        validation_mse = if (success) {
          mean((fold_data$validation$y - result$pred)^2)
        } else Inf,
        success = success,
        error = if (success) "" else conditionMessage(result),
        runtime_seconds = proc.time()[3] - fit_started,
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }
  fold_results <- do.call(rbind, rows)
  candidate_results <- do.call(rbind, lapply(
    split(fold_results, fold_results$candidate_id),
    function(x) {
      data.frame(
        candidate_id = x$candidate_id[1], k = x$k[1],
        cv_mse = if (all(x$success)) mean(x$validation_mse) else Inf,
        cv_mse_sd = if (all(x$success)) sd(x$validation_mse) else NA_real_,
        successful_folds = sum(x$success),
        total_runtime_seconds = sum(x$runtime_seconds),
        stringsAsFactors = FALSE
      )
    }
  ))
  candidate_results <- candidate_results[
    order(candidate_results$candidate_id), ]
  best_index <- which.min(candidate_results$cv_mse)
  list(
    folds = folds, fold_results = fold_results,
    candidate_results = candidate_results,
    best = candidate_results[best_index, , drop = FALSE],
    runtime_seconds = proc.time()[3] - started
  )
}

make_fixed_cv_folds <- function(n, k_folds = 5L, seed = 123L) {
  set.seed(seed)
  sample(rep(seq_len(k_folds), length.out = n))
}

fold_assignment_checksum <- function(folds) {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  write.csv(
    data.frame(row_id = seq_along(folds), fold = as.integer(folds)),
    tmp, row.names = FALSE
  )
  unname(tools::md5sum(tmp))
}

standardize_with_training_only <- function(x_train, x_validation) {
  x_mean <- colMeans(x_train)
  x_sd <- apply(x_train, 2, stats::sd)
  if (any(!is.finite(x_sd)) || any(x_sd <= 0)) {
    stop("Non-positive or non-finite fold-training predictor SD.")
  }
  list(
    train = sweep(sweep(x_train, 2, x_mean, "-"), 2, x_sd, "/"),
    validation = sweep(sweep(x_validation, 2, x_mean, "-"), 2, x_sd, "/"),
    x_mean = x_mean,
    x_sd = x_sd
  )
}

cv_data_objects_leakage_free <- function(data, folds, fold) {
  validation_idx <- which(folds == fold)
  train_idx <- which(folds != fold)
  scaled <- standardize_with_training_only(
    data$train$X[train_idx, , drop = FALSE],
    data$train$X[validation_idx, , drop = FALSE]
  )
  list(
    train = list(
      y = data$train$y[train_idx],
      X = scaled$train,
      coords = data$train$coords[train_idx, , drop = FALSE],
      row_id = train_idx
    ),
    validation = list(
      y = data$train$y[validation_idx],
      X = scaled$validation,
      coords = data$train$coords[validation_idx, , drop = FALSE],
      row_id = validation_idx
    ),
    preprocessing = list(
      x_mean = scaled$x_mean,
      x_sd = scaled$x_sd,
      estimated_from_rows = train_idx,
      applied_to_validation_rows = validation_idx
    )
  )
}

make_fold_basis_metadata <- function(fold_data, K, centered) {
  meta <- fit_basis_metadata(
    fold_data$train$coords,
    n_basis = K,
    centered = centered
  )
  phi_train <- apply_basis_metadata(fold_data$train$coords, meta)
  phi_validation <- apply_basis_metadata(fold_data$validation$coords, meta)
  list(
    metadata = meta,
    train_center_error = if (centered) max(abs(colMeans(phi_train))) else NA_real_,
    validation_nrow = nrow(phi_validation),
    estimated_from_rows = fold_data$train$row_id,
    applied_to_validation_rows = fold_data$validation$row_id
  )
}

make_external_cv_preprocessing_metadata <- function(data, folds, K_values,
                                                    methods) {
  rows <- list()
  metadata <- list()
  row_id <- 1L
  checksum <- fold_assignment_checksum(folds)
  for (method in names(methods)) {
    centered <- isTRUE(methods[[method]]$centered_basis)
    uses_basis <- isTRUE(methods[[method]]$uses_b_spline_basis)
    candidate_K <- if (uses_basis) K_values else NA_integer_
    for (fold in sort(unique(folds))) {
      fold_data <- cv_data_objects_leakage_free(data, folds, fold)
      for (K in candidate_K) {
        basis <- if (uses_basis) {
          make_fold_basis_metadata(fold_data, K, centered)
        } else {
          NULL
        }
        key <- paste(method, fold, ifelse(is.na(K), "no_basis", K), sep = "__")
        metadata[[key]] <- list(
          method = method,
          fold = fold,
          K = K,
          H = if (is.na(K)) NA_integer_ else K^2L,
          centered_basis = centered,
          train_rows = fold_data$train$row_id,
          validation_rows = fold_data$validation$row_id,
          preprocessing = fold_data$preprocessing,
          basis = basis
        )
        rows[[row_id]] <- data.frame(
          method = method,
          fold = fold,
          n_train = nrow(fold_data$train$X),
          n_validation = nrow(fold_data$validation$X),
          fold_assignment_checksum = checksum,
          training_predictor_means = paste(
            format(fold_data$preprocessing$x_mean, digits = 10),
            collapse = ";"
          ),
          training_predictor_sds = paste(
            format(fold_data$preprocessing$x_sd, digits = 10),
            collapse = ";"
          ),
          K = K,
          H = if (is.na(K)) NA_integer_ else K^2L,
          basis_centering_used = centered,
          basis_train_center_error = if (uses_basis && centered) {
            basis$train_center_error
          } else {
            NA_real_
          },
          validation_rows = paste(fold_data$validation$row_id, collapse = ";"),
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1L
      }
    }
  }
  list(audit_table = do.call(rbind, rows), metadata = metadata)
}

check_external_cv_preprocessing <- function(metadata_bundle, folds, methods,
                                            tolerance = 1e-8) {
  audit <- metadata_bundle$audit_table
  metadata <- metadata_bundle$metadata
  checks <- list()
  add_check <- function(name, passed, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name,
      passed = isTRUE(passed),
      detail = detail,
      stringsAsFactors = FALSE
    )
  }

  expected_checksum <- fold_assignment_checksum(folds)
  add_check(
    "all_rows_have_identical_fold_checksum",
    all(audit$fold_assignment_checksum == expected_checksum),
    expected_checksum
  )

  fold_validation_by_method <- unique(audit[, c("method", "fold", "validation_rows")])
  for (fold in sort(unique(folds))) {
    z <- fold_validation_by_method[fold_validation_by_method$fold == fold, ]
    add_check(
      paste0("same_validation_rows_fold_", fold),
      length(unique(z$validation_rows)) == 1L,
      unique(z$validation_rows)[1]
    )
  }

  for (key in names(metadata)) {
    item <- metadata[[key]]
    overlap <- intersect(
      item$preprocessing$estimated_from_rows,
      item$preprocessing$applied_to_validation_rows
    )
    add_check(
      paste0("no_scaling_validation_rows__", key),
      length(overlap) == 0L,
      if (length(overlap)) paste(overlap, collapse = ",") else ""
    )
    if (!is.null(item$basis)) {
      basis_overlap <- intersect(
        item$basis$estimated_from_rows,
        item$basis$applied_to_validation_rows
      )
      add_check(
        paste0("no_basis_validation_rows__", key),
        length(basis_overlap) == 0L,
        if (length(basis_overlap)) paste(basis_overlap, collapse = ",") else ""
      )
      if (isTRUE(item$centered_basis)) {
        add_check(
          paste0("centered_basis_train_mean_zero__", key),
          is.finite(item$basis$train_center_error) &&
            item$basis$train_center_error < tolerance,
          format(item$basis$train_center_error, digits = 12)
        )
      }
    }
  }

  method_names <- names(methods)
  add_check(
    "all_externally_tuned_methods_present",
    all(method_names %in% unique(audit$method)),
    paste(method_names, collapse = ",")
  )

  do.call(rbind, checks)
}

expand_basis_hyperparameter_grid <- function(K_values, hyper_grid) {
  out <- do.call(rbind, lapply(K_values, function(K) {
    cbind(
      data.frame(K = as.integer(K), H = as.integer(K)^2L),
      hyper_grid,
      row.names = NULL
    )
  }))
  rownames(out) <- NULL
  out
}

select_best_cv_candidate <- function(candidate_results) {
  ok <- candidate_results[is.finite(candidate_results$cv_mse), , drop = FALSE]
  if (!nrow(ok)) return(candidate_results[which.min(candidate_results$cv_mse), , drop = FALSE])
  min_mse <- min(ok$cv_mse)
  tied <- ok[abs(ok$cv_mse - min_mse) <=
               max(1e-12, 1e-10 * max(abs(min_mse), 1)), , drop = FALSE]
  if ("H" %in% names(tied)) {
    tied <- tied[order(tied$H, tied$candidate_id), , drop = FALSE]
  } else if ("K" %in% names(tied)) {
    tied <- tied[order(tied$K, tied$candidate_id), , drop = FALSE]
  } else {
    tied <- tied[order(tied$candidate_id), , drop = FALSE]
  }
  tied[1, , drop = FALSE]
}

cv_wrapper_grid_leakage_free <- function(data, method_function, param_grid,
                                         folds, cv_n_iter = 500L,
                                         cv_burn_in = 100L,
                                         fixed_model_config = list(),
                                         seed_start = 740000L,
                                         fold_metadata = NULL) {
  if (!("K" %in% names(param_grid))) {
    stop("param_grid must include K for leakage-free spatial CV.")
  }
  if (!identical(sort(unique(folds)), seq_len(length(unique(folds))))) {
    stop("Fold IDs must be consecutive positive integers.")
  }
  results <- list()
  row_id <- 1L
  started <- proc.time()[3]
  for (candidate_id in seq_len(nrow(param_grid))) {
    params_all <- as.list(param_grid[candidate_id, , drop = FALSE])
    K <- as.integer(params_all$K)
    model_params <- params_all[setdiff(names(params_all), c("K", "H"))]
    for (fold in sort(unique(folds))) {
      fold_data <- cv_data_objects_leakage_free(data, folds, fold)
      model_config <- modifyList(fixed_model_config, model_params)
      model_config$b_sigma <- stats::var(fold_data$train$y) / 2
      fit_started <- proc.time()[3]
      fit <- tryCatch(
        method_function(
          fold_data$train,
          fold_data$validation,
          list(coords = fold_data$validation$coords),
          list(n_basis = K),
          model_config,
          list(n_iter = cv_n_iter, burn_in = cv_burn_in),
          seed = as.integer(seed_start + 1000L * candidate_id + fold)
        ),
        error = identity
      )
      success <- !inherits(fit, "error")
      results[[row_id]] <- data.frame(
        candidate_id = candidate_id,
        fold = fold,
        validation_mse = if (success) {
          mean((fold_data$validation$y - fit$pred_test)^2)
        } else {
          Inf
        },
        success = success,
        error = if (success) "" else conditionMessage(fit),
        runtime_seconds = proc.time()[3] - fit_started,
        preprocessing_checksum = fold_assignment_checksum(folds),
        stringsAsFactors = FALSE
      )
      for (name in names(params_all)) results[[row_id]][[name]] <- params_all[[name]]
      row_id <- row_id + 1L
    }
  }
  fold_results <- do.call(rbind, results)
  candidate_results <- do.call(rbind, lapply(
    split(fold_results, fold_results$candidate_id),
    function(x) {
      parameter_columns <- setdiff(
        names(x),
        c("candidate_id", "fold", "validation_mse", "success",
          "error", "runtime_seconds", "preprocessing_checksum")
      )
      out <- data.frame(
        candidate_id = x$candidate_id[1],
        cv_mse = if (all(x$success)) mean(x$validation_mse) else Inf,
        cv_mse_sd = if (all(x$success)) sd(x$validation_mse) else NA_real_,
        successful_folds = sum(x$success),
        total_runtime_seconds = sum(x$runtime_seconds),
        preprocessing_checksum = x$preprocessing_checksum[1],
        stringsAsFactors = FALSE
      )
      for (name in parameter_columns) out[[name]] <- x[[name]][1]
      out
    }
  ))
  candidate_results <- candidate_results[
    order(candidate_results$candidate_id), ]
  list(
    folds = folds,
    fold_results = fold_results,
    candidate_results = candidate_results,
    best = select_best_cv_candidate(candidate_results),
    runtime_seconds = proc.time()[3] - started,
    fold_metadata = fold_metadata
  )
}

make_gam_reml_default_formula <- function(p) {
  smooth_terms <- paste0(
    "s(s1, s2, by = X", seq_len(p), ", bs = 'tp')",
    collapse = " + "
  )
  as.formula(paste0("y ~ 0 + ", smooth_terms))
}

fit_gam_reml_default_fold <- function(train, validation) {
  p <- ncol(train$X)
  train_df <- as.data.frame(train$X)
  names(train_df) <- paste0("X", seq_len(p))
  train_df$y <- train$y
  train_df$s1 <- train$coords[, 1]
  train_df$s2 <- train$coords[, 2]
  fit <- mgcv::gam(
    make_gam_reml_default_formula(p),
    data = train_df,
    method = "REML"
  )
  validation_df <- as.data.frame(validation$X)
  names(validation_df) <- paste0("X", seq_len(p))
  validation_df$s1 <- validation$coords[, 1]
  validation_df$s2 <- validation$coords[, 2]
  list(
    fit = fit,
    pred = as.numeric(stats::predict(fit, newdata = validation_df)),
    k_check = suppressWarnings(mgcv::k.check(fit))
  )
}

cv_gam_reml_default_leakage_free <- function(data, folds) {
  rows <- list()
  row_id <- 1L
  started <- proc.time()[3]
  for (fold in sort(unique(folds))) {
    fold_data <- cv_data_objects_leakage_free(data, folds, fold)
    fit_started <- proc.time()[3]
    result <- tryCatch(
      fit_gam_reml_default_fold(fold_data$train, fold_data$validation),
      error = identity
    )
    success <- !inherits(result, "error")
    rows[[row_id]] <- data.frame(
      candidate_id = 1L,
      fold = fold,
      validation_mse = if (success) {
        mean((fold_data$validation$y - result$pred)^2)
      } else {
        Inf
      },
      success = success,
      error = if (success) "" else conditionMessage(result),
      runtime_seconds = proc.time()[3] - fit_started,
      preprocessing_checksum = fold_assignment_checksum(folds),
      stringsAsFactors = FALSE
    )
    row_id <- row_id + 1L
  }
  fold_results <- do.call(rbind, rows)
  candidate_results <- data.frame(
    candidate_id = 1L,
    cv_mse = if (all(fold_results$success)) {
      mean(fold_results$validation_mse)
    } else {
      Inf
    },
    cv_mse_sd = if (all(fold_results$success)) {
      sd(fold_results$validation_mse)
    } else {
      NA_real_
    },
    successful_folds = sum(fold_results$success),
    total_runtime_seconds = sum(fold_results$runtime_seconds),
    preprocessing_checksum = fold_results$preprocessing_checksum[1],
    method = "REML_default_k",
    stringsAsFactors = FALSE
  )
  list(
    folds = folds,
    fold_results = fold_results,
    candidate_results = candidate_results,
    best = candidate_results,
    runtime_seconds = proc.time()[3] - started
  )
}
