# Concise public entry points.
#
# The validated implementations remain in wrappers.R so existing scripts and
# frozen results stay reproducible.

find_ssgl_root <- function() {
  candidates <- c(getwd(), file.path(getwd(), "new_method_comparison"))
  hit <- candidates[file.exists(file.path(
    candidates, "R", "methods", "wrappers.R"
  ))]
  if (!length(hit)) stop("Cannot locate new_method_comparison.")
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

.ssgl_root <- find_ssgl_root()
source(file.path(.ssgl_root, "R", "methods", "wrappers.R"))

fit_newssgl <- fit_proposed_ssgl
fit_oldssgl <- fit_original_ssgl

load_newssgl_intercept_entry <- function(rebuild = FALSE) {
  source(file.path(.ssgl_root, "R", "methods", "accelerated",
                   "newssgl_intercept.R"))
  load_newssgl_intercept_fast(rebuild = rebuild)
  invisible(TRUE)
}

load_matched_wsssgl_intercept_entry <- function(rebuild = FALSE) {
  source(file.path(.ssgl_root, "R", "methods", "accelerated",
                   "matched_wsssgl_intercept.R"))
  load_matched_wsssgl_intercept_fast(rebuild = rebuild)
  invisible(TRUE)
}
