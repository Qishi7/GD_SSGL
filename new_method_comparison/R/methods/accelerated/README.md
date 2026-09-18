# newssgl

Public names:

- `fit_newssgl()`: validated RcppArmadillo implementation, current default.
- `fit_newssgl_reference()`: retained original R reference implementation.
- `fit_newssgl_fast()`: validated RcppArmadillo implementation.
- `fit_oldssgl()`: old SSGL.

The original R reference is retained for audit/reproducibility. After the
validation in `results_debug/newssgl_validation/`, future runs should use the
RcppArmadillo implementation through `fit_newssgl()` or `fit_newssgl_fast()`.

- `newssgl.cpp`: RcppArmadillo Gibbs loop and deterministic
  conditional-parameter helper.
- `newssgl.R`: loader, shared input preparation,
  accelerated wrapper, and an R diagnostic reference sampler.

Validation is run by:

```r
source("new_method_comparison/tests/test_newssgl.R")
```

Results are written to:

`new_method_comparison/results_debug/newssgl_validation/`
