# GD-SSGL

Clean reproducible code for the **global + spatial-deviation spike-and-slab group lasso (GD-SSGL)** model.

The current baseline protocol uses:

- response model with an explicit constant intercept;
- true simulation intercept `beta0 = 0.5`;
- SVD full-rank centered tensor-product B-spline basis for GD-SSGL and Gaussian SVC;
- original uncentered whole-surface basis for WS-SSGL;
- per-replicate cross-validation for GD-SSGL, WS-SSGL, and Gaussian SVC in the main formal simulations;
- grouped reporting for land-cover variables in the MODIS application.

Large simulation outputs, fitted objects, MODIS data, figures, and temporary diagnostics are intentionally excluded from this repository.

## Repository Layout

```text
new_method_comparison/
  R/
    common/                         shared utilities
    methods/
      accelerated/                  Rcpp/RcppArmadillo samplers
      wrappers.R                    method wrappers
      ssgl.R                        reference/helper code
  formal_simulation/
    R/                              simulation framework, tuning, reporting helpers
    configs/                        tuning grids and protocol settings
    main_simulation/                reproducible experiment runners
    scripts/                        legacy/general simulation scripts
    tests/                          formal simulation tests
  real_data_analysis/modis/scripts/ MODIS analysis runners and plotting scripts
  tests/                            sampler and implementation validation tests
```

## Requirements

Core R packages:

```r
install.packages(c(
  "Rcpp", "RcppArmadillo", "GIGrvg", "mgcv", "MASS",
  "Matrix", "splines", "sp", "sf", "ggplot2"
))
```

Optional real-data comparators and figure scripts may also use:

```r
install.packages(c("GWmodel", "spData"))
```

Python MGWR real-data follow-up uses:

```bash
python -m pip install mgwr numpy pandas
```

## Quick Smoke Tests

From the repository root:

```bash
Rscript new_method_comparison/tests/test_newssgl_intercept_minimal.R
Rscript new_method_comparison/tests/compare_newssgl_intercept_samplers.R
```

The scripts compile the Rcpp samplers on demand using `Rcpp::sourceCpp()`.

To run a complete one-replicate beta0 = 0.5 intercept/SVD example with
GD-SSGL, WS-SSGL, Gaussian SVC, Bayesian Lasso, and GAM:

```bash
Rscript new_method_comparison/formal_simulation/main_simulation/four_function_beta0_0p5_one_rep_direct_intercept_svd/scripts/run_beta0_0p5_one_rep.R
```

This demo writes metrics and PIP summaries under
`new_method_comparison/formal_simulation/main_simulation/four_function_beta0_0p5_one_rep_direct_intercept_svd/`.
It is a runnable verification example, not the full per-replicate-tuned
100-replicate analysis.

## Main Formal Simulation Runners

Run commands from the repository root. The scripts write outputs under their corresponding `new_method_comparison/formal_simulation/main_simulation/...` folders.

### Main four-function simulation

```bash
Rscript new_method_comparison/formal_simulation/main_simulation/four_function_beta0_0p5_100rep_perrep_tuned_intercept_svd/scripts/run_beta0_0p5_100rep_perrep_tuned.R --workers=1
```

### No-spatial-deviation diagnostic

```bash
Rscript new_method_comparison/formal_simulation/main_simulation/no_spatial_deviation_beta0_0p5_100rep_perrep_tuned_intercept_svd/scripts/run_no_spatial_beta0_0p5_100rep_perrep_tuned.R --workers=1
```

### Correlated-predictor stress test, rho = 0.7

```bash
Rscript new_method_comparison/formal_simulation/main_simulation/correlated_predictor_rho07_beta0_0p5_50rep_perrep_tuned_intercept_svd/scripts/run_rho07_beta0_0p5_50rep_perrep_tuned.R --workers=1
```

### Weak spatial-signal gamma grid

```bash
Rscript new_method_comparison/formal_simulation/main_simulation/weak_spatial_gamma_0p15_0p25_beta0_0p5_intercept_svd_GDSSGL/perrep_tuned_R100/scripts/run_gamma_0p15_0p25_perrep_tuned_R100.R --workers=1
```

### W6 on-basis calibration

```bash
Rscript new_method_comparison/formal_simulation/main_simulation/w6_onbasis_K6_beta0_0p5_intercept_svd_one_rep_pilot/scripts/run_w6_beta0_0p5_intercept_svd_one_rep.R
Rscript new_method_comparison/formal_simulation/main_simulation/w6_onbasis_K6_beta0_0p5_intercept_svd_R100/scripts/run_w6_beta0_0p5_intercept_svd_R100.R
```

### MGWR simulation extensions

MGWR is fit as a post-hoc extension on frozen replicate datasets produced by
the main simulation runners. The R implementation uses `GWmodel::gwr.multiscale()`
with the default formula intercept.

```bash
Rscript new_method_comparison/formal_simulation/main_simulation/four_function_beta0_0p5_100rep_perrep_tuned_intercept_svd/scripts/run_mgwr_extension_100rep.R --workers=1
Rscript new_method_comparison/formal_simulation/main_simulation/no_spatial_deviation_beta0_0p5_100rep_perrep_tuned_intercept_svd/scripts/run_mgwr_extension_100rep.R --workers=1
Rscript new_method_comparison/formal_simulation/main_simulation/correlated_predictor_rho07_beta0_0p5_50rep_perrep_tuned_intercept_svd/scripts/run_mgwr_extension_R50.R --workers=1
```

## MODIS Real-data Application

The MODIS data are not included. Set `MODIS_RDATA` to the cleaned `.RData` file before running:

```bash
export MODIS_RDATA=/path/to/data_cleaned_small_expanded.RData
export GDSSGL_ROOT=$(pwd)
Rscript new_method_comparison/real_data_analysis/modis/scripts/run_modis_intercept_svd_centered_y_cv_tuned_correct.R
```

MGWR real-data follow-ups are separate because the full multiscale search can be
slow. The R/GWmodel runner supports `MODIS_MGWR_MODE=auto`,
`practical_fixedbw`, or `practical_localwls`:

```bash
export MODIS_RDATA=/path/to/data_cleaned_small_expanded.RData
export GDSSGL_ROOT=$(pwd)
MODIS_MGWR_MODE=practical_localwls Rscript new_method_comparison/real_data_analysis/modis/scripts/run_modis_centered_y_mgwr_only_evi_metrics.R
```

The Python `mgwr` grouped-LC runner expects CSV inputs written by the MODIS
analysis/preparation pipeline and uses `constant=True` for an explicit
intercept:

```bash
export MODIS_RDATA=/path/to/data_cleaned_small_expanded.RData
export GDSSGL_ROOT=$(pwd)
PY_MGWR_INPUT_TAG=grouped_lc MODIS_PY_MGWR_LC_ENCODING=grouped_numeric \
  Rscript new_method_comparison/real_data_analysis/modis/scripts/prepare_modis_centered_y_python_mgwr_inputs.R
PY_MGWR_INPUT_TAG=grouped_lc PY_MGWR_MIN_BW=4000 PY_MGWR_MAX_ITER_MULTI=1 \
  python new_method_comparison/real_data_analysis/modis/scripts/run_modis_centered_y_python_mgwr_evi_metrics.py
```

Optional map/figure scripts may use Natural Earth shapefiles. Set:

```bash
export NATURALEARTH_DIR=/path/to/geo_naturalearth_10m
```

## Path Handling

Most scripts assume they are launched from the repository root. Alternatively, set:

```bash
export GDSSGL_ROOT=/path/to/GD_SSGL
```

## Notes on Inclusion Probabilities

- GD-SSGL PIPs target **nonconstant spatial-deviation inclusion**.
- WS-SSGL PIPs target **whole coefficient-surface inclusion**.
- These two PIP definitions should not be interpreted as the same scientific target.
