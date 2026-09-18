# Reproducibility Notes

This repository is a cleaned code release. It intentionally excludes:

- MCMC fit objects;
- raw simulation output folders;
- MODIS `.RData` files;
- generated figures and tables;
- intermediate exploratory diagnostics.

The main simulation protocol is the current intercept/SVD version:

- `beta0 = 0.5` in the data-generating process;
- explicit constant-intercept samplers for GD-SSGL, WS-SSGL, Gaussian SVC, and Bayesian Lasso;
- SVD full-rank centered basis for GD-SSGL and Gaussian SVC;
- uncentered whole-surface basis for WS-SSGL;
- per-replicate CV for GD-SSGL, WS-SSGL, and Gaussian SVC in the main formal simulation scripts.

For MODIS, provide the data path via:

```bash
export MODIS_RDATA=/path/to/data_cleaned_small_expanded.RData
```

For map plotting, provide Natural Earth data via:

```bash
export NATURALEARTH_DIR=/path/to/geo_naturalearth_10m
```

