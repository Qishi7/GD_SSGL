#!/usr/bin/env python3

import json
import os
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd
from mgwr.gwr import MGWR
from mgwr.sel_bw import Sel_BW


PROJECT_ROOT = Path(os.environ.get("GDSSGL_ROOT", Path.cwd())).resolve()
OUT_ROOT = PROJECT_ROOT / "new_method_comparison" / "real_data_analysis" / "modis_n10000_intercept_svd_centered_y_check"
DATA_DIR = OUT_ROOT / "data"
RESULTS_DIR = OUT_ROOT / "results"
REPORTS_DIR = OUT_ROOT / "reports"
LOGS_DIR = OUT_ROOT / "logs"
for d in (DATA_DIR, RESULTS_DIR, REPORTS_DIR, LOGS_DIR):
    d.mkdir(parents=True, exist_ok=True)


def metric_row(method, observed, predicted, runtime_sec, scale_name, beta0_hat):
    resid = observed - predicted
    return {
        "method": method,
        "scale": scale_name,
        "MSPE": float(np.mean(resid ** 2)),
        "RMSE": float(np.sqrt(np.mean(resid ** 2))),
        "MAE": float(np.mean(np.abs(resid))),
        "bias": float(np.mean(resid)),
        "test_R2": float(1.0 - np.sum(resid ** 2) / np.sum((observed - np.mean(observed)) ** 2)),
        "correlation": float(np.corrcoef(observed, predicted)[0, 1]),
        "beta0_hat": float(beta0_hat),
        "runtime_sec": float(runtime_sec),
    }


def idw_coefficients(train_coords, coef_mat, new_coords, power=2.0):
    train_coords = np.asarray(train_coords, dtype=float)
    new_coords = np.asarray(new_coords, dtype=float)
    coef_mat = np.asarray(coef_mat, dtype=float)
    out = np.empty((new_coords.shape[0], coef_mat.shape[1]), dtype=float)
    for i, point in enumerate(new_coords):
        d = np.sqrt(np.sum((train_coords - point[None, :]) ** 2, axis=1))
        exact = np.where(d < 1e-12)[0]
        if exact.size:
            out[i, :] = coef_mat[exact[0], :]
        else:
            w = 1.0 / np.maximum(d, 1e-12) ** power
            w /= np.sum(w)
            out[i, :] = w @ coef_mat
    return out


def main():
    input_tag = os.environ.get("PY_MGWR_INPUT_TAG", "").strip()
    input_mid = f"_{input_tag}" if input_tag else ""
    train = pd.read_csv(DATA_DIR / f"modis_centered_y_mgwr_python{input_mid}_train.csv")
    test = pd.read_csv(DATA_DIR / f"modis_centered_y_mgwr_python{input_mid}_test.csv")
    metadata = pd.read_csv(DATA_DIR / f"modis_centered_y_mgwr_python{input_mid}_metadata.csv")
    y_train_mean = float(metadata.loc[metadata["key"] == "y_train_mean", "value"].iloc[0])

    x_cols = [c for c in train.columns if c not in ("y_centered", "y_log", "EVI", "s1", "s2")]
    coords_train = train[["s1", "s2"]].to_numpy(float)
    coords_test = test[["s1", "s2"]].to_numpy(float)
    y = train[["y_centered"]].to_numpy(float)
    x_train = train[x_cols].to_numpy(float)
    x_test = test[x_cols].to_numpy(float)

    method = "MGWR (Python mgwr)" if not input_tag else f"MGWR (Python mgwr, {input_tag})"
    started = time.time()
    selector = Sel_BW(
        coords_train,
        y,
        x_train,
        multi=True,
        constant=True,
        fixed=False,
        kernel="bisquare",
        n_jobs=1,
    )
    # Use moderate iteration limits to keep the real-data run feasible while
    # retaining true multiscale bandwidth selection.
    min_bw = int(os.environ.get("PY_MGWR_MIN_BW", "1000"))
    max_iter_multi = int(os.environ.get("PY_MGWR_MAX_ITER_MULTI", "50"))
    bws = selector.search(
        criterion="AICc",
        max_iter=200,
        max_iter_multi=max_iter_multi,
        tol=1e-6,
        tol_multi=1e-5,
        multi_bw_min=[min_bw] * (x_train.shape[1] + 1),
        verbose=True,
    )
    model = MGWR(
        coords_train,
        y,
        x_train,
        selector,
        fixed=False,
        kernel="bisquare",
        constant=True,
        hat_matrix=False,
        n_jobs=1,
        name_x=x_cols,
    )
    results = model.fit()
    runtime_sec = time.time() - started

    coef_train = np.asarray(results.params, dtype=float)
    coef_test = idw_coefficients(coords_train, coef_train, coords_test)
    pred_log_centered = coef_test[:, 0] + np.sum(x_test * coef_test[:, 1:], axis=1)
    pred_log = pred_log_centered + y_train_mean
    pred_evi = np.exp(pred_log) - 1.0
    obs_log = test["y_log"].to_numpy(float)
    obs_evi = test["EVI"].to_numpy(float)
    beta0_hat = float(np.mean(coef_test[:, 0]) + y_train_mean)

    pred_df = pd.DataFrame({
        "method": method,
        "observed_log": obs_log,
        "predicted_log": pred_log,
        "residual_log": obs_log - pred_log,
        "observed_evi": obs_evi,
        "predicted_evi": pred_evi,
        "residual_evi": obs_evi - pred_evi,
        "s1": coords_test[:, 0],
        "s2": coords_test[:, 1],
    })
    pred_path = DATA_DIR / f"modis_centered_y_mgwr_python{input_mid}_predictions.csv"
    pred_df.to_csv(pred_path, index=False)

    metrics = pd.DataFrame([
        metric_row(method, obs_log, pred_log, runtime_sec, "log(EVI+1)", beta0_hat),
        metric_row(method, obs_evi, pred_evi, runtime_sec, "EVI", beta0_hat),
    ])
    metrics_path = DATA_DIR / f"modis_centered_y_mgwr_python{input_mid}_prediction_metrics_log_and_evi_scale.csv"
    metrics.to_csv(metrics_path, index=False)

    primary_metrics_path = DATA_DIR / f"modis_centered_y_mgwr{input_mid}_prediction_metrics_log_and_evi_scale.csv"
    metrics.to_csv(primary_metrics_path, index=False)

    fit_path = RESULTS_DIR / f"mgwr_intercept_centered_y_python_mgwr{input_mid}_fit.pkl"
    with fit_path.open("wb") as f:
        pickle.dump({
            "method": method,
            "bws": np.asarray(bws).tolist(),
            "params": coef_train,
            "x_cols": x_cols,
            "y_train_mean": y_train_mean,
            "runtime_sec": runtime_sec,
            "selector": selector,
            "summary": {
                "aicc": float(getattr(results, "aicc", np.nan)),
                "bic": float(getattr(results, "bic", np.nan)),
                "sigma2": float(getattr(results, "sigma2", np.nan)),
                "min_bw": min_bw,
                "max_iter_multi": max_iter_multi,
            },
        }, f)

    diag_path = DATA_DIR / f"modis_centered_y_mgwr_python{input_mid}_bandwidths.csv"
    pd.DataFrame({"term": ["Intercept"] + x_cols, "bandwidth": np.asarray(bws, dtype=float)}).to_csv(diag_path, index=False)

    all_evi_path = DATA_DIR / "modis_n10000_centered_y_check_prediction_metrics_evi_scale.csv"
    if all_evi_path.exists():
        all_evi = pd.read_csv(all_evi_path)
        mgwr_evi = metrics.loc[metrics["scale"] == "EVI"].iloc[0]
        mgwr_evi_out = pd.DataFrame([{
            "method": method,
            "MSPE_EVI": mgwr_evi["MSPE"],
            "RMSE_EVI": mgwr_evi["RMSE"],
            "MAE_EVI": mgwr_evi["MAE"],
            "bias_EVI": mgwr_evi["bias"],
            "test_R2_EVI": mgwr_evi["test_R2"],
            "correlation_EVI": mgwr_evi["correlation"],
            "min_pred_EVI": float(np.min(pred_evi)),
            "max_pred_EVI": float(np.max(pred_evi)),
        }])
        combined = pd.concat([all_evi, mgwr_evi_out[all_evi.columns]], ignore_index=True)
        combined = combined.sort_values("MSPE_EVI")
        combined.to_csv(DATA_DIR / f"modis_n10000_centered_y_check_prediction_metrics_evi_scale_with_mgwr_python{input_mid}.csv", index=False)

    report_path = REPORTS_DIR / f"modis_centered_y_mgwr_python{input_mid}_evi_metrics_report.md"
    report_path.write_text(
        "# MODIS centered-y MGWR Python EVI-scale metrics\n\n"
        f"- Output root: `{OUT_ROOT}`.\n"
        f"- Input tag: `{input_tag if input_tag else 'default_p15_dummy_lc'}`.\n"
        "- Implementation: Python `mgwr` package with multiscale bandwidth selection.\n"
        "- Fit response: centered `log(EVI+1)`.\n"
        "- Intercept: included via `constant=True`.\n"
        "- Prediction: IDW interpolation of fitted training-location local coefficients to test locations.\n"
        f"- Minimum adaptive bandwidth allowed during selection: {min_bw} nearest neighbors.\n"
        f"- Maximum multiscale iterations: {max_iter_multi}.\n"
        f"- Training response mean added back before EVI transform: {y_train_mean:.12f}.\n"
        f"- Runtime sec: {runtime_sec:.3f}.\n\n"
        "## Bandwidths\n\n"
        + pd.DataFrame({"term": ["Intercept"] + x_cols, "bandwidth": np.asarray(bws, dtype=float)}).to_markdown(index=False)
        + "\n\n## Metrics\n\n"
        + metrics.to_markdown(index=False)
        + "\n"
    )

    status = {
        "fit_path": str(fit_path),
        "metrics_path": str(metrics_path),
        "primary_metrics_path": str(primary_metrics_path),
        "predictions_path": str(pred_path),
        "bandwidths_path": str(diag_path),
        "report_path": str(report_path),
        "runtime_sec": runtime_sec,
    }
    (RESULTS_DIR / f"mgwr_intercept_centered_y_python_mgwr{input_mid}_status.json").write_text(json.dumps(status, indent=2))
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
