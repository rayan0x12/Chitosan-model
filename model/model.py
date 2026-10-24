"""One-command GMDH pipeline: train / validate / sensitivity / dose-sweep.

Run from anywhere with a single command::

    python model/model.py

The pipeline is a *framework*.  It reads the curated dataset
(``data/sheet1_cleaned.csv``), builds the three nested input sets (A/B/C) and
the two normalised targets (``ratio`` = UCS_treated/q0 and its log), then runs:

  1. Training with four validation protocols:
       - random 80/20 split
       - k-fold cross-validation (k = 5)
       - leave-one-study-out (LOSO) — the headline protocol
       - bootstrap resampling
  2. Sensitivity / importance:
       - Liong et al. (2000) one-at-a-time perturbation
       - SHAP values (if the ``shap`` package is installed)
       - gradient (finite-difference) importance of the white-box polynomial
  3. A dose-sweep: predict UCS/q0 across a chitosan-dose grid at fixed values
     of the other inputs.

All random steps use the fixed seeds in ``config.py``.  Every artefact is
written to ``output/`` (gitignored) so the repository never ships with
calibration numbers baked into the source tree.

Status: model calibration and validation — in progress.
"""

from __future__ import annotations

import csv
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Make this file importable as a module and as a script from the repo root.
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from config import (  # noqa: E402
    DATA_DIR, SHEET1, FIG_DIR, OUT_DIR, FEATURE_LABELS, INPUT_SETS, TARGETS,
    TARGET_LABELS, FEATURE_KEY, GMDH_KWARGS, SEEDS, KFOLD, N_BOOTSTRAP,
    CHANGE_FACTORS, PLOT_RC, PALETTE, DPI,
)
from gmdh import GMDH  # noqa: E402

plt.rcParams.update(PLOT_RC)
os.makedirs(OUT_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# Numeric parsing and data loading (read-only on the source CSV)
# ---------------------------------------------------------------------------
_NUM = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?")


def parse_num(cell):
    """Parse a cell that may be ``"13.16 (derived ...)"`` -> 13.16, else None."""
    if cell is None:
        return None
    s = str(cell).strip()
    if s == "":
        return None
    try:
        return float(s)
    except ValueError:
        pass
    m = _NUM.match(s)
    if not m:
        return None
    rest = s[m.end():].strip()
    if rest == "" or rest.startswith("(") or rest.startswith("["):
        return float(m.group())
    return None


def load_rows():
    """Return one dict per eligible row, with the v2 inclusion rules applied."""
    with open(SHEET1, newline="", encoding="utf-8-sig") as fh:
        rd = list(csv.reader(fh))
    hdr = rd[0]
    idx = {h: i for i, h in enumerate(hdr)}
    rows = []
    for r in rd[1:]:
        ref = r[idx["Ref_ID"]].strip()
        status = r[idx["Row_Status"]].strip()
        if ref in ("R08", "R11"):      # R08 = secondary values; R11 dropped
            continue
        if status.lower() == "incomplete":
            continue
        cyc = parse_num(r[idx["WetDry_FT_Cycles_Number"]])
        if cyc is not None and cyc > 0:  # no wet-dry / freeze-thaw cycled rows
            continue
        rows.append({
            "ref": ref,
            "pi": parse_num(r[idx["PI_untreated"]]),
            "cc": parse_num(r[idx["Chitosan_Content_%"]]),
            "cti": parse_num(r[idx["Curing_Time_days"]]),
            "wmc": parse_num(r[idx["Initial_Moisture_Content_%"]]),
            "q0": parse_num(r[idx["UCS_untreated_kPa_q0"]]),   # kPa
            "ucs": parse_num(r[idx["UCS_treated_kPa"]]),       # kPa
            "ctype": r[idx["Chitosan_Type"]].strip(),
            "solv": r[idx["Solvent"]].strip(),
            "dda": r[idx["DDA_%"]].strip(),
        })
    return rows


def is_calsil(r: dict) -> bool:
    t = r["ctype"].lower()
    return ("calcium metasilicate" in t) or ("calsil" in t)


def study_label(r: dict) -> str:
    if r["ref"] == "R13":
        return "R13 (CTS + calsil)" if is_calsil(r) else "R13 (CTS)"
    return r["ref"]


def build_df(rows, feats):
    """Keep rows with every feature + q0 + ucs present and positive."""
    out = []
    for r in rows:
        vals = [r["q0"], r["ucs"]] + [r[FEATURE_KEY[f]] for f in feats]
        if any(v is None for v in vals):
            continue
        if r["q0"] <= 0 or r["ucs"] <= 0:
            continue
        out.append(r)
    df = pd.DataFrame(out)
    df["ratio"] = df["ucs"] / df["q0"]
    df["logratio"] = np.log(df["ucs"] / df["q0"])
    df["treat"] = [study_label(r) for r in out]
    return df


def Xy(df, feats, target):
    cols = []
    for f in feats:
        v = df[FEATURE_KEY[f]].astype(float).values
        if f == "q0":
            v = v / 1000.0  # kPa -> MPa (feature convention)
        cols.append(v)
    X = np.column_stack(cols)
    y = df[target].astype(float).values
    return X, y


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def metrics(y, yp):
    y = np.asarray(y, float)
    yp = np.asarray(yp, float)
    ss_res = float(np.sum((y - yp) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else np.nan
    rmse = float(np.sqrt(np.mean((y - yp) ** 2)))
    mae = float(np.mean(np.abs(y - yp)))
    return r2, rmse, mae


def _fit(Xtr, ytr):
    m = GMDH(random_state=SEEDS["gmdh"], **GMDH_KWARGS)
    m.fit(Xtr, ytr)
    return m


# ---------------------------------------------------------------------------
# Validation protocols
# ---------------------------------------------------------------------------
def split_80_20(df, feats, target):
    """One deterministic 80/20 train/test split."""
    X, y = Xy(df, feats, target)
    rng = np.random.RandomState(SEEDS["split"])
    perm = rng.permutation(len(y))
    n_tr = int(round(0.8 * len(y)))
    tr, te = perm[:n_tr], perm[n_tr:]
    m = _fit(X[tr], y[tr])
    return metrics(y[tr], m.predict(X[tr])), metrics(y[te], m.predict(X[te]))


def kfold_cv(df, feats, target, k=KFOLD):
    """Deterministic k-fold cross-validation; returns fold test metrics."""
    X, y = Xy(df, feats, target)
    rng = np.random.RandomState(SEEDS["cv"])
    perm = rng.permutation(len(y))
    folds = np.array_split(perm, k)
    out = []
    for f in folds:
        mask = np.ones(len(y), dtype=bool)
        mask[f] = False
        m = _fit(X[mask], y[mask])
        out.append(metrics(y[f], m.predict(X[f])))
    return out


def leave_one_study_out(df, feats, target):
    """Leave-one-study-out: train on all studies but one, predict the held-out."""
    y_all, yp_all, per = [], [], {}
    for study, sub in df.groupby("ref"):
        tr = df[df["ref"] != study]
        if len(tr) < 8:
            continue
        Xtr, ytr = Xy(tr, feats, target)
        Xte, yte = Xy(sub, feats, target)
        m = _fit(Xtr, ytr)
        yp = m.predict(Xte)
        y_all += list(yte)
        yp_all += list(yp)
        per[study] = metrics(yte, yp)
    return metrics(np.array(y_all), np.array(yp_all)), per


def bootstrap(df, feats, target, n_boot=N_BOOTSTRAP):
    """Bootstrap resample the full dataset; returns the distribution of R2/RMSE."""
    X, y = Xy(df, feats, target)
    rng = np.random.RandomState(SEEDS["bootstrap"])
    n = len(y)
    r2s, rmses = [], []
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        m = _fit(X[idx], y[idx])
        r2, rmse, _ = metrics(y[idx], m.predict(X[idx]))
        r2s.append(r2)
        rmses.append(rmse)
    return np.array(r2s), np.array(rmses)


# ---------------------------------------------------------------------------
# Sensitivity / importance
# ---------------------------------------------------------------------------
def liong_sensitivity(model, X, feats):
    """Liong et al. (2000): percent change in output under per-input perturbation."""
    base = model.predict(X)
    out = {}
    for k, f in enumerate(feats):
        deltas = []
        for cf in CHANGE_FACTORS:
            Xp = X.copy()
            Xp[:, k] = Xp[:, k] * cf
            deltas.append(float(np.mean(np.abs((model.predict(Xp) - base) / (base + 1e-12)))))
        out[f] = float(np.mean(deltas))
    return out


def gradient_importance(model, X, feats):
    """Mean |d(output)/d(input)| via central finite differences on raw inputs."""
    out = {}
    for k, f in enumerate(feats):
        s = X[:, k].std() * 1e-3 + 1e-9
        Xp, Xm = X.copy(), X.copy()
        Xp[:, k] += s
        Xm[:, k] -= s
        grad = (model.predict(Xp) - model.predict(Xm)) / (2 * s)
        out[f] = float(np.mean(np.abs(grad)))
    return out


def shap_importance(model, X, feats):
    """SHAP values for the GMDH (callable), if the ``shap`` package is present."""
    try:
        import shap  # type: ignore
    except Exception:
        return None
    # KernelSHAP needs a predict that returns a scalar per row.
    def f(x):
        return model.predict(x)
    X_bg = X[: min(len(X), 50)]
    explainer = shap.KernelExplainer(f, X_bg)
    sh = explainer.shap_values(X[: min(len(X), 200)], nsamples=50)
    if isinstance(sh, list):
        sh = sh[0]
    return {f: float(np.mean(np.abs(sh[:, k]))) for k, f in enumerate(feats)}


# ---------------------------------------------------------------------------
# Dose sweep (framework: predict UCS/q0 over a dose grid)
# ---------------------------------------------------------------------------
def dose_sweep(df, feats, target, dose_grid=None):
    """Sweep the chitosan dose at median values of the other inputs.

    Returns a DataFrame with the swept dose and the predicted target.  This is
    a *procedure*, not a reported optimum: any optimum that sits on the grid
    boundary is an artefact of the sweep range, not a physical optimum.
    """
    if dose_grid is None:
        dose_grid = np.linspace(0.1, 2.5, 49)
    X, y = Xy(df, feats, target)
    m = _fit(X, y)

    med = np.median(X, axis=0)
    # locate the dose (CC) column
    cc_idx = feats.index("CC")
    rows = []
    for cc in dose_grid:
        x = med.copy()
        x[cc_idx] = cc
        rows.append([cc, float(m.predict(x.reshape(1, -1))[0])])
    return pd.DataFrame(rows, columns=["Chitosan_dose_%", "predicted_" + target])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    rows = load_rows()
    dfs = {name: build_df(rows, feats) for name, feats in INPUT_SETS.items()}

    # 1) validation table
    table = []
    for name, feats in INPUT_SETS.items():
        df = dfs[name]
        for tgt in TARGETS:
            tr_met, te_met = split_80_20(df, feats, tgt)
            folds = kfold_cv(df, feats, tgt)
            loso_met, _ = leave_one_study_out(df, feats, tgt)
            r2_boot, rmse_boot = bootstrap(df, feats, tgt)
            table.append({
                "set": name, "target": tgt, "n": len(df),
                "studies": df["ref"].nunique(),
                "train_r2": tr_met[0], "test_r2": te_met[0], "test_rmse": te_met[1],
                "cv_r2_mean": float(np.mean([m[0] for m in folds])),
                "cv_rmse_mean": float(np.mean([m[1] for m in folds])),
                "loso_r2": loso_met[0], "loso_rmse": loso_met[1],
                "boot_r2_med": float(np.median(r2_boot)),
                "boot_rmse_med": float(np.median(rmse_boot)),
            })
    pd.DataFrame(table).to_csv(os.path.join(OUT_DIR, "validation.csv"), index=False)

    # 2) sensitivity / importance (on set A, ratio target — the largest set)
    dfA = dfs["A"]
    X, y = Xy(dfA, INPUT_SETS["A"], "ratio")
    m = _fit(X, y)
    sens = {
        "liong": liong_sensitivity(m, X, INPUT_SETS["A"]),
        "gradient": gradient_importance(m, X, INPUT_SETS["A"]),
        "shap": shap_importance(m, X, INPUT_SETS["A"]),
    }
    with open(os.path.join(OUT_DIR, "sensitivity.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["method"] + INPUT_SETS["A"])
        for method, vals in sens.items():
            if vals is None:
                continue
            w.writerow([method] + [f"{vals[f]:.4g}" for f in INPUT_SETS["A"]])

    # 3) dose sweep (procedure only)
    for name, feats in INPUT_SETS.items():
        df = dfs[name]
        for tgt in TARGETS:
            sweep = dose_sweep(df, feats, tgt)
            sweep.to_csv(os.path.join(OUT_DIR, f"dose_sweep_{name}_{tgt}.csv"), index=False)

    # 4) a small model-figure for the output folder (parity on set A, ratio)
    _parity_fig(dfA, INPUT_SETS["A"], "ratio", os.path.join(OUT_DIR, "parity_A_ratio.png"))

    print(f"rows loaded: {len(rows)}")
    for name, df in dfs.items():
        print(f"  set {name}: n={len(df)} studies={df['ref'].nunique()}")
    print("wrote validation.csv, sensitivity.csv, dose_sweep_*.csv -> output/")


def _parity_fig(df, feats, target, fname):
    y_all, yp_all = [], []
    for study, sub in df.groupby("ref"):
        tr = df[df["ref"] != study]
        if len(tr) < 8:
            continue
        Xtr, ytr = Xy(tr, feats, target)
        Xte, yte = Xy(sub, feats, target)
        m = _fit(Xtr, ytr)
        yp = m.predict(Xte)
        y_all += list(yte)
        yp_all += list(yp)
    y = np.array(y_all)
    yp = np.array(yp_all)
    fig, ax = plt.subplots(figsize=(4, 4))
    lo, hi = min(y.min(), yp.min()), max(y.max(), yp.max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=1, alpha=0.5)
    ax.scatter(y, yp, s=18, alpha=0.75, color="tab:blue")
    ax.set_xlabel(TARGET_LABELS[target] + " (measured)")
    ax.set_ylabel(TARGET_LABELS[target] + " (predicted, LOSO)")
    ax.set_title("Set A, ratio target — LOSO parity")
    fig.tight_layout()
    fig.savefig(fname, dpi=DPI)
    plt.close(fig)


if __name__ == "__main__":
    main()
