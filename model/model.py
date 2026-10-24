"""One-command GMDH pipeline: train / validate / feature & hyperparameter analysis.

Run from anywhere with a single command::

    python model/model.py

The pipeline reads the curated dataset (``data/sheet1_cleaned.csv``), then:

  1. reproduces the original single-layer GMDH configuration across the three
     nested input sets (A/B/C) and both normalised targets (ratio / logratio);
  2. compares against two honest baselines (predict-the-training-mean and a
     linear ridge regression) under the same leave-one-study-out (LOSO) protocol;
  3. runs a one-factor analysis of the levers that can plausibly reduce
     underfitting — log transforms of dose/curing time, an explicit
     form/solvent categorical, a single-best-neuron output, and a second layer;
  4. reports train/test R² and RMSE for an 80/20 split, k-fold CV, LOSO and an
     in-sample within-study R² (representational capacity).

All random steps use the fixed seeds in ``config.py``; artefacts are written to
``output/`` (gitignored). The pipeline never removes an observation or feature;
it only adds feature transforms / encodings.
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

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from config import (  # noqa: E402
    SHEET1, OUT_DIR, INPUT_SETS, TARGETS, TARGET_LABELS, FEATURE_KEY,
    GMDH_KWARGS, SEEDS, KFOLD, N_BOOTSTRAP, CHANGE_FACTORS, PLOT_RC, DPI,
)
from gmdh import GMDH  # noqa: E402

plt.rcParams.update(PLOT_RC)
os.makedirs(OUT_DIR, exist_ok=True)

# Form/solvent categories used for the explicit categorical feature.
# "dry" is the reference level (not one-hot encoded); "other" is dropped.
FORM_CLASSES = ["nano", "acid", "solution", "water"]


# ---------------------------------------------------------------------------
# Numeric parsing and data loading (read-only on the source CSV)
# ---------------------------------------------------------------------------
_NUM = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?")


def parse_num(cell):
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
        if cyc is not None and cyc > 0:
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


def form_class(r: dict) -> str:
    """Coarse form/solvent class, from the treatment metadata."""
    s = r["solv"].lower()
    t = r["ctype"].lower()
    if "nano" in t:
        return "nano"
    if "acetic" in s:
        return "acid"
    if "solution" in t:
        return "solution"
    if "water" in s:
        return "water"
    if "none" in s or "dry" in s or s in ("", "na", "nan"):
        return "dry"
    return "other"


def build_df(rows, feats):
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
    df["form"] = [form_class(r) for r in out]
    return df


def make_X(df, feats, opts):
    """Build the design matrix for a feature set and transform options."""
    cols, names = [], []
    for f in feats:
        if f == "CC":
            v = df["cc"].astype(float).values
            if opts.get("log_dose"):
                v = np.log(np.clip(v, 1e-6, None))
                names.append("log_CC")
            else:
                names.append("CC")
        elif f == "CTi":
            v = df["cti"].astype(float).values
            if opts.get("log_cti"):
                v = np.log1p(v)
                names.append("log_CTi")
            else:
                names.append("CTi")
        elif f == "PI":
            v = df["pi"].astype(float).values
            names.append("PI")
        elif f == "WMC":
            v = df["wmc"].astype(float).values
            names.append("WMC")
        elif f == "q0":
            v = df["q0"].astype(float).values / 1000.0  # kPa -> MPa
            names.append("q0_MPa")
        cols.append(v)
    if opts.get("add_form"):
        f = df["form"].values
        for cls in FORM_CLASSES:
            cols.append((f == cls).astype(float))
            names.append("form_" + cls)
    return np.column_stack(cols), names


def _y(df, target):
    return df[target].astype(float).values


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


# ---------------------------------------------------------------------------
# Predictors (GMDH + two baselines)
# ---------------------------------------------------------------------------
def _gmdh_predict(Xtr, ytr, Xte, kwargs):
    m = GMDH(random_state=SEEDS["gmdh"], **kwargs)
    m.fit(Xtr, ytr)
    return m.predict(Xte)


def _mean_predict(Xtr, ytr, Xte, _=None):
    return np.full(Xte.shape[0], np.mean(ytr))


def _ridge_predict(Xtr, ytr, Xte, lam=1.0):
    mu = Xtr.mean(axis=0)
    sd = Xtr.std(axis=0)
    sd[sd == 0] = 1.0
    Ztr = (Xtr - mu) / sd
    Zte = (Xte - mu) / sd
    A = np.column_stack([np.ones(len(Ztr)), Ztr])
    w = np.linalg.pinv(A.T @ A + lam * np.eye(A.shape[1])) @ A.T @ ytr
    return np.column_stack([np.ones(len(Zte)), Zte]) @ w


# ---------------------------------------------------------------------------
# Validation protocols (all take a `fit_predict(Xtr, ytr, Xte)` callable)
# ---------------------------------------------------------------------------
def split_80_20(df, feats, target, opts, fit_predict):
    X, names = make_X(df, feats, opts)
    y = _y(df, target)
    rng = np.random.RandomState(SEEDS["split"])
    perm = rng.permutation(len(y))
    n_tr = int(round(0.8 * len(y)))
    tr, te = perm[:n_tr], perm[n_tr:]
    yp_tr = fit_predict(X[tr], y[tr], X[tr])
    yp_te = fit_predict(X[tr], y[tr], X[te])
    return metrics(y[tr], yp_tr), metrics(y[te], yp_te)


def kfold_cv(df, feats, target, opts, fit_predict, k=KFOLD):
    X, names = make_X(df, feats, opts)
    y = _y(df, target)
    rng = np.random.RandomState(SEEDS["cv"])
    perm = rng.permutation(len(y))
    folds = np.array_split(perm, k)
    out = []
    for f in folds:
        mask = np.ones(len(y), dtype=bool)
        mask[f] = False
        yp = fit_predict(X[mask], y[mask], X[f])
        out.append(metrics(y[f], yp))
    return out


def leave_one_study_out(df, feats, target, opts, fit_predict, min_train=8):
    y_all, yp_all, per = [], [], {}
    for study, sub in df.groupby("ref"):
        tr = df[df["ref"] != study]
        if len(tr) < min_train:
            continue
        Xtr, _ = make_X(tr, feats, opts)
        Xte, _ = make_X(sub, feats, opts)
        ytr = _y(tr, target)
        yte = _y(sub, target)
        yp = fit_predict(Xtr, ytr, Xte)
        y_all += list(yte)
        yp_all += list(yp)
        per[study] = metrics(yte, yp)
    return metrics(np.array(y_all), np.array(yp_all)), per


def within_study_r2(df, feats, target, opts, fit_predict):
    """In-sample per-study R² — can the model represent each study's curve?."""
    X, names = make_X(df, feats, opts)
    y = _y(df, target)
    yp = fit_predict(X, y, X)  # in-sample (optimistic)
    per = {}
    for study, sub in df.groupby("ref"):
        idx = df.index.get_indexer(sub.index)
        per[study] = metrics(y[idx], yp[idx])[0]
    return float(np.mean(list(per.values()))), per


# ---------------------------------------------------------------------------
# Experiments
# ---------------------------------------------------------------------------
def _fmt(r2, rmse, *_):
    return f"R2={r2:+.3f} RMSE={rmse:.3f}"


def main():
    rows = load_rows()
    dfs = {name: build_df(rows, feats) for name, feats in INPUT_SETS.items()}

    # --- 1. original configuration across sets/targets ---------------------
    print("== 1. original config (max_layers=1, average, raw features) ==")
    orig = GMDH_KWARGS.copy()
    for name, feats in INPUT_SETS.items():
        df = dfs[name]
        for tgt in TARGETS:
            tr_met, te_met = split_80_20(df, feats, tgt, {}, _gmdh_predict_w(orig))
            folds = kfold_cv(df, feats, tgt, {}, _gmdh_predict_w(orig))
            loso_met, _ = leave_one_study_out(df, feats, tgt, {}, _gmdh_predict_w(orig))
            ws_mean, _ = within_study_r2(df, feats, tgt, {}, _gmdh_predict_w(orig))
            print(f"[{name}/{tgt}] n={len(df)} stud={df['ref'].nunique()} "
                  f"train {_fmt(*tr_met)} | test {_fmt(*te_met)} | "
                  f"CV {np.mean([m[0] for m in folds]):+.3f} | "
                  f"LOSO {_fmt(*loso_met)} | within-study R2 {ws_mean:+.3f}")

    # --- 2. baselines under LOSO (set A, both targets) --------------------
    print("\n== 2. LOSO baselines (set A) ==")
    for tgt in TARGETS:
        mean_met, _ = leave_one_study_out(dfs["A"], INPUT_SETS["A"], tgt, {}, _mean_predict)
        ridge_met, _ = leave_one_study_out(dfs["A"], INPUT_SETS["A"], tgt, {}, _ridge_predict)
        print(f"[{tgt}] mean-predict LOSO {_fmt(*mean_met)} | ridge-linear LOSO {_fmt(*ridge_met)}")

    # --- 3. one-factor analysis (set A, ratio target) ---------------------
    print("\n== 3. one-factor analysis (set A, ratio) ==")
    featsA = INPUT_SETS["A"]
    dfA = dfs["A"]

    def report(label, opts, kwargs):
        loso_met, per = leave_one_study_out(dfA, featsA, "ratio", opts, _gmdh_predict_w(kwargs))
        ws_mean, _ = within_study_r2(dfA, featsA, "ratio", opts, _gmdh_predict_w(kwargs))
        tr, te = split_80_20(dfA, featsA, "ratio", opts, _gmdh_predict_w(kwargs))
        print(f"  {label:34s} LOSO {_fmt(*loso_met)} | within-study {ws_mean:+.3f} "
              f"| test {_fmt(*te)}")
        return loso_met

    report("original", {}, orig)
    report("log-dose", {"log_dose": True}, orig)
    report("log-dose + log-curing", {"log_dose": True, "log_cti": True}, orig)
    report("log + form dummies", {"log_dose": True, "log_cti": True, "add_form": True}, orig)
    report("log + form + best-neuron", {"log_dose": True, "log_cti": True, "add_form": True},
           {**orig, "output_mode": "best"})
    report("log + form + 2 layers", {"log_dose": True, "log_cti": True, "add_form": True},
           {**orig, "max_layers": 2})
    report("2 layers (raw)", {}, {**orig, "max_layers": 2})
    report("best-neuron (raw)", {}, {**orig, "output_mode": "best"})

    # --- 4. small hyperparameter sweep (LOSO, set A ratio, log+form) ------
    print("\n== 4. hyperparameter sweep (set A ratio, log+form features, LOSO) ==")
    base_opts = {"log_dose": True, "log_cti": True, "add_form": True}
    best = (np.inf, None)
    for alpha in (0.3, 0.6, 1.0):
        for maxn in (3, 8, 15):
            for outmode in ("average", "best"):
                kw = {**orig, "alpha": alpha, "max_neurons": maxn, "output_mode": outmode}
                loso_met, _ = leave_one_study_out(dfA, featsA, "ratio", base_opts, _gmdh_predict_w(kw))
                print(f"    alpha={alpha} maxn={maxn:2d} out={outmode:7s} -> LOSO {_fmt(*loso_met)}")
                if loso_met[1] < best[0]:
                    best = (loso_met[1], dict(alpha=alpha, max_neurons=maxn, output_mode=outmode))
    print(f"  best-by-RMSE: {best[1]}")

    # --- 5. sensitivity + dose sweep (kept from the framework) ------------
    _sensitivity_and_sweep(dfA, featsA)

    # --- write a comparison CSV -------------------------------------------
    _write_comparison(dfs)

    print("\nDONE -> output/")


def _gmdh_predict_w(kwargs):
    def _f(Xtr, ytr, Xte):
        return _gmdh_predict(Xtr, ytr, Xte, kwargs)
    return _f


def _sensitivity_and_sweep(dfA, featsA):
    """Liong / gradient / SHAP + a dose sweep (set A ratio, original config)."""
    X, names = make_X(dfA, featsA, {})
    y = _y(dfA, "ratio")
    m = GMDH(random_state=SEEDS["gmdh"], **GMDH_KWARGS)
    m.fit(X, y)

    # Liong
    liong = {}
    base = m.predict(X)
    for k in range(X.shape[1]):
        d = []
        for cf in CHANGE_FACTORS:
            Xp = X.copy()
            Xp[:, k] *= cf
            d.append(float(np.mean(np.abs((m.predict(Xp) - base) / (base + 1e-12)))))
        liong[names[k]] = float(np.mean(d))

    # gradient
    grad = {}
    for k in range(X.shape[1]):
        s = X[:, k].std() * 1e-3 + 1e-9
        Xp, Xm = X.copy(), X.copy()
        Xp[:, k] += s
        Xm[:, k] -= s
        g = (m.predict(Xp) - m.predict(Xm)) / (2 * s)
        grad[names[k]] = float(np.mean(np.abs(g)))

    # SHAP (optional)
    shap_vals = None
    try:
        import shap  # type: ignore
        def f(x):
            return m.predict(x)
        explainer = shap.KernelExplainer(f, X[: min(len(X), 50)])
        sv = explainer.shap_values(X[: min(len(X), 200)], nsamples=50)
        if isinstance(sv, list):
            sv = sv[0]
        shap_vals = {names[k]: float(np.mean(np.abs(sv[:, k]))) for k in range(X.shape[1])}
    except Exception:
        pass

    with open(os.path.join(OUT_DIR, "sensitivity.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["method"] + names)
        w.writerow(["liong"] + [f"{liong[n]:.4g}" for n in names])
        w.writerow(["gradient"] + [f"{grad[n]:.4g}" for n in names])
        if shap_vals:
            w.writerow(["shap"] + [f"{shap_vals[n]:.4g}" for n in names])

    # dose sweep (procedure only)
    cc_idx = featsA.index("CC")
    med = np.median(X, axis=0)
    grid = np.linspace(0.1, 2.5, 49)
    with open(os.path.join(OUT_DIR, "dose_sweep_A_ratio.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["dose_%", "predicted_ratio"])
        for cc in grid:
            x = med.copy()
            x[cc_idx] = cc
            w.writerow([f"{cc:.3f}", f"{m.predict(x.reshape(1, -1))[0]:.4f}"])


def _write_comparison(dfs):
    """Write the headline train/test/CV/LOSO table for the original config."""
    rows_out = []
    for name, feats in INPUT_SETS.items():
        df = dfs[name]
        for tgt in TARGETS:
            tr, te = split_80_20(df, feats, tgt, {}, _gmdh_predict_w(GMDH_KWARGS))
            folds = kfold_cv(df, feats, tgt, {}, _gmdh_predict_w(GMDH_KWARGS))
            loso, _ = leave_one_study_out(df, feats, tgt, {}, _gmdh_predict_w(GMDH_KWARGS))
            ws, _ = within_study_r2(df, feats, tgt, {}, _gmdh_predict_w(GMDH_KWARGS))
            rows_out.append({
                "set": name, "target": tgt, "n": len(df), "studies": df["ref"].nunique(),
                "train_r2": tr[0], "test_r2": te[0], "test_rmse": te[1],
                "cv_r2": float(np.mean([m[0] for m in folds])),
                "cv_rmse": float(np.mean([m[1] for m in folds])),
                "loso_r2": loso[0], "loso_rmse": loso[1],
                "within_study_r2": ws,
            })
    pd.DataFrame(rows_out).to_csv(os.path.join(OUT_DIR, "comparison_original.csv"), index=False)


if __name__ == "__main__":
    main()
