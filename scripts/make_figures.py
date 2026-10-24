"""Generate the repository's data figures.

Run from the repo root with one command::

    python scripts/make_figures.py

Every figure is a *data* figure (no model output): row counts, distributions,
correlations, VIF, dose/curing curves, UCS/q0 relationships, boxplots, a
coverage heatmap and an input-space scatter.  All numbers are computed from
``data/sheet1_cleaned.csv``; nothing is hard-coded.  Output is written to
``figures/`` at 300 dpi, white background, labelled axes with units, one
categorical palette (tab10) throughout.
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
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(HERE)
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from model.config import PLOT_RC, PALETTE, DPI, FIG_DIR  # noqa: E402

plt.rcParams.update(PLOT_RC)
os.makedirs(FIG_DIR, exist_ok=True)

SHEET1 = os.path.join(REPO_ROOT, "data", "sheet1_cleaned.csv")

# Numeric columns we care about, in display order.
NUM_COLS = [
    "PI_untreated", "Chitosan_Content_%", "Curing_Time_days",
    "Initial_Moisture_Content_%", "UCS_untreated_kPa_q0", "UCS_treated_kPa",
]

# Additional numeric columns parsed for the coverage heatmap.
COVERAGE_COLS = [
    "LL_untreated", "PL_untreated", "FSI_untreated_%", "DDA_%",
    "CBR_treated_unsoaked_%", "Secant_Modulus_MPa",
]

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


def load_all():
    """Load every row of sheet1 with numeric columns parsed."""
    with open(SHEET1, newline="", encoding="utf-8-sig") as fh:
        rd = list(csv.reader(fh))
    hdr = rd[0]
    idx = {h: i for i, h in enumerate(hdr)}
    rows = []
    for r in rd[1:]:
        rec = {
            "Ref_ID": r[idx["Ref_ID"]].strip(),
            "Soil_Name": r[idx["Soil_Name"]].strip(),
            "USCS": r[idx["USCS_Class"]].strip(),
            "Country": r[idx["Country"]].strip(),
            "Chitosan_Type": r[idx["Chitosan_Type"]].strip(),
            "Solvent": r[idx["Solvent"]].strip(),
            "DDA": r[idx["DDA_%"]].strip(),
            "Value_Source": r[idx["Value_Source"]].strip(),
            "Row_Status": r[idx["Row_Status"]].strip(),
            "Notes": r[idx["Notes"]].strip(),
        }
        for c in NUM_COLS:
            rec[c] = parse_num(r[idx[c]])
        for c in COVERAGE_COLS:
            rec[c] = parse_num(r[idx[c]])
        rec["q0"] = rec["UCS_untreated_kPa_q0"]
        rec["ucs"] = rec["UCS_treated_kPa"]
        rows.append(rec)
    return rows


def provenance(vs: str) -> str:
    """Normalise a free-text Value_Source to one of four provenance labels."""
    s = vs.lower()
    if s.startswith("derived"):
        return "derived"
    if "graph" in s:
        return "graph-estimated"
    if "table" in s:
        return "table"
    return "text"


def _save(fig, name):
    path = os.path.join(FIG_DIR, name)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------------------
# 1. rows per study
# ---------------------------------------------------------------------------
def fig_rows_per_study(rows):
    df = pd.DataFrame(rows)
    counts = df.groupby("Ref_ID").size()
    counts = counts.reindex(sorted(counts.index, key=lambda k: int(k[1:])))
    fig, ax = plt.subplots(figsize=(7, 3.4))
    cols = [plt.get_cmap(PALETTE)(i % 10) for i in range(len(counts))]
    ax.bar(counts.index, counts.values, color=cols)
    ax.set_xlabel("Study (Ref_ID)")
    ax.set_ylabel("Number of rows")
    ax.set_title("Rows per study in sheet1 (n = %d)" % len(df))
    ax.tick_params(axis="x", rotation=45)
    for i, v in enumerate(counts.values):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=7)
    fig.tight_layout()
    _save(fig, "fig01_rows_per_study.png")


# ---------------------------------------------------------------------------
# 2. histograms with normal fits
# ---------------------------------------------------------------------------
def fig_histograms(rows):
    df = pd.DataFrame(rows)
    fig, axes = plt.subplots(2, 3, figsize=(9, 5.2))
    axes = axes.ravel()
    for ax, c in zip(axes, NUM_COLS):
        v = df[c].dropna().astype(float)
        if len(v) == 0:
            ax.axis("off")
            continue
        ax.hist(v, bins=20, density=True, color="tab:blue", alpha=0.55,
                edgecolor="white")
        mu, sd = v.mean(), v.std()
        x = np.linspace(v.min(), v.max(), 200)
        ax.plot(x, stats.norm.pdf(x, mu, sd), "r-", lw=1.2,
                label="normal fit")
        ax.set_title(c, fontsize=8)
        ax.legend(fontsize=6)
        ax.set_xlabel("value")
        ax.set_ylabel("density")
    fig.suptitle("Distributions of key variables with normal fits")
    fig.tight_layout()
    _save(fig, "fig02_histograms.png")


# ---------------------------------------------------------------------------
# 3. correlation heatmap
# ---------------------------------------------------------------------------
def fig_correlation(rows):
    df = pd.DataFrame(rows)
    cols = ["PI_untreated", "Chitosan_Content_%", "Curing_Time_days",
            "Initial_Moisture_Content_%", "UCS_untreated_kPa_q0", "UCS_treated_kPa"]
    labels = ["PI", "Dose", "Curing", "WMC", "q0", "UCS treated"]
    sub = df[cols].dropna().astype(float)
    corr = sub.corr()
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(corr.values, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)
    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(j, i, f"{corr.values[i, j]:.2f}", ha="center", va="center",
                    fontsize=7, color="black")
    fig.colorbar(im, ax=ax, label="Pearson r")
    ax.set_title(f"Pearson correlation (n = {len(sub)} complete rows)")
    fig.tight_layout()
    _save(fig, "fig03_correlation_heatmap.png")


# ---------------------------------------------------------------------------
# 4. VIF
# ---------------------------------------------------------------------------
def _vif(df, cols):
    X = df[cols].dropna().astype(float).values
    vifs = {}
    for k, name in enumerate(cols):
        y = X[:, k]
        rest = np.delete(X, k, axis=1)
        rest = np.column_stack([np.ones(len(y)), rest])
        coef, *_ = np.linalg.lstsq(rest, y, rcond=None)
        yp = rest @ coef
        ss_res = float(np.sum((y - yp) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
        vifs[name] = float("inf") if r2 >= 1.0 else 1.0 / (1.0 - r2)
    return vifs


def fig_vif(rows):
    df = pd.DataFrame(rows)
    cols = ["PI_untreated", "Chitosan_Content_%", "Curing_Time_days",
            "Initial_Moisture_Content_%"]
    labels = ["PI", "Dose", "Curing", "WMC"]
    vifs = _vif(df, cols)
    fig, ax = plt.subplots(figsize=(5, 3.4))
    vals = [vifs[c] for c in cols]
    ax.bar(labels, vals, color="tab:blue", alpha=0.8)
    ax.axhline(5, color="r", ls="--", lw=1, label="threshold = 5")
    ax.set_ylabel("Variance inflation factor (VIF)")
    ax.set_title("Multicollinearity of the four predictors")
    ax.legend()
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    _save(fig, "fig04_vif.png")


# ---------------------------------------------------------------------------
# 5. pairplot (input space)
# ---------------------------------------------------------------------------
def fig_pairplot(rows):
    df = pd.DataFrame(rows)
    cols = ["PI_untreated", "Chitosan_Content_%", "Curing_Time_days",
            "Initial_Moisture_Content_%"]
    labels = ["PI", "Dose (%)", "Curing (d)", "WMC (%)"]
    sub = df[["Ref_ID"] + cols].dropna()
    refs = sorted(sub["Ref_ID"].unique(), key=lambda k: int(k[1:]))
    cmap = plt.get_cmap(PALETTE)
    color_of = {r: cmap(i % 10) for i, r in enumerate(refs)}
    n = len(cols)
    fig, axes = plt.subplots(n, n, figsize=(9, 8))
    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            if i == j:
                ax.hist(sub[cols[i]], bins=20, color="tab:blue", alpha=0.6)
                ax.set_xlabel(labels[i], fontsize=7)
                ax.set_ylabel("count", fontsize=7)
            else:
                for r in refs:
                    d = sub[sub["Ref_ID"] == r]
                    ax.scatter(d[cols[j]], d[cols[i]], s=10, alpha=0.6,
                               color=color_of[r], label=r if (i == 1 and j == 0) else None)
                ax.set_xlabel(labels[j], fontsize=7)
                ax.set_ylabel(labels[i], fontsize=7)
            ax.tick_params(labelsize=6)
    handles = [plt.Line2D([0], [0], marker="o", ls="", color=color_of[r], label=r)
               for r in refs]
    fig.legend(handles=handles, loc="lower right", fontsize=6, ncol=2,
               title="Study", title_fontsize=7)
    fig.suptitle("Pairwise input space (rows with all four inputs present)")
    fig.tight_layout()
    _save(fig, "fig05_pairplot.png")


# ---------------------------------------------------------------------------
# 6-9. dose / curing curves and ratios
# ---------------------------------------------------------------------------
def _study_groups(rows):
    """Studies with dose + UCS (raw treated), grouped by study and curing time."""
    from collections import defaultdict
    g = defaultdict(list)
    for r in rows:
        if r["Ref_ID"] in ("R08", "R11"):
            continue
        if r["Row_Status"].lower() == "incomplete":
            continue
        if None in (r["Chitosan_Content_%"], r["UCS_treated_kPa"]):
            continue
        g[r["Ref_ID"]].append(r)
    return g


def fig_ucs_vs_dose(rows):
    g = _study_groups(rows)
    refs = sorted(g, key=lambda k: int(k[1:]))
    fig, axes = plt.subplots(3, 3, figsize=(10.5, 8))
    axes = axes.ravel()
    cmap = plt.get_cmap("viridis")
    for ax, ref in zip(axes, refs):
        rs = g[ref]
        ctis = sorted(set(r["Curing_Time_days"] for r in rs if r["Curing_Time_days"] is not None))
        for k, cti in enumerate(ctis):
            pts = sorted((r["Chitosan_Content_%"], r["UCS_treated_kPa"])
                         for r in rs if r["Curing_Time_days"] == cti)
            cc = [p[0] for p in pts]
            u = [p[1] for p in pts]
            col = cmap(k / max(1, len(ctis) - 1)) if len(ctis) > 1 else "tab:blue"
            ax.plot(cc, u, "-o", ms=3.5, lw=1, color=col, label=f"{cti:g} d")
        ax.set_title(ref, fontsize=8)
        ax.set_xlabel("Chitosan dose (%)")
        ax.set_ylabel("UCS treated (kPa)")
        ax.tick_params(labelsize=7)
    for ax in axes[len(refs):]:
        ax.axis("off")
    axes[0].legend(fontsize=6, ncol=2, title="Curing time")
    fig.suptitle("UCS treated vs chitosan dose, per study (lines = curing time)")
    fig.tight_layout()
    _save(fig, "fig06_ucs_vs_dose.png")


def fig_ucs_vs_curing(rows):
    g = _study_groups(rows)
    refs = sorted(g, key=lambda k: int(k[1:]))
    fig, axes = plt.subplots(3, 3, figsize=(10.5, 8))
    axes = axes.ravel()
    cmap = plt.get_cmap("viridis")
    for ax, ref in zip(axes, refs):
        rs = g[ref]
        doses = sorted(set(r["Chitosan_Content_%"] for r in rs))
        for k, dose in enumerate(doses):
            pts = sorted((r["Curing_Time_days"], r["UCS_treated_kPa"])
                         for r in rs if r["Chitosan_Content_%"] == dose
                         and r["Curing_Time_days"] is not None)
            ct = [p[0] for p in pts]
            u = [p[1] for p in pts]
            col = cmap(k / max(1, len(doses) - 1)) if len(doses) > 1 else "tab:blue"
            ax.plot(ct, u, "-o", ms=3.5, lw=1, color=col, label=f"{dose:g} %")
        ax.set_title(ref, fontsize=8)
        ax.set_xlabel("Curing time (d)")
        ax.set_ylabel("UCS treated (kPa)")
        ax.tick_params(labelsize=7)
    for ax in axes[len(refs):]:
        ax.axis("off")
    axes[0].legend(fontsize=6, ncol=2, title="Dose")
    fig.suptitle("UCS treated vs curing time, per study (lines = dose)")
    fig.tight_layout()
    _save(fig, "fig07_ucs_vs_curing.png")


def fig_ratio_dose_pi(rows):
    """UCS/q0 vs dose (left) and vs PI (right), coloured by study."""
    # model-eligible rows (set A subset) to compute the ratio
    from collections import defaultdict
    def label(r):
        return r["Ref_ID"]
    g = defaultdict(list)
    for r in rows:
        if r["Ref_ID"] in ("R08", "R11"):
            continue
        if r["Row_Status"].lower() == "incomplete":
            continue
        if None in (r["Chitosan_Content_%"], r["UCS_treated_kPa"],
                    r["UCS_untreated_kPa_q0"], r["PI_untreated"]):
            continue
        if r["UCS_untreated_kPa_q0"] <= 0 or r["UCS_treated_kPa"] <= 0:
            continue
        g[r["Ref_ID"]].append(r)
    refs = sorted(g, key=lambda k: int(k[1:]))
    cmap = plt.get_cmap(PALETTE)
    color_of = {r: cmap(i % 10) for i, r in enumerate(refs)}

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for r in refs:
        d = g[r]
        dose = [x["Chitosan_Content_%"] for x in d]
        ratio = [x["UCS_treated_kPa"] / x["UCS_untreated_kPa_q0"] for x in d]
        axes[0].scatter(dose, ratio, s=16, alpha=0.7, color=color_of[r], label=r)
        pi = [x["PI_untreated"] for x in d]
        axes[1].scatter(pi, ratio, s=16, alpha=0.7, color=color_of[r])
    axes[0].set_xlabel("Chitosan dose (%)")
    axes[0].set_ylabel("UCS treated / q0 (-)")
    axes[0].set_title("UCS/q0 vs dose")
    axes[1].set_xlabel("Plasticity index (-)")
    axes[1].set_ylabel("UCS treated / q0 (-)")
    axes[1].set_title("UCS/q0 vs PI")
    axes[0].legend(fontsize=6, ncol=2, title="Study", title_fontsize=7)
    fig.tight_layout()
    _save(fig, "fig08_ratio_dose_pi.png")


# ---------------------------------------------------------------------------
# 10. boxplots by chitosan form and curing regime
# ---------------------------------------------------------------------------
def _form(r):
    s = r["Solvent"].lower()
    t = r["Chitosan_Type"].lower()
    if "nano" in t:
        return "nano"
    if "acetic" in s:
        return "acid-dissolved"
    if "solution" in t:
        return "solution"
    if "water" in s:
        return "water-dissolved"
    if "none" in s or "dry" in s or s in ("", "na"):
        return "dry-mixed"
    return "other"


def _regime(cti):
    if cti is None:
        return "NA"
    if cti <= 1:
        return "0–1 d"
    if cti <= 7:
        return "2–7 d"
    if cti <= 28:
        return "8–28 d"
    if cti <= 90:
        return "29–90 d"
    return ">90 d"


def fig_boxplots(rows):
    recs = []
    for r in rows:
        if r["Ref_ID"] in ("R08", "R11"):
            continue
        if r["Row_Status"].lower() == "incomplete":
            continue
        if None in (r["UCS_treated_kPa"], r["UCS_untreated_kPa_q0"]):
            continue
        if r["UCS_untreated_kPa_q0"] <= 0:
            continue
        recs.append({
            "ratio": r["UCS_treated_kPa"] / r["UCS_untreated_kPa_q0"],
            "form": _form(r),
            "regime": _regime(r["Curing_Time_days"]),
        })
    df = pd.DataFrame(recs)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))

    forms = sorted(df["form"].unique())
    data = [df[df["form"] == f]["ratio"] for f in forms]
    bp = axes[0].boxplot(data, tick_labels=forms, patch_artist=True)
    for patch, i in zip(bp["boxes"], range(len(forms))):
        patch.set_facecolor(plt.get_cmap(PALETTE)(i % 10))
        patch.set_alpha(0.6)
    axes[0].set_xlabel("Chitosan form / solvent")
    axes[0].set_ylabel("UCS treated / q0 (-)")
    axes[0].set_title("UCS/q0 by chitosan form")
    axes[0].tick_params(axis="x", rotation=30, labelsize=7)

    regimes = ["0–1 d", "2–7 d", "8–28 d", "29–90 d", ">90 d", "NA"]
    data = [df[df["regime"] == rg]["ratio"] for rg in regimes]
    data = [d for d in data if len(d) > 0]
    regs = [rg for rg, d in zip(regimes, data) if len(d) > 0]
    bp = axes[1].boxplot(data, tick_labels=regs, patch_artist=True)
    for patch, i in zip(bp["boxes"], range(len(regs))):
        patch.set_facecolor(plt.get_cmap(PALETTE)(i % 10))
        patch.set_alpha(0.6)
    axes[1].set_xlabel("Curing regime")
    axes[1].set_ylabel("UCS treated / q0 (-)")
    axes[1].set_title("UCS/q0 by curing regime")
    axes[1].tick_params(axis="x", rotation=30, labelsize=7)

    fig.tight_layout()
    _save(fig, "fig09_boxplots_form_curing.png")


# ---------------------------------------------------------------------------
# 11. coverage heatmap (study x property)
# ---------------------------------------------------------------------------
def fig_coverage(rows):
    df = pd.DataFrame(rows)
    props = ["LL_untreated", "PL_untreated", "PI_untreated", "FSI_untreated_%",
             "UCS_untreated_kPa_q0", "Chitosan_Content_%", "DDA_%",
             "Initial_Moisture_Content_%", "Curing_Time_days", "UCS_treated_kPa",
             "CBR_treated_unsoaked_%", "Secant_Modulus_MPa"]
    short = ["LL", "PL", "PI", "FSI", "q0", "Dose", "DDA", "WMC", "Curing",
             "UCS_t", "CBR_t", "E50"]
    # use untreated + treated labels precisely
    full = ["LL_untreated", "PL_untreated", "PI_untreated", "FSI_untreated_%",
            "UCS_untreated_kPa_q0", "Chitosan_Content_%", "DDA_%",
            "Initial_Moisture_Content_%", "Curing_Time_days", "UCS_treated_kPa",
            "CBR_treated_unsoaked_%", "Secant_Modulus_MPa"]
    refs = sorted(df["Ref_ID"].unique(), key=lambda k: int(k[1:]))
    mat = np.zeros((len(refs), len(full)))
    for i, ref in enumerate(refs):
        sub = df[df["Ref_ID"] == ref]
        for j, c in enumerate(full):
            mat[i, j] = sub[c].notna().mean()
    fig, ax = plt.subplots(figsize=(9, 5))
    im = ax.imshow(mat, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(short)))
    ax.set_xticklabels(short, rotation=45, ha="right")
    ax.set_yticks(range(len(refs)))
    ax.set_yticklabels(refs)
    for i in range(len(refs)):
        for j in range(len(short)):
            v = mat[i, j]
            ax.text(j, i, f"{v:.0%}", ha="center", va="center", fontsize=6,
                    color="white" if v < 0.5 else "black")
    fig.colorbar(im, ax=ax, label="fraction of rows with a value")
    ax.set_title("Data coverage: fraction of non-missing values per study and property")
    fig.tight_layout()
    _save(fig, "fig10_coverage_heatmap.png")


# ---------------------------------------------------------------------------
# 12. input-space scatter by study (PI vs dose, sized/coloured)
# ---------------------------------------------------------------------------
def fig_input_space(rows):
    from collections import defaultdict
    g = defaultdict(list)
    for r in rows:
        if r["Ref_ID"] in ("R08", "R11"):
            continue
        if r["Row_Status"].lower() == "incomplete":
            continue
        if None in (r["PI_untreated"], r["Chitosan_Content_%"]):
            continue
        g[r["Ref_ID"]].append(r)
    refs = sorted(g, key=lambda k: int(k[1:]))
    cmap = plt.get_cmap(PALETTE)
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for i, ref in enumerate(refs):
        d = g[ref]
        pi = [x["PI_untreated"] for x in d]
        cc = [x["Chitosan_Content_%"] for x in d]
        ax.scatter(pi, cc, s=18, alpha=0.7, color=cmap(i % 10), label=ref)
    ax.set_xlabel("Plasticity index (-)")
    ax.set_ylabel("Chitosan dose (%)")
    ax.set_title("Input-space coverage by study (PI vs dose)")
    ax.legend(fontsize=6, ncol=2, title="Study", title_fontsize=7)
    fig.tight_layout()
    _save(fig, "fig11_input_space_by_study.png")


def main():
    rows = load_all()
    fig_rows_per_study(rows)
    fig_histograms(rows)
    fig_correlation(rows)
    fig_vif(rows)
    fig_pairplot(rows)
    fig_ucs_vs_dose(rows)
    fig_ucs_vs_curing(rows)
    fig_ratio_dose_pi(rows)
    fig_boxplots(rows)
    fig_coverage(rows)
    fig_input_space(rows)
    print("done ->", FIG_DIR)


if __name__ == "__main__":
    main()
