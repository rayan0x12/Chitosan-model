"""Central configuration for the GMDH modelling framework.

Everything the pipeline needs — input sets, targets, seeds, paths, units and
plotting style — lives here so ``model.py`` and ``make_figures.py`` share one
source of truth.  Edit this file to add a study, change a seed or re-define an
input set; nothing else needs to change.
"""

from __future__ import annotations

import os

# ---------------------------------------------------------------------------
# Paths (relative to this file -> repo root)
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(HERE)
DATA_DIR = os.path.join(REPO_ROOT, "data")
FIG_DIR = os.path.join(REPO_ROOT, "figures")
OUT_DIR = os.path.join(REPO_ROOT, "output")  # gitignored — run artefacts only

SHEET1 = os.path.join(DATA_DIR, "sheet1_cleaned.csv")

# ---------------------------------------------------------------------------
# Input sets and targets
# ---------------------------------------------------------------------------
# Feature short-name -> the human label + unit used on every axis/table.
FEATURE_LABELS = {
    "CC": "Chitosan content (%)",
    "CTi": "Curing time (d)",
    "PI": "Plasticity index (-)",
    "WMC": "Initial moisture content (%)",
    "q0": "Untreated UCS q0 (MPa)",
}

# Nested input sets.  A is the largest row set; B adds WMC; C adds q0 as a
# predictor (C uses the same rows as B because q0 is already required for the
# normalised target).
INPUT_SETS = {
    "A": ["CC", "CTi", "PI"],
    "B": ["CC", "CTi", "PI", "WMC"],
    "C": ["CC", "CTi", "PI", "WMC", "q0"],
}

# Targets.  "ratio" = UCS_treated / q0 (dimensionless); "logratio" is its log.
TARGETS = ["ratio", "logratio"]
TARGET_LABELS = {
    "ratio": "UCS treated / q0 (-)",
    "logratio": "ln(UCS treated / q0) (-)",
}

# Mapping from a feature short-name to the column key produced by the loader.
FEATURE_KEY = {"CC": "cc", "CTi": "cti", "PI": "pi", "WMC": "wmc", "q0": "q0"}

# ---------------------------------------------------------------------------
# GMDH hyper-parameters
# ---------------------------------------------------------------------------
GMDH_KWARGS = {
    "alpha": 0.6,          # selection threshold (1 + alpha)
    "ridge": 10.0,         # Tikhonov ridge on standardised features
    "max_layers": 1,       # single-layer GMDH (stable on small data)
    "val_frac": 0.30,      # internal neuron-selection split
    "max_neurons": 8,      # cap on neurons carried forward
}

# ---------------------------------------------------------------------------
# Fixed seeds (reproducibility)
# ---------------------------------------------------------------------------
SEEDS = {
    "gmdh": 1,        # internal neuron-selection split
    "split": 42,      # 80/20 random split
    "cv": 0,          # k-fold cross-validation (k = 5)
    "bootstrap": 123, # bootstrap resampling
}

KFOLD = 5
N_BOOTSTRAP = 1000

# Liong et al. (2000) perturbation factors for classical sensitivity analysis.
CHANGE_FACTORS = [0.8, 0.9, 0.95, 0.99]

# ---------------------------------------------------------------------------
# Plotting style (shared by every figure)
# ---------------------------------------------------------------------------
PLOT_RC = {
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "font.size": 9,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
}
PALETTE = "tab10"
DPI = 300
