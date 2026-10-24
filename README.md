# Chitosan-stabilized expansive clay — dataset and white-box GMDH framework

A documented dataset and modelling framework for a graduation project on
**chitosan as a green stabilizer for expansive/cohesive clay**. The project
replicates the methodology of a published white-box *Group Method of Data
Handling* (GMDH) study on xanthan-gum-stabilized soils, but with **chitosan**
as the biopolymer.

The repository contains three things:

1. **A curated dataset** — 167 rows transcribed from 17 peer-reviewed studies
   (`R01`–`R17`) of chitosan-treated cohesive/expansive soils, with a
   column-level data dictionary and a per-value provenance label
   (table / text / graph-estimated / derived).
2. **A data-extraction protocol** — search strategy, inclusion/exclusion
   criteria, paywalled-paper handling, duplicate detection, figure-digitization
   QA and derived-value rules (`docs/technical_report.md` §2).
3. **A white-box GMDH modelling framework** — a self-contained, documented GMDH
   regressor and a one-command train / validate / sensitivity / dose-sweep
   pipeline with fixed seeds.

> **Status: model calibration and validation — in progress.** The model is
> provided as a reproducible *framework*; this repository intentionally ships
> no fitted equations, performance scores or optimum doses.

## Quick start

```bash
# 1. install dependencies
pip install -r requirements.txt

# 2. integrity-check the dataset (exits 0 only if all checks pass)
python scripts/dataset_checks.py

# 3. regenerate the data figures -> figures/*.png
python scripts/make_figures.py

# 4. run the modelling pipeline -> output/ (gitignored)
python model/model.py
```

All scripts run from the repository root with a single command and use the
fixed seeds in `model/config.py`.

## Repository layout

```
├── README.md
├── LICENSE                     # MIT (code); data = extracted facts — cite sources
├── CITATION.cff                # citation metadata (fill in the TODOs)
├── requirements.txt
├── .gitignore
├── data/
│   ├── sheet1_cleaned.csv          # 167 rows × 48 cols (primary dataset)
│   ├── sheet3_modified_chitosan.csv# modified-chitosan studies (M01–M05)
│   ├── study_summary.csv           # one row per study, row counts
│   ├── sources.csv                 # Ref_ID, authors, year, title, journal, DOI, access, rows
│   ├── data_used.csv               # 139 rows fed to the model
│   └── data_dictionary.md          # columns, units, provenance labels
├── model/
│   ├── gmdh.py                     # self-contained GMDH regressor (documented)
│   ├── model.py                    # one-command pipeline (train/validate/sensitivity/dose-sweep)
│   └── config.py                   # input sets A/B/C, targets, seeds, style
├── scripts/
│   ├── make_figures.py             # regenerates every figure in figures/
│   └── dataset_checks.py           # integrity + reproducibility checks
├── figures/                        # 300 dpi PNG data figures (see docs §4)
├── docs/
│   └── technical_report.md         # journal-style report
└── output/                         # gitignored run artefacts
```

## Input sets and targets

The model uses three nested input sets and a normalized target that removes the
per-study baseline strength:

| Set | Inputs |
|---|---|
| A | chitosan dose `CC`, curing time `CTi`, plasticity index `PI` |
| B | A + compaction moisture `WMC` |
| C | B + untreated UCS `q0` (as a predictor) |

Targets: `ratio = UCS_treated / q0` and `logratio = ln(UCS_treated / q0)`.

## Data and citation

The code is MIT-licensed. The files under `data/` contain **extracted facts**
transcribed from published papers; if you reuse them, cite the original studies
listed in `data/sources.csv` in addition to this repository. See `CITATION.cff`
for the repository citation (replace the author/contact placeholders).

## Note on scope

This is an open-source data + modelling template for a student project. The
modelling sections are written as a protocol; calibration numbers are produced
locally by `model/model.py` but are not part of the repository.
