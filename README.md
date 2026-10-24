# Chitosan-Based Expansive Soil Stabilization | AI & White-Box GMDH

An independent research project integrating **sustainable geotechnical engineering, chitosan-based soil stabilization, and interpretable artificial intelligence** to investigate and predict the unconfined compressive strength (UCS) of treated expansive and cohesive soils.

The project combines a curated literature-derived dataset with a dedicated **Group Method of Data Handling (GMDH)** modeling framework to analyze the nonlinear relationships between chitosan dosage, soil properties, curing conditions, and strength development.

## Research Objectives

- **Dataset Development:** Consolidate experimental findings from published studies on chitosan-treated cohesive and expansive soils into a structured, traceable dataset.
- **AI-Based Prediction:** Develop an interpretable GMDH regression framework for modeling the mechanical response of treated soils.
- **Parameter Analysis:** Investigate the relationships between treatment conditions, intrinsic soil properties, and UCS.
- **Dosage Investigation:** Use model-based analysis to explore chitosan dosage ranges and their relationship with predicted soil strength.
- **Reproducible Research:** Establish a documented workflow for data processing, validation, visualization, and model evaluation.

## Key Features

- 167 records compiled from 17 peer-reviewed studies.
- 48-column primary dataset with variable definitions and units.
- Value-level data provenance and source traceability.
- Dedicated GMDH regression implementation.
- Configurable predictor sets for comparative modeling.
- Strength normalization relative to untreated soil.
- Automated model training, validation, sensitivity analysis, and dosage sweeps.
- Reproducible dataset checks and figure generation.

## AI Modeling Framework

The project employs **white-box machine learning** through the Group Method of Data Handling (GMDH), an approach designed to construct explicit mathematical relationships from experimental data.

Unlike models that provide predictions without a directly accessible analytical expression, a GMDH framework can represent learned relationships through polynomial equations, supporting mathematical inspection and investigation of variable interactions.

### Predictor Configurations

| Configuration | Input Variables |
|---|---|
| **A — Core Features** | Chitosan content (`CC`), curing time (`CTi`), plasticity index (`PI`) |
| **B — Extended Features** | Configuration A + compaction moisture content (`WMC`) |
| **C — Baseline Strength** | Configuration B + untreated UCS (`q0`) |

### Prediction Targets

The framework defines normalized strength targets:

\[
\mathrm{ratio}=\frac{\mathrm{UCS}_{\mathrm{treated}}}{q_0}
\]

\[
\mathrm{logratio}=\ln\left(\frac{\mathrm{UCS}_{\mathrm{treated}}}{q_0}\right)
\]

Here, \(\mathrm{UCS}_{\mathrm{treated}}\) denotes the treated soil's unconfined compressive strength, while \(q_0\) represents the untreated soil's baseline UCS.

These normalized targets support analysis of strength changes relative to the initial soil condition.

## Dataset & Research Documentation

The primary dataset contains 167 records from 17 peer-reviewed studies. Each record is associated with documented data-processing conventions and provenance labels identifying values extracted from tables, extracted from text, estimated from figures, or derived through calculations.

The repository also includes an additional dataset covering modified-chitosan studies.

| File | Purpose |
|---|---|
| `data/sheet1_cleaned.csv` | Primary chitosan dataset |
| `data/sheet3_modified_chitosan.csv` | Modified-chitosan study records |
| `data/study_summary.csv` | Study-level summary and record counts |
| `data/sources.csv` | Bibliographic metadata and source identifiers |
| `data/data_used.csv` | Records designated for modeling |
| `data/data_dictionary.md` | Variable definitions, units, and provenance |
| `docs/technical_report.md` | Detailed methodology and technical documentation |

## Reproducible Workflow

### Installation

```bash
pip install -r requirements.txt
```

### Dataset Validation

```bash
python scripts/dataset_checks.py
```

### Generate Research Figures

```bash
python scripts/make_figures.py
```

### Execute the Modeling Pipeline

```bash
python model/model.py
```

The modeling pipeline is configured through `model/config.py`. Generated outputs are written to `output/`, which is excluded from version control.

## Repository Structure

```text
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
├── data/
│   ├── sheet1_cleaned.csv
│   ├── sheet3_modified_chitosan.csv
│   ├── study_summary.csv
│   ├── sources.csv
│   ├── data_used.csv
│   └── data_dictionary.md
├── model/
│   ├── gmdh.py
│   ├── model.py
│   └── config.py
├── scripts/
│   ├── make_figures.py
│   └── dataset_checks.py
├── figures/
├── docs/
│   └── technical_report.md
└── output/
```

## Scientific Context

This project investigates the application of interpretable AI to a sustainable soil stabilization problem using chitosan as the treatment material.

The broader research context includes published work on white-box machine learning for biopolymer-stabilized soils. The present project focuses on its own chitosan-derived dataset, modeling implementation, and analysis.

The aim is to investigate whether interpretable mathematical models can capture the relationships between soil characteristics, treatment parameters, and strength development, while providing a reproducible basis for further geotechnical research.

## License & Citation

The project code is distributed under the MIT License.

The dataset contains information extracted from published scientific literature. Users reusing the data should cite the original studies listed in `data/sources.csv`, alongside the repository where appropriate.

Repository citation metadata is maintained in `CITATION.cff`.

---

**Research Areas:** Artificial Intelligence · White-Box Machine Learning · GMDH Regression · Chitosan-Based Soil Stabilization · Expansive Soils · Sustainable Geotechnical Engineering

