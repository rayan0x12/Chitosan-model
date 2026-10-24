# Data dictionary

The primary dataset is `data/sheet1_cleaned.csv` (167 rows × 48 columns), one row
per reported measurement. Derived, model-ready views are `data/data_used.csv`
(139 rows fed to the model), `data/sources.csv` (one row per study) and
`data/study_summary.csv`. `data/sheet3_modified_chitosan.csv` records the
modified-chitosan studies that inform the planned laboratory programme but are
out of scope for the UCS model.

## Provenance labels

Every row carries a `Value_Source` that states *where the number came from*. For
summaries and filtering, each row is classified into exactly one of four labels
(see `scripts/dataset_checks.py`):

| Label | Meaning | Rows |
|---|---|---|
| `table` | read directly from a table in the paper | 101 |
| `text` | read from the article text (sometimes marked "approximate") | 30 |
| `graph-estimated` | read from a figure (axis/data labels) | 28 |
| `derived` | computed arithmetically from other values (e.g. `UCS_treated = q0 × (1 + ΔUCS/100)`) | 8 |

## Unit and value conventions

- **Units are per-column** and given in the table below. Stress is kPa for the
  raw dataset and MPa only inside the model feature `q0` (see `model/config.py`).
- **Missing values** are the empty string (rendered `NA` where a value is
  genuinely absent/not reported). A value is *never* imputed or invented.
- **Derived values** are flagged as such in `Value_Source` and in `Notes`; the
  arithmetic is written out so it can be re-checked.
- **Ranges vs points**: cells that print a range (e.g. `24-25`) are left
  non-numeric and are *not* averaged; the model's `parse_num` returns `None`
  for them.

## `sheet1_cleaned.csv` columns

### Identification

| Column | Description |
|---|---|
| `Ref_ID` | Study key, `R01`–`R17`. |
| `Authors` | Authors as printed in the paper. |
| `Year` | Publication year. |
| `Title` | Article title. |
| `Journal` | Journal / proceedings and volume (issue):pages where available. |
| `DOI_URL` | DOI or stable URL; `(NOT OPEN)` marks an inaccessible full text. |
| `Soil_Name` | Soil description as reported (e.g. "Kaolinite silt", "Bhavnagar black cotton soil"). |
| `USCS_Class` | Unified Soil Classification System class (e.g. `CL`, `CH`, `MH`). |
| `Clay_Mineralogy` | Reported mineralogy (e.g. montmorillonite, quartz-dominant silt). |
| `Country` | Where the soil was sampled / tests performed. |

### Untreated soil properties

| Column | Unit | Description |
|---|---|---|
| `LL_untreated` | % | Liquid limit. |
| `PL_untreated` | % | Plastic limit. |
| `PI_untreated` | — | Plasticity index (LL − PL). |
| `FSI_untreated_%` | % | Free swell index. |
| `Swell_Percent_untreated` | % | Swell (volume change) percent. |
| `Swell_Pressure_untreated_kPa` | kPa | Swell pressure. |
| `UCS_untreated_kPa_q0` | kPa | Untreated UCS; the model's `q0`. |
| `CBR_untreated_soaked_%` | % | California bearing ratio (soaked). |
| `CBR_untreated_unsoaked_%` | % | CBR (unsoaked). |
| `MDD_g_cm3` | g/cm³ | Maximum dry density. |
| `OMC_%` | % | Optimum moisture content. |

### Treatment (chitosan)

| Column | Unit | Description |
|---|---|---|
| `Chitosan_Content_%` | % | Chitosan dose by dry weight of soil. |
| `Chitosan_Type` | — | Raw / nano / quaternized / composite, source, CAS where given. |
| `DDA_%` | % | Degree of deacetylation (often not reported → `NA`). |
| `Molecular_Weight` | — | Molecular weight, when quantified. |
| `Solvent` | — | Dissolution medium (e.g. acetic acid, distilled water, none). |
| `Form` | — | Physical form (powder, solution). |
| `Mixing_Method` | — | How the chitosan was combined with the soil. |
| `Initial_Moisture_Content_%` | % | Compaction moisture (the model's `WMC`). |
| `Curing_Time_days` | d | Curing period (the model's `CTi`). |
| `Curing_Temperature_C` | °C | Curing temperature, when reported numerically. |
| `Curing_Condition` | — | Sealed / humid / air-dried, etc. |

### Treated soil properties

| Column | Unit | Description |
|---|---|---|
| `LL_treated`, `PL_treated`, `PI_treated` | %, %, — | Post-treatment Atterberg limits. |
| `FSI_treated_%` | % | Free swell index after treatment. |
| `Swell_Percent_treated` | % | Swell percent after treatment. |
| `Swell_Pressure_treated_kPa` | kPa | Swell pressure after treatment. |
| `UCS_treated_kPa` | kPa | Treated UCS (the model output). |
| `CBR_treated_soaked_%` / `CBR_treated_unsoaked_%` | % | Treated CBR. |
| `Secant_Modulus_MPa` | MPa | Secant (Young's) modulus. |
| `WetDry_FT_Cycles_Number` | — | Number of wet-dry / freeze-thaw cycles (0 = un-cycled). |
| `Properties_after_cycles` | — | Reported post-cycle property change. |

### Provenance / QA

| Column | Description |
|---|---|
| `Value_Source` | Free-text provenance (see labels above). |
| `Page_Table_Figure` | Location in the paper (e.g. "Table 2", "Figure 6b"). |
| `Row_Status` | `complete` or `incomplete`. |
| `Notes` | Extraction caveats, derivation formulae, flags. |

## `data_used.csv` columns (model input)

| Column | Unit | Description |
|---|---|---|
| `Ref_ID` | — | Study key. |
| `Study` | — | Treatment label (R13 is split into `R13 (CTS)` and `R13 (CTS + calsil)`). |
| `PI` | — | Plasticity index. |
| `CC` | % | Chitosan dose. |
| `CTi` | d | Curing time. |
| `WMC` | % | Compaction moisture (may be `NA` for set A rows). |
| `q0_kPa` | kPa | Untreated UCS. |
| `UCS_treated_kPa` | kPa | Treated UCS. |
| `UCS_q0_ratio` | — | `UCS_treated / q0` (model target "ratio"). |
| `log_UCS_q0` | — | `ln(UCS_treated / q0)` (model target "logratio"). |
| `Chitosan_Type`, `Solvent`, `DDA_%` | — | Treatment metadata. |
| `Value_Source`, `Row_Status` | — | Provenance and QA. |

## `sources.csv` columns

`Ref_ID, Authors, Year, Title, Journal, DOI, Access, Rows` — one row per study
(R01–R17). `Access` is `full text` or `not accessed`; `Rows` is the number of
rows that study contributes to `sheet1_cleaned.csv`.
