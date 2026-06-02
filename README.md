# Adolescent Exposure to Deadly Neighborhood Gun Violence and Asthma Burden in Early Adulthood

**Martz, C.D., Navas Nazario, A., & Gaydosh, L.**  

---

## Overview

This repository contains the R analysis code for a prospective study examining whether residential exposure to deadly neighborhood gun violence during adolescence is associated with incident asthma diagnosis, hospitalization, and functional limitation in early adulthood, and whether these associations differ by sex. Analyses use data from the Future of Families and Child Wellbeing Study (FFCWS), linked to geocoded gun violence incidence data from the Gun Violence Archive (GVA) and county violent crime data from the Uniform Crime Reports (UCR).

---

## Data Access

**This repository contains code only. No data are included.**

The FFCWS data used in this study are restricted-use and must be obtained directly from the data distributor:

- **FFCWS core and Wave 7 data**: Available via restricted-use contract from [Princeton's Center for Research on Child Wellbeing](https://ffcws.princeton.edu/data-and-documentation/data)
- **FFCWS–GVA linked data** (`ff_gva_15y_res1`): Available as a restricted-use appendage from Princeton/Columbia; see [documentation](https://ffcws.princeton.edu/sites/g/files/toruqf4356/files/documents/ff_gva_15y_res1_20190603.pdf)
- **Census and UCR data**: Linked administrative files distributed with the FFCWS restricted-use package

Secondary analyses were determined to be Not Human Subjects Research by the Institutional Review Board at the University of Central Florida.

---

## Script Structure

All analyses are contained in a single script: `asthma_gv_JAH_040226.R`

### Data Cleaning (Sections 1–13)

| Section | Description |
|---------|-------------|
| [1] | Load libraries |
| [2] | Data import |
| [3] | Merge FFCWS core, Wave 7, GVA, UCR, and census files |
| [4] | Date and exposure window variables (GVA coverage flags) |
| [5] | Gun violence exposure variables (binary, count, categorical, annualized rate) |
| [6] | Distance-banded counts (100m intervals, 0–1,600m) |
| [7] | Distance-weighted exposure (quadratic decay) |
| [8] | Time-windowed exposure counts (30-day increments) |
| [9] | Individual-level covariates (race/ethnicity, SES, smoking/vaping, BMI) |
| [10] | Neighborhood and household covariates (tract poverty, violent crime, household composition) |
| [11] | Housing quality and asthma environmental risk index (observer-rated, Wave 6) |
| [12] | Asthma outcome variables (incident diagnosis, hospitalization, functional limitation) |
| [13] | Standardize continuous predictors across full dataset |

### Analysis (Sections 14–26)

| Section | Description |
|---------|-------------|
| [14] | Helper functions for HC1-robust CI extraction (single term and linear combinations) |
| [15] | Create analysis-ready `df_reduced` with factor coding |
| [16] | Analytic sample preparation (`df_main_nomiss`; n = 1,936; 109 incident cases) |
| [16.1] | Attrition analysis (loss to follow-up at age 22) |
| [17] | ICC clustering checks (state- and census tract-level) |
| [18] | Main models: incident asthma (logistic, Models 1–6; primary and sex-interaction specifications) |
| [19] | Secondary analytic datasets (incident Dx and any Dx subsamples for morbidity outcomes) |
| [20] | Hospitalization models (logistic; incident Dx n = 103; any Dx n = 527) |
| [21] | Functional limitation models (OLS; incident Dx n = 106; any Dx n = 533) |
| [22] | Supplemental: spatial sensitivity analyses (1,000m and 500m buffers) |
| [23] | Supplemental: null findings at age 15 (falsification tests) |
| [24] | Supplemental: built environment sensitivity (observer-rated subsample; n = 512) |
| [25] | Supplemental: additional sensitivity checks (prenatal smoking, BMI, caregiver asthma) |
| [26] | Supplemental: propensity score matched analysis (full matching via MatchIt; n = 1,229) |

### Outputs (Sections 27–28)

| Section | Description |
|---------|-------------|
| [27] | Descriptive statistics table (overall and stratified by sex) |
| [28] | Figures 1–4 (forest plots, predicted probability curves; saved as PNG) |

---

## Primary Exposure

**Deadly gun violence at age 15**: Annual count of deadly shooting incidents within 1,600 meters of participants' residential addresses at the age 15 assessment, drawn from GVA data linked to FFCWS geocoded residential locations. Participants with less than one full year of GVA coverage prior to their assessment had counts annualized by dividing by the proportion of the year covered. The count was standardized (mean = 0, SD = 1) within each analytic sample.

---

## Key Outcomes

All outcomes are self-reported at age 22 (youth or primary caregiver proxy):

- **Incident asthma diagnosis**: New physician diagnosis between ages 15 and 22 among participants without a prior diagnosis at age 15 (primary outcome; n = 1,936)
- **Asthma-related hospitalization**: Any hospitalization in the past 12 months, among incident cases (n = 103)
- **Functional limitation**: 4-point ordinal scale ("How much does asthma limit your normal daily activities?"), among incident cases (n = 106)

---

## Statistical Approach

- Logistic regression (incident asthma, hospitalization) and OLS regression (functional limitation)
- Sequential covariate adjustment across 5 models; primary specification is Model 5 (fully adjusted)
- Model 6 tests sex-by-exposure interaction
- All models use HC1 heteroskedasticity-robust standard errors via the `sandwich` package
- Supplemental analyses: spatial buffers (500m, 1,000m), falsification tests, built environment adjustment, additional covariate sensitivity, propensity score matching (MatchIt, full matching, ATT estimand)

---

## R Dependencies

R version 4.0.3. Install all required packages with:

```r
install.packages(c(
  "lubridate", "lmtest", "sandwich", "sjPlot", "dplyr", "corrplot",
  "visdat", "purrr", "broom", "forcats", "stringr", "scales", "cowplot",
  "marginaleffects", "ggplot2", "patchwork", "ggeffects", "haven",
  "lme4", "performance", "tidyr", "DT", "MatchIt", "cobalt",
  "WeightIt", "tibble", "flextable"
))
```

---

## File Paths

The script uses absolute file paths specific to the original author's system (OneDrive/UCF). Before running, update all paths in **Section [2] (Data Import)** and **Section [28] (Figure export)** to point to your local data directory and output folder.

---

## Citation

Martz, C.D., Navas Nazario, A., & Gaydosh, L. Adolescent exposure to deadly neighborhood gun violence and asthma burden in early adulthood. (under review).

---

## Contact

Connor D. Martz, PhD  
Assistant Professor, Department of Population Health Sciences  
UCF College of Medicine  
connor.martz@ucf.edu
