# Adolescent Exposure to Deadly Neighborhood Gun Violence and Asthma Burden in Early Adulthood

**Martz, C.D., Navas Nazario, A.A., & Gaydosh, L.**  
*Journal of Urban Health*

---

## Overview

This repository contains the R analysis code for a prospective study examining whether residential exposure to deadly neighborhood gun violence during adolescence is associated with incident asthma diagnosis, hospitalization, and functional limitation in early adulthood, and whether these associations differ by sex. Analyses use data from the Future of Families and Child Wellbeing Study (FFCWS), linked to geocoded gun violence incidence data from the Gun Violence Archive (GVA) and county violent crime data from the Uniform Crime Reports (UCR).

---

## Data Access

**This repository contains code only. No data are included.**

The FFCWS data used in this study are restricted-use and must be obtained directly from the data distributor:

- **FFCWS core and Wave 7 data**: Available from [Princeton's Center for Research on Child Wellbeing](https://ffcws.princeton.edu/data-and-documentation/data)
- **FFCWS–GVA linked data** (`ff_gva_15y`): Available as a restricted-use appendage; see [documentation](https://ffcws.princeton.edu/sites/g/files/toruqf4356/files/documents/ff_gva_15y_res1_20190603.pdf)
- **Census, geographic identifier, and UCR data**: Restricted-use appendages distributed with the FFCWS restricted-use package

The script reads six Stata files: `FF_wave7_2024v2.dta`, `FF_allwaves_2020v2.dta`, `ffgeo7_all_pub.dta`, `contractcity7pub_121824.dta`, `ff_gva_15y_pub1_20190426.dta`, and `ff_UCR_pub1.dta`.

Secondary analyses were determined to be Not Human Subjects Research by the Institutional Review Board at the University of Central Florida.

---

## Script Structure

All analyses are contained in a single script: `asthma_gv_JUH_final.R`

### Data Cleaning (Sections 1–13)

| Section | Description |
|---------|-------------|
| [1] | Load libraries |
| [2] | Data import (set `data_dir` and `out_dir` here) |
| [3] | Merge FFCWS core, Wave 7, geographic, census, GVA, and UCR files |
| [4] | Date and exposure window variables (GVA coverage flags) |
| [5] | Gun violence exposure variables (binary, count, categorical, annualized rate) |
| [6] | Distance-banded counts (100m intervals, 0–1,600m) |
| [7] | Distance-weighted exposure (quadratic decay) |
| [8] | Time-windowed exposure counts (30-day increments) |
| [9] | Individual-level covariates (race/ethnicity, SES, smoking/vaping, BMI) |
| [10] | Neighborhood and household covariates (tract poverty, violent crime, household composition, traffic density) |
| [11] | Housing quality and built environment index (12 interviewer-rated items, Wave 6) |
| [12] | Asthma outcome variables (incident diagnosis, hospitalization, functional limitation) |
| [13] | Standardize continuous predictors across full dataset |

Sections [6]–[8] and [13] construct variables retained from earlier versions of the analysis that are not used in the reported models.

### Analysis (Sections 14–26)

| Section | Description |
|---------|-------------|
| [14] | Helper functions (robust CIs, sex-specific ORs, table construction) |
| [15] | Create analysis-ready `df_reduced` with factor coding |
| [16] | Analytic sample preparation (`df_main_nomiss`) |
| [16.1] | Attrition analysis (loss to follow-up at age 22) |
| [16.2] | Sample flowchart counts (Figure S1) |
| [17] | ICC clustering checks (state and census tract) |
| [18] | Main models: incident asthma (logistic, Models 1–6; Table S1) |
| [19] | Secondary analytic datasets (incident diagnosis and any diagnosis samples) |
| [20] | Hospitalization models (Firth penalized logistic regression; Tables S2, S4) |
| [21] | Functional limitation models (OLS; Tables S3, S5) |
| [22] | Spatial sensitivity analyses (1,600m, 1,000m, 500m; Figure S2) |
| [23] | Falsification tests (incident asthma ages 9–15; prevalent asthma at age 15) |
| [24] | Built environment sensitivity (home-visit subsample) |
| [25] | Additional confounders (prenatal smoking, BMI, parent asthma) |
| [26] | Propensity score analyses: full matching (ATT) and subclassification (Figures S3–S4) |

### Outputs (Sections 27–28)

| Section | Description |
|---------|-------------|
| [27] | Descriptive statistics (Table 1; overall and by sex) |
| [28] | Figures 1–3 (forest plot, predicted probabilities, secondary outcomes) |

### Additional Sensitivity Analyses (Sections 29–31; Table S6)

| Section | Description |
|---------|-------------|
| [29] | Exposure data quality: GVA–UCR correlation; excluding Michigan and California; restricting to ≥12 months of GVA data |
| [30] | Exposure distribution: categorical exposure with linear trend; log-transformed exposure; excluding values above the 99th and 95th percentiles |
| [31] | Youth-reported neighborhood traffic density (adjustment and interaction) |

Supplemental tables (S1–S5) are saved as .docx files and figures as .png files to `out_dir`.

---

## Primary Exposure

**Deadly gun violence at age 15**: Count of deadly shooting incidents within 1,600 meters of participants' residential addresses in the 12 months before the age 15 assessment, drawn from GVA data linked to FFCWS geocoded residential locations. Participants with less than one full year of GVA coverage before their assessment had counts annualized by dividing by the proportion of the year covered. The annualized rate was standardized (mean = 0, SD = 1) within each analytic sample; spatial buffer and trimmed-sample models use the unstandardized rate (per additional incident per year).

---

## Key Outcomes

All outcomes are self-reported at age 22 (youth report, with primary caregiver proxy report when missing):

- **Incident asthma diagnosis**: New physician diagnosis between ages 15 and 22 among participants without a prior diagnosis at age 15 (primary outcome)
- **Asthma-related hospitalization**: Any hospitalization in the past 12 months, among incident cases and, as a sensitivity analysis, among all participants with asthma by age 22
- **Functional limitation**: 4-point ordinal scale ("How much does asthma limit your normal daily activities?"), among incident cases and among all participants with asthma by age 22

---

## Statistical Approach

- Logistic regression (incident asthma), Firth penalized logistic regression (hospitalization), and OLS regression (functional limitation)
- Sequential covariate adjustment across five models; Model 5 (fully adjusted) is the primary specification
- Model 6 adds an exposure × sex interaction to Model 5
- Heteroskedasticity-robust standard errors via the `sandwich` package (HC type noted in each section); Firth models use profile penalized-likelihood CIs; the full-matching model uses cluster-robust SEs by matched subclass
- Sensitivity analyses: spatial buffers, falsification tests, attrition, GVA data quality, exposure distribution, built environment, prenatal smoking, BMI, parent asthma, traffic density, and propensity score methods (MatchIt full matching and subclassification, ATT estimand)

---

## R Dependencies

Analyses were conducted in R version 4.0.3. Install required packages with:

```r
install.packages(c(
  "lubridate", "haven", "dplyr", "tidyr", "tibble", "lmtest", "sandwich",
  "lme4", "performance", "logistf", "MatchIt", "optmatch", "cobalt",
  "ggplot2", "ggeffects", "scales", "patchwork", "flextable", "Hmisc"
))
```

`optmatch` is required by `MatchIt` for full matching.

---

## File Paths

Before running, set `data_dir` (folder containing the six FFCWS .dta files) and `out_dir` (folder for tables and figures) in **Section [2]**.

---

## Version History

- **Revision (accepted version)**: Hospitalization models re-estimated with Firth penalized logistic regression; added sample flowchart, propensity score subclassification, GVA data quality checks, exposure distribution checks, and traffic density adjustment; expanded Table 1; added code for Tables S1–S5 and Figures S2 and S4; combined secondary outcome figures into Figure 3; replaced local file paths with `data_dir`/`out_dir`.
- **Initial submission**: Original analysis code.

---

## Citation

Martz, C.D., Navas Nazario, A.A., & Gaydosh, L. Adolescent exposure to deadly neighborhood gun violence and asthma burden in early adulthood. *Journal of Urban Health* (in press).

---

## Contact

Connor D. Martz, PhD  
Assistant Professor, Department of Population Health Sciences  
UCF College of Medicine  
connor.martz@ucf.edu
