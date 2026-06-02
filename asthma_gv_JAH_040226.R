# =============================================================================
# ADOLESCENT EXPOSURE TO DEADLY GUN VIOLENCE & ASTHMA BURDEN IN EARLY ADULTHOOD
# Martz, Navas Nazario & Gaydosh | Journal of Adolescent Health Submission
# =============================================================================
#
# SCRIPT STRUCTURE:
#
#   DATA CLEANING (Sections 1–13)
#     [1]  Load Libraries
#     [2]  Data Import
#     [3]  Merge Data
#     [4]  Date & Exposure Window Variables
#     [5]  Gun Violence Exposure Variables
#     [6]  Distance-Banded Counts
#     [7]  Distance-Weighted Exposure
#     [8]  Time-Windowed Exposure Counts
#     [9]  Individual-Level Covariates
#     [10] Neighborhood & Household Covariates
#     [11] Housing Quality & Asthma Environmental Risk
#     [12] Asthma Outcomes
#     [13] Standardize Key Continuous Predictors (full df)
#
#   ANALYSIS (Sections 14–26)
#     [14] Helper Functions (robust CI extraction)
#     [15] Create df_reduced
#     [16] Analytic Sample Preparation (df_main_nomiss; n=1,936)
#     [16.1] Attrition Analysis
#     [17] Check Clustering (ICC)
#     [18] Main Models: Incident Asthma (n=1,936)
#     [19] Secondary Analytic Datasets
#     [20] Hospitalization Models
#     [21] Severity Models
#     [22] Supplemental: Spatial Sensitivity
#     [23] Supplemental: Null Findings at Age 15
#     [24] Supplemental: Built Environment Sensitivity
#     [25] Supplemental: Additional Sensitivity Checks
#     [26] Supplemental: Propensity Score Matched Analysis
#
#   OUTPUTS (Sections 27–28)
#     [27] Descriptive Statistics Table
#     [28] Figures
#
#   EXPLORATORY (Section 29 — not reported in manuscript)
#     [29] Policy Analysis (exploratory; retained for reference)
#
# =============================================================================


# =============================================================================
# [1] LOAD LIBRARIES ----
# =============================================================================

library(lubridate)
library(lmtest)
library(sandwich)
library(sjPlot)
library(dplyr)
library(corrplot)
library(visdat)
library(purrr)
library(broom)
library(forcats)
library(stringr)
library(scales)
library(cowplot)
library(marginaleffects)
library(ggplot2)
library(patchwork)
library(ggeffects)
library(haven)
library(lme4)
library(performance)
library(tidyr)
library(DT)
library(MatchIt)
library(cobalt)
library(WeightIt)
library(tibble)
library(flextable)
library(tibble)


# =============================================================================
# [2] DATA IMPORT ----
# =============================================================================

y22     <- read_dta(
  "/Users/co453621/Library/CloudStorage/OneDrive-UniversityofCentralFlorida/2. MANUSCRIPTS/MARTZ GUN VIOLENCE ASTHMA/data/FF_wave7_2024v2.dta"
)
core    <- read_dta(
  "/Users/co453621/Library/CloudStorage/OneDrive-UniversityofCentralFlorida/2. MANUSCRIPTS/MARTZ GUN VIOLENCE ASTHMA/data/FF_allwaves_2020v2.dta"
)
census00 <- read_csv("/Users/co453621/Desktop/census.csv")
geoid   <- read_csv("/Users/co453621/Desktop/geoid.csv")
FF_gva  <- read_csv("/Users/co453621/Desktop/gva.csv")
FF_ucr  <- read_csv("/Users/co453621/Desktop/ucr.csv")


# =============================================================================
# [3] MERGE DATA ----
# =============================================================================

df <- core %>%
  inner_join(y22, by = "idnum") %>%
  left_join(geoid, by = "idnum") %>%
  left_join(census00, by = "idnum") %>%
  left_join(FF_gva, by = "idnum") %>%
  left_join(FF_ucr, by = "idnum")


# =============================================================================
# [4] DATE & EXPOSURE WINDOW VARIABLES ----
# =============================================================================

df <- df %>%
  mutate(
    # Date of Year 15 assessment (estimated from birth + age in months)
    date_Y15 = case_when(
      !is.na(ck6intyr) & !is.na(ck6intmon) & ck6intyr > 0 & ck6intmon > 0
      ~ as.Date(paste(ck6intyr, ck6intmon, "15", sep = "-")),
      TRUE ~ NA_Date_
    ),
    
    # GVA data begins January 1, 2014; compute years of available exposure prior to Y15
    gva_start_date = as.Date("2014-01-01"),
    years_dgv_exposure_prior_to_Y15 =
      as.numeric(difftime(date_Y15, gva_start_date, units = "days")) / 365,
    
    # Flag exposure coverage:
    #   1 = at least one full year of GVA data available before Y15 assessment
    #   0 = partial year only (assessment within first year of GVA coverage)
    #   NA = assessment predates GVA launch (negative window; no valid data)
    atleast_1year_exposure = case_when(
      years_dgv_exposure_prior_to_Y15 <= 0 ~ NA_real_,  # predates GVA; invalid
      years_dgv_exposure_prior_to_Y15 <  1 ~ 0,          # partial year
      years_dgv_exposure_prior_to_Y15 >= 1 ~ 1           # full year or more
    )
  )

# Verify exposure window distribution
summary(df$years_dgv_exposure_prior_to_Y15)
table(df$atleast_1year_exposure, useNA = "ifany")
# Expected: most participants have >= 1 year; small number with partial coverage;
# participants with negative values (pre-2014 assessments) should be NA


# =============================================================================
# [5] GUN VIOLENCE EXPOSURE VARIABLES ----
# =============================================================================

df <- df %>%
  mutate(
    
    # --- 5.1 Binary: Any deadly gun violence (DGV) within buffer ---
    gvh1600_any = case_when(
      rg6gva_totl_1600m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_1600m_home_365d_k6 >  0 ~ 1,
      TRUE ~ NA_real_
    ),
    gvh1000_any = case_when(
      rg6gva_totl_1000m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_1000m_home_365d_k6 >  0 ~ 1,
      TRUE ~ NA_real_
    ),
    gvh500_any = case_when(
      rg6gva_totl_500m_home_365d_k6 == 0 ~ 0,
      rg6gva_totl_500m_home_365d_k6 >  0 ~ 1,
      TRUE ~ NA_real_
    ),
    
    # --- 5.2 Continuous count within primary buffer (1600m) ---
    gvh1600count = ifelse(rg6gva_totl_1600m_home_365d_k6 >= 0,
                          rg6gva_totl_1600m_home_365d_k6, NA_real_),
    gvh1600count_log = log(gvh1600count + 1),
    
    # Sensitivity buffers
    gvh1000count = ifelse(rg6gva_totl_1000m_home_365d_k6 >= 0,
                          rg6gva_totl_1000m_home_365d_k6, NA_real_),
    gvh500count  = ifelse(rg6gva_totl_500m_home_365d_k6  >= 0,
                          rg6gva_totl_500m_home_365d_k6,  NA_real_),
    
    # --- 5.3 Categorical count (1600m): 0=unexposed, 1=low, 2=moderate, 3=high ---
    gvh1600cat = case_when(
      gvh1600count == 0                    ~ 0,
      gvh1600count == 1                    ~ 1,
      gvh1600count > 1 & gvh1600count <= 5 ~ 2,
      gvh1600count > 5                     ~ 3,
      TRUE ~ NA_real_
    ),
    gvh1600count_cat = case_when(
      gvh1600count == 0        ~ 0,
      gvh1600count %in% 1:5   ~ 1,
      gvh1600count %in% 6:10  ~ 2,
      gvh1600count > 10        ~ 3,
      TRUE ~ NA_real_
    ),
    
    # --- 5.4 Annualized rate (incidents per year) ----
    gvh1600_rate = case_when(
      is.na(atleast_1year_exposure)  ~ NA_real_,
      atleast_1year_exposure == 1    ~ gvh1600count,
      atleast_1year_exposure == 0    ~ gvh1600count / years_dgv_exposure_prior_to_Y15
    ),
    gvh1000_rate = case_when(
      is.na(atleast_1year_exposure)  ~ NA_real_,
      atleast_1year_exposure == 1    ~ gvh1000count,
      atleast_1year_exposure == 0    ~ gvh1000count / years_dgv_exposure_prior_to_Y15
    ),
    gvh500_rate = case_when(
      is.na(atleast_1year_exposure)  ~ NA_real_,
      atleast_1year_exposure == 1    ~ gvh500count,
      atleast_1year_exposure == 0    ~ gvh500count  / years_dgv_exposure_prior_to_Y15
    )
  )


# =============================================================================
# [6] DISTANCE-BANDED COUNTS (100m intervals + irregular intervals, 0–1600m) ----
# Used for distance-decay weighting and spatial sensitivity analyses
# =============================================================================

# Helper to compute incidents in each distance band (outer - inner buffer)
make_band <- function(outer, inner) {
  case_when(
    outer == 0 & inner == 0 ~ 0,
    outer >= 0 & inner >= 0 ~ outer - inner,
    TRUE ~ NA_real_
  )
}

df <- df %>%
  mutate(
    # Primary 100m bands (0–1600m)
    gvh0_100count     = ifelse(rg6gva_totl_100m_home_365d_k6 >= 0,
                               rg6gva_totl_100m_home_365d_k6, NA_real_),
    gvh100_200count   = make_band(rg6gva_totl_200m_home_365d_k6,  rg6gva_totl_100m_home_365d_k6),
    gvh200_300count   = make_band(rg6gva_totl_300m_home_365d_k6,  rg6gva_totl_200m_home_365d_k6),
    gvh300_400count   = make_band(rg6gva_totl_400m_home_365d_k6,  rg6gva_totl_300m_home_365d_k6),
    gvh400_500count   = make_band(rg6gva_totl_500m_home_365d_k6,  rg6gva_totl_400m_home_365d_k6),
    gvh500_600count   = make_band(rg6gva_totl_600m_home_365d_k6,  rg6gva_totl_500m_home_365d_k6),
    gvh600_700count   = make_band(rg6gva_totl_700m_home_365d_k6,  rg6gva_totl_600m_home_365d_k6),
    gvh700_800count   = make_band(rg6gva_totl_800m_home_365d_k6,  rg6gva_totl_700m_home_365d_k6),
    gvh800_900count   = make_band(rg6gva_totl_900m_home_365d_k6,  rg6gva_totl_800m_home_365d_k6),
    gvh900_1000count  = make_band(rg6gva_totl_1000m_home_365d_k6, rg6gva_totl_900m_home_365d_k6),
    gvh1000_1100count = make_band(rg6gva_totl_1100m_home_365d_k6, rg6gva_totl_1000m_home_365d_k6),
    gvh1100_1200count = make_band(rg6gva_totl_1200m_home_365d_k6, rg6gva_totl_1100m_home_365d_k6),
    gvh1200_1300count = make_band(rg6gva_totl_1300m_home_365d_k6, rg6gva_totl_1200m_home_365d_k6),
    gvh1300_1400count = make_band(rg6gva_totl_1400m_home_365d_k6, rg6gva_totl_1300m_home_365d_k6),
    gvh1400_1500count = make_band(rg6gva_totl_1500m_home_365d_k6, rg6gva_totl_1400m_home_365d_k6),
    gvh1500_1600count = make_band(rg6gva_totl_1600m_home_365d_k6, rg6gva_totl_1500m_home_365d_k6),
    
    # Irregular intervals (retained for sensitivity analyses)
    gvh200_250count   = make_band(rg6gva_totl_250m_home_365d_k6,  rg6gva_totl_200m_home_365d_k6),
    gvh250_300count   = make_band(rg6gva_totl_300m_home_365d_k6,  rg6gva_totl_250m_home_365d_k6),
    gvh700_750count   = make_band(rg6gva_totl_750m_home_365d_k6,  rg6gva_totl_700m_home_365d_k6),
    gvh750_800count   = make_band(rg6gva_totl_800m_home_365d_k6,  rg6gva_totl_750m_home_365d_k6),
    gvh1200_1250count = make_band(rg6gva_totl_1250m_home_365d_k6, rg6gva_totl_1200m_home_365d_k6)
  )


# =============================================================================
# [7] DISTANCE-WEIGHTED EXPOSURE (quadratic decay) ----
# Bands closer to home receive greater weight; decay from ring 16 inward
# =============================================================================

df <- df %>%
  mutate(
    dist_weight_count =
      (gvh0_100count     * (17 - 1)^2) +
      (gvh100_200count   * (17 - 2)^2) +
      (gvh200_300count   * (17 - 3)^2) +
      (gvh300_400count   * (17 - 4)^2) +
      (gvh400_500count   * (17 - 5)^2) +
      (gvh500_600count   * (17 - 6)^2) +
      (gvh600_700count   * (17 - 7)^2) +
      (gvh700_800count   * (17 - 8)^2) +
      (gvh800_900count   * (17 - 9)^2) +
      (gvh900_1000count  * (17 - 10)^2) +
      (gvh1000_1100count * (17 - 11)^2) +
      (gvh1100_1200count * (17 - 12)^2) +
      (gvh1200_1300count * (17 - 13)^2) +
      (gvh1300_1400count * (17 - 14)^2) +
      (gvh1400_1500count * (17 - 15)^2) +
      (gvh1500_1600count * (17 - 16)^2),
    dist_weight_count_log = log(dist_weight_count + 1)
  )


# =============================================================================
# [8] TIME-WINDOWED EXPOSURE COUNTS (30-day increments, 1600m) ----
# For temporal sensitivity analyses
# =============================================================================

df <- df %>%
  mutate(
    # Raw counts at each time window
    gvh1600count_365 = ifelse(rg6gva_totl_1600m_home_365d_k6 >= 0, rg6gva_totl_1600m_home_365d_k6, NA_real_),
    gvh1600count_360 = ifelse(rg6gva_totl_1600m_home_360d_k6 >= 0, rg6gva_totl_1600m_home_360d_k6, NA_real_),
    gvh1600count_330 = ifelse(rg6gva_totl_1600m_home_330d_k6 >= 0, rg6gva_totl_1600m_home_330d_k6, NA_real_),
    gvh1600count_300 = ifelse(rg6gva_totl_1600m_home_300d_k6 >= 0, rg6gva_totl_1600m_home_300d_k6, NA_real_),
    gvh1600count_270 = ifelse(rg6gva_totl_1600m_home_270d_k6 >= 0, rg6gva_totl_1600m_home_270d_k6, NA_real_),
    gvh1600count_240 = ifelse(rg6gva_totl_1600m_home_240d_k6 >= 0, rg6gva_totl_1600m_home_240d_k6, NA_real_),
    gvh1600count_210 = ifelse(rg6gva_totl_1600m_home_210d_k6 >= 0, rg6gva_totl_1600m_home_210d_k6, NA_real_),
    gvh1600count_180 = ifelse(rg6gva_totl_1600m_home_180d_k6 >= 0, rg6gva_totl_1600m_home_180d_k6, NA_real_),
    gvh1600count_150 = ifelse(rg6gva_totl_1600m_home_150d_k6 >= 0, rg6gva_totl_1600m_home_150d_k6, NA_real_),
    gvh1600count_120 = ifelse(rg6gva_totl_1600m_home_120d_k6 >= 0, rg6gva_totl_1600m_home_120d_k6, NA_real_),
    gvh1600count_90  = ifelse(rg6gva_totl_1600m_home_90d_k6  >= 0, rg6gva_totl_1600m_home_90d_k6,  NA_real_),
    gvh1600count_60  = ifelse(rg6gva_totl_1600m_home_60d_k6  >= 0, rg6gva_totl_1600m_home_60d_k6,  NA_real_),
    gvh1600count_30  = ifelse(rg6gva_totl_1600m_home_30d_k6  >= 0, rg6gva_totl_1600m_home_30d_k6,  NA_real_),
    gvh1600count_14  = ifelse(rg6gva_totl_1600m_home_14d_k6  >= 0, rg6gva_totl_1600m_home_14d_k6,  NA_real_),
    gvh1600count_7   = ifelse(rg6gva_totl_1600m_home_7d_k6   >= 0, rg6gva_totl_1600m_home_7d_k6,   NA_real_),
    
    # Incremental time bands (events within each 30-day window)
    gvh1600count_365_360 = gvh1600count_365 - gvh1600count_360,
    gvh1600count_360_330 = gvh1600count_360 - gvh1600count_330,
    gvh1600count_330_300 = gvh1600count_330 - gvh1600count_300,
    gvh1600count_300_270 = gvh1600count_300 - gvh1600count_270,
    gvh1600count_270_240 = gvh1600count_270 - gvh1600count_240,
    gvh1600count_240_210 = gvh1600count_240 - gvh1600count_210,
    gvh1600count_210_180 = gvh1600count_210 - gvh1600count_180,
    gvh1600count_180_150 = gvh1600count_180 - gvh1600count_150,
    gvh1600count_150_120 = gvh1600count_150 - gvh1600count_120,
    gvh1600count_120_90  = gvh1600count_120 - gvh1600count_90,
    gvh1600count_90_60   = gvh1600count_90  - gvh1600count_60,
    gvh1600count_60_30   = gvh1600count_60  - gvh1600count_30,
    gvh1600count_30_14   = gvh1600count_30  - gvh1600count_14,
    gvh1600count_14_7    = gvh1600count_14  - gvh1600count_7,
    gvh1600count_7_0     = gvh1600count_7
  )


# =============================================================================
# [9] INDIVIDUAL-LEVEL COVARIATES ----
# =============================================================================

df <- df %>%
  mutate(
    
    # --- 9.1 Race/ethnicity (time-invariant; prioritize youth self-report at Y22) ---
    race = case_when(
      ck7ethrace == 1 ~ 1,                 # Non-Hispanic White
      ck7ethrace == 2 ~ 2,                 # Non-Hispanic Black
      ck7ethrace == 3 ~ 3,                 # Hispanic
      ck7ethrace %in% c(4, 5) ~ 4,         # Multiracial / Other Non-Hispanic
      ck7ethrace %in% c(-1,-2,-3) & cm1ethrace == 2 & (cf1ethrace == 2 | (m1i4 == 2 & m1i4a != 1)) ~ 2,
      ck7ethrace %in% c(-1,-2,-3) & cm1ethrace == 1 & (cf1ethrace == 1 | (m1i4 == 1 & m1i4a != 1)) ~ 1,
      ck7ethrace %in% c(-1,-2,-3) & cm1ethrace == 3 & (cf1ethrace == 3 | (m1i4 == 5 & m1i4a == 1)) ~ 3,
      ck7ethrace %in% c(-1,-2,-3) & (cm1ethrace %in% c(4,5) | cf1ethrace %in% c(4,5)) ~ 4,
      ck7ethrace %in% c(-1,-2,-3) & cm1ethrace == 1 & cf1ethrace %in% c(2,3) ~ 4,
      ck7ethrace %in% c(-1,-2,-3) & cf1ethrace == 1 & cm1ethrace %in% c(2,3) ~ 4,
      ck7ethrace %in% c(-1,-2,-3) & cm1ethrace == 2 & cf1ethrace %in% c(1,3) ~ 4,
      ck7ethrace %in% c(-1,-2,-3) & cf1ethrace == 2 & cm1ethrace %in% c(1,3) ~ 4,
      ck7ethrace %in% c(-1,-2,-3) & cm1ethrace == 3 & cf1ethrace %in% c(1,2) ~ 4,
      ck7ethrace %in% c(-1,-2,-3) & cf1ethrace == 3 & cm1ethrace %in% c(1,2) ~ 4,
      ck7ethrace == -9 & cm1ethrace == 1 & cf1ethrace == 1 ~ 1,
      ck7ethrace == -9 & cm1ethrace == 2 & cf1ethrace == 2 ~ 2,
      ck7ethrace == -9 & cm1ethrace == 3 & cf1ethrace == 3 ~ 3,
      ck7ethrace == -9 & cm1ethrace == 1 & cf1ethrace %in% c(2,3,4,5) ~ 4,
      ck7ethrace == -9 & cf1ethrace == 1 & cm1ethrace %in% c(2,3,4,5) ~ 4,
      TRUE ~ NA_real_
    ),
    # Race/ethnicity indicator variables
    black  = case_when(race == 2 ~ 1, race %in% c(1,3,4) ~ 0, TRUE ~ NA_real_),
    white  = case_when(race == 1 ~ 1, race %in% c(2,3,4) ~ 0, TRUE ~ NA_real_),
    hisp   = case_when(race == 3 ~ 1, race %in% c(1,2,4) ~ 0, TRUE ~ NA_real_),
    other  = case_when(race == 4 ~ 1, race %in% c(1,2,3) ~ 0, TRUE ~ NA_real_),
    
    # --- 9.2 Sex (time-invariant) ---
    female = case_when(cm1bsex == 1 ~ 0, cm1bsex == 2 ~ 1, TRUE ~ NA_real_),
    
    # --- 9.3 Maternal smoking at birth ---
    momsmk = ifelse(is.na(m1g4), NA, ifelse(m1g4 == 4, 0, 1)),
    
    # --- 9.4 Age at each wave ---
    age_t5 = case_when(ch5agem > 0 ~ ch5agem / 12, TRUE ~ NA_real_),
    age_t6 = case_when(ck6yagem > 0 ~ ck6yagem / 12, TRUE ~ NA_real_),
    age_t7 = case_when(ck7yagem > 0 ~ ck7yagem / 12, TRUE ~ NA_real_),
    
    # --- 9.5 Caregiver education ---
    edu_t5 = case_when(cm5edu > 0 ~ cm5edu, TRUE ~ NA_real_),
    edu_t6 = case_when(cp6edu > 0 ~ cp6edu, TRUE ~ NA_real_),
    edu_t7 = case_when(ck7edu > 0 ~ ck7edu, TRUE ~ NA_real_),
    
    # --- 9.6 Income-to-poverty ratio ---
    pov_t5 = case_when(cm5povco %in% c(-1,-3,-9) ~ NA_real_, TRUE ~ cm5povco),
    pov_t6 = case_when(cp6povco %in% c(-1,-3,-9) ~ NA_real_, TRUE ~ cp6povco),
    pov_t7 = case_when(ck7povco %in% c(-1,-3,-9) ~ NA_real_, TRUE ~ ck7povco),
    
    # --- 9.7 BMI z-score ---
    bmiz_t5 = case_when(is.na(ch5bmiz) ~ NA_real_, TRUE ~ ch5bmiz),
    bmiz_t6 = case_when(is.na(ck6bmiz) ~ NA_real_, TRUE ~ ck6bmiz),
    bmiz_t7 = case_when(is.na(ck7bmiz) ~ NA_real_, TRUE ~ ck7bmiz),
    
    # --- 9.8 Smoking and vaping ---
    smoking_t6 = case_when(      # Any cigarette use in past month (Y15)
      k6d42 == 1 | k6d42 == -6 ~ 0,
      k6d42 > 1               ~ 1,
      TRUE ~ NA_real_
    ),
    smoking_t7 = case_when(      # Any cigarette use in past 30 days (Y22)
      k7j4 == 0 ~ 0,
      k7j4 > 0  ~ 1,
      TRUE ~ NA_real_
    ),
    vaping_t7 = case_when(       # Vaped nicotine in past 12 months (Y22)
      k7j2a == 2 | k7j2a == -6 ~ 0,
      k7j2a == 1               ~ 1,
      TRUE ~ NA_real_
    )
  )


# =============================================================================
# [10] NEIGHBORHOOD & HOUSEHOLD COVARIATES ----
# =============================================================================

df <- df %>%
  mutate(
    
    # --- 10.1 State of residence ---
    state_t4 = case_when(m4stfips %in% c(-3,-7,-9) ~ NA_real_, TRUE ~ m4stfips),
    state_t5 = case_when(m5stfips %in% c(-3,-7,-9) ~ NA_real_, TRUE ~ m5stfips),
    state_t6 = case_when(p6state_n %in% c(-3,-9)   ~ NA_real_, TRUE ~ p6state_n),
    
    # --- 10.2 Census tract poverty (%) ---
    trct_pct_pov_t5 = case_when(
      tm5pfbpl_cen00 %in% c(-9,-7,-3) ~ NA_real_,
      TRUE ~ tm5pfbpl_cen00 * 100
    ),
    trct_pct_pov_t6 = case_when(
      tp6pfbpl_acs15 %in% c(-9,-7,-3) ~ NA_real_,
      TRUE ~ tp6pfbpl_acs15 * 100
    ),
    
    # --- 10.3 Residential mobility ---
    moved_tracts_9_15 = case_when(
      cp6_moved_9_15 == 1 ~ 1,
      cp6_moved_9_15 == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # --- 10.4 County violent crime rate (UCR) ---
    vc_rate0  = case_when(rg1ucr_mviort >= 0 ~ as.numeric(rg1ucr_mviort), TRUE ~ NA_real_),
    vc_rate1  = case_when(rg2ucr_mviort >= 0 ~ as.numeric(rg2ucr_mviort), TRUE ~ NA_real_),
    vc_rate3  = case_when(rg3ucr_mviort >= 0 ~ as.numeric(rg3ucr_mviort), TRUE ~ NA_real_),
    vc_rate5  = case_when(rg4ucr_mviort >= 0 ~ as.numeric(rg4ucr_mviort), TRUE ~ NA_real_),
    vc_rate9  = case_when(rg5ucr_mviort >= 0 ~ as.numeric(rg5ucr_mviort), TRUE ~ NA_real_),
    vc_rate15 = case_when(rg6ucr_pviort >= 0 ~ as.numeric(rg6ucr_pviort), TRUE ~ NA_real_),
    
    # --- 10.5 Household composition ---
    household_smokers_t6 = case_when(
      p6h78 == 0 ~ 0,
      p6h78 >  0 ~ 1,
      TRUE ~ NA_real_
    ),
    household_size_t6_raw = ifelse(
      cp6adult %in% c(-9,-3) | cp6kids %in% c(-9,-3), NA_real_,
      cp6adult + cp6kids
    ),
    household_size_t7_raw = ifelse(ck7hhsize %in% c(-9,-3,-1), NA_real_, ck7hhsize),
    household_size_t6 = case_when(household_size_t6_raw > 10 ~ 11, TRUE ~ household_size_t6_raw),
    household_size_t7 = case_when(household_size_t7_raw > 10 ~ 11, TRUE ~ household_size_t7_raw),
    
    # --- 10.6 Tract ID (Wave 6 / Year 15) ---
    trct_id6 = case_when(tp6tract_cen00 >= 0 ~ tp6tract_cen00, TRUE ~ NA)
  )


# =============================================================================
# [11] HOUSING QUALITY & ASTHMA ENVIRONMENTAL RISK ----
# =============================================================================
df <- df %>%
  mutate(
    # --- 11.1 Individual housing condition indicators (Y15 observer-rated) ---
    broken_windows        = case_when(o6c1     == 1 ~ 1, o6c1     == 2 ~ 0, o6c1     %in% c(-3,-9,-1) ~ NA_real_),
    exposed_wires         = case_when(o6c2     == 1 ~ 1, o6c2     == 2 ~ 0, o6c2     %in% c(-3,-9,-1) ~ NA_real_),
    wall_cracks           = case_when(o6c3     == 1 ~ 1, o6c3     == 2 ~ 0, o6c3     %in% c(-3,-9,-1) ~ NA_real_),
    floor_holes           = case_when(o6c4     == 1 ~ 1, o6c4     == 2 ~ 0, o6c4     %in% c(-3,-9,-1) ~ NA_real_),
    peeling_paint         = case_when(o6c5     == 1 ~ 1, o6c5     == 2 ~ 0, o6c5     %in% c(-3,-9,-1) ~ NA_real_),
    darkness              = case_when(o6c6     == 1 ~ 1, o6c6     == 2 ~ 0, o6c6     %in% c(-3,-9,-1) ~ NA_real_),
    crowding              = case_when(o6c7     == 1 ~ 1, o6c7     == 2 ~ 0, o6c7     %in% c(-3,-9,-1) ~ NA_real_),
    clutter               = case_when(o6c8     == 1 ~ 1, o6c8     == 2 ~ 0, o6c8     %in% c(-3,-9,-1) ~ NA_real_),
    dirtiness             = case_when(o6c9     == 1 ~ 1, o6c9     == 2 ~ 0, o6c9     %in% c(-3,-9,-1) ~ NA_real_),
    mice_rats             = case_when(o6c10a_2 == 1 ~ 1, TRUE ~ 0),
    hazardness            = case_when(o6c10    == 1 ~ 1, o6c10    == 2 ~ 0, o6c10    %in% c(-3,-9,-1) ~ NA_real_),
    
    # --- 11.2 Neighborhood and home trash ---
    neighborhood_trash_yn = case_when(
      o6a1 %in% c(-3,-9) ~ NA_real_,
      o6a1 == 1          ~ 0,
      o6a1 > 1           ~ 1
    ),
    home_trash_yn = case_when(
      o6a6f %in% c(-3,-9) ~ NA_real_,
      o6a6f == 1 ~ 1,
      o6a6f == 2 ~ 0
    )
  ) %>%
  mutate(
    # --- 11.3 Asthma-relevant environmental risk index ---
    # Components: peeling paint, dirtiness, mice/rats, neighborhood trash, home trash
    # Excludes household_smokers_t6 (entered separately as primary model covariate)
    # Range: 0-5; no capping needed given maximum is 5
    asthma_env_risk = case_when(
      rowSums(is.na(select(., broken_windows:mice_rats,
                           neighborhood_trash_yn, home_trash_yn))) > 0 ~ NA_real_,
      TRUE ~ rowSums(select(., broken_windows:mice_rats,
                            neighborhood_trash_yn, home_trash_yn))
    )
  )

table(df$asthma_env_risk)

# =============================================================================
# [12] ASTHMA OUTCOMES ----
# =============================================================================

df <- df %>%
  mutate(
    
    # --- 12.1 Incident asthma diagnosis (Y9 to Y15) ---
    newasthmadx_15 = case_when(
      p5h1b %in% c(1, -2, -9)                                       ~ NA_real_,
      p6b2  %in% c(-9, -2)                                          ~ NA_real_,
      p5h1b == 2 & p6b2 == 2                                        ~ 0,
      p5h1b == 2 & p6b2 %in% c(-9,-3,-2,-1) & p6b27_3 == 0         ~ 0,
      p5h1b == 2 & p6b2 == 1                                        ~ 1,
      p5h1b == 2 & p6b2 %in% c(-9,-3,-2,-1) & p6b27_3 == 1         ~ 1,
      TRUE ~ NA_real_
    ),
    
    # --- 12.2 Incident asthma diagnosis (Y15 to Y22) | PRIMARY OUTCOME ---
    newasthmadx_22 = case_when(
      p6b2  %in% c(1, -2, -9)                                       ~ NA_real_,
      k7i4_1 %in% c(-9,-3,-2,-1) & p7i3_1 %in% c(-9,-6,-3,-1)      ~ NA_real_,
      p6b2 == 2 & k7i4_1 == 0                                       ~ 0,
      p6b2 == 2 & k7i4_1 %in% c(-9,-3,-2,-1) & p7i3_1 == 0         ~ 0,
      p6b2 == 2 & k7i4_1 == 1                                       ~ 1,
      p6b2 == 2 & k7i4_1 %in% c(-9,-3,-2,-1) & p7i3_1 == 1         ~ 1,
      TRUE ~ NA_real_
    ),
    
    # --- 12.3 Any asthma diagnosis at Y15 or Y22 (prevalence) ---
    anyasthmadx_15 = case_when(p6b2 == 1 ~ 1, p6b2 == 2 ~ 0, TRUE ~ NA_real_),
    anyasthmadx_22 = case_when(k7i4_1 == 1 ~ 1, k7i4_1 == 0 ~ 0, TRUE ~ NA_real_),
    
    # --- 12.4 Age at asthma diagnosis (Y22 self-report) ---
    age_asthma_dx = case_when(k7i9 %in% c(-9,-6,-3,-1) ~ NA_real_, TRUE ~ k7i9),
    
    # --- 12.5 Asthma medication use at Y22 ---
    asthma_meds = case_when(
      k7i11 %in% c(-9,-6,-3,-1) ~ NA_real_,
      k7i11 == 1 ~ 1,
      k7i11 == 2 ~ 0
    ),
    
    # --- 12.6 Asthma-related hospitalizations at Y22 (ordinal: 0, 1, 2, 3+) ---
    hosp_asthma = case_when(
      k7i10 == 0  ~ 0,
      k7i10 == 1  ~ 1,
      k7i10 == 2  ~ 2,
      k7i10 >  2  ~ 3,
      TRUE ~ NA_real_
    ),
    hosp_asthma_yn = case_when(
      k7i10 == 0 ~ 0,
      k7i10 >  0 ~ 1,
      TRUE ~ NA_real_
    ),
    
    # --- 12.7 Functional limitation due to asthma (reverse-coded; higher = more limited) ---
    limit_asthma = case_when(
      k7i29_1 == 0 ~ 0,
      k7i29_1 == 1 ~ 4,
      k7i29_1 == 2 ~ 3,
      k7i29_1 == 3 ~ 2,
      k7i29_1 == 4 ~ 1,
      TRUE ~ NA_real_
    ),
    
    # --- 12.8 PCG-reported asthma (Y22 and Y15; used as proxy when youth report missing) ---
    pcg_asthma = case_when(
      p7l2a == 1 | p6h4_2 == 1                                   ~ 1,
      p7l2a == 2                                                  ~ 0,
      p7l2a %in% c(-9,-3,-1) & p6h4_2 %in% c(0, -6)             ~ 0,
      p7l2a %in% c(-9,-3,-1) & p6h4_2 %in% c(-9,-3)             ~ NA_real_
    )
  )


# =============================================================================
# [13] STANDARDIZE KEY CONTINUOUS PREDICTORS (full df) ----
# Note: variables are re-standardized within each analytic subsample below
# =============================================================================

df <- df %>%
  mutate(
    gvh1600count_std      = as.numeric(scale(gvh1600count)),
    gvh1000count_std      = as.numeric(scale(gvh1000count)),
    gvh500count_std       = as.numeric(scale(gvh500count)),
    dist_weight_count_std = as.numeric(scale(dist_weight_count)),
    gvh1600rate_std       = as.numeric(scale(gvh1600_rate)),
    gvh1000rate_std       = as.numeric(scale(gvh1000_rate)),
    gvh500rate_std        = as.numeric(scale(gvh500_rate)),
    vc_rate15_std         = as.numeric(scale(vc_rate15)),
    vc_rate9_std          = as.numeric(scale(vc_rate9)),
    pov_t6_std            = as.numeric(scale(pov_t6)),
    pov_t5_std            = as.numeric(scale(pov_t5)),
    trct_pct_pov_t6_std   = as.numeric(scale(trct_pct_pov_t6)),
    trct_pct_pov_t5_std   = as.numeric(scale(trct_pct_pov_t5))
  )


# =============================================================================
# ---- ANALYSIS SCRIPT BEGINS HERE ----
# =============================================================================


# =============================================================================
# [14] HELPER FUNCTIONS FOR ROBUST CI EXTRACTION ----
# Defined at top of analysis section so available to all downstream sections
# =============================================================================

# Single coefficient: HC1 robust CI, with optional transformation (default: exp for OR)
robust_ci_term <- function(fit, term, transform = exp, type = "HC1") {
  V  <- sandwich::vcovHC(fit, type = type)
  b  <- coef(fit)[[term]]
  se <- sqrt(diag(V))[term]
  ci <- b + c(-1, 1) * 1.96 * se
  tibble(estimate = transform(b), l95 = transform(ci[1]), u95 = transform(ci[2]))
}

# Linear combination of coefficients (e.g., main + interaction term for female OR)
robust_ci_lincomb <- function(fit, coefs, transform = exp, type = "HC1") {
  V   <- sandwich::vcovHC(fit, type = type)
  b   <- coef(fit)[names(coefs)]
  cv  <- matrix(unname(coefs), ncol = 1)
  est <- sum(b * coefs)
  se  <- sqrt(drop(t(cv) %*% V[names(coefs), names(coefs), drop = FALSE] %*% cv))
  ci  <- est + c(-1, 1) * 1.96 * se
  tibble(estimate = transform(est), l95 = transform(ci[1]), u95 = transform(ci[2]))
}


# =============================================================================
# [15] CREATE df_reduced ----
# Selects all variables needed across main and supplemental analyses
# Input: df (from data cleaning, Sections 1–13)
# Output: df_reduced (analysis-ready; factor coding applied)
# =============================================================================

df_reduced <- df %>%
  dplyr::select(
    idnum, ck7kint, state_t6, state_t5, state_t4, trct_id6,
    # Gun violence exposure
    gvh1600_any, gvh1600count, gvh1600count_log, gvh1600count_cat, gvh1600cat,
    gvh1000_any, gvh1000count,
    gvh500_any,  gvh500count,
    gvh1600_rate, gvh1000_rate, gvh500_rate,
    # Time-varying individual covariates
    age_t7, age_t6, age_t5,
    edu_t7, edu_t6, edu_t5,
    pov_t7, pov_t6, pov_t5,
    race, female, momsmk,
    smoking_t6, smoking_t7, vaping_t7,
    bmiz_t5, bmiz_t6, bmiz_t7,
    # Neighborhood covariates
    trct_pct_pov_t5, trct_pct_pov_t6,
    vc_rate15, vc_rate9, vc_rate5, vc_rate3, vc_rate1, vc_rate0,
    # Household/housing
    household_smokers_t6, household_size_t6, household_size_t7,
    asthma_env_risk,
    broken_windows, exposed_wires, wall_cracks, floor_holes, peeling_paint,
    darkness, crowding, clutter, dirtiness, mice_rats,
    #high_asthma_risk, asthma_risk,
    neighborhood_trash_yn, home_trash_yn,
    # Asthma outcomes
    newasthmadx_22, newasthmadx_15, anyasthmadx_22, anyasthmadx_15,
    age_asthma_dx, asthma_meds, hosp_asthma, hosp_asthma_yn,
    limit_asthma, pcg_asthma
  )

# Examine missingness
#vis_miss(df_reduced)

# Factor coding
relevant_factor_vars <- c("edu_t7", "edu_t6", "edu_t5", "gvh1600cat", "gvh1600count_cat")
df_reduced[relevant_factor_vars] <- lapply(df_reduced[relevant_factor_vars], as.factor)

df_reduced$race <- factor(
  df_reduced$race,
  levels = c(1, 2, 3, 4),
  labels = c("White", "Black", "Hispanic", "Other")
)
df_reduced$race <- relevel(df_reduced$race, ref = "Black")


# =============================================================================
# [16] ANALYTIC SAMPLE PREPARATION ----
# Restricted to participants without asthma at age 15 (anyasthmadx_15 == 0)
# Listwise deletion on all primary model covariates
# Primary analytic sample: n=1,936; incident asthma cases: n=109 (5.6%)
# =============================================================================

df_reduced_main <- df_reduced %>%
  dplyr::select(
    idnum,
    age_t7, age_t6, female, race,
    smoking_t7, smoking_t6, vaping_t7, household_smokers_t6, household_size_t6,
    newasthmadx_22,
    gvh1600count, gvh1600count_log, gvh1600_any, gvh1600count_cat,
    gvh1000count, gvh500count, gvh500_any,
    gvh1600_rate, gvh1000_rate, gvh500_rate,
    edu_t7, edu_t6, pov_t7, pov_t6, trct_pct_pov_t6, vc_rate15
  )

df_main_nomiss <- df_reduced_main %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(
    gvh1600count_std    = as.numeric(scale(gvh1600count)),
    gvh1000count_std    = as.numeric(scale(gvh1000count)),
    gvh500count_std     = as.numeric(scale(gvh500count)),
    gvh1600rate_std     = as.numeric(scale(gvh1600_rate)),
    gvh1000rate_std     = as.numeric(scale(gvh1000_rate)),
    gvh500rate_std      = as.numeric(scale(gvh500_rate)),
    vc_rate15_std       = as.numeric(scale(vc_rate15)),
    pov_t7_std          = as.numeric(scale(pov_t7)),
    pov_t6_std          = as.numeric(scale(pov_t6)),
    trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6)),
    gvh1600count_top    = ifelse(gvh1600count > 15, 15, gvh1600count)
  )

cat("Primary analytic sample n =", nrow(df_main_nomiss), "\n")                         # n=1,936
cat("Incident asthma cases n =", sum(df_main_nomiss$newasthmadx_22, na.rm = TRUE), "\n") # n=109

# --- [16.1] Attrition Analysis ----
# Tests whether loss to follow-up at Y22 is associated with exposure
# Eligible: participants without asthma at Y15 (n=2,690)
# Lost: n=441 (16.4%) — those missing both youth and PCG report at Y22

df_eligible <- df_reduced %>%
  filter(anyasthmadx_15 == 0) %>%
  mutate(
    attrited    = ifelse(is.na(newasthmadx_22), 1, 0),
    gvh1600count_std = as.numeric(scale(gvh1600count))
  ) %>%
  filter(!is.na(gvh1600count_std) & !is.na(female))

m_attrition <- glm(
  attrited ~ gvh1600count_std + female,
  data = df_eligible, family = binomial(link = "logit")
)
coeftest(m_attrition, vcov. = vcovHC, type = "HC1")
exp(coef(m_attrition))
exp(coefci(m_attrition, vcov. = vcovHC, type = "HC1"))
# Result: exposure OR=1.11 (p=.044); female OR=0.62 (p<.001)
# Direction: higher-exposed participants more likely to attrite
# → Primary estimates are conservative (downward bias)


# =============================================================================
# [17] CHECK CLUSTERING ----
# Tests whether multilevel modeling is warranted
# =============================================================================

# State-level clustering (**NEED TO BRING IN STATE AND TRACT IDENTIFIERS**)
m_null_state <- glmer(
  newasthmadx_22 ~ 1 + (1 | state_t6),
  data   = df_main_nomiss,
  family = binomial(link = "logit")
)
icc(m_null_state)
# ICC = 0.039 — minimal state-level clustering; standard regression sufficient

# Census tract-level clustering
m_null_tract <- glmer(
  newasthmadx_22 ~ 1 + (1 | trct_id6),
  data   = df_main_nomiss,
  family = binomial(link = "logit")
)
icc(m_null_tract)
# ICC ≈ 0 / singular fit — no meaningful tract-level clustering


# =============================================================================
# [18] MAIN MODELS: INCIDENT ASTHMA (n=1,936) ----
# Sequential logistic regression; M5 = primary fully adjusted model
# All models use HC1 robust standard errors
# =============================================================================

out_main_list <- list(
  
  # M1: demographics + smoking/vaping
  m1 = glm(newasthmadx_22 ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7,
           data = df_main_nomiss, family = binomial(link = "logit")),
  
  # M2: + individual SES (Y15 and Y22)
  m2 = glm(newasthmadx_22 ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 + 
             pov_t7_std + edu_t7 + edu_t6 + pov_t6_std,
           data = df_main_nomiss, family = binomial(link = "logit")),
  
  # M3: + race/ethnicity
  m3 = glm(newasthmadx_22 ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 + 
             pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race,
           data = df_main_nomiss, family = binomial(link = "logit")),
  
  # M4: + neighborhood SES and county violent crime
  m4 = glm(newasthmadx_22 ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 + 
             pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std,
           data = df_main_nomiss, family = binomial(link = "logit")),
  
  # M5: fully adjusted [PRIMARY MODEL]
  # OR = 1.22 (95% CI: 1.03–1.45, p=.022)
  m5 = glm(newasthmadx_22 ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 + 
             pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
             household_smokers_t6 + household_size_t6,
           data = df_main_nomiss, family = binomial(link = "logit")),
  
  # M6: sex interaction [PRIMARY MODERATION MODEL]
  # Male OR = 1.50 (95% CI: 1.18–1.91); Female OR = 1.03 (NS); interaction p=.006
  m6_sex = glm(newasthmadx_22 ~ gvh1600rate_std * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 + 
                 pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                 household_smokers_t6 + household_size_t6,
               data = df_main_nomiss, family = binomial(link = "logit"))
)

# Results with robust SEs
lapply(out_main_list, function(m) {
  list(
    n       = nobs(m),
    robust  = coeftest(m, vcov. = vcovHC, type = "HC1"),
    ci      = coefci(m,   vcov. = vcovHC, type = "HC1"),
    OR      = exp(coef(m)["gvh1600rate_std"]),
    OR_ci   = exp(coefci(m, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])
  )
})

# Sex-specific ORs from M6 (linear combination: main + interaction)
int_term <- grep("gvh1600rate_std:female|female:gvh1600rate_std",
                 names(coef(out_main_list$m6_sex)), value = TRUE)
# Females
robust_ci_lincomb(out_main_list$m6_sex, setNames(c(1, 1), c("gvh1600rate_std", int_term)))
# Males
robust_ci_term(out_main_list$m6_sex, "gvh1600rate_std")


# M5 main effect
exp(coef(out_main_list$m5)["gvh1600rate_std"])
exp(coefci(out_main_list$m5, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])

# Main results summary:
# M5 overall:  OR=1.25 (95% CI: 1.06–1.49), p=.010
# M6 males:    OR=1.55 (95% CI: 1.24–1.94)
# M6 females:  OR=1.06 (95% CI: 0.87–1.30), NS
# Interaction: p=.006

# output for supplemental tables 

# ── Helper: extract all estimates from one model ──────────────────────────────
extract_model <- function(mod, label) {
  ct  <- coeftest(mod, vcov. = vcovHC, type = "HC1")
  ci  <- coefci(mod,   vcov. = vcovHC, type = "HC1")
  
  terms <- rownames(ct)
  
  tibble(
    term  = terms,
    OR    = exp(ct[, "Estimate"]),
    lower = exp(ci[, 1]),
    upper = exp(ci[, 2]),
    pval  = ct[, "Pr(>|z|)"],
    model = label
  )
}

# ── Format a single cell: "OR (LL–UL)" + p ───────────────────────────────────
fmt_cell <- function(OR, lower, upper, pval) {
  p_str <- ifelse(pval < .001, "< .001",
                  ifelse(pval < .01,  formatC(pval, digits = 3, format = "f"),
                         formatC(pval, digits = 3, format = "f")))
  sprintf("%.2f (%.2f–%.2f)\np = %s", OR, lower, upper, p_str)
}

# ── Extract all models ────────────────────────────────────────────────────────
results <- bind_rows(
  extract_model(out_main_list$m1,     "M1"),
  extract_model(out_main_list$m2,     "M2"),
  extract_model(out_main_list$m3,     "M3"),
  extract_model(out_main_list$m4,     "M4"),
  extract_model(out_main_list$m5,     "M5"),
  extract_model(out_main_list$m6_sex, "M6")
)

# ── Build wide table (terms as rows, models as columns) ──────────────────────
results_wide <- results %>%
  mutate(cell = mapply(fmt_cell, OR, lower, upper, pval)) %>%
  select(term, model, cell) %>%
  tidyr::pivot_wider(names_from = model, values_from = cell, values_fill = "—")

# ── Define covariate display order and labels ─────────────────────────────────
term_labels <- c(
  "(Intercept)"                    = "Intercept",
  "gvh1600rate_std"                = "Gun Violence Exposure (1,600m), std.",
  "gvh1600rate_std:female"         = "Gun Violence × Female (interaction)",
  "female"                         = "Female (ref: Male)",
  "age_t7"                         = "Age at Age 22 Assessment, std.",
  "smoking_t6"                     = "Smoking, Past Month (Age 15)",
  "smoking_t7"                     = "Smoking, Past Month (Age 22)",
  "vaping_t7"                      = "Vaping, Past Year (Age 22)",
  "pov_t6_std"                     = "Income-to-Poverty Ratio (Age 15), std.",
  "pov_t7_std"                     = "Income-to-Poverty Ratio (Age 22), std.",
  "edu_t62"                        = "Parental Education (Age 15): HS or Equivalent",
  "edu_t63"                        = "Parental Education (Age 15): Some College",
  "edu_t64"                        = "Parental Education (Age 15): College+",
  "edu_t72"                        = "Own Education (Age 22): HS or Equivalent",
  "edu_t73"                        = "Own Education (Age 22): Some College",
  "edu_t74"                        = "Own Education (Age 22): College+",
  "raceWhite"                      = "Race/Ethnicity: Non-Hispanic White (ref: NH Black)",
  "raceHispanic"                   = "Race/Ethnicity: Hispanic",
  "raceOther"                      = "Race/Ethnicity: Multiracial/Other",
  "trct_pct_pov_t6_std"            = "Census Tract Poverty Rate (Age 15), std.",
  "vc_rate15_std"                  = "County Violent Crime Rate (Age 15), std.",
  "household_smokers_t6"           = "Household Smokers (Age 15)",
  "household_size_t6"              = "Household Size (Age 15), std."
)

# Reorder rows to match term_labels order; keep any unlabeled terms at end
results_wide <- results_wide %>%
  mutate(
    order = match(term, names(term_labels)),
    label = ifelse(term %in% names(term_labels), term_labels[term], term)
  ) %>%
  arrange(order) %>%
  select(label, M1, M2, M3, M4, M5, M6)

# ── N per model ───────────────────────────────────────────────────────────────
n_row <- tibble(
  label = "N",
  M1 = as.character(nobs(out_main_list$m1)),
  M2 = as.character(nobs(out_main_list$m2)),
  M3 = as.character(nobs(out_main_list$m3)),
  M4 = as.character(nobs(out_main_list$m4)),
  M5 = as.character(nobs(out_main_list$m5)),
  M6 = as.character(nobs(out_main_list$m6_sex))
)

results_wide <- bind_rows(n_row, results_wide)

# ── flextable ─────────────────────────────────────────────────────────────────
ft <- flextable(results_wide) %>%
  set_header_labels(
    label = "Covariate",
    M1 = "Model 1",
    M2 = "Model 2",
    M3 = "Model 3",
    M4 = "Model 4",
    M5 = "Model 5\n(Primary)",
    M6 = "Model 6\n(Sex Interaction)"
  ) %>%
  bold(i = 1, part = "header") %>%
  bold(j = 1) %>%
  # Highlight exposure row
  bold(i = ~ label == "Gun Violence Exposure (1,600m), std.") %>%
  bg(i = ~ label == "Gun Violence Exposure (1,600m), std.", bg = "#f0f0f0") %>%
  # Highlight interaction row
  bold(i = ~ label == "Gun Violence × Female (interaction)") %>%
  bg(i = ~ label == "Gun Violence × Female (interaction)", bg = "#f0f0f0") %>%
  fontsize(size = 9, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  align(j = 2:7, align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  width(j = 1, width = 2.2) %>%
  width(j = 2:7, width = 1.15) %>%
  add_footer_lines("Note: Estimates are odds ratios (ORs) with 95% confidence intervals and p-values from logistic regression models with heteroskedasticity-robust standard errors (HC1). All continuous covariates standardized (mean = 0, SD = 1). Model 5 is the primary fully adjusted specification. Model 6 replaces the main effect of sex with a Gun Violence × Female interaction term; sex-specific ORs for Model 6 are reported in the main text. — = covariate not included in model.") %>%
  fontsize(size = 8, part = "footer") %>%
  set_table_properties(layout = "autofit")

# ── Print / export ────────────────────────────────────────────────────────────
ft  # renders in RStudio Viewer; copy from there into Word



# =============================================================================
# [19] SECONDARY ANALYTIC DATASETS ----
# Separate complete-case datasets per outcome to avoid unnecessary case loss
# from unrelated outcome missingness. Re-standardized within each subsample.
# =============================================================================

# Helper: re-standardize within analytic subsample
add_std_vars <- function(data) {
  data %>% mutate(
    gvh1600count_std    = as.numeric(scale(gvh1600count)),
    gvh1000count_std    = as.numeric(scale(gvh1000count)),
    gvh500count_std     = as.numeric(scale(gvh500count)),
    gvh1600rate_std     = as.numeric(scale(gvh1600_rate)),
    gvh1000rate_std     = as.numeric(scale(gvh1000_rate)),
    gvh500rate_std      = as.numeric(scale(gvh500_rate)),
    vc_rate15_std       = as.numeric(scale(vc_rate15)),
    pov_t7_std          = as.numeric(scale(pov_t7)),
    pov_t6_std          = as.numeric(scale(pov_t6)),
    trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6))
  )
}

# Shared covariate columns used across all secondary datasets
# Note: gvh1600_rate corrected from earlier gvh1600rate (no underscore)
secondary_covariates <- c(
  "idnum", "age_t7", "age_t6", "female", "race",
  "smoking_t7", "smoking_t6", "vaping_t7", "household_smokers_t6", "household_size_t6",
  "gvh1600count", "gvh1000count", "gvh500count",
  "gvh1600_any", "gvh500_any", "gvh1600count_cat",
  "gvh1600_rate", "gvh1000_rate", "gvh500_rate",
  "edu_t7", "edu_t6", "pov_t7", "pov_t6", "trct_pct_pov_t6", "vc_rate15"
)

# --- [19.1] Incident Dx datasets (primary secondary outcomes) ----
# Source: df_main_nomiss filtered to incident cases; outcomes joined from df_reduced
# Only outcome variable drives case loss

# [19.1a] Hospitalization among incident Dx (n=103)
df_newdx_hosp <- df_main_nomiss %>%
  filter(newasthmadx_22 == 1) %>%
  left_join(dplyr::select(df_reduced, idnum, hosp_asthma_yn), by = "idnum") %>%
  filter(!is.na(hosp_asthma_yn)) %>%
  as.data.frame()
cat("Incident Dx | Hospitalization n =", nrow(df_newdx_hosp), "\n")  # n=103

# [19.1b] Severity among incident Dx (n=106)
df_newdx_sev <- df_main_nomiss %>%
  filter(newasthmadx_22 == 1) %>%
  left_join(dplyr::select(df_reduced, idnum, limit_asthma), by = "idnum") %>%
  filter(!is.na(limit_asthma)) %>%
  as.data.frame()
cat("Incident Dx | Severity n =", nrow(df_newdx_sev), "\n")  # n=106

# --- [19.2] Any Dx datasets (supplemental secondary outcomes) ----
# Source: df_reduced filtered to anyasthmadx_22 == 1

# [19.2a] Hospitalization among any Dx (n=527)
df_anydx_hosp <- df_reduced %>%
  filter(anyasthmadx_22 == 1) %>%
  dplyr::select(all_of(secondary_covariates), hosp_asthma_yn) %>%
  na.omit() %>%
  as.data.frame() %>%
  add_std_vars()
cat("Any Dx | Hospitalization n =", nrow(df_anydx_hosp), "\n")  # n=527

# [19.2b] Severity among any Dx (n=533)
df_anydx_sev <- df_reduced %>%
  filter(anyasthmadx_22 == 1) %>%
  dplyr::select(all_of(secondary_covariates), limit_asthma) %>%
  na.omit() %>%
  as.data.frame() %>%
  add_std_vars()
cat("Any Dx | Severity n =", nrow(df_anydx_sev), "\n")  # n=533


# =============================================================================
# [20] HOSPITALIZATION MODELS ----
# =============================================================================

# Primary: incident Dx sample (n=103)
out_hosp_new_list <- list(
  m1 = glm(hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7,
           data = df_newdx_hosp, family = binomial(link = "logit")),
  m2 = glm(hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
             pov_t7_std + edu_t7 + edu_t6 + pov_t6_std,
           data = df_newdx_hosp, family = binomial(link = "logit")),
  m3 = glm(hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
             pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race,
           data = df_newdx_hosp, family = binomial(link = "logit")),
  m4 = glm(hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
             pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std,
           data = df_newdx_hosp, family = binomial(link = "logit")),
  # M5: fully adjusted [PRIMARY]
  m5 = glm(hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
             pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
             household_smokers_t6 + household_size_t6,
           data = df_newdx_hosp, family = binomial(link = "logit")),
  # M6: sex interaction — NS (p=.631)
  m6_sex = glm(hosp_asthma_yn ~ gvh1600rate_std * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
                 pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                 household_smokers_t6 + household_size_t6,
               data = df_newdx_hosp, family = binomial(link = "logit"))
)
# Note: warnings re fitted probabilities of 0/1 expected given small n and binary outcome.
# edu_t74 shows large SE due to quasi-complete separation — does not affect primary estimate.
# Unable to test for moderation by race due to small cell sizes. 

lapply(out_hosp_new_list, function(m) summary(m))
lapply(out_hosp_new_list, function(m) nobs(m))

# Robust SEs and CIs: M5 primary result — OR = 4.59 (95% CI: 1.56–13.47, p=.006)
exp(coef(out_hosp_new_list$m5)["gvh1600rate_std"])
exp(coefci(out_hosp_new_list$m5, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])

exp(coef(out_hosp_new_list$m4)["gvh1600rate_std"])
exp(coefci(out_hosp_new_list$m4, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])

exp(coef(out_hosp_new_list$m3)["gvh1600rate_std"])
exp(coefci(out_hosp_new_list$m3, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])

exp(coef(out_hosp_new_list$m2)["gvh1600rate_std"])
exp(coefci(out_hosp_new_list$m2, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])

exp(coef(out_hosp_new_list$m1)["gvh1600rate_std"])
exp(coefci(out_hosp_new_list$m1, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])

# sex interaaction: p = .862
coeftest(out_hosp_new_list$m6, vcov. = vcovHC, type = "HC1")
coefci(out_hosp_new_list$m6,   vcov. = vcovHC, type = "HC1", level = 0.95)

# Supplemental: any Dx sample (n=527)
out_hosp_any_list <- list(
  m1 = glm(hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7,
           data = df_anydx_hosp, family = binomial(link = "logit")),
  m2 = glm(hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
             pov_t7_std + edu_t7 + edu_t6 + pov_t6_std,
           data = df_anydx_hosp, family = binomial(link = "logit")),
  m3 = glm(hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
             pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race,
           data = df_anydx_hosp, family = binomial(link = "logit")),
  m4 = glm(hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
             pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std,
           data = df_anydx_hosp, family = binomial(link = "logit")),
  # M5: fully adjusted — OR = 1.53 (95% CI: 1.19–1.96, p=.001)
  m5 = glm(hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
             pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
             household_smokers_t6 + household_size_t6,
           data = df_anydx_hosp, family = binomial(link = "logit")),
  # M6: sex interaction — NS (interaction p=.101)
  m6_sex = glm(hosp_asthma_yn ~ gvh1600rate_std * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
                 pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                 household_smokers_t6 + household_size_t6,
               data = df_anydx_hosp, family = binomial(link = "logit"))
)

lapply(out_hosp_any_list, summary)
lapply(out_hosp_any_list, function(m) nobs(m))

exp(coef(out_hosp_any_list$m6)["gvh1600rate_std"])
exp(coefci(out_hosp_any_list$m6, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])
exp(coef(out_hosp_any_list$m6)["gvh1600rate_std:female"])
exp(coefci(out_hosp_any_list$m6, vcov. = vcovHC, type = "HC1")["gvh1600rate_std:female", ])

exp(coef(out_hosp_any_list$m5)["gvh1600rate_std"])
exp(coefci(out_hosp_any_list$m5, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])

exp(coef(out_hosp_any_list$m4)["gvh1600rate_std"])
exp(coefci(out_hosp_any_list$m4, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])

exp(coef(out_hosp_any_list$m3)["gvh1600rate_std"])
exp(coefci(out_hosp_any_list$m3, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])

exp(coef(out_hosp_any_list$m2)["gvh1600rate_std"])
exp(coefci(out_hosp_any_list$m2, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])

exp(coef(out_hosp_any_list$m1)["gvh1600rate_std"])
exp(coefci(out_hosp_any_list$m1, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])

# generate table for supplemental materials 
# ── Helper: extract all estimates from one model ──────────────────────────────
extract_model <- function(mod, label) {
  ct  <- coeftest(mod, vcov. = vcovHC, type = "HC1")
  ci  <- coefci(mod,   vcov. = vcovHC, type = "HC1")
  
  tibble(
    term  = rownames(ct),
    OR    = exp(ct[, "Estimate"]),
    lower = exp(ci[, 1]),
    upper = exp(ci[, 2]),
    pval  = ct[, "Pr(>|z|)"],
    model = label
  )
}

# ── Format cell: OR (LL–UL) only, p on second line ───────────────────────────
fmt_cell <- function(OR, lower, upper, pval) {
  p_str <- ifelse(pval < .001, "< .001", formatC(pval, digits = 3, format = "f"))
  sprintf("%.2f (%.2f\u2013%.2f)", OR, lower, upper)
}

# ── Term labels (adjust suffixes to match your actual coef names) ─────────────
term_labels <- c(
  "(Intercept)"          = "Intercept",
  "gvh1600rate_std"      = "Gun Violence Exposure (1,600m), std.",
  "gvh1600rate_std:female" = "Gun Violence \u00d7 Female (interaction)",
  "female"               = "Female (ref: Male)",
  "age_t7"               = "Age at Age 22 Assessment, std.",
  "smoking_t6"           = "Smoking, Past Month (Age 15)",
  "smoking_t7"           = "Smoking, Past Month (Age 22)",
  "vaping_t7"            = "Vaping, Past Year (Age 22)",
  "pov_t6_std"           = "Income-to-Poverty Ratio (Age 15), std.",
  "pov_t7_std"           = "Income-to-Poverty Ratio (Age 22), std.",
  "edu_t62"              = "Parental Education (Age 15): HS or Equivalent",
  "edu_t63"              = "Parental Education (Age 15): Some College",
  "edu_t64"              = "Parental Education (Age 15): College+",
  "edu_t72"              = "Own Education (Age 22): HS or Equivalent",
  "edu_t73"              = "Own Education (Age 22): Some College",
  "edu_t74"              = "Own Education (Age 22): College+",
  "raceWhite"            = "Race/Ethnicity: Non-Hispanic White (ref: NH Black)",
  "raceHispanic"         = "Race/Ethnicity: Hispanic",
  "raceMultiracial"      = "Race/Ethnicity: Multiracial/Other",
  "trct_pct_pov_t6_std"  = "Census Tract Poverty Rate (Age 15), std.",
  "vc_rate15_std"        = "County Violent Crime Rate (Age 15), std.",
  "household_smokers_t6" = "Household Smokers (Age 15)",
  "household_size_t6"    = "Household Size (Age 15), std."
)

shared_footer <- "Note: Estimates are odds ratios (ORs) with 95% confidence intervals from logistic regression models with heteroskedasticity-robust standard errors (HC1). All continuous covariates standardized (mean = 0, SD = 1). Model 5 is the primary fully adjusted specification. Model 6 replaces the main effect of sex with a Gun Violence \u00d7 Female interaction term. \u2014 = covariate not included in model."

# ── Generic table-builder function ───────────────────────────────────────────
build_hosp_table <- function(model_list, model_names, col_labels, n_row_vals,
                             highlight_interaction = TRUE, footer_note) {
  results <- bind_rows(
    mapply(extract_model, model_list, model_names, SIMPLIFY = FALSE)
  )
  
  results_wide <- results %>%
    mutate(cell = mapply(fmt_cell, OR, lower, upper, pval)) %>%
    select(term, model, cell) %>%
    pivot_wider(names_from = model, values_from = cell, values_fill = "\u2014") %>%
    mutate(
      order = match(term, names(term_labels)),
      label = ifelse(term %in% names(term_labels), term_labels[term], term)
    ) %>%
    arrange(order) %>%
    select(label, all_of(model_names))
  
  n_row <- as_tibble(
    c(list(label = "N"), setNames(as.list(as.character(n_row_vals)), model_names))
  )
  results_wide <- bind_rows(n_row, results_wide)
  
  # Build flextable with do.call to avoid !!! splicing issue
  ft <- flextable(results_wide)
  ft <- do.call(set_header_labels,
                c(list(x = ft), list(label = "Covariate"), as.list(col_labels)))
  ft <- ft %>%
    bold(i = 1, part = "header") %>%
    bold(j = 1) %>%
    bg(i = ~ label == "Gun Violence Exposure (1,600m), std.", bg = "#f0f0f0") %>%
    bold(i = ~ label == "Gun Violence Exposure (1,600m), std.")
  
  if (highlight_interaction) {
    ft <- ft %>%
      bg(i = ~ label == "Gun Violence \u00d7 Female (interaction)", bg = "#f0f0f0") %>%
      bold(i = ~ label == "Gun Violence \u00d7 Female (interaction)")
  }
  
  ft %>%
    fontsize(size = 9, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    align(j = 2:(length(model_names) + 1), align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "all") %>%
    width(j = 1, width = 2.2) %>%
    width(j = 2:(length(model_names) + 1), width = 1.1) %>%
    add_footer_lines(footer_note) %>%
    fontsize(size = 8, part = "footer") %>%
    set_table_properties(layout = "autofit")
}

# ── Table S2a: Incident Diagnosis Sample (n = 103) ───────────────────────────
ft_hosp_new <- build_hosp_table(
  model_list = list(
    out_hosp_new_list$m1,
    out_hosp_new_list$m2,
    out_hosp_new_list$m3,
    out_hosp_new_list$m4,
    out_hosp_new_list$m5,
    out_hosp_new_list$m6_sex
  ),
  model_names = c("M1", "M2", "M3", "M4", "M5", "M6"),
  col_labels  = c(
    M1 = "Model 1", M2 = "Model 2", M3 = "Model 3",
    M4 = "Model 4", M5 = "Model 5\n(Primary)", M6 = "Model 6\n(Sex Interaction)"
  ),
  n_row_vals  = c(
    nobs(out_hosp_new_list$m1), nobs(out_hosp_new_list$m2),
    nobs(out_hosp_new_list$m3), nobs(out_hosp_new_list$m4),
    nobs(out_hosp_new_list$m5), nobs(out_hosp_new_list$m6_sex)
  ),
  footer_note = paste(shared_footer,
                      "Separation warnings for edu_t74 reflect small cell sizes and do not affect the primary exposure estimate.")
)

ft_hosp_new  # renders in RStudio Viewer

# ── Table S2b: Any Diagnosis Sample (n = 527) ────────────────────────────────
ft_hosp_any <- build_hosp_table(
  model_list = list(
    out_hosp_any_list$m1,
    out_hosp_any_list$m2,
    out_hosp_any_list$m3,
    out_hosp_any_list$m4,
    out_hosp_any_list$m5,
    out_hosp_any_list$m6_sex
  ),
  model_names = c("M1", "M2", "M3", "M4", "M5", "M6"),
  col_labels  = c(
    M1 = "Model 1", M2 = "Model 2", M3 = "Model 3",
    M4 = "Model 4", M5 = "Model 5\n(Primary)", M6 = "Model 6\n(Sex Interaction)"
  ),
  n_row_vals  = c(
    nobs(out_hosp_any_list$m1), nobs(out_hosp_any_list$m2),
    nobs(out_hosp_any_list$m3), nobs(out_hosp_any_list$m4),
    nobs(out_hosp_any_list$m5), nobs(out_hosp_any_list$m6_sex)
  ),
  footer_note = shared_footer
)

ft_hosp_any 



# =============================================================================
# [21] SEVERITY MODELS (FUNCTIONAL LIMITATION) ----
# =============================================================================

df_newdx_sev$limit_asthma_std <- scale(df_newdx_sev$limit_asthma)

# Primary: incident Dx sample (n=106)
out_sev_new_list <- list(
  m1 = lm(limit_asthma_std ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7,
          data = df_newdx_sev),
  m2 = lm(limit_asthma_std ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
            pov_t7_std + edu_t7 + edu_t6 + pov_t6_std,
          data = df_newdx_sev),
  m3 = lm(limit_asthma_std ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
            pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race,
          data = df_newdx_sev),
  m4 = lm(limit_asthma_std ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
            pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std,
          data = df_newdx_sev),
  # M5: fully adjusted [PRIMARY]
  m5 = lm(limit_asthma_std ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
            pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
            household_smokers_t6 + household_size_t6,
          data = df_newdx_sev),
  # M6: sex interaction — NS (interaction p=.834)
  m6_sex = lm(limit_asthma_std ~ gvh1600rate_std * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
                pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                household_smokers_t6 + household_size_t6,
              data = df_newdx_sev)
)

lapply(out_sev_new_list, function(m) summary(m))
lapply(out_sev_new_list, function(m) nobs(m))

# Robust SEs and CIs: M6 interaction: B = -0.03 (95% CI: -0.36–0.29, p=.834)
coeftest(out_sev_new_list$m6_sex, vcov. = vcovHC, type = "HC1")
coefci(out_sev_new_list$m6_sex,   vcov. = vcovHC, type = "HC1", level = 0.95)

# Robust SEs and CIs: M5 primary result — B = 0.25 (95% CI: 0.06–0.44, p=.010)
coeftest(out_sev_new_list$m5, vcov. = vcovHC, type = "HC1")
coefci(out_sev_new_list$m5,   vcov. = vcovHC, type = "HC1", level = 0.95)

coeftest(out_sev_new_list$m4, vcov. = vcovHC, type = "HC1")
coefci(out_sev_new_list$m4,   vcov. = vcovHC, type = "HC1", level = 0.95)

coeftest(out_sev_new_list$m3, vcov. = vcovHC, type = "HC1")
coefci(out_sev_new_list$m3,   vcov. = vcovHC, type = "HC1", level = 0.95)

coeftest(out_sev_new_list$m2, vcov. = vcovHC, type = "HC1")
coefci(out_sev_new_list$m2,   vcov. = vcovHC, type = "HC1", level = 0.95)

coeftest(out_sev_new_list$m1, vcov. = vcovHC, type = "HC1")
coefci(out_sev_new_list$m1,   vcov. = vcovHC, type = "HC1", level = 0.95)


# Supplemental: any Dx sample (n=533)
df_anydx_sev$limit_asthma_std <- scale(df_anydx_sev$limit_asthma)
out_sev_any_list <- list(
  m1 = lm(limit_asthma_std ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7,
          data = df_anydx_sev),
  m2 = lm(limit_asthma_std ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
            pov_t7_std + edu_t7 + edu_t6 + pov_t6_std,
          data = df_anydx_sev),
  m3 = lm(limit_asthma_std ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
            pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race,
          data = df_anydx_sev),
  m4 = lm(limit_asthma_std ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
            pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std,
          data = df_anydx_sev),
  # M5: fully adjusted — B = 0.10 (95% CI: 0.02–0.17, p=.017)
  m5 = lm(limit_asthma_std ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
            pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + trct_pct_pov_t6_std + vc_rate15_std +
            race + household_smokers_t6 + household_size_t6,
          data = df_anydx_sev),
  # M6: sex interaction — NS (interaction p=.930)
  m6_sex = lm(limit_asthma_std ~ gvh1600rate_std * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
                pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + trct_pct_pov_t6_std + vc_rate15_std +
                race + household_smokers_t6 + household_size_t6,
              data = df_anydx_sev)
)

lapply(out_sev_any_list, summary)
lapply(out_sev_any_list, function(m) nobs(m))

coeftest(out_sev_any_list$m6_sex, vcov. = vcovHC, type = "HC1")
coefci(out_sev_any_list$m6_sex,   vcov. = vcovHC, type = "HC1", level = 0.95)

coeftest(out_sev_any_list$m5, vcov. = vcovHC, type = "HC1")
coefci(out_sev_any_list$m5,   vcov. = vcovHC, type = "HC1", level = 0.95)

coeftest(out_sev_any_list$m4, vcov. = vcovHC, type = "HC1")
coefci(out_sev_any_list$m4,   vcov. = vcovHC, type = "HC1", level = 0.95)

coeftest(out_sev_any_list$m3, vcov. = vcovHC, type = "HC1")
coefci(out_sev_any_list$m3,   vcov. = vcovHC, type = "HC1", level = 0.95)

coeftest(out_sev_any_list$m2, vcov. = vcovHC, type = "HC1")
coefci(out_sev_any_list$m2,   vcov. = vcovHC, type = "HC1", level = 0.95)

coeftest(out_sev_any_list$m1, vcov. = vcovHC, type = "HC1")
coefci(out_sev_any_list$m1,   vcov. = vcovHC, type = "HC1", level = 0.95)



# =============================================================================
# [22] SUPPLEMENTAL: SPATIAL SENSITIVITY ANALYSES ----
# Tests whether associations hold at closer residential proximities (1000m, 500m)
# All models use fully adjusted M5 covariate specification
# Sample: df_main_nomiss (n=1,936)
# =============================================================================

out_spatial_list <- list(
  
  # --- Main effects ---
  
  # 1600m [primary; included for comparison]
  # OR = 1.22 (95% CI: 1.03–1.45, p=.022)
  m_1600 = glm(newasthmadx_22 ~ gvh1600_rate + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
                 pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                 household_smokers_t6 + household_size_t6,
               data = df_main_nomiss, family = binomial(link = "logit")),
  
  # 1000m — significant; slightly attenuated vs. 1600m
  # OR = 1.21 (95% CI: 1.02–1.44, p=.031)
  m_1000 = glm(newasthmadx_22 ~ gvh1000_rate + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
                 pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                 household_smokers_t6 + household_size_t6,
               data = df_main_nomiss, family = binomial(link = "logit")),
  
  # 500m — NS in fully adjusted model (p=.166)
  # OR = 1.13 (95% CI: 0.95–1.34)
  m_500  = glm(newasthmadx_22 ~ gvh500_rate + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
                 pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                 household_smokers_t6 + household_size_t6,
               data = df_main_nomiss, family = binomial(link = "logit")),
  
  # --- Sex interaction models ---
  
  # 1600m × sex: significant interaction (p=.006)
  # Male OR = 1.50 (95% CI: 1.21–1.85); Female OR = 1.03 (NS)
  m_1600_sex = glm(newasthmadx_22 ~ gvh1600_rate * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
                     pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                     household_smokers_t6 + household_size_t6,
                   data = df_main_nomiss, family = binomial(link = "logit")),
  
  # 1000m × sex: significant interaction (p=.041)
  # Male OR = 1.42 (95% CI: 1.15–1.76); Female OR = 1.06 (NS)
  m_1000_sex = glm(newasthmadx_22 ~ gvh1000_rate * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
                     pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                     household_smokers_t6 + household_size_t6,
                   data = df_main_nomiss, family = binomial(link = "logit")),
  
  # 500m × sex: marginal interaction (p=.092)
  # Male OR = 1.36 (95% CI: 1.07–1.73); Female OR = 1.07 (NS)
  m_500_sex  = glm(newasthmadx_22 ~ gvh500_rate * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
                     pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                     household_smokers_t6 + household_size_t6,
                   data = df_main_nomiss, family = binomial(link = "logit"))
)

# Robust SEs, CIs, and n for all spatial models
lapply(out_spatial_list, function(m) {
  list(coeftest = coeftest(m, vcov. = vcovHC, type = "HC1"),
       coefci   = coefci(m,   vcov. = vcovHC, type = "HC1", level = 0.95),
       n        = nobs(m))
})

# Summary table of primary exposure ORs across distances
spatial_summary <- lapply(names(out_spatial_list), function(nm) {
  m       <- out_spatial_list[[nm]]
  exp_var <- grep("rate$", names(coef(m)), value = TRUE)[1]
  ci      <- coefci(m, vcov. = vcovHC, type = "HC1")
  data.frame(
    model    = nm,
    exposure = exp_var,
    OR       = round(exp(coef(m)[exp_var]), 3),
    lower95  = round(exp(ci[exp_var, 1]), 3),
    upper95  = round(exp(ci[exp_var, 2]), 3)
  )
})
do.call(rbind, spatial_summary)

# Results summary (main effects):
# 1600m: OR=1.22 (1.03–1.45), p=.022 [primary]
# 1000m: OR=1.21 (1.02–1.44), p=.031 [consistent]
# 500m:  OR=1.13 (0.95–1.34), p=.166 [NS; attenuated at closer proximity]
#
# Results summary (sex interactions):
# 1600m × sex: interaction p=.006 — males OR=1.50 (1.21–1.85); females NS
# 1000m × sex: interaction p=.041 — males OR=1.42 (1.15–1.76); females NS
# 500m  × sex: interaction p=.092 — marginal; males OR=1.36 (1.07–1.73); females NS


# =============================================================================
# [23] SUPPLEMENTAL: NULL FINDINGS AT AGE 15 ----
# Tests whether DGV at age 15 was associated with asthma at age 15
# Two outcomes: incident dx Y9→Y15 (newasthmadx_15) and prevalent dx at Y15
# (anyasthmadx_15). Separate datasets per outcome to avoid unnecessary case loss.
# Expected: NS for gun violence exposure — confirms no pre-existing association
# that could confound primary Y22 findings.
# =============================================================================

y15_covariates <- c(
  "idnum", "age_t6", "female", "race", "smoking_t6",
  "household_smokers_t6", "household_size_t6",
  "gvh1600count", "gvh1600_rate", "edu_t6", "pov_t6", "trct_pct_pov_t6", "vc_rate15"
)

# --- [23.1] Incident asthma Y9→Y15 (n=2,372) ----
df_newdx_y15 <- df_reduced %>%
  dplyr::select(all_of(y15_covariates), newasthmadx_15) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(
    gvh1600count_std    = as.numeric(scale(gvh1600count)),
    gvh1600rate_std    = as.numeric(scale(gvh1600_rate)),
    vc_rate15_std       = as.numeric(scale(vc_rate15)),
    pov_t6_std          = as.numeric(scale(pov_t6)),
    trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6))
  )
df_newdx_y15$edu_t6 <- factor(df_newdx_y15$edu_t6)
cat("Incident dx Y9→Y15 n =", nrow(df_newdx_y15), "\n")  # n=2,372

# --- [23.2] Prevalent asthma at Y15 (n=3,279) ----
df_anydx_y15 <- df_reduced %>%
  dplyr::select(all_of(y15_covariates), anyasthmadx_15) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(
    gvh1600count_std    = as.numeric(scale(gvh1600count)),
    gvh1600rate_std    = as.numeric(scale(gvh1600_rate)),
    vc_rate15_std       = as.numeric(scale(vc_rate15)),
    pov_t6_std          = as.numeric(scale(pov_t6)),
    trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6))
  )
df_anydx_y15$edu_t6 <- factor(df_anydx_y15$edu_t6)
cat("Prevalent dx at Y15 n =", nrow(df_anydx_y15), "\n")  # n=3,279

# --- [23.3] Models ----

# Incident asthma Y9→Y15 — Result: OR=1.00 (95% CI: 0.87–1.14, p=.953) — NS
m_newdx_y15 <- glm(
  newasthmadx_15 ~ gvh1600rate_std + age_t6 + female + smoking_t6 +
    edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
    household_smokers_t6 + household_size_t6,
  data = df_newdx_y15, family = binomial(link = "logit"))

summary(m_newdx_y15)
coeftest(m_newdx_y15, vcov. = vcovHC, type = "HC1")
coefci(m_newdx_y15,   vcov. = vcovHC, type = "HC1", level = 0.95)
exp(coef(m_newdx_y15)["gvh1600rate_std"])
exp(coefci(m_newdx_y15, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])
# OR=1.02 (95% CI: 0.89–1.17), p=.743

# Prevalent asthma at Y15 — Result: OR=1.00 (95% CI: 0.92–1.09, p=.991) — NS
m_anydx_y15 <- glm(
  anyasthmadx_15 ~ gvh1600rate_std + age_t6 + female + smoking_t6 +
    edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
    household_smokers_t6 + household_size_t6,
  data = df_anydx_y15, family = binomial(link = "logit"))

summary(m_anydx_y15)
coeftest(m_anydx_y15, vcov. = vcovHC, type = "HC1")
coefci(m_anydx_y15,   vcov. = vcovHC, type = "HC1", level = 0.95)
exp(coef(m_anydx_y15)["gvh1600rate_std"])
exp(coefci(m_anydx_y15, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])
# OR=1.01 (95% CI: 0.92–1.10), p=.870

# Summary of null Y15 findings:
# Incident dx Y9→Y15: OR=1.02 (95% CI: 0.89–1.17), p=.743 — NS
# Prevalent dx at Y15: OR=1.01 (95% CI: 0.92–1.10), p=.870 — NS
# Interpretation: gun violence exposure at Y15 was not associated with concurrent
# or prior asthma, supporting a prospective relationship with Y22 outcomes
# and ruling out reverse causation.


# =============================================================================
# [24] SUPPLEMENTAL: BUILT ENVIRONMENT SENSITIVITY ----
# Adjusts for interviewer-rated housing/neighborhood conditions
# In-home observations conducted only for hard-to-reach/non-respondent
# participants — subsample is approximately half the full sample (n=529)
# Only 31 incident asthma cases; interpret with caution
# Primary question: does the association persist after adjusting for
# environmental asthma triggers that co-occur with gun violence?
# =============================================================================

# --- [24.1] Analytic dataset (n=529; incident cases n=31) ----
be_covariates <- c(
  "idnum", "age_t7", "female", "race",
  "smoking_t7", "smoking_t6", "vaping_t7", "household_smokers_t6", "household_size_t6",
  "newasthmadx_22",
  "gvh1600count", "gvh1000count", "gvh500count", "gvh1600_rate",
  "gvh1600_any", "gvh1600count_cat",
  "edu_t7", "edu_t6", "pov_t7", "pov_t6", "trct_pct_pov_t6", "vc_rate15",
  "asthma_env_risk"
)

df_be <- df_reduced %>%
  dplyr::select(all_of(be_covariates)) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(
    gvh1600count_std    = as.numeric(scale(gvh1600count)),
    gvh1600rate_std    = as.numeric(scale(gvh1600_rate)),
    gvh1000count_std    = as.numeric(scale(gvh1000count)),
    gvh500count_std     = as.numeric(scale(gvh500count)),
    vc_rate15_std       = as.numeric(scale(vc_rate15)),
    pov_t7_std          = as.numeric(scale(pov_t7)),
    pov_t6_std          = as.numeric(scale(pov_t6)),
    trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6))
  )

cat("Built environment subsample n =", nrow(df_be), "\n")                              # n=512
cat("Incident asthma cases n =", sum(df_be$newasthmadx_22, na.rm = TRUE), "\n")       # n=28
# Note: small case count (n=31) limits power; interpret with caution
# edu_t7 shows quasi-complete separation in this subsample — robust SEs used throughout

# --- [24.2] Models ----
out_be_list <- list(
  
  # M1: asthma-specific environmental risk index (peeling_paint, dirtiness, mice_rats, neighborhood_trash_yn)
  # Result: OR=1.43 (95% CI: 1.04–1.97, p=.029) — exposure remains significant
  m1_risk = glm(
    newasthmadx_22 ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
      household_smokers_t6 + household_size_t6 + asthma_env_risk,
    data = df_be, family = binomial(link = "logit")),
  
  # M2: sex interaction + built environment indicators
  # Result: interaction NS (p=.415); low power given n=31 cases
  m2_sex = glm(
    newasthmadx_22 ~ gvh1600rate_std * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
      household_smokers_t6 + household_size_t6 + asthma_env_risk,
    data = df_be, family = binomial(link = "logit"))
)

# Results with robust SEs, CIs, and ORs
lapply(out_be_list, function(m) {
  list(
    summary = summary(m),
    n       = nobs(m),
    robust  = coeftest(m, vcov. = vcovHC, type = "HC1"),
    ci      = coefci(m,   vcov. = vcovHC, type = "HC1"),
    OR      = exp(coef(m)["gvh1600rate_std"]),
    OR_ci   = exp(coefci(m, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])
  )
})

# Summary:
# M1 (adjust asthma-specific environmental risk index):  OR=1.60 (95% CI: 1.08–2.39), p=.025
# M2 (sex interaction):  interaction NS (p=.217); unreliable given n=31 cases
# Interpretation: association persists after built environment adjustment;
# effect not explained by correlated environmental asthma triggers.

exp(coef(out_be_list$m1)["asthma_env_risk"])
exp(coefci(out_be_list$m1, vcov. = vcovHC, type = "HC1")["asthma_env_risk", ])


# Sex-specific ORs from (linear combination: main + interaction)
int_term <- grep("gvh1600rate_std:female|female:gvh1600rate_std",
                 names(coef(out_be_list$m2)), value = TRUE)
# Females
robust_ci_lincomb(out_be_list$m2, setNames(c(1, 1), c("gvh1600rate_std", int_term)))
# Males
robust_ci_term(out_be_list$m2, "gvh1600rate_std")


# =============================================================================
# [25] SUPPLEMENTAL: ADDITIONAL SENSITIVITY CHECKS ----
# Tests robustness to maternal prenatal smoking, youth BMI, PCG-reported asthma
# Separate datasets per model; all yield n=1,938 (no additional case loss)
# Key finding: OR range 1.22–1.23 — highly stable across all specifications
# =============================================================================

# --- [25.1] Analytic datasets ----

# Maternal prenatal smoking (n=1,938)
df_sens_momsmk <- df_reduced %>%
  dplyr::select(idnum, age_t7, female, race, smoking_t7, smoking_t6, vaping_t7, 
                household_smokers_t6, household_size_t6, newasthmadx_22,
                gvh1600count, gvh1600_rate, edu_t7, edu_t6, pov_t7, pov_t6, trct_pct_pov_t6, vc_rate15, momsmk) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(gvh1600count_std = as.numeric(scale(gvh1600count)), gvh1600rate_std = as.numeric(scale(gvh1600_rate)), 
         vc_rate15_std = as.numeric(scale(vc_rate15)),
         pov_t7_std = as.numeric(scale(pov_t7)), pov_t6_std = as.numeric(scale(pov_t6)),
         trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6)))
cat("momsmk sensitivity n =", nrow(df_sens_momsmk), "\n")  # n=1,938

# Youth BMI at Y22 (n=1,938)
df_sens_bmi <- df_reduced %>%
  dplyr::select(idnum, age_t7, female, race, smoking_t7, smoking_t6, vaping_t7, 
                household_smokers_t6, household_size_t6, newasthmadx_22,
                gvh1600count, gvh1600_rate,edu_t7, edu_t6, pov_t7, pov_t6, trct_pct_pov_t6, vc_rate15, bmiz_t7) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(gvh1600count_std = as.numeric(scale(gvh1600count)), gvh1600rate_std = as.numeric(scale(gvh1600_rate)), 
         vc_rate15_std = as.numeric(scale(vc_rate15)),
         pov_t7_std = as.numeric(scale(pov_t7)), pov_t6_std = as.numeric(scale(pov_t6)),
         trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6)))
cat("BMI sensitivity n =", nrow(df_sens_bmi), "\n")  # n=1,938

# PCG-reported asthma dx (n=1,938)
df_sens_pcg <- df_reduced %>%
  dplyr::select(idnum, age_t7, female, race, smoking_t7, smoking_t6, vaping_t7, 
                household_smokers_t6, household_size_t6, newasthmadx_22,
                gvh1600count, gvh1600_rate,edu_t7, edu_t6, pov_t7, pov_t6, trct_pct_pov_t6, vc_rate15, pcg_asthma) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(gvh1600count_std = as.numeric(scale(gvh1600count)), gvh1600rate_std = as.numeric(scale(gvh1600_rate)), 
         vc_rate15_std = as.numeric(scale(vc_rate15)),
         pov_t7_std = as.numeric(scale(pov_t7)), pov_t6_std = as.numeric(scale(pov_t6)),
         trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6)))
cat("PCG asthma sensitivity n =", nrow(df_sens_pcg), "\n")  # n=1,938

# All three sensitivity covariates combined (n=1,938)
df_sens_all <- df_reduced %>%
  dplyr::select(idnum, age_t7, female, race, smoking_t7, smoking_t6, vaping_t7, 
                household_smokers_t6, household_size_t6, newasthmadx_22,
                gvh1600count, gvh1600_rate, edu_t7, edu_t6, pov_t7, pov_t6, trct_pct_pov_t6, vc_rate15,
                momsmk, bmiz_t7, pcg_asthma) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(gvh1600count_std = as.numeric(scale(gvh1600count)), gvh1600rate_std = as.numeric(scale(gvh1600_rate)), 
         vc_rate15_std = as.numeric(scale(vc_rate15)),
         pov_t7_std = as.numeric(scale(pov_t7)), pov_t6_std = as.numeric(scale(pov_t6)),
         trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6)))
cat("All sensitivity covariates n =", nrow(df_sens_all), "\n")  # n=1,938

# --- [25.2] Models ----
out_sens_list <- list(
  
  # M1: + maternal prenatal smoking — OR=1.22 (1.03–1.45, p=.021); momsmk NS (p=.733)
  m1_momsmk = glm(
    newasthmadx_22 ~ gvh1600rate_std + momsmk + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
      household_smokers_t6 + household_size_t6,
    data = df_sens_momsmk, family = binomial(link = "logit")),
  
  # M2: + youth BMI at Y22 — OR=1.22 (1.03–1.45, p=.021); BMI NS (p=.372)
  m2_bmi = glm(
    newasthmadx_22 ~ gvh1600rate_std + bmiz_t7 + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
      household_smokers_t6 + household_size_t6,
    data = df_sens_bmi, family = binomial(link = "logit")),
  
  # M3: + PCG-reported asthma — OR=1.23 (1.03–1.46, p=.021); PCG asthma NS (p=.586)
  m3_pcg = glm(
    newasthmadx_22 ~ gvh1600rate_std + pcg_asthma + age_t7 + female + smoking_t6 + smoking_t7 +
      vaping_t7 + pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
      household_smokers_t6 + household_size_t6,
    data = df_sens_pcg, family = binomial(link = "logit")),
  
  # M4: all three simultaneously [most conservative] — OR=1.23 (1.03–1.46, p=.020)
  m4_all = glm(
    newasthmadx_22 ~ gvh1600rate_std + momsmk + bmiz_t7 + pcg_asthma + age_t7 + female +
      smoking_t6 + smoking_t7 + vaping_t7 + pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race +
      trct_pct_pov_t6_std + vc_rate15_std + household_smokers_t6 + household_size_t6,
    data = df_sens_all, family = binomial(link = "logit")),
  
  # M5: sex interaction + all — interaction p=.006; male OR=1.51 (1.22–1.87)
  m5_sex = glm(
    newasthmadx_22 ~ gvh1600rate_std * female + momsmk + bmiz_t7 + pcg_asthma + age_t7 +
      smoking_t6 + smoking_t7 + vaping_t7 + pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race +
      trct_pct_pov_t6_std + vc_rate15_std + household_smokers_t6 + household_size_t6,
    data = df_sens_all, family = binomial(link = "logit"))
)

# Results with robust SEs, CIs, and ORs
lapply(names(out_sens_list), function(nm) {
  m <- out_sens_list[[nm]]
  
  # Identify which exposure term is present in this model
  exp_var <- intersect(
    c("gvh1600rate_std", "gvh1600rate_std"),
    names(coef(m))
  )[1]
  
  list(
    model  = nm,
    n      = nobs(m),
    robust = coeftest(m, vcov. = vcovHC, type = "HC1"),
    ci     = coefci(m,   vcov. = vcovHC, type = "HC1"),
    OR     = exp(coef(m)[exp_var]),
    OR_ci  = exp(coefci(m, vcov. = vcovHC, type = "HC1")[exp_var, ])
  )
})

# Summary:
# M1 (+ prenatal smoking): OR=1.26 (1.06–1.49), p=.009
# M2 (+ BMI):              OR=1.26 (1.06–1.49), p=.009
# M3 (+ PCG asthma):       OR=1.26 (1.06–1.50), p=.009
# M4 (+ all three):        OR=1.26 (1.06–1.49), p=.009
# M5 (sex interaction):    male OR=1.56 (1.25–1.95), interaction p=.007


# =============================================================================
# [26] SUPPLEMENTAL: PROPENSITY SCORE MATCHED ANALYSIS ----
# Full matching on Y15 neighborhood characteristics to balance exposed vs.
# unexposed participants on pre-exposure socioeconomic factors.
# Restricted to participants without asthma at Y15.
# =============================================================================

# --- [26.1] PSM analytic dataset ----
# Pre-match: n=1,640 (exposed n=791, unexposed n=849; incident cases n=94)
df_psm <- df_reduced %>%
  dplyr::select(
    idnum, gvh1600_any, gvh1600count, gvh1000_any, gvh1000count, gvh500_any, gvh500count, gvh1600_rate,
    age_t7, age_t6, edu_t6, pov_t7, pov_t6, pov_t5, race, female,
    smoking_t6, smoking_t7, vaping_t7, bmiz_t6, bmiz_t7,
    trct_pct_pov_t6, trct_pct_pov_t5, vc_rate15, vc_rate9,
    household_size_t6, pcg_asthma, household_smokers_t6,
    newasthmadx_22, anyasthmadx_15
  ) %>%
  filter(anyasthmadx_15 == 0) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(
    gvh1600count_std    = as.numeric(scale(gvh1600count)),
    gvh1600rate_std    = as.numeric(scale(gvh1600_rate)),
    pov_t6_std          = as.numeric(scale(pov_t6)),
    pov_t5_std          = as.numeric(scale(pov_t5)),
    pov_t7_std          = as.numeric(scale(pov_t7)),
    trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6)),
    trct_pct_pov_t5_std = as.numeric(scale(trct_pct_pov_t5)),
    vc_rate15_std       = as.numeric(scale(vc_rate15)),
    vc_rate9_std        = as.numeric(scale(vc_rate9))
  )

cat("PSM analytic sample n =", nrow(df_psm), "\n")               # n=1,640
cat("Exposed (any DGV) n =", sum(df_psm$gvh1600_any == 1), "\n") # n=791
cat("Unexposed n =", sum(df_psm$gvh1600_any == 0), "\n")          # n=849
cat("Incident asthma cases n =", sum(df_psm$newasthmadx_22, na.rm = TRUE), "\n")  # n=94

# --- [26.2] Full matching on Y15 neighborhood characteristics ----
# Exact matching on race; calipers on continuous variables to prevent poor matches
m_out_psm <- matchit(
  gvh1600_any ~ race + pov_t6_std + edu_t6 + trct_pct_pov_t6_std +
    trct_pct_pov_t5_std + pov_t5_std + vc_rate9_std,
  data        = df_psm,
  method      = "full",
  distance    = "glm",
  link        = "logit",
  estimand    = "att",
  caliper     = c(pov_t6_std = .5, trct_pct_pov_t6_std = .5,
                  trct_pct_pov_t5_std = .5, pov_t5_std = .5, 
                  vc_rate9_std = 1),
  std.caliper = TRUE,
  exact       = ~ race
)

# Balance diagnostics — all SMDs < 0.10 except edu_t63 (SMD=0.102, marginally above threshold)
# Race perfectly balanced (SMD=0.000) due to exact matching
# 411 participants unmatched due to caliper restrictions (209 control, 202 treated)
summary(m_out_psm, un = FALSE)
love.plot(m_out_psm, threshold = 0.10)
plot(m.out.main, type = "jitter", interactive = FALSE)

# --- [26.3] Matched dataset ----
# Matched sample: n=1,229 (treated n=589, control ESS=182)
matched_data_psm <- match.data(m_out_psm, subclass = "subclass")
cat("Matched sample n =", nrow(matched_data_psm), "\n")  # n=1,229
table(matched_data_psm$newasthmadx_22)                   # 0=1,162; 1=67

# re-standardize exposure within treated group 
library(Hmisc)
standardize_att <- function(data, var_name, treated_idx, weights_var = "weights") {
  # Get the variable and weights
  var <- data[[var_name]]
  weights <- data[[weights_var]]
  
  # Calculate weighted mean and SD for treated group
  mean_treated <- wtd.mean(var[treated_idx], weights = weights[treated_idx])
  sd_treated <- sqrt(wtd.var(var[treated_idx], weights = weights[treated_idx]))
  
  # Standardize the entire variable
  standardized <- (var - mean_treated) / sd_treated
  
  return(standardized)
}

# Identify treated observations once
treated_idx <- matched_data_psm$gvh1600_any == 1

# Apply standardization to all variables
vars_to_standardize <- c("gvh1600_rate")

for (var in vars_to_standardize) {
  new_var_name <- paste0(var, "_att_std")
  matched_data_psm[[new_var_name]] <- standardize_att(
    matched_data_psm, 
    var, 
    treated_idx
  )
}
summary(matched_data_psm$gvh1600_rate_att_std)

# --- [26.4] Weighted regression on matched sample ----
# Cluster-robust SEs clustered on matched subclass
# Result: OR=1.22 (95% CI: 1.00–1.48, p=.047) — consistent with primary model
m_psm_weighted <- glm(
  newasthmadx_22 ~ gvh1600_rate_att_std  + female + age_t7 + pov_t7_std +
    smoking_t6 + smoking_t7 + vaping_t7 +
    household_smokers_t6 + household_size_t6 + vc_rate15_std,
  family  = quasibinomial(link = "logit"),
  data    = matched_data_psm,
  weights = weights
)
nobs(m_psm_weighted)  # n=1,229
coeftest(m_psm_weighted, vcov. = vcovCL, cluster = ~subclass)
coefci(m_psm_weighted,   vcov. = vcovCL, cluster = ~subclass, level = 0.95)
exp(coef(m_psm_weighted)["gvh1600_rate_att_std"])
exp(coefci(m_psm_weighted, vcov. = vcovCL, cluster = ~subclass)["gvh1600_rate_att_std", ])
# OR=1.25 (95% CI: 0.99–1.59), p=.061

m_psm_weighted_int <- glm(
  newasthmadx_22 ~ gvh1600_rate_att_std * female + age_t7 + pov_t7_std +
    smoking_t6 + smoking_t7 + vaping_t7 +
    household_smokers_t6 + household_size_t6 + vc_rate15_std,
  family  = quasibinomial(link = "logit"),
  data    = matched_data_psm,
  weights = weights
)

ct_int  <- coeftest(m_psm_weighted_int, vcov. = vcovCL, cluster = ~subclass)
ci_int  <- coefci(m_psm_weighted_int,   vcov. = vcovCL, cluster = ~subclass, level = 0.95)

# Interaction term OR = 0.75 (0.46 - 1.23), p=.256
int_term <- grep("gvh1600_rate_att_std:female|female:gvh1600_rate_att_std",
                 rownames(ct_int), value = TRUE)
cat("Interaction term OR:\n")
cat(sprintf("OR = %.2f (95%% CI: %.2f-%.2f), p = %.3f\n",
            exp(ct_int[int_term, "Estimate"]),
            exp(ci_int[int_term, 1]),
            exp(ci_int[int_term, 2]),
            ct_int[int_term, "Pr(>|z|)"]))

# Summary: matched OR virtually identical to primary (OR=1.22 vs. 1.22),
# confirming association is not explained by pre-exposure neighborhood confounding.

# --- Sample size comparison: primary vs. matched, overall and by sex ---

# Primary analytic sample (df_main_nomiss or equivalent — adjust name as needed)
cat("=== PRIMARY SAMPLE ===\n")
cat("Overall n =", nrow(df_main_nomiss), "\n")
cat("Incident asthma cases n =", sum(df_main_nomiss$newasthmadx_22, na.rm = TRUE), "\n")
table(Sex = df_main_nomiss$female)
df_main_nomiss %>%
  dplyr::group_by(female) %>%
  dplyr::summarise(total = dplyr::n(), cases = sum(newasthmadx_22, na.rm = TRUE))

# PSM pre-match sample
cat("\n=== PSM PRE-MATCH SAMPLE ===\n")
cat("Overall n =", nrow(df_psm), "\n")
cat("Incident asthma cases n =", sum(df_psm$newasthmadx_22, na.rm = TRUE), "\n")
df_psm %>%
  dplyr::group_by(female) %>%
  dplyr::summarise(total = dplyr::n(), cases = sum(newasthmadx_22, na.rm = TRUE))

# PSM matched sample
cat("\n=== PSM MATCHED SAMPLE ===\n")
cat("Overall n =", nrow(matched_data_psm), "\n")
cat("Incident asthma cases n =", sum(matched_data_psm$newasthmadx_22, na.rm = TRUE), "\n")
matched_data_psm %>%
  dplyr::group_by(female) %>%
  dplyr::summarise(
    total     = dplyr::n(),
    cases     = sum(newasthmadx_22, na.rm = TRUE),
    exposed   = sum(gvh1600_any == 1),
    pct_cases = round(cases / total * 100, 1)
  )

# Effective sample size reduction
cat("\n=== REDUCTION ===\n")
cat("Primary -> matched, overall:", nrow(df_main_nomiss), "->", nrow(matched_data_psm), "\n")
cat("% retained:", round(nrow(matched_data_psm) / nrow(df_main_nomiss) * 100, 1), "%\n")



# =============================================================================
# [27] DESCRIPTIVE STATISTICS TABLE ----
# =============================================================================

mean_sd <- function(x) sprintf("%.1f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
n_pct   <- function(x) {
  n <- sum(x == 1, na.rm = TRUE); d <- sum(!is.na(x))
  sprintf("%d (%.1f%%)", n, 100 * n / d)
}

# Join hospitalization and severity outcomes from their respective subsamples
# hosp_asthma_yn: df_newdx_hosp (incident Dx; n=103)
# limit_asthma:   df_newdx_sev  (incident Dx; n=106)
df_desc <- df_main_nomiss %>%
  dplyr::select(idnum, female, newasthmadx_22, gvh1600_any, gvh1600_rate,
                pov_t6, trct_pct_pov_t6) %>%
  left_join(dplyr::select(df_newdx_hosp, idnum, hosp_asthma_yn), by = "idnum") %>%
  left_join(dplyr::select(df_newdx_sev,  idnum, limit_asthma),   by = "idnum")

desc_summary <- function(data, group_label) {
  data %>% summarise(
    Group                    = group_label,
    `Incident Asthma N (%)`  = n_pct(newasthmadx_22),
    `Hospitalization N (%)`  = n_pct(hosp_asthma_yn),
    `Severity M (SD)`        = mean_sd(limit_asthma),
    `Any GV Exposure N (%)`  = n_pct(gvh1600_any),
    `GV Rate M (SD)`        = mean_sd(gvh1600_rate),
    `Poverty Ratio M (SD)`   = mean_sd(pov_t6),
    `Tract Poverty % M (SD)` = mean_sd(trct_pct_pov_t6)
  )
}

table1 <- bind_rows(
  desc_summary(df_desc, "Overall"),
  desc_summary(df_desc %>% filter(female == 0), "Male"),
  desc_summary(df_desc %>% filter(female == 1), "Female")
)
print(table1)

table(df_main_nomiss$female)
# =============================================================================
# [28] FIGURES ----
# =============================================================================

# Shared theme
theme_poster <- theme_classic(base_size = 14) +
  theme(
    axis.text         = element_text(colour = "black"),
    axis.title        = element_text(colour = "black"),
    plot.title        = element_text(face = "bold", size = 13),
    legend.background = element_rect(fill = NA, colour = NA)
  )
TICK_TEXT <- 12

# --- Figure 1: Forest plot — overall and sex-specific ORs (incident asthma) ---
gvh_term <- "gvh1600rate_std"
int_term  <- grep("gvh1600rate_std:female|female:gvh1600rate_std",
                  names(coef(out_main_list$m6_sex)), value = TRUE)

plot_tbl_f1 <- bind_rows(
  robust_ci_term(out_main_list$m5,     gvh_term) %>% mutate(group = "Overall"),
  robust_ci_term(out_main_list$m6_sex, gvh_term) %>% mutate(group = "Male"),
  robust_ci_lincomb(out_main_list$m6_sex,
                    setNames(c(1, 1), c(gvh_term, int_term))) %>% mutate(group = "Female")
) %>%
  mutate(
    # Step 1: relabel to display strings
    group = dplyr::recode(group,
                          "Overall" = "Overall (N = 1,936)",
                          "Male"    = "Male (n = 872)",
                          "Female"  = "Female (n = 1,064)"
    ),
    # Step 2: order for plot (bottom to top on y-axis)
    group = factor(group, levels = c(
      "Overall (N = 1,936)",
      "Male (n = 872)",
      "Female (n = 1,064)"
    )),
    label = sprintf("%.2f (%.2f\u2013%.2f)", estimate, l95, u95)
  )

fig1_forest <- ggplot(plot_tbl_f1, aes(x = estimate, y = group)) +
  geom_errorbar(aes(xmin = l95, xmax = u95), width = 0, linewidth = 1, alpha = 0.9,
                orientation = "y") +
  geom_point(size = 4) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 1.2) +
  geom_text(aes(label = label), hjust = -0.10, vjust = 1.5, size = 5.5) +
  scale_x_continuous(trans = "log10", limits = c(0.8, 2),
                     breaks = c(0.8, 1, 2), labels = c("0.8", "1.0", "2.0"),
                     position = "bottom") +
  labs(title = "Odds of incident Asthma by Age 22",
       x = "Odds Ratio (per 1 SD exposure)", y = NULL) +
  theme_poster +
  theme(axis.text.y  = element_text(size = TICK_TEXT, face = "bold", colour = "black"),
        axis.title.y = element_blank())

# --- Figure 2: Predicted probabilities by sex ---
# Uses unstandardized count (gvh1600count) for interpretable x-axis
new_dx_1600_int <- glm(
  newasthmadx_22 ~ gvh1600_rate * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 + 
    pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + trct_pct_pov_t6_std + vc_rate15_std +
    race + household_smokers_t6 + household_size_t6,
  data = df_main_nomiss, family = binomial(link = "logit"))

pred <- ggpredict(new_dx_1600_int, terms = c("gvh1600_rate [all]", "female")) %>%
  as_tibble() %>%
  mutate(sex = ifelse(group == "0", "Male", "Female"))

fig2_pred <- ggplot(pred, aes(x = x, y = predicted, colour = sex, fill = sex)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 1.4) +
  scale_colour_manual(values = c(Male = "#0072B2", Female = "#D55E00")) +
  scale_fill_manual(values   = c(Male = "#0072B2", Female = "#D55E00")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Predicted Probability of Incident Asthma by Age 22",
       x = "Annual Exposure Rate (Age 15)", y = NULL, colour = NULL, fill = NULL) +
  theme_poster +
  theme(legend.position      = c(0.01, 0.99),
        legend.justification = c(0, 1))

# --- Figure 3: Hospitalization ORs (incident and any Dx) ---
plot_tbl_hosp <- bind_rows(
  robust_ci_term(out_hosp_new_list$m5, gvh_term) %>% mutate(group = "Incident Dx\n(n=103)   "),
  robust_ci_term(out_hosp_any_list$m5, gvh_term) %>% mutate(group = "Any Dx\n(n=527)")
) %>%
  mutate(group = factor(group, levels = c("Any Dx\n(n=527)", "Incident Dx\n(n=103)   ")),
         label = sprintf("%.2f (%.2f\u2013%.2f)", estimate, l95, u95))

fig3_hosp <- ggplot(plot_tbl_hosp, aes(x = estimate, y = group)) +
  geom_errorbar(aes(xmin = l95, xmax = u95), width = 0, linewidth = 1, alpha = 0.9,
                orientation = "y") +
  geom_point(size = 4) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 1.2) +
  geom_text(aes(label = label), hjust = -0.10, vjust = 1.5, size = 5.5) +
  scale_x_continuous(trans = "log10", limits = c(0.8, 30),
                     breaks = c(1, 2, 5, 25), labels = c("1", "2", "5", "25"),
                     position = "bottom") +
  labs(title = "Odds of Asthma Hospitalization by Age 22",
       x = "Odds Ratio (per 1 SD exposure)", y = NULL) +
  theme_poster +
  theme(axis.text.y = element_text(face = "bold"))

# --- Figure 4: Severity betas (incident and any Dx) ---
# Models already fitted in Sections 20–21 — no refitting needed
plot_tbl_sev <- bind_rows(
  robust_ci_term(out_sev_new_list$m5, gvh_term, transform = identity) %>%
    mutate(group = "Incident Dx\n(n=106)   "),
  robust_ci_term(out_sev_any_list$m5, gvh_term, transform = identity) %>%
    mutate(group = "Any Dx\n(n=533)")
) %>%
  mutate(group = factor(group, levels = c("Any Dx\n(n=533)", "Incident Dx\n(n=106)   ")),
         label = sprintf("%.2f (%.2f\u2013%.2f)", estimate, l95, u95))

fig4_sev <- ggplot(plot_tbl_sev, aes(x = estimate, y = group)) +
  geom_errorbar(aes(xmin = l95, xmax = u95), width = 0, linewidth = 1, alpha = 0.9,
                orientation = "y") +
  geom_point(size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 1.2) +
  geom_text(aes(label = label), hjust = -0.10, vjust = 1.5, size = 5.5) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), position = "bottom") +
  labs(title = "Associations with Asthma Severity at Age 22",
       x = "\u03b2 (per 1 SD exposure)", y = NULL) +
  theme_poster +
  theme(axis.text.y = element_text(face = "bold"))

# --- Print and save ---
print(fig1_forest)
print(fig2_pred)
print(fig3_hosp)
print(fig4_sev)

ggsave("/Users/co453621/Library/CloudStorage/OneDrive-UniversityofCentralFlorida/2. MANUSCRIPTS/MARTZ GUN VIOLENCE ASTHMA/plots/fig1_forest_overall_sex.png", fig1_forest, width = 7.5, height = 6, dpi = 500)
ggsave("/Users/co453621/Library/CloudStorage/OneDrive-UniversityofCentralFlorida/2. MANUSCRIPTS/MARTZ GUN VIOLENCE ASTHMA/plots/fig2_pred_prob_sex.png",      fig2_pred,   width = 7.5, height = 6, dpi = 500)
ggsave("/Users/co453621/Library/CloudStorage/OneDrive-UniversityofCentralFlorida/2. MANUSCRIPTS/MARTZ GUN VIOLENCE ASTHMA/plots/fig3_hosp_OR.png",            fig3_hosp,   width = 7.5, height = 6, dpi = 500)
ggsave("/Users/co453621/Library/CloudStorage/OneDrive-UniversityofCentralFlorida/2. MANUSCRIPTS/MARTZ GUN VIOLENCE ASTHMA/plots/fig4_severity_B.png",         fig4_sev,    width = 7.5, height = 6, dpi = 500)


# =============================================================================