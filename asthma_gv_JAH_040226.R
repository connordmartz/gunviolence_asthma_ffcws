# =============================================================================
# ADOLESCENT EXPOSURE TO DEADLY GUN VIOLENCE & ASTHMA BURDEN IN EARLY ADULTHOOD
# Martz, Navas Nazario & Gaydosh | Journal of Urban Health
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
#     [14]   Helper Functions
#     [15]   Create df_reduced
#     [16]   Analytic Sample Preparation (df_main_nomiss)
#     [16.1] Attrition Analysis
#     [16.2] Sample Flowchart Counts (Figure S1)
#     [17]   Check Clustering (ICC)
#     [18]   Main Models: Incident Asthma (Table S1)
#     [19]   Secondary Analytic Datasets
#     [20]   Hospitalization Models: Firth Logistic Regression (Tables S2, S4)
#     [21]   Functional Limitation Models (Tables S3, S5)
#     [22]   Supplemental: Spatial Sensitivity (Figure S2)
#     [23]   Supplemental: Falsification Tests at Age 15
#     [24]   Supplemental: Built Environment Sensitivity
#     [25]   Supplemental: Additional Confounders
#     [26]   Supplemental: Propensity Score Analyses (Figures S3–S4)
#
#   OUTPUTS (Sections 27–28)
#     [27] Descriptive Statistics (Table 1)
#     [28] Figures 1–3
#
#   ADDITIONAL SENSITIVITY ANALYSES (Sections 29–31; Table S6)
#     [29] Exposure Data Quality
#     [30] Exposure Distribution
#     [31] Neighborhood Traffic Density
#
# =============================================================================


# =============================================================================
# [1] LOAD LIBRARIES ----
# =============================================================================

library(lubridate) 
library(haven)     
library(dplyr)
library(tidyr)
library(tibble)
library(lmtest)    
library(sandwich)   
library(lme4)        
library(performance) 
library(logistf)    
library(MatchIt)   
library(cobalt)     
library(ggplot2)
library(ggeffects)  
library(scales)
library(patchwork)
library(flextable)

# =============================================================================
# [2] DATA IMPORT ----
# =============================================================================

core    <- read_dta("/Users/co453621/Library/CloudStorage/OneDrive-UniversityofCentralFlorida/15. DATA/FFCWS ALL WAVES Y0-Y22/Raw Data/ffcws_core_y0_y22.dta")
census00 <- read_csv("/Users/co453621/Library/CloudStorage/OneDrive-UniversityofCentralFlorida/15. DATA/FFCWS ALL WAVES Y0-Y22/Restricted Data/FFCWS Restricted Use Data_GRP - Documents/Restricted Data/ffgeo7_all_pub.dta")
geoid   <- read_csv("/Users/co453621/Library/CloudStorage/OneDrive-UniversityofCentralFlorida/15. DATA/FFCWS ALL WAVES Y0-Y22/Restricted Data/FFCWS Restricted Use Data_GRP - Documents/Restricted Data/contractcity7pub_121824.dta")
FF_gva  <- read_csv("/Users/co453621/Library/CloudStorage/OneDrive-UniversityofCentralFlorida/15. DATA/FFCWS ALL WAVES Y0-Y22/Restricted Data/FFCWS Restricted Use Data_GRP - Documents/Restricted Data/ff_gva_15y_pub1_20190426.dta")
FF_ucr  <- read_csv("/Users/co453621/Library/CloudStorage/OneDrive-UniversityofCentralFlorida/15. DATA/FFCWS ALL WAVES Y0-Y22/Restricted Data/FFCWS Restricted Use Data_GRP - Documents/Restricted Data/ff_UCR_pub1.dta")

# =============================================================================
# [3] MERGE DATA ----
# =============================================================================

df <- core %>%
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

    # --- 10.7 Youth-reported neighborhood traffic, age 15 (higher = more traffic) ---
    traffic_t6 = ifelse(k6e4a < 0, NA, 5 - k6e4a)
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
      p5h1b == 2 & p6b2 %in% c(-9,-3,-2,-1) & p6b27_3 == 0          ~ 0,
      p5h1b == 2 & p6b2 == 1                                        ~ 1,
      p5h1b == 2 & p6b2 %in% c(-9,-3,-2,-1) & p6b27_3 == 1          ~ 1,
      TRUE ~ NA_real_
    ),
    
    # --- 12.2 Incident asthma diagnosis (Y15 to Y22) | PRIMARY OUTCOME ---
    newasthmadx_22 = case_when(
      p6b2  %in% c(1, -2, -9)                                       ~ NA_real_,
      k7i4_1 %in% c(-9,-3,-2,-1) & p7i3_1 %in% c(-9,-6,-3,-1)       ~ NA_real_,
      p6b2 == 2 & k7i4_1 == 0                                       ~ 0,
      p6b2 == 2 & k7i4_1 %in% c(-9,-3,-2,-1) & p7i3_1 == 0          ~ 0,
      p6b2 == 2 & k7i4_1 == 1                                       ~ 1,
      p6b2 == 2 & k7i4_1 %in% c(-9,-3,-2,-1) & p7i3_1 == 1          ~ 1,
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

# Female-specific OR (HC1) from an exposure x female model
female_or <- function(m, exp_var) {
  int <- paste0(exp_var, ":female")
  V   <- vcovHC(m, type = "HC1")
  b   <- coef(m)[exp_var] + coef(m)[int]
  se  <- sqrt(V[exp_var, exp_var] + V[int, int] + 2 * V[exp_var, int])
  round(exp(c(OR = b, lower = b - 1.96 * se, upper = b + 1.96 * se)), 2)
}
 
# All coefficients with HC1 robust 95% CIs (exp = OR; identity = beta)
extract_robust <- function(mod, label, transform = exp) {
  ct <- coeftest(mod, vcov. = vcovHC, type = "HC1")
  ci <- coefci(mod,   vcov. = vcovHC, type = "HC1")
  tibble(term  = rownames(ct),
         est   = transform(ct[, 1]),
         lower = transform(ci[, 1]),
         upper = transform(ci[, 2]),
         pval  = ct[, 4],
         model = label)
}
 
# All coefficients from a Firth model (profile penalized-likelihood CIs)
extract_firth <- function(mod, label) {
  tibble(term  = names(coef(mod)),
         est   = exp(coef(mod)),
         lower = exp(mod$ci.lower),
         upper = exp(mod$ci.upper),
         pval  = mod$prob,
         model = label)
}

# Row labels for supplemental Tables S1-S5
term_labels <- c(
  "gvh1600rate_std"                = "Gun Violence Exposure (1,600m), std.",
  "gvh1600rate_std:female"         = "Gun Violence x Female (interaction)",
  "female"                         = "Female (ref: Male)",
  "age_t7"                         = "Age at Age 22 Assessment, std.",
  "smoking_t6"                     = "Smoking, Past Month (Age 15)",
  "smoking_t7"                     = "Smoking, Past Month (Age 22)",
  "vaping_t7"                      = "Vaping, Past Year (Age 22)",
  "pov_t6_std"                     = "Income-Poverty Ratio (Age 15), std.",
  "pov_t7_std"                     = "Income-Poverty Ratio (Age 22), std.",
  "edu_t62"                        = "Parent Education (Age 15): HS or Equivalent",
  "edu_t63"                        = "Parent Education (Age 15): Some College",
  "edu_t64"                        = "Parent Education (Age 15): College/Graduate Degree",
  "edu_t72"                        = "Own Education (Age 22): HS or Equivalent",
  "edu_t73"                        = "Own Education (Age 22): Some College",
  "edu_t74"                        = "Own Education (Age 22): College/Graduate Degree",
  "edu_t7_binSome College or More" = "Own Education (Age 22): Some College or More (ref: <= HS)",
  "raceWhite"                      = "Race/Ethnicity: Non-Hispanic White (ref: Non-Hispanic Black)",
  "raceHispanic"                   = "Race/Ethnicity: Hispanic",
  "raceOther"                      = "Race/Ethnicity: Multiracial/Other",
  "trct_pct_pov_t6_std"            = "Census Tract % Poverty (Age 15), std.",
  "vc_rate15_std"                  = "County Violent Crime Rate (Age 15), std.",
  "household_smokers_t6"           = "Household Smokers (Age 15)",
  "household_size_t6"              = "Household Size (Age 15), std."
)

# Wide table: terms as rows, Models 1-6 as columns, "est (LL-UL)" cells
make_wide_table <- function(results) {
  results %>%
    filter(term != "(Intercept)") %>%
    mutate(cell = sprintf("%.2f (%.2f-%.2f)", est, lower, upper)) %>%
    dplyr::select(term, model, cell) %>%
    tidyr::pivot_wider(names_from = model, values_from = cell, values_fill = "-") %>%
    mutate(order = match(term, names(term_labels)),
           label = ifelse(term %in% names(term_labels), term_labels[term], term)) %>%
    arrange(order) %>%
    dplyr::select(label, M1, M2, M3, M4, M5, M6)
}
 
# Format a wide table as a flextable and save it to out_dir as .docx
save_table <- function(wide, file) {
  ft <- flextable(wide) %>%
    set_header_labels(label = "Covariate",
                      M1 = "Model 1", M2 = "Model 2", M3 = "Model 3", M4 = "Model 4",
                      M5 = "Model 5 (Primary)", M6 = "Model 6 (Sex Interaction)") %>%
    fontsize(size = 8, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    align(j = 2:7, align = "center", part = "all") %>%
    set_table_properties(layout = "autofit")
  save_as_docx(ft, path = file.path(out_dir, file))
  ft
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
  ) %>%
  # Variables for sensitivity analyses and ICCs, added after listwise
  # deletion so they do not change the analytic sample
  left_join(dplyr::select(df_reduced, idnum, pcg_asthma, state_t6, trct_id6), by = "idnum") %>%
  left_join(dplyr::select(df, idnum, atleast_1year_exposure, traffic_t6),     by = "idnum")

cat("Primary analytic sample n =", nrow(df_main_nomiss), "\n")                         
cat("Incident asthma cases n =", sum(df_main_nomiss$newasthmadx_22, na.rm = TRUE), "\n") 

# --- [16.1] Attrition Analysis ----
# Tests whether loss to follow-up at Y22 is associated with exposure
# Eligible: participants without asthma at Y15 (n=2,690)
# Lost: n=441 (16.4%) — those missing both youth and PCG report at Y22

df_eligible <- df_reduced %>%
  filter(anyasthmadx_15 == 0) %>%
  mutate(attrited = ifelse(is.na(newasthmadx_22), 1, 0)) %>%
  filter(!is.na(gvh1600count) & !is.na(female))
 
mu  <- mean(df_main_nomiss$gvh1600_rate)
sdv <- sd(df_main_nomiss$gvh1600_rate)
df_eligible$gvh1600rate_std_a <- (df_eligible$gvh1600_rate - mu) / sdv
 
m_attrition <- glm(attrited ~ gvh1600rate_std_a + female,
                   family = binomial, data = df_eligible)
nobs(m_attrition)
table(model.frame(m_attrition)$attrited)
coeftest(m_attrition, vcov. = vcovHC, type = "HC0")
exp(coef(m_attrition))
exp(coefci(m_attrition, vcov. = vcovHC, type = "HC0"))

# --- [16.2] Sample Flowchart Counts (Figure S1) ----
step <- function(x, label) cat(sprintf("%-55s n = %5d\n", label, nrow(x)))
 
d0 <- df
step(d0, "0. FFCWS sample")
 
d1 <- d0 %>% filter(ck6kint == 1)
step(d1, "1. Participated at age 15")
 
d2 <- d1 %>% filter(!is.na(gvh1600_rate))
step(d2, "2. Valid age 15 GVA exposure data")
 
d3 <- d2 %>% filter(ck7kint == 1)
step(d3, "3. Participated at age 22")
 
d4 <- d3 %>% filter(!is.na(anyasthmadx_15))
cat("Missing age 15 asthma data:", nrow(d3) - nrow(d4), "\n")
 
d5 <- d4 %>% filter(anyasthmadx_15 == 0)
cat("Asthma diagnosis by age 15:", nrow(d4) - nrow(d5), "\n")
 
d6 <- d5 %>% filter(!is.na(newasthmadx_22))
cat("Missing age 22 asthma outcome:", nrow(d5) - nrow(d6), "\n")
 
d7 <- d6 %>% tidyr::drop_na(age_t7, age_t6, female, smoking_t6, smoking_t7, vaping_t7,
                            edu_t6, edu_t7, pov_t6, pov_t7, race, trct_pct_pov_t6, vc_rate15,
                            household_smokers_t6, household_size_t6)
cat("Missing covariate data:", nrow(d6) - nrow(d7), "\n")
cat("Analytic sample:", nrow(d7), "| cases:", sum(d7$newasthmadx_22), "\n")

# =============================================================================
# [17] CHECK CLUSTERING ----
# Tests whether multilevel modeling is warranted
# =============================================================================

m_null_state <- glmer(newasthmadx_22 ~ 1 + (1 | state_t6),
                      data = df_main_nomiss, family = binomial(link = "logit"))
icc(m_null_state)
 
m_null_tract <- glmer(newasthmadx_22 ~ 1 + (1 | trct_id6),
                      data = df_main_nomiss, family = binomial(link = "logit"))
icc(m_null_tract)

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

# Robust coefficient tests for all models
lapply(out_main_list, coeftest, vcov. = vcovHC, type = "HC1")
 
# Table S1
results_main <- bind_rows(
  extract_robust(out_main_list$m1,     "M1"),
  extract_robust(out_main_list$m2,     "M2"),
  extract_robust(out_main_list$m3,     "M3"),
  extract_robust(out_main_list$m4,     "M4"),
  extract_robust(out_main_list$m5,     "M5"),
  extract_robust(out_main_list$m6_sex, "M6")
)
results_main %>% filter(term %in% c("gvh1600rate_std", "gvh1600rate_std:female"))
table_s1 <- make_wide_table(results_main)
save_table(table_s1, "table_s1_incident_asthma.docx")
 
# Sex-specific ORs from M6 (female = main effect + interaction)
int_term <- "gvh1600rate_std:female"
robust_ci_term(out_main_list$m6_sex, "gvh1600rate_std")                                    # Males
robust_ci_lincomb(out_main_list$m6_sex, setNames(c(1, 1), c("gvh1600rate_std", int_term))) # Females

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

# [19.1a] Hospitalization among incident Dx
df_newdx_hosp <- df_main_nomiss %>%
  filter(newasthmadx_22 == 1) %>%
  left_join(dplyr::select(df_reduced, idnum, hosp_asthma_yn), by = "idnum") %>%
  filter(!is.na(hosp_asthma_yn)) %>%
  as.data.frame()
cat("Incident Dx | Hospitalization n =", nrow(df_newdx_hosp), "\n")
 
# [19.1b] Functional limitation among incident Dx
df_newdx_sev <- df_main_nomiss %>%
  filter(newasthmadx_22 == 1) %>%
  left_join(dplyr::select(df_reduced, idnum, limit_asthma), by = "idnum") %>%
  filter(!is.na(limit_asthma)) %>%
  as.data.frame()
cat("Incident Dx | Functional limitation n =", nrow(df_newdx_sev), "\n")

# --- [19.2] Any Dx datasets (supplemental secondary outcomes) ----

# [19.2a] Hospitalization among any Dx
df_anydx_hosp <- df_reduced %>%
  filter(anyasthmadx_22 == 1) %>%
  dplyr::select(all_of(secondary_covariates), hosp_asthma_yn) %>%
  na.omit() %>%
  as.data.frame() %>%
  add_std_vars()
cat("Any Dx | Hospitalization n =", nrow(df_anydx_hosp), "\n")
 
# [19.2b] Functional limitation among any Dx
df_anydx_sev <- df_reduced %>%
  filter(anyasthmadx_22 == 1) %>%
  dplyr::select(all_of(secondary_covariates), limit_asthma) %>%
  na.omit() %>%
  as.data.frame() %>%
  add_std_vars()
cat("Any Dx | Functional limitation n =", nrow(df_anydx_sev), "\n")

# =============================================================================
# [20] HOSPITALIZATION MODELS: FIRTH PENALIZED LOGISTIC REGRESSION (Tables S2, S4) ----
# Profile penalized-likelihood CIs (logistf default).
# =============================================================================
 
# --- [20.1] Incident Dx sample (Table S2) ----
# Sparse cells: smoking at age 15 omitted; own education at age 22 collapsed
# to binary (some college or more vs. high school or less)
table(df_newdx_hosp$smoking_t6, df_newdx_hosp$hosp_asthma_yn)
table(df_newdx_hosp$edu_t7,     df_newdx_hosp$hosp_asthma_yn)
 
df_newdx_hosp <- df_newdx_hosp %>%
  mutate(
    edu_t7_bin = case_when(
      edu_t7 %in% c("1", "2") ~ "HS or Less",
      edu_t7 %in% c("3", "4") ~ "Some College or More",
      TRUE ~ NA_character_
    ),
    edu_t7_bin = factor(edu_t7_bin, levels = c("HS or Less", "Some College or More"))
  )
table(df_newdx_hosp$edu_t7, df_newdx_hosp$edu_t7_bin, useNA = "always")
 
out_hosp_firth_list <- list(
  m1 = logistf(
    hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t7 + vaping_t7,
    data = df_newdx_hosp),
 
  m2 = logistf(
    hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7_bin + edu_t6 + pov_t6_std,
    data = df_newdx_hosp),
 
  m3 = logistf(
    hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7_bin + edu_t6 + pov_t6_std + race,
    data = df_newdx_hosp),
 
  m4 = logistf(
    hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7_bin + edu_t6 + pov_t6_std + race +
      trct_pct_pov_t6_std + vc_rate15_std,
    data = df_newdx_hosp),
 
  m5 = logistf(
    hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7_bin + edu_t6 + pov_t6_std + race +
      trct_pct_pov_t6_std + vc_rate15_std + household_smokers_t6 + household_size_t6,
    data = df_newdx_hosp),
 
  m6_sex = logistf(
    hosp_asthma_yn ~ gvh1600rate_std * female + age_t7 + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7_bin + edu_t6 + pov_t6_std + race +
      trct_pct_pov_t6_std + vc_rate15_std + household_smokers_t6 + household_size_t6,
    data = df_newdx_hosp)
)
 
lapply(out_hosp_firth_list, summary)
 
results_hosp_new <- bind_rows(
  extract_firth(out_hosp_firth_list$m1,     "M1"),
  extract_firth(out_hosp_firth_list$m2,     "M2"),
  extract_firth(out_hosp_firth_list$m3,     "M3"),
  extract_firth(out_hosp_firth_list$m4,     "M4"),
  extract_firth(out_hosp_firth_list$m5,     "M5"),
  extract_firth(out_hosp_firth_list$m6_sex, "M6")
)
results_hosp_new %>% filter(term %in% c("gvh1600rate_std", "gvh1600rate_std:female"))
table_s2 <- make_wide_table(results_hosp_new)
save_table(table_s2, "table_s2_hosp_incident.docx")
 
# --- [20.2] Any Dx sample (Table S4) ----
# Primary covariate specification (no sparse-cell modifications needed)
table(df_anydx_hosp$smoking_t6, df_anydx_hosp$hosp_asthma_yn)
table(df_anydx_hosp$edu_t7,     df_anydx_hosp$hosp_asthma_yn)
 
out_hosp_firth_anydx_list <- list(
  m1 = logistf(
    hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7,
    data = df_anydx_hosp),
 
  m2 = logistf(
    hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7 + edu_t6 + pov_t6_std,
    data = df_anydx_hosp),
 
  m3 = logistf(
    hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race,
    data = df_anydx_hosp),
 
  m4 = logistf(
    hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race +
      trct_pct_pov_t6_std + vc_rate15_std,
    data = df_anydx_hosp),
 
  m5 = logistf(
    hosp_asthma_yn ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race +
      trct_pct_pov_t6_std + vc_rate15_std + household_smokers_t6 + household_size_t6,
    data = df_anydx_hosp),
 
  m6_sex = logistf(
    hosp_asthma_yn ~ gvh1600rate_std * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race +
      trct_pct_pov_t6_std + vc_rate15_std + household_smokers_t6 + household_size_t6,
    data = df_anydx_hosp)
)
 
lapply(out_hosp_firth_anydx_list, summary)
 
results_hosp_any <- bind_rows(
  extract_firth(out_hosp_firth_anydx_list$m1,     "M1"),
  extract_firth(out_hosp_firth_anydx_list$m2,     "M2"),
  extract_firth(out_hosp_firth_anydx_list$m3,     "M3"),
  extract_firth(out_hosp_firth_anydx_list$m4,     "M4"),
  extract_firth(out_hosp_firth_anydx_list$m5,     "M5"),
  extract_firth(out_hosp_firth_anydx_list$m6_sex, "M6")
)
results_hosp_any %>% filter(term %in% c("gvh1600rate_std", "gvh1600rate_std:female"))
table_s4 <- make_wide_table(results_hosp_any)
save_table(table_s4, "table_s4_hosp_anydx.docx")
 
 
# =============================================================================
# [21] FUNCTIONAL LIMITATION MODELS (Tables S3, S5) ----
# Linear regression; outcome standardized within each sample. HC1 SEs.
# =============================================================================
 
# --- [21.1] Incident Dx sample (Table S3) ----
df_newdx_sev$limit_asthma_std <- as.numeric(scale(df_newdx_sev$limit_asthma))
 
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
  # M6: sex interaction
  m6_sex = lm(limit_asthma_std ~ gvh1600rate_std * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
                pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                household_smokers_t6 + household_size_t6,
              data = df_newdx_sev)
)
 
lapply(out_sev_new_list, coeftest, vcov. = vcovHC, type = "HC1")
 
results_sev_new <- bind_rows(
  extract_robust(out_sev_new_list$m1,     "M1", transform = identity),
  extract_robust(out_sev_new_list$m2,     "M2", transform = identity),
  extract_robust(out_sev_new_list$m3,     "M3", transform = identity),
  extract_robust(out_sev_new_list$m4,     "M4", transform = identity),
  extract_robust(out_sev_new_list$m5,     "M5", transform = identity),
  extract_robust(out_sev_new_list$m6_sex, "M6", transform = identity)
)
results_sev_new %>% filter(term %in% c("gvh1600rate_std", "gvh1600rate_std:female"))
table_s3 <- make_wide_table(results_sev_new)
save_table(table_s3, "table_s3_limitation_incident.docx")
 
# --- [21.2] Any Dx sample (Table S5) ----
df_anydx_sev$limit_asthma_std <- as.numeric(scale(df_anydx_sev$limit_asthma))
 
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
  # M5: fully adjusted
  m5 = lm(limit_asthma_std ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
            pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + trct_pct_pov_t6_std + vc_rate15_std +
            race + household_smokers_t6 + household_size_t6,
          data = df_anydx_sev),
  # M6: sex interaction
  m6_sex = lm(limit_asthma_std ~ gvh1600rate_std * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
                pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + trct_pct_pov_t6_std + vc_rate15_std +
                race + household_smokers_t6 + household_size_t6,
              data = df_anydx_sev)
)
 
lapply(out_sev_any_list, coeftest, vcov. = vcovHC, type = "HC1")
 
results_sev_any <- bind_rows(
  extract_robust(out_sev_any_list$m1,     "M1", transform = identity),
  extract_robust(out_sev_any_list$m2,     "M2", transform = identity),
  extract_robust(out_sev_any_list$m3,     "M3", transform = identity),
  extract_robust(out_sev_any_list$m4,     "M4", transform = identity),
  extract_robust(out_sev_any_list$m5,     "M5", transform = identity),
  extract_robust(out_sev_any_list$m6_sex, "M6", transform = identity)
)
results_sev_any %>% filter(term %in% c("gvh1600rate_std", "gvh1600rate_std:female"))
table_s5 <- make_wide_table(results_sev_any)
save_table(table_s5, "table_s5_limitation_anydx.docx")

# =============================================================================
# [22] SUPPLEMENTAL: SPATIAL SENSITIVITY ANALYSES (Figure S2) ----
# Exposure entered as the unstandardized annualized rate (OR per additional
# incident per year) at 1,600m, 1,000m, and 500m. M5 covariates. HC1 SEs.
# =============================================================================
 
out_spatial_list <- list(
 
  # --- Main effects ---
  m_1600 = glm(newasthmadx_22 ~ gvh1600_rate + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
                 pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                 household_smokers_t6 + household_size_t6,
               data = df_main_nomiss, family = binomial(link = "logit")),
 
  m_1000 = glm(newasthmadx_22 ~ gvh1000_rate + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
                 pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                 household_smokers_t6 + household_size_t6,
               data = df_main_nomiss, family = binomial(link = "logit")),
 
  m_500  = glm(newasthmadx_22 ~ gvh500_rate + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
                 pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                 household_smokers_t6 + household_size_t6,
               data = df_main_nomiss, family = binomial(link = "logit")),
 
  # --- Sex interaction models ---
  m_1600_sex = glm(newasthmadx_22 ~ gvh1600_rate * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
                     pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                     household_smokers_t6 + household_size_t6,
                   data = df_main_nomiss, family = binomial(link = "logit")),
 
  m_1000_sex = glm(newasthmadx_22 ~ gvh1000_rate * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
                     pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                     household_smokers_t6 + household_size_t6,
                   data = df_main_nomiss, family = binomial(link = "logit")),
 
  m_500_sex  = glm(newasthmadx_22 ~ gvh500_rate * female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
                     pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
                     household_smokers_t6 + household_size_t6,
                   data = df_main_nomiss, family = binomial(link = "logit"))
)
 
lapply(out_spatial_list, coeftest, vcov. = vcovHC, type = "HC1")
 
# Overall, male, and female ORs per additional incident at each buffer
buffers <- c("1,600m" = "gvh1600_rate", "1,000m" = "gvh1000_rate", "500m" = "gvh500_rate")
models_main <- list(out_spatial_list$m_1600,     out_spatial_list$m_1000,     out_spatial_list$m_500)
models_sex  <- list(out_spatial_list$m_1600_sex, out_spatial_list$m_1000_sex, out_spatial_list$m_500_sex)
 
fp <- bind_rows(lapply(seq_along(buffers), function(i) {
  v  <- buffers[[i]]
  ov <- robust_ci_term(models_main[[i]], v)
  ml <- robust_ci_term(models_sex[[i]],  v)
  fm <- unname(female_or(models_sex[[i]], v))   # OR, lower, upper
  tibble(buffer = names(buffers)[i],
         group  = c("Overall", "Males", "Females"),
         OR     = c(ov$estimate, ml$estimate, fm[1]),
         lo     = c(ov$l95,      ml$l95,      fm[2]),
         hi     = c(ov$u95,      ml$u95,      fm[3]))
}))
fp
fp$buffer <- factor(fp$buffer, levels = names(buffers))
fp$group  <- factor(fp$group,  levels = c("Females", "Males", "Overall"))
 
# --- Figure S2: spatial buffer forest plot ----
figS2_spatial <- ggplot(fp, aes(x = OR, y = group)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.15, orientation = "y") +
  geom_point(size = 2) +
  facet_grid(buffer ~ ., switch = "y") +
  scale_x_log10() +
  labs(x = "Odds ratio per additional deadly gun violence incident (log scale)", y = NULL) +
  theme_classic() +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0))
 
print(figS2_spatial)
ggsave(file.path(out_dir, "figS2_spatial_buffers.png"), figS2_spatial,
       width = 6, height = 5, dpi = 500)
 
 
# =============================================================================
# [23] SUPPLEMENTAL: FALSIFICATION TESTS AT AGE 15 ----
# Tests whether DGV at age 15 was associated with asthma at or before age 15.
# Two outcomes: incident dx Y9→Y15 (newasthmadx_15) and prevalent dx at Y15
# (anyasthmadx_15). Separate datasets per outcome to avoid unnecessary case loss.
# HC1 SEs.
# =============================================================================
 
y15_covariates <- c(
  "idnum", "age_t6", "female", "race", "smoking_t6",
  "household_smokers_t6", "household_size_t6",
  "gvh1600count", "gvh1600_rate", "edu_t6", "pov_t6", "trct_pct_pov_t6", "vc_rate15"
)
 
# --- [23.1] Incident asthma Y9→Y15 ----
df_newdx_y15 <- df_reduced %>%
  dplyr::select(all_of(y15_covariates), newasthmadx_15) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(
    gvh1600count_std    = as.numeric(scale(gvh1600count)),
    gvh1600rate_std     = as.numeric(scale(gvh1600_rate)),
    vc_rate15_std       = as.numeric(scale(vc_rate15)),
    pov_t6_std          = as.numeric(scale(pov_t6)),
    trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6))
  )
df_newdx_y15$edu_t6 <- factor(df_newdx_y15$edu_t6)
cat("Incident dx Y9→Y15 n =", nrow(df_newdx_y15), "\n")
 
# --- [23.2] Prevalent asthma at Y15 ----
df_anydx_y15 <- df_reduced %>%
  dplyr::select(all_of(y15_covariates), anyasthmadx_15) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(
    gvh1600count_std    = as.numeric(scale(gvh1600count)),
    gvh1600rate_std     = as.numeric(scale(gvh1600_rate)),
    vc_rate15_std       = as.numeric(scale(vc_rate15)),
    pov_t6_std          = as.numeric(scale(pov_t6)),
    trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6))
  )
df_anydx_y15$edu_t6 <- factor(df_anydx_y15$edu_t6)
cat("Prevalent dx at Y15 n =", nrow(df_anydx_y15), "\n")
 
# --- [23.3] Models ----
m_newdx_y15 <- glm(
  newasthmadx_15 ~ gvh1600rate_std + age_t6 + female + smoking_t6 +
    edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
    household_smokers_t6 + household_size_t6,
  data = df_newdx_y15, family = binomial(link = "logit"))
 
nobs(m_newdx_y15); table(model.frame(m_newdx_y15)[[1]])
coeftest(m_newdx_y15, vcov. = vcovHC, type = "HC1")
exp(coef(m_newdx_y15)["gvh1600rate_std"])
exp(coefci(m_newdx_y15, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])
 
m_anydx_y15 <- glm(
  anyasthmadx_15 ~ gvh1600rate_std + age_t6 + female + smoking_t6 +
    edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
    household_smokers_t6 + household_size_t6,
  data = df_anydx_y15, family = binomial(link = "logit"))
 
nobs(m_anydx_y15); table(model.frame(m_anydx_y15)[[1]])
coeftest(m_anydx_y15, vcov. = vcovHC, type = "HC1")
exp(coef(m_anydx_y15)["gvh1600rate_std"])
exp(coefci(m_anydx_y15, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])
 
 
# =============================================================================
# [24] SUPPLEMENTAL: BUILT ENVIRONMENT SENSITIVITY ----
# Participants with age 15 in-home assessments; adjusts for the 12-item
# interviewer-rated built environment index (Section 11). HC1 SEs.
# =============================================================================
 
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
    gvh1600rate_std     = as.numeric(scale(gvh1600_rate)),
    gvh1000count_std    = as.numeric(scale(gvh1000count)),
    gvh500count_std     = as.numeric(scale(gvh500count)),
    vc_rate15_std       = as.numeric(scale(vc_rate15)),
    pov_t7_std          = as.numeric(scale(pov_t7)),
    pov_t6_std          = as.numeric(scale(pov_t6)),
    trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6))
  )
 
cat("Built environment subsample n =", nrow(df_be), "\n")
cat("Incident asthma cases n =", sum(df_be$newasthmadx_22), "\n")
 
m_be <- glm(
  newasthmadx_22 ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
    pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
    household_smokers_t6 + household_size_t6 + asthma_env_risk,
  data = df_be, family = binomial(link = "logit"))
 
coeftest(m_be, vcov. = vcovHC, type = "HC1")
exp(coef(m_be)["gvh1600rate_std"])
exp(coefci(m_be, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ])
 
 
# =============================================================================
# [25] SUPPLEMENTAL: ADDITIONAL CONFOUNDERS ----
# Maternal prenatal smoking, youth BMI (age 22), and parent asthma diagnosis.
# Separate complete-case datasets per added covariate. HC1 SEs.
# =============================================================================
 
# --- [25.1] Analytic datasets ----
sens_vars <- c("idnum", "age_t7", "female", "race", "smoking_t7", "smoking_t6", "vaping_t7",
               "household_smokers_t6", "household_size_t6", "newasthmadx_22",
               "gvh1600count", "gvh1600_rate", "edu_t7", "edu_t6", "pov_t7", "pov_t6",
               "trct_pct_pov_t6", "vc_rate15")
 
df_sens_momsmk <- df_reduced %>%
  dplyr::select(all_of(sens_vars), momsmk) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(gvh1600count_std = as.numeric(scale(gvh1600count)), gvh1600rate_std = as.numeric(scale(gvh1600_rate)),
         vc_rate15_std = as.numeric(scale(vc_rate15)),
         pov_t7_std = as.numeric(scale(pov_t7)), pov_t6_std = as.numeric(scale(pov_t6)),
         trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6)))
 
df_sens_bmi <- df_reduced %>%
  dplyr::select(all_of(sens_vars), bmiz_t7) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(gvh1600count_std = as.numeric(scale(gvh1600count)), gvh1600rate_std = as.numeric(scale(gvh1600_rate)),
         vc_rate15_std = as.numeric(scale(vc_rate15)),
         pov_t7_std = as.numeric(scale(pov_t7)), pov_t6_std = as.numeric(scale(pov_t6)),
         trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6)))
 
df_sens_pcg <- df_reduced %>%
  dplyr::select(all_of(sens_vars), pcg_asthma) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(gvh1600count_std = as.numeric(scale(gvh1600count)), gvh1600rate_std = as.numeric(scale(gvh1600_rate)),
         vc_rate15_std = as.numeric(scale(vc_rate15)),
         pov_t7_std = as.numeric(scale(pov_t7)), pov_t6_std = as.numeric(scale(pov_t6)),
         trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6)))
 
df_sens_all <- df_reduced %>%
  dplyr::select(all_of(sens_vars), momsmk, bmiz_t7, pcg_asthma) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(gvh1600count_std = as.numeric(scale(gvh1600count)), gvh1600rate_std = as.numeric(scale(gvh1600_rate)),
         vc_rate15_std = as.numeric(scale(vc_rate15)),
         pov_t7_std = as.numeric(scale(pov_t7)), pov_t6_std = as.numeric(scale(pov_t6)),
         trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6)))
 
sapply(list(momsmk = df_sens_momsmk, bmi = df_sens_bmi, pcg = df_sens_pcg, all = df_sens_all), nrow)
 
# --- [25.2] Models ----
out_sens_list <- list(
  m1_momsmk = glm(
    newasthmadx_22 ~ gvh1600rate_std + momsmk + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
      household_smokers_t6 + household_size_t6,
    data = df_sens_momsmk, family = binomial(link = "logit")),
 
  m2_bmi = glm(
    newasthmadx_22 ~ gvh1600rate_std + bmiz_t7 + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
      pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
      household_smokers_t6 + household_size_t6,
    data = df_sens_bmi, family = binomial(link = "logit")),
 
  m3_pcg = glm(
    newasthmadx_22 ~ gvh1600rate_std + pcg_asthma + age_t7 + female + smoking_t6 + smoking_t7 +
      vaping_t7 + pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race + trct_pct_pov_t6_std + vc_rate15_std +
      household_smokers_t6 + household_size_t6,
    data = df_sens_pcg, family = binomial(link = "logit")),
 
  m4_all = glm(
    newasthmadx_22 ~ gvh1600rate_std + momsmk + bmiz_t7 + pcg_asthma + age_t7 + female +
      smoking_t6 + smoking_t7 + vaping_t7 + pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + race +
      trct_pct_pov_t6_std + vc_rate15_std + household_smokers_t6 + household_size_t6,
    data = df_sens_all, family = binomial(link = "logit"))
)
 
lapply(out_sens_list, function(m) {
  list(n     = nobs(m),
       OR    = exp(coef(m)["gvh1600rate_std"]),
       OR_ci = exp(coefci(m, vcov. = vcovHC, type = "HC1")["gvh1600rate_std", ]))
})
 
# Parent asthma OR
exp(c(OR = coef(out_sens_list$m3_pcg)[["pcg_asthma"]],
      coefci(out_sens_list$m3_pcg, vcov. = vcovHC, type = "HC1")["pcg_asthma", ]))
 
 
# =============================================================================
# [26] SUPPLEMENTAL: PROPENSITY SCORE ANALYSES (Figures S3–S4) ----
# Treatment: any deadly gun violence within 1,600m at age 15 (binary).
# Restricted to participants without asthma at Y15. Outcome models use the
# exposure rate standardized among treated (exposed) participants, weighted
# by matching weights.
# =============================================================================
 
# --- [26.1] PSM analytic dataset ----
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
    gvh1600rate_std     = as.numeric(scale(gvh1600_rate)),
    pov_t6_std          = as.numeric(scale(pov_t6)),
    pov_t5_std          = as.numeric(scale(pov_t5)),
    pov_t7_std          = as.numeric(scale(pov_t7)),
    trct_pct_pov_t6_std = as.numeric(scale(trct_pct_pov_t6)),
    trct_pct_pov_t5_std = as.numeric(scale(trct_pct_pov_t5)),
    vc_rate15_std       = as.numeric(scale(vc_rate15)),
    vc_rate9_std        = as.numeric(scale(vc_rate9))
  )
 
cat("PSM pre-match n =", nrow(df_psm), "| exposed =", sum(df_psm$gvh1600_any == 1),
    "| unexposed =", sum(df_psm$gvh1600_any == 0),
    "| cases =", sum(df_psm$newasthmadx_22), "\n")
 
# Weighted mean/SD among treated, used to standardize exposure
standardize_att <- function(data, var_name, treated_idx, weights_var = "weights") {
  var     <- data[[var_name]]
  weights <- data[[weights_var]]
  mean_treated <- Hmisc::wtd.mean(var[treated_idx], weights = weights[treated_idx])
  sd_treated   <- sqrt(Hmisc::wtd.var(var[treated_idx], weights = weights[treated_idx]))
  (var - mean_treated) / sd_treated
}
 
# --- [26.2] Full matching (ATT) ----
# Exact matching on race/ethnicity; calipers of 0.5 SD on poverty variables and
# 1.0 SD on county violent crime at age 9
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
 
summary(m_out_psm, un = FALSE)
 
# Figure S3: propensity score distributions
plot(m_out_psm, type = "jitter", interactive = FALSE)
 
# Figure S4: covariate balance (SMDs)
figS4_balance <- love.plot(m_out_psm, threshold = 0.10)
print(figS4_balance)
ggsave(file.path(out_dir, "figS4_love_plot.png"), figS4_balance,
       width = 7, height = 5, dpi = 500)
 
# --- [26.3] Matched dataset ----
matched_data_psm <- match.data(m_out_psm, subclass = "subclass")
cat("Matched sample n =", nrow(matched_data_psm), "\n")
table(matched_data_psm$newasthmadx_22)
 
matched_data_psm$gvh1600_rate_att_std <- standardize_att(
  matched_data_psm, "gvh1600_rate", matched_data_psm$gvh1600_any == 1
)
summary(matched_data_psm$gvh1600_rate_att_std)
 
# --- [26.4] Weighted regression on matched sample ----
# Cluster-robust SEs clustered on matched subclass
m_psm_weighted <- glm(
  newasthmadx_22 ~ gvh1600_rate_att_std + female + age_t7 + pov_t7_std +
    smoking_t6 + smoking_t7 + vaping_t7 +
    household_smokers_t6 + household_size_t6 + vc_rate15_std,
  family  = quasibinomial(link = "logit"),
  data    = matched_data_psm,
  weights = weights
)
nobs(m_psm_weighted)
coeftest(m_psm_weighted, vcov. = vcovCL, cluster = ~subclass)
exp(coef(m_psm_weighted)["gvh1600_rate_att_std"])
exp(coefci(m_psm_weighted, vcov. = vcovCL, cluster = ~subclass)["gvh1600_rate_att_std", ])
 
m_psm_weighted_int <- glm(
  newasthmadx_22 ~ gvh1600_rate_att_std * female + age_t7 + pov_t7_std +
    smoking_t6 + smoking_t7 + vaping_t7 +
    household_smokers_t6 + household_size_t6 + vc_rate15_std,
  family  = quasibinomial(link = "logit"),
  data    = matched_data_psm,
  weights = weights
)
ct_int <- coeftest(m_psm_weighted_int, vcov. = vcovCL, cluster = ~subclass)
ci_int <- coefci(m_psm_weighted_int,   vcov. = vcovCL, cluster = ~subclass)
cat(sprintf("Interaction OR = %.2f (95%% CI: %.2f-%.2f), p = %.3f\n",
            exp(ct_int["gvh1600_rate_att_std:female", "Estimate"]),
            exp(ci_int["gvh1600_rate_att_std:female", 1]),
            exp(ci_int["gvh1600_rate_att_std:female", 2]),
            ct_int["gvh1600_rate_att_std:female", "Pr(>|z|)"]))
 
# --- [26.5] Propensity score subclassification (10 strata) ----
# Sensitivity check retaining all pre-match participants. HC0 SEs.
m_out_sub <- matchit(
  gvh1600_any ~ race + pov_t6_std + edu_t6 + trct_pct_pov_t6_std +
    trct_pct_pov_t5_std + pov_t5_std + vc_rate9_std,
  data     = df_psm,
  method   = "subclass",
  subclass = 10,
  distance = "glm",
  link     = "logit",
  estimand = "ATT"
)
 
summary(m_out_sub, un = FALSE)                   # pooled balance
summary(m_out_sub, subclass = TRUE, un = FALSE)  # within-subclass balance
 
matched_data_sub <- match.data(m_out_sub, subclass = "subclass")
cat("Subclassified sample n =", nrow(matched_data_sub), "\n")
table(matched_data_sub$newasthmadx_22)
table(matched_data_sub$subclass, matched_data_sub$gvh1600_any)   # controls per stratum
 
matched_data_sub$gvh1600_rate_att_std <- standardize_att(
  matched_data_sub, "gvh1600_rate", matched_data_sub$gvh1600_any == 1
)
 
m_sub_weighted <- glm(
  newasthmadx_22 ~ gvh1600_rate_att_std + female + age_t7 + pov_t7_std +
    smoking_t6 + smoking_t7 + vaping_t7 +
    household_smokers_t6 + household_size_t6 + vc_rate15_std,
  family  = quasibinomial(link = "logit"),
  data    = matched_data_sub,
  weights = weights
)
ct_sub <- coeftest(m_sub_weighted, vcov. = vcovHC, type = "HC0")
ci_sub <- coefci(m_sub_weighted,   vcov. = vcovHC, type = "HC0")
cat(sprintf("OR = %.2f (95%% CI: %.2f-%.2f), p = %.3f\n",
            exp(ct_sub["gvh1600_rate_att_std", "Estimate"]),
            exp(ci_sub["gvh1600_rate_att_std", 1]),
            exp(ci_sub["gvh1600_rate_att_std", 2]),
            ct_sub["gvh1600_rate_att_std", "Pr(>|z|)"]))
 
 
# =============================================================================
# [27] DESCRIPTIVE STATISTICS (Table 1) ----
# Raw (unstandardized) values; overall and by sex. Hospitalization and
# functional limitation are among incident cases.
# =============================================================================
 
mean_sd <- function(x) sprintf("%.1f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
 
n_pct <- function(x) {
  n <- sum(x == 1, na.rm = TRUE); d <- sum(!is.na(x))
  sprintf("%d (%.1f)", n, 100 * n / d)
}
 
n_pct_level <- function(x, level) {
  n <- sum(x == level, na.rm = TRUE); d <- sum(!is.na(x))
  sprintf("%d (%.1f)", n, 100 * n / d)
}
 
df_desc <- df_main_nomiss %>%
  left_join(dplyr::select(df_newdx_hosp, idnum, hosp_asthma_yn), by = "idnum") %>%
  left_join(dplyr::select(df_newdx_sev,  idnum, limit_asthma),   by = "idnum")
 
desc_column <- function(d) {
  c(
    "Incident Asthma (age 15 to 22), N (%)"              = n_pct(d$newasthmadx_22),
    "Hospitalization, Past 12 Months (age 22), N (%)"    = n_pct(d$hosp_asthma_yn),
    "Functional Limitation (age 22), M (SD)"             = mean_sd(d$limit_asthma),
    "Any Gun Violence Exposure (age 15), N (%)"          = n_pct(d$gvh1600_any),
    "Gun Violence Exposure Rate (age 15), M (SD)"        = mean_sd(d$gvh1600_rate),
    "Age at Exposure Assessment (age 15), M (SD)"        = mean_sd(d$age_t6),
    "Age at Follow-Up Assessment (age 22), M (SD)"       = mean_sd(d$age_t7),
    "Smoking, Past Month (age 15), N (%)"                = n_pct(d$smoking_t6),
    "Smoking, Past Month (age 22), N (%)"                = n_pct(d$smoking_t7),
    "Vaping, Past Year (age 22), N (%)"                  = n_pct(d$vaping_t7),
    "Income-Poverty Ratio (age 15), M (SD)"              = mean_sd(d$pov_t6),
    "Income-Poverty Ratio (age 22), M (SD)"              = mean_sd(d$pov_t7),
    "Parental Education (age 15), N (%)"                 = "",
    "  Less than HS"                                     = n_pct_level(d$edu_t6, 1),
    "  HS or Equivalent"                                 = n_pct_level(d$edu_t6, 2),
    "  Some College"                                     = n_pct_level(d$edu_t6, 3),
    "  College/Graduate Degree"                          = n_pct_level(d$edu_t6, 4),
    "Own Education (age 22), N (%)"                      = "",
    "  Less than HS "                                    = n_pct_level(d$edu_t7, 1),
    "  HS or Equivalent "                                = n_pct_level(d$edu_t7, 2),
    "  Some College "                                    = n_pct_level(d$edu_t7, 3),
    "  College/Graduate Degree "                         = n_pct_level(d$edu_t7, 4),
    "Self-Reported Race/Ethnicity, N (%)"                = "",
    "  Non-Hispanic Black/African American"              = n_pct_level(d$race, "Black"),
    "  Non-Hispanic White"                               = n_pct_level(d$race, "White"),
    "  Hispanic"                                         = n_pct_level(d$race, "Hispanic"),
    "  Multiracial or Some Other Race"                   = n_pct_level(d$race, "Other"),
    "Census Tract Poverty % (age 15), M (SD)"            = mean_sd(d$trct_pct_pov_t6),
    "County Violent Crime Rate (age 15), M (SD)"         = mean_sd(d$vc_rate15),
    "Household Smokers (age 15), N (%)"                  = n_pct(d$household_smokers_t6),
    "Household Size (age 15), M (SD)"                    = mean_sd(d$household_size_t6),
    "Parental Asthma Diagnosis, N (%)"                   = n_pct(d$pcg_asthma)
  )
}
 
table1 <- tibble(
  Measure = names(desc_column(df_desc)),
  Overall = desc_column(df_desc),
  Male    = desc_column(df_desc %>% filter(female == 0)),
  Female  = desc_column(df_desc %>% filter(female == 1))
)
print(table1, n = Inf)
 
# Sample sizes by sex
table(df_desc$female)
 
 
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
gvh_term  <- "gvh1600rate_std"
 
# --- Figure 1: overall (M5) and sex-specific (M6) ORs, incident asthma ----
n_all <- nrow(df_main_nomiss)
n_m   <- sum(df_main_nomiss$female == 0)
n_f   <- sum(df_main_nomiss$female == 1)
lab_overall <- sprintf("Overall (N = %s)", format(n_all, big.mark = ","))
lab_male    <- sprintf("Male (n = %s)",    format(n_m,   big.mark = ","))
lab_female  <- sprintf("Female (n = %s)",  format(n_f,   big.mark = ","))
 
plot_tbl_f1 <- bind_rows(
  robust_ci_term(out_main_list$m5,     gvh_term) %>% mutate(group = lab_overall),
  robust_ci_term(out_main_list$m6_sex, gvh_term) %>% mutate(group = lab_male),
  robust_ci_lincomb(out_main_list$m6_sex,
                    setNames(c(1, 1), c(gvh_term, int_term))) %>% mutate(group = lab_female)
) %>%
  mutate(group = factor(group, levels = c(lab_female, lab_male, lab_overall)))
 
fig1_forest <- ggplot(plot_tbl_f1, aes(x = estimate, y = group)) +
  geom_errorbar(aes(xmin = l95, xmax = u95), width = 0, linewidth = 1, alpha = 0.9,
                orientation = "y") +
  geom_point(size = 4) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 1.2) +
  scale_x_continuous(trans = "log10", limits = c(0.8, 2),
                     breaks = c(0.8, 1, 2), labels = c("0.8", "1.0", "2.0")) +
  labs(x = "Odds Ratio (per 1 SD exposure) for Incident Asthma Diagnosis", y = NULL) +
  theme_poster +
  theme(axis.text.y = element_text(size = TICK_TEXT, face = "bold", colour = "black"))
 
# --- Figure 2: predicted probability of incident asthma by sex ----
# Annualized exposure rate (incidents per year) for an interpretable x-axis
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
  geom_rug(data = df_main_nomiss, aes(x = gvh1600_rate), inherit.aes = FALSE,
           sides = "b", alpha = 0.15, length = unit(0.025, "npc")) +
  scale_colour_manual(values = c(Male = "#0072B2", Female = "#D55E00")) +
  scale_fill_manual(values   = c(Male = "#0072B2", Female = "#D55E00")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Annual Exposure Rate at Age 15 (Incidents per Year within 1,600m of Home)",
       y = "Predicted Probability of Incident Asthma",
       colour = NULL, fill = NULL) +
  theme_poster +
  theme(legend.position = c(0.01, 0.99), legend.justification = c(0, 1))
 
# --- Figure 3a: hospitalization ORs (Firth, profile-likelihood CIs) ----
m_inc <- out_hosp_firth_list$m5
m_any <- out_hosp_firth_anydx_list$m5
lab_hosp_inc <- sprintf("Incident Dx\n(n=%d)   ", nrow(df_newdx_hosp))
lab_hosp_any <- sprintf("Any Dx\n(n=%d)",         nrow(df_anydx_hosp))
 
plot_tbl_hosp <- bind_rows(
  tibble(estimate = exp(m_inc$coefficients[gvh_term]),
         l95      = exp(m_inc$ci.lower[gvh_term]),
         u95      = exp(m_inc$ci.upper[gvh_term]),
         group    = lab_hosp_inc),
  tibble(estimate = exp(m_any$coefficients[gvh_term]),
         l95      = exp(m_any$ci.lower[gvh_term]),
         u95      = exp(m_any$ci.upper[gvh_term]),
         group    = lab_hosp_any)
) %>%
  mutate(group = factor(group, levels = c(lab_hosp_any, lab_hosp_inc)))
 
fig3_hosp <- ggplot(plot_tbl_hosp, aes(x = estimate, y = group)) +
  geom_errorbar(aes(xmin = l95, xmax = u95), width = 0, linewidth = 1, alpha = 0.9,
                orientation = "y") +
  geom_point(size = 4) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 1.2) +
  scale_x_continuous(trans = "log10", limits = c(0.8, 30),
                     breaks = c(1, 2, 5, 25), labels = c("1", "2", "5", "25")) +
  labs(x = "Odds Ratio (per 1 SD Exposure) for Asthma Hospitalization", y = NULL) +
  theme_poster +
  theme(axis.text.y = element_text(face = "bold"))
 
# --- Figure 3b: functional limitation betas (HC1) ----
lab_sev_inc <- sprintf("Incident Dx\n(n=%d)   ", nrow(df_newdx_sev))
lab_sev_any <- sprintf("Any Dx\n(n=%d)",         nrow(df_anydx_sev))
 
plot_tbl_sev <- bind_rows(
  robust_ci_term(out_sev_new_list$m5, gvh_term, transform = identity) %>%
    mutate(group = lab_sev_inc),
  robust_ci_term(out_sev_any_list$m5, gvh_term, transform = identity) %>%
    mutate(group = lab_sev_any)
) %>%
  mutate(group = factor(group, levels = c(lab_sev_any, lab_sev_inc)))
 
fig3_sev <- ggplot(plot_tbl_sev, aes(x = estimate, y = group)) +
  geom_errorbar(aes(xmin = l95, xmax = u95), width = 0, linewidth = 1, alpha = 0.9,
                orientation = "y") +
  geom_point(size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 1.2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  labs(x = "β (per 1 SD Exposure) for Functional Limitation", y = NULL) +
  theme_poster +
  theme(axis.text.y = element_text(face = "bold"))
 
# --- Figure 3: combined panels ----
fig3_combined <- (fig3_hosp / fig3_sev) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 16))
 
# --- Print and save ---
print(fig1_forest)
print(fig2_pred)
print(fig3_combined)
 
ggsave(file.path(out_dir, "fig1_forest_overall_sex.png"),  fig1_forest,   width = 7,   height = 5, dpi = 500)
ggsave(file.path(out_dir, "fig2_pred_prob_sex.png"),       fig2_pred,     width = 7.5, height = 6, dpi = 500)
ggsave(file.path(out_dir, "fig3_secondary_combined.png"),  fig3_combined, width = 7,   height = 7, dpi = 500)
 
 
# =============================================================================
# [29] SUPPLEMENTAL: EXPOSURE DATA QUALITY ----
# HC0 robust SEs.
# =============================================================================
 
# --- [29.1] Convergent validity: GVA exposure vs. county violent crime (UCR) ----
cor.test(df_main_nomiss$gvh1600_rate, df_main_nomiss$vc_rate15, method = "spearman")
 
# --- [29.2] Excluding California (FIPS 6) and Michigan (FIPS 26) residents ----
df_sens_state <- df_main_nomiss %>% filter(!state_t6 %in% c(6, 26))
nrow(df_sens_state); sum(df_sens_state$newasthmadx_22)
 
m_sens_state <- glm(
  newasthmadx_22 ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
    edu_t6 + edu_t7 + pov_t6_std + pov_t7_std + race +
    trct_pct_pov_t6_std + vc_rate15_std + household_smokers_t6 + household_size_t6,
  family = binomial(link = "logit"), data = df_sens_state)
ct_state <- coeftest(m_sens_state, vcov. = vcovHC, type = "HC0")
ci_state <- coefci(m_sens_state,   vcov. = vcovHC, type = "HC0")
exp(c(OR = ct_state["gvh1600rate_std", "Estimate"], ci_state["gvh1600rate_std", ]))
ct_state["gvh1600rate_std", "Pr(>|z|)"]
 
# --- [29.3] Restricted to participants with >= 12 months of GVA data ----
# (exposure counts not annualized)
df_sens_full <- df_main_nomiss %>% filter(atleast_1year_exposure == 1)
nrow(df_sens_full); sum(df_sens_full$newasthmadx_22)
 
m_sens_full <- glm(
  newasthmadx_22 ~ gvh1600rate_std + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
    edu_t6 + edu_t7 + pov_t6_std + pov_t7_std + race +
    trct_pct_pov_t6_std + vc_rate15_std + household_smokers_t6 + household_size_t6,
  family = binomial(link = "logit"), data = df_sens_full)
ct_full <- coeftest(m_sens_full, vcov. = vcovHC, type = "HC0")
ci_full <- coefci(m_sens_full,   vcov. = vcovHC, type = "HC0")
exp(c(OR = ct_full["gvh1600rate_std", "Estimate"], ci_full["gvh1600rate_std", ]))
ct_full["gvh1600rate_std", "Pr(>|z|)"]
 
 
# =============================================================================
# [30] SUPPLEMENTAL: EXPOSURE DISTRIBUTION ----
# HC0 robust SEs, except female-specific ORs (female_or, HC1).
# =============================================================================
 
# --- [30.1] Categorical exposure: none, low (<= median among exposed), high ----
mean(df_main_nomiss$gvh1600_rate == 0)
med_exp <- median(df_main_nomiss$gvh1600_rate[df_main_nomiss$gvh1600_rate > 0])
med_exp
 
df_main_nomiss <- df_main_nomiss %>%
  mutate(gv_cat = case_when(gvh1600_rate == 0       ~ "None",
                            gvh1600_rate <= med_exp ~ "Low",
                            TRUE                    ~ "High"),
         gv_cat = factor(gv_cat, levels = c("None", "Low", "High")),
         gv_cat_trend = as.numeric(gv_cat) - 1)
 
table(df_main_nomiss$gv_cat, df_main_nomiss$newasthmadx_22)
range(df_main_nomiss$gvh1600_rate[df_main_nomiss$gv_cat == "Low"])
range(df_main_nomiss$gvh1600_rate[df_main_nomiss$gv_cat == "High"])
 
m_cat <- glm(
  newasthmadx_22 ~ gv_cat + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
    edu_t6 + edu_t7 + pov_t6_std + pov_t7_std + race +
    trct_pct_pov_t6_std + vc_rate15_std + household_smokers_t6 + household_size_t6,
  family = binomial(link = "logit"), data = df_main_nomiss)
ct_cat <- coeftest(m_cat, vcov. = vcovHC, type = "HC0")
ci_cat <- coefci(m_cat,   vcov. = vcovHC, type = "HC0")
round(exp(cbind(OR = coef(m_cat), ci_cat))[c("gv_catLow", "gv_catHigh"), ], 2)
ct_cat[c("gv_catLow", "gv_catHigh"), ]
 
# Linear trend across categories
m_trend <- update(m_cat, . ~ . - gv_cat + gv_cat_trend)
coeftest(m_trend, vcov. = vcovHC, type = "HC0")["gv_cat_trend", ]
 
# --- [30.2] Log-transformed exposure, log(rate + 1) ----
m_log_main <- glm(
  newasthmadx_22 ~ log1p(gvh1600_rate) + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
    pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + trct_pct_pov_t6_std + vc_rate15_std +
    race + household_smokers_t6 + household_size_t6,
  data = df_main_nomiss, family = binomial)
ci_log <- coefci(m_log_main, vcov. = vcovHC, type = "HC0")["log1p(gvh1600_rate)", ]
exp(c(OR = coef(m_log_main)[["log1p(gvh1600_rate)"]], ci_log))
 
m_log_sex <- update(new_dx_1600_int, . ~ . - gvh1600_rate * female + log1p(gvh1600_rate) * female)
coeftest(m_log_sex, vcov. = vcovHC, type = "HC0")[c("log1p(gvh1600_rate)", "log1p(gvh1600_rate):female"), ]
 
# --- [30.3] Excluding participants above the 99th and 95th exposure percentiles ----
# OR per additional incident per year
for (p in c(0.99, 0.95)) {
  cut  <- quantile(df_main_nomiss$gvh1600_rate, p)
  dsub <- subset(df_main_nomiss, gvh1600_rate <= cut)
 
  m_main <- glm(newasthmadx_22 ~ gvh1600_rate + female + age_t7 + smoking_t6 + smoking_t7 + vaping_t7 +
                  pov_t7_std + edu_t7 + edu_t6 + pov_t6_std + trct_pct_pov_t6_std + vc_rate15_std +
                  race + household_smokers_t6 + household_size_t6,
                data = dsub, family = binomial)
  m_int  <- update(new_dx_1600_int, data = dsub)
 
  ci_main <- coefci(m_main, vcov. = vcovHC, type = "HC0")["gvh1600_rate", ]
  ci_male <- coefci(m_int,  vcov. = vcovHC, type = "HC0")["gvh1600_rate", ]
  ct_int  <- coeftest(m_int, vcov. = vcovHC, type = "HC0")
 
  cat(sprintf("\n<= %.0fth percentile (rate <= %.1f): n = %d, cases = %d\n",
              p * 100, cut, nobs(m_main), sum(dsub$newasthmadx_22)))
  cat(sprintf("  Overall OR = %.2f (%.2f-%.2f)\n",
              exp(coef(m_main)["gvh1600_rate"]), exp(ci_main[1]), exp(ci_main[2])))
  cat(sprintf("  Males OR   = %.2f (%.2f-%.2f)\n",
              exp(coef(m_int)["gvh1600_rate"]), exp(ci_male[1]), exp(ci_male[2])))
  cat("  Females OR:\n"); print(female_or(m_int, "gvh1600_rate"))
  cat(sprintf("  Interaction p = %.3f\n", ct_int["gvh1600_rate:female", 4]))
}
 
 
# =============================================================================
# [31] SUPPLEMENTAL: NEIGHBORHOOD TRAFFIC DENSITY ----
# Youth-reported traffic density at age 15 (proxy for traffic-related air
# pollution). HC0 robust SEs.
# =============================================================================
 
table(df_main_nomiss$traffic_t6, useNA = "ifany")
 
df_traffic <- df_main_nomiss %>% filter(!is.na(traffic_t6))
nrow(df_traffic); sum(df_traffic$newasthmadx_22)
 
m_traffic <- glm(
  newasthmadx_22 ~ gvh1600rate_std + traffic_t6 + age_t7 + female + smoking_t6 + smoking_t7 + vaping_t7 +
    edu_t6 + edu_t7 + pov_t6_std + pov_t7_std + race +
    trct_pct_pov_t6_std + vc_rate15_std + household_smokers_t6 + household_size_t6,
  family = binomial(link = "logit"), data = df_traffic)
ct_tr <- coeftest(m_traffic, vcov. = vcovHC, type = "HC0")
ci_tr <- coefci(m_traffic,   vcov. = vcovHC, type = "HC0")
exp(c(OR = ct_tr["gvh1600rate_std", "Estimate"], ci_tr["gvh1600rate_std", ]))
 
# Exposure x traffic interaction
m_traffic_int <- update(m_traffic, . ~ . + gvh1600rate_std:traffic_t6)
coeftest(m_traffic_int, vcov. = vcovHC, type = "HC0")["gvh1600rate_std:traffic_t6", ]
