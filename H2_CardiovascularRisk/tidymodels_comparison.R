### CoDiet CVD Prediction Model Comparison
### Author: Luisa Delius

### Maybe add an explaination here of what I've been doing

# ============================================
# 1.1 LOAD PACKAGES
# ============================================
library(tidyverse)
library(tidymodels)
library(vip)
library(future)
library(janitor)
library(tidytext)
library(qs2)
library(plsmod)
library(ggplot2)
library(writexl)
library(flextable)
library(officer)
library(patchwork)

# ============================================
# 1.2 GLOBAL SETTINGS
# ============================================
set.seed(42)
options(scipen = 999, expressions = 500000)

# Parallel processing
CPUS <- parallel::detectCores() - 1
cat("Using", CPUS, "CPU cores for parallel processing\n")

# Resolve common function conflicts
select <- dplyr::select
slice <- dplyr::slice
rename <- dplyr::rename
filter <- dplyr::filter

wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files/processed_data"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

# ============================================
# 2. Define CVD Scores and Data Sets
# ============================================
# ---- 2.1 Load raw datasets ----
df_fatty_acids_statin_suppl <- readRDS("df_fatty_acids_predictor_statin_suppl.rds") 

df_fatty_acids <- df_fatty_acids_statin_suppl %>%
  select(-starts_with("z_"), -QRISK3_risk, -Statins, -Supplements)

df_lipidomics <- readRDS("df_lipidomics_predictor_statin_suppl.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk, -Statins, -Supplements) %>%
  filter(!is.na(Sample_ID))

df_urine_nmr <- readRDS("df_urine_NMR_data.rds") %>%
  select(-starts_with("z_"))

df_body_composition <- readRDS("df_body_composition_metrics.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk) %>%
  filter(!is.na(Sample_ID))

df_REDcap_demographics <- readRDS("df_REDcap_demographics_ElaNet.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk, -recruitment_site)

df_risk_factors <- readRDS("df_risk_factor_predictors.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk,
         -Heart.Rate, -Age, -Total.Cholesterol.mg.dl, -HDL.mg.dl,
         -Body.Weight, -Height, -BMI, -Gender,
         -stress_resilience_status, -stress_index_status, -Age.Risk) %>%
  full_join(df_fatty_acids_statin_suppl %>% select(Sample_ID, Statins, Supplements), by = "Sample_ID")

df_all_cvd_risk_scores <- readRDS("df_all_risk_scores.rds") %>%
  rename(Sample_ID = PatientID) %>%
  select(-SCORE2_strat)

df_ascvd_frs_score2_input <- readRDS("ASCVD_SCORE2_Framingham_input.rds") %>%
  rename(Sample_ID = PatientID)

df_QRISK3_input <- readRDS("QRISK3_calculation_input.rds") %>%
  rename(Sample_ID = PatientID)

CVD_scores <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y", "mean_risk")

# ---- 2.2 Create datasets that shall be tested ----
df_fatty_acids_model <- df_fatty_acids %>%
  full_join(df_all_cvd_risk_scores, by = "Sample_ID")

df_lipidomics_model <- df_lipidomics %>%
  full_join(df_all_cvd_risk_scores, by = "Sample_ID")

df_urine_nmr_model <- df_urine_nmr %>%
  full_join(df_all_cvd_risk_scores, by = "Sample_ID")

df_body_composition_model <- df_body_composition %>%
  full_join(df_all_cvd_risk_scores, by = "Sample_ID")

df_REDcap_demographics_model <- df_REDcap_demographics %>%
  full_join(df_all_cvd_risk_scores, by = "Sample_ID")

df_risk_factors_model <- df_risk_factors %>%
  full_join(df_all_cvd_risk_scores, by = "Sample_ID")

df_all_data_model <- df_all_cvd_risk_scores %>%
  full_join(df_fatty_acids, by = "Sample_ID") %>%
  full_join(df_lipidomics, by = "Sample_ID") %>%
  full_join(df_urine_nmr, by = "Sample_ID") %>%
  full_join(df_body_composition, by = "Sample_ID") %>%
  full_join(df_REDcap_demographics, by = "Sample_ID") %>%
  full_join(df_risk_factors, by = "Sample_ID")

# # Significant predictors
# BH_sig_predictors <- c("ecw_tbw", "mean_hrt", "ecm_bcm", "rbc_22_4n_6",
#                        "naps_during_day", "Living_Status")
# 
# df_BH_sig_predictors <- df_all_data_model %>%
#   select(Sample_ID, all_of(BH_sig_predictors), all_of(CVD_scores))
# 
# sig_predictors <- c("AGE.reader", "Trigonelline", "Hba1C", "ALT.unit.L",
#                     "rbc_dpa_22_5n3", "rbc_eicosadienoic_20_2n6",
#                     "rbc_epa_20_5n3", "pe_o_19_1_20_5",
#                     "3_hydroxybutyric_acid", "hippuric_acid",
#                     "Employment_Status", "Education_Level")
# 
# df_sig_predictors <- df_all_data_model %>%
#   select(Sample_ID, all_of(BH_sig_predictors), all_of(sig_predictors), all_of(CVD_scores))

datasets_list <- list(
  all_data = df_all_data_model,
  fatty_acids = df_fatty_acids_model,
  lipids = df_lipidomics_model,
  urine_NMR = df_urine_nmr_model,
  body_composition = df_body_composition_model,
  sociodemographics_lifestyle = df_REDcap_demographics_model,
  clinical_risk_factors = df_risk_factors_model,
  # mtc_sign = df_BH_sig_predictors,
  # sign_predictors = df_sig_predictors
)

# Score-specific input datasets
df_model0_QRISK3 <- df_QRISK3_input %>%
  mutate(townsend = replace_na(townsend, 0)) %>%
  full_join(df_all_cvd_risk_scores %>% select(Sample_ID, QRISK3_risk), by = "Sample_ID") %>%
  filter(!is.na(QRISK3_risk))

df_model0_ascvd <- df_ascvd_frs_score2_input %>%
  select(-Risk.region, -mean_LDL_mg_dl) %>%
  full_join(df_all_cvd_risk_scores %>% select(Sample_ID, ascvd_10y), by = "Sample_ID") %>%
  filter(!is.na(ascvd_10y))

df_model0_score2 <- df_ascvd_frs_score2_input %>%
  select(-mean_LDL_mg_dl, -blood_pressure_treatment, -race_ascvd) %>%
  full_join(df_all_cvd_risk_scores %>% select(Sample_ID, SCORE2_score), by = "Sample_ID") %>%
  filter(!is.na(SCORE2_score))

df_model0_frs <- df_ascvd_frs_score2_input %>%
  select(-Risk.region, -mean_LDL_mg_dl, -race_ascvd) %>%
  full_join(df_all_cvd_risk_scores %>% select(Sample_ID, frs_10y), by = "Sample_ID") %>%
  filter(!is.na(frs_10y))

df_model0_composite <- df_QRISK3_input %>%
  mutate(townsend = replace_na(townsend, 0)) %>%
  full_join(df_ascvd_frs_score2_input %>% select(Sample_ID, Risk.region), by = "Sample_ID") %>%
  full_join(df_all_cvd_risk_scores %>% select(Sample_ID, mean_risk), by = "Sample_ID") %>%
  filter(!is.na(mean_risk))

score_specific_datasets <- list(
  QRISK3_risk = df_model0_QRISK3,
  SCORE2_score = df_model0_score2,
  frs_10y = df_model0_frs,
  ascvd_10y = df_model0_ascvd,
  mean_risk = df_model0_composite
)

# # Body Composition Extended
# height_weight_vars <- df_QRISK3_input %>%
#   select(Sample_ID, Height_cm, Weight_kg)
# 
# df_body_comp_extended <- df_body_composition_model %>%
#   full_join(height_weight_vars, by = "Sample_ID")
# 
# scores_without_ht_wt <- c("SCORE2_score", "ascvd_10y", "frs_10y")

# ============================================
# 3. HELPER FUNCTIONS
# ============================================

prep_extra_predictors <- function(extra_df) {
  keep_cols <- setdiff(names(extra_df), CVD_scores)
  extra_df[, keep_cols, drop = FALSE]
}

build_score_plus_block <- function(score, block_name, score_specific_datasets, datasets_list) {
  base_df <- score_specific_datasets[[score]]
  extra_df <- datasets_list[[block_name]]
  extra_predictors <- prep_extra_predictors(extra_df)
  dplyr::left_join(base_df, extra_predictors, by = "Sample_ID")
}

# # ============================================
# # 4. EXCLUDE PARTICIPANTS WITH >80% MISSING
# # ============================================
# 
# cat("\n=== EXCLUDING HIGH-MISSINGNESS PARTICIPANTS ===\n")
# 
# # Calculate row-level missingness
# predictor_cols <- setdiff(names(df_all_data_model), c("Sample_ID", CVD_scores))
# row_missing_pct <- rowMeans(is.na(df_all_data_model[predictor_cols])) * 100
# 
# df_all_data_model$row_missing_pct <- row_missing_pct
# 
# # Identify participants with >80% missing
# EXCLUDE_HIGH_MISSING <- df_all_data_model %>%
#   filter(row_missing_pct > 80) %>%
#   pull(Sample_ID)
# 
# cat(sprintf("Participants with >80%% missing: %d\n", length(EXCLUDE_HIGH_MISSING)))
# cat("IDs:", paste(EXCLUDE_HIGH_MISSING, collapse = ", "), "\n")
# 
# # Remove from df_all_data_model (will propagate to datasets_list later)
# df_all_data_model <- df_all_data_model %>%
#   filter(!Sample_ID %in% EXCLUDE_HIGH_MISSING) %>%
#   select(-row_missing_pct)
# 
#qs2::qs_save(EXCLUDE_HIGH_MISSING, "excluded_high_missing_ids.qs2")
# cat("✓ Excluded participant IDs saved\n")
EXCLUDE_HIGH_MISSING <- character(0)

# ============================================
# 5. SELECT 20 HOLD-OUT PARTICIPANTS
# ============================================

cat("\n=== SELECTING HOLD-OUT PARTICIPANTS ===\n")

# Build participant info table
participant_info <- df_all_data_model %>%
  select(Sample_ID) %>%
  mutate(
    missing_pct = rowMeans(is.na(df_all_data_model %>% select(-Sample_ID, -all_of(CVD_scores)))) * 100
  ) %>%
  left_join(
    df_all_cvd_risk_scores %>% select(Sample_ID, mean_risk),
    by = "Sample_ID"
  ) %>%
  left_join(
    df_QRISK3_input %>% select(Sample_ID, Sex),
    by = "Sample_ID"
  ) %>%
  left_join(
    df_risk_factors %>% select(Sample_ID, Country),
    by = "Sample_ID"
  ) %>%
  filter(
    !is.na(mean_risk),
    !is.na(Sex),
    !is.na(Country)
  )

cat(sprintf("Total eligible participants: %d\n", nrow(participant_info)))

# Cohort distribution
cat("\n--- Full Cohort Distribution ---\n")

cat("\nBy Country:\n")
country_dist <- table(participant_info$Country)
print(country_dist)
print(round(prop.table(country_dist) * 100, 1))

cat("\nBy Sex:\n")
sex_dist <- table(participant_info$Sex)
print(sex_dist)
print(round(prop.table(sex_dist) * 100, 1))

# Create risk tertiles based on mean_risk
participant_info <- participant_info %>%
  mutate(
    risk_tertile = ntile(mean_risk, 3),
    risk_category = factor(
      case_when(
        risk_tertile == 1 ~ "Low",
        risk_tertile == 2 ~ "Medium",
        risk_tertile == 3 ~ "High"
      ),
      levels = c("Low", "Medium", "High")
    )
  )

cat("\nBy Risk Category (mean_risk tertiles):\n")
risk_dist <- table(participant_info$risk_category)
print(risk_dist)

cat("\nmean_risk ranges per tertile:\n")
participant_info %>%
  group_by(risk_category) %>%
  summarise(
    min = round(min(mean_risk), 2),
    max = round(max(mean_risk), 2),
    .groups = "drop"
  ) %>%
  print()

# ----- STRATIFIED SELECTION -----
set.seed(42)
n_holdout <- 20
min_per_country <- 2

cat(sprintf("\n--- Selection Strategy ---\n"))
cat(sprintf("Target: %d participants\n", n_holdout))
cat(sprintf("Minimum %d per country\n", min_per_country))

# Prioritize low missingness
holdout_candidates <- participant_info %>%
  arrange(missing_pct)

# Step 1: Select 2 per country (lowest missingness, balanced by sex if possible)
holdout_ids <- c()

for (ctry in names(country_dist)) {
  country_candidates <- holdout_candidates %>%
    filter(Country == ctry, !Sample_ID %in% holdout_ids)
  
  # Try to get 1 male and 1 female if possible
  males <- country_candidates %>% filter(Sex == "Male") %>% slice_head(n = 1)
  females <- country_candidates %>% filter(Sex == "Female") %>% slice_head(n = 1)
  
  if (nrow(males) > 0 && nrow(females) > 0) {
    selected <- bind_rows(males, females)
  } else {
    selected <- country_candidates %>% slice_head(n = 2)
  }
  
  holdout_ids <- c(holdout_ids, selected$Sample_ID)
}

cat(sprintf("\nAfter 2 per country: %d selected\n", length(holdout_ids)))

# Step 2: Fill remaining slots to match sex and risk distribution
remaining_slots <- n_holdout - length(holdout_ids)

if (remaining_slots > 0) {
  remaining_candidates <- holdout_candidates %>%
    filter(!Sample_ID %in% holdout_ids)
  
  current_selection <- participant_info %>% filter(Sample_ID %in% holdout_ids)
  current_sex <- table(current_selection$Sex)
  current_risk <- table(current_selection$risk_category)
  
  target_sex <- round(prop.table(sex_dist) * n_holdout)
  target_risk <- round(prop.table(risk_dist) * n_holdout)
  
  need_sex <- pmax(0, target_sex - current_sex)
  need_risk <- pmax(0, target_risk - current_risk)
  
  cat("\nCurrent vs Target:\n")
  cat("Sex - Current:", paste(names(current_sex), current_sex, collapse = ", "), "\n")
  cat("Sex - Target:", paste(names(target_sex), target_sex, collapse = ", "), "\n")
  cat("Risk - Current:", paste(names(current_risk), current_risk, collapse = ", "), "\n")
  cat("Risk - Target:", paste(names(target_risk), target_risk, collapse = ", "), "\n")
  
  additional <- remaining_candidates %>%
    mutate(
      sex_need = case_when(
        Sex == "Male" & need_sex["Male"] > 0 ~ 2,
        Sex == "Female" & need_sex["Female"] > 0 ~ 2,
        TRUE ~ 0
      ),
      risk_need = case_when(
        risk_category == "Low" & need_risk["Low"] > 0 ~ 1,
        risk_category == "Medium" & need_risk["Medium"] > 0 ~ 1,
        risk_category == "High" & need_risk["High"] > 0 ~ 1,
        TRUE ~ 0
      ),
      priority_score = sex_need + risk_need
    ) %>%
    arrange(desc(priority_score), missing_pct) %>%
    slice_head(n = remaining_slots)
  
  holdout_ids <- c(holdout_ids, additional$Sample_ID)
}

# Final hold-out set
holdout_participants <- participant_info %>%
  filter(Sample_ID %in% holdout_ids)

# ----- VERIFY DISTRIBUTION -----
cat(sprintf("\n=== FINAL HOLD-OUT SET: %d participants ===\n", nrow(holdout_participants)))

cat("\nHold-out vs Cohort - Country:\n")
comparison_country <- data.frame(
  Country = names(country_dist),
  Cohort_N = as.numeric(country_dist),
  Cohort_Pct = round(as.numeric(prop.table(country_dist)) * 100, 1),
  Holdout_N = sapply(names(country_dist), function(x) sum(holdout_participants$Country == x)),
  stringsAsFactors = FALSE
)
comparison_country$Holdout_Pct <- round(comparison_country$Holdout_N / n_holdout * 100, 1)
print(comparison_country)

cat("\nHold-out vs Cohort - Sex:\n")
comparison_sex <- data.frame(
  Sex = names(sex_dist),
  Cohort_N = as.numeric(sex_dist),
  Cohort_Pct = round(as.numeric(prop.table(sex_dist)) * 100, 1),
  Holdout_N = sapply(names(sex_dist), function(x) sum(holdout_participants$Sex == x)),
  stringsAsFactors = FALSE
)
comparison_sex$Holdout_Pct <- round(comparison_sex$Holdout_N / n_holdout * 100, 1)
print(comparison_sex)

cat("\nHold-out vs Cohort - Risk (mean_risk tertiles):\n")
comparison_risk <- data.frame(
  Risk = c("Low", "Medium", "High"),
  Cohort_N = as.numeric(risk_dist),
  Cohort_Pct = round(as.numeric(prop.table(risk_dist)) * 100, 1),
  Holdout_N = sapply(c("Low", "Medium", "High"), function(x) sum(holdout_participants$risk_category == x)),
  stringsAsFactors = FALSE
)
comparison_risk$Holdout_Pct <- round(comparison_risk$Holdout_N / n_holdout * 100, 1)
print(comparison_risk)

cat("\nMissingness in Hold-out:\n")
cat(sprintf("  Mean: %.1f%%\n", mean(holdout_participants$missing_pct)))
cat(sprintf("  Max:  %.1f%%\n", max(holdout_participants$missing_pct)))

# Save hold-out IDs
HOLDOUT_IDS <- holdout_participants$Sample_ID
qs2::qs_save(HOLDOUT_IDS, "holdout_participant_ids.qs2")
cat("\n✓ Hold-out IDs saved to holdout_participant_ids.qs2\n")

cat("\nSelected participants:\n")
print(holdout_participants %>% 
        select(Sample_ID, Country, Sex, risk_category, missing_pct, mean_risk) %>%
        arrange(Country, Sex) %>%
        as.data.frame())

# ============================================
# 6. APPLY 40% COLUMN MISSINGNESS EXCLUSION
# ============================================

MISSINGNESS_THRESHOLD <- 0.40

cat("\n=== APPLYING 40% COLUMN MISSINGNESS EXCLUSION ===\n")

remove_high_missing <- function(df, threshold = MISSINGNESS_THRESHOLD,
                                exclude_cols = c("Sample_ID", CVD_scores)) {
  cols_to_check <- setdiff(names(df), exclude_cols)
  if (length(cols_to_check) == 0) return(df)
  
  missing_pct <- sapply(df[cols_to_check], function(x) mean(is.na(x)))
  cols_to_remove <- names(missing_pct)[missing_pct >= threshold]
  
  if (length(cols_to_remove) > 0) {
    df <- df %>% select(-all_of(cols_to_remove))
  }
  df
}

# Show what will be removed
predictor_cols <- setdiff(names(df_all_data_model), c("Sample_ID", CVD_scores))
missing_pct <- sapply(df_all_data_model[predictor_cols], function(x) mean(is.na(x)) * 100)
vars_excluded <- names(missing_pct)[missing_pct >= MISSINGNESS_THRESHOLD * 100]

if (length(vars_excluded) > 0) {
  cat("Variables excluded:\n")
  for (v in vars_excluded) {
    cat(sprintf("  %.1f%% - %s\n", missing_pct[v], v))
  }
}

cat("\nNote: MAR assumption supported for all variables.\n")

# ============================================
# 7. REBUILD DATASETS WITH EXCLUSIONS
# ============================================
# KEY: apply_exclusions() ONLY removes high-missing participants
# Holdout participants are KEPT - the modeling function handles the split

cat("\n=== REBUILDING DATASETS WITH EXCLUSIONS ===\n")

# ONLY remove high-missing participants, NOT holdout
apply_exclusions <- function(df, high_missing_ids = EXCLUDE_HIGH_MISSING) {
  df %>%
    filter(!Sample_ID %in% high_missing_ids)
}

# Rebuild individual datasets (INCLUDING holdout participants)
df_fatty_acids_model <- df_fatty_acids %>%
  full_join(df_all_cvd_risk_scores, by = "Sample_ID") %>%
  apply_exclusions()

df_lipidomics_model <- df_lipidomics %>%
  full_join(df_all_cvd_risk_scores, by = "Sample_ID") %>%
  apply_exclusions()

df_urine_nmr_model <- df_urine_nmr %>%
  full_join(df_all_cvd_risk_scores, by = "Sample_ID") %>%
  apply_exclusions()

df_body_composition_model <- df_body_composition %>%
  full_join(df_all_cvd_risk_scores, by = "Sample_ID") %>%
  apply_exclusions()

df_REDcap_demographics_model <- df_REDcap_demographics %>%
  full_join(df_all_cvd_risk_scores, by = "Sample_ID") %>%
  apply_exclusions()

df_risk_factors_model <- df_risk_factors %>%
  full_join(df_all_cvd_risk_scores, by = "Sample_ID") %>%
  apply_exclusions()

df_all_data_model <- df_all_cvd_risk_scores %>%
  full_join(df_fatty_acids, by = "Sample_ID") %>%
  full_join(df_lipidomics, by = "Sample_ID") %>%
  full_join(df_urine_nmr, by = "Sample_ID") %>%
  full_join(df_body_composition, by = "Sample_ID") %>%
  full_join(df_REDcap_demographics, by = "Sample_ID") %>%
  full_join(df_risk_factors, by = "Sample_ID") %>%
  apply_exclusions()

df_BH_sig_predictors <- df_all_data_model %>%
  select(Sample_ID, all_of(BH_sig_predictors), all_of(CVD_scores))

df_sig_predictors <- df_all_data_model %>%
  select(Sample_ID, any_of(c(BH_sig_predictors, sig_predictors)), all_of(CVD_scores))

# Rebuild datasets_list (all contain holdout participants)
datasets_list <- list(
  all_data = df_all_data_model,
  fatty_acids = df_fatty_acids_model,
  lipids = df_lipidomics_model,
  urine_NMR = df_urine_nmr_model,
  body_composition = df_body_composition_model,
  sociodemographics_lifestyle = df_REDcap_demographics_model,
  clinical_risk_factors = df_risk_factors_model,
  mtc_sign = df_BH_sig_predictors,
  sign_predictors = df_sig_predictors
)

# Apply column missingness exclusion
datasets_list <- lapply(datasets_list, remove_high_missing)

# Rebuild score-specific datasets (also INCLUDING holdout)
df_model0_QRISK3 <- df_QRISK3_input %>%
  mutate(townsend = replace_na(townsend, 0)) %>%
  full_join(df_all_cvd_risk_scores %>% select(Sample_ID, QRISK3_risk), by = "Sample_ID") %>%
  filter(!is.na(QRISK3_risk)) %>%
  apply_exclusions()

df_model0_ascvd <- df_ascvd_frs_score2_input %>%
  select(-Risk.region, -mean_LDL_mg_dl) %>%
  full_join(df_all_cvd_risk_scores %>% select(Sample_ID, ascvd_10y), by = "Sample_ID") %>%
  filter(!is.na(ascvd_10y)) %>%
  apply_exclusions()

df_model0_score2 <- df_ascvd_frs_score2_input %>%
  select(-mean_LDL_mg_dl, -blood_pressure_treatment, -race_ascvd) %>%
  full_join(df_all_cvd_risk_scores %>% select(Sample_ID, SCORE2_score), by = "Sample_ID") %>%
  filter(!is.na(SCORE2_score)) %>%
  apply_exclusions()

df_model0_frs <- df_ascvd_frs_score2_input %>%
  select(-Risk.region, -mean_LDL_mg_dl, -race_ascvd) %>%
  full_join(df_all_cvd_risk_scores %>% select(Sample_ID, frs_10y), by = "Sample_ID") %>%
  filter(!is.na(frs_10y)) %>%
  apply_exclusions()

df_model0_composite <- df_QRISK3_input %>%
  mutate(townsend = replace_na(townsend, 0)) %>%
  full_join(df_ascvd_frs_score2_input %>% select(Sample_ID, Risk.region), by = "Sample_ID") %>%
  full_join(df_all_cvd_risk_scores %>% select(Sample_ID, mean_risk), by = "Sample_ID") %>%
  filter(!is.na(mean_risk)) %>%
  apply_exclusions()

score_specific_datasets <- list(
  QRISK3_risk = df_model0_QRISK3,
  SCORE2_score = df_model0_score2,
  frs_10y = df_model0_frs,
  ascvd_10y = df_model0_ascvd,
  mean_risk = df_model0_composite
)

score_specific_datasets <- lapply(score_specific_datasets, remove_high_missing)

# Body comp extended
df_body_comp_extended <- df_body_composition %>%
  full_join(df_all_cvd_risk_scores, by = "Sample_ID") %>%
  full_join(height_weight_vars, by = "Sample_ID") %>%
  apply_exclusions() %>%
  remove_high_missing()

# Verify holdout participants are present
cat("\n--- Verification: Holdout participants in datasets ---\n")
for (ds_name in names(datasets_list)) {
  n_holdout_in_ds <- sum(HOLDOUT_IDS %in% datasets_list[[ds_name]]$Sample_ID)
  cat(sprintf("  %s: %d/%d holdout participants present\n", 
              ds_name, n_holdout_in_ds, length(HOLDOUT_IDS)))
}

n_train_expected <- nrow(df_all_data_model) - length(HOLDOUT_IDS)
cat(sprintf("\n✓ Total participants: %d\n", nrow(df_all_data_model)))
cat(sprintf("✓ Training set size (expected): %d participants\n", n_train_expected))
cat(sprintf("✓ Hold-out set size: %d participants\n", length(HOLDOUT_IDS)))

# ============================================
# 9. DATA SANITIZATION
# ============================================

sanitize_df <- function(df) {
  for (col in names(df)) {
    if (col == "Sample_ID") next
    if (col %in% CVD_scores) next
    
    x <- df[[col]]
    
    if (is.logical(x)) {
      df[[col]] <- as.integer(x)
    } else if (inherits(x, "Date") || inherits(x, "POSIXct") || inherits(x, "POSIXlt")) {
      df[[col]] <- as.numeric(x)
    } else if (!is.numeric(x) && !is.factor(x)) {
      df[[col]] <- as.factor(x)
    }
  }
  df
}

cat("\n=== SANITIZING DATASETS ===\n")
datasets_list <- lapply(datasets_list, sanitize_df)
score_specific_datasets <- lapply(score_specific_datasets, sanitize_df)
df_body_comp_extended <- sanitize_df(df_body_comp_extended)
cat("✓ All datasets sanitized\n")

# ============================================
# 10. CREATE SCORE-SPECIFIC + BLOCK COMBINATIONS
# ============================================

cat("\n=== CREATING SCORE-SPECIFIC + BLOCK COMBINATIONS ===\n")

score_plus_block_datasets <- list()

for (score in names(score_specific_datasets)) {
  for (block in names(datasets_list)) {
    combo_name <- paste0(score, "_plus_", block)
    score_plus_block_datasets[[combo_name]] <- build_score_plus_block(
      score = score,
      block_name = block,
      score_specific_datasets = score_specific_datasets,
      datasets_list = datasets_list
    )
  }
}

cat(sprintf("✓ Created %d score-specific + block combinations\n", length(score_plus_block_datasets)))

# ============================================
# 12. CV AND TUNING SETTINGS
# ============================================

CV_FOLDS <- 5
CV_REPEATS <- 5
GRID_SIZE <- 20
N_IMPORTANCE_SIMS <- 3

cat("\n=== SETTINGS ===\n")
cat(sprintf("CV: %d-fold × %d repeats = %d resamples\n", CV_FOLDS, CV_REPEATS, CV_FOLDS * CV_REPEATS))
cat(sprintf("Grid size: %d random combinations per model\n", GRID_SIZE))
cat(sprintf("Importance simulations: %d\n", N_IMPORTANCE_SIMS))

# ============================================
# 13. DEFINE MODELS
# ============================================

cat("\n=== MODEL SPECIFICATIONS ===\n")

elastic_net_spec <- linear_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>%
  set_mode("regression")

pls_spec <- pls(num_comp = tune(), predictor_prop = tune()) %>%
  set_engine("mixOmics") %>%
  set_mode("regression")

rf_spec <- rand_forest(trees = tune(), mtry = tune(), min_n = tune()) %>%
  set_engine("ranger", importance = "permutation") %>%
  set_mode("regression")

xgb_spec <- boost_tree(
  trees = tune(), tree_depth = tune(), min_n = tune(),
  loss_reduction = tune(), sample_size = tune(), mtry = tune(), learn_rate = tune()
) %>%
  set_engine("xgboost") %>%
  set_mode("regression")

svm_spec <- svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
  set_engine("kernlab") %>%
  set_mode("regression")

knn_spec <- nearest_neighbor(neighbors = tune(), weight_func = tune(), dist_power = tune()) %>%
  set_engine("kknn") %>%
  set_mode("regression")

model_specs <- list(
  elastic_net = elastic_net_spec,
  pls = pls_spec,
  random_forest = rf_spec,
  xgboost = xgb_spec,
  svm = svm_spec,
  knn = knn_spec
)

cat("Models: Elastic Net, PLS, Random Forest, XGBoost, SVM-RBF, kNN\n")

# ============================================
# 14. MAIN MODELING FUNCTION
# ============================================

run_models_for_dataset <- function(df, outcome, dataset_name, 
                                   model_specs, cv_folds = 5, cv_repeats = 3,
                                   grid_size = 10, seed = 42,
                                   calculate_importance = FALSE, # switch to true for importance calculation of predictors
                                   n_importance_sims = 3,
                                   holdout_ids = NULL) {
  
  set.seed(seed)
  
  cvd_scores_local <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y", "mean_risk")
  
  if (!outcome %in% names(df)) {
    warning(sprintf("Outcome '%s' not found in dataset '%s'", outcome, dataset_name))
    return(NULL)
  }
  
  # Prepare data - df contains ALL participants (including holdout)
  other_scores <- setdiff(cvd_scores_local, outcome)
  
  df_full <- df %>% 
    select(-any_of(other_scores)) %>%
    filter(!is.na(.data[[outcome]]))
  
  # Split into train and holdout INTERNALLY
  if (!is.null(holdout_ids) && length(holdout_ids) > 0) {
    df_train <- df_full %>% filter(!Sample_ID %in% holdout_ids)
    holdout_clean <- df_full %>% filter(Sample_ID %in% holdout_ids)
    n_holdout <- nrow(holdout_clean)
    
    cat(sprintf("    Split: %d train, %d holdout (of %d IDs requested)\n", 
                nrow(df_train), n_holdout, length(holdout_ids)))
  } else {
    df_train <- df_full
    holdout_clean <- NULL
    n_holdout <- NA
  }
  
  # ============================================
  # PER-DATASET PARTICIPANT EXCLUSION (>80% missing)
  # ============================================
  predictor_cols_here <- setdiff(names(df_train), c("Sample_ID", outcome))
  
  if (length(predictor_cols_here) > 0) {
    # Check train set
    train_row_missing <- rowMeans(is.na(df_train[predictor_cols_here])) * 100
    high_missing_rows <- train_row_missing > 80
    
    if (any(high_missing_rows)) {
      removed_ids <- df_train$Sample_ID[high_missing_rows]
      n_removed <- length(removed_ids)
      
      cat(sprintf("    Removing %d participant(s) with >80%% missing in %s:\n", 
                  n_removed, dataset_name))
      cat(sprintf("      IDs: %s\n", paste(removed_ids, collapse = ", ")))
      
      # Show their missingness percentages
      for (id in removed_ids) {
        pct <- train_row_missing[df_train$Sample_ID == id]
        cat(sprintf("        %s: %.1f%% missing\n", id, pct))
      }
      
      df_train <- df_train[!high_missing_rows, ]
    }
    
    # Also check holdout set
    if (!is.null(holdout_clean) && nrow(holdout_clean) > 0) {
      holdout_row_missing <- rowMeans(is.na(holdout_clean[predictor_cols_here])) * 100
      holdout_high_missing <- holdout_row_missing > 80
      
      if (any(holdout_high_missing)) {
        removed_holdout_ids <- holdout_clean$Sample_ID[holdout_high_missing]
        cat(sprintf("    Also removing %d from holdout: %s\n", 
                    length(removed_holdout_ids), 
                    paste(removed_holdout_ids, collapse = ", ")))
        holdout_clean <- holdout_clean[!holdout_high_missing, ]
      }
    }
  }
  # ============================================
  
  n_train <- nrow(df_train)
  
  if (n_train < 30) {
    warning(sprintf("Too few training observations (%d) for %s ~ %s", n_train, outcome, dataset_name))
    return(NULL)
  }
  
  n_predictors <- ncol(df_train) - 2  # Minus Sample_ID and outcome
  
  cat(sprintf("  %s ~ %s: n_train=%d, n_holdout=%d, p=%d\n", 
              outcome, dataset_name, n_train, n_holdout, n_predictors))
  
  # CV splits with strata
  cv_splits <- vfold_cv(df_train, v = cv_folds, repeats = cv_repeats, 
                        strata = !!sym(outcome))
  
  # Recipe with imputation
  rec <- recipe(as.formula(paste(outcome, "~ .")), data = df_train) %>%
    update_role(Sample_ID, new_role = "ID") %>%
    step_zv(all_predictors()) %>%
    step_impute_mode(all_nominal_predictors()) %>%
    step_impute_knn(all_numeric_predictors(), neighbors = 5) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_zv(all_predictors()) %>%
    step_YeoJohnson(all_numeric_predictors()) %>%
    step_normalize(all_numeric_predictors()) %>%
    step_dummy(all_nominal_predictors())
  
  wf_set <- workflow_set(
    preproc = list(recipe = rec),
    models = model_specs,
    cross = TRUE
  )
  
  wf_results <- wf_set %>%
    workflow_map(
      fn = "tune_grid",
      resamples = cv_splits,
      grid = grid_size,
      metrics = metric_set(mae, rmse, rsq),
      control = control_grid(
        save_pred = TRUE,
        save_workflow = FALSE,
        verbose = FALSE,
        parallel_over = "everything"
      ),
      seed = seed,
      verbose = TRUE
    )
  
  # ---- EXTRACT COMPREHENSIVE METRICS ----
  all_model_metrics <- tibble()
  
  for (model_name in wf_results$wflow_id) {
    
    best_result <- wf_results %>%
      extract_workflow_set_result(model_name)
    
    best_params <- best_result %>%
      select_best(metric = "rmse")
    
    cv_per_fold <- best_result %>%
      collect_metrics(summarize = FALSE) %>%
      filter(.config == best_params$.config)
    
    # CV summary
    cv_summary <- cv_per_fold %>%
      group_by(.metric) %>%
      summarise(
        cv_mean = mean(.estimate, na.rm = TRUE),
        cv_sd = sd(.estimate, na.rm = TRUE),
        cv_se = sd(.estimate, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )
    
    # Fit final model on training data
    final_wf <- wf_results %>%
      extract_workflow(model_name) %>%
      finalize_workflow(best_params)
    
    final_fit <- tryCatch({
      fit(final_wf, data = df_train)
    }, error = function(e) NULL)
    
    # ---- TRAINING METRICS ----
    train_mae <- NA; train_rmse <- NA; train_rsq <- NA
    train_cal_slope <- NA; train_cal_intercept <- NA
    
    if (!is.null(final_fit)) {
      train_preds <- predict(final_fit, new_data = df_train) %>%
        bind_cols(df_train %>% select(all_of(outcome)))
      
      train_mae <- mean(abs(train_preds$.pred - train_preds[[outcome]]), na.rm = TRUE)
      train_rmse <- sqrt(mean((train_preds$.pred - train_preds[[outcome]])^2, na.rm = TRUE))
      train_rsq <- cor(train_preds$.pred, train_preds[[outcome]], use = "complete.obs")^2
      
      cal_model <- lm(as.formula(paste(outcome, "~ .pred")), data = train_preds)
      train_cal_intercept <- coef(cal_model)[1]
      train_cal_slope <- coef(cal_model)[2]
    }
    
    # ---- HOLDOUT METRICS ----
    holdout_mae <- NA; holdout_rmse <- NA; holdout_rsq <- NA
    holdout_cal_slope <- NA; holdout_cal_intercept <- NA
    
    if (!is.null(holdout_clean) && nrow(holdout_clean) >= 5 && !is.null(final_fit)) {
      holdout_preds <- tryCatch({
        predict(final_fit, new_data = holdout_clean) %>%
          bind_cols(holdout_clean %>% select(all_of(outcome)))
      }, error = function(e) NULL)
      
      if (!is.null(holdout_preds) && nrow(holdout_preds) > 0) {
        holdout_mae <- mean(abs(holdout_preds$.pred - holdout_preds[[outcome]]), na.rm = TRUE)
        holdout_rmse <- sqrt(mean((holdout_preds$.pred - holdout_preds[[outcome]])^2, na.rm = TRUE))
        
        if (sd(holdout_preds$.pred, na.rm = TRUE) > 0.001) {
          holdout_rsq <- cor(holdout_preds$.pred, holdout_preds[[outcome]], use = "complete.obs")^2
          
          cal_model <- lm(as.formula(paste(outcome, "~ .pred")), data = holdout_preds)
          holdout_cal_intercept <- coef(cal_model)[1]
          holdout_cal_slope <- coef(cal_model)[2]
        }
      }
    }
    
    # Extract CV metrics
    cv_mae <- cv_summary %>% filter(.metric == "mae") %>% pull(cv_mean)
    cv_mae_sd <- cv_summary %>% filter(.metric == "mae") %>% pull(cv_sd)
    cv_rmse <- cv_summary %>% filter(.metric == "rmse") %>% pull(cv_mean)
    cv_rmse_sd <- cv_summary %>% filter(.metric == "rmse") %>% pull(cv_sd)
    cv_rsq <- cv_summary %>% filter(.metric == "rsq") %>% pull(cv_mean)
    cv_rsq_sd <- cv_summary %>% filter(.metric == "rsq") %>% pull(cv_sd)
    
    # Calculate gaps and ratio
    # rsq_gap: positive = overfitting (train better than holdout)
    rsq_gap <- train_rsq - holdout_rsq
    mae_gap <- holdout_mae - train_mae
    rmse_gap <- holdout_rmse - train_rmse
    cv_rsq_gap <- train_rsq - cv_rsq
    
    # rsq_ratio: >1 means train better than holdout (overfitting)
    rsq_ratio <- ifelse(!is.na(holdout_rsq) & holdout_rsq > 0.01, 
                        train_rsq / holdout_rsq, NA)
    
    model_metrics <- tibble(
      dataset = dataset_name,
      outcome = outcome,
      model = str_remove(model_name, "recipe_"),
      n_train = n_train,
      n_holdout = n_holdout,
      n_predictors = n_predictors,
      
      # Training metrics
      train_mae = round(train_mae, 3),
      train_rmse = round(train_rmse, 3),
      train_rsq = round(train_rsq, 3),
      train_cal_slope = round(train_cal_slope, 3),
      train_cal_intercept = round(train_cal_intercept, 3),
      
      # CV metrics
      cv_mae = round(cv_mae, 3),
      cv_mae_sd = round(cv_mae_sd, 3),
      cv_mae_display = sprintf("%.2f ± %.2f", cv_mae, cv_mae_sd),
      cv_rmse = round(cv_rmse, 3),
      cv_rmse_sd = round(cv_rmse_sd, 3),
      cv_rmse_display = sprintf("%.2f ± %.2f", cv_rmse, cv_rmse_sd),
      cv_rsq = round(cv_rsq, 3),
      cv_rsq_sd = round(cv_rsq_sd, 3),
      cv_rsq_display = sprintf("%.2f ± %.2f", cv_rsq, cv_rsq_sd),
      
      # Holdout metrics
      holdout_mae = round(holdout_mae, 3),
      holdout_rmse = round(holdout_rmse, 3),
      holdout_rsq = round(holdout_rsq, 3),
      holdout_cal_slope = round(holdout_cal_slope, 3),
      holdout_cal_intercept = round(holdout_cal_intercept, 3),
      
      # Overfitting indicators
      rsq_gap = round(rsq_gap, 3),
      rsq_ratio = round(rsq_ratio, 3),
      mae_gap = round(mae_gap, 3),
      rmse_gap = round(rmse_gap, 3),
      cv_rsq_gap = round(cv_rsq_gap, 3),
      
      # Overfit flag: catches memorization (train_rsq ≈ 1) OR large gaps
      overfit_flag = case_when(
        is.na(holdout_rsq) ~ NA,
        train_rsq > 0.95 ~ TRUE,       # Near-perfect training fit = memorization
        rsq_gap > 0.3 ~ TRUE,          # Large gap = poor generalization
        !is.na(rsq_ratio) & rsq_ratio > 1.5 ~ TRUE,  # Train >> holdout
        TRUE ~ FALSE
      )
    )
    
    all_model_metrics <- bind_rows(all_model_metrics, model_metrics)
  }
  
  # ---- SELECT BEST MODEL (OBJECTIVE) ----
  best_model_name <- select_best_model(all_model_metrics)
  
  cat(sprintf("    Best model: %s\n", best_model_name))
  
  # ---- VARIABLE IMPORTANCE FOR BEST MODEL ----
  importance_results <- NULL
  
  if (calculate_importance) {
    cat(sprintf("    Calculating importance for: %s\n", best_model_name))
    
    best_params <- wf_results %>%
      extract_workflow_set_result(paste0("recipe_", best_model_name)) %>%
      select_best(metric = "rmse")
    
    final_wf <- wf_results %>%
      extract_workflow(paste0("recipe_", best_model_name)) %>%
      finalize_workflow(best_params)
    
    final_fit <- fit(final_wf, data = df_train)
    
    importance_results <- tryCatch({
      vi_permute(
        object = final_fit,
        train = df_train,
        target = outcome,
        metric = "rmse",
        nsim = n_importance_sims,
        pred_wrapper = function(object, newdata) {
          predict(object, new_data = newdata)$.pred
        }
      ) %>%
        mutate(
          dataset = dataset_name,
          outcome = outcome,
          best_model = best_model_name
        ) %>%
        filter(Importance > 0)
    }, error = function(e) {
      cat(sprintf("    Warning: Importance failed: %s\n", e$message))
      NULL
    })
  }
  
  # ---- CV PREDICTIONS FOR BEST MODEL ----
  best_result <- wf_results %>%
    extract_workflow_set_result(paste0("recipe_", best_model_name))
  
  best_config <- best_result %>%
    select_best(metric = "rmse")
  
  cv_predictions <- tryCatch({
    best_result %>%
      collect_predictions() %>%
      filter(.config == best_config$.config) %>%
      mutate(
        dataset = dataset_name,
        outcome_name = outcome,
        best_model = best_model_name
      )
  }, error = function(e) NULL)
  
  # ---- HOLDOUT PREDICTIONS FOR BEST MODEL ----
  holdout_predictions <- NULL
  if (!is.null(holdout_clean) && nrow(holdout_clean) >= 5) {
    final_wf <- wf_results %>%
      extract_workflow(paste0("recipe_", best_model_name)) %>%
      finalize_workflow(best_config)
    
    final_fit <- fit(final_wf, data = df_train)
    
    holdout_predictions <- tryCatch({
      predict(final_fit, new_data = holdout_clean) %>%
        bind_cols(holdout_clean %>% select(Sample_ID, all_of(outcome))) %>%
        mutate(
          dataset = dataset_name,
          outcome_name = outcome,
          best_model = best_model_name
        )
    }, error = function(e) NULL)
  }
  
  return(list(
    metrics = all_model_metrics,
    importance = importance_results,
    cv_predictions = cv_predictions,
    holdout_predictions = holdout_predictions,
    wf_results = wf_results,
    best_model = best_model_name
  ))
}


# ============================================
# BEST MODEL SELECTION - OBJECTIVE VERSION
# ============================================

select_best_model <- function(metrics_df) {
  # Objective selection: best holdout performance with acceptable overfitting
  # NO preference for interpretable models - justify choice in discussion
  
  candidates <- metrics_df %>%
    filter(!is.na(holdout_mae), !is.na(rsq_gap))
  
  # If no holdout metrics, fall back to CV
  if (nrow(candidates) == 0) {
    return(metrics_df %>% 
             slice_min(cv_rmse, n = 1) %>% 
             pull(model) %>% 
             `[`(1))
  }
  
  # Exclude severely overfit models:
  # - train_rsq > 0.95 (memorization)
  # - rsq_gap > 0.4 (large generalization gap)
  reasonable_models <- candidates %>%
    filter(train_rsq < 0.95, rsq_gap < 0.4)
  
  # If all models are overfit, pick the least overfit
  if (nrow(reasonable_models) == 0) {
    # Create composite overfit score: penalize high train_rsq AND high gap
    return(candidates %>%
             mutate(overfit_score = (train_rsq * 0.5) + (rsq_gap * 0.5)) %>%
             slice_min(overfit_score, n = 1) %>%
             pull(model) %>%
             `[`(1))
  }
  
  # Among reasonable models, pick best holdout MAE
  reasonable_models %>%
    slice_min(holdout_mae, n = 1) %>%
    pull(model) %>%
    `[`(1)
}

# ============================================
# 15. BUILD TASK LIST
# ============================================

cat("\n=== BUILDING TASK LIST ===\n")

task_list <- list()
task_counter <- 1

# A) Predictor blocks × CVD scores
for (dataset_name in names(datasets_list)) {
  for (outcome in CVD_scores) {
    task_list[[task_counter]] <- list(
      df = datasets_list[[dataset_name]],
      holdout_ids = HOLDOUT_IDS,
      outcome = outcome,
      dataset_name = dataset_name,
      category = "predictor_block"
    )
    task_counter <- task_counter + 1
  }
}

# B) Score-specific
for (score_name in names(score_specific_datasets)) {
  task_list[[task_counter]] <- list(
    df = score_specific_datasets[[score_name]],
    holdout_ids = HOLDOUT_IDS,
    outcome = score_name,
    dataset_name = paste0(score_name, "_specific"),
    category = "score_specific"
  )
  task_counter <- task_counter + 1
}

# C) Score + block combos
for (combo_name in names(score_plus_block_datasets)) {
  outcome <- sub("_plus_.*", "", combo_name)
  task_list[[task_counter]] <- list(
    df = score_plus_block_datasets[[combo_name]],
    holdout_ids = HOLDOUT_IDS,
    outcome = outcome,
    dataset_name = combo_name,
    category = "score_plus_block"
  )
  task_counter <- task_counter + 1
}

cat(sprintf("Total tasks: %d\n", length(task_list)))
cat(sprintf("  Predictor blocks: %d\n", length(datasets_list) * length(CVD_scores)))
cat(sprintf("  Score-specific: %d\n", length(score_specific_datasets)))
cat(sprintf("  Score + block: %d\n", length(score_plus_block_datasets)))


# # ============================================
# # TEST RUN: LIPIDS ~ QRISK3 ONLY
# # ============================================
# 
# cat("\n=== TEST RUN: Lipids ~ QRISK3_risk ===\n")
# 
# # Check the dataset first
# cat("\nDataset diagnostics:\n")
# cat(sprintf("  Rows: %d\n", nrow(datasets_list$lipids)))
# cat(sprintf("  Columns: %d\n", ncol(datasets_list$lipids)))
# cat(sprintf("  Predictors: %d\n", ncol(datasets_list$lipids) - length(CVD_scores) - 1))
# 
# # Check missingness in lipids dataset
# lipids_predictors <- setdiff(names(datasets_list$lipids), c("Sample_ID", CVD_scores))
# lipids_missing <- sapply(datasets_list$lipids[lipids_predictors], function(x) mean(is.na(x)) * 100)
# 
# cat("\nMissingness in lipids predictors:\n")
# cat(sprintf("  Mean: %.1f%%\n", mean(lipids_missing)))
# cat(sprintf("  Median: %.1f%%\n", median(lipids_missing)))
# cat(sprintf("  Max: %.1f%%\n", max(lipids_missing)))
# cat(sprintf("  Variables with >30%% missing: %d\n", sum(lipids_missing > 30)))
# 
# if (any(lipids_missing > 30)) {
#   cat("\n  Variables with >30% missing:\n")
#   high_missing_vars <- names(lipids_missing)[lipids_missing > 30]
#   for (v in high_missing_vars) {
#     cat(sprintf("    %.1f%% - %s\n", lipids_missing[v], v))
#   }
# }
# 
# # Row-level missingness
# lipids_row_missing <- rowMeans(is.na(datasets_list$lipids[lipids_predictors])) * 100
# cat(sprintf("\nParticipant missingness in lipids:\n"))
# cat(sprintf("  Mean: %.1f%%\n", mean(lipids_row_missing)))
# cat(sprintf("  Max: %.1f%%\n", max(lipids_row_missing)))
# 
# # Start parallel
# plan(multisession, workers = CPUS)
# cat(sprintf("\nParallel workers: %d\n", future::nbrOfWorkers()))
# 
# # Run the model
# cat("\n--- Starting model run ---\n")
# start_time <- Sys.time()
# 
# test_result <- tryCatch({
#   run_models_for_dataset(
#     df = datasets_list$lipids,
#     outcome = "QRISK3_risk",
#     dataset_name = "lipids",
#     model_specs = model_specs,
#     cv_folds = CV_FOLDS,
#     cv_repeats = CV_REPEATS,
#     grid_size = GRID_SIZE,
#     seed = 42,
#     calculate_importance = FALSE,  # Skip importance for speed
#     n_importance_sims = N_IMPORTANCE_SIMS,
#     holdout_ids = HOLDOUT_IDS
#   )
# }, error = function(e) {
#   cat(sprintf("\nERROR: %s\n", e$message))
#   cat("\nFull error:\n")
#   print(e)
#   return(NULL)
# })
# 
# plan(sequential)
# 
# elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
# cat(sprintf("\n--- Test completed in %.2f minutes ---\n", elapsed))
# 
# if (!is.null(test_result)) {
#   cat("\n✓ SUCCESS!\n")
#   cat(sprintf("  Best model: %s\n", test_result$best_model))
#   
#   cat("\n--- Metrics ---\n")
#   print(test_result$metrics %>% 
#           select(model, train_rsq, cv_rsq, holdout_rsq, 
#                  holdout_mae, rsq_gap, overfit_flag))
# } else {
#   cat("\n✗ FAILED - see error above\n")
# }
# 
# # Save result if successful
# if (!is.null(test_result)) {
#   qs2::qs_save(test_result, "test_lipids_qrisk3.qs2")
#   cat("\n✓ Saved to test_lipids_qrisk3.qs2\n")
# }


# # ============================================
# # 16. TEST RUN — SCORE-SPECIFIC ONLY (5 tasks)
# # ============================================
# 
# cat("\n=== TEST RUN: SCORE-SPECIFIC DATASETS ONLY ===\n")
# 
# # Filter to score-specific tasks only
# score_specific_tasks <- task_list[sapply(task_list, function(x) x$category == "score_specific")]
# 
# cat(sprintf("Tasks to run: %d\n", length(score_specific_tasks)))
# for (t in score_specific_tasks) {
#   cat(sprintf("  - %s ~ %s\n", t$outcome, t$dataset_name))
# }
# 
# plan(multisession, workers = CPUS)
# cat(sprintf("\nParallel workers: %d\n", future::nbrOfWorkers()))
# cat(sprintf("Started at: %s\n", Sys.time()))
# 
# test_results <- list()
# failed_tasks <- c()
# start_time <- Sys.time()
# 
# for (i in seq_along(score_specific_tasks)) {
#   task <- score_specific_tasks[[i]]
#   task_name <- paste0(task$outcome, " ~ ", task$dataset_name)
# 
#   cat(sprintf("\n[%d/%d] %s\n", i, length(score_specific_tasks), task_name))
# 
#   result <- tryCatch({
#     run_models_for_dataset(
#       df = task$df,
#       outcome = task$outcome,
#       dataset_name = task$dataset_name,
#       model_specs = model_specs,
#       cv_folds = CV_FOLDS,
#       cv_repeats = CV_REPEATS,
#       grid_size = GRID_SIZE,
#       seed = 42,
#       calculate_importance = FALSE,
#       n_importance_sims = N_IMPORTANCE_SIMS,
#       holdout_ids = task$holdout_ids
#     )
#   }, error = function(e) {
#     cat(sprintf("  ERROR: %s\n", e$message))
#     return(NULL)
#   })
# 
#   if (!is.null(result)) {
#     test_results[[task_name]] <- result
#     cat(sprintf("  ✓ Completed | Best: %s\n", result$best_model))
#   } else {
#     failed_tasks <- c(failed_tasks, task_name)
#     cat("  ✗ Failed\n")
#   }
# }
# 
# plan(sequential)
# 
# total_time <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
# cat(sprintf("\n=== TEST RUN COMPLETED ===\n"))
# cat(sprintf("Total time: %.1f minutes\n", total_time))
# cat(sprintf("Successful: %d | Failed: %d\n", length(test_results), length(failed_tasks)))
# 
# # ---- REVIEW RESULTS ----
# 
# cat("\n=== RESULTS REVIEW ===\n")
# 
# # Compile metrics
# test_metrics <- bind_rows(lapply(test_results, function(x) x$metrics))
# 
# # Verify all models present
# cat("\n--- Models per task ---\n")
# print(table(test_metrics$model))
# 
# # Best model selected for each
# cat("\n--- Best Model Selected Per Task ---\n")
# best_selected_test <- tibble(
#   task = names(test_results),
#   best_model = sapply(test_results, function(x) x$best_model)
# )
# print(best_selected_test)
# 
# # Selection counts
# cat("\n--- Selection Counts ---\n")
# print(table(best_selected_test$best_model))
# 
# # Full metrics table with new columns
# cat("\n--- All Models Comparison (with rsq_ratio) ---\n")
# print(test_metrics %>%
#         select(dataset, outcome, model,
#                train_rsq, holdout_rsq, rsq_gap, rsq_ratio,
#                holdout_mae, holdout_rmse, holdout_cal_slope, overfit_flag) %>%
#         arrange(dataset, desc(overfit_flag), rsq_gap))
# 
# # Summary by model
# cat("\n--- Model Summary Across All Tasks ---\n")
# model_summary <- test_metrics %>%
#   group_by(model) %>%
#   summarise(
#     n_selected = sum(model == best_selected_test$best_model[match(paste(outcome, "~", dataset),
#                                                                   best_selected_test$task)], na.rm = TRUE),
#     mean_train_rsq = round(mean(train_rsq, na.rm = TRUE), 3),
#     mean_holdout_rsq = round(mean(holdout_rsq, na.rm = TRUE), 3),
#     mean_rsq_gap = round(mean(rsq_gap, na.rm = TRUE), 3),
#     mean_rsq_ratio = round(mean(rsq_ratio, na.rm = TRUE), 3),
#     mean_holdout_mae = round(mean(holdout_mae, na.rm = TRUE), 3),
#     mean_cal_slope = round(mean(holdout_cal_slope, na.rm = TRUE), 3),
#     pct_overfit = round(mean(overfit_flag, na.rm = TRUE) * 100, 1),
#     .groups = "drop"
#   ) %>%
#   arrange(desc(n_selected), mean_rsq_gap)
# 
# print(model_summary)
# 
# # Flag summary
# cat("\n--- Overfitting Summary ---\n")
# cat("Overfit criteria: train_rsq > 0.95 OR rsq_gap > 0.3 OR rsq_ratio > 1.5\n")
# overfit_summary <- test_metrics %>%
#   group_by(model) %>%
#   summarise(
#     n_overfit = sum(overfit_flag, na.rm = TRUE),
#     n_total = n(),
#     pct_overfit = round(mean(overfit_flag, na.rm = TRUE) * 100, 1),
#     .groups = "drop"
#   )
# print(overfit_summary)
# 
# # Save test results
# qs2::qs_save(test_results, "ml_score_specific_test_results.qs2")
# qs2::qs_save(test_metrics, "ml_score_specific_test_metrics.qs2")
# 
# cat("\n✓ Test results saved to ml_score_specific_test_metrics.qs2\n")
# cat("\n=== REVIEW COMPLETE — Ready for full run if satisfied ===\n")

# ============================================
# 17. FULL RUN
# ============================================

plan(multisession, workers = CPUS)

all_results <- list()
failed_tasks <- c()
start_time <- Sys.time()

for (i in seq_along(task_list)) {
  task <- task_list[[i]]
  task_name <- paste0(task$outcome, " ~ ", task$dataset_name)

  cat(sprintf("\n[%d/%d] %s\n", i, length(task_list), task_name))

  result <- tryCatch({
    run_models_for_dataset(
      df = task$df,
      outcome = task$outcome,
      dataset_name = task$dataset_name,
      model_specs = model_specs,
      cv_folds = CV_FOLDS,
      cv_repeats = CV_REPEATS,
      grid_size = GRID_SIZE,
      seed = 42,
      calculate_importance = FALSE,
      n_importance_sims = N_IMPORTANCE_SIMS,
      holdout_ids = task$holdout_ids
    )
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    return(NULL)
  })

  if (!is.null(result)) {
    all_results[[task_name]] <- result
    cat("  ✓ Completed\n")
  } else {
    failed_tasks <- c(failed_tasks, task_name)
    cat("  ✗ Failed\n")
  }

  # Intermediate save every 10 tasks
  if (i %% 10 == 0) {
    intermediate_metrics <- bind_rows(lapply(all_results, function(x) x$metrics))
    qs2::qs_save(intermediate_metrics, "ml_comparison_intermediate.qs2")
    qs2::qs_save(all_results, "ml_comparison_all_results_partial.qs2")
    cat(sprintf("  [Intermediate save: %d tasks]\n", i))
  }
}

plan(sequential)

total_time <- round(difftime(Sys.time(), start_time, units = "hours"), 2)
cat(sprintf("\n✓ Completed in %.2f hours\n", total_time))
cat(sprintf("  Successful: %d | Failed: %d\n", length(all_results), length(failed_tasks)))

qs2::qs_save(all_results, "ml_comparison_all_results.qs2")
if (length(failed_tasks) > 0) {
  writeLines(failed_tasks, "ml_comparison_failed_tasks.txt")
}

# ============================================
# 18. COMPILE RESULTS
# ============================================

all_metrics <- bind_rows(lapply(all_results, function(x) x$metrics))
all_importance <- bind_rows(lapply(all_results, function(x) x$importance))
all_cv_preds <- bind_rows(lapply(all_results, function(x) x$cv_predictions))
all_holdout_preds <- bind_rows(lapply(all_results, function(x) x$holdout_predictions))

# Best models summary - get the selected best model for each task
best_models_summary <- tibble()
for (task_name in names(all_results)) {
  result <- all_results[[task_name]]
  best_row <- result$metrics %>%
    filter(model == result$best_model)
  best_models_summary <- bind_rows(best_models_summary, best_row)
}

# Model rankings based on holdout MAE
model_rankings <- all_metrics %>%
  filter(!is.na(holdout_mae)) %>%
  group_by(dataset, outcome) %>%
  mutate(rank = rank(holdout_mae)) %>%
  ungroup() %>%
  group_by(model) %>%
  summarise(
    mean_rank = round(mean(rank), 2),
    times_best = sum(rank == 1),
    times_top3 = sum(rank <= 3),
    mean_train_rsq = round(mean(train_rsq, na.rm = TRUE), 3),
    mean_cv_rsq = round(mean(cv_rsq, na.rm = TRUE), 3),
    mean_holdout_rsq = round(mean(holdout_rsq, na.rm = TRUE), 3),
    mean_rsq_gap = round(mean(rsq_gap, na.rm = TRUE), 3),
    mean_rsq_ratio = round(mean(rsq_ratio, na.rm = TRUE), 3),
    mean_holdout_mae = round(mean(holdout_mae, na.rm = TRUE), 3),
    mean_holdout_rmse = round(mean(holdout_rmse, na.rm = TRUE), 3),
    mean_cal_slope = round(mean(holdout_cal_slope, na.rm = TRUE), 3),
    pct_overfit = round(mean(overfit_flag, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(mean_rank)

cat("\n--- Model Rankings ---\n")
print(model_rankings)

# Times each model was selected as best
cat("\n--- Best Model Selection Counts ---\n")
selection_counts <- table(sapply(all_results, function(x) x$best_model))
print(selection_counts)

# Overfitting summary by model
cat("\n--- Overfitting Summary by Model ---\n")
overfit_by_model <- all_metrics %>%
  group_by(model) %>%
  summarise(
    n_tasks = n(),
    n_overfit = sum(overfit_flag, na.rm = TRUE),
    pct_overfit = round(mean(overfit_flag, na.rm = TRUE) * 100, 1),
    mean_train_rsq = round(mean(train_rsq, na.rm = TRUE), 3),
    n_train_rsq_gt_95 = sum(train_rsq > 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_overfit))
print(overfit_by_model)

# Performance by dataset category
cat("\n--- Performance by Dataset Category ---\n")
performance_by_category <- all_metrics %>%
  left_join(
    tibble(
      dataset = sapply(task_list, function(x) x$dataset_name),
      category = sapply(task_list, function(x) x$category)
    ) %>% distinct(),
    by = "dataset"
  ) %>%
  group_by(category) %>%
  summarise(
    n_tasks = n_distinct(paste(dataset, outcome)),
    mean_holdout_rsq = round(mean(holdout_rsq, na.rm = TRUE), 3),
    mean_holdout_mae = round(mean(holdout_mae, na.rm = TRUE), 3),
    mean_rsq_gap = round(mean(rsq_gap, na.rm = TRUE), 3),
    pct_overfit = round(mean(overfit_flag, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )
print(performance_by_category)

# Performance by outcome
cat("\n--- Performance by Outcome ---\n")
performance_by_outcome <- all_metrics %>%
  group_by(outcome) %>%
  summarise(
    n_tasks = n_distinct(dataset),
    mean_holdout_rsq = round(mean(holdout_rsq, na.rm = TRUE), 3),
    mean_holdout_mae = round(mean(holdout_mae, na.rm = TRUE), 3),
    best_holdout_rsq = round(max(holdout_rsq, na.rm = TRUE), 3),
    best_model_for_outcome = model[which.max(holdout_rsq)],
    .groups = "drop"
  )
print(performance_by_outcome)

# Top variable importance across all tasks
cat("\n--- Top 20 Most Important Variables (across all tasks) ---\n")
if (nrow(all_importance) > 0) {
  top_importance <- all_importance %>%
    group_by(Variable) %>%
    summarise(
      mean_importance = round(mean(Importance, na.rm = TRUE), 4),
      times_appeared = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_importance)) %>%
    head(20)
  print(top_importance)
}

# Excel output
excel_output <- list(
  "Best_Models" = best_models_summary,
  "All_Metrics" = all_metrics,
  "Model_Rankings" = model_rankings,
  "Variable_Importance" = all_importance,
  "Overfitting_Check" = all_metrics %>%
    select(dataset, outcome, model, train_rsq, cv_rsq, holdout_rsq,
           rsq_gap, rsq_ratio, overfit_flag) %>%
    arrange(desc(overfit_flag), desc(rsq_gap)),
  "Performance_by_Category" = performance_by_category,
  "Performance_by_Outcome" = performance_by_outcome
)

write_xlsx(excel_output, "ml_model_comparison_results.xlsx")

# Save compiled R objects
qs2::qs_save(list(
  metrics = all_metrics,
  importance = all_importance,
  cv_preds = all_cv_preds,
  holdout_preds = all_holdout_preds,
  best_models = best_models_summary,
  model_rankings = model_rankings
), "ml_comparison_compiled_objects.qs2")

cat("\n✓ Results saved to ml_model_comparison_results.xlsx\n")
cat("✓ R objects saved to ml_comparison_compiled_objects.qs2\n")

# Print failed tasks if any
if (length(failed_tasks) > 0) {
  cat("\n--- Failed Tasks ---\n")
  for (ft in failed_tasks) {
    cat(sprintf("  - %s\n", ft))
  }
}

cat("\n=== ANALYSIS COMPLETE ===\n")


################################################################################
# Creating results table
################################################################################

# Load results
all_metrics <- qs2::qs_read("ml_comparison_compiled_objects.qs2")$metrics

# Dataset and model definitions (same as before)
dataset_order <- c(
  "QRISK3_risk_specific", "SCORE2_score_specific", "frs_10y_specific", 
  "ascvd_10y_specific", "mean_risk_specific",
  "all_data",
  "fatty_acids",
  "lipids",
  "urine_NMR",
  "body_composition",
  "sociodemographics_lifestyle",
  "clinical_risk_factors"
)

dataset_labels <- c(
  "QRISK3_risk_specific" = "Score-specific inputs",
  "SCORE2_score_specific" = "Score-specific inputs",
  "frs_10y_specific" = "Score-specific inputs",
  "ascvd_10y_specific" = "Score-specific inputs",
  "mean_risk_specific" = "Score-specific inputs",
  "all_data" = "All data",
  "fatty_acids" = "Fatty acids",
  "lipids" = "Lipidomics",
  "urine_NMR" = "Urinary NMR",
  "body_composition" = "Body composition",
  "sociodemographics_lifestyle" = "Sociodemographics/Lifestyle",
  "clinical_risk_factors" = "Clinical risk factors"
)

model_order <- c("elastic_net", "pls", "random_forest", "xgboost", "svm", "knn")
model_labels <- c(
  "elastic_net" = "Elastic Net",
  "pls" = "PLS",
  "random_forest" = "Random Forest",
  "xgboost" = "XGBoost",
  "svm" = "SVM-RBF",
  "knn" = "k-NN"
)

outcomes <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y", "mean_risk")
outcome_labels <- c(
  "QRISK3_risk" = "QRISK3 Risk",
  "SCORE2_score" = "SCORE2 Risk",
  "frs_10y" = "Framingham 10-year Risk",
  "ascvd_10y" = "ASCVD 10-year Risk",
  "mean_risk" = "Composite Score Risk"
)

# ============================================
# CREATE FORMATTED TABLE FUNCTION - COMPACT
# ============================================

create_publication_table <- function(outcome_name, metrics_df) {
  
  # Prepare data
  table_data <- metrics_df %>%
    filter(
      outcome == outcome_name,
      dataset %in% dataset_order
    ) %>%
    mutate(
      model = factor(model, levels = model_order),
      dataset = factor(dataset, levels = dataset_order)
    ) %>%
    arrange(dataset, model)
  
  # Format metrics
  formatted_table <- table_data %>%
    mutate(
      Dataset = dataset_labels[as.character(dataset)],
      Model = model_labels[as.character(model)],
      
      `Train R²` = sprintf("%.3f", train_rsq),
      
      `CV R²` = ifelse(!is.na(cv_rsq_sd) & cv_rsq_sd > 0,
                       sprintf("%.3f ± %.3f", cv_rsq, cv_rsq_sd),
                       sprintf("%.3f", cv_rsq)),
      
      `Holdout R²` = ifelse(!is.na(holdout_rsq), 
                            sprintf("%.3f", holdout_rsq), 
                            "—"),
      
      `R² Gap` = ifelse(!is.na(rsq_gap), 
                        sprintf("%.3f", rsq_gap), 
                        "—"),
      
      `R² Ratio` = ifelse(!is.na(rsq_ratio), 
                          sprintf("%.2f", rsq_ratio), 
                          "—"),
      
      `CV MAE` = ifelse(!is.na(cv_mae_sd) & cv_mae_sd > 0,
                        sprintf("%.3f ± %.3f", cv_mae, cv_mae_sd),
                        sprintf("%.3f", cv_mae)),
      
      `Holdout MAE` = ifelse(!is.na(holdout_mae), 
                             sprintf("%.3f", holdout_mae), 
                             "—"),
      
      `CV RMSE` = ifelse(!is.na(cv_rmse_sd) & cv_rmse_sd > 0,
                         sprintf("%.3f ± %.3f", cv_rmse, cv_rmse_sd),
                         sprintf("%.3f", cv_rmse)),
      
      `Holdout RMSE` = ifelse(!is.na(holdout_rmse), 
                              sprintf("%.3f", holdout_rmse), 
                              "—"),
      
      `Cal. Slope` = ifelse(!is.na(holdout_cal_slope), 
                            sprintf("%.3f", holdout_cal_slope), 
                            "—")
    ) %>%
    select(Dataset, Model, `Train R²`, `CV R²`, `Holdout R²`, 
           `R² Gap`, `R² Ratio`,
           `CV MAE`, `Holdout MAE`, 
           `CV RMSE`, `Holdout RMSE`, 
           `Cal. Slope`)
  
  # Create flextable with TIGHT widths
  ft <- formatted_table %>%
    flextable() %>%
    merge_v(j = "Dataset") %>%
    set_caption(caption = sprintf("Table: Model Performance for %s", 
                                  outcome_labels[outcome_name])) %>%
    align(j = 1:2, align = "left", part = "all") %>%
    align(j = 3:12, align = "center", part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    fontsize(size = 8, part = "body") %>%     # Size 8
    fontsize(size = 9, part = "header") %>%   # Size 9
    bold(part = "header") %>%
    bold(j = 1, part = "body") %>%
    # TIGHT column widths (in inches)
    width(j = 1, width = 1.5) %>%   # Dataset
    width(j = 2, width = 1.0) %>%   # Model
    width(j = 3, width = 0.6) %>%   # Train R²
    width(j = 4, width = 1.0) %>%   # CV R² (has ± SD, needs more space)
    width(j = 5, width = 0.65) %>%  # Holdout R²
    width(j = 6, width = 0.6) %>%   # R² Gap
    width(j = 7, width = 0.6) %>%   # R² Ratio
    width(j = 8, width = 1.0) %>%   # CV MAE (has ± SD)
    width(j = 9, width = 0.7) %>%   # Holdout MAE
    width(j = 10, width = 1.0) %>%  # CV RMSE (has ± SD)
    width(j = 11, width = 0.7) %>%  # Holdout RMSE
    width(j = 12, width = 0.65) %>% # Cal. Slope
    border_remove() %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1.5), part = "header") %>%
    hline(i = seq(6, nrow(formatted_table), 6), 
          border = fp_border(width = 1), part = "body") %>%
    vline(border = fp_border(width = 0.5, color = "#CCCCCC"), part = "all") %>%
    padding(padding = 2, part = "all") %>%   # Reduced padding
    hrule(rule = "exact")  # Force exact widths, no autofit
  
  return(ft)
}

# ============================================
# GENERATE ALL TABLES WITH LANDSCAPE ORIENTATION
# ============================================

cat("\n=== GENERATING PUBLICATION TABLES (LANDSCAPE, COMPACT) ===\n")

all_tables <- list()

for (outcome in outcomes) {
  cat(sprintf("Creating table for: %s\n", outcome_labels[outcome]))
  
  ft <- create_publication_table(outcome, all_metrics)
  all_tables[[outcome]] <- ft
  
  # Create Word document with landscape orientation
  doc <- read_docx() %>%
    body_add_par("", style = "Normal") %>%
    body_end_section_landscape() %>%
    body_add_flextable(value = ft)
  
  print(doc, target = sprintf("Table_%s.docx", outcome))
  
  # Also save as HTML
  save_as_html(ft, path = sprintf("Table_%s.html", outcome))
}

cat("\n✓ All tables saved as .docx (LANDSCAPE, COMPACT) and .html files\n")

# ============================================
# CREATE COMBINED DOCUMENT WITH ALL TABLES (LANDSCAPE)
# ============================================

doc <- read_docx()

for (i in seq_along(outcomes)) {
  outcome <- outcomes[i]
  
  # Add landscape section and table
  doc <- doc %>%
    body_add_par("", style = "Normal") %>%
    body_end_section_landscape() %>%
    body_add_flextable(value = all_tables[[outcome]])
  
  # Add page break if not last table
  if (i < length(outcomes)) {
    doc <- doc %>% body_add_break()
  }
}

print(doc, target = "All_Tables_Combined.docx")

cat("✓ Combined document saved as All_Tables_Combined.docx (LANDSCAPE, COMPACT)\n")


################################################################################
# Creating summary results table
################################################################################

# Load results
all_metrics <- qs2::qs_read("ml_comparison_compiled_objects.qs2")$metrics

model_order <- c("elastic_net", "pls", "random_forest", "xgboost", "svm", "knn")
model_labels <- c(
  "elastic_net" = "Elastic Net",
  "pls" = "PLS",
  "random_forest" = "Random Forest",
  "xgboost" = "XGBoost",
  "svm" = "SVM-RBF",
  "knn" = "k-NN"
)

# ============================================
# FIRST: Check what datasets we actually have
# ============================================

cat("\n=== Checking available datasets ===\n")
cat("Score-specific datasets:\n")
print(unique(all_metrics$dataset[str_detect(all_metrics$dataset, "_specific")]))

cat("\nGeneral datasets:\n")
general_check <- unique(all_metrics$dataset[!str_detect(all_metrics$dataset, "_specific") & 
                                              !str_detect(all_metrics$dataset, "_plus_")])
print(general_check)

cat("\nScore + block datasets:\n")
plus_check <- unique(all_metrics$dataset[str_detect(all_metrics$dataset, "_plus_")])
print(head(plus_check, 20))

# ============================================
# TABLE A: SCORE-SPECIFIC ONLY (MUST BE 6 ROWS)
# ============================================

cat("\n=== Preparing Table A: Score-Specific Inputs ===\n")

score_specific_datasets <- c("QRISK3_risk_specific", "SCORE2_score_specific", 
                             "frs_10y_specific", "ascvd_10y_specific", "mean_risk_specific")

# Debug: check how many rows before grouping
temp_a <- all_metrics %>%
  filter(dataset %in% score_specific_datasets)
cat(sprintf("Rows before grouping: %d\n", nrow(temp_a)))
cat(sprintf("Unique models: %d\n", n_distinct(temp_a$model)))
cat(sprintf("Unique datasets: %d\n", n_distinct(temp_a$dataset)))

# Average across ALL 5 score-specific datasets
table_a_data <- all_metrics %>%
  filter(dataset %in% score_specific_datasets) %>%
  group_by(model) %>%  # CRITICAL: Only group by model
  summarise(
    n_datasets = n(),  # Should be 5 (one per CVD score)
    train_rsq = mean(train_rsq, na.rm = TRUE),
    cv_rsq = mean(cv_rsq, na.rm = TRUE),
    cv_rsq_sd = mean(cv_rsq_sd, na.rm = TRUE),
    holdout_rsq = mean(holdout_rsq, na.rm = TRUE),
    rsq_gap = mean(rsq_gap, na.rm = TRUE),
    rsq_ratio = mean(rsq_ratio, na.rm = TRUE),
    cv_mae = mean(cv_mae, na.rm = TRUE),
    cv_mae_sd = mean(cv_mae_sd, na.rm = TRUE),
    holdout_mae = mean(holdout_mae, na.rm = TRUE),
    cv_rmse = mean(cv_rmse, na.rm = TRUE),
    cv_rmse_sd = mean(cv_rmse_sd, na.rm = TRUE),
    holdout_rmse = mean(holdout_rmse, na.rm = TRUE),
    holdout_cal_slope = mean(holdout_cal_slope, na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf("Table A rows after grouping: %d (should be 6)\n", nrow(table_a_data)))
print(table_a_data %>% select(model, n_datasets))

# Add dataset label AFTER grouping
table_a_data <- table_a_data %>%
  mutate(dataset_label = "Score-specific inputs") %>%
  select(-n_datasets)  # Remove debug column

# ============================================
# TABLE B: GENERAL DATASETS
# ============================================

cat("\n=== Preparing Table B: General Datasets ===\n")

general_datasets <- c("all_data", "fatty_acids", "lipids", "urine_NMR", 
                      "body_composition", "sociodemographics_lifestyle", 
                      "clinical_risk_factors")

dataset_labels_b <- c(
  "all_data" = "All data",
  "fatty_acids" = "Fatty acids",
  "lipids" = "Lipidomics",
  "urine_NMR" = "Urinary NMR",
  "body_composition" = "Body composition",
  "sociodemographics_lifestyle" = "Sociodemographics/Lifestyle",
  "clinical_risk_factors" = "Clinical risk factors"
)

table_b_data <- all_metrics %>%
  filter(dataset %in% general_datasets) %>%
  group_by(dataset, model) %>%
  summarise(
    train_rsq = mean(train_rsq, na.rm = TRUE),
    cv_rsq = mean(cv_rsq, na.rm = TRUE),
    cv_rsq_sd = mean(cv_rsq_sd, na.rm = TRUE),
    holdout_rsq = mean(holdout_rsq, na.rm = TRUE),
    rsq_gap = mean(rsq_gap, na.rm = TRUE),
    rsq_ratio = mean(rsq_ratio, na.rm = TRUE),
    cv_mae = mean(cv_mae, na.rm = TRUE),
    cv_mae_sd = mean(cv_mae_sd, na.rm = TRUE),
    holdout_mae = mean(holdout_mae, na.rm = TRUE),
    cv_rmse = mean(cv_rmse, na.rm = TRUE),
    cv_rmse_sd = mean(cv_rmse_sd, na.rm = TRUE),
    holdout_rmse = mean(holdout_rmse, na.rm = TRUE),
    holdout_cal_slope = mean(holdout_cal_slope, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(dataset_label = dataset_labels_b[dataset])

cat(sprintf("Table B rows: %d (should be ~42)\n", nrow(table_b_data)))

# ============================================
# TABLE C: SCORE-SPECIFIC + BLOCK
# ============================================

cat("\n=== Preparing Table C: Score-Specific + Block ===\n")

# Extract all blocks that exist
all_blocks <- all_metrics %>%
  filter(str_detect(dataset, "_plus_")) %>%
  mutate(block = str_extract(dataset, "(?<=_plus_).*")) %>%
  pull(block) %>%
  unique() %>%
  sort()

cat("Blocks found:\n")
print(all_blocks)

# Create comprehensive block labels
block_labels <- c(
  "all_data" = "All data",
  "fatty_acids" = "Fatty acids",
  "lipids" = "Lipidomics",
  "urine_NMR" = "Urinary NMR",
  "body_composition" = "Body composition",
  "sociodemographics_lifestyle" = "Sociodemographics/Lifestyle",
  "clinical_risk_factors" = "Clinical risk factors",
  "mtc_sign" = "MTC significant predictors",
  "sign_predictors" = "Significant predictors"
)

# Process Table C
score_plus_data <- all_metrics %>%
  filter(str_detect(dataset, "_plus_")) %>%
  mutate(
    score = str_extract(dataset, "^[^_]+_[^_]+(?=_plus)"),
    block = str_extract(dataset, "(?<=_plus_).*")
  ) %>%
  filter(!is.na(block), !is.na(score))

cat(sprintf("Score+block combinations before grouping: %d\n", nrow(score_plus_data)))

table_c_data <- score_plus_data %>%
  group_by(block, model) %>%
  summarise(
    n_scores = n(),
    train_rsq = mean(train_rsq, na.rm = TRUE),
    cv_rsq = mean(cv_rsq, na.rm = TRUE),
    cv_rsq_sd = mean(cv_rsq_sd, na.rm = TRUE),
    holdout_rsq = mean(holdout_rsq, na.rm = TRUE),
    rsq_gap = mean(rsq_gap, na.rm = TRUE),
    rsq_ratio = mean(rsq_ratio, na.rm = TRUE),
    cv_mae = mean(cv_mae, na.rm = TRUE),
    cv_mae_sd = mean(cv_mae_sd, na.rm = TRUE),
    holdout_mae = mean(holdout_mae, na.rm = TRUE),
    cv_rmse = mean(cv_rmse, na.rm = TRUE),
    cv_rmse_sd = mean(cv_rmse_sd, na.rm = TRUE),
    holdout_rmse = mean(holdout_rmse, na.rm = TRUE),
    holdout_cal_slope = mean(holdout_cal_slope, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    dataset_label = paste("Score-specific +", block_labels[block]),
    # Check for NA labels
    has_na = is.na(block_labels[block])
  )

cat(sprintf("Table C rows: %d\n", nrow(table_c_data)))

# Check for NA labels
if (any(table_c_data$has_na)) {
  cat("\nWARNING: Some blocks have no label:\n")
  print(table_c_data %>% filter(has_na) %>% select(block) %>% distinct())
}

table_c_data <- table_c_data %>%
  rename(dataset = block) %>%
  select(-n_scores, -has_na)

# ============================================
# FORMATTING FUNCTION
# ============================================

create_summary_table <- function(data, table_title) {
  
  # Get unique dataset labels in order they appear
  dataset_levels <- unique(data$dataset_label)
  
  # Format data
  formatted_table <- data %>%
    mutate(
      model = factor(model, levels = model_order),
      dataset_label = factor(dataset_label, levels = dataset_levels)
    ) %>%
    arrange(dataset_label, model) %>%
    mutate(
      Dataset = as.character(dataset_label),
      Model = model_labels[as.character(model)],
      
      `Train R²` = sprintf("%.3f", train_rsq),
      
      `CV R²` = ifelse(
        !is.na(cv_rsq_sd) & cv_rsq_sd > 0,
        sprintf("%.3f ± %.3f", cv_rsq, cv_rsq_sd),
        sprintf("%.3f", cv_rsq)
      ),
      
      `Holdout R²` = ifelse(!is.na(holdout_rsq), sprintf("%.3f", holdout_rsq), "—"),
      `R² Gap` = ifelse(!is.na(rsq_gap), sprintf("%.3f", rsq_gap), "—"),
      `R² Ratio` = ifelse(!is.na(rsq_ratio), sprintf("%.2f", rsq_ratio), "—"),
      
      `CV MAE` = ifelse(
        !is.na(cv_mae_sd) & cv_mae_sd > 0,
        sprintf("%.3f ± %.3f", cv_mae, cv_mae_sd),
        sprintf("%.3f", cv_mae)
      ),
      
      `Holdout MAE` = ifelse(!is.na(holdout_mae), sprintf("%.3f", holdout_mae), "—"),
      
      `CV RMSE` = ifelse(
        !is.na(cv_rmse_sd) & cv_rmse_sd > 0,
        sprintf("%.3f ± %.3f", cv_rmse, cv_rmse_sd),
        sprintf("%.3f", cv_rmse)
      ),
      
      `Holdout RMSE` = ifelse(!is.na(holdout_rmse), sprintf("%.3f", holdout_rmse), "—"),
      `Cal. Slope` = ifelse(!is.na(holdout_cal_slope), sprintf("%.3f", holdout_cal_slope), "—")
    ) %>%
    select(Dataset, Model, `Train R²`, `CV R²`, `Holdout R²`, 
           `R² Gap`, `R² Ratio`, `CV MAE`, `Holdout MAE`, 
           `CV RMSE`, `Holdout RMSE`, `Cal. Slope`)
  
  # Create flextable
  ft <- formatted_table %>%
    flextable() %>%
    merge_v(j = "Dataset") %>%
    set_caption(caption = table_title) %>%
    align(j = 1:2, align = "left", part = "all") %>%
    align(j = 3:12, align = "center", part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    fontsize(size = 8, part = "body") %>%
    fontsize(size = 9, part = "header") %>%
    bold(part = "header") %>%
    bold(j = 1, part = "body") %>%
    width(j = 1, width = 1.5) %>%
    width(j = 2, width = 1.0) %>%
    width(j = 3, width = 0.6) %>%
    width(j = 4, width = 1.0) %>%
    width(j = 5, width = 0.65) %>%
    width(j = 6, width = 0.6) %>%
    width(j = 7, width = 0.6) %>%
    width(j = 8, width = 1.0) %>%
    width(j = 9, width = 0.7) %>%
    width(j = 10, width = 1.0) %>%
    width(j = 11, width = 0.7) %>%
    width(j = 12, width = 0.65) %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1.5), part = "header") %>%
    hline(i = seq(6, nrow(formatted_table), 6), 
          border = fp_border(width = 1), part = "body") %>%
    vline(border = fp_border(width = 0.5, color = "#CCCCCC"), part = "all") %>%
    padding(padding = 2, part = "all") %>%
    hrule(rule = "exact")
  
  return(ft)
}

# ============================================
# GENERATE TABLES
# ============================================

cat("\n=== Generating Summary Tables ===\n")

# Table A
ft_a <- create_summary_table(
  table_a_data,
  "Table A: Mean Performance Across CVD Scores - Score-Specific Inputs"
)

doc_a <- read_docx() %>%
  body_add_par("", style = "Normal") %>%
  body_end_section_landscape() %>%
  body_add_flextable(value = ft_a)
print(doc_a, target = "Summary_Table_A_ScoreSpecific.docx")
save_as_html(ft_a, path = "Summary_Table_A_ScoreSpecific.html")

cat(sprintf("✓ Table A saved (%d rows)\n", nrow(table_a_data)))

# Table B
ft_b <- create_summary_table(
  table_b_data,
  "Table B: Mean Performance Across CVD Scores - General Datasets"
)

doc_b <- read_docx() %>%
  body_add_par("", style = "Normal") %>%
  body_end_section_landscape() %>%
  body_add_flextable(value = ft_b)
print(doc_b, target = "Summary_Table_B_GeneralDatasets.docx")
save_as_html(ft_b, path = "Summary_Table_B_GeneralDatasets.html")

cat(sprintf("✓ Table B saved (%d rows)\n", nrow(table_b_data)))

# Table C
ft_c <- create_summary_table(
  table_c_data,
  "Table C: Mean Performance Across CVD Scores - Score-Specific + Block Combinations"
)

doc_c <- read_docx() %>%
  body_add_par("", style = "Normal") %>%
  body_end_section_landscape() %>%
  body_add_flextable(value = ft_c)
print(doc_c, target = "Summary_Table_C_ScorePlusBlock.docx")
save_as_html(ft_c, path = "Summary_Table_C_ScorePlusBlock.html")

cat(sprintf("✓ Table C saved (%d rows)\n", nrow(table_c_data)))

# Combined
doc_combined <- read_docx()

doc_combined <- doc_combined %>%
  body_add_par("", style = "Normal") %>%
  body_end_section_landscape() %>%
  body_add_flextable(value = ft_a) %>%
  body_add_break() %>%
  body_add_par("", style = "Normal") %>%
  body_end_section_landscape() %>%
  body_add_flextable(value = ft_b) %>%
  body_add_break() %>%
  body_add_par("", style = "Normal") %>%
  body_end_section_landscape() %>%
  body_add_flextable(value = ft_c)

print(doc_combined, target = "Summary_Tables_All_Combined.docx")

cat("\n✓ All summary tables saved\n")









################################################################################
# Creating bar plot for R2 and MAE hold-out and CV data
################################################################################


# Load results
all_metrics <- qs2::qs_read("ml_comparison_compiled_objects.qs2")$metrics

# ============================================
# PREPARE DATA FOR SCORE-SPECIFIC ONLY
# ============================================

# Filter for score-specific datasets only
score_specific_data <- all_metrics %>%
  filter(
    dataset %in% c("QRISK3_risk_specific", "SCORE2_score_specific", 
                   "frs_10y_specific", "ascvd_10y_specific", "mean_risk_specific")
  )

# Define order and labels
outcome_order <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y", "mean_risk")
outcome_labels <- c(
  "QRISK3_risk" = "QRISK3",
  "SCORE2_score" = "SCORE2",
  "frs_10y" = "Framingham",
  "ascvd_10y" = "ASCVD",
  "mean_risk" = "Composite"
)

model_order <- c("elastic_net", "pls", "random_forest", "xgboost", "svm", "knn")
model_labels <- c(
  "elastic_net" = "Elastic Net",
  "pls" = "PLS",
  "random_forest" = "Random Forest",
  "xgboost" = "XGBoost",
  "svm" = "SVM-RBF",
  "knn" = "k-NN"
)

# Define color palette for models
model_colors <- c(
  "Elastic Net" = "#E69F00",
  "PLS" = "#56B4E9",
  "Random Forest" = "#009E73",
  "XGBoost" = "#F0E442",
  "SVM-RBF" = "#0072B2",
  "k-NN" = "#D55E00"
)

# Prepare data with proper ordering
plot_data <- score_specific_data %>%
  mutate(
    outcome = factor(outcome, levels = outcome_order),
    outcome_label = outcome_labels[as.character(outcome)],
    outcome_label = factor(outcome_label, levels = outcome_labels),
    model = factor(model, levels = model_order),
    model_label = model_labels[as.character(model)],
    model_label = factor(model_label, levels = model_labels)
  )

# ============================================
# PLOT 1: CV R²
# ============================================

p1_cv_rsq <- plot_data %>%
  ggplot(aes(x = cv_rsq, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(
    aes(xmin = cv_rsq - cv_rsq_sd, xmax = cv_rsq + cv_rsq_sd),
    position = position_dodge(width = 0.8),
    height = 0.3,
    linewidth = 0.4
  ) +
  geom_text(
    aes(x = cv_rsq + cv_rsq_sd, label = sprintf("%.2f", cv_rsq)),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 2
  ) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(
    limits = c(0, 1.0),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "Cross-Validation",
    x = expression(R^2),
    y = NULL,
    fill = "Model"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = element_text(size = 9),
    axis.title.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# ============================================
# PLOT 2: Holdout R²
# ============================================

p2_holdout_rsq <- plot_data %>%
  filter(!is.na(holdout_rsq)) %>%
  ggplot(aes(x = holdout_rsq, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(
    aes(label = sprintf("%.2f", holdout_rsq)),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 2
  ) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(
    limits = c(0, 1.0),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "Holdout Set",
    x = expression(R^2),
    y = NULL,
    fill = "Model"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = element_text(size = 9),
    axis.title.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold")
  )

# ============================================
# PLOT 3: CV MAE
# ============================================

p3_cv_mae <- plot_data %>%
  ggplot(aes(x = cv_mae, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(
    aes(xmin = cv_mae - cv_mae_sd, xmax = cv_mae + cv_mae_sd),
    position = position_dodge(width = 0.8),
    height = 0.3,
    linewidth = 0.4
  ) +
  geom_text(
    aes(x = cv_mae + cv_mae_sd, label = sprintf("%.2f", cv_mae)),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 2
  ) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "",
    x = "MAE",
    y = NULL,
    fill = "Model"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 9),
    axis.title.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# ============================================
# PLOT 4: Holdout MAE
# ============================================

p4_holdout_mae <- plot_data %>%
  filter(!is.na(holdout_mae)) %>%
  ggplot(aes(x = holdout_mae, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(
    aes(label = sprintf("%.2f", holdout_mae)),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 2
  ) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(
    expand = expansion(mult = c( 0, 0.15))
  ) +
  labs(
    title = "",
    x = "MAE",
    y = NULL,
    fill = "Model"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 9),
    axis.title.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold")
  )

# ============================================
# COMBINE WITH PATCHWORK
# ============================================

combined_plot <- (p1_cv_rsq | p2_holdout_rsq) / 
  (p3_cv_mae | p4_holdout_mae) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Model Performance Comparison: Score-Specific Inputs",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    )
  )

# Display
print(combined_plot)

# ============================================
# SAVE
# ============================================

# Save as high-resolution image
ggsave(
  "model_comparison_score_specific.png",
  plot = combined_plot,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  "model_comparison_score_specific.pdf",
  plot = combined_plot,
  width = 12,
  height = 8,
  device = "pdf"
)

cat("\n✓ Plot saved as:\n")
cat("  - model_comparison_score_specific.png\n")
cat("  - model_comparison_score_specific.pdf\n")