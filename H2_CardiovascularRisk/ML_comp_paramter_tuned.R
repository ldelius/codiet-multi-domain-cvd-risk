################################################################################
# CoDiet CVD Prediction Model Comparison
# Author: Luisa Delius

# This script compares multiple ML algorithms (Elastic Net, PLS, Random Forest,
# XGBoost, SVM-RBF, k-NN) for predicting CVD risk scores using various datasets
# from the CoDiet cohort.
################################################################################

# ============================================
# 1. Setup
# ============================================

# ---- 1.1 Load packages ----
library(tidyverse)
library(tidymodels)
library(vip)
library(future)
library(qs2)
library(plsmod)
library(writexl)
library(flextable)
library(officer)
library(patchwork)

# ---- 1.2 Global settings ----
set.seed(42)
options(scipen = 999, expressions = 500000)

CPUS <- parallel::detectCores() - 3
cat("Using", CPUS, "CPU cores for parallel processing\n")

select <- dplyr::select
slice <- dplyr::slice
rename <- dplyr::rename
filter <- dplyr::filter

wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files/processed_data"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

# Create output folder for all results
base_name <- "ml_comparison_results_v_"
existing_dirs <- list.dirs(wkdir, full.names = FALSE, recursive = FALSE)
existing_versions <- existing_dirs[grepl(paste0("^", base_name), existing_dirs)]

if (length(existing_versions) == 0) {
  version_num <- 1
} else {
  # Extract version numbers and find max
  version_nums <- as.numeric(sub(paste0(base_name, "(\\d+)"), "\\1", existing_versions))
  version_num <- max(version_nums, na.rm = TRUE) + 1
}

# Create new versioned folder
output_dir <- file.path(wkdir, paste0(base_name, sprintf("%02d", version_num)))
dir.create(output_dir)
cat("Output directory:", output_dir, "(version", version_num, ")\n")

# ============================================
# 2. Load and prepare data
# ============================================

# ---- 2.1 Load raw datasets ----
df_fatty_acids <- readRDS("df_fatty_acids_predictor_statin_suppl.rds") %>%
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

df_fatty_acids_statin_suppl <- readRDS("df_fatty_acids_predictor_statin_suppl.rds")

df_risk_factors <- readRDS("df_risk_factor_predictors.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk,
         -Heart.Rate, -Age, -Total.Cholesterol.mg.dl, -HDL.mg.dl,
         -Body.Weight, -Height, -BMI, -Gender,
         -stress_resilience_status, -stress_index_status, -Age.Risk,
         -Body.fat, -Fat.free.mass) %>%
  full_join(df_fatty_acids_statin_suppl %>% select(Sample_ID, Statins, Supplements), 
            by = "Sample_ID")

# Add Body.fat and Fat.free.mass to body composition
df_risk_factors_raw <- readRDS("df_risk_factor_predictors.rds")
body_comp_extras <- df_risk_factors_raw %>%
  select(Sample_ID, Body.fat, Fat.free.mass) %>%
  filter(!is.na(Sample_ID))

df_body_composition <- df_body_composition %>%
  left_join(body_comp_extras, by = "Sample_ID")

df_all_cvd_risk_scores <- readRDS("df_all_risk_scores.rds") %>%
  rename(Sample_ID = PatientID) %>%
  select(-SCORE2_strat)

df_ascvd_frs_score2_input <- readRDS("ASCVD_SCORE2_Framingham_input.rds") %>%
  rename(Sample_ID = PatientID)

df_QRISK3_input <- readRDS("QRISK3_calculation_input.rds") %>%
  rename(Sample_ID = PatientID)

CVD_scores <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y", "mean_risk")

# ---- 2.2 Create predictor block datasets ----
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

datasets_list <- list(
  all_data = df_all_data_model,
  fatty_acids = df_fatty_acids_model,
  lipids = df_lipidomics_model,
  urine_NMR = df_urine_nmr_model,
  body_composition = df_body_composition_model,
  sociodemographics_lifestyle = df_REDcap_demographics_model,
  clinical_risk_factors = df_risk_factors_model
)

# ---- 2.3 Score-specific input datasets ----
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

# ============================================
# 3. Helper functions
# ============================================

# Removes CVD score columns from a dataframe, returning only Sample_ID and predictor variables
prep_extra_predictors <- function(extra_df) {
  keep_cols <- setdiff(names(extra_df), CVD_scores)
  extra_df[, keep_cols, drop = FALSE]
}

# Joins score-specific input variables with an additional predictor block
build_score_plus_block <- function(score, block_name, score_specific_datasets, datasets_list) {
  base_df <- score_specific_datasets[[score]]
  extra_df <- datasets_list[[block_name]]
  extra_predictors <- prep_extra_predictors(extra_df)
  dplyr::left_join(base_df, extra_predictors, by = "Sample_ID")
}

# ============================================
# 4. Select hold-out participants
# ============================================
# Creates a representative hold-out set (n=20) for final model evaluation.
# 
# Priorities (in order):
# 1. CVD risk distribution - match cohort's risk tertile proportions
# 2. Sex distribution - match cohort's overall sex ratio
# 3. Low missingness - prefer participants with complete data
# 4. Country representation - ensure at least 1 per country
#
# Method: Stratified sampling within risk tertiles, then adjust for sex/country

cat("\n=== SELECTING HOLD-OUT PARTICIPANTS ===\n")

# Use all_data as reference for holdout selection (covers all participants)
reference_df <- df_all_data_model

participant_info <- reference_df %>%
  select(Sample_ID) %>%
  mutate(
    missing_pct = rowMeans(is.na(reference_df %>% select(-Sample_ID, -all_of(CVD_scores)))) * 100
  ) %>%
  left_join(df_all_cvd_risk_scores %>% select(Sample_ID, mean_risk), by = "Sample_ID") %>%
  left_join(df_QRISK3_input %>% select(Sample_ID, Sex), by = "Sample_ID") %>%
  left_join(df_risk_factors %>% select(Sample_ID, Country), by = "Sample_ID") %>%
  filter(!is.na(mean_risk), !is.na(Sex), !is.na(Country)) %>%
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

cat(sprintf("Total eligible participants: %d\n", nrow(participant_info)))

# Cohort distributions
country_dist <- table(participant_info$Country)
sex_dist <- table(participant_info$Sex)
risk_dist <- table(participant_info$risk_category)

set.seed(42)
n_holdout <- 20

# Calculate target numbers per stratum
target_risk <- round(prop.table(risk_dist) * n_holdout)
target_sex <- round(prop.table(sex_dist) * n_holdout)

cat(sprintf("Target risk distribution: Low=%d, Medium=%d, High=%d\n", 
            target_risk["Low"], target_risk["Medium"], target_risk["High"]))
cat(sprintf("Target sex distribution: Female=%d, Male=%d\n", 
            target_sex["Female"], target_sex["Male"]))

# Stratified selection: sample from each risk tertile proportionally
# Within each tertile, prioritize by: (1) low missingness, (2) sex balance needed
holdout_ids <- c()

for (risk_cat in c("Low", "Medium", "High")) {
  n_from_tertile <- target_risk[risk_cat]
  
  tertile_candidates <- participant_info %>%
    filter(risk_category == risk_cat) %>%
    arrange(missing_pct)  # Prioritize low missingness
  
  # Check current sex balance and prioritize underrepresented sex
  if (length(holdout_ids) > 0) {
    current_sex <- table(factor(participant_info$Sex[participant_info$Sample_ID %in% holdout_ids], 
                                levels = c("Female", "Male")))
    sex_need <- target_sex - current_sex
    
    # Sort by: underrepresented sex first, then missingness
    tertile_candidates <- tertile_candidates %>%
      mutate(sex_priority = case_when(
        Sex == "Female" & sex_need["Female"] > 0 ~ 1,
        Sex == "Male" & sex_need["Male"] > 0 ~ 1,
        TRUE ~ 2
      )) %>%
      arrange(sex_priority, missing_pct) %>%
      select(-sex_priority)
  }
  
  selected <- tertile_candidates %>%
    filter(!Sample_ID %in% holdout_ids) %>%
    slice_head(n = n_from_tertile)
  
  holdout_ids <- c(holdout_ids, selected$Sample_ID)
}

# Ensure at least 1 participant per country (swap if needed)
for (ctry in names(country_dist)) {
  if (!any(participant_info$Country[participant_info$Sample_ID %in% holdout_ids] == ctry)) {
    # Find best candidate from this country (lowest missingness)
    country_candidate <- participant_info %>%
      filter(Country == ctry, !Sample_ID %in% holdout_ids) %>%
      arrange(missing_pct) %>%
      slice_head(n = 1)
    
    if (nrow(country_candidate) > 0) {
      # Swap with highest-missingness participant from same risk category
      swap_out <- participant_info %>%
        filter(Sample_ID %in% holdout_ids, 
               risk_category == country_candidate$risk_category) %>%
        arrange(desc(missing_pct)) %>%
        slice_head(n = 1)
      
      if (nrow(swap_out) > 0) {
        holdout_ids <- setdiff(holdout_ids, swap_out$Sample_ID)
        holdout_ids <- c(holdout_ids, country_candidate$Sample_ID)
        cat(sprintf("  Swapped %s for %s to include country %s\n", 
                    swap_out$Sample_ID, country_candidate$Sample_ID, ctry))
      }
    }
  }
}

HOLDOUT_IDS <- holdout_ids
holdout_participants <- participant_info %>% filter(Sample_ID %in% HOLDOUT_IDS)

# Report final distribution
cat(sprintf("\nHold-out set: %d participants\n", nrow(holdout_participants)))
cat(sprintf("  Mean missingness: %.1f%%\n", mean(holdout_participants$missing_pct)))
cat(sprintf("  Risk: Low=%d, Medium=%d, High=%d\n", 
            sum(holdout_participants$risk_category == "Low"),
            sum(holdout_participants$risk_category == "Medium"),
            sum(holdout_participants$risk_category == "High")))
cat(sprintf("  Sex: Female=%d, Male=%d\n",
            sum(holdout_participants$Sex == "Female"),
            sum(holdout_participants$Sex == "Male")))
cat(sprintf("  Countries represented: %d/%d\n",
            n_distinct(holdout_participants$Country), length(country_dist)))

qs2::qs_save(HOLDOUT_IDS, file.path(output_dir, "holdout_participant_ids.qs2"))

# ============================================
# 5. Apply column missingness exclusion
# ============================================

MISSINGNESS_THRESHOLD <- 0.40

remove_high_missing_cols <- function(df, threshold = MISSINGNESS_THRESHOLD,
                                     exclude_cols = c("Sample_ID", CVD_scores)) {
  
  cols_to_check <- setdiff(names(df), exclude_cols)
  if (length(cols_to_check) == 0) return(df)
  
  missing_pct <- sapply(df[cols_to_check], function(x) mean(is.na(x)))
  cols_to_remove <- names(missing_pct)[missing_pct >= threshold]
  
  if (length(cols_to_remove) > 0) {
    cat(sprintf("  Removing %d columns with >%.0f%% missing\n", 
                length(cols_to_remove), threshold * 100))
    df <- df %>% select(-all_of(cols_to_remove))
  }
  df
}

# Apply to all dataset lists
datasets_list <- lapply(datasets_list, remove_high_missing_cols)
score_specific_datasets <- lapply(score_specific_datasets, remove_high_missing_cols)

# ============================================
# 6. Data sanitisation
# ============================================
# Convert character columns to factors (required by recipes for nominal predictors).
# Unexpected types (dates, logicals) are NOT silently converted — they will cause
# an informative error during model fitting if present.

sanitize_df <- function(df) {
  for (col in names(df)) {
    if (col %in% c("Sample_ID", CVD_scores)) next
    x <- df[[col]]
    if (!is.numeric(x) && !is.factor(x)) df[[col]] <- as.factor(x)
  }
  df
}

datasets_list <- lapply(datasets_list, sanitize_df)
score_specific_datasets <- lapply(score_specific_datasets, sanitize_df)

# ============================================
# 7. Create Score-specific + Block combinations
# ============================================

cat("\n=== CREATING SCORE-SPECIFIC + BLOCK COMBINATIONS ===\n")

score_plus_block_datasets <- list()

for (score in names(score_specific_datasets)) {
  for (block in names(datasets_list)) {
    combo_name <- paste0(score, "_plus_", block)
    score_plus_block_datasets[[combo_name]] <- build_score_plus_block(
      score, block, score_specific_datasets, datasets_list
    )
  }
}

cat(sprintf("Created %d combinations\n", length(score_plus_block_datasets)))

# Verify holdout participants present
cat("\n--- Holdout verification ---\n")
for (ds_name in names(datasets_list)) {
  n_present <- sum(HOLDOUT_IDS %in% datasets_list[[ds_name]]$Sample_ID)
  cat(sprintf("  %s: %d/%d\n", ds_name, n_present, length(HOLDOUT_IDS)))
}

# ============================================
# 8. CV and tuning settings
# ============================================

CV_FOLDS <- 5
CV_REPEATS <- 5
GRID_SIZE <- 20
N_IMPORTANCE_SIMS <- 3

# Overfitting threshold: absolute difference between train and holdout Q²
# E.g., 0.25 means train R²=0.80 and holdout Q²=0.55 is acceptable (gap=0.25)
RSQ_GAP_THRESHOLD <- 0.25

cat("\n=== SETTINGS ===\n")
cat(sprintf("CV: %d-fold x %d repeats\n", CV_FOLDS, CV_REPEATS))
cat(sprintf("Grid size: %d\n", GRID_SIZE))
cat(sprintf("Q² gap threshold: %.2f (absolute difference)\n", RSQ_GAP_THRESHOLD))

# ============================================
# 9. Model specifications
# ============================================

elastic_net_spec <- linear_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>%
  set_mode("regression")

sparse_pls_spec <- pls(num_comp = tune(), predictor_prop = tune()) %>%
  set_engine("mixOmics") %>%
  set_mode("regression")

# Override default num_comp range (default caps at 4)
sparse_pls_num_comp_range <- num_comp(range = c(1L, 10L))

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
  sparse_pls = sparse_pls_spec,
  random_forest = rf_spec,
  xgboost = xgb_spec,
  svm = svm_spec,
  knn = knn_spec
)

# ============================================
# 10. Best model selection function
# ============================================
# Selects best model based on holdout MAE, excluding models with excessive overfitting.
# rsq_gap = train_rsq - holdout Q² (absolute difference, not ratio)

select_best_model <- function(metrics_df, rsq_gap_threshold = RSQ_GAP_THRESHOLD) {
  candidates <- metrics_df %>%
    filter(!is.na(holdout_mae), !is.na(rsq_gap))
  
  # Fallback to CV if no holdout metrics
  if (nrow(candidates) == 0) {
    return(metrics_df %>% slice_min(cv_rmse, n = 1) %>% pull(model) %>% `[`(1))
  }
  
  # Exclude overfit models (R² gap > threshold)
  reasonable <- candidates %>% filter(rsq_gap <= rsq_gap_threshold)
  
  # If all models exceed threshold, pick the one with smallest gap
  if (nrow(reasonable) == 0) {
    return(candidates %>% slice_min(rsq_gap, n = 1) %>% pull(model) %>% `[`(1))
  }
  
  # Among acceptable models, pick best holdout MAE
  reasonable %>% slice_min(holdout_mae, n = 1) %>% pull(model) %>% `[`(1)
}

# ============================================
# 11. Main modelling function
# ============================================

run_models_for_dataset <- function(df, outcome, dataset_name, 
                                   model_specs, cv_folds = 5, cv_repeats = 3,
                                   grid_size = 10, seed = 42,
                                   calculate_importance = FALSE,
                                   n_importance_sims = 3,
                                   holdout_ids = NULL) {
  set.seed(seed)
  cvd_scores_local <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y", "mean_risk")
  
  if (!outcome %in% names(df)) {
    warning(sprintf("Outcome '%s' not found in '%s'", outcome, dataset_name))
    return(NULL)
  }
  
  # Prepare data
  other_scores <- setdiff(cvd_scores_local, outcome)
  df_full <- df %>% 
    select(-any_of(other_scores)) %>%
    filter(!is.na(.data[[outcome]]))
  
  # Split into train and holdout
  if (!is.null(holdout_ids) && length(holdout_ids) > 0) {
    df_train <- df_full %>% filter(!Sample_ID %in% holdout_ids)
    holdout_clean <- df_full %>% filter(Sample_ID %in% holdout_ids)
    n_holdout <- nrow(holdout_clean)
    cat(sprintf("    Split: %d train, %d holdout\n", nrow(df_train), n_holdout))
  } else {
    df_train <- df_full
    holdout_clean <- NULL
    n_holdout <- NA
  }
  
  # Per-dataset >80% row missingness exclusion
  predictor_cols <- setdiff(names(df_train), c("Sample_ID", outcome))
  if (length(predictor_cols) > 0) {
    train_row_missing <- rowMeans(is.na(df_train[predictor_cols])) * 100
    high_missing <- train_row_missing > 80
    
    if (any(high_missing)) {
      cat(sprintf("    Removing %d participants with >80%% missing\n", sum(high_missing)))
      df_train <- df_train[!high_missing, ]
    }
    
    if (!is.null(holdout_clean) && nrow(holdout_clean) > 0) {
      holdout_missing <- rowMeans(is.na(holdout_clean[predictor_cols])) * 100
      holdout_clean <- holdout_clean[holdout_missing <= 80, ]
    }
  }
  
  n_train <- nrow(df_train)
  if (n_train < 30) {
    warning(sprintf("Too few observations (%d) for %s ~ %s", n_train, outcome, dataset_name))
    return(NULL)
  }
  
  n_predictors <- ncol(df_train) - 2
  cat(sprintf("  %s ~ %s: n=%d, p=%d\n", outcome, dataset_name, n_train, n_predictors))
  
  # CV setup
  cv_splits <- vfold_cv(df_train, v = cv_folds, repeats = cv_repeats, strata = !!sym(outcome))
  
  # Recipe
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
  
  # Fit all models
  wf_set <- workflow_set(preproc = list(recipe = rec), models = model_specs, cross = TRUE)
  
  # Set default grid size for all models
  wf_set <- wf_set %>%
    option_add(grid = grid_size)
  
  # Override sparse PLS num_comp range (extend beyond default of 4)
  max_comp <- min(8L, max(4L, n_predictors - 4L))
  spls_params <- wf_set %>%
    extract_parameter_set_dials("recipe_sparse_pls") %>%
    update(num_comp = num_comp(range = c(1L, max_comp)))
  wf_set <- wf_set %>%
    option_remove(id = "recipe_sparse_pls") %>%
    option_add(param_info = spls_params, id = "recipe_sparse_pls")
  cat(sprintf("    Sparse PLS: num_comp range 1-%d\n", max_comp))
  
  # Override SVM with extended cost range (default caps at ~32, extend to ~1000)
  svm_grid <- grid_latin_hypercube(
    cost(range = c(-10, 10)),
    rbf_sigma(range = c(-10, 0)),
    size = grid_size
  )
  wf_set <- wf_set %>%
    option_remove(id = "recipe_svm") %>%
    option_add(grid = svm_grid, id = "recipe_svm")
  
  wf_results <- wf_set %>%
    workflow_map(
      fn = "tune_grid",
      resamples = cv_splits,
      metrics = metric_set(mae, rmse, rsq),
      control = control_grid(save_pred = TRUE, save_workflow = FALSE, 
                             verbose = FALSE, parallel_over = "everything"),
      seed = seed,
      verbose = TRUE
    )
  
  # Extract metrics for all models
  all_model_metrics <- tibble()
  
  # One-SE rule helper: select simplest model within 1 SE of best RMSE.
  # Each model type has a different "simplicity" direction.
  select_one_se <- function(tune_result, wf_id) {
    model_type <- str_remove(wf_id, "recipe_")
    tryCatch({
      switch(model_type,
             elastic_net = tune_result %>% select_by_one_std_err(metric = "rmse", desc(penalty)),
             sparse_pls  = tune_result %>% select_by_one_std_err(metric = "rmse", num_comp),
             random_forest = tune_result %>% select_by_one_std_err(metric = "rmse", min_n),
             xgboost     = tune_result %>% select_by_one_std_err(metric = "rmse", min_n),
             svm         = tune_result %>% select_by_one_std_err(metric = "rmse", cost),
             knn         = tune_result %>% select_by_one_std_err(metric = "rmse", desc(neighbors)),
             tune_result %>% select_best(metric = "rmse")
      )
    }, error = function(e) {
      tune_result %>% select_best(metric = "rmse")
    })
  }
  
  for (model_name in wf_results$wflow_id) {
    best_result <- wf_results %>% extract_workflow_set_result(model_name)
    best_params <- select_one_se(best_result, model_name)
    
    cv_per_fold <- best_result %>%
      collect_metrics(summarize = FALSE) %>%
      filter(.config == best_params$.config)
    
    cv_summary <- cv_per_fold %>%
      group_by(.metric) %>%
      summarise(cv_mean = mean(.estimate, na.rm = TRUE),
                cv_sd = sd(.estimate, na.rm = TRUE), .groups = "drop")
    
    # Fit final model
    final_wf <- wf_results %>%
      extract_workflow(model_name) %>%
      finalize_workflow(best_params)
    
    final_fit <- tryCatch(fit(final_wf, data = df_train), error = function(e) NULL)
    
    # Training metrics
    train_mae <- train_rmse <- train_rsq <- train_cal_slope <- train_cal_intercept <- NA
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
    
    # Holdout metrics (Q²)
    q2 <- holdout_mae <- holdout_rmse <- holdout_cal_slope <- holdout_cal_intercept <- NA
    if (!is.null(holdout_clean) && nrow(holdout_clean) >= 5 && !is.null(final_fit)) {
      holdout_preds <- tryCatch({
        predict(final_fit, new_data = holdout_clean) %>%
          bind_cols(holdout_clean %>% select(all_of(outcome)))
      }, error = function(e) NULL)
      
      if (!is.null(holdout_preds) && nrow(holdout_preds) > 0) {
        holdout_mae <- mean(abs(holdout_preds$.pred - holdout_preds[[outcome]]), na.rm = TRUE)
        holdout_rmse <- sqrt(mean((holdout_preds$.pred - holdout_preds[[outcome]])^2, na.rm = TRUE))
        
        if (sd(holdout_preds$.pred, na.rm = TRUE) > 0.001) {
          q2 <- cor(holdout_preds$.pred, holdout_preds[[outcome]], use = "complete.obs")^2
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
    
    # Overfitting indicators
    rsq_gap <- train_rsq - q2
    rsq_ratio <- ifelse(!is.na(q2) & q2 > 0.01, train_rsq / q2, NA)
    
    model_metrics <- tibble(
      dataset = dataset_name,
      outcome = outcome,
      model = str_remove(model_name, "recipe_"),
      n_train = n_train,
      n_holdout = n_holdout,
      n_predictors = n_predictors,
      train_mae = round(train_mae, 3),
      train_rmse = round(train_rmse, 3),
      train_rsq = round(train_rsq, 3),
      train_cal_slope = round(train_cal_slope, 3),
      train_cal_intercept = round(train_cal_intercept, 3),
      cv_mae = round(cv_mae, 3),
      cv_mae_sd = round(cv_mae_sd, 3),
      cv_rmse = round(cv_rmse, 3),
      cv_rmse_sd = round(cv_rmse_sd, 3),
      cv_rsq = round(cv_rsq, 3),
      cv_rsq_sd = round(cv_rsq_sd, 3),
      holdout_mae = round(holdout_mae, 3),
      holdout_rmse = round(holdout_rmse, 3),
      q2 = round(q2, 3),
      holdout_cal_slope = round(holdout_cal_slope, 3),
      holdout_cal_intercept = round(holdout_cal_intercept, 3),
      rsq_gap = round(rsq_gap, 3),
      rsq_ratio = round(rsq_ratio, 3),
      overfit_flag = case_when(
        is.na(q2) ~ NA,
        rsq_gap > RSQ_GAP_THRESHOLD ~ TRUE,
        TRUE ~ FALSE
      )
    )
    
    all_model_metrics <- bind_rows(all_model_metrics, model_metrics)
  }
  
  # Select best model
  best_model_name <- select_best_model(all_model_metrics)
  cat(sprintf("    Best model: %s\n", best_model_name))
  
  # Variable importance (optional)
  importance_results <- NULL
  if (calculate_importance) {
    best_wf_id <- paste0("recipe_", best_model_name)
    best_params <- select_one_se(
      wf_results %>% extract_workflow_set_result(best_wf_id), best_wf_id
    )
    
    final_wf <- wf_results %>%
      extract_workflow(paste0("recipe_", best_model_name)) %>%
      finalize_workflow(best_params)
    
    final_fit <- fit(final_wf, data = df_train)
    
    importance_results <- tryCatch({
      vi_permute(
        object = final_fit, train = df_train, target = outcome,
        metric = "rmse", nsim = n_importance_sims,
        pred_wrapper = function(object, newdata) predict(object, new_data = newdata)$.pred
      ) %>%
        mutate(dataset = dataset_name, outcome = outcome, best_model = best_model_name) %>%
        filter(Importance > 0)
    }, error = function(e) NULL)
  }
  
  # CV and holdout predictions for best model
  best_result <- wf_results %>% extract_workflow_set_result(paste0("recipe_", best_model_name))
  best_config <- select_one_se(best_result, paste0("recipe_", best_model_name))
  
  cv_predictions <- tryCatch({
    best_result %>%
      collect_predictions() %>%
      filter(.config == best_config$.config) %>%
      mutate(dataset = dataset_name, outcome_name = outcome, best_model = best_model_name)
  }, error = function(e) NULL)
  
  holdout_predictions <- NULL
  if (!is.null(holdout_clean) && nrow(holdout_clean) >= 5) {
    final_wf <- wf_results %>%
      extract_workflow(paste0("recipe_", best_model_name)) %>%
      finalize_workflow(best_config)
    final_fit <- fit(final_wf, data = df_train)
    
    holdout_predictions <- tryCatch({
      predict(final_fit, new_data = holdout_clean) %>%
        bind_cols(holdout_clean %>% select(Sample_ID, all_of(outcome))) %>%
        mutate(dataset = dataset_name, outcome_name = outcome, best_model = best_model_name)
    }, error = function(e) NULL)
  }
  
  list(
    metrics = all_model_metrics,
    importance = importance_results,
    cv_predictions = cv_predictions,
    holdout_predictions = holdout_predictions,
    wf_results = wf_results,
    best_model = best_model_name
  )
}

# ============================================
# UPDATED SECTIONS 12–19
# Changes:
#   - Composite/mean_risk excluded
#   - "all_data" renamed to "full_predictor_set" throughout
#   - Combined table: one table, two score pairs, two dataset rows each
#   - Tuning plots: per-dataset overview + per-score detailed
# ============================================

# ============================================
# 12. Build task list (score-specific + full predictor set, no composite)
# ============================================

cat("\n=== BUILDING TASK LIST ===\n")

CVD_scores_final <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y")

task_list <- list()
task_counter <- 1

# A) Score-specific
for (score_name in CVD_scores_final) {
  task_list[[task_counter]] <- list(
    df = score_specific_datasets[[score_name]],
    holdout_ids = HOLDOUT_IDS,
    outcome = score_name,
    dataset_name = paste0(score_name, "_specific"),
    category = "score_specific"
  )
  task_counter <- task_counter + 1
}

# B) Full predictor set (renamed from all_data)
for (outcome in CVD_scores_final) {
  task_list[[task_counter]] <- list(
    df = datasets_list[["all_data"]],
    holdout_ids = HOLDOUT_IDS,
    outcome = outcome,
    dataset_name = "full_predictor_set",
    category = "predictor_block"
  )
  task_counter <- task_counter + 1
}

cat(sprintf("Total tasks: %d\n", length(task_list)))

# ============================================
# 13. Test run
# ============================================

# cat("\n=== TEST RUN ===\n")
# plan(multisession, workers = CPUS)
# 
# test_result <- run_models_for_dataset(
#   df = score_specific_datasets[["QRISK3_risk"]],
#   outcome = "QRISK3_risk",
#   dataset_name = "QRISK3_risk_specific",
#   model_specs = model_specs,
#   cv_folds = CV_FOLDS,
#   cv_repeats = CV_REPEATS,
#   grid_size = GRID_SIZE,
#   seed = 42,
#   calculate_importance = FALSE,
#   holdout_ids = HOLDOUT_IDS
# )
# 
# plan(sequential)
# 
# if (!is.null(test_result)) {
#   cat("Test successful. Best model:", test_result$best_model, "\n")
#   print(test_result$metrics %>% select(model, train_rsq, cv_rsq, q2, rsq_gap))
# }

# ============================================
# 14. Full run
# ============================================

plan(multisession, workers = CPUS)

checkpoint_file <- file.path(output_dir, "ml_comparison_checkpoint.qs2")
if (file.exists(checkpoint_file)) {
  checkpoint <- qs2::qs_read(checkpoint_file)
  all_results <- checkpoint$all_results
  failed_tasks <- checkpoint$failed_tasks
  start_idx <- checkpoint$last_completed + 1
  cat(sprintf("Resuming from checkpoint: task %d/%d\n", start_idx, length(task_list)))
} else {
  all_results <- list()
  failed_tasks <- c()
  start_idx <- 1
}

start_time <- Sys.time()

for (i in start_idx:length(task_list)) {
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
      holdout_ids = task$holdout_ids
    )
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    return(NULL)
  })
  
  if (!is.null(result)) {
    all_results[[task_name]] <- result
    cat("  Done\n")
  } else {
    failed_tasks <- c(failed_tasks, task_name)
    cat("  Failed\n")
  }
  
  if (i %% 5 == 0 || i == length(task_list)) {
    qs2::qs_save(
      list(all_results = all_results, failed_tasks = failed_tasks, last_completed = i),
      checkpoint_file
    )
    cat(sprintf("  [Checkpoint saved: %d/%d tasks]\n", i, length(task_list)))
  }
}

plan(sequential)

total_time <- round(difftime(Sys.time(), start_time, units = "hours"), 2)
cat(sprintf("\nCompleted in %.2f hours\n", total_time))
cat(sprintf("Successful: %d | Failed: %d\n", length(all_results), length(failed_tasks)))

qs2::qs_save(all_results, file.path(output_dir, "ml_comparison_all_results.qs2"))
if (file.exists(checkpoint_file)) file.remove(checkpoint_file)
if (length(failed_tasks) > 0) writeLines(failed_tasks, file.path(output_dir, "ml_comparison_failed_tasks.txt"))

# ============================================
# 15. Compile results
# ============================================

all_metrics <- bind_rows(lapply(all_results, function(x) x$metrics))
all_importance <- bind_rows(lapply(all_results, function(x) x$importance))
all_cv_preds <- bind_rows(lapply(all_results, function(x) x$cv_predictions))
all_holdout_preds <- bind_rows(lapply(all_results, function(x) x$holdout_predictions))

best_models_summary <- map_dfr(names(all_results), function(task_name) {
  result <- all_results[[task_name]]
  result$metrics %>% filter(model == result$best_model)
})

model_rankings <- all_metrics %>%
  filter(!is.na(holdout_mae)) %>%
  group_by(dataset, outcome) %>%
  mutate(rank = rank(holdout_mae)) %>%
  ungroup() %>%
  group_by(model) %>%
  summarise(
    mean_rank = round(mean(rank), 2),
    times_best = sum(rank == 1),
    mean_q2 = round(mean(q2, na.rm = TRUE), 3),
    mean_holdout_mae = round(mean(holdout_mae, na.rm = TRUE), 3),
    mean_rsq_gap = round(mean(rsq_gap, na.rm = TRUE), 3),
    pct_overfit = round(mean(overfit_flag, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(mean_rank)

cat("\n--- Model Rankings ---\n")
print(model_rankings)

# Save
excel_output <- list(
  Best_Models = best_models_summary,
  All_Metrics = all_metrics,
  Model_Rankings = model_rankings,
  Variable_Importance = all_importance
)
write_xlsx(excel_output, file.path(output_dir, "ml_model_comparison_results.xlsx"))

qs2::qs_save(list(
  metrics = all_metrics,
  importance = all_importance,
  cv_preds = all_cv_preds,
  holdout_preds = all_holdout_preds,
  best_models = best_models_summary,
  model_rankings = model_rankings
), file.path(output_dir, "ml_comparison_compiled_objects.qs2"))

# ============================================
# 16. Combined comparison table
# ============================================
# One table with:
#   Left pair: QRISK3 | SCORE2    (side by side)
#   Right pair: Framingham | ASCVD (below)
# Each pair shows score-specific rows then full predictor set rows

cat("\n=== GENERATING COMBINED TABLE ===\n")

model_order <- c("elastic_net", "sparse_pls", "random_forest", "xgboost", "svm", "knn")
model_labels_table <- c(
  elastic_net = "Elastic Net", sparse_pls = "Sparse PLS",
  random_forest = "RF", xgboost = "XGBoost",
  svm = "SVM-RBF", knn = "k-NN"
)

outcome_labels <- c(
  QRISK3_risk = "QRISK3", SCORE2_score = "SCORE2",
  frs_10y = "Framingham 10-year", ascvd_10y = "ASCVD 10-year"
)

dataset_display <- c(
  "score_specific" = "Score-specific inputs",
  "full_predictor_set" = "Full predictor set"
)

# Helper: format one block
format_block <- function(metrics_df, score, dataset_name) {
  metrics_df %>%
    filter(outcome == score, dataset == dataset_name) %>%
    mutate(model = factor(model, levels = model_order)) %>%
    arrange(model) %>%
    transmute(
      Model = model_labels_table[as.character(model)],
      `CV R²` = sprintf("%.3f ± %.3f", cv_rsq, cv_rsq_sd),
      `Q²` = ifelse(!is.na(q2), sprintf("%.3f", q2), "–"),
      `R²-Q² Gap` = ifelse(!is.na(rsq_gap), sprintf("%.3f", rsq_gap), "–"),
      `Holdout MAE` = ifelse(!is.na(holdout_mae), sprintf("%.3f", holdout_mae), "–")
    )
}

# Score pairs
score_pairs <- list(
  c("QRISK3_risk", "SCORE2_score"),
  c("frs_10y", "ascvd_10y")
)

# Build all body rows with label rows
all_rows <- tibble()
label_row_indices <- c()
dataset_label_indices <- c()

for (pair in score_pairs) {
  left_score <- pair[1]
  right_score <- pair[2]
  
  # Score label row
  label_row <- tibble(
    Model_L = outcome_labels[left_score],
    CVR2_L = "", Q2_L = "", Gap_L = "", MAE_L = "",
    Spacer = "",
    Model_R = outcome_labels[right_score],
    CVR2_R = "", Q2_R = "", Gap_R = "", MAE_R = ""
  )
  label_row_indices <- c(label_row_indices, nrow(all_rows) + 1)
  all_rows <- bind_rows(all_rows, label_row)
  
  # For each dataset (score-specific first, then full predictor set)
  for (ds_type in c("score_specific", "full_predictor_set")) {
    
    # Dataset sub-label row
    ds_label <- dataset_display[ds_type]
    ds_label_row <- tibble(
      Model_L = ds_label,
      CVR2_L = "", Q2_L = "", Gap_L = "", MAE_L = "",
      Spacer = "",
      Model_R = ds_label,
      CVR2_R = "", Q2_R = "", Gap_R = "", MAE_R = ""
    )
    dataset_label_indices <- c(dataset_label_indices, nrow(all_rows) + 1)
    all_rows <- bind_rows(all_rows, ds_label_row)
    
    # Get actual dataset names
    if (ds_type == "score_specific") {
      left_ds <- paste0(left_score, "_specific")
      right_ds <- paste0(right_score, "_specific")
    } else {
      left_ds <- "full_predictor_set"
      right_ds <- "full_predictor_set"
    }
    
    df_left <- format_block(all_metrics, left_score, left_ds)
    df_right <- format_block(all_metrics, right_score, right_ds)
    
    data_rows <- tibble(
      Model_L = df_left$Model,
      CVR2_L = df_left$`CV R²`,
      Q2_L = df_left$`Q²`,
      Gap_L = df_left$`R²-Q² Gap`,
      MAE_L = df_left$`Holdout MAE`,
      Spacer = rep("", nrow(df_left)),
      Model_R = df_right$Model,
      CVR2_R = df_right$`CV R²`,
      Q2_R = df_right$`Q²`,
      Gap_R = df_right$`R²-Q² Gap`,
      MAE_R = df_right$`Holdout MAE`
    )
    
    all_rows <- bind_rows(all_rows, data_rows)
  }
}

# Build flextable
ft <- flextable(all_rows) %>%
  set_header_labels(
    Model_L = "Model", CVR2_L = "CV R²", Q2_L = "Q²",
    Gap_L = "R²–Q² Gap", MAE_L = "Holdout MAE",
    Spacer = "",
    Model_R = "Model", CVR2_R = "CV R²", Q2_R = "Q²",
    Gap_R = "R²–Q² Gap", MAE_R = "Holdout MAE"
  ) %>%
  # Typography
  font(fontname = "Arial", part = "all") %>%
  fontsize(size = 9, part = "body") %>%
  fontsize(size = 9, part = "header") %>%
  bold(part = "header") %>%
  # Score label rows: bold, larger
  bold(i = label_row_indices, j = c(1, 7), part = "body") %>%
  fontsize(i = label_row_indices, size = 11, part = "body") %>%
  # Dataset label rows: italic
  italic(i = dataset_label_indices, j = c(1, 7), part = "body") %>%
  fontsize(i = dataset_label_indices, size = 9, part = "body") %>%
  # Alignment
  align(j = c(1, 7), align = "left", part = "all") %>%
  align(j = c(2:5, 8:11), align = "center", part = "all") %>%
  align(j = 6, align = "center", part = "all") %>%
  # Column widths
  width(j = 6, width = 0.15) %>%
  width(j = c(1, 7), width = 0.9) %>%
  width(j = c(2, 8), width = 1.1) %>%
  width(j = c(3, 9), width = 0.5) %>%
  width(j = c(4, 10), width = 0.7) %>%
  width(j = c(5, 11), width = 0.8) %>%
  # Borders
  border_remove() %>%
  hline_top(border = fp_border(width = 2), part = "all") %>%
  hline_bottom(border = fp_border(width = 2), part = "body") %>%
  hline(i = 1, border = fp_border(width = 1), part = "header")

# Add lines above each score label row
for (idx in label_row_indices) {
  ft <- ft %>%
    hline(i = idx, border = fp_border(width = 1.5), part = "body")
}

# Add subtle lines above each dataset label row
for (idx in dataset_label_indices) {
  ft <- ft %>%
    hline(i = idx - 1, border = fp_border(width = 0.5, color = "grey60"), part = "body")
}

# Save
doc <- read_docx() %>%
  body_end_section_landscape() %>%
  body_add_flextable(value = ft)
print(doc, target = file.path(output_dir, "Table_Model_Comparison.docx"))
cat("Saved Table_Model_Comparison.docx\n")

# ============================================
# 17. Visualisation — bar plots
# ============================================

model_order <- c("elastic_net", "sparse_pls", "random_forest", "xgboost", "svm", "knn")
model_labels <- c(elastic_net = "Elastic Net", sparse_pls = "Sparse PLS",
                  random_forest = "Random Forest", xgboost = "XGBoost",
                  svm = "SVM-RBF", knn = "k-NN")

outcome_labels <- c(QRISK3_risk = "QRISK3 Risk", SCORE2_score = "SCORE2 Risk",
                    frs_10y = "Framingham 10-year Risk", ascvd_10y = "ASCVD 10-year Risk")

model_colors <- c("Elastic Net" = "#E69F00", "Sparse PLS" = "#56B4E9",
                  "Random Forest" = "#009E73", "XGBoost" = "#F0E442",
                  "SVM-RBF" = "#0072B2", "k-NN" = "#D55E00")

score_specific_names <- c("QRISK3_risk_specific", "SCORE2_score_specific",
                          "frs_10y_specific", "ascvd_10y_specific")

# Score-specific bar plots
score_specific_data <- all_metrics %>%
  filter(dataset %in% score_specific_names,
         outcome %in% names(outcome_labels)) %>%
  mutate(
    outcome_label = factor(outcome_labels[outcome], levels = outcome_labels),
    model_label = factor(model_labels[model], levels = model_labels)
  )

p1 <- score_specific_data %>%
  ggplot(aes(x = cv_rsq, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(aes(xmin = cv_rsq - cv_rsq_sd, xmax = cv_rsq + cv_rsq_sd),
                 position = position_dodge(width = 0.8), height = 0.3) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Cross-Validation", x = expression(R^2), y = NULL, fill = "Model") +
  theme_minimal() +
  theme(legend.position = "none", panel.grid.major.y = element_blank())

p2 <- score_specific_data %>%
  filter(!is.na(q2)) %>%
  ggplot(aes(x = q2, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Holdout Set", x = expression(Q^2), y = NULL, fill = "Model") +
  theme_minimal() +
  theme(legend.position = "right", panel.grid.major.y = element_blank())

combined_ss <- (p1 | p2) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Model Performance: Score-Specific Inputs",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

ggsave(file.path(output_dir, "barplot_score_specific.png"), combined_ss, width = 12, height = 6, dpi = 300)
ggsave(file.path(output_dir, "barplot_score_specific.pdf"), combined_ss, width = 12, height = 6)

# Full predictor set bar plots
fps_data <- all_metrics %>%
  filter(dataset == "full_predictor_set",
         outcome %in% names(outcome_labels)) %>%
  mutate(
    outcome_label = factor(outcome_labels[outcome], levels = outcome_labels),
    model_label = factor(model_labels[model], levels = model_labels)
  )

p3 <- fps_data %>%
  ggplot(aes(x = cv_rsq, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(aes(xmin = cv_rsq - cv_rsq_sd, xmax = cv_rsq + cv_rsq_sd),
                 position = position_dodge(width = 0.8), height = 0.3) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Cross-Validation", x = expression(R^2), y = NULL, fill = "Model") +
  theme_minimal() +
  theme(legend.position = "none", panel.grid.major.y = element_blank())

p4 <- fps_data %>%
  filter(!is.na(q2)) %>%
  ggplot(aes(x = q2, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Holdout Set", x = expression(Q^2), y = NULL, fill = "Model") +
  theme_minimal() +
  theme(legend.position = "right", panel.grid.major.y = element_blank())

combined_fps <- (p3 | p4) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Model Performance: Full Predictor Set",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

ggsave(file.path(output_dir, "barplot_full_predictor_set.png"), combined_fps, width = 12, height = 6, dpi = 300)
ggsave(file.path(output_dir, "barplot_full_predictor_set.pdf"), combined_fps, width = 12, height = 6)

cat("Bar plots saved\n")

# ============================================
# 18. Tuning parameter diagnostics
# ============================================

cat("\n=== TUNING PARAMETER DIAGNOSTICS ===\n")

primary_tuning_params <- c(
  recipe_elastic_net = "penalty",
  recipe_sparse_pls = "num_comp",
  recipe_random_forest = "mtry",
  recipe_xgboost = "learn_rate",
  recipe_svm = "cost",
  recipe_knn = "neighbors"
)

model_labels <- c(elastic_net = "Elastic Net", sparse_pls = "Sparse PLS",
                  random_forest = "Random Forest", xgboost = "XGBoost",
                  svm = "SVM-RBF", knn = "k-NN")

outcome_labels <- c(QRISK3_risk = "QRISK3 Risk", SCORE2_score = "SCORE2 Risk",
                    frs_10y = "Framingham 10-year Risk", ascvd_10y = "ASCVD 10-year Risk")

dataset_display_labels <- c(
  "QRISK3_risk_specific" = "Score-specific",
  "SCORE2_score_specific" = "Score-specific",
  "frs_10y_specific" = "Score-specific",
  "ascvd_10y_specific" = "Score-specific",
  "full_predictor_set" = "Full predictor set"
)

param_display_names <- c(
  penalty = "Penalty (\u03bb)", num_comp = "Number of Components",
  mtry = "Variables per Split (mtry)", learn_rate = "Learning Rate",
  cost = "Cost (C)", neighbors = "Neighbors (k)"
)

# ---- Extract tuning data ----
extract_tuning_data_with_best <- function(result, dataset_name, outcome_name) {
  wf_res <- result$wf_results
  
  map_dfr(wf_res$wflow_id, function(wf_id) {
    primary_param <- primary_tuning_params[wf_id]
    if (is.na(primary_param)) return(NULL)
    
    tryCatch({
      res <- wf_res %>% extract_workflow_set_result(wf_id)
      metrics <- res %>% collect_metrics() %>% filter(.metric == "rmse")
      if (!primary_param %in% names(metrics)) return(NULL)
      
      model_name <- str_remove(wf_id, "recipe_")
      best_config <- tryCatch({
        switch(model_name,
               elastic_net   = res %>% select_by_one_std_err(metric = "rmse", desc(penalty)),
               sparse_pls    = res %>% select_by_one_std_err(metric = "rmse", num_comp),
               random_forest = res %>% select_by_one_std_err(metric = "rmse", min_n),
               xgboost       = res %>% select_by_one_std_err(metric = "rmse", min_n),
               svm           = res %>% select_by_one_std_err(metric = "rmse", cost),
               knn           = res %>% select_by_one_std_err(metric = "rmse", desc(neighbors)),
               res %>% select_best(metric = "rmse")
        )
      }, error = function(e) res %>% select_best(metric = "rmse"))
      
      metrics %>%
        mutate(
          model = model_name,
          model_label = model_labels[model_name],
          dataset = dataset_name,
          outcome = outcome_name,
          tuning_param_name = primary_param,
          tuning_param_value = .data[[primary_param]],
          is_selected = (.config == best_config$.config)
        ) %>%
        select(model, model_label, dataset, outcome,
               tuning_param_name, tuning_param_value,
               cv_rmse = mean, cv_rmse_se = std_err, is_selected)
    }, error = function(e) NULL)
  })
}

# Collect from all results
all_tuning_data <- map_dfr(names(all_results), function(task_name) {
  result <- all_results[[task_name]]
  parts <- str_split(task_name, " ~ ", simplify = TRUE)
  outcome_name <- parts[1]
  dataset_name <- parts[2]
  extract_tuning_data_with_best(result, dataset_name, outcome_name)
})

if (nrow(all_tuning_data) > 0) {
  
  all_tuning_data <- all_tuning_data %>%
    filter(outcome %in% names(outcome_labels)) %>%
    mutate(
      outcome_label = factor(outcome_labels[outcome], levels = outcome_labels),
      model_label = factor(model_label, levels = model_labels),
      dataset_label = factor(
        ifelse(str_detect(dataset, "_specific"), "Score-specific inputs", "Full predictor set"),
        levels = c("Score-specific inputs", "Full predictor set")
      )
    )
  
  # ---- 18a. Per-model tuning plot function ----
  create_model_tuning_plot <- function(data, model_name, plot_subtitle) {
    model_data <- data %>% filter(model == model_name)
    if (nrow(model_data) == 0) return(NULL)
    
    param_name <- unique(model_data$tuning_param_name)
    x_label <- param_display_names[param_name]
    selected_pts <- model_data %>% filter(is_selected)
    
    p <- ggplot(model_data, aes(x = tuning_param_value, y = cv_rmse)) +
      geom_line(colour = "grey60", linewidth = 0.4) +
      geom_point(colour = "grey40", size = 1.5, alpha = 0.6) +
      geom_errorbar(aes(ymin = cv_rmse - cv_rmse_se, ymax = cv_rmse + cv_rmse_se),
                    colour = "grey70", alpha = 0.4, width = 0) +
      geom_point(data = selected_pts, colour = "#E41A1C", size = 3, shape = 18) +
      geom_vline(data = selected_pts, aes(xintercept = tuning_param_value),
                 colour = "#E41A1C", linetype = "dashed", alpha = 0.5) +
      facet_wrap(~ outcome_label, scales = "free_y", ncol = 1) +
      labs(title = model_labels[model_name], subtitle = plot_subtitle,
           x = x_label, y = "CV RMSE") +
      theme_minimal(base_size = 10) +
      theme(strip.text = element_text(face = "bold", size = 9),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 9, colour = "grey40"))
    
    if (param_name %in% c("penalty", "cost")) p <- p + scale_x_log10()
    p
  }
  
  # ---- 18b. Per-dataset overview plots (existing style) ----
  model_order <- c("elastic_net", "sparse_pls", "random_forest", "xgboost", "svm", "knn")
  
  for (ds_pattern in c("_specific", "full_predictor_set")) {
    if (ds_pattern == "_specific") {
      ds_data <- all_tuning_data %>% filter(str_detect(dataset, "_specific"))
      ds_suffix <- "score_specific"
      ds_title <- "Score-Specific Inputs"
    } else {
      ds_data <- all_tuning_data %>% filter(dataset == "full_predictor_set")
      ds_suffix <- "full_predictor_set"
      ds_title <- "Full Predictor Set"
    }
    
    if (nrow(ds_data) > 0) {
      plots_list <- list()
      for (m in model_order) {
        p <- create_model_tuning_plot(ds_data, m, ds_title)
        if (!is.null(p)) plots_list[[m]] <- p
      }
      
      if (length(plots_list) > 0) {
        combined <- wrap_plots(plots_list, ncol = 6) +
          plot_annotation(
            title = paste("Tuning Parameter Diagnostics:", ds_title),
            subtitle = "Red diamond = selected configuration (one-SE rule on RMSE)",
            theme = theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = 10)
            )
          )
        
        ggsave(file.path(output_dir, paste0("tuning_diagnostics_", ds_suffix, ".png")),
               combined, width = 20, height = 12, dpi = 300)
        ggsave(file.path(output_dir, paste0("tuning_diagnostics_", ds_suffix, ".pdf")),
               combined, width = 20, height = 12)
        cat(sprintf("  Saved tuning_diagnostics_%s\n", ds_suffix))
      }
    }
  }
  
  # ---- 18c. Per-score detailed tuning plots ----
  # One figure per CVD score: Row 1 = score-specific, Row 2 = full predictor set
  # 6 columns = 6 models
  
  create_single_score_tuning_plot <- function(data, model_name, dataset_lab) {
    model_data <- data %>% filter(model == model_name)
    if (nrow(model_data) == 0) return(NULL)
    
    param_name <- unique(model_data$tuning_param_name)
    x_label <- param_display_names[param_name]
    selected_pts <- model_data %>% filter(is_selected)
    
    p <- ggplot(model_data, aes(x = tuning_param_value, y = cv_rmse)) +
      geom_line(colour = "grey60", linewidth = 0.4) +
      geom_point(colour = "grey40", size = 1.5, alpha = 0.6) +
      geom_errorbar(aes(ymin = cv_rmse - cv_rmse_se, ymax = cv_rmse + cv_rmse_se),
                    colour = "grey70", alpha = 0.4, width = 0) +
      geom_point(data = selected_pts, colour = "#E41A1C", size = 3, shape = 18) +
      geom_vline(data = selected_pts, aes(xintercept = tuning_param_value),
                 colour = "#E41A1C", linetype = "dashed", alpha = 0.5) +
      labs(title = model_labels[model_name], x = x_label, y = "CV RMSE") +
      theme_minimal(base_size = 9) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 10, hjust = 0.5))
    
    if (param_name %in% c("penalty", "cost")) p <- p + scale_x_log10()
    p
  }
  
  for (score in names(outcome_labels)) {
    score_label <- outcome_labels[score]
    score_specific_ds <- paste0(score, "_specific")
    
    # Get data for this score
    ss_data <- all_tuning_data %>%
      filter(outcome == score, dataset == score_specific_ds)
    fps_data <- all_tuning_data %>%
      filter(outcome == score, dataset == "full_predictor_set")
    
    if (nrow(ss_data) == 0 && nrow(fps_data) == 0) next
    
    # Row 1: score-specific (6 models)
    row1_plots <- list()
    for (m in model_order) {
      p <- create_single_score_tuning_plot(ss_data, m, "Score-specific")
      if (!is.null(p)) {
        row1_plots[[m]] <- p
      } else {
        row1_plots[[m]] <- ggplot() + theme_void()
      }
    }
    
    # Row 2: full predictor set (6 models)
    row2_plots <- list()
    for (m in model_order) {
      p <- create_single_score_tuning_plot(fps_data, m, "Full predictor set")
      if (!is.null(p)) {
        row2_plots[[m]] <- p
      } else {
        row2_plots[[m]] <- ggplot() + theme_void()
      }
    }
    
    # Combine: row 1 on top, row 2 below
    top_row <- wrap_plots(row1_plots, ncol = 6)
    bottom_row <- wrap_plots(row2_plots, ncol = 6)
    
    combined_score <- top_row / bottom_row +
      plot_annotation(
        title = paste("Tuning Parameter Diagnostics:", score_label),
        subtitle = "Top row: Score-specific inputs | Bottom row: Full predictor set\nRed diamond = selected configuration (one-SE rule on RMSE)",
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = 10)
        )
      )
    
    ggsave(file.path(output_dir, paste0("tuning_detailed_", score, ".png")),
           combined_score, width = 20, height = 8, dpi = 300)
    ggsave(file.path(output_dir, paste0("tuning_detailed_", score, ".pdf")),
           combined_score, width = 20, height = 8)
    cat(sprintf("  Saved tuning_detailed_%s\n", score))
  }
  
} else {
  cat("  No tuning data available\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")