################################################################################
# CoDiet CVD Prediction Model Comparison
# Author: Luisa Delius
#
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

CPUS <- parallel::detectCores() - 1
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
         -stress_resilience_status, -stress_index_status, -Age.Risk) %>%
  full_join(df_fatty_acids_statin_suppl %>% select(Sample_ID, Statins, Supplements), 
            by = "Sample_ID")

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

participant_info <- df_all_data_model %>%
  select(Sample_ID) %>%
  mutate(
    missing_pct = rowMeans(is.na(df_all_data_model %>% select(-Sample_ID, -all_of(CVD_scores)))) * 100
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

sanitize_df <- function(df) {
  for (col in names(df)) {
    if (col %in% c("Sample_ID", CVD_scores)) next
    x <- df[[col]]
    
    if (is.logical(x)) df[[col]] <- as.integer(x)
    else if (inherits(x, c("Date", "POSIXct", "POSIXlt"))) df[[col]] <- as.numeric(x)
    else if (!is.numeric(x) && !is.factor(x)) df[[col]] <- as.factor(x)
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

# Overfitting threshold: absolute difference between train and holdout R²
# E.g., 0.25 means train R²=0.80 and holdout R²=0.55 is acceptable (gap=0.25)
RSQ_GAP_THRESHOLD <- 0.25

cat("\n=== SETTINGS ===\n")
cat(sprintf("CV: %d-fold x %d repeats\n", CV_FOLDS, CV_REPEATS))
cat(sprintf("Grid size: %d\n", GRID_SIZE))
cat(sprintf("R² gap threshold: %.2f (absolute difference)\n", RSQ_GAP_THRESHOLD))

# ============================================
# 9. Model specifications
# ============================================

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

# ============================================
# 10. Best model selection function
# ============================================
# Selects best model based on holdout MAE, excluding models with excessive overfitting.
# rsq_gap = train_rsq - holdout_rsq (absolute difference, not ratio)

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
  
  wf_results <- wf_set %>%
    workflow_map(
      fn = "tune_grid",
      resamples = cv_splits,
      grid = grid_size,
      metrics = metric_set(mae, rmse, rsq),
      control = control_grid(save_pred = TRUE, save_workflow = FALSE, 
                             verbose = FALSE, parallel_over = "everything"),
      seed = seed,
      verbose = TRUE
    )
  
  # Extract metrics for all models
  all_model_metrics <- tibble()
  
  for (model_name in wf_results$wflow_id) {
    best_result <- wf_results %>% extract_workflow_set_result(model_name)
    best_params <- best_result %>% select_best(metric = "rmse")
    
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
    
    # Holdout metrics
    holdout_mae <- holdout_rmse <- holdout_rsq <- holdout_cal_slope <- holdout_cal_intercept <- NA
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
    
    # Overfitting indicators
    rsq_gap <- train_rsq - holdout_rsq
    rsq_ratio <- ifelse(!is.na(holdout_rsq) & holdout_rsq > 0.01, train_rsq / holdout_rsq, NA)
    
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
      holdout_rsq = round(holdout_rsq, 3),
      holdout_cal_slope = round(holdout_cal_slope, 3),
      holdout_cal_intercept = round(holdout_cal_intercept, 3),
      rsq_gap = round(rsq_gap, 3),
      rsq_ratio = round(rsq_ratio, 3),
      overfit_flag = case_when(
        is.na(holdout_rsq) ~ NA,
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
    best_params <- wf_results %>%
      extract_workflow_set_result(paste0("recipe_", best_model_name)) %>%
      select_best(metric = "rmse")
    
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
  best_config <- best_result %>% select_best(metric = "rmse")
  
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
# 12. Build task list
# ============================================

cat("\n=== BUILDING TASK LIST ===\n")

task_list <- list()
task_counter <- 1

# A) Predictor blocks x CVD scores
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

# ============================================
# 13. Test run
# ============================================

cat("\n=== TEST RUN ===\n")
plan(multisession, workers = CPUS)

test_result <- run_models_for_dataset(
  df = datasets_list$lipids,
  outcome = "QRISK3_risk",
  dataset_name = "lipids",
  model_specs = model_specs,
  cv_folds = CV_FOLDS,
  cv_repeats = CV_REPEATS,
  grid_size = GRID_SIZE,
  seed = 42,
  calculate_importance = FALSE,
  holdout_ids = HOLDOUT_IDS
)

plan(sequential)

if (!is.null(test_result)) {
  cat("Test successful. Best model:", test_result$best_model, "\n")
  print(test_result$metrics %>% select(model, train_rsq, cv_rsq, holdout_rsq, rsq_gap))
}

# ============================================
# 14. Full run
# ============================================

# plan(multisession, workers = CPUS)
# 
# all_results <- list()
# failed_tasks <- c()
# start_time <- Sys.time()
# 
# for (i in seq_along(task_list)) {
#   task <- task_list[[i]]
#   task_name <- paste0(task$outcome, " ~ ", task$dataset_name)
#   
#   cat(sprintf("\n[%d/%d] %s\n", i, length(task_list), task_name))
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
#       holdout_ids = task$holdout_ids
#     )
#   }, error = function(e) {
#     cat(sprintf("  ERROR: %s\n", e$message))
#     return(NULL)
#   })
#   
#   if (!is.null(result)) {
#     all_results[[task_name]] <- result
#     cat("  Done\n")
#   } else {
#     failed_tasks <- c(failed_tasks, task_name)
#     cat("  Failed\n")
#   }
#   
#   # Intermediate save every 10 tasks
#   if (i %% 10 == 0) {
#     qs2::qs_save(all_results, file.path(output_dir, "ml_comparison_all_results_partial.qs2"))
#     cat(sprintf("  [Saved checkpoint: %d tasks]\n", i))
#   }
# }
# 
# plan(sequential)
# 
# total_time <- round(difftime(Sys.time(), start_time, units = "hours"), 2)
# cat(sprintf("\nCompleted in %.2f hours\n", total_time))
# cat(sprintf("Successful: %d | Failed: %d\n", length(all_results), length(failed_tasks)))
# 
# qs2::qs_save(all_results, file.path(output_dir, "ml_comparison_all_results.qs2"))
# if (length(failed_tasks) > 0) writeLines(failed_tasks, file.path(output_dir, "ml_comparison_failed_tasks.txt"))

# ============================================
# 15. Compile and summarise results
# ============================================

all_metrics <- bind_rows(lapply(all_results, function(x) x$metrics))
all_importance <- bind_rows(lapply(all_results, function(x) x$importance))
all_cv_preds <- bind_rows(lapply(all_results, function(x) x$cv_predictions))
all_holdout_preds <- bind_rows(lapply(all_results, function(x) x$holdout_predictions))

# Best models summary
best_models_summary <- map_dfr(names(all_results), function(task_name) {
  result <- all_results[[task_name]]
  result$metrics %>% filter(model == result$best_model)
})

# Model rankings
model_rankings <- all_metrics %>%
  filter(!is.na(holdout_mae)) %>%
  group_by(dataset, outcome) %>%
  mutate(rank = rank(holdout_mae)) %>%
  ungroup() %>%
  group_by(model) %>%
  summarise(
    mean_rank = round(mean(rank), 2),
    times_best = sum(rank == 1),
    mean_holdout_rsq = round(mean(holdout_rsq, na.rm = TRUE), 3),
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

cat("\nResults saved to", file.path(output_dir, "ml_model_comparison_results.xlsx"), "\n")

# ============================================
# 16. Create tables
# ============================================

# Definitions
dataset_order <- c(
  "QRISK3_risk_specific", "SCORE2_score_specific", "frs_10y_specific", 
  "ascvd_10y_specific", "mean_risk_specific", "all_data", "fatty_acids",
  "lipids", "urine_NMR", "body_composition", "sociodemographics_lifestyle",
  "clinical_risk_factors"
)

dataset_labels <- c(
  QRISK3_risk_specific = "Score-specific inputs", SCORE2_score_specific = "Score-specific inputs",
  frs_10y_specific = "Score-specific inputs", ascvd_10y_specific = "Score-specific inputs",
  mean_risk_specific = "Score-specific inputs", all_data = "All data",
  fatty_acids = "Fatty acids", lipids = "Lipidomics", urine_NMR = "Urinary NMR",
  body_composition = "Body composition", sociodemographics_lifestyle = "Sociodemographics/Lifestyle",
  clinical_risk_factors = "Clinical risk factors"
)

model_order <- c("elastic_net", "pls", "random_forest", "xgboost", "svm", "knn")
model_labels <- c(elastic_net = "Elastic Net", pls = "PLS", random_forest = "Random Forest",
                  xgboost = "XGBoost", svm = "SVM-RBF", knn = "k-NN")

outcome_labels <- c(QRISK3_risk = "QRISK3 Risk", SCORE2_score = "SCORE2 Risk",
                    frs_10y = "Framingham 10-year Risk", ascvd_10y = "ASCVD 10-year Risk",
                    mean_risk = "Composite Score Risk")

create_publication_table <- function(outcome_name, metrics_df) {
  table_data <- metrics_df %>%
    filter(outcome == outcome_name, dataset %in% dataset_order) %>%
    mutate(
      model = factor(model, levels = model_order),
      dataset = factor(dataset, levels = dataset_order)
    ) %>%
    arrange(dataset, model)
  
  formatted <- table_data %>%
    mutate(
      Dataset = dataset_labels[as.character(dataset)],
      Model = model_labels[as.character(model)],
      `Train R2` = sprintf("%.3f", train_rsq),
      `CV R2` = sprintf("%.3f +/- %.3f", cv_rsq, cv_rsq_sd),
      `Holdout R2` = ifelse(!is.na(holdout_rsq), sprintf("%.3f", holdout_rsq), "-"),
      `R2 Gap` = ifelse(!is.na(rsq_gap), sprintf("%.3f", rsq_gap), "-"),
      `CV MAE` = sprintf("%.3f +/- %.3f", cv_mae, cv_mae_sd),
      `Holdout MAE` = ifelse(!is.na(holdout_mae), sprintf("%.3f", holdout_mae), "-"),
      `Cal. Slope` = ifelse(!is.na(holdout_cal_slope), sprintf("%.3f", holdout_cal_slope), "-")
    ) %>%
    select(Dataset, Model, `Train R2`, `CV R2`, `Holdout R2`, `R2 Gap`, 
           `CV MAE`, `Holdout MAE`, `Cal. Slope`)
  
  flextable(formatted) %>%
    merge_v(j = "Dataset") %>%
    set_caption(sprintf("Model Performance for %s", outcome_labels[outcome_name])) %>%
    align(j = 1:2, align = "left", part = "all") %>%
    align(j = 3:9, align = "center", part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    fontsize(size = 8, part = "body") %>%
    fontsize(size = 9, part = "header") %>%
    bold(part = "header") %>%
    autofit() %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1.5), part = "header") %>%
    hline(i = seq(6, nrow(formatted), 6), border = fp_border(width = 1), part = "body")
}

# Generate tables for each outcome
cat("\n=== GENERATING TABLES ===\n")

for (outcome in names(outcome_labels)) {
  ft <- create_publication_table(outcome, all_metrics)
  
  doc <- read_docx() %>%
    body_end_section_landscape() %>%
    body_add_flextable(value = ft)
  print(doc, target = file.path(output_dir, sprintf("Table_%s.docx", outcome)))
  
  cat(sprintf("  Saved Table_%s.docx\n", outcome))
}

# ============================================
# 17. Summary tables
# ============================================

create_summary_table <- function(data, title) {
  dataset_levels <- unique(data$dataset_label)
  
  formatted <- data %>%
    mutate(
      model = factor(model, levels = model_order),
      dataset_label = factor(dataset_label, levels = dataset_levels)
    ) %>%
    arrange(dataset_label, model) %>%
    mutate(
      Dataset = as.character(dataset_label),
      Model = model_labels[as.character(model)],
      `Train R2` = sprintf("%.3f", train_rsq),
      `CV R2` = sprintf("%.3f +/- %.3f", cv_rsq, cv_rsq_sd),
      `Holdout R2` = ifelse(!is.na(holdout_rsq), sprintf("%.3f", holdout_rsq), "-"),
      `R2 Gap` = ifelse(!is.na(rsq_gap), sprintf("%.3f", rsq_gap), "-"),
      `CV MAE` = sprintf("%.3f +/- %.3f", cv_mae, cv_mae_sd),
      `Holdout MAE` = ifelse(!is.na(holdout_mae), sprintf("%.3f", holdout_mae), "-")
    ) %>%
    select(Dataset, Model, `Train R2`, `CV R2`, `Holdout R2`, `R2 Gap`, `CV MAE`, `Holdout MAE`)
  
  flextable(formatted) %>%
    merge_v(j = "Dataset") %>%
    set_caption(title) %>%
    align(j = 1:2, align = "left", part = "all") %>%
    align(j = 3:8, align = "center", part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    fontsize(size = 8, part = "body") %>%
    fontsize(size = 9, part = "header") %>%
    bold(part = "header") %>%
    autofit() %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1.5), part = "header")
}

# Table A: Score-specific (averaged across all 5 CVD scores)
score_specific_datasets_names <- c("QRISK3_risk_specific", "SCORE2_score_specific", 
                                   "frs_10y_specific", "ascvd_10y_specific", "mean_risk_specific")

table_a_data <- all_metrics %>%
  filter(dataset %in% score_specific_datasets_names) %>%
  group_by(model) %>%
  summarise(across(c(train_rsq, cv_rsq, cv_rsq_sd, holdout_rsq, rsq_gap, cv_mae, cv_mae_sd, holdout_mae), 
                   ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(dataset_label = "Score-specific inputs")

ft_a <- create_summary_table(table_a_data, "Table A: Score-Specific Inputs (averaged across CVD scores)")

# Table B: General datasets
general_datasets <- c("all_data", "fatty_acids", "lipids", "urine_NMR", 
                      "body_composition", "sociodemographics_lifestyle", "clinical_risk_factors")

table_b_data <- all_metrics %>%
  filter(dataset %in% general_datasets) %>%
  group_by(dataset, model) %>%
  summarise(across(c(train_rsq, cv_rsq, cv_rsq_sd, holdout_rsq, rsq_gap, cv_mae, cv_mae_sd, holdout_mae), 
                   ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(dataset_label = dataset_labels[dataset])

ft_b <- create_summary_table(table_b_data, "Table B: General Datasets (averaged across CVD scores)")

# Table C: Score-specific + block
block_labels <- c(
  all_data = "All data", fatty_acids = "Fatty acids", lipids = "Lipidomics",
  urine_NMR = "Urinary NMR", body_composition = "Body composition",
  sociodemographics_lifestyle = "Sociodemographics/Lifestyle",
  clinical_risk_factors = "Clinical risk factors"
)

table_c_data <- all_metrics %>%
  filter(str_detect(dataset, "_plus_")) %>%
  mutate(block = str_extract(dataset, "(?<=_plus_).*")) %>%
  filter(!is.na(block), block %in% names(block_labels)) %>%
  group_by(block, model) %>%
  summarise(across(c(train_rsq, cv_rsq, cv_rsq_sd, holdout_rsq, rsq_gap, cv_mae, cv_mae_sd, holdout_mae), 
                   ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(dataset_label = paste("Score-specific +", block_labels[block]))

ft_c <- create_summary_table(table_c_data, "Table C: Score-Specific + Block (averaged across CVD scores)")

# Save summary tables
for (tbl in list(list(ft_a, "Summary_Table_A"), list(ft_b, "Summary_Table_B"), list(ft_c, "Summary_Table_C"))) {
  doc <- read_docx() %>% body_end_section_landscape() %>% body_add_flextable(value = tbl[[1]])
  print(doc, target = file.path(output_dir, paste0(tbl[[2]], ".docx")))
}

cat("Summary tables saved\n")

# ============================================
# 18. Visualisation
# ============================================

score_specific_data <- all_metrics %>%
  filter(dataset %in% score_specific_datasets_names) %>%
  mutate(
    outcome = factor(outcome, levels = names(outcome_labels)),
    outcome_label = factor(outcome_labels[outcome], levels = outcome_labels),
    model = factor(model, levels = model_order),
    model_label = factor(model_labels[model], levels = model_labels)
  )

model_colors <- c("Elastic Net" = "#E69F00", "PLS" = "#56B4E9", "Random Forest" = "#009E73",
                  "XGBoost" = "#F0E442", "SVM-RBF" = "#0072B2", "k-NN" = "#D55E00")

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
  filter(!is.na(holdout_rsq)) %>%
  ggplot(aes(x = holdout_rsq, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Holdout Set", x = expression(R^2), y = NULL, fill = "Model") +
  theme_minimal() +
  theme(legend.position = "right", panel.grid.major.y = element_blank())

p3 <- score_specific_data %>%
  ggplot(aes(x = cv_mae, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(aes(xmin = cv_mae - cv_mae_sd, xmax = cv_mae + cv_mae_sd),
                 position = position_dodge(width = 0.8), height = 0.3) +
  scale_fill_manual(values = model_colors) +
  labs(x = "MAE", y = NULL, fill = "Model") +
  theme_minimal() +
  theme(legend.position = "none", panel.grid.major.y = element_blank())

p4 <- score_specific_data %>%
  filter(!is.na(holdout_mae)) %>%
  ggplot(aes(x = holdout_mae, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = model_colors) +
  labs(x = "MAE", y = NULL, fill = "Model") +
  theme_minimal() +
  theme(legend.position = "right", panel.grid.major.y = element_blank())

combined_plot <- (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Model Performance: Score-Specific Inputs",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

ggsave(file.path(output_dir, "model_comparison_score_specific.png"), combined_plot, width = 12, height = 8, dpi = 300)
ggsave(file.path(output_dir, "model_comparison_score_specific.pdf"), combined_plot, width = 12, height = 8)

cat("\nPlots saved\n")
cat("\n=== ANALYSIS COMPLETE ===\n")