################################################################################
# CoDiet CVD Prediction Model Comparison — Nested Cross-Validation
# Author: Luisa Delius
#
# Architecture:
#   Outer loop: 5-fold CV x 5 repeats (unbiased performance estimation)
#   Inner loop: 5-fold CV within each outer training fold (hyperparameter tuning)
#   Hyperparameter selection: one-SE rule on RMSE
#   Best model selection: lowest mean outer-fold CV MAE
################################################################################

# ============================================
# 1. Setup
# ============================================

# ---- 1.1 Load packages ----
library(tidyverse)
library(tidymodels)
library(future)
library(qs2)
library(plsmod)
library(openxlsx)
library(flextable)
library(officer)
library(patchwork)

# ---- 1.2 Global settings ----
set.seed(42)
options(scipen = 999, expressions = 500000)

CPUS <- parallel::detectCores() - 3
cat("Using", CPUS, "CPU cores for parallel processing\n")

select <- dplyr::select
slice  <- dplyr::slice
rename <- dplyr::rename
filter <- dplyr::filter

wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files/processed_data"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

# Create output folder for all results
base_name <- "ml_nested_cv_results_v_"
existing_dirs <- list.dirs(wkdir, full.names = FALSE, recursive = FALSE)
existing_versions <- existing_dirs[grepl(paste0("^", base_name), existing_dirs)]

if (length(existing_versions) == 0) {
  version_num <- 1
} else {
  version_nums <- as.numeric(sub(paste0(base_name, "(\\d+)"), "\\1", existing_versions))
  version_num <- max(version_nums, na.rm = TRUE) + 1
}

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
  all_data               = df_all_data_model,
  fatty_acids            = df_fatty_acids_model,
  lipids                 = df_lipidomics_model,
  urine_NMR              = df_urine_nmr_model,
  body_composition       = df_body_composition_model,
  sociodemographics_lifestyle = df_REDcap_demographics_model,
  clinical_risk_factors  = df_risk_factors_model
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
  QRISK3_risk  = df_model0_QRISK3,
  SCORE2_score = df_model0_score2,
  frs_10y      = df_model0_frs,
  ascvd_10y    = df_model0_ascvd,
  mean_risk    = df_model0_composite
)

# ============================================
# 3. Helper functions
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

# ============================================
# 4. Apply column missingness exclusion
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

datasets_list <- lapply(datasets_list, remove_high_missing_cols)
score_specific_datasets <- lapply(score_specific_datasets, remove_high_missing_cols)

# ============================================
# 5. Data sanitisation
# ============================================

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
# 6. Score-specific + Block combinations (UNUSED — kept for reference)
# ============================================

# cat("\n=== CREATING SCORE-SPECIFIC + BLOCK COMBINATIONS ===\n")
# 
# score_plus_block_datasets <- list()
# for (score in names(score_specific_datasets)) {
#   for (block in names(datasets_list)) {
#     combo_name <- paste0(score, "_plus_", block)
#     score_plus_block_datasets[[combo_name]] <- build_score_plus_block(
#       score, block, score_specific_datasets, datasets_list
#     )
#   }
# }
# cat(sprintf("Created %d combinations\n", length(score_plus_block_datasets)))


# ============================================
# 7. CV and tuning settings
# ============================================

OUTER_FOLDS   <- 5
OUTER_REPEATS <- 5
INNER_FOLDS   <- 5
GRID_SIZE     <- 20

cat("\n=== SETTINGS ===\n")
cat(sprintf("Outer CV: %d-fold x %d repeats\n", OUTER_FOLDS, OUTER_REPEATS))
cat(sprintf("Inner CV: %d-fold (for hyperparameter tuning)\n", INNER_FOLDS))
cat(sprintf("Grid size: %d\n", GRID_SIZE))

# ============================================
# 8. Model specifications
# ============================================

elastic_net_spec <- linear_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>% set_mode("regression")

sparse_pls_spec <- pls(num_comp = tune(), predictor_prop = tune()) %>%
  set_engine("mixOmics") %>% set_mode("regression")

rf_spec <- rand_forest(trees = tune(), mtry = tune(), min_n = tune()) %>%
  set_engine("ranger", importance = "permutation") %>% set_mode("regression")

xgb_spec <- boost_tree(
  trees = tune(), tree_depth = tune(), min_n = tune(),
  loss_reduction = tune(), sample_size = tune(), mtry = tune(), learn_rate = tune()
) %>% set_engine("xgboost") %>% set_mode("regression")

svm_spec <- svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
  set_engine("kernlab") %>% set_mode("regression")

knn_spec <- nearest_neighbor(neighbors = tune(), weight_func = tune(), dist_power = tune()) %>%
  set_engine("kknn") %>% set_mode("regression")

model_specs <- list(
  elastic_net = elastic_net_spec, sparse_pls = sparse_pls_spec,
  random_forest = rf_spec, xgboost = xgb_spec,
  svm = svm_spec, knn = knn_spec
)

# ============================================
# 9. One-SE rule helper
# ============================================

select_one_se <- function(tune_result, wf_id) {
  model_type <- str_remove(wf_id, "recipe_")
  tryCatch({
    switch(model_type,
           elastic_net   = tune_result %>% select_by_one_std_err(metric = "rmse", desc(penalty)),
           sparse_pls    = tune_result %>% select_by_one_std_err(metric = "rmse", num_comp),
           random_forest = tune_result %>% select_by_one_std_err(metric = "rmse", min_n),
           xgboost       = tune_result %>% select_by_one_std_err(metric = "rmse", min_n),
           svm           = tune_result %>% select_by_one_std_err(metric = "rmse", cost),
           knn           = tune_result %>% select_by_one_std_err(metric = "rmse", desc(neighbors)),
           tune_result %>% select_best(metric = "rmse")
    )
  }, error = function(e) tune_result %>% select_best(metric = "rmse"))
}

# ============================================
# 10. Main nested CV function
# ============================================

run_nested_cv <- function(df, outcome, dataset_name,
                          model_specs,
                          outer_folds = 5, outer_repeats = 5,
                          inner_folds = 5, grid_size = 20,
                          seed = 42) {
  
  set.seed(seed)
  cvd_scores_local <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y", "mean_risk")
  
  if (!outcome %in% names(df)) {
    warning(sprintf("Outcome '%s' not found in '%s'", outcome, dataset_name))
    return(NULL)
  }
  
  # ---- Prepare data ----
  other_scores <- setdiff(cvd_scores_local, outcome)
  df_clean <- df %>%
    select(-any_of(other_scores)) %>%
    filter(!is.na(.data[[outcome]]))
  
  predictor_cols <- setdiff(names(df_clean), c("Sample_ID", outcome))
  if (length(predictor_cols) > 0) {
    row_missing <- rowMeans(is.na(df_clean[predictor_cols])) * 100
    high_missing <- row_missing > 80
    if (any(high_missing)) {
      cat(sprintf("    Removing %d participants with >80%% missing\n", sum(high_missing)))
      df_clean <- df_clean[!high_missing, ]
    }
  }
  
  n_total <- nrow(df_clean)
  n_predictors <- length(predictor_cols)
  
  if (n_total < 30) {
    warning(sprintf("Too few observations (%d) for %s ~ %s", n_total, outcome, dataset_name))
    return(NULL)
  }
  
  cat(sprintf("  %s ~ %s: n=%d, p=%d\n", outcome, dataset_name, n_total, n_predictors))
  
  # ---- Outer CV splits ----
  outer_splits <- vfold_cv(df_clean, v = outer_folds, repeats = outer_repeats,
                           strata = !!sym(outcome))
  
  # ---- Recipe builder ----
  build_recipe <- function(train_data) {
    recipe(as.formula(paste(outcome, "~ .")), data = train_data) %>%
      update_role(Sample_ID, new_role = "ID") %>%
      step_zv(all_predictors()) %>%
      step_impute_mode(all_nominal_predictors()) %>%
      step_impute_knn(all_numeric_predictors(), neighbors = 5) %>%
      step_impute_median(all_numeric_predictors()) %>%
      step_zv(all_predictors()) %>%
      step_YeoJohnson(all_numeric_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_dummy(all_nominal_predictors())
  }
  
  # ---- Inner tuning function ----
  run_inner_tuning <- function(train_data) {
    rec <- build_recipe(train_data)
    inner_cv <- vfold_cv(train_data, v = inner_folds, strata = !!sym(outcome))
    
    wf_set <- workflow_set(preproc = list(recipe = rec), models = model_specs, cross = TRUE)
    wf_set <- wf_set %>% option_add(grid = grid_size)
    
    n_pred_local <- ncol(train_data) - 2
    max_comp <- min(8L, max(4L, n_pred_local - 4L))
    spls_params <- wf_set %>%
      extract_parameter_set_dials("recipe_sparse_pls") %>%
      update(num_comp = num_comp(range = c(1L, max_comp)))
    wf_set <- wf_set %>%
      option_remove(id = "recipe_sparse_pls") %>%
      option_add(param_info = spls_params, id = "recipe_sparse_pls")
    
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
        resamples = inner_cv,
        metrics = metric_set(mae, rmse, rsq),
        control = control_grid(save_pred = FALSE, save_workflow = FALSE,
                               verbose = FALSE, parallel_over = "everything"),
        seed = seed, verbose = FALSE
      )
    
    best_params_list <- list()
    for (wf_id in wf_results$wflow_id) {
      model_name <- str_remove(wf_id, "recipe_")
      tune_res <- wf_results %>% extract_workflow_set_result(wf_id)
      best_params_list[[model_name]] <- select_one_se(tune_res, wf_id)
    }
    
    list(wf_results = wf_results, best_params = best_params_list)
  }
  
  # ---- Outer loop ----
  outer_fold_results <- list()
  
  for (i in seq_len(nrow(outer_splits))) {
    split_obj  <- outer_splits$splits[[i]]
    fold_id    <- outer_splits$id[[i]]
    repeat_id  <- if ("id2" %in% names(outer_splits)) outer_splits$id2[[i]] else "Repeat1"
    fold_label <- paste(repeat_id, fold_id, sep = "_")
    
    cat(sprintf("    Outer fold %d/%d (%s) ... ", i, nrow(outer_splits), fold_label))
    
    outer_train <- analysis(split_obj)
    outer_test  <- assessment(split_obj)
    
    inner_result <- tryCatch(
      run_inner_tuning(outer_train),
      error = function(e) { cat(sprintf("INNER ERROR: %s\n", e$message)); NULL }
    )
    
    if (is.null(inner_result)) { cat("skipped\n"); next }
    
    fold_metrics <- list()
    for (model_name in names(model_specs)) {
      wf_id <- paste0("recipe_", model_name)
      best_params <- inner_result$best_params[[model_name]]
      
      fit_result <- tryCatch({
        final_wf <- inner_result$wf_results %>%
          extract_workflow(wf_id) %>%
          finalize_workflow(best_params)
        final_fit <- fit(final_wf, data = outer_train)
        
        test_preds <- predict(final_fit, new_data = outer_test) %>%
          bind_cols(outer_test %>% select(all_of(outcome)))
        train_preds <- predict(final_fit, new_data = outer_train) %>%
          bind_cols(outer_train %>% select(all_of(outcome)))
        
        # ---- FIX: Classic R² and Q² (1 - SS_res/SS_tot) ----
        train_mean <- mean(train_preds[[outcome]], na.rm = TRUE)
        
        # Training R²
        ss_res_train <- sum((train_preds[[outcome]] - train_preds$.pred)^2, na.rm = TRUE)
        ss_tot_train <- sum((train_preds[[outcome]] - train_mean)^2, na.rm = TRUE)
        train_rsq <- 1 - ss_res_train / ss_tot_train
        
        # Q² (held-out test data, using TRAINING mean)
        ss_res_test <- sum((test_preds[[outcome]] - test_preds$.pred)^2, na.rm = TRUE)
        ss_tot_test <- sum((test_preds[[outcome]] - train_mean)^2, na.rm = TRUE)
        test_rsq <- 1 - ss_res_test / ss_tot_test
        
        test_mae  <- mean(abs(test_preds$.pred - test_preds[[outcome]]), na.rm = TRUE)
        test_rmse <- sqrt(mean((test_preds$.pred - test_preds[[outcome]])^2, na.rm = TRUE))
        
        train_mae  <- mean(abs(train_preds$.pred - train_preds[[outcome]]), na.rm = TRUE)
        train_rmse <- sqrt(mean((train_preds$.pred - train_preds[[outcome]])^2, na.rm = TRUE))
        
        test_cal_slope <- test_cal_intercept <- NA
        if (!is.na(test_rsq) && is.finite(test_rsq)) {
          cal_model <- lm(as.formula(paste(outcome, "~ .pred")), data = test_preds)
          test_cal_intercept <- coef(cal_model)[1]
          test_cal_slope     <- coef(cal_model)[2]
        }
        
        tibble(
          model = model_name, outer_fold = fold_label,
          n_train = nrow(outer_train), n_test = nrow(outer_test),
          train_mae = train_mae, train_rmse = train_rmse, train_rsq = train_rsq,
          test_mae = test_mae, test_rmse = test_rmse, test_rsq = test_rsq,
          test_cal_slope = test_cal_slope, test_cal_intercept = test_cal_intercept
        )
      }, error = function(e) {
        tibble(
          model = model_name, outer_fold = fold_label,
          n_train = nrow(outer_train), n_test = nrow(outer_test),
          train_mae = NA, train_rmse = NA, train_rsq = NA,
          test_mae = NA, test_rmse = NA, test_rsq = NA,
          test_cal_slope = NA, test_cal_intercept = NA
        )
      })
      
      fold_metrics[[model_name]] <- fit_result
    }
    
    outer_fold_results[[fold_label]] <- bind_rows(fold_metrics)
    cat("done\n")
  }
  
  # ---- Aggregate outer-fold metrics ----
  all_outer_metrics <- bind_rows(outer_fold_results)
  
  summary_metrics <- all_outer_metrics %>%
    group_by(model) %>%
    summarise(
      n_folds    = n(),
      n_train    = round(mean(n_train, na.rm = TRUE)),
      n_test     = round(mean(n_test, na.rm = TRUE)),
      train_mae  = round(mean(train_mae, na.rm = TRUE), 3),
      train_rmse = round(mean(train_rmse, na.rm = TRUE), 3),
      train_rsq  = round(mean(train_rsq, na.rm = TRUE), 3),
      q2         = round(mean(test_rsq, na.rm = TRUE), 3),
      q2_sd      = round(sd(test_rsq, na.rm = TRUE), 3),
      cv_mae     = round(mean(test_mae, na.rm = TRUE), 3),
      cv_mae_sd  = round(sd(test_mae, na.rm = TRUE), 3),
      cv_rmse    = round(mean(test_rmse, na.rm = TRUE), 3),
      cv_rmse_sd = round(sd(test_rmse, na.rm = TRUE), 3),
      cv_cal_slope     = round(mean(test_cal_slope, na.rm = TRUE), 3),
      cv_cal_intercept = round(mean(test_cal_intercept, na.rm = TRUE), 3),
      # FIX: Q²/R² ratio (was R²/Q²)
      q2_r2_ratio  = round(
        ifelse(mean(train_rsq, na.rm = TRUE) > 0.01,
               mean(test_rsq, na.rm = TRUE) / mean(train_rsq, na.rm = TRUE), NA), 2),
      .groups = "drop"
    ) %>%
    mutate(dataset = dataset_name, outcome = outcome,
           n_predictors = n_predictors, n_total = n_total) %>%
    select(dataset, outcome, model, n_total, n_train, n_test, n_predictors, everything())
  
  # ---- Final refit on full data for tuning diagnostics ----
  cat("    Final refit on full data for tuning plots...\n")
  final_inner <- tryCatch(run_inner_tuning(df_clean), error = function(e) NULL)
  
  list(
    summary_metrics    = summary_metrics,
    outer_fold_metrics = all_outer_metrics,
    final_tuning       = final_inner,
    dataset_name       = dataset_name,
    outcome            = outcome
  )
}



# ============================================
# 11. Build task list (no composite)
# ============================================

cat("\n=== BUILDING TASK LIST ===\n")

CVD_scores_final <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y")

task_list <- list()
task_counter <- 1

for (score_name in CVD_scores_final) {
  task_list[[task_counter]] <- list(
    df = score_specific_datasets[[score_name]],
    outcome = score_name,
    dataset_name = paste0(score_name, "_specific"),
    category = "score_specific"
  )
  task_counter <- task_counter + 1
}

for (outcome in CVD_scores_final) {
  task_list[[task_counter]] <- list(
    df = datasets_list[["all_data"]],
    outcome = outcome,
    dataset_name = "full_predictor_set",
    category = "predictor_block"
  )
  task_counter <- task_counter + 1
}

cat(sprintf("Total tasks: %d\n", length(task_list)))

# ============================================
# 12. Full run with incremental saving
# ============================================

plan(multisession, workers = CPUS)

checkpoint_file <- file.path(output_dir, "nested_cv_checkpoint.qs2")
if (file.exists(checkpoint_file)) {
  checkpoint <- qs2::qs_read(checkpoint_file)
  all_results  <- checkpoint$all_results
  failed_tasks <- checkpoint$failed_tasks
  start_idx    <- checkpoint$last_completed + 1
  cat(sprintf("Resuming from checkpoint: task %d/%d\n", start_idx, length(task_list)))
} else {
  all_results  <- list()
  failed_tasks <- c()
  start_idx    <- 1
}

start_time <- Sys.time()

for (i in start_idx:length(task_list)) {
  task <- task_list[[i]]
  task_name <- paste0(task$outcome, " ~ ", task$dataset_name)
  
  cat(sprintf("\n[%d/%d] %s\n", i, length(task_list), task_name))
  
  result <- tryCatch({
    run_nested_cv(
      df = task$df, outcome = task$outcome, dataset_name = task$dataset_name,
      model_specs = model_specs,
      outer_folds = OUTER_FOLDS, outer_repeats = OUTER_REPEATS,
      inner_folds = INNER_FOLDS, grid_size = GRID_SIZE, seed = 42
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
  
  qs2::qs_save(
    list(all_results = all_results, failed_tasks = failed_tasks, last_completed = i),
    checkpoint_file
  )
  cat(sprintf("  [Checkpoint saved: %d/%d tasks]\n", i, length(task_list)))
}

plan(sequential)

total_time <- round(difftime(Sys.time(), start_time, units = "hours"), 2)
cat(sprintf("\nCompleted in %.2f hours\n", total_time))
cat(sprintf("Successful: %d | Failed: %d\n", length(all_results), length(failed_tasks)))

qs2::qs_save(all_results, file.path(output_dir, "nested_cv_all_results.qs2"))
if (file.exists(checkpoint_file)) file.remove(checkpoint_file)
if (length(failed_tasks) > 0) {
  writeLines(failed_tasks, file.path(output_dir, "nested_cv_failed_tasks.txt"))
}

# ============================================
# 13. Compile results
# ============================================

all_metrics <- bind_rows(lapply(all_results, function(x) x$summary_metrics))

all_outer_folds <- bind_rows(lapply(names(all_results), function(task_name) {
  x <- all_results[[task_name]]
  x$outer_fold_metrics %>%
    mutate(dataset = x$dataset_name, outcome_name = x$outcome)
}))


model_rankings <- all_metrics %>%
  group_by(dataset, outcome) %>%
  mutate(rank = rank(cv_mae)) %>%
  ungroup() %>%
  group_by(model) %>%
  summarise(
    mean_rank      = round(mean(rank), 2),
    times_best     = sum(rank == 1),
    mean_q2        = round(mean(q2, na.rm = TRUE), 3),
    mean_cv_mae    = round(mean(cv_mae, na.rm = TRUE), 3),
    mean_q2_r2_ratio = round(mean(q2_r2_ratio, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(mean_rank)

cat("\n--- Model Rankings ---\n")
print(model_rankings)

excel_output <- list(
  All_Metrics       = all_metrics,
  Model_Rankings    = model_rankings,
  Outer_Fold_Detail = all_outer_folds
)

wb <- createWorkbook()
for (sheet_name in names(excel_output)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, excel_output[[sheet_name]])
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = 1:ncol(excel_output[[sheet_name]]), widths = "auto")
}

saveWorkbook(wb, file.path(output_dir, "nested_cv_model_comparison_results.xlsx"), overwrite = TRUE)
qs2::qs_save(list(
  metrics        = all_metrics,
  outer_folds    = all_outer_folds,
  model_rankings = model_rankings
), file.path(output_dir, "nested_cv_compiled_objects.qs2"))



# ============================================
# 14. Combined comparison table (with Q²/R² ratio)
# ============================================

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
  "score_specific"     = "Score-specific inputs",
  "full_predictor_set" = "Full predictor set"
)

format_block <- function(metrics_df, score, dataset_name) {
  metrics_df %>%
    filter(outcome == score, dataset == dataset_name) %>%
    mutate(model = factor(model, levels = model_order)) %>%
    arrange(model) %>%
    transmute(
      Model      = model_labels_table[as.character(model)],
      `Train R²` = sprintf("%.3f", train_rsq),
      `Q²`       = sprintf("%.3f \u00b1 %.3f", q2, q2_sd),
      `CV MAE`   = sprintf("%.3f \u00b1 %.3f", cv_mae, cv_mae_sd),
      # FIX: Q²/R² ratio instead of R²/Q²
      `Q²/R² Ratio` = ifelse(!is.na(q2_r2_ratio), sprintf("%.2f", q2_r2_ratio), "\u2013")
    )
}

score_pairs <- list(
  c("QRISK3_risk", "SCORE2_score"),
  c("frs_10y", "ascvd_10y")
)

all_rows <- tibble()
label_row_indices <- c()
dataset_label_indices <- c()

for (pair in score_pairs) {
  left_score  <- pair[1]
  right_score <- pair[2]
  
  label_row <- tibble(
    Model_L = outcome_labels[left_score],
    TrainR2_L = "", Q2_L = "", CVMAE_L = "", Ratio_L = "",
    Spacer = "",
    Model_R = outcome_labels[right_score],
    TrainR2_R = "", Q2_R = "", CVMAE_R = "", Ratio_R = ""
  )
  label_row_indices <- c(label_row_indices, nrow(all_rows) + 1)
  all_rows <- bind_rows(all_rows, label_row)
  
  for (ds_type in c("score_specific", "full_predictor_set")) {
    ds_label <- dataset_display[ds_type]
    ds_label_row <- tibble(
      Model_L = ds_label,
      TrainR2_L = "", Q2_L = "", CVMAE_L = "", Ratio_L = "",
      Spacer = "",
      Model_R = ds_label,
      TrainR2_R = "", Q2_R = "", CVMAE_R = "", Ratio_R = ""
    )
    dataset_label_indices <- c(dataset_label_indices, nrow(all_rows) + 1)
    all_rows <- bind_rows(all_rows, ds_label_row)
    
    if (ds_type == "score_specific") {
      left_ds  <- paste0(left_score, "_specific")
      right_ds <- paste0(right_score, "_specific")
    } else {
      left_ds  <- "full_predictor_set"
      right_ds <- "full_predictor_set"
    }
    
    df_left  <- format_block(all_metrics, left_score, left_ds)
    df_right <- format_block(all_metrics, right_score, right_ds)
    
    data_rows <- tibble(
      Model_L   = df_left$Model,
      TrainR2_L = df_left$`Train R²`,
      Q2_L      = df_left$`Q²`,
      CVMAE_L   = df_left$`CV MAE`,
      Ratio_L   = df_left$`Q²/R² Ratio`,
      Spacer    = rep("", nrow(df_left)),
      Model_R   = df_right$Model,
      TrainR2_R = df_right$`Train R²`,
      Q2_R      = df_right$`Q²`,
      CVMAE_R   = df_right$`CV MAE`,
      Ratio_R   = df_right$`Q²/R² Ratio`
    )
    all_rows <- bind_rows(all_rows, data_rows)
  }
}

ft <- flextable(all_rows) %>%
  set_header_labels(
    Model_L = "Model", TrainR2_L = "Train R\u00b2", Q2_L = "Q\u00b2",
    CVMAE_L = "CV MAE", Ratio_L = "Q\u00b2/R\u00b2",
    Spacer = "",
    Model_R = "Model", TrainR2_R = "Train R\u00b2", Q2_R = "Q\u00b2",
    CVMAE_R = "CV MAE", Ratio_R = "Q\u00b2/R\u00b2"
  ) %>%
  style(part = "body",
        pr_t = fp_text(font.size = 9, font.family = "Arial")) %>%
  style(part = "header",
        pr_t = fp_text(font.size = 9, font.family = "Arial", bold = TRUE)) %>%
  style(i = label_row_indices, j = c(1, 7), part = "body",
        pr_t = fp_text(font.size = 11, font.family = "Arial", bold = TRUE)) %>%
  style(i = dataset_label_indices, j = c(1, 7), part = "body",
        pr_t = fp_text(font.size = 9, font.family = "Arial", italic = TRUE)) %>%
  align(j = c(1, 7), align = "left", part = "all") %>%
  align(j = c(2:5, 8:11), align = "center", part = "all") %>%
  align(j = 6, align = "center", part = "all") %>%
  width(j = 6, width = 0.15) %>%
  width(j = c(1, 7), width = 0.9) %>%
  width(j = c(2, 8), width = 0.7) %>%
  width(j = c(3, 9), width = 1.1) %>%
  width(j = c(4, 10), width = 1.1) %>%
  width(j = c(5, 11), width = 0.7) %>%
  border_remove() %>%
  hline_top(border = fp_border(width = 2), part = "all") %>%
  hline_bottom(border = fp_border(width = 2), part = "body") %>%
  hline(i = 1, border = fp_border(width = 1), part = "header")

for (idx in label_row_indices) {
  ft <- ft %>% hline(i = idx, border = fp_border(width = 1.5), part = "body")
}
for (idx in dataset_label_indices) {
  ft <- ft %>% hline(i = idx - 1, border = fp_border(width = 0.5, color = "grey60"), part = "body")
}

doc <- read_docx() %>%
  body_end_section_landscape() %>%
  body_add_flextable(value = ft)
print(doc, target = file.path(output_dir, "Table_Model_Comparison_NestedCV.docx"))
cat("Saved Table_Model_Comparison_NestedCV.docx\n")

# ============================================
# 15. Bar plots and tuning diagnostics
# ============================================

# ---- Shared labels and colours (defined once) ----
model_labels <- c(elastic_net = "Elastic Net", sparse_pls = "Sparse PLS",
                  random_forest = "Random Forest", xgboost = "XGBoost",
                  svm = "SVM-RBF", knn = "k-NN")

outcome_labels <- c(QRISK3_risk = "QRISK3 Risk", SCORE2_score = "SCORE2 Risk",
                    frs_10y = "Framingham 10-year Risk", ascvd_10y = "ASCVD 10-year Risk")

model_colors <- c("Elastic Net" = "#E69F00", "Sparse PLS" = "#56B4E9",
                  "Random Forest" = "#009E73", "XGBoost" = "#F0E442",
                  "SVM-RBF" = "#0072B2", "k-NN" = "#D55E00")

param_display_names <- c(
  penalty = "Penalty (\u03bb)", num_comp = "Number of Components",
  mtry = "Variables per Split (mtry)", learn_rate = "Learning Rate",
  cost = "Cost (C)", neighbors = "Neighbors (k)"
)

model_plot_order <- c("elastic_net", "sparse_pls", "random_forest",
                      "xgboost", "svm", "knn")

score_specific_names <- c("QRISK3_risk_specific", "SCORE2_score_specific",
                          "frs_10y_specific", "ascvd_10y_specific")

# ---- Bar plots: Score-specific ----
ss_data <- all_metrics %>%
  filter(dataset %in% score_specific_names, outcome %in% names(outcome_labels)) %>%
  mutate(outcome_label = factor(outcome_labels[outcome], levels = outcome_labels),
         model_label = factor(model_labels[model], levels = model_labels))

p1 <- ss_data %>%
  ggplot(aes(x = train_rsq, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Training R\u00b2", x = expression(R^2), y = NULL, fill = "Model") +
  theme_minimal() +
  theme(legend.position = "none", panel.grid.major.y = element_blank())

p2 <- ss_data %>%
  ggplot(aes(x = q2, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(aes(xmin = q2 - q2_sd, xmax = q2 + q2_sd),
                 position = position_dodge(width = 0.8), height = 0.3) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Q\u00b2", x = expression(Q^2), y = NULL, fill = "Model") +
  theme_minimal() +
  theme(legend.position = "right", panel.grid.major.y = element_blank())

combined_ss <- (p1 | p2) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Model Performance: Score-Specific Inputs (Nested CV)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

ggsave(file.path(output_dir, "barplot_score_specific.png"), combined_ss, width = 12, height = 6, dpi = 300)
ggsave(file.path(output_dir, "barplot_score_specific.pdf"), combined_ss, width = 12, height = 6)

# ---- Bar plots: Full predictor set ----
fps_data <- all_metrics %>%
  filter(dataset == "full_predictor_set", outcome %in% names(outcome_labels)) %>%
  mutate(outcome_label = factor(outcome_labels[outcome], levels = outcome_labels),
         model_label = factor(model_labels[model], levels = model_labels))

p3 <- fps_data %>%
  ggplot(aes(x = train_rsq, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Training R\u00b2", x = expression(R^2), y = NULL, fill = "Model") +
  theme_minimal() +
  theme(legend.position = "none", panel.grid.major.y = element_blank())

p4 <- fps_data %>%
  ggplot(aes(x = q2, y = outcome_label, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(aes(xmin = q2 - q2_sd, xmax = q2 + q2_sd),
                 position = position_dodge(width = 0.8), height = 0.3) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Q\u00b2", x = expression(Q^2), y = NULL, fill = "Model") +
  theme_minimal() +
  theme(legend.position = "right", panel.grid.major.y = element_blank())

combined_fps <- (p3 | p4) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Model Performance: Full Predictor Set (Nested CV)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

ggsave(file.path(output_dir, "barplot_full_predictor_set.png"), combined_fps, width = 12, height = 6, dpi = 300)
ggsave(file.path(output_dir, "barplot_full_predictor_set.pdf"), combined_fps, width = 12, height = 6)

cat("Bar plots saved\n")

# ============================================
# 16. Tuning parameter diagnostics
# ============================================

cat("\n=== TUNING PARAMETER DIAGNOSTICS ===\n")

primary_tuning_params <- c(
  recipe_elastic_net   = "penalty",
  recipe_sparse_pls    = "num_comp",
  recipe_random_forest = "mtry",
  recipe_xgboost       = "learn_rate",
  recipe_svm           = "cost",
  recipe_knn           = "neighbors"
)

# ---- Extract tuning data from final full-data refit ----
extract_tuning_data <- function(result) {
  if (is.null(result$final_tuning)) return(NULL)
  wf_res <- result$final_tuning$wf_results
  
  map_dfr(wf_res$wflow_id, function(wf_id) {
    primary_param <- primary_tuning_params[wf_id]
    if (is.na(primary_param)) return(NULL)
    
    tryCatch({
      res     <- wf_res %>% extract_workflow_set_result(wf_id)
      metrics <- res %>% collect_metrics() %>% filter(.metric == "rmse")
      if (!primary_param %in% names(metrics)) return(NULL)
      
      model_name  <- str_remove(wf_id, "recipe_")
      best_config <- select_one_se(res, wf_id)
      
      metrics %>%
        mutate(
          model              = model_name,
          model_label        = model_labels[model_name],
          dataset            = result$dataset_name,
          outcome            = result$outcome,
          tuning_param_name  = primary_param,
          tuning_param_value = .data[[primary_param]],
          is_selected        = (.config == best_config$.config)
        ) %>%
        select(model, model_label, dataset, outcome,
               tuning_param_name, tuning_param_value,
               cv_rmse = mean, cv_rmse_se = std_err, is_selected)
    }, error = function(e) NULL)
  })
}

all_tuning_data <- map_dfr(all_results, extract_tuning_data)

if (nrow(all_tuning_data) > 0) {
  
  all_tuning_data <- all_tuning_data %>%
    filter(outcome %in% names(outcome_labels)) %>%
    mutate(
      outcome_label = factor(outcome_labels[outcome], levels = outcome_labels),
      model_label   = factor(model_label, levels = model_labels)
    )
  
  # ---- Overview plot helper (multi-outcome faceted) ----
  create_model_tuning_plot <- function(data, model_name, plot_subtitle) {
    model_data <- data %>% filter(model == model_name)
    if (nrow(model_data) == 0) return(NULL)
    
    param_name   <- unique(model_data$tuning_param_name)
    x_label      <- param_display_names[param_name]
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
           x = x_label, y = "Inner CV RMSE") +
      theme_minimal(base_size = 10) +
      theme(strip.text    = element_text(face = "bold", size = 9),
            panel.grid.minor = element_blank(),
            plot.title    = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 9, colour = "grey40"))
    
    if (param_name %in% c("penalty", "cost")) {
      p <- p + scale_x_continuous(
        trans = "log10",
        breaks = 10^seq(-10, 10, by = 2),
        labels = function(x) parse(text = paste0("10^", log10(x)))
      )
    }
    p
  }
  
  # ---- Thesis-style single-score tuning plot (defined once) ----
  create_thesis_tuning_plot <- function(data, model_name) {
    model_data <- data %>% filter(model == model_name)
    if (nrow(model_data) == 0) return(ggplot() + theme_void())
    
    param_name   <- unique(model_data$tuning_param_name)
    x_label      <- param_display_names[param_name]
    selected_pts <- model_data %>% filter(is_selected)
    
    p <- ggplot(model_data, aes(x = tuning_param_value, y = cv_rmse)) +
      geom_line(colour = "grey60", linewidth = 0.5) +
      geom_point(colour = "grey40", size = 2, alpha = 0.7) +
      geom_errorbar(aes(ymin = cv_rmse - cv_rmse_se, ymax = cv_rmse + cv_rmse_se),
                    colour = "grey70", alpha = 0.5, width = 0) +
      geom_point(data = selected_pts, colour = "#E41A1C", size = 3.5, shape = 18) +
      geom_vline(data = selected_pts, aes(xintercept = tuning_param_value),
                 colour = "#E41A1C", linetype = "dashed", alpha = 0.5) +
      labs(title = model_labels[model_name], x = x_label, y = "Inner CV RMSE") +
      theme_minimal(base_size = 16) +
      theme(
        text             = element_text(family = "Arial", size = 16),
        plot.title       = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title       = element_text(size = 15),
        axis.text        = element_text(size = 13),
        strip.text       = element_text(size = 14, face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    if (param_name %in% c("penalty", "cost")) {
      p <- p + scale_x_log10(
        breaks = 10^seq(-10, 10, by = 2),
        labels = function(x) parse(text = paste0("10^", log10(x)))
      )
    }
    p
  }
  
  # ---- 16a. Per-dataset overview plots ----
  for (ds_pattern in c("_specific", "full_predictor_set")) {
    if (ds_pattern == "_specific") {
      ds_data   <- all_tuning_data %>% filter(str_detect(dataset, "_specific"))
      ds_suffix <- "score_specific"
      ds_title  <- "Score-Specific Inputs"
    } else {
      ds_data   <- all_tuning_data %>% filter(dataset == "full_predictor_set")
      ds_suffix <- "full_predictor_set"
      ds_title  <- "Full Predictor Set"
    }
    
    if (nrow(ds_data) > 0) {
      plots_list <- list()
      for (m in model_plot_order) {
        p <- create_model_tuning_plot(ds_data, m, ds_title)
        if (!is.null(p)) plots_list[[m]] <- p
      }
      
      if (length(plots_list) > 0) {
        combined <- wrap_plots(plots_list, ncol = 6) +
          plot_annotation(
            title    = paste("Tuning Parameter Diagnostics:", ds_title),
            subtitle = "Red diamond = selected configuration (one-SE rule on inner CV RMSE)",
            theme = theme(
              plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
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
  
  # ---- 16b. Per-score detailed tuning plots ----
  for (score in names(outcome_labels)) {
    score_label       <- outcome_labels[score]
    score_specific_ds <- paste0(score, "_specific")
    
    ss_tune  <- all_tuning_data %>% filter(outcome == score, dataset == score_specific_ds)
    fps_tune <- all_tuning_data %>% filter(outcome == score, dataset == "full_predictor_set")
    
    if (nrow(ss_tune) == 0 && nrow(fps_tune) == 0) next
    
    row1_plots <- list()
    row2_plots <- list()
    for (m in model_plot_order) {
      p1 <- create_thesis_tuning_plot(ss_tune, m)
      row1_plots[[m]] <- if (!is.null(p1)) p1 else ggplot() + theme_void()
      
      p2 <- create_thesis_tuning_plot(fps_tune, m)
      row2_plots[[m]] <- if (!is.null(p2)) p2 else ggplot() + theme_void()
    }
    
    # Add A/B tags to first plot in each row
    row1_plots[[1]] <- row1_plots[[1]] +
      labs(tag = "A") +
      theme(plot.tag = element_text(face = "bold", family = "Arial", size = 18))
    
    row2_plots[[1]] <- row2_plots[[1]] +
      labs(tag = "B") +
      theme(plot.tag = element_text(face = "bold", family = "Arial", size = 18))
    
    top_row    <- wrap_plots(row1_plots, ncol = 6)
    bottom_row <- wrap_plots(row2_plots, ncol = 6)
    
    combined_score <- top_row / bottom_row +
      plot_annotation(
        title    = paste("Tuning Parameter Diagnostics:", score_label),
        subtitle = "A: Score-specific inputs | B: Full predictor set\nRed diamond = selected (one-SE rule on inner CV RMSE)",
        theme = theme(
          plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = 10)
        )
      )
    
    ggsave(file.path(output_dir, paste0("tuning_detailed_", score, ".png")),
           combined_score, width = 26, height = 11, dpi = 300)
    ggsave(file.path(output_dir, paste0("tuning_detailed_", score, ".pdf")),
           combined_score, width = 26, height = 11)
    cat(sprintf("  Saved tuning_detailed_%s\n", score))
  }
  
} else {
  cat("  No tuning data available\n")
}

# ============================================
# 17. Selected hyperparameter table
# ============================================

cat("\n=== SELECTED HYPERPARAMETER TABLE ===\n")

# Parameter display names for the table
param_labels <- c(
  penalty        = "\u03bb (penalty)",
  mixture        = "\u03b1 (mixture)",
  num_comp       = "Components",
  predictor_prop = "Predictor prop.",
  trees          = "Trees",
  mtry           = "mtry",
  min_n          = "Min node size",
  tree_depth     = "Tree depth",
  loss_reduction = "Loss reduction",
  sample_size    = "Subsample prop.",
  learn_rate     = "Learning rate",
  cost           = "Cost (C)",
  rbf_sigma      = "\u03c3 (RBF)",
  neighbors      = "Neighbours (k)",
  weight_func    = "Weight function",
  dist_power     = "Minkowski p"
)

# Define which parameters each model has (in display order)
model_param_map <- list(
  elastic_net   = c("penalty", "mixture"),
  sparse_pls    = c("num_comp", "predictor_prop"),
  random_forest = c("trees", "mtry", "min_n"),
  xgboost       = c("trees", "tree_depth", "min_n", "loss_reduction",
                    "sample_size", "mtry", "learn_rate"),
  svm           = c("cost", "rbf_sigma"),
  knn           = c("neighbors", "weight_func", "dist_power")
)

# Extract best params from final_tuning (full-data refit) for each task
extract_best_params <- function(result) {
  if (is.null(result$final_tuning)) return(NULL)
  best_params <- result$final_tuning$best_params
  
  map_dfr(names(best_params), function(model_name) {
    bp <- best_params[[model_name]]
    if (is.null(bp)) return(NULL)
    
    expected_params <- model_param_map[[model_name]]
    
    param_values <- sapply(expected_params, function(p) {
      if (p %in% names(bp)) {
        val <- bp[[p]]
        if (is.numeric(val)) {
          if (abs(val) < 0.001 & val != 0) {
            sprintf("%.2e", val)
          } else if (val == round(val)) {
            sprintf("%d", as.integer(val))
          } else {
            sprintf("%.4f", val)
          }
        } else {
          as.character(val)
        }
      } else {
        NA_character_
      }
    })
    
    tibble(
      dataset  = result$dataset_name,
      outcome  = result$outcome,
      model    = model_name,
      param    = expected_params,
      param_label = param_labels[expected_params],
      value    = param_values
    )
  })
}

all_best_params <- map_dfr(all_results, extract_best_params)

if (nrow(all_best_params) > 0) {
  
  # ---- Build one table per dataset type ----
  for (ds_type in c("score_specific", "full_predictor_set")) {
    
    if (ds_type == "score_specific") {
      ds_data <- all_best_params %>% filter(str_detect(dataset, "_specific"))
      ds_title <- "Score-Specific Inputs"
      ds_suffix <- "score_specific"
    } else {
      ds_data <- all_best_params %>% filter(dataset == "full_predictor_set")
      ds_title <- "Full Predictor Set"
      ds_suffix <- "full_predictor_set"
    }
    
    if (nrow(ds_data) == 0) next
    
    # Pivot: rows = model + parameter, columns = outcome scores
    score_order <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y")
    score_labels_short <- c(
      QRISK3_risk = "QRISK3", SCORE2_score = "SCORE2",
      frs_10y = "Framingham", ascvd_10y = "ASCVD"
    )
    
    model_order_local <- c("elastic_net", "sparse_pls", "random_forest",
                           "xgboost", "svm", "knn")
    model_labels_local <- c(
      elastic_net = "Elastic Net", sparse_pls = "Sparse PLS",
      random_forest = "Random Forest", xgboost = "XGBoost",
      svm = "SVM-RBF", knn = "k-NN"
    )
    
    # Build table rows
    table_rows <- tibble()
    model_header_indices <- c()
    
    for (m in model_order_local) {
      m_data <- ds_data %>% filter(model == m)
      if (nrow(m_data) == 0) next
      
      # Model header row
      header_row <- tibble(Parameter = model_labels_local[m])
      for (s in score_order) {
        header_row[[score_labels_short[s]]] <- ""
      }
      model_header_indices <- c(model_header_indices, nrow(table_rows) + 1)
      table_rows <- bind_rows(table_rows, header_row)
      
      # Parameter rows
      params_for_model <- model_param_map[[m]]
      for (p in params_for_model) {
        param_row <- tibble(Parameter = paste0("    ", param_labels[p]))
        for (s in score_order) {
          val <- m_data %>% filter(outcome == s, param == p) %>% pull(value)
          param_row[[score_labels_short[s]]] <- if (length(val) == 1 && !is.na(val)) val else "\u2013"
        }
        table_rows <- bind_rows(table_rows, param_row)
      }
    }
    
    # Create flextable
    ft_params <- flextable(table_rows) %>%
      font(fontname = "Arial", part = "all") %>%
      style(part = "body",
            pr_t = fp_text(font.size = 9, font.family = "Arial")) %>%
      style(part = "header",
            pr_t = fp_text(font.size = 9, font.family = "Arial", bold = TRUE)) %>%
      style(i = model_header_indices, j = 1, part = "body",
            pr_t = fp_text(font.size = 9, font.family = "Arial", bold = TRUE)) %>%
      style(i = setdiff(seq_len(nrow(table_rows)), model_header_indices),
            j = 1, part = "body",
            pr_t = fp_text(font.size = 9, font.family = "Arial", italic = TRUE)) %>%
      align(j = 1, align = "left", part = "all") %>%
      align(j = 2:5, align = "center", part = "all") %>%
      width(j = 1, width = 1.8) %>%
      width(j = 2:5, width = 1.2) %>%
      border_remove() %>%
      hline_top(border = fp_border(width = 2), part = "all") %>%
      hline_bottom(border = fp_border(width = 2), part = "body") %>%
      hline(i = 1, border = fp_border(width = 1), part = "header")
    
    for (idx in model_header_indices) {
      ft_params <- ft_params %>%
        hline(i = idx, border = fp_border(width = 0.5, color = "grey60"), part = "body")
    }
    
    # Save
    doc_params <- read_docx() %>%
      body_add_par(paste("Selected Hyperparameters:", ds_title),
                   style = "heading 2") %>%
      body_add_flextable(value = ft_params)
    print(doc_params, target = file.path(output_dir,
                                         paste0("Table_Selected_Hyperparameters_", ds_suffix, ".docx")))
    cat(sprintf("  Saved Table_Selected_Hyperparameters_%s.docx\n", ds_suffix))
  }
  
  # ---- Also save as Excel for reference ----
  params_wide <- all_best_params %>%
    filter(outcome %in% c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y")) %>%
    mutate(
      model_label = model_labels[model],
      outcome_label = c(QRISK3_risk = "QRISK3", SCORE2_score = "SCORE2",
                        frs_10y = "Framingham", ascvd_10y = "ASCVD")[outcome]
    ) %>%
    select(dataset, outcome_label, model_label, param_label, value) %>%
    pivot_wider(names_from = outcome_label, values_from = value)
  
  wb_hp <- createWorkbook()
  addWorksheet(wb_hp, "Selected_Hyperparameters")
  writeData(wb_hp, "Selected_Hyperparameters", params_wide)
  freezePane(wb_hp, "Selected_Hyperparameters", firstRow = TRUE)
  setColWidths(wb_hp, "Selected_Hyperparameters", cols = 1:ncol(params_wide), widths = "auto")
  saveWorkbook(wb_hp, file.path(output_dir, "selected_hyperparameters.xlsx"), overwrite = TRUE)
  
  cat("  Saved selected_hyperparameters.xlsx\n")
  
} else {
  cat("  No tuning results available for hyperparameter table\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")