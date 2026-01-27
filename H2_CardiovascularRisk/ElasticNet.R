### Elastic Net Analysis - CoDiet CVD Prediction
### Author: Luisa Delius

### Table of Content
# Step 1: Load data
# Step 2: create Data Sets
# 3 Missing Data Imputation


# Load packages
library(glmnet)      # Elastic net regression
library(dplyr)       # Data manipulation
library(tidyr)       # Data reshaping
library(ggplot2)     # Visualization
library(openxlsx)

set.seed(42)

# Set working directory
wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files/processed_data"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

# ============================================
## 1. Uploading the data and keeping only the needed columns 
# ============================================
df_fatty_acids_ElaNet <- readRDS("df_fatty_acids_predictor_statin_suppl.rds") %>%
  select(Sample_ID, Statins, Supplements, starts_with("z_"))

df_lipidomics_ElaNet <- readRDS("df_lipidomics_predictor_statin_suppl.rds") %>%
  select(Sample_ID, starts_with("z_"))

df_urine_nmr_ElaNet <- readRDS("df_urine_NMR_data.rds") %>%
  select(Sample_ID, starts_with("z_"))

df_body_composition_ElaNet <- readRDS("df_body_composition_metrics.rds") %>%
  select(Sample_ID, starts_with("z_")) %>%
  filter(!is.na(Sample_ID))

df_REDcap_demographics_ElaNet <- readRDS("df_REDcap_demographics_ElaNet.rds") %>%
  select(Sample_ID, starts_with("z_"), Marital_status, Living_Status, Employment_Status,
         annual_net_salary, Working_time, Education_Level, naps_during_day,
         olive_oil_as_main_culinary_fat, olive_oil_given_day, wine_per_week, servings_nuts_per_week)
   
df_risk_factors_ElaNet <- readRDS("df_risk_factor_predictors.rds") %>%  
  select(Sample_ID, Country, starts_with("z_")) %>% #HERE KEEP THE COUNTRY!!
  select(-z_Heart.Rate, -z_Age, -z_Total.Cholesterol.mg.dl, -z_HDL.mg.dl,
         -z_Body.Weight, -z_Height, -z_BMI) # remove the factors we dont want to include. I removed heart rate here as its also included in body composition metrics and i dont want to use it twice.
  
df_all_cvd_risk_scores_ElaNet <- readRDS("df_all_risk_scores.rds") %>% 
  rename(Sample_ID = PatientID) %>%
  select(-SCORE2_strat)

df_ascvd_frs_score2_input <- readRDS("ASCVD_SCORE2_Framingham_input.rds") %>%
  rename(Sample_ID = PatientID)

df_QRISK3_input <-  readRDS("QRISK3_calculation_input.rds") %>%
  rename(Sample_ID = PatientID)

CVD_scores <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y", "mean_risk")

# ============================================
## 2. Create Data Sets
# ============================================
# 2.1 Model 1: Metabolites only
df_model1 <- df_fatty_acids_ElaNet %>% 
  select(-Statins, -Supplements) %>%
  full_join(df_lipidomics_ElaNet, by = "Sample_ID") %>%
  full_join(df_urine_nmr_ElaNet, by = "Sample_ID") %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID") %>%
  select(where(~!all(is.na(.)))) # Remove completely empty columns (100% missing)
    
# 2.2: Model 2: Metabolites + sociodemographics and clinical risk factors (=ALL data)
df_model2 <- df_model1 %>%
  full_join(df_fatty_acids_ElaNet %>% select(Sample_ID, Statins, Supplements), by = "Sample_ID") %>%
  full_join(df_body_composition_ElaNet, by = "Sample_ID") %>%
  full_join(df_risk_factors_ElaNet, by = "Sample_ID") %>%
  full_join(df_REDcap_demographics_ElaNet, by = "Sample_ID")

# 2.3: Model 3: Significant after multiple testing correction
BH_sig_predictors <- c("z_ecw_tbw", "z_mean_hrt", "z_ecm_bcm", "z_rbc_22_4n_6", 
                       "naps_during_day", "Living_Status")

df_model3 <- df_model2 %>%
  select(Sample_ID, all_of(BH_sig_predictors), all_of(CVD_scores))

# 2.4: Model 4: Significant after multiple testing correction (mtc) + min. 3x w/o mtc
sig_predictors <- c("z_AGE.reader", "z_Trigonelline", "z_Hba1C", "z_ALT.unit.L", 
                    "z_rbc_dpa_22_5n3", "z_rbc_eicosadienoic_20_2n6", 
                    "z_rbc_epa_20_5n3", "z_pe_o_19_1_20_5",
                    "z_3_hydroxybutyric_acid", "z_hippuric_acid", 
                    "Employment_Status", "Education_Level")

df_model4 <- df_model2 %>%
  select(Sample_ID,all_of(BH_sig_predictors), all_of(sig_predictors), all_of(CVD_scores))

# 2.5: Model 0: All columns that were used for the respective score calculation
df_model0_QRISK3 <- df_QRISK3_input %>%
  mutate(townsend = replace_na(townsend, 0)) %>% # set 0 for everyone we dont know the score cause thats how we calcualted QRISK3.
  full_join(df_all_cvd_risk_scores_ElaNet %>% select(Sample_ID, QRISK3_risk), 
            by = "Sample_ID") %>%
  filter(!is.na(QRISK3_risk))

df_model0_ascvd <- df_ascvd_frs_score2_input %>%
  select(-Risk.region, -mean_LDL_mg_dl) %>%
  full_join(df_all_cvd_risk_scores_ElaNet %>% select(Sample_ID, ascvd_10y), 
            by = "Sample_ID") %>%
  filter(!is.na(ascvd_10y))

df_model0_score2 <- df_ascvd_frs_score2_input %>%
  select(-mean_LDL_mg_dl, -blood_pressure_treatment, -race_ascvd) %>%
  full_join(df_all_cvd_risk_scores_ElaNet %>% select(Sample_ID, SCORE2_score), 
            by = "Sample_ID") %>%
  filter(!is.na(SCORE2_score))
  
df_model0_frs <- df_ascvd_frs_score2_input %>%
  select(-Risk.region, -mean_LDL_mg_dl, -race_ascvd) %>%
  full_join(df_all_cvd_risk_scores_ElaNet %>% select(Sample_ID, frs_10y), 
            by = "Sample_ID") %>%
  filter(!is.na(frs_10y))

df_model0_composite <-df_QRISK3_input %>%
  mutate(townsend = replace_na(townsend, 0)) %>%
  full_join(df_ascvd_frs_score2_input %>% select(Sample_ID, Risk.region),
            by = "Sample_ID") %>%
  full_join(df_all_cvd_risk_scores_ElaNet %>% select(Sample_ID, mean_risk), 
            by = "Sample_ID") %>%
  filter(!is.na(mean_risk))

# Score-specific datasets
score_specific_datasets <- list(
  QRISK3_risk = df_model0_QRISK3,
  SCORE2_score = df_model0_score2,
  frs_10y = df_model0_frs,
  ascvd_10y = df_model0_ascvd,
  mean_risk = df_model0_composite
)

# ============================================
# 3. Missing Data Imputation using Mean of each column (numeric) or mode (factor)
# ============================================
# DATA LEAKAGE! Imputation happenms before corss-validation splits. this causes
# minor data leakage cause mean/mode is calculated using all samples, including 
# the test fold. --> I have to impute inside each CV-fold.

impute_mean_mode <- function(df) {
    df %>%
        # Mean imputation for z-scored numeric columns
      mutate(across(starts_with("z_"), ~ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
        # Mode imputation for factor/character columns (except Country and Sample_ID)
      mutate(across(
          where(~(is.factor(.) | is.character(.)) & !all(is.na(.))),
          ~{
              col_name <- cur_column()
              if (col_name %in% c("Country", "Sample_ID")) {
                  .  # Don't impute these
                } else {
                    mode_val <- names(sort(table(., useNA = "ifany"), decreasing = TRUE))[1]
                    x_imputed <- as.character(.)
                    x_imputed[is.na(x_imputed)] <- mode_val
                    if (is.factor(.)) factor(x_imputed, levels = levels(.)) else x_imputed
                  }
            }
        ))
  }



# Apply imputation to all models
df_model1_imputed <- impute_mean_mode(df_model1)
df_model2_imputed <- impute_mean_mode(df_model2)
df_model3_imputed <- impute_mean_mode(df_model3)
df_model4_imputed <- impute_mean_mode(df_model4)

df_body_comp_imputed <- impute_mean_mode(df_body_composition_ElaNet) %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")

df_demographics_imputed <- df_REDcap_demographics_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID") %>%
  impute_mean_mode()

df_risk_factors_imputed <- df_risk_factors_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID") %>%
  select(where(~mean(is.na(.)) <= 0.5)) %>% # excluded columns with more then 50% missingness (could be even more strict I guess...)
  impute_mean_mode()

# list of all predictor block datasets
datasets_list <-  list(all_data = df_model2_imputed,
                       metabolites = df_model1_imputed,
                       mtc_sign = df_model3_imputed,
                       sign_wo_mtc = df_model4_imputed,
                       body_composition = df_body_comp_imputed,
                       sociodemographics_lifestyle = df_demographics_imputed,
                       clinical_risk_factors = df_risk_factors_imputed)

# ============================================
# 4. Helper function to build score-specific + 1 further predictor-block datasets
# ============================================
# Common blocks contain all CVD scores (must be dropped before join to avoid duplicates)
prep_extra_predictors <- function(extra_df) {
  keep_cols <- setdiff(names(extra_df), CVD_scores)
  extra_df[, keep_cols, drop = FALSE]
}

build_score_plus_block <- function(score, block_name,
                                   score_specific_datasets, datasets_list) {
  base_df  <- score_specific_datasets[[score]]
  extra_df <- datasets_list[[block_name]]
  extra_predictors <- prep_extra_predictors(extra_df)
  dplyr::left_join(base_df, extra_predictors, by = "Sample_ID")
}

# ============================================
# 5. CORE ELASTIC NET FUNCTIONS
# ============================================
get_oof_predictions_glmnet_nested <- function(X, y, alpha_grid, nfolds = 10,
                                              type_measure = "mse", foldid_outer = NULL,
                                              seed = 42) {
  set.seed(seed) 
  n <- length(y)
  
  if (is.null(foldid_outer)) {
    foldid_outer <- sample(rep(1:nfolds, length.out = n))
  }
  
  y_oof <- rep(NA_real_, n)
  alpha_selected_per_fold <- numeric(nfolds)  # ADD THIS LINE
  
  
  for (k in 1:nfolds) {
    idx_test <- which(foldid_outer == k)
    idx_train <- setdiff(1:n, idx_test)
    
    # Inner CV: grid search over alpha on training fold only
    cv_results_fold <- lapply(alpha_grid, function(a) {
      cv_fit <- cv.glmnet(
        x = X[idx_train, , drop = FALSE],
        y = y[idx_train],
        alpha = a,
        nfolds = max(3, nfolds - 1),
        type.measure = type_measure,
        standardize = TRUE,
        family = "gaussian"
      )
      list(
        alpha = a,
        cv_mse = min(cv_fit$cvm),
        cv_fit = cv_fit
      )
    })
    
    # Select best alpha for this fold
    cv_mse_vals <- sapply(cv_results_fold, function(x) x$cv_mse)
    best_idx <- which.min(cv_mse_vals)
    alpha_selected_per_fold[k] <- cv_results_fold[[best_idx]]$alpha  # ADD THIS LINE
    best_cv_fit_fold <- cv_results_fold[[best_idx]]$cv_fit
    
    # Predict on held-out fold
    y_oof[idx_test] <- as.numeric(
      predict(best_cv_fit_fold, newx = X[idx_test, , drop = FALSE], s = "lambda.min")
    )
  }
  
  return(list(predictions = y_oof, alpha_per_fold = alpha_selected_per_fold))
}




run_elastic_net_model <- function(
    X,                          # Predictor matrix (after imputation)
    y,                          # Outcome vector
    cvd_score_name,            # Name of CVD score
    dataset_name,              # Name of predictor dataset
    imputation_method,         # "mean_mode" 
    alpha_grid = seq(0, 1, by = 0.1),  # Alpha values to test
    nfolds = 10,
    type_measure = "mse",
    seed = 42
) {
  
  set.seed(seed)
  
  # Remove any rows with missing outcomes
  complete_cases <- complete.cases(y)
  X_clean <- X[complete_cases, ]
  y_clean <- y[complete_cases]
  
  n_obs <- length(y_clean)
  n_pred <- ncol(X_clean)
  predictor_names <- colnames(X_clean)
  
  cat(sprintf("\n  Sample size: %d, Predictors: %d\n", n_obs, n_pred))
  cat("  Testing alpha values:", paste(alpha_grid, collapse = ", "), "\n")
  
  # Fixed fold assignments for reproducibility
  set.seed(seed)
  foldid_inner <- sample(rep(1:nfolds, length.out = length(y_clean)))
  set.seed(seed + 1)
  foldid_outer <- sample(rep(1:nfolds, length.out = length(y_clean)))
  
  # Grid search over alpha values
  cv_results <- lapply(alpha_grid, function(a) {
    cv_fit <- cv.glmnet(
      x = X_clean,
      y = y_clean,
      alpha = a,
      nfolds = nfolds,
      foldid = foldid_inner,
      type.measure = type_measure,
      standardize = TRUE,
      family = "gaussian"
    )
    
    # Return CV error at lambda.min for this alpha
    min_cvm_idx <- which(cv_fit$lambda == cv_fit$lambda.min)
    list(
      alpha = a,
      cv_mse_min = cv_fit$cvm[min_cvm_idx],
      cv_fit = cv_fit
    )
  })
  
  # Find best alpha (lowest CV MSE)
  cv_mse_values <- sapply(cv_results, function(x) x$cv_mse_min)
  best_alpha_idx <- which.min(cv_mse_values)
  best_alpha_fulldata <- alpha_grid[best_alpha_idx]
  best_cv_fit <- cv_results[[best_alpha_idx]]$cv_fit
  
  # Extract lambda.min and lambda.1se
  lambda_min <- best_cv_fit$lambda.min
  lambda_1se <- best_cv_fit$lambda.1se
  cv_mse_min <- min(best_cv_fit$cvm)
  cv_mse_1se <- best_cv_fit$cvm[which(best_cv_fit$lambda == lambda_1se)]
  
  # Fit final model with optimal alpha at lambda.min
  final_coef <- coef(best_cv_fit, s = "lambda.min")
  intercept <- as.numeric(final_coef[1])
  coef_vector <- as.numeric(final_coef[-1])
  names(coef_vector) <- predictor_names
  
  # Count non-zero coefficients
  n_nonzero <- sum(coef_vector != 0)
  
  # ---- PART B: In-sample predictions (optimistic) ----
  y_pred_in <- as.numeric(predict(best_cv_fit, newx = X_clean, s = "lambda.min"))
  res_in <- y_clean - y_pred_in
  mse_in <- mean(res_in^2)
  rmse_in <- sqrt(mse_in)
  mae_in <- mean(abs(res_in))
  ss_total_in <- sum((y_clean - mean(y_clean))^2)
  r2_in <- 1 - sum(res_in^2) / ss_total_in
  cor_in <- cor(y_pred_in, y_clean)
  
  # ---- PART C: Out-of-fold predictions (honest, nested CV) ----
  oof_results <- get_oof_predictions_glmnet_nested(
    X = X_clean,
    y = y_clean,
    alpha_grid = alpha_grid,
    nfolds = nfolds,
    type_measure = type_measure,
    foldid_outer = foldid_outer,
    seed = seed + 100
  )
  
  y_pred_oof <- oof_results$predictions
  alpha_per_fold <- oof_results$alpha_per_fold
  
  res_oof <- y_clean - y_pred_oof
  mse_oof <- mean(res_oof^2)
  rmse_oof <- sqrt(mse_oof)
  mae_oof <- mean(abs(res_oof))
  ss_total_oof <- sum((y_clean - mean(y_clean))^2)
  r2_oof <- 1 - sum(res_oof^2) / ss_total_oof
  cor_oof <- cor(y_pred_oof, y_clean)
  
  
  # Deviance explained (from glmnet)
  dev_explained <- best_cv_fit$glmnet.fit$dev.ratio[which.min(abs(best_cv_fit$glmnet.fit$lambda - lambda_min))]
  
  # Extract top predictors (by absolute coefficient value)
  nonzero_coefs <- coef_vector[coef_vector != 0]
  if (length(nonzero_coefs) > 0) {
    sorted_coefs <- sort(abs(nonzero_coefs), decreasing = TRUE)
    top_n <- min(10, length(sorted_coefs))
    top_predictor_names <- names(sorted_coefs)[1:top_n]
    top_predictor_coefs <- nonzero_coefs[top_predictor_names]
  } else {
    top_predictor_names <- NA
    top_predictor_coefs <- NA
  }
  
  # Compile results
  result_row <- tibble(
    model_name = "elastic_net",
    cvd_score = cvd_score_name,
    dataset_name = dataset_name,
    
    # Sample characteristics
    n_observations = n_obs,
    n_predictors = n_pred,
    n_nonzero_coefs = n_nonzero,
    
    # Hyperparameters
    alpha_optimal_fulldata = best_alpha_fulldata,
    alpha_mean_nested = mean(alpha_per_fold),
    lambda_min = lambda_min,
    lambda_1se = lambda_1se,
    n_folds = nfolds,
    type_measure = type_measure,
    family = "gaussian",
    standardize = TRUE,
    
    # Imputation
    imputation_method = imputation_method,
    
    # Cross-validation performance
    cv_mse_min = cv_mse_min,
    cv_mse_1se = cv_mse_1se,
    
    # Full data performance (at lambda.min)
    # In-sample (optimistic)
    in_sample_mse = mse_in,
    in_sample_rmse = rmse_in,
    in_sample_mae = mae_in,
    in_sample_r2 = r2_in,
    in_sample_cor = cor_in,
    
    # Out-of-fold (honest)
    oof_mse = mse_oof,
    oof_rmse = rmse_oof,
    oof_mae = mae_oof,
    oof_r2 = r2_oof,
    oof_cor = cor_oof,
    
    pct_deviance_explained = dev_explained * 100,
    
    # Model characteristics
    intercept = intercept,
    top_10_predictors = list(top_predictor_names),
    top_10_coefficients = list(top_predictor_coefs),
    all_predictor_names = list(predictor_names),
    all_coefficients = list(coef_vector),
    alpha_per_fold_nested = list(alpha_per_fold),
    
    # Reproducibility
    seed = seed,
    date_run = as.character(Sys.time()),
    glmnet_version = as.character(packageVersion("glmnet")),
    
    # Store the fitted model object (optional, for later use)
    cv_fit_object = list(best_cv_fit)
  )
  
  return(result_row)
}


# ============================================
# 6. MODEL RUNNER FUNCTIONS
# ============================================
# ---------- 6.1 Run score specific models ----------
run_score_specific_models <- function(
    score_specific_datasets,
    alpha_grid = seq(0, 1, by = 0.1),
    nfolds = 10,
    seed = 123,
    output_file = "elastic_net_score_specific_results.rds"
) {
  results_list <- list()
  n_models <- length(score_specific_datasets)
  
  cat("==========================================================\n")
  cat("ELASTIC NET: SCORE-SPECIFIC DATASETS\n")
  cat(sprintf("Total models to run: %d\n", n_models))
  cat("==========================================================\n")
  
  for (i in seq_along(score_specific_datasets)) {
    cvd_score <- names(score_specific_datasets)[i]
    df <- score_specific_datasets[[i]]
    dataset_name <- paste0(cvd_score, "_specific")
    
    cat(sprintf("\n[Model %d/%d] %s ~ %s\n", i, n_models, cvd_score, dataset_name))
    cat("----------------------------------------------------------\n")
    
    # Extract outcome (y) and remove from predictors (X)
    if (!cvd_score %in% colnames(df)) {
      cat(sprintf("  ✗ ERROR: Column '%s' not found in dataset!\n", cvd_score))
      next
    }
    
    y <- df[[cvd_score]]
    X <- df %>% select(-all_of(cvd_score), -Sample_ID) %>% as.matrix()
    
    tryCatch({
      result <- run_elastic_net_model(
        X = X,
        y = y,
        cvd_score_name = cvd_score,
        dataset_name = dataset_name,
        imputation_method = "mean_mode",
        alpha_grid = alpha_grid,
        nfolds = nfolds,
        seed = seed
      )
      results_list[[i]] <- result
      cat("  ✓ Model completed successfully!\n")
    }, error = function(e) {
      cat(sprintf("  ✗ ERROR: %s\n", e$message))
      results_list[[i]] <- NULL
    })
  }
  
  results_df <- bind_rows(results_list)   # Combine results
  saveRDS(results_df, file = output_file)
  cat(sprintf("✓ Complete! %d/%d models successful\n", nrow(results_df), n_models))
  cat(sprintf("✓ Results saved to: %s\n", output_file))
  return(results_df)
}


# ---------- 6.2 Run created data set models ----------
run_common_datasets_models <- function(
    datasets_list,             # Named list of datasets (each has Sample_ID + all CVD scores)
    cvd_score_names,           # Vector of CVD score column names
    alpha_grid = seq(0, 1, by = 0.1),
    nfolds = 10,
    seed = 42,
    output_file = "elastic_net_common_datasets_results.rds"
) {
  results_list <- list()
  counter <- 1
  n_models <- length(datasets_list) * length(cvd_score_names)
  
  cat("==========================================================\n")
  cat("ELASTIC NET: COMMON DATASETS × CVD SCORES\n")
  cat(sprintf("Total models to run: %d\n", n_models))
  cat(sprintf("  %d datasets × %d CVD scores\n", length(datasets_list), length(cvd_score_names)))
  cat("==========================================================\n")
  
  for (dataset_name in names(datasets_list)) {
    df <- datasets_list[[dataset_name]]
    
    for (cvd_score in cvd_score_names) {
      
      cat(sprintf("\n[Model %d/%d] %s ~ %s\n", counter, n_models, cvd_score, dataset_name))
      cat("----------------------------------------------------------\n")
      
      # Check CVD score exists
      if (!cvd_score %in% colnames(df)) {
        cat(sprintf("  ✗ ERROR: Column '%s' not found in dataset!\n", cvd_score))
        counter <- counter + 1
        next
      }
      
      # Filter to rows with non-missing CVD score
      df_complete <- df %>% filter(!is.na(.data[[cvd_score]]))
      
      cat(sprintf("  Rows with valid %s: %d (from %d total)\n", 
                  cvd_score, nrow(df_complete), nrow(df)))
      
      # Extract y and X
      y <- df_complete[[cvd_score]]
      X <- df_complete %>% 
        select(-Sample_ID, -all_of(cvd_score_names)) %>%  # Remove Sample_ID and ALL CVD scores
        as.matrix()
      
      cat(sprintf("  Predictors: %d columns\n", ncol(X)))
      
      tryCatch({
        result <- run_elastic_net_model(
          X = X,
          y = y,
          cvd_score_name = cvd_score,
          dataset_name = dataset_name,
          imputation_method = "mean_mode",
          alpha_grid = alpha_grid,
          nfolds = nfolds,
          seed = seed
        )
        
        results_list[[counter]] <- result
        cat("  ✓ Model completed successfully!\n")
        
      }, error = function(e) {
        cat(sprintf("  ✗ ERROR: %s\n", e$message))
        results_list[[counter]] <- NULL
      })
      
      counter <- counter + 1
    }
  }
  
  results_df <- bind_rows(results_list) # Combine results
  saveRDS(results_df, file = output_file)
  cat(sprintf("✓ Complete! %d/%d models successful\n", nrow(results_df), n_models))
  cat(sprintf("✓ Results saved to: %s\n", output_file))
  return(results_df)
}

# ---------- 6.3 Run score specific + 1 further dataset models ----------
run_score_specific_plus_blocks_models <- function(
    score_specific_datasets,     # named list: each includes Sample_ID + its outcome + score-specific predictors
    datasets_list,               # named list: each includes Sample_ID + ALL outcomes + predictors
    alpha_grid = seq(0, 1, by = 0.1),
    nfolds = 10,
    seed = 42,
    output_file = "elastic_net_score_specific_plus_blocks_results.rds"
) {
  results_list <- list()
  counter <- 1
  n_models <- length(score_specific_datasets) * length(datasets_list)
  
  cat("==========================================================\n")
  cat("ELASTIC NET: SCORE-SPECIFIC + ONE COMMON BLOCK\n")
  cat(sprintf("Total models to run: %d\n", n_models))
  cat(sprintf("  %d scores × %d blocks\n", length(score_specific_datasets), length(datasets_list)))
  cat("==========================================================\n")
  
  for (cvd_score in names(score_specific_datasets)) {
    
    base_df <- score_specific_datasets[[cvd_score]]
    
    for (block_name in names(datasets_list)) {
      
      dataset_name <- paste0(cvd_score, "_specific+", block_name)
      
      cat(sprintf("\n[Model %d/%d] %s ~ %s\n", counter, n_models, cvd_score, dataset_name))
      cat("----------------------------------------------------------\n")
      
      # Build combined df (base + predictors from block)
      df_combined <- build_score_plus_block(
        score = cvd_score,
        block_name = block_name,
        score_specific_datasets = score_specific_datasets,
        datasets_list = datasets_list
      )
      
      
      # Need outcome column present
      if (!cvd_score %in% colnames(df_combined)) {
        cat(sprintf("  ✗ ERROR: Outcome '%s' not found after join!\n", cvd_score))
        results_list[[counter]] <- NULL
        counter <- counter + 1
        next
      }
      
      # Filter to valid outcomes
      df_complete <- df_combined %>% dplyr::filter(!is.na(.data[[cvd_score]]))
      
      cat(sprintf("  Rows with valid %s: %d (from %d total)\n",
                  cvd_score, nrow(df_complete), nrow(df_combined)))
      
      # Build X/y (IMPORTANT: here we remove ONLY the current outcome + Sample_ID
      # because the joined df should not contain other outcomes anymore)
      y <- df_complete[[cvd_score]]
      X <- df_complete %>%
        dplyr::select(-Sample_ID, -dplyr::all_of(cvd_score)) %>%
        as.matrix()
      
      cat(sprintf("  Predictors: %d columns\n", ncol(X)))
      
      tryCatch({
        result <- run_elastic_net_model(
          X = X,
          y = y,
          cvd_score_name = cvd_score,
          dataset_name = dataset_name,
          imputation_method = "mean_mode",
          alpha_grid = alpha_grid,
          nfolds = nfolds,
          seed = seed
        )
        
        results_list[[counter]] <- result
        cat("  ✓ Model completed successfully!\n")
        
      }, error = function(e) {
        cat(sprintf("  ✗ ERROR: %s\n", e$message))
        results_list[[counter]] <- NULL
      })
      
      counter <- counter + 1
    }
  }
  
  results_df <- dplyr::bind_rows(results_list)
  saveRDS(results_df, file = output_file)
  cat(sprintf("✓ Complete! %d/%d models successful\n", nrow(results_df), n_models))
  cat(sprintf("✓ Results saved to: %s\n", output_file))
  return(results_df)
}


# ============================================
# 7. RESULTS VIEWER
# ============================================
view_results_summary <- function(results_df) {
  cat("\n=== RESULTS SUMMARY ===\n\n")
  
  summary_table <- results_df %>%
    select(cvd_score, dataset_name, n_observations, n_predictors, 
           n_nonzero_coefs, alpha_optimal_fulldata, oof_r2, oof_rmse) %>%
    mutate(across(where(is.numeric), ~round(., 3)))
  
  print(summary_table, n = Inf)
  
  cat("\n=== TOP PREDICTORS PER MODEL ===\n\n")
  for (i in 1:nrow(results_df)) {
    cat(sprintf("%s:\n", results_df$cvd_score[i]))
    top_preds <- results_df$top_10_predictors[[i]]
    top_coefs <- results_df$top_10_coefficients[[i]]
    if (!is.na(top_preds[1])) {
      for (j in seq_along(top_preds)) {
        cat(sprintf("  %d. %s (coef = %.4f)\n", j, top_preds[j], top_coefs[j]))
      }
    } else {
      cat("  No predictors selected\n")
    }
    cat("\n")
  }
}

# ============================================
# 8. EXECUTE PIPELINE
# ============================================

results_score_specific <- run_score_specific_models(
  score_specific_datasets = score_specific_datasets,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 10,
  seed = 42,
  output_file = "elastic_net_score_specific_results.rds"
)


# RUN THE 20 MODELS
results_common <- run_common_datasets_models(
  datasets_list = datasets_list,
  cvd_score_names = CVD_scores,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 10,
  seed = 42,
  output_file = "elastic_net_common_datasets_results.rds"
)

# RUN SCORE-SPECIFIC + ONE BLOCK (5 × 7 = 35 MODELS)
results_score_plus_blocks <- run_score_specific_plus_blocks_models(
  score_specific_datasets = score_specific_datasets,
  datasets_list = datasets_list,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 10,
  seed = 42,
  output_file = "elastic_net_score_specific_plus_blocks_results.rds"
)

# Combine all results
results_all <- bind_rows(results_score_specific, results_common, results_score_plus_blocks)

# View summary
View(results_all)
view_results_summary(results_all)

# ============================================
# 9. EXPORT TO EXCEL
# ============================================

# ---------- Sheet 1: Main Results ----------
results_main <- results_all %>%
  select(
    cvd_score, dataset_name, model_name,
    n_observations, n_predictors, n_nonzero_coefs,
    alpha_optimal_fulldata, alpha_mean_nested, lambda_min, lambda_1se,
    in_sample_r2, in_sample_rmse, in_sample_mae, in_sample_cor,
    oof_r2, oof_rmse, oof_mae, oof_cor,
    cv_mse_min, cv_mse_1se, pct_deviance_explained,
    imputation_method, n_folds, seed, date_run
  )


# ---------- Sheet 2: Top Predictors per Model ----------
predictor_rows <- list()
for (i in 1:nrow(results_all)) {
  cvd <- results_all$cvd_score[i]
  dataset <- results_all$dataset_name[i]
  preds <- results_all$top_10_predictors[[i]]
  coefs <- results_all$top_10_coefficients[[i]]
  
  # Check if predictors exist and are not NA
  if (!is.null(preds) && length(preds) > 0 && !all(is.na(preds))) {
    for (j in seq_along(preds)) {
      if (!is.na(preds[j])) {
        predictor_rows[[length(predictor_rows) + 1]] <- tibble(
          cvd_score = cvd,
          dataset_name = dataset,
          rank = j,
          predictor_name = preds[j],
          coefficient = round(coefs[j], 4)
        )
      }
    }
  }
}

top_predictors_sheet <- bind_rows(predictor_rows) # Combine all rows


# ---------- Sheet 3: OOF R² overview tables
# ---- Table A: ALL MODELS oof_r2 (rows = dataset_name, cols = CVD scores) ----
table_A_all_models <- results_all %>%
  select(dataset_name, cvd_score, oof_r2) %>%
  mutate(oof_r2 = round(oof_r2, 3)) %>%
  distinct(dataset_name, cvd_score, .keep_all = TRUE) %>%  # safety
  pivot_wider(names_from = cvd_score, values_from = oof_r2) %>%
  mutate(
    dataset_group = case_when(
      grepl("_specific\\+", dataset_name) ~ "specific+block",
      grepl("_specific$", dataset_name)   ~ "specific",
      TRUE                                ~ "block_only"
    )
  ) %>%
  arrange(
    factor(dataset_group, levels = c("specific", "block_only", "specific+block")),
    dataset_name
  ) %>%
  select(-dataset_group)

# ---- Table B: delta oof_r2 = (score-specific + block) − (score-specific baseline) ----
baseline_specific <- results_all %>%
  filter(dataset_name == paste0(cvd_score, "_specific")) %>%
  select(cvd_score, baseline_oof_r2 = oof_r2)

table_B_delta <- results_all %>%
  filter(grepl("_specific\\+", dataset_name)) %>%
  mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
  select(cvd_score, block, oof_r2) %>%
  left_join(baseline_specific, by = "cvd_score") %>%
  mutate(delta_oof_r2 = round(oof_r2 - baseline_oof_r2, 3)) %>%
  select(block, cvd_score, delta_oof_r2) %>%
  pivot_wider(names_from = cvd_score, values_from = delta_oof_r2) %>%
  arrange(match(block, names(datasets_list)))


# ---------- Create Excel Workbook ----------
headerStyle <- createStyle(textDecoration = "bold", fgFill = "#D3D3D3")

wb <- createWorkbook()

addWorksheet(wb, "All_Results")
writeData(wb, "All_Results", results_main)
addStyle(wb, "All_Results", headerStyle, rows = 1, cols = 1:ncol(results_main), gridExpand = TRUE)
freezePane(wb, "All_Results", firstRow = TRUE)
setColWidths(wb, "All_Results", cols = 1:ncol(results_main), widths = "auto")

addWorksheet(wb, "Top_Predictors")
writeData(wb, "Top_Predictors", top_predictors_sheet)
addStyle(wb, "Top_Predictors", headerStyle, rows = 1, cols = 1:ncol(top_predictors_sheet), gridExpand = TRUE)
freezePane(wb, "Top_Predictors", firstRow = TRUE)
setColWidths(wb, "Top_Predictors", cols = 1:ncol(top_predictors_sheet), widths = "auto")

addWorksheet(wb, "OOF_R2_Overview")
writeData(wb, "OOF_R2_Overview", "Table A: OOF R² (all models)", startRow = 1, startCol = 1)
writeData(wb, "OOF_R2_Overview", table_A_all_models, startRow = 2, startCol = 1)

start_row_B <- nrow(table_A_all_models) + 5
writeData(wb, "OOF_R2_Overview", 
                       "Table B: Δ OOF R² = (score-specific + block) − (score-specific baseline)",
                       startRow = start_row_B, startCol = 1)
writeData(wb, "OOF_R2_Overview", table_B_delta, startRow = start_row_B + 1, startCol = 1)

addStyle(wb, "OOF_R2_Overview", headerStyle,
                     rows = c(2, start_row_B + 1),
                     cols = 1:max(ncol(table_A_all_models), ncol(table_B_delta)),
                     gridExpand = TRUE)
setColWidths(wb, "OOF_R2_Overview",
                             cols = 1:max(ncol(table_A_all_models), ncol(table_B_delta)),
                             widths = "auto")

saveWorkbook(wb, "elastic_net_results.xlsx", overwrite = TRUE)

cat("\n✓ Excel file saved: elastic_net_results.xlsx\n")
cat("  - Sheet 1: All_Results (", nrow(results_main), " models)\n")
cat("  - Sheet 2: Top_Predictors (", nrow(top_predictors_sheet), " rows)\n")
cat("  - Sheet 3: OOF_R2_Overview\n")