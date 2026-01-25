### Elastic Net Analysis - CoDiet CVD Prediction
### Author: Luisa Delius

### Table of Content
# Data preparation: Step 1 and 2

# Load packages
# Load required packages
library(glmnet)      # Elastic net regression
library(caret)       # Cross-validation framework
library(dplyr)       # Data manipulation
library(tidyr)       # Data reshaping
library(ggplot2)     # Visualization
library(missForest)  # Imputation
library(parallel)    # Parallel processing
library(doParallel)

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

CVD_scores <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y", "mean_risk")

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
  mutate(townsend = replace_na(townsend, 0)) %>%
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
# 3.1 Missing Data Imputation using Mean of each column (numeric) or mode (factor)
# ============================================
# Function for mean imputation of z_ columns only
impute_mean_z <- function(df) {
  df %>%
    mutate(across(starts_with("z_"), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))
}

# Function for mode imputation of factor columns (excluding Country)
impute_mode_factors <- function(df) {
  df %>%
    mutate(across(
      where(~(is.factor(.) | is.character(.)) & !all(is.na(.))),
      ~{
        col_name <- cur_column()
        if(col_name %in% c("Country", "Sample_ID")) {
          .  # Don't impute
        } else {
          # Get mode value
          mode_val <- names(sort(table(., useNA="ifany"), decreasing = TRUE))[1]
          # Preserve factor/character type
          x_imputed <- as.character(.)
          x_imputed[is.na(x_imputed)] <- mode_val
          # Convert back to factor if original was factor
          if(is.factor(.)) {
            factor(x_imputed, levels = levels(.))
          } else {
            x_imputed
          }
        }
      }
    ))
}

# Apply imputation to all models
df_model1_imputed <- impute_mean_z(df_model1)

df_model2_imputed <- df_model2 %>%
  impute_mean_z() %>%
  impute_mode_factors()

df_model3_imputed <- df_model3 %>%
  impute_mean_z() %>%
  impute_mode_factors()

df_model4_imputed <- df_model4 %>%
  impute_mean_z() %>%
  impute_mode_factors()

df_body_comp_imputed <- impute_mean_z(df_body_composition_ElaNet) %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")

df_demographics_imputed <- df_REDcap_demographics_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID") %>%
  impute_mean_z() %>%
  impute_mode_factors()

df_risk_factors_imputed <-  df_risk_factors_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID") %>%
  select(where(~mean(is.na(.)) <= 0.5)) %>% # excluded columns with more then 50% missingness (could be even more strict I guess...)
  impute_mean_z() %>%
  impute_mode_factors()

datasets_list <-  list(all_data = df_model2_imputed,
                       metabolites = df_model1_imputed,
                       mtc_sign = df_model3_imputed,
                       sign_wo_mtc = df_model4_imputed,
                       body_composition = df_body_comp_imputed,
                       sociodemographics_lifestyle = df_demographics_imputed,
                       clinical_risk_factors = df_risk_factors_imputed)

# ============================================
## 3.2 Missing Data Imputation - missForest

# # --- STEP 1: Set up parallel processing ---
# n_cores <- detectCores() - 1  # Use all cores except 1 (leave one for system)
# cl <- makeCluster(n_cores)
# registerDoParallel(cl)
# cat("Using", n_cores, "cores\n")
# 
# # --- STEP 2: Prepare data for imputation ---
# # Extract only predictor columns (no Sample_ID, Country, or CVD outcomes)
# df_to_impute <- df_model1 %>%
#   select(-Sample_ID, -Country, -QRISK3_risk, -SCORE2_score, -ascvd_10y, -frs_10y, -mean_risk)  # predictor_cols_m1 was defined earlier
# 
# # --- STEP 3: Run missForest imputation ---
# imputation_result <- missForest(
#   df_to_impute,        # Data to impute
#   maxiter = 10,              # Maximum iterations (usually converges before this)
#   ntree = 500,                # Number of trees per random forest
#   parallelize = "variables",    # Use parallel processing across trees
#   verbose = TRUE              # Show progress during imputation
# )
# 
# # --- STEP 4: Stop parallel processing ---
# stopCluster(cl)
# 
# # --- STEP 5: Extract and verify imputed data ---
# df_imputed <- imputation_result$ximp
# 
# # Check imputation quality using Out-of-Bag (OOB) error
# # Lower = better (want <0.3 for excellent, <0.5 for acceptable)
# cat("Imputation Quality (NRMSE - Normalized Root Mean Square Error):\n")
# cat("  OOB error:", round(imputation_result$OOBerror, 4), "\n")


# ============================================
## 4. ELASTIC NET REGRESSION PIPELINE FOR CVD RISK PREDICTION
# ============================================
# ----------------------------------------------------------------------------
# 4.1 FUNCTION: Run single elastic net model with alpha grid search
# ----------------------------------------------------------------------------

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
  
  # Grid search over alpha values
  cv_results <- lapply(alpha_grid, function(a) {
    cv_fit <- cv.glmnet(
      x = X_clean,
      y = y_clean,
      alpha = a,
      nfolds = nfolds,
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
  best_alpha <- alpha_grid[best_alpha_idx]
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
  
  # Get predictions on full data
  y_pred <- predict(best_cv_fit, newx = X_clean, s = "lambda.min")[,1]
  
  # Calculate performance metrics
  residuals <- y_clean - y_pred
  mse <- mean(residuals^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(residuals))
  
  # R-squared
  ss_total <- sum((y_clean - mean(y_clean))^2)
  ss_residual <- sum(residuals^2)
  r2 <- 1 - (ss_residual / ss_total)
  
  # Correlation between predicted and observed
  cor_pred_obs <- cor(y_pred, y_clean)
  
  # Deviance explained (from glmnet)
  dev_explained <- best_cv_fit$glmnet.fit$dev.ratio[which(best_cv_fit$glmnet.fit$lambda == lambda_min)]
  
  # Extract top predictors (by absolute coefficient value)
  nonzero_coefs <- coef_vector[coef_vector != 0]
  if (length(nonzero_coefs) > 0) {
    sorted_coefs <- sort(abs(nonzero_coefs), decreasing = TRUE)
    top_n <- min(5, length(sorted_coefs))
    top_predictor_names <- names(sorted_coefs)[1:top_n]
    top_predictor_coefs <- nonzero_coefs[top_predictor_names]
  } else {
    top_predictor_names <- NA
    top_predictor_coefs <- NA
  }
  
  # Compile results
  result_row <- tibble(
    # Model identifiers
    model_name = "elastic_net",
    cvd_score = cvd_score_name,
    dataset_name = dataset_name,
    
    # Sample characteristics
    n_observations = n_obs,
    n_predictors = n_pred,
    n_nonzero_coefs = n_nonzero,
    
    # Hyperparameters
    alpha_optimal = best_alpha,
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
    full_data_mse = mse,
    full_data_rmse = rmse,
    full_data_mae = mae,
    full_data_r2 = r2,
    cor_pred_obs = cor_pred_obs,
    pct_deviance_explained = dev_explained * 100,
    
    # Model characteristics
    intercept = intercept,
    top_5_predictors = list(top_predictor_names),
    top_5_coefficients = list(top_predictor_coefs),
    all_predictor_names = list(predictor_names),
    all_coefficients = list(coef_vector),
    
    # Reproducibility
    seed = seed,
    date_run = as.character(Sys.time()),
    glmnet_version = as.character(packageVersion("glmnet")),
    
    # Store the fitted model object (optional, for later use)
    cv_fit_object = list(best_cv_fit)
  )
  
  return(result_row)
}


# ----------------------------------------------------------------------------
# FUNCTION 4.2: Run score-specific models
# ----------------------------------------------------------------------------

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
  cat("==========================================================\n")
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
  
  # Combine results
  results_df <- bind_rows(results_list)
  
  # Save
  saveRDS(results_df, file = output_file)
  
  cat("\n==========================================================\n")
  cat(sprintf("✓ Complete! %d/%d models successful\n", nrow(results_df), n_models))
  cat(sprintf("✓ Results saved to: %s\n", output_file))
  cat("==========================================================\n\n")
  
  return(results_df)
}


# ----------------------------------------------------------------------------
# FUNCTION 4.3: Run common datasets across all CVD scores
# ----------------------------------------------------------------------------

run_common_datasets_models <- function(
    datasets_list,              # Named list of datasets (each has Sample_ID + all CVD scores)
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
  cat("==========================================================\n")
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
  
  # Combine results
  results_df <- bind_rows(results_list)
  
  # Save
  saveRDS(results_df, file = output_file)
  
  cat("\n==========================================================\n")
  cat(sprintf("✓ Complete! %d/%d models successful\n", nrow(results_df), n_models))
  cat(sprintf("✓ Results saved to: %s\n", output_file))
  cat("==========================================================\n\n")
  
  return(results_df)
}


# ----------------------------------------------------------------------------
# FUNCTION 3: View results summary
# ----------------------------------------------------------------------------

view_results_summary <- function(results_df) {
  
  cat("\n=== RESULTS SUMMARY ===\n\n")
  
  summary_table <- results_df %>%
    select(cvd_score, dataset_name, n_observations, n_predictors, 
           n_nonzero_coefs, alpha_optimal, full_data_r2, full_data_rmse) %>%
    mutate(across(where(is.numeric), ~round(., 3)))
  
  print(summary_table, n = Inf)
  
  cat("\n=== TOP PREDICTORS PER MODEL ===\n\n")
  for (i in 1:nrow(results_df)) {
    cat(sprintf("%s:\n", results_df$cvd_score[i]))
    top_preds <- results_df$top_5_predictors[[i]]
    top_coefs <- results_df$top_5_coefficients[[i]]
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

# ----------------------------------------------------------------------------
### Execute the pipeline and view results
# ----------------------------------------------------------------------------
# RUN THIS CODE:
results_score_specific <- run_score_specific_models(
  score_specific_datasets = score_specific_datasets,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 10,
  seed = 42,
  output_file = "elastic_net_score_specific_results.rds"
)

View(results_score_specific)
view_results_summary(results_score_specific)



# RUN THE 20 MODELS
results_common <- run_common_datasets_models(
  datasets_list = datasets_list,
  cvd_score_names = CVD_scores,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 10,
  seed = 42,
  output_file = "elastic_net_common_datasets_results.rds"
)

# COMBINE WITH YOUR 5 SCORE-SPECIFIC RESULTS
results_all <- bind_rows(results_score_specific, results_common)

# VIEW
View(results_all)
view_results_summary(results_all)

# ============================================
## save as excel
# ============================================
library(openxlsx)
# ============================================
# Prepare Sheet 1: Main Results
# ============================================
results_main <- results_all %>%
  select(
    # Model identifiers
    cvd_score, dataset_name, model_name,
    # Sample info
    n_observations, n_predictors, n_nonzero_coefs,
    # Hyperparameters
    alpha_optimal, lambda_min, lambda_1se,
    # Performance metrics
    full_data_r2, full_data_rmse, full_data_mae, 
    cv_mse_min, cv_mse_1se, pct_deviance_explained, cor_pred_obs,
    # Metadata
    imputation_method, n_folds, seed, date_run
  )

# ============================================
# Prepare Sheet 2: Top Predictors per Model (ROBUST VERSION)
# ============================================
# Build predictor list row by row
predictor_rows <- list()

for (i in 1:nrow(results_all)) {
  cvd <- results_all$cvd_score[i]
  dataset <- results_all$dataset_name[i]
  
  preds <- results_all$top_5_predictors[[i]]
  coefs <- results_all$top_5_coefficients[[i]]
  
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

# Combine all rows
top_predictors_sheet <- bind_rows(predictor_rows)

# ============================================
# Create Excel Workbook
# ============================================
wb <- createWorkbook()

# Add Sheet 1: Main Results
addWorksheet(wb, "All_Results")
writeData(wb, "All_Results", results_main)

# Add Sheet 2: Top Predictors
addWorksheet(wb, "Top_Predictors")
writeData(wb, "Top_Predictors", top_predictors_sheet)

# Optional: Add formatting
headerStyle <- createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
addStyle(wb, "All_Results", headerStyle, rows = 1, cols = 1:ncol(results_main), gridExpand = TRUE)
addStyle(wb, "Top_Predictors", headerStyle, rows = 1, cols = 1:ncol(top_predictors_sheet), gridExpand = TRUE)

# Freeze first row
freezePane(wb, "All_Results", firstRow = TRUE)
freezePane(wb, "Top_Predictors", firstRow = TRUE)

# Auto-adjust column widths
setColWidths(wb, "All_Results", cols = 1:ncol(results_main), widths = "auto")
setColWidths(wb, "Top_Predictors", cols = 1:ncol(top_predictors_sheet), widths = "auto")

# Save
saveWorkbook(wb, "elastic_net_results.xlsx", overwrite = TRUE)

cat("\n✓ Excel file saved: elastic_net_results.xlsx\n")
cat("  - Sheet 1: All_Results (main summary with", nrow(results_main), "models)\n")
cat("  - Sheet 2: Top_Predictors (", nrow(top_predictors_sheet), "predictor rows)\n")
