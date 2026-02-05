### Elastic Net Regression with Preprocessing Inside Each CV-Fold (No Data Leakage)
### Author: Luisa Delius
 
### KEY DIFFERENCE FROM PREVIOUS VERSION:
### - Scaling (z-scoring) happens INSIDE each CV fold
### - Imputation happens INSIDE each CV fold  
### - Dummy coding happens INSIDE each CV fold
### - Parameters computed from TRAINING data only, then applied to test data
### - This eliminates data leakage and provides honest performance estimates

# ============================================
# SETUP
# ============================================
library(glmnet)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(parallel)
library(flextable)
library(officer)
library(patchwork)

set.seed(42)

# Set working directory
wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files/processed_data"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

# Parallel processing setup
CPUS <- parallel::detectCores() - 3 # change to -1 if running over night
cat(sprintf("Using %d CPU cores for parallel processing\n", CPUS))

# ============================================
# 1. DATA LOADING (RAW, UNSCALED DATA)
# ============================================
df_fatty_acids_ElaNet <- readRDS("df_fatty_acids_predictor_statin_suppl.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk)

df_lipidomics_ElaNet <- readRDS("df_lipidomics_predictor_statin_suppl.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk, -Statins, -Supplements) %>%
  filter(!is.na(Sample_ID))

df_urine_nmr_ElaNet <- readRDS("df_urine_NMR_data.rds") %>%
  select(-starts_with("z_"))

df_body_composition_ElaNet <- readRDS("df_body_composition_metrics.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk) %>%
  filter(!is.na(Sample_ID))

df_REDcap_demographics_ElaNet <- readRDS("df_REDcap_demographics_ElaNet.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk, -recruitment_site)

df_risk_factors_ElaNet <- readRDS("df_risk_factor_predictors.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk,
         -Heart.Rate, -Age, -Total.Cholesterol.mg.dl, -HDL.mg.dl,
         -Body.Weight, -Height, -BMI, -Gender,
         -stress_resilience_status, -stress_index_status, -Age.Risk)

df_all_cvd_risk_scores_ElaNet <- readRDS("df_all_risk_scores.rds") %>%
  rename(Sample_ID = PatientID) %>%
  select(-SCORE2_strat)

df_ascvd_frs_score2_input <- readRDS("ASCVD_SCORE2_Framingham_input.rds") %>%
  rename(Sample_ID = PatientID)

df_QRISK3_input <- readRDS("QRISK3_calculation_input.rds") %>%
  rename(Sample_ID = PatientID)

CVD_scores <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y", "mean_risk")

# ============================================
# 2. CREATE DATASETS
# ============================================
# Model 1: Metabolites only
df_model1 <- df_fatty_acids_ElaNet %>%
  select(-Statins, -Supplements) %>%
  full_join(df_lipidomics_ElaNet, by = "Sample_ID") %>%
  full_join(df_urine_nmr_ElaNet, by = "Sample_ID") %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID") %>%
  select(where(~!all(is.na(.))))

# Model 2: All predictors
df_model2 <- df_model1 %>%
  full_join(df_fatty_acids_ElaNet %>% select(Sample_ID, Statins, Supplements), by = "Sample_ID") %>%
  full_join(df_body_composition_ElaNet, by = "Sample_ID") %>%
  full_join(df_risk_factors_ElaNet, by = "Sample_ID") %>%
  full_join(df_REDcap_demographics_ElaNet, by = "Sample_ID")

# Model 3: BH-significant predictors
BH_sig_predictors <- c("ecw_tbw", "mean_hrt", "ecm_bcm", "rbc_22_4n_6",
                       "naps_during_day", "Living_Status")

df_model3 <- df_model2 %>%
  select(Sample_ID, all_of(BH_sig_predictors), all_of(CVD_scores))

# Model 4: BH-significant + nominally significant
sig_predictors <- c("AGE.reader", "Trigonelline", "Hba1C", "ALT.unit.L",
                    "rbc_dpa_22_5n3", "rbc_eicosadienoic_20_2n6",
                    "rbc_epa_20_5n3", "pe_o_19_1_20_5",
                    "3_hydroxybutyric_acid", "hippuric_acid",
                    "Employment_Status", "Education_Level")

df_model4 <- df_model2 %>%
  select(Sample_ID, all_of(BH_sig_predictors), all_of(sig_predictors), all_of(CVD_scores))

# Model 0: Score-specific input datasets
df_model0_QRISK3 <- df_QRISK3_input %>%
  mutate(townsend = replace_na(townsend, 0)) %>%
  full_join(df_all_cvd_risk_scores_ElaNet %>% select(Sample_ID, QRISK3_risk), by = "Sample_ID") %>%
  filter(!is.na(QRISK3_risk))

df_model0_ascvd <- df_ascvd_frs_score2_input %>%
  select(-Risk.region, -mean_LDL_mg_dl) %>%
  full_join(df_all_cvd_risk_scores_ElaNet %>% select(Sample_ID, ascvd_10y), by = "Sample_ID") %>%
  filter(!is.na(ascvd_10y))

df_model0_score2 <- df_ascvd_frs_score2_input %>%
  select(-mean_LDL_mg_dl, -blood_pressure_treatment, -race_ascvd) %>%
  full_join(df_all_cvd_risk_scores_ElaNet %>% select(Sample_ID, SCORE2_score), by = "Sample_ID") %>%
  filter(!is.na(SCORE2_score))

df_model0_frs <- df_ascvd_frs_score2_input %>%
  select(-Risk.region, -mean_LDL_mg_dl, -race_ascvd) %>%
  full_join(df_all_cvd_risk_scores_ElaNet %>% select(Sample_ID, frs_10y), by = "Sample_ID") %>%
  filter(!is.na(frs_10y))

df_model0_composite <- df_QRISK3_input %>%
  mutate(townsend = replace_na(townsend, 0)) %>%
  full_join(df_ascvd_frs_score2_input %>% select(Sample_ID, Risk.region), by = "Sample_ID") %>%
  full_join(df_all_cvd_risk_scores_ElaNet %>% select(Sample_ID, mean_risk), by = "Sample_ID") %>%
  filter(!is.na(mean_risk))

score_specific_datasets <- list(
  QRISK3_risk = df_model0_QRISK3,
  SCORE2_score = df_model0_score2,
  frs_10y = df_model0_frs,
  ascvd_10y = df_model0_ascvd,
  mean_risk = df_model0_composite
)

# Predictor block datasets (with CVD scores for later filtering)
df_body_comp_raw <- df_body_composition_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")

df_demographics_raw <- df_REDcap_demographics_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")

df_risk_factors_raw <- df_risk_factors_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID") %>%
  select(where(~mean(is.na(.)) <= 0.5))  # Exclude columns with >50% missing

df_lipids_raw <- df_lipidomics_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")
  
df_fatty_Acids_raw <- df_fatty_acids_ElaNet %>%
  select(-Statins, -Supplements) %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")
  
df_urine_nmr_raw <- df_urine_nmr_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")

  
datasets_list <- list(
  all_data = df_model2,
  metabolites = df_model1,
  mtc_sign = df_model3,
  sign_wo_mtc = df_model4,
  body_composition = df_body_comp_raw,
  sociodemographics_lifestyle = df_demographics_raw,
  clinical_risk_factors = df_risk_factors_raw,
  lipids = df_lipids_raw,
  fatty_acids = df_fatty_Acids_raw,
  urine_nmr = df_urine_nmr_raw
)

# ============================================
# ADDITIONAL MODELS: Body Composition + Height/Weight
# ============================================

# Get height and weight from QRISK3 input (where they exist)
height_weight_vars <- df_QRISK3_input %>%
  select(Sample_ID, Height_cm, Weight_kg)

# Create extended body composition dataset
df_body_comp_extended <- df_body_composition_ElaNet %>%
  full_join(height_weight_vars, by = "Sample_ID") %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")

# Scores that DON'T already include height/weight
scores_without_ht_wt <- c("SCORE2_score", "ascvd_10y", "frs_10y")



# ============================================
# 2b. SANITIZE DATA TYPES
# ============================================
# Force all columns to numeric or factor to prevent preprocessing issues
sanitize_df <- function(df) {
    for (col in names(df)) {
        # Skip ID and outcome columns
        if (col == "Sample_ID") next
        if (col %in% CVD_scores) next
        
        x <- df[[col]]
        
        if (is.logical(x)) {
            # Convert logical to 0/1 integer
            df[[col]] <- as.integer(x)
          } else if (inherits(x, "Date") || inherits(x, "POSIXct") || inherits(x, "POSIXlt")) {
              # Convert dates to numeric (days since epoch)
              df[[col]] <- as.numeric(x)
            } else if (!is.numeric(x) && !is.factor(x)) {
                # Force character/other types to factor
                df[[col]] <- as.factor(x)
              }
      }
    df
  }

# Apply sanitization to all datasets
cat("Sanitizing data types...\n")
datasets_list <- lapply(datasets_list, sanitize_df)
score_specific_datasets <- lapply(score_specific_datasets, sanitize_df)
df_body_comp_extended <- sanitize_df(df_body_comp_extended)  # ADD THIS LINE
cat("✓ All datasets sanitized (logical → integer, character → factor)\n")


# ============================================
# 3. PREPROCESSING FUNCTIONS (FOLD-WISE)
# ============================================

# Compute preprocessing parameters from training data
compute_preprocess_params <- function(df_train, exclude_cols = NULL) {
  
  # Identify column types
  cols_to_process <- setdiff(names(df_train), exclude_cols)
  
  numeric_cols <- cols_to_process[sapply(df_train[cols_to_process], is.numeric)]
  factor_cols <- cols_to_process[sapply(df_train[cols_to_process], function(x) is.factor(x) | is.character(x))]
  logical_cols <- cols_to_process[sapply(df_train[cols_to_process], is.logical)]
  
  # Compute means and SDs for numeric columns
  means <- sapply(df_train[numeric_cols], mean, na.rm = TRUE)
  sds <- sapply(df_train[numeric_cols], sd, na.rm = TRUE)
  # Handle zero SD (constant columns) - set to 1 to avoid division by zero
  sds[sds == 0 | is.na(sds)] <- 1
  
  # Compute modes for imputation (numeric and factor)
  compute_mode <- function(x) {
    if (all(is.na(x))) return(NA)
    ux <- unique(na.omit(x))
    ux[which.max(tabulate(match(na.omit(x), ux)))]
  }
  
  numeric_modes <- sapply(df_train[numeric_cols], function(x) mean(x, na.rm = TRUE))  # Use mean for numeric
  factor_modes <- sapply(df_train[factor_cols], compute_mode)
  logical_modes <- sapply(df_train[logical_cols], function(x) {
    if (all(is.na(x))) return(FALSE)
    as.logical(compute_mode(as.integer(x)))
  })
  
  # Get factor levels from training data
  factor_levels <- lapply(df_train[factor_cols], function(x) {
    if (is.factor(x)) levels(x) else unique(na.omit(x))
  })
  
  list(
    numeric_cols = numeric_cols,
    factor_cols = factor_cols,
    logical_cols = logical_cols,
    means = means,
    sds = sds,
    numeric_modes = numeric_modes,
    factor_modes = factor_modes,
    logical_modes = logical_modes,
    factor_levels = factor_levels
  )
}


# Apply preprocessing parameters to data (train or test)
apply_preprocess <- function(df, params, exclude_cols = NULL) {
  
  df_processed <- df
  
  
  # --- 1. Impute and scale numeric columns ---
  for (col in params$numeric_cols) {
    if (col %in% names(df_processed)) {
      # Impute with training mean
      df_processed[[col]][is.na(df_processed[[col]])] <- params$numeric_modes[col]
      # Z-score using training mean and SD
      df_processed[[col]] <- (df_processed[[col]] - params$means[col]) / params$sds[col]
    }
  }
  
  # --- 2. Impute and process factor columns ---
  for (col in params$factor_cols) {
    if (col %in% names(df_processed)) {
      x <- df_processed[[col]]
      
      # Convert to character for processing
      x <- as.character(x)
      
      # Impute missing with training mode
      x[is.na(x)] <- as.character(params$factor_modes[col])
      
      # Handle unseen levels: replace with mode
      known_levels <- params$factor_levels[[col]]
      x[!x %in% known_levels] <- as.character(params$factor_modes[col])
      
      # Convert to factor with training levels
      df_processed[[col]] <- factor(x, levels = known_levels)
    }
  }
  
  # --- 3. Impute logical columns ---
  for (col in params$logical_cols) {
    if (col %in% names(df_processed)) {
      df_processed[[col]][is.na(df_processed[[col]])] <- params$logical_modes[col]
    }
  }
  
  # --- 4. Create model matrix (dummy coding) ---
  # Keep only processed columns
  cols_for_matrix <- c(params$numeric_cols, params$factor_cols, params$logical_cols)
  cols_for_matrix <- intersect(cols_for_matrix, names(df_processed))
  
  df_for_matrix <- df_processed[, cols_for_matrix, drop = FALSE]
  
  # Create formula for model.matrix
  # model.matrix automatically creates dummies for factors
  X <- model.matrix(~ . - 1, data = df_for_matrix)  # -1 removes intercept (glmnet adds its own)
  
  return(X)
}


# Preprocess training and test data together (ensures consistent columns)
preprocess_train_test <- function(df_train, df_test, exclude_cols = NULL) {
  
  # Compute parameters from training data only
  params <- compute_preprocess_params(df_train, exclude_cols)
  
  # Apply to both datasets
  X_train <- apply_preprocess(df_train, params, exclude_cols)
  X_test <- apply_preprocess(df_test, params, exclude_cols)
  
  # Ensure same columns in both matrices
  common_cols <- intersect(colnames(X_train), colnames(X_test))
  X_train <- X_train[, common_cols, drop = FALSE]
  X_test <- X_test[, common_cols, drop = FALSE]
  
  list(
    X_train = X_train,
    X_test = X_test,
    params = params,
    feature_names = common_cols
  )
}

# ============================================
# 4. PERMUTATION TEST FUNCTION
# ============================================

#' Compute Q²_Y for a single permutation (shuffled y)
#' Uses same preprocessing and CV structure as main analysis

compute_permutation_Q2 <- function(df_raw, y_shuffled, exclude_cols,
                                   alpha_grid = seq(0, 1, by = 0.1),
                                   nfolds = 10, seed = 42) {
  set.seed(seed)
  n <- length(y_shuffled)
  foldid <- sample(rep(1:nfolds, length.out = n))
  
  y_oof <- rep(NA_real_, n)
  
  for (k in 1:nfolds) {
    idx_test <- which(foldid == k)
    idx_train <- setdiff(1:n, idx_test)
    
    df_train <- df_raw[idx_train, , drop = FALSE]
    df_test <- df_raw[idx_test, , drop = FALSE]
    y_train <- y_shuffled[idx_train]
    
    preprocessed <- tryCatch({
      preprocess_train_test(df_train, df_test, exclude_cols)
    }, error = function(e) NULL)
    
    if (is.null(preprocessed)) return(NA_real_)
    
    X_train <- preprocessed$X_train
    X_test <- preprocessed$X_test
    
    # Quick alpha selection (use fewer folds for speed)
    best_cv_fit <- NULL
    best_mse <- Inf
    
    for (a in alpha_grid) {
      cv_fit <- tryCatch({
        cv.glmnet(x = X_train, y = y_train, alpha = a,
                  nfolds = min(5, nfolds - 1),  # Faster inner CV
                  type.measure = "mse", standardize = FALSE, family = "gaussian")
      }, error = function(e) NULL)
      
      if (!is.null(cv_fit) && min(cv_fit$cvm) < best_mse) {
        best_mse <- min(cv_fit$cvm)
        best_cv_fit <- cv_fit
      }
    }
    
    if (is.null(best_cv_fit)) return(NA_real_)
    
    y_oof[idx_test] <- as.numeric(
      predict(best_cv_fit, newx = X_test, s = "lambda.min")
    )
  }
  
  # Compute Q²_Y
  ss_res <- sum((y_shuffled - y_oof)^2)
  ss_tot <- sum((y_shuffled - mean(y_shuffled))^2)
  Q2_Y <- 1 - ss_res / ss_tot
  
  return(Q2_Y)
}


#' Run permutation test to compute empirical p-value

run_permutation_test <- function(df_raw, y, exclude_cols, observed_Q2,
                                 n_permutations = 100,
                                 alpha_grid = seq(0, 1, by = 0.1),
                                 nfolds = 10, seed = 42) {
  set.seed(seed)
  
  cat(sprintf("  Running permutation test (%d permutations) on %d cores...\n", 
              n_permutations, CPUS))
  
  # Pre-generate all shuffled y vectors for reproducibility
  set.seed(seed)
  shuffled_y_list <- lapply(1:n_permutations, function(p) sample(y))
  
  # Run permutations in parallel
  null_Q2_values <- parallel::mclapply(1:n_permutations, function(p) {
    compute_permutation_Q2(
      df_raw = df_raw,
      y_shuffled = shuffled_y_list[[p]],
      exclude_cols = exclude_cols,
      alpha_grid = alpha_grid,
      nfolds = nfolds,
      seed = seed + p
    )
  }, mc.cores = CPUS)
  
  # Convert list to numeric vector
  null_Q2_values <- unlist(null_Q2_values)
  
  # Remove NA values (failed permutations)
  null_Q2_values <- null_Q2_values[!is.na(null_Q2_values)]
  
  # Empirical p-value: proportion of permuted Q² >= observed Q²
  p_value <- (sum(null_Q2_values >= observed_Q2) + 1) / (length(null_Q2_values) + 1)
  
  cat(sprintf("  Permutation test complete: p = %.4f\n", p_value))
  
  list(
    p_value = p_value,
    observed_Q2 = observed_Q2,
    null_distribution = null_Q2_values,
    n_permutations_completed = length(null_Q2_values)
  )
}


# ============================================
# 5. HELPER FUNCTIONS
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
# 6. CORE ELASTIC NET FUNCTION (WITH FOLD-WISE PREPROCESSING)
# ============================================

#' Nested CV with preprocessing inside each fold (NO DATA LEAKAGE)
#' 
#' For each outer fold:
#'   1. Split into train/test
#'   2. Compute preprocessing params from TRAINING ONLY
#'   3. Apply preprocessing to both train and test
#'   4. Grid search alpha via inner CV on processed training data
#'   5. Predict on processed test data

get_oof_predictions_clean <- function(df_raw, y, exclude_cols,
                                      alpha_grid = seq(0, 1, by = 0.1),
                                      nfolds = 10, seed = 42) {
  set.seed(seed)
  n <- length(y)
  y_mean_global <- mean(y)
  
  foldid_outer <- sample(rep(1:nfolds, length.out = n))
  
  y_oof <- rep(NA_real_, n)
  alpha_selected_per_fold <- numeric(nfolds)
  n_features_per_fold <- numeric(nfolds)
  Q2_per_fold <- numeric(nfolds)
  MAE_per_fold <- numeric(nfolds)
  RMSE_per_fold <- numeric(nfolds)
  
  for (k in 1:nfolds) {
    idx_test <- which(foldid_outer == k)
    idx_train <- setdiff(1:n, idx_test)
    
    df_train <- df_raw[idx_train, , drop = FALSE]
    df_test <- df_raw[idx_test, , drop = FALSE]
    y_train <- y[idx_train]
    y_test <- y[idx_test]
    
    preprocessed <- preprocess_train_test(df_train, df_test, exclude_cols)
    X_train <- preprocessed$X_train
    X_test <- preprocessed$X_test
    
    n_features_per_fold[k] <- ncol(X_train)
    
    cv_results_fold <- lapply(alpha_grid, function(a) {
      cv_fit <- cv.glmnet(
        x = X_train, y = y_train, alpha = a,
        nfolds = max(3, nfolds - 1),
        type.measure = "mse", standardize = FALSE, family = "gaussian"
      )
      list(alpha = a, cv_mse = min(cv_fit$cvm), cv_fit = cv_fit)
    })
    
    cv_mse_vals <- sapply(cv_results_fold, function(x) x$cv_mse)
    best_idx <- which.min(cv_mse_vals)
    alpha_selected_per_fold[k] <- cv_results_fold[[best_idx]]$alpha
    best_cv_fit <- cv_results_fold[[best_idx]]$cv_fit
    
    y_pred_k <- as.numeric(predict(best_cv_fit, newx = X_test, s = "lambda.min"))
    y_oof[idx_test] <- y_pred_k
    
    # Per-fold metrics
    ss_res_k <- sum((y_test - y_pred_k)^2)
    ss_tot_k <- sum((y_test - y_mean_global)^2)
    Q2_per_fold[k] <- 1 - ss_res_k / ss_tot_k
    MAE_per_fold[k] <- mean(abs(y_test - y_pred_k))
    RMSE_per_fold[k] <- sqrt(mean((y_test - y_pred_k)^2))
  }
  
  list(
    predictions = y_oof,
    alpha_per_fold = alpha_selected_per_fold,
    n_features_per_fold = n_features_per_fold,
    foldid = foldid_outer,
    Q2_per_fold = Q2_per_fold,
    MAE_per_fold = MAE_per_fold,
    RMSE_per_fold = RMSE_per_fold
  )
}


#' Run elastic net model with fold-wise preprocessing

run_elastic_net_model <- function(df_raw, outcome_col, cvd_score_name, dataset_name,
                                  exclude_cols = "Sample_ID",
                                  alpha_grid = seq(0, 1, by = 0.1),
                                  nfolds = 10, 
                                  n_permutations = 100,
                                  seed = 42) {
  set.seed(seed)
  
  # Remove rows with missing outcome
  df_complete <- df_raw %>% filter(!is.na(.data[[outcome_col]]))
  y <- df_complete[[outcome_col]]
  
  # Columns to exclude from predictors
  all_exclude <- unique(c(exclude_cols, outcome_col, CVD_scores))
  
  n_obs <- nrow(df_complete)
  cat(sprintf("\n  Sample size: %d\n", n_obs))
  cat("  Testing alpha values:", paste(alpha_grid, collapse = ", "), "\n")
  
  # ---- PART A: Full-data model (for coefficients) ----
  # Preprocess full dataset
  params_full <- compute_preprocess_params(df_complete, all_exclude)
  X_full <- apply_preprocess(df_complete, params_full, all_exclude)
  n_pred <- ncol(X_full)
  predictor_names <- colnames(X_full)
  
  cat(sprintf("  Predictors after processing: %d\n", n_pred))
  
  # Grid search on full data (still has minor leakage for alpha selection, but
  # this is only for final coefficient estimation, not performance reporting)
  set.seed(seed)
  foldid_inner <- sample(rep(1:nfolds, length.out = n_obs))
  
  cv_results <- lapply(alpha_grid, function(a) {
    cv_fit <- cv.glmnet(
      x = X_full, y = y, alpha = a,
      nfolds = nfolds, foldid = foldid_inner,
      type.measure = "mse", standardize = FALSE, family = "gaussian"
    )
    min_cvm_idx <- which(cv_fit$lambda == cv_fit$lambda.min)
    list(alpha = a, cv_mse_min = cv_fit$cvm[min_cvm_idx], cv_fit = cv_fit)
  })
  
  cv_mse_values <- sapply(cv_results, function(x) x$cv_mse_min)
  best_alpha_idx <- which.min(cv_mse_values)
  best_alpha_fulldata <- alpha_grid[best_alpha_idx]
  best_cv_fit <- cv_results[[best_alpha_idx]]$cv_fit
  
  lambda_min <- best_cv_fit$lambda.min
  lambda_1se <- best_cv_fit$lambda.1se
  cv_mse_min <- min(best_cv_fit$cvm)
  cv_mse_1se <- best_cv_fit$cvm[which(best_cv_fit$lambda == lambda_1se)]
  
  final_coef <- coef(best_cv_fit, s = "lambda.min")
  intercept <- as.numeric(final_coef[1])
  coef_vector <- as.numeric(final_coef[-1])
  names(coef_vector) <- predictor_names
  n_nonzero <- sum(coef_vector != 0)
  
  # ---- PART B: In-sample predictions (R²_Y - optimistic) ----
  y_pred_in <- as.numeric(predict(best_cv_fit, newx = X_full, s = "lambda.min"))
  res_in <- y - y_pred_in
  mse_in <- mean(res_in^2)
  rmse_in <- sqrt(mse_in)
  mae_in <- mean(abs(res_in))
  R2_Y <- 1 - sum(res_in^2) / sum((y - mean(y))^2)
  cor_in <- cor(y_pred_in, y)
  
  # ---- PART C: Out-of-fold predictions (Q²_Y - honest, no leakage) ----
  oof_results <- get_oof_predictions_clean(
    df_raw = df_complete,
    y = y,
    exclude_cols = all_exclude,
    alpha_grid = alpha_grid,
    nfolds = nfolds,
    seed = seed + 100
  )
  
  y_pred_oof <- oof_results$predictions
  alpha_per_fold <- oof_results$alpha_per_fold
  Q2_per_fold <- oof_results$Q2_per_fold  # NEW
  
  res_oof <- y - y_pred_oof
  mse_oof <- mean(res_oof^2)
  rmse_oof <- sqrt(mse_oof)
  mae_oof <- mean(abs(res_oof))
  Q2_Y <- 1 - sum(res_oof^2) / sum((y - mean(y))^2)
  cor_oof <- cor(y_pred_oof, y)
  
  # NEW: Per-fold Q² statistics
  Q2_fold_mean <- mean(Q2_per_fold)
  Q2_fold_sd <- sd(Q2_per_fold)
  
  MAE_per_fold <- oof_results$MAE_per_fold
  RMSE_per_fold <- oof_results$RMSE_per_fold
  MAE_fold_mean <- mean(MAE_per_fold)
  MAE_fold_sd <- sd(MAE_per_fold)
  RMSE_fold_mean <- mean(RMSE_per_fold)
  RMSE_fold_sd <- sd(RMSE_per_fold)
  
  # ---- PART D: Model quality metrics ----
  Q2_R2_ratio <- Q2_Y / R2_Y
  R2_Q2_gap <- R2_Y - Q2_Y
  
  # Deviance explained
  dev_explained <- best_cv_fit$glmnet.fit$dev.ratio[
    which.min(abs(best_cv_fit$glmnet.fit$lambda - lambda_min))
  ]
  
  # ---- PART E: Permutation test (empirical p-value) ----
  if (n_permutations > 0) {
    perm_results <- run_permutation_test(
      df_raw = df_complete,
      y = y,
      exclude_cols = all_exclude,
      observed_Q2 = Q2_Y,
      n_permutations = n_permutations,
      alpha_grid = alpha_grid,
      nfolds = nfolds,
      seed = seed + 200
    )
    perm_p_value <- perm_results$p_value
    perm_n_completed <- perm_results$n_permutations_completed
    perm_null_mean <- mean(perm_results$null_distribution)
  } else {
    perm_p_value <- NA_real_
    perm_n_completed <- 0
    perm_null_mean <- NA_real_
    cat("  Permutation test skipped (n_permutations = 0)\n")
  }
  
  # Top predictors
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
  tibble(
    model_name = "elastic_net",
    cvd_score = cvd_score_name,
    dataset_name = dataset_name,
    n_observations = n_obs,
    n_predictors = n_pred,
    n_nonzero_coefs = n_nonzero,
    alpha_optimal_fulldata = best_alpha_fulldata,
    alpha_mean_nested = mean(alpha_per_fold),
    lambda_min = lambda_min,
    lambda_1se = lambda_1se,
    n_folds = nfolds,
    type_measure = "mse",
    family = "gaussian",
    standardize = FALSE,
    imputation_method = "fold_wise_mean_mode",
    preprocessing = "fold_wise_scaling_and_imputation",
    cv_mse_min = cv_mse_min,
    cv_mse_1se = cv_mse_1se,
    # In-sample metrics (R²_Y)
    R2_Y = R2_Y,
    in_sample_rmse = rmse_in,
    in_sample_mae = mae_in,
    in_sample_cor = cor_in,
    # Cross-validated metrics (Q²_Y)
    Q2_Y = Q2_Y,
    Q2_Y_rmse = rmse_oof,
    Q2_Y_mae = mae_oof,
    MAE_fold_mean = MAE_fold_mean,
    MAE_fold_sd = MAE_fold_sd,
    RMSE_fold_mean = RMSE_fold_mean,
    RMSE_fold_sd = RMSE_fold_sd,
    MAE_per_fold = list(MAE_per_fold),
    RMSE_per_fold = list(RMSE_per_fold),
    Q2_Y_cor = cor_oof,
    Q2_fold_mean = Q2_fold_mean,        # NEW
    Q2_fold_sd = Q2_fold_sd,            # NEW
    Q2_fold_summary = sprintf("%.3f ± %.3f", Q2_fold_mean, Q2_fold_sd),  # NEW: nice display format
    Q2_per_fold = list(Q2_per_fold),    # NEW: store all 10 values
    # Model quality indicators
    Q2_R2_ratio = Q2_R2_ratio,
    R2_Q2_gap = R2_Q2_gap,
    # Permutation test
    permutation_p_value = perm_p_value,
    permutation_n = perm_n_completed,
    permutation_null_mean_Q2 = perm_null_mean,
    pct_deviance_explained = dev_explained * 100,
    intercept = intercept,
    top_10_predictors = list(top_predictor_names),
    top_10_coefficients = list(top_predictor_coefs),
    all_predictor_names = list(predictor_names),
    all_coefficients = list(coef_vector),
    alpha_per_fold_nested = list(alpha_per_fold),
    n_features_per_fold = list(oof_results$n_features_per_fold),
    seed = seed,
    date_run = as.character(Sys.time()),
    glmnet_version = as.character(packageVersion("glmnet")),
    cv_fit_object = list(best_cv_fit)
  )
}


# ============================================
# 7. MODEL RUNNER FUNCTIONS
# ============================================

#' Run score-specific models
run_score_specific_models <- function(score_specific_datasets,
                                      alpha_grid = seq(0, 1, by = 0.1),
                                      nfolds = 10, 
                                      n_permutations = 100,
                                      seed = 42,
                                      output_file = "elastic_net_foldwise_preproc_score_specific.rds") {
  results_list <- list()
  n_models <- length(score_specific_datasets)
  
  cat("==========================================================\n")
  cat("ELASTIC NET (FOLD-WISE PREPROCESSING): SCORE-SPECIFIC DATASETS\n")
  cat(sprintf("Total models: %d\n", n_models))
  cat(sprintf("Permutations per model: %d\n", n_permutations))
  cat("==========================================================\n")
  
  for (i in seq_along(score_specific_datasets)) {
    cvd_score <- names(score_specific_datasets)[i]
    df <- score_specific_datasets[[i]]
    dataset_name <- paste0(cvd_score, "_specific")
    
    cat(sprintf("\n[%d/%d] %s ~ %s\n", i, n_models, cvd_score, dataset_name))
    cat("----------------------------------------------------------\n")
    
    tryCatch({
      results_list[[i]] <- run_elastic_net_model(
        df_raw = df,
        outcome_col = cvd_score,
        cvd_score_name = cvd_score,
        dataset_name = dataset_name,
        exclude_cols = "Sample_ID",
        alpha_grid = alpha_grid,
        nfolds = nfolds,
        n_permutations = n_permutations,
        seed = seed
      )
      cat("  ✓ Success\n")
    }, error = function(e) {
      cat(sprintf("  ✗ ERROR: %s\n", e$message))
      results_list[[i]] <- NULL
    })
  }
  
  results_df <- bind_rows(results_list)
  saveRDS(results_df, file = output_file)
  cat(sprintf("\n✓ Saved %d models to %s\n", nrow(results_df), output_file))
  results_df
}


#' Run common dataset models
run_common_datasets_models <- function(datasets_list, cvd_score_names,
                                       alpha_grid = seq(0, 1, by = 0.1),
                                       nfolds = 10, 
                                       n_permutations = 100,
                                       seed = 42,
                                       output_file = "elastic_net_foldwise_preproc_common.rds") {
  results_list <- list()
  counter <- 1
  n_models <- length(datasets_list) * length(cvd_score_names)
  
  cat("==========================================================\n")
  cat("ELASTIC NET (FOLD-WISE PREPROCESSING): COMMON DATASETS × CVD SCORES\n")
  cat(sprintf("Total models: %d (%d datasets × %d scores)\n",
              n_models, length(datasets_list), length(cvd_score_names)))
  cat(sprintf("Permutations per model: %d\n", n_permutations))
  cat("==========================================================\n")
  
  for (dataset_name in names(datasets_list)) {
    df <- datasets_list[[dataset_name]]
    
    for (cvd_score in cvd_score_names) {
      cat(sprintf("\n[%d/%d] %s ~ %s\n", counter, n_models, cvd_score, dataset_name))
      cat("----------------------------------------------------------\n")
      
      if (!cvd_score %in% colnames(df)) {
        cat(sprintf("  ✗ Column '%s' not found\n", cvd_score))
        counter <- counter + 1
        next
      }
      
      tryCatch({
        results_list[[counter]] <- run_elastic_net_model(
          df_raw = df,
          outcome_col = cvd_score,
          cvd_score_name = cvd_score,
          dataset_name = dataset_name,
          exclude_cols = "Sample_ID",
          alpha_grid = alpha_grid,
          nfolds = nfolds,
          n_permutations = n_permutations,
          seed = seed
        )
        cat("  ✓ Success\n")
      }, error = function(e) {
        cat(sprintf("  ✗ ERROR: %s\n", e$message))
        results_list[[counter]] <- NULL
      })
      counter <- counter + 1
    }
  }
  
  results_df <- bind_rows(results_list)
  saveRDS(results_df, file = output_file)
  cat(sprintf("\n✓ Saved %d models to %s\n", nrow(results_df), output_file))
  results_df
}


#' Run score-specific + block models
run_score_specific_plus_blocks_models <- function(score_specific_datasets, datasets_list,
                                                  alpha_grid = seq(0, 1, by = 0.1),
                                                  nfolds = 10, 
                                                  n_permutations = 100,
                                                  seed = 42,
                                                  output_file = "elastic_net_foldwise_preproc_score_plus_blocks.rds") {
  results_list <- list()
  counter <- 1
  n_models <- length(score_specific_datasets) * length(datasets_list)
  
  cat("==========================================================\n")
  cat("ELASTIC NET (FOLD-WISE PREPROCESSING): SCORE-SPECIFIC + PREDICTOR BLOCK\n")
  cat(sprintf("Total models: %d (%d scores × %d blocks)\n",
              n_models, length(score_specific_datasets), length(datasets_list)))
  cat(sprintf("Permutations per model: %d\n", n_permutations))
  cat("==========================================================\n")
  
  for (cvd_score in names(score_specific_datasets)) {
    for (block_name in names(datasets_list)) {
      dataset_name <- paste0(cvd_score, "_specific+", block_name)
      
      cat(sprintf("\n[%d/%d] %s ~ %s\n", counter, n_models, cvd_score, dataset_name))
      cat("----------------------------------------------------------\n")
      
      df_combined <- build_score_plus_block(
        score = cvd_score, block_name = block_name,
        score_specific_datasets = score_specific_datasets,
        datasets_list = datasets_list
      )
      
      if (!cvd_score %in% colnames(df_combined)) {
        cat(sprintf("  ✗ Outcome '%s' not found\n", cvd_score))
        counter <- counter + 1
        next
      }
      
      tryCatch({
        results_list[[counter]] <- run_elastic_net_model(
          df_raw = df_combined,
          outcome_col = cvd_score,
          cvd_score_name = cvd_score,
          dataset_name = dataset_name,
          exclude_cols = "Sample_ID",
          alpha_grid = alpha_grid,
          nfolds = nfolds,
          n_permutations = n_permutations,
          seed = seed
        )
        cat("  ✓ Success\n")
      }, error = function(e) {
        cat(sprintf("  ✗ ERROR: %s\n", e$message))
        results_list[[counter]] <- NULL
      })
      counter <- counter + 1
    }
  }
  
  results_df <- bind_rows(results_list)
  saveRDS(results_df, file = output_file)
  cat(sprintf("\n✓ Saved %d models to %s\n", nrow(results_df), output_file))
  results_df
}


# ============================================
# 8. RESULTS VIEWER
# ============================================
view_results_summary <- function(results_df) {
  cat("\n=== RESULTS SUMMARY (FOLD-WISE PREPROCESSING) ===\n\n")
  
  summary_table <- results_df %>%
    select(cvd_score, dataset_name, n_observations, n_predictors,
           n_nonzero_coefs, alpha_optimal_fulldata, R2_Y, Q2_Y, 
           Q2_R2_ratio, permutation_p_value) %>%
    mutate(across(where(is.numeric), ~round(., 3)))
  
  print(summary_table, n = Inf)
  
  cat("\n=== TOP PREDICTORS PER MODEL ===\n\n")
  for (i in 1:nrow(results_df)) {
    cat(sprintf("%s ~ %s:\n", results_df$cvd_score[i], results_df$dataset_name[i]))
    top_preds <- results_df$top_10_predictors[[i]]
    top_coefs <- results_df$top_10_coefficients[[i]]
    if (!is.na(top_preds[1])) {
      for (j in seq_along(top_preds)) {
        cat(sprintf("  %d. %s (β = %.4f)\n", j, top_preds[j], top_coefs[j]))
      }
    } else {
      cat("  No predictors selected (null model)\n")
    }
    cat("\n")
  }
}


# ============================================
# 9. EXECUTE PIPELINE
# ============================================
# NOTE: With 100 permutations per model, this will take ~2-3 min per model
# Set n_permutations = 0 to skip permutation testing for faster runs

results_score_specific <- run_score_specific_models(
  score_specific_datasets = score_specific_datasets,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 5, 
  n_permutations = 0,
  seed = 42,
  output_file = "elastic_net_foldwise_preproc_score_specific_5fold.rds"
)

results_common <- run_common_datasets_models(
  datasets_list = datasets_list,
  cvd_score_names = CVD_scores,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 5, 
  n_permutations = 0,
  seed = 42,
  output_file = "elastic_net_foldwise_preproc_common_5fold.rds"
)

results_score_plus_blocks <- run_score_specific_plus_blocks_models(
  score_specific_datasets = score_specific_datasets,
  datasets_list = datasets_list,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 5, 
  n_permutations = 0,
  seed = 42,
  output_file = "elastic_net_foldwise_preproc_score_plus_blocks_5fold.rds"
)


# ADDITIONAL: Body Composition + Height/Weight Models
# For scores that don't already include height/weight (SCORE2, ASCVD, Framingham)

# Option A: Body composition + height/weight alone (3 models)
results_body_comp_extended <- list()

for (i in seq_along(scores_without_ht_wt)) {
  score <- scores_without_ht_wt[i]
  
  cat(sprintf("\n[%d/3] %s ~ body_comp_with_ht_wt\n", i, score))

  df_model <- df_body_comp_extended %>%
    filter(!is.na(.data[[score]]))
  
  tryCatch({
    results_body_comp_extended[[i]] <- run_elastic_net_model(
      df_raw = df_model,
      outcome_col = score,
      cvd_score_name = score,
      dataset_name = "body_comp_with_ht_wt",
      exclude_cols = "Sample_ID",
      alpha_grid = seq(0, 1, by = 0.1),
      nfolds = 5,
      n_permutations = 0,
      seed = 42
    )
    cat("  ✓ Success\n")
  }, error = function(e) {
    cat(sprintf("  ✗ ERROR: %s\n", e$message))
    results_body_comp_extended[[i]] <- NULL
  })
}

results_body_comp_extended_df <- bind_rows(results_body_comp_extended)
saveRDS(results_body_comp_extended_df, "elastic_net_body_comp_extended_5fold.rds")

# Option B: Score-specific + body composition with height/weight (3 models)
results_score_plus_body_extended <- list()

for (i in seq_along(scores_without_ht_wt)) {
  score <- scores_without_ht_wt[i]
  dataset_name <- paste0(score, "_specific+body_comp_with_ht_wt")
  
  cat(sprintf("\n[%d/3] %s ~ %s\n", i, score, dataset_name))

  # Get score-specific base
  base_df <- score_specific_datasets[[score]]
  
  # Get body comp predictors (excluding CVD scores)
  body_comp_predictors <- df_body_comp_extended %>%
    select(-all_of(CVD_scores))
  
  # Combine
  df_combined <- base_df %>%
    left_join(body_comp_predictors, by = "Sample_ID") %>%
    filter(!is.na(.data[[score]]))
  
  tryCatch({
    results_score_plus_body_extended[[i]] <- run_elastic_net_model(
      df_raw = df_combined,
      outcome_col = score,
      cvd_score_name = score,
      dataset_name = dataset_name,
      exclude_cols = "Sample_ID",
      alpha_grid = seq(0, 1, by = 0.1),
      nfolds = 5,
      n_permutations = 0,
      seed = 42
    )
    cat("  ✓ Success\n")
  }, error = function(e) {
    cat(sprintf("  ✗ ERROR: %s\n", e$message))
    results_score_plus_body_extended[[i]] <- NULL
  })
}

results_score_plus_body_extended_df <- bind_rows(results_score_plus_body_extended)
saveRDS(results_score_plus_body_extended_df, "elastic_net_score_plus_body_extended_5fold.rds")


# Combine all results
results_all <- bind_rows(results_score_specific, results_common, results_score_plus_blocks,
                         results_body_comp_extended_df, results_score_plus_body_extended_df)

# View summary
view_results_summary(results_all)


# ============================================
# 10. EXPORT TO EXCEL
# ============================================
# Sheet 1: Main results
results_main <- results_all %>%
  select(
    cvd_score, dataset_name, model_name,
    n_observations, n_predictors, n_nonzero_coefs,
    alpha_optimal_fulldata, alpha_mean_nested, lambda_min, lambda_1se,
    R2_Y, in_sample_rmse, in_sample_mae, in_sample_cor,
    Q2_Y, Q2_fold_mean, Q2_fold_sd, Q2_fold_summary,
    Q2_Y_rmse, Q2_Y_mae, Q2_Y_cor,
    Q2_R2_ratio, R2_Q2_gap,
    permutation_p_value, permutation_n, permutation_null_mean_Q2,
    cv_mse_min, cv_mse_1se, pct_deviance_explained,
    imputation_method, preprocessing, n_folds, seed, date_run
  )

# Sheet 2: Top predictors
predictor_rows <- list()
for (i in 1:nrow(results_all)) {
  preds <- results_all$top_10_predictors[[i]]
  coefs <- results_all$top_10_coefficients[[i]]
  if (!is.null(preds) && length(preds) > 0 && !all(is.na(preds))) {
    for (j in seq_along(preds)) {
      if (!is.na(preds[j])) {
        predictor_rows[[length(predictor_rows) + 1]] <- tibble(
          cvd_score = results_all$cvd_score[i],
          dataset_name = results_all$dataset_name[i],
          rank = j,
          predictor_name = preds[j],
          coefficient = round(coefs[j], 4)
        )
      }
    }
  }
}
top_predictors_sheet <- bind_rows(predictor_rows)

# Sheet 3: Q²_Y overview
table_A_all_models <- results_all %>%
  select(dataset_name, cvd_score, Q2_Y) %>%
  mutate(Q2_Y = round(Q2_Y, 3)) %>%
  distinct(dataset_name, cvd_score, .keep_all = TRUE) %>%
  pivot_wider(names_from = cvd_score, values_from = Q2_Y) %>%
  mutate(
    dataset_group = case_when(
      grepl("_specific\\+", dataset_name) ~ "specific+block",
      grepl("_specific$", dataset_name) ~ "specific",
      TRUE ~ "block_only"
    )
  ) %>%
  arrange(factor(dataset_group, levels = c("specific", "block_only", "specific+block")), dataset_name) %>%
  select(-dataset_group)

baseline_specific <- results_all %>%
  filter(dataset_name == paste0(cvd_score, "_specific")) %>%
  select(cvd_score, baseline_Q2_Y = Q2_Y)

table_B_delta <- results_all %>%
  filter(grepl("_specific\\+", dataset_name)) %>%
  mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
  select(cvd_score, block, Q2_Y) %>%
  left_join(baseline_specific, by = "cvd_score") %>%
  mutate(delta_Q2_Y = round(Q2_Y - baseline_Q2_Y, 3)) %>%
  select(block, cvd_score, delta_Q2_Y) %>%
  pivot_wider(names_from = cvd_score, values_from = delta_Q2_Y) %>%
  arrange(match(block, names(datasets_list)))

# Create workbook
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

addWorksheet(wb, "Q2_Y_Overview")
writeData(wb, "Q2_Y_Overview", "Table A: Q²_Y (all models) - Fold-wise Preprocessing", startRow = 1, startCol = 1)
writeData(wb, "Q2_Y_Overview", table_A_all_models, startRow = 2, startCol = 1)

start_row_B <- nrow(table_A_all_models) + 5
writeData(wb, "Q2_Y_Overview",
          "Table B: Δ Q²_Y = (score-specific + block) − (score-specific baseline)",
          startRow = start_row_B, startCol = 1)
writeData(wb, "Q2_Y_Overview", table_B_delta, startRow = start_row_B + 1, startCol = 1)

addStyle(wb, "Q2_Y_Overview", headerStyle,
         rows = c(2, start_row_B + 1),
         cols = 1:max(ncol(table_A_all_models), ncol(table_B_delta)),
         gridExpand = TRUE)
setColWidths(wb, "Q2_Y_Overview",
             cols = 1:max(ncol(table_A_all_models), ncol(table_B_delta)),
             widths = "auto")

saveWorkbook(wb, "elastic_net_foldwise_preproc_results_5fold.xlsx", overwrite = TRUE)

cat("\n✓ Excel saved: elastic_net_foldwise_preproc_results.xlsx\n")
cat("  - Sheet 1: All_Results (", nrow(results_main), " models)\n")
cat("  - Sheet 2: Top_Predictors (", nrow(top_predictors_sheet), " rows)\n")
cat("  - Sheet 3: Q2_Y_Overview (Fold-wise Preprocessing)\n")



# ============================================
# FIXED TABLE AND PLOT CODE FOR ELASTIC NET RESULTS
# Run this AFTER the main elastic net analysis completes
# ============================================

# ============================================
# REQUIRED VARIABLE DEFINITIONS (were missing)
# ============================================

# CVD score ordering and labels
cvd_score_order <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y", "mean_risk")

cvd_score_labels <- c(
  "QRISK3_risk" = "QRISK3",
  "SCORE2_score" = "SCORE2", 
  "frs_10y" = "Framingham",
  "ascvd_10y" = "ASCVD",
  "mean_risk" = "Composite"
)

# Dataset labels for tables
dataset_labels_table <- c(
  "all_data" = "All predictors",
  "metabolites" = "All metabolites",
  
  "mtc_sign" = "BH-significant",
  "sign_wo_mtc" = "Nominally significant",
  "body_composition" = "Body composition",
  "sociodemographics_lifestyle" = "Sociodemographics & lifestyle",
  "clinical_risk_factors" = "Clinical risk factors",
  "lipids" = "Lipids",
  "fatty_acids" = "Fatty acids",
  "urine_nmr" = "Urine NMR",
  "body_comp_with_ht_wt" = "Body composition + Ht/Wt"
)

# ============================================
# TABLE CREATION FUNCTIONS
# ============================================

#' Create table for score-specific models (5 rows, one per CVD score)
create_score_specific_table <- function(results_df, table_title) {
  
  # Prepare data
  table_data <- results_df %>%
    filter(cvd_score %in% cvd_score_order) %>%
    mutate(
      cvd_score = factor(cvd_score, levels = cvd_score_order),
      CVD_Score = cvd_score_labels[as.character(cvd_score)]
    ) %>%
    arrange(cvd_score)
  
  # Check if fold-wise SDs exist
  has_mae_sd <- "MAE_fold_sd" %in% names(table_data) && any(!is.na(table_data$MAE_fold_sd))
  has_rmse_sd <- "RMSE_fold_sd" %in% names(table_data) && any(!is.na(table_data$RMSE_fold_sd))
  
  # Format metrics
  formatted_table <- table_data %>%
    mutate(
      # R² (in-sample, no SD available)
      `R²` = sprintf("%.3f", R2_Y),
      
      # Q² with fold SD
      `Q²` = ifelse(!is.na(Q2_fold_sd) & Q2_fold_sd > 0,
                    sprintf("%.3f ± %.3f", Q2_Y, Q2_fold_sd),
                    sprintf("%.3f", Q2_Y)),
      
      # MAE with fold SD if available
      MAE = if (has_mae_sd) {
        ifelse(!is.na(MAE_fold_sd) & MAE_fold_sd > 0,
               sprintf("%.3f ± %.3f", Q2_Y_mae, MAE_fold_sd),
               sprintf("%.3f", Q2_Y_mae))
      } else {
        sprintf("%.3f", Q2_Y_mae)
      },
      
      # RMSE with fold SD if available
      RMSE = if (has_rmse_sd) {
        ifelse(!is.na(RMSE_fold_sd) & RMSE_fold_sd > 0,
               sprintf("%.3f ± %.3f", Q2_Y_rmse, RMSE_fold_sd),
               sprintf("%.3f", Q2_Y_rmse))
      } else {
        sprintf("%.3f", Q2_Y_rmse)
      },
      
      # Permutation p-value
      `p-value` = ifelse(is.na(permutation_p_value), 
                         "—",
                         ifelse(permutation_p_value < 0.001, 
                                "<0.001",
                                sprintf("%.3f", permutation_p_value)))
    ) %>%
    select(CVD_Score, `R²`, `Q²`, MAE, RMSE, `p-value`)
  
  # Create flextable
  ft <- formatted_table %>%
    flextable() %>%
    set_caption(caption = table_title) %>%
    set_header_labels(CVD_Score = "CVD Score") %>%
    align(j = 1, align = "left", part = "all") %>%
    align(j = 2:6, align = "center", part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    fontsize(size = 9, part = "body") %>%
    fontsize(size = 10, part = "header") %>%
    bold(part = "header") %>%
    width(j = 1, width = 1.2) %>%
    width(j = 2, width = 0.8) %>%
    width(j = 3, width = 1.2) %>%
    width(j = 4, width = 1.2) %>%
    width(j = 5, width = 1.2) %>%
    width(j = 6, width = 0.8) %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1.5), part = "header") %>%
    padding(padding = 3, part = "all")
  
  return(ft)
}


#' Create table for multiple CVD scores across datasets
create_multi_score_table <- function(results_df, table_title, 
                                     dataset_order = NULL,
                                     dataset_labels = dataset_labels_table) {
  
  # Prepare data
  table_data <- results_df %>%
    filter(cvd_score %in% cvd_score_order) %>%
    mutate(
      cvd_score = factor(cvd_score, levels = cvd_score_order),
      CVD_Score = cvd_score_labels[as.character(cvd_score)]
    )
  
  # Apply dataset ordering if provided
  if (!is.null(dataset_order)) {
    table_data <- table_data %>%
      filter(dataset_name %in% dataset_order) %>%
      mutate(dataset_name = factor(dataset_name, levels = dataset_order))
  }
  
  # Check if we have data
  if (nrow(table_data) == 0) {
    warning("No data remaining after filtering!")
    return(NULL)
  }
  
  table_data <- table_data %>%
    arrange(dataset_name, cvd_score)
  
  # Check if fold-wise SDs exist
  has_mae_sd <- "MAE_fold_sd" %in% names(table_data) && any(!is.na(table_data$MAE_fold_sd))
  has_rmse_sd <- "RMSE_fold_sd" %in% names(table_data) && any(!is.na(table_data$RMSE_fold_sd))
  
  # Format metrics
  formatted_table <- table_data %>%
    mutate(
      Dataset = ifelse(dataset_name %in% names(dataset_labels),
                       dataset_labels[as.character(dataset_name)],
                       as.character(dataset_name)),
      
      # R² (in-sample)
      `R²` = sprintf("%.3f", R2_Y),
      
      # Q² with fold SD
      `Q²` = ifelse(!is.na(Q2_fold_sd) & Q2_fold_sd > 0,
                    sprintf("%.3f ± %.3f", Q2_Y, Q2_fold_sd),
                    sprintf("%.3f", Q2_Y)),
      
      # MAE with fold SD if available
      MAE = if (has_mae_sd) {
        ifelse(!is.na(MAE_fold_sd) & MAE_fold_sd > 0,
               sprintf("%.3f ± %.3f", Q2_Y_mae, MAE_fold_sd),
               sprintf("%.3f", Q2_Y_mae))
      } else {
        sprintf("%.3f", Q2_Y_mae)
      },
      
      # RMSE with fold SD if available
      RMSE = if (has_rmse_sd) {
        ifelse(!is.na(RMSE_fold_sd) & RMSE_fold_sd > 0,
               sprintf("%.3f ± %.3f", Q2_Y_rmse, RMSE_fold_sd),
               sprintf("%.3f", Q2_Y_rmse))
      } else {
        sprintf("%.3f", Q2_Y_rmse)
      },
      
      # Permutation p-value
      `p-value` = ifelse(is.na(permutation_p_value), 
                         "—",
                         ifelse(permutation_p_value < 0.001, 
                                "<0.001",
                                sprintf("%.3f", permutation_p_value)))
    ) %>%
    select(Dataset, CVD_Score, `R²`, `Q²`, MAE, RMSE, `p-value`)
  
  n_rows <- nrow(formatted_table)
  
  # Create flextable
  ft <- formatted_table %>%
    flextable() %>%
    merge_v(j = "Dataset") %>%
    set_caption(caption = table_title) %>%
    set_header_labels(CVD_Score = "CVD Score") %>%
    align(j = 1:2, align = "left", part = "all") %>%
    align(j = 3:7, align = "center", part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    fontsize(size = 9, part = "body") %>%
    fontsize(size = 10, part = "header") %>%
    bold(part = "header") %>%
    bold(j = 1, part = "body") %>%
    width(j = 1, width = 1.8) %>%
    width(j = 2, width = 1.0) %>%
    width(j = 3, width = 0.7) %>%
    width(j = 4, width = 1.2) %>%
    width(j = 5, width = 1.2) %>%
    width(j = 6, width = 1.2) %>%
    width(j = 7, width = 0.7) %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1.5), part = "header") %>%
    padding(padding = 3, part = "all")
  
  # Add horizontal lines between datasets (every 5 rows for 5 CVD scores)
  if (n_rows > 5) {
    hline_rows <- seq(5, n_rows - 1, by = 5)
    if (length(hline_rows) > 0) {
      ft <- ft %>%
        hline(i = hline_rows, 
              border = fp_border(width = 0.5, color = "#CCCCCC"), part = "body")
    }
  }
  
  return(ft)
}


# ============================================
# CREATE THE 4 TABLES
# ============================================

cat("\n==========================================================\n")
cat("CREATING PUBLICATION TABLES\n")
cat("==========================================================\n")

# TABLE 1: Score-specific only (simple table, 5 rows)
cat("\nTable 1: Score-specific models...\n")
table1 <- create_score_specific_table(
  results_df = results_score_specific,
  table_title = "Table 1: Elastic Net Performance – Score-Specific Input Variables"
)
cat("  ✓ Created\n")

# TABLE 2: Common datasets (excluding mtc_sign and sign_wo_mtc)
cat("\nTable 2: Common datasets...\n")
datasets_table2 <- c("all_data", "metabolites", "lipids", "fatty_acids", "urine_nmr", 
                     "body_composition", "sociodemographics_lifestyle", "clinical_risk_factors")

# Filter to only include datasets that exist in results
datasets_table2 <- datasets_table2[datasets_table2 %in% unique(results_common$dataset_name)]

table2 <- create_multi_score_table(
  results_df = results_common %>% 
    filter(dataset_name %in% datasets_table2),
  table_title = "Table 2: Elastic Net Performance – Predictor Blocks",
  dataset_order = datasets_table2
)
cat(sprintf("  ✓ Created (%d datasets)\n", length(datasets_table2)))

# TABLE 3: Score-specific + block (excluding mtc_sign and sign_wo_mtc)
cat("\nTable 3: Score-specific + predictor blocks...\n")

results_table3 <- results_score_plus_blocks %>%
  mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
  filter(!block %in% c("mtc_sign", "sign_wo_mtc"))

blocks_table3 <- c("all_data", "metabolites", "lipids", "fatty_acids", "urine_nmr", 
                   "body_composition", "sociodemographics_lifestyle", "clinical_risk_factors")

dataset_order_table3 <- results_table3 %>%
  filter(sub("^.*_specific\\+", "", dataset_name) %in% blocks_table3) %>%
  arrange(match(sub("^.*_specific\\+", "", dataset_name), blocks_table3)) %>%
  pull(dataset_name) %>%
  unique()

# Create labels for score+block datasets
dataset_labels_table3 <- dataset_labels_table
for (ds in dataset_order_table3) {
  block <- sub("^.*_specific\\+", "", ds)
  block_label <- ifelse(block %in% names(dataset_labels_table), 
                        dataset_labels_table[block], block)
  dataset_labels_table3[ds] <- paste0("Score inputs + ", block_label)
}

table3 <- create_multi_score_table(
  results_df = results_table3,
  table_title = "Table 3: Elastic Net Performance – Score-Specific + Predictor Blocks",
  dataset_order = dataset_order_table3,
  dataset_labels = dataset_labels_table3
)
cat(sprintf("  ✓ Created (%d dataset combinations)\n", length(dataset_order_table3)))

# TABLE 4: Significance-based predictors
cat("\nTable 4: Significance-based predictor sets...\n")

results_sig_common <- results_common %>%
  filter(dataset_name %in% c("mtc_sign", "sign_wo_mtc"))

results_sig_plus <- results_score_plus_blocks %>%
  mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
  filter(block %in% c("mtc_sign", "sign_wo_mtc"))

results_table4 <- bind_rows(results_sig_common, results_sig_plus)

if (nrow(results_table4) > 0) {
  dataset_order_table4 <- c(
    intersect(c("mtc_sign", "sign_wo_mtc"), unique(results_sig_common$dataset_name)),
    results_sig_plus %>% 
      arrange(match(sub("^.*_specific\\+", "", dataset_name), c("mtc_sign", "sign_wo_mtc"))) %>%
      pull(dataset_name) %>% unique()
  )
  
  dataset_labels_table4 <- dataset_labels_table
  for (ds in dataset_order_table4) {
    if (!ds %in% names(dataset_labels_table4)) {
      block <- sub("^.*_specific\\+", "", ds)
      block_label <- ifelse(block == "mtc_sign", "BH-significant", "Nominally significant")
      dataset_labels_table4[ds] <- paste0("Score inputs + ", block_label)
    }
  }
  
  table4 <- create_multi_score_table(
    results_df = results_table4,
    table_title = "Table 4: Elastic Net Performance – Significance-Based Predictor Sets",
    dataset_order = dataset_order_table4,
    dataset_labels = dataset_labels_table4
  )
  cat(sprintf("  ✓ Created (%d rows)\n", nrow(results_table4)))
} else {
  cat("  ⚠ No significance-based results found, skipping Table 4\n")
  table4 <- NULL
}


# ============================================
# SAVE TABLES
# ============================================

cat("\nSaving tables...\n")

save_as_docx(table1, path = "Table1_score_specific.docx")
save_as_docx(table2, path = "Table2_predictor_blocks.docx")
save_as_docx(table3, path = "Table3_score_plus_blocks.docx")
if (!is.null(table4)) {
  save_as_docx(table4, path = "Table4_significance_based.docx")
}

save_as_html(table1, path = "Table1_score_specific.html")
save_as_html(table2, path = "Table2_predictor_blocks.html")
save_as_html(table3, path = "Table3_score_plus_blocks.html")
if (!is.null(table4)) {
  save_as_html(table4, path = "Table4_significance_based.html")
}

# Combined document
doc <- read_docx() %>%
  body_add_par("Elastic Net Regression Results", style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(value = table1) %>%
  body_add_break() %>%
  body_add_flextable(value = table2) %>%
  body_add_break() %>%
  body_add_flextable(value = table3)

if (!is.null(table4)) {
  doc <- doc %>%
    body_add_break() %>%
    body_add_flextable(value = table4)
}

print(doc, target = "Elastic_Net_All_Tables.docx")

cat("\n✓ Tables saved:\n")
cat("  - Table1_score_specific.docx/.html\n")
cat("  - Table2_predictor_blocks.docx/.html\n")
cat("  - Table3_score_plus_blocks.docx/.html\n")
if (!is.null(table4)) {
  cat("  - Table4_significance_based.docx/.html\n")
}
cat("  - Elastic_Net_All_Tables.docx (combined)\n")


# ============================================
# VISUALIZATION: ELASTIC NET PERFORMANCE
# ============================================

cat("\n==========================================================\n")
cat("CREATING VISUALIZATION PLOTS\n")
cat("==========================================================\n")

# CVD score color palette
cvd_colors <- c(
  "QRISK3" = "#E69F00",
  "SCORE2" = "#56B4E9",
  "Framingham" = "#009E73",
  "ASCVD" = "#F0E442",
  "Composite" = "#D55E00"
)

# ============================================
# PREPARE DATA FOR PANELS 1 & 3 (LEFT COLUMN)
# Score-specific (grouped as one) + common blocks (all scores)
# ============================================

# Score-specific: group all 5 CVD scores into ONE "Score-Specific" row
score_specific_grouped <- results_score_specific %>%
  filter(cvd_score %in% cvd_score_order) %>%
  mutate(
    dataset_name = "score_specific",
    dataset_label = "Score-Specific"
  )

# Common blocks: show all 5 CVD scores for each dataset
# EXCLUDE "all_data" (All metabolites)
datasets_for_left <- c("metabolites", "lipids", "fatty_acids", "urine_nmr", 
                       "body_composition", "clinical_risk_factors", 
                       "sociodemographics_lifestyle")

# Filter to existing datasets
datasets_for_left <- datasets_for_left[datasets_for_left %in% unique(results_common$dataset_name)]

# Build dataset labels for left panel
dataset_labels_left <- c(
  "score_specific" = "Score-Specific",
  dataset_labels_table[datasets_for_left]
)

common_blocks_all_scores <- results_common %>%
  filter(
    dataset_name %in% datasets_for_left,
    cvd_score %in% cvd_score_order
  )

# Define dataset order for left panels (REVERSED - bottom to top)
dataset_order_left <- rev(c("score_specific", datasets_for_left))

plot_data_left <- bind_rows(
  score_specific_grouped,
  common_blocks_all_scores
) %>%
  mutate(
    cvd_score = factor(cvd_score, levels = cvd_score_order),
    cvd_label = cvd_score_labels[as.character(cvd_score)],
    cvd_label = factor(cvd_label, levels = cvd_score_labels),
    dataset_name = factor(dataset_name, levels = dataset_order_left),
    dataset_label = ifelse(dataset_name %in% names(dataset_labels_left),
                           dataset_labels_left[as.character(dataset_name)],
                           as.character(dataset_name)),
    dataset_label = factor(dataset_label, levels = unique(dataset_labels_left[dataset_order_left]))
  )

# ============================================
# PREPARE DATA FOR PANELS 2 & 4 (RIGHT COLUMN)
# Score-specific (baseline) + grouped by block
# ============================================

# Get score-specific baseline grouped as ONE row
baseline_grouped <- results_score_specific %>%
  filter(cvd_score %in% cvd_score_order) %>%
  mutate(
    dataset_name = "score_specific",
    dataset_label = "Score-Specific",
    block = "score_specific"
  )

# Get score+block combinations - EXCLUDE "all_data"
blocks_for_plus <- c("metabolites", "lipids", "fatty_acids", "urine_nmr", 
                     "body_composition", "clinical_risk_factors", 
                     "sociodemographics_lifestyle")

score_plus_data <- results_score_plus_blocks %>%
  mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
  filter(
    block %in% blocks_for_plus,
    cvd_score %in% cvd_score_order
  ) %>%
  mutate(
    dataset_name = block,  # Group by block, not by score+block combination
    dataset_label = ifelse(block %in% names(dataset_labels_table),
                           dataset_labels_table[block], block)
  )

# Define dataset order for right panels (REVERSED - bottom to top)
dataset_order_right <- rev(c("score_specific", blocks_for_plus))

# Create labels for right panel
dataset_labels_right <- c(
  "score_specific" = "Score-Specific",
  dataset_labels_table[blocks_for_plus]
)

plot_data_right <- bind_rows(baseline_grouped, score_plus_data) %>%
  filter(dataset_name %in% dataset_order_right) %>%
  mutate(
    cvd_score = factor(cvd_score, levels = cvd_score_order),
    cvd_label = cvd_score_labels[as.character(cvd_score)],
    cvd_label = factor(cvd_label, levels = cvd_score_labels),
    dataset_name = factor(dataset_name, levels = dataset_order_right),
    dataset_label = ifelse(dataset_name %in% names(dataset_labels_right),
                           dataset_labels_right[as.character(dataset_name)],
                           as.character(dataset_name)),
    dataset_label = factor(dataset_label, levels = unique(dataset_labels_right[dataset_order_right]))
  )

# ============================================
# PANEL 1: Left column - Q²
# ============================================

cat("Creating Q² panel (left)...\n")

p1_left_q2 <- plot_data_left %>%
  ggplot(aes(x = Q2_Y, y = dataset_label, fill = cvd_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(
    aes(xmin = pmax(0, Q2_Y - Q2_fold_sd), xmax = Q2_Y + Q2_fold_sd),
    position = position_dodge(width = 0.8),
    height = 0.3,
    linewidth = 0.4
  ) +
  geom_text(
    aes(x = Q2_Y + Q2_fold_sd, label = sprintf("%.2f", Q2_Y)),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 2.2
  ) +
  scale_fill_manual(values = cvd_colors) +
  scale_x_continuous(
    limits = c(0, 1.0),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "Score-Specific & Predictor Blocks",
    x = expression(Q^2),
    y = NULL,
    fill = "CVD Score"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# ============================================
# PANEL 2: Right column - Q²
# ============================================

cat("Creating Q² panel (right)...\n")

p2_right_q2 <- plot_data_right %>%
  ggplot(aes(x = Q2_Y, y = dataset_label, fill = cvd_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(
    aes(xmin = pmax(0, Q2_Y - Q2_fold_sd), xmax = Q2_Y + Q2_fold_sd),
    position = position_dodge(width = 0.8),
    height = 0.3,
    linewidth = 0.4
  ) +
  geom_text(
    aes(x = Q2_Y + Q2_fold_sd, label = sprintf("%.2f", Q2_Y)),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 2.2
  ) +
  scale_fill_manual(values = cvd_colors) +
  scale_x_continuous(
    limits = c(0, 1.0),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "Score-Specific + Predictor Blocks",
    x = expression(Q^2),
    y = NULL,
    fill = "CVD Score"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold")
  )

# ============================================
# PANEL 3: Left column - MAE
# ============================================

cat("Creating MAE panel (left)...\n")

# Handle potential missing MAE_fold_sd column
if (!"MAE_fold_sd" %in% names(plot_data_left)) {
  plot_data_left$MAE_fold_sd <- 0
}

p3_left_mae <- plot_data_left %>%
  ggplot(aes(x = Q2_Y_mae, y = dataset_label, fill = cvd_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(
    aes(xmin = pmax(0, Q2_Y_mae - MAE_fold_sd), xmax = Q2_Y_mae + MAE_fold_sd),
    position = position_dodge(width = 0.8),
    height = 0.3,
    linewidth = 0.4
  ) +
  geom_text(
    aes(x = Q2_Y_mae + MAE_fold_sd, label = sprintf("%.2f", Q2_Y_mae)),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 2.2
  ) +
  scale_fill_manual(values = cvd_colors) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "",
    x = "MAE",
    y = NULL,
    fill = "CVD Score"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# ============================================
# PANEL 4: Right column - MAE
# ============================================

cat("Creating MAE panel (right)...\n")

# Handle potential missing MAE_fold_sd column
if (!"MAE_fold_sd" %in% names(plot_data_right)) {
  plot_data_right$MAE_fold_sd <- 0
}

p4_right_mae <- plot_data_right %>%
  ggplot(aes(x = Q2_Y_mae, y = dataset_label, fill = cvd_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(
    aes(xmin = pmax(0, Q2_Y_mae - MAE_fold_sd), xmax = Q2_Y_mae + MAE_fold_sd),
    position = position_dodge(width = 0.8),
    height = 0.3,
    linewidth = 0.4
  ) +
  geom_text(
    aes(x = Q2_Y_mae + MAE_fold_sd, label = sprintf("%.2f", Q2_Y_mae)),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 2.2
  ) +
  scale_fill_manual(values = cvd_colors) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "",
    x = "MAE",
    y = NULL,
    fill = "CVD Score"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold")
  )

# ============================================
# COMBINE WITH PATCHWORK
# ============================================

cat("Combining panels...\n")

combined_elastic_plot <- (p1_left_q2 | p2_right_q2) / 
  (p3_left_mae | p4_right_mae) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Elastic Net Performance: Cross-Validated Predictive Ability",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  )

# Display
print(combined_elastic_plot)

# ============================================
# SAVE PLOTS
# ============================================

ggsave(
  "elastic_net_performance_comparison.png",
  plot = combined_elastic_plot,
  width = 16,
  height = 12,
  dpi = 300,
  bg = "white"
)

ggsave(
  "elastic_net_performance_comparison.pdf",
  plot = combined_elastic_plot,
  width = 16,
  height = 12,
  device = "pdf"
)

cat("\n✓ Elastic Net performance plot saved as:\n")
cat("  - elastic_net_performance_comparison.png\n")
cat("  - elastic_net_performance_comparison.pdf\n")

cat("\n==========================================================\n")
cat("TABLE AND PLOT CREATION COMPLETE\n")
cat("==========================================================\n")


#' # ============================================
#' # XY. VISUALIZE RESIDUALS
#' # ============================================
# ============================================
# OBSERVED VS. PREDICTED: THREE FIGURE SET
# ============================================

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

cat("\n==========================================================\n")
cat("CREATING OBSERVED VS. PREDICTED FIGURES\n")
cat("==========================================================\n")

# CVD score color palette
cvd_colors <- c(
  "QRISK3_risk" = "#E69F00",
  "SCORE2_score" = "#56B4E9",
  "frs_10y" = "#009E73",
  "ascvd_10y" = "#F0E442",
  "mean_risk" = "#D55E00"
)

cvd_score_order <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y", "mean_risk")

cvd_score_labels <- c(
  "QRISK3_risk" = "QRISK3",
  "SCORE2_score" = "SCORE2",
  "frs_10y" = "Framingham",
  "ascvd_10y" = "ASCVD",
  "mean_risk" = "Composite"
)

# ============================================
# FUNCTION: Create single observed vs predicted panel
# ============================================

create_obs_pred_panel <- function(residuals_subset, cvd_score, color) {
  
  # Calculate metrics
  R2 <- cor(residuals_subset$observed, residuals_subset$predicted)^2
  RMSE <- sqrt(mean(residuals_subset$residual^2))
  MAE <- mean(abs(residuals_subset$residual))
  
  # Get label
  score_label <- cvd_score_labels[cvd_score]
  
  # Create plot
  p <- ggplot(residuals_subset, aes(x = predicted, y = observed)) +
    geom_point(alpha = 0.6, size = 1, color = color) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                linewidth = 0.8, color = "black") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, 
                alpha = 0.2, color = color, fill = color) +
    annotate("text",
             x = -Inf, y = Inf,
             label = sprintf("R² = %.3f\nRMSE = %.2f\nMAE = %.2f", R2, RMSE, MAE),
             hjust = -0.1, vjust = 1.1,
             size = 3, color = "black") +
    labs(
      title = score_label,
      x = "Predicted",
      y = "Observed"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      aspect.ratio = 1
    )
  
  return(p)
}

# ============================================
# FIGURE 1: SCORE-SPECIFIC MODELS (1 row × 5 columns)
# ============================================

cat("\n--- Creating Figure 1: Score-Specific Models ---\n")

# Get all score-specific rows
score_specific_rows <- which(
  results_all$dataset_name %in% paste0(cvd_score_order, "_specific")
)

# Extract residuals for all score-specific models
residuals_score_specific <- extract_multiple_residuals(
  results_df = results_all,
  row_indices = score_specific_rows,
  datasets_list = datasets_list,
  score_specific_datasets = score_specific_datasets
)

# Create individual plots
plot_list_fig1 <- list()
for (i in seq_along(cvd_score_order)) {
  score <- cvd_score_order[i]
  color <- cvd_colors[score]
  
  # Filter data for this score
  data_subset <- residuals_score_specific %>%
    filter(cvd_score == score)
  
  if (nrow(data_subset) > 0) {
    plot_list_fig1[[i]] <- create_obs_pred_panel(data_subset, score, color)
  }
}

# Combine into single row
fig1_combined <- wrap_plots(plot_list_fig1, nrow = 1) +
  plot_annotation(
    title = "Observed vs. Predicted: Score-Specific Inputs",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

print(fig1_combined)

ggsave(
  "figure1_observed_vs_predicted_score_specific.png",
  plot = fig1_combined,
  width = 18,
  height = 4,
  dpi = 300,
  bg = "white"
)

cat("✓ Figure 1 saved\n")

# ============================================
# FIGURE 2: BLOCK PREDICTORS (rows × 5 columns)
# ============================================

cat("\n--- Creating Figure 2: Block Predictors ---\n")

# Define blocks to include (exclude all_data, mtc_sign, sign_wo_mtc)
blocks_to_plot <- c(
  "metabolites",
  "lipids", 
  "fatty_acids",
  "urine_nmr",
  "body_composition",
  "clinical_risk_factors",
  "sociodemographics_lifestyle"
)

# Filter to blocks that exist
blocks_to_plot <- blocks_to_plot[blocks_to_plot %in% unique(results_all$dataset_name)]

# Get all block predictor rows
block_rows <- which(
  results_all$dataset_name %in% blocks_to_plot
)

# Extract residuals
residuals_blocks <- extract_multiple_residuals(
  results_df = results_all,
  row_indices = block_rows,
  datasets_list = datasets_list,
  score_specific_datasets = score_specific_datasets
)

# Create grid of plots (rows = blocks, columns = scores)
plot_list_fig2 <- list()
plot_index <- 1

for (block in blocks_to_plot) {
  for (score in cvd_score_order) {
    color <- cvd_colors[score]
    
    # Filter data
    data_subset <- residuals_blocks %>%
      filter(dataset == block, cvd_score == score)
    
    if (nrow(data_subset) > 0) {
      p <- create_obs_pred_panel(data_subset, score, color)
      
      # Add y-axis label only for first column
      if (score == cvd_score_order[1]) {
        block_label <- dataset_labels_table[block]
        p <- p + labs(y = paste0(block_label, "\n\nObserved"))
      } else {
        p <- p + labs(y = "")
      }
      
      # Remove title for all except top row
      if (block != blocks_to_plot[1]) {
        p <- p + labs(title = "")
      }
      
      plot_list_fig2[[plot_index]] <- p
      plot_index <- plot_index + 1
    } else {
      # Add empty plot if no data
      plot_list_fig2[[plot_index]] <- ggplot() + theme_void()
      plot_index <- plot_index + 1
    }
  }
}

# Combine into grid
fig2_combined <- wrap_plots(plot_list_fig2, ncol = 5, byrow = TRUE) +
  plot_annotation(
    title = "Observed vs. Predicted: Predictor Blocks Across CVD Scores",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

print(fig2_combined)

ggsave(
  "figure2_observed_vs_predicted_blocks.png",
  plot = fig2_combined,
  width = 18,
  height = length(blocks_to_plot) * 3.5,
  dpi = 300,
  bg = "white"
)

cat("✓ Figure 2 saved\n")

# ============================================
# FIGURE 3: SCORE-SPECIFIC + BLOCK (rows × 5 columns)
# ============================================

cat("\n--- Creating Figure 3: Score-Specific + Block Predictors ---\n")

# Get all score + block rows
score_plus_block_rows <- which(
  grepl("_specific\\+", results_all$dataset_name) &
    grepl(paste0("\\+(", paste(blocks_to_plot, collapse = "|"), ")$"), results_all$dataset_name)
)

# Extract residuals
residuals_score_plus <- extract_multiple_residuals(
  results_df = results_all,
  row_indices = score_plus_block_rows,
  datasets_list = datasets_list,
  score_specific_datasets = score_specific_datasets
)

# Add block column for grouping
residuals_score_plus <- residuals_score_plus %>%
  mutate(block = sub("^.*_specific\\+", "", dataset))

# Create grid of plots (rows = blocks, columns = scores)
plot_list_fig3 <- list()
plot_index <- 1

for (block in blocks_to_plot) {
  for (score in cvd_score_order) {
    color <- cvd_colors[score]
    
    # Filter data
    data_subset <- residuals_score_plus %>%
      filter(block == !!block, cvd_score == score)
    
    if (nrow(data_subset) > 0) {
      p <- create_obs_pred_panel(data_subset, score, color)
      
      # Add y-axis label only for first column
      if (score == cvd_score_order[1]) {
        block_label <- dataset_labels_table[block]
        p <- p + labs(y = paste0("+ ", block_label, "\n\nObserved"))
      } else {
        p <- p + labs(y = "")
      }
      
      # Remove title for all except top row
      if (block != blocks_to_plot[1]) {
        p <- p + labs(title = "")
      }
      
      plot_list_fig3[[plot_index]] <- p
      plot_index <- plot_index + 1
    } else {
      # Add empty plot if no data
      plot_list_fig3[[plot_index]] <- ggplot() + theme_void()
      plot_index <- plot_index + 1
    }
  }
}

# Combine into grid
fig3_combined <- wrap_plots(plot_list_fig3, ncol = 5, byrow = TRUE) +
  plot_annotation(
    title = "Observed vs. Predicted: Score-Specific + Predictor Blocks",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

print(fig3_combined)

ggsave(
  "figure3_observed_vs_predicted_score_plus_blocks.png",
  plot = fig3_combined,
  width = 18,
  height = length(blocks_to_plot) * 3.5,
  dpi = 300,
  bg = "white"
)

cat("✓ Figure 3 saved\n")

cat("\n==========================================================\n")
cat("ALL FIGURES COMPLETE\n")
cat("==========================================================\n")
cat("Figure 1: Score-Specific (1 row × 5 CVD scores)\n")
cat("Figure 2: Predictor Blocks (", length(blocks_to_plot), " rows × 5 CVD scores)\n")
cat("Figure 3: Score-Specific + Blocks (", length(blocks_to_plot), " rows × 5 CVD scores)\n")
cat("==========================================================\n")




#' # ============================================
#' # 11. VISUALIZE NON-ZERO COEFFICIENTS
#' # ============================================
#'
#' #' Extract all non-zero coefficients from model results
#' extract_nonzero_coefs <- function(results_df) {
#'   coef_list <- list()
#'
#'   for (i in 1:nrow(results_df)) {
#'     predictor_names <- results_df$all_predictor_names[[i]]
#'     coefficients <- results_df$all_coefficients[[i]]
#'
#'     # Get non-zero coefficients
#'     nonzero_mask <- coefficients != 0
#'
#'     if (sum(nonzero_mask) > 0) {
#'       coef_list[[i]] <- tibble(
#'         cvd_score = results_df$cvd_score[i],
#'         dataset_name = results_df$dataset_name[i],
#'         predictor = predictor_names[nonzero_mask],
#'         coefficient = coefficients[nonzero_mask],
#'         abs_coefficient = abs(coefficients[nonzero_mask]),
#'         direction = ifelse(coefficients[nonzero_mask] > 0, "Positive", "Negative")
#'       )
#'     }
#'   }
#'
#'   bind_rows(coef_list)
#' }
#'
#'
#' #' Create bar plot for a single model's non-zero coefficients
#' plot_nonzero_coefs <- function(coef_data, cvd_score, dataset_name) {
#'   # Order predictors by absolute coefficient value
#'   coef_data <- coef_data %>%
#'     arrange(abs_coefficient) %>%
#'     mutate(predictor = factor(predictor, levels = predictor))
#'
#'   # Create plot
#'   p <- ggplot(coef_data, aes(x = predictor, y = abs_coefficient, fill = direction)) +
#'     geom_bar(stat = "identity", width = 0.7) +
#'     scale_fill_manual(values = c("Positive" = "#FF8C00", "Negative" = "#4169E1"),
#'                       name = "Direction") +
#'     coord_flip() +
#'     labs(
#'       title = paste0("Non-Zero Coefficients: ", cvd_score),
#'       subtitle = paste0("Dataset: ", dataset_name, " (n = ", nrow(coef_data), " predictors)"),
#'       x = "Predictor",
#'       y = "Absolute Coefficient Value"
#'     ) +
#'     theme_minimal(base_size = 11) +
#'     theme(
#'       plot.title = element_text(face = "bold", size = 14),
#'       plot.subtitle = element_text(size = 11, color = "gray40"),
#'       axis.title = element_text(face = "bold"),
#'       axis.text.y = element_text(size = 8),
#'       legend.position = "top",
#'       panel.grid.major.y = element_blank(),
#'       panel.grid.minor = element_blank()
#'     )
#'
#'   return(p)
#' }
#'
#'
#' #' Create coefficient plots for all score+block models
#' create_coefficient_plots <- function(results_df, output_dir = "coefficient_plots") {
#'   # Create output directory if it doesn't exist
#'   if (!dir.exists(output_dir)) {
#'     dir.create(output_dir, recursive = TRUE)
#'   }
#'
#'   # Extract all non-zero coefficients
#'   cat("\nExtracting non-zero coefficients...\n")
#'   all_coefs <- extract_nonzero_coefs(results_df)
#'
#'   # Get unique model combinations
#'   models <- results_df %>%
#'     select(cvd_score, dataset_name) %>%
#'     distinct()
#'
#'   cat(sprintf("Creating %d coefficient plots...\n", nrow(models)))
#'
#'   plot_list <- list()
#'
#'   for (i in 1:nrow(models)) {
#'     cvd <- models$cvd_score[i]
#'     dataset <- models$dataset_name[i]
#'
#'     # Filter coefficients for this model
#'     model_coefs <- all_coefs %>%
#'       filter(cvd_score == cvd, dataset_name == dataset)
#'
#'     if (nrow(model_coefs) == 0) {
#'       cat(sprintf("  [%d/%d] %s ~ %s: No non-zero coefficients\n", i, nrow(models), cvd, dataset))
#'       next
#'     }
#'
#'     cat(sprintf("  [%d/%d] %s ~ %s: %d non-zero coefficients\n",
#'                 i, nrow(models), cvd, dataset, nrow(model_coefs)))
#'
#'     # Create plot
#'     p <- plot_nonzero_coefs(model_coefs, cvd, dataset)
#'     plot_list[[i]] <- p
#'
#'     # Save plot
#'     filename <- paste0(output_dir, "/coef_", cvd, "_",
#'                        gsub("[^A-Za-z0-9_]", "_", dataset), ".pdf")
#'
#'     # Adjust plot height based on number of predictors
#'     plot_height <- max(6, min(20, 3 + 0.3 * nrow(model_coefs)))
#'
#'     ggsave(filename, plot = p, width = 10, height = plot_height, device = "pdf")
#'   }
#'
#'   cat(sprintf("\n✓ Saved %d plots to %s/\n", length(plot_list), output_dir))
#'
#'   return(list(plots = plot_list, coefficients = all_coefs))
#' }
#'
#'
#' # Run visualization for score+block models
#' cat("\n==========================================================\n")
#' cat("VISUALIZING NON-ZERO COEFFICIENTS (SCORE+BLOCK MODELS)\n")
#' cat("==========================================================\n")
#'
#' coef_plots_results <- create_coefficient_plots(
#'   results_df = results_score_plus_blocks,
#'   output_dir = "coefficient_plots_score_plus_blocks"
#' )
#'
#' # Optional: Create summary statistics
#' coef_summary <- coef_plots_results$coefficients %>%
#'   group_by(cvd_score, dataset_name) %>%
#'   summarise(
#'     n_nonzero = n(),
#'     n_positive = sum(direction == "Positive"),
#'     n_negative = sum(direction == "Negative"),
#'     mean_abs_coef = mean(abs_coefficient),
#'     max_abs_coef = max(abs_coefficient),
#'     .groups = "drop"
#'   ) %>%
#'   arrange(cvd_score, desc(n_nonzero))
#'
#' print(coef_summary)
#'
# write.xlsx(coef_summary, "coefficient_summary.xlsx")
#'
#'
#' # ============================================
#' # 11b. BLOCK-LEVEL PLOTS (FIXED TRUNCATION)
#' # ============================================
#'
#' #' Smart truncation: keep start and end of long names
#' smart_truncate <- function(x, max_len = 30) {
#'   ifelse(
#'     nchar(x) > max_len,
#'     paste0(substr(x, 1, 12), "..", substr(x, nchar(x) - 12, nchar(x))),
#'     x
#'   )
#' }
#'
#' #' Create faceted plot for one predictor block across all CVD scores
#' plot_block_coefficients <- function(coef_data, block_name) {
#'
#'   # Adaptive text size for dense plots
#'   n_predictors <- n_distinct(coef_data$predictor)
#'   x_text_size <- case_when(
#'     block_name %in% c("all_data", "metabolites") ~ 5,
#'     n_predictors > 30 ~ 6,
#'     TRUE ~ 7
#'   )
#'
#'   # Recode CVD score names for display
#'   coef_data <- coef_data %>%
#'     mutate(
#'       cvd_score_label = recode(cvd_score,
#'                                "ascvd_10y" = "ASCVD",
#'                                "frs_10y" = "Framingham",
#'                                "QRISK3_risk" = "QRISK3",
#'                                "SCORE2_score" = "SCORE2"
#'       )
#'     )
#'
#'   # Prepare data with smart truncation
#'   coef_data <- coef_data %>%
#'     mutate(
#'       direction = ifelse(coefficient > 0, "Positive", "Negative"),
#'       abs_coef = abs(coefficient),
#'       predictor_short = smart_truncate(predictor, max_len = 30),
#'       predictor_id = paste(predictor_short, cvd_score_label, sep = "___")
#'     ) %>%
#'     arrange(cvd_score_label, desc(abs_coef))
#'
#'   # Set factor levels in sorted order
#'   coef_data$predictor_id <- factor(coef_data$predictor_id, levels = unique(coef_data$predictor_id))
#'
#'   # Order facets consistently
#'   coef_data$cvd_score_label <- factor(coef_data$cvd_score_label,
#'                                       levels = c("ASCVD", "Framingham", "QRISK3", "SCORE2"))
#'
#'   ggplot(coef_data, aes(x = predictor_id, y = abs_coef, fill = direction)) +
#'     geom_col(width = 0.7) +
#'     facet_wrap(~ cvd_score_label, scales = "free_x", nrow = 1) +
#'     scale_x_discrete(labels = function(x) gsub("___.*$", "", x)) +
#'     scale_fill_manual(
#'       values = c("Positive" = "#E07B39", "Negative" = "#3B7EA1"),
#'       name = "Effect Direction"
#'     ) +
#'     labs(
#'       title = sprintf("Non-Zero Coefficients by CVD Score: %s", block_name),
#'       x = NULL,
#'       y = "Absolute Coefficient"
#'     ) +
#'     theme_minimal(base_size = 10) +
#'     theme(
#'       plot.title = element_text(face = "bold", hjust = 0.5),
#'       axis.text.x = element_text(size = x_text_size, angle = 45, hjust = 1),
#'       strip.text = element_text(face = "bold"),
#'       legend.position = "bottom",
#'       panel.grid.major.x = element_blank()
#'     )
#' }
#'
#'
#'
#' #' Create separate plots for each predictor block
#' create_block_plots <- function(results_df,
#'                                blocks = c("all_data", "metabolites", "lipids", "fatty_acids", "urine_nmr", "body_composition",
#'                                           "clinical_risk_factors", "sociodemographics_lifestyle"),
#'                                output_dir = "coefficient_plots_by_block") {
#'
#'   dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
#'
#'   all_coefs <- extract_nonzero_coefs(results_df) %>%
#'     mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
#'     filter(block %in% blocks) %>%
#'     filter(cvd_score != "mean_risk")  # Exclude composite score
#'
#'   cat(sprintf("\nCreating plots for %d blocks (excluding mean_risk)...\n", length(blocks)))
#'
#'   for (current_block in blocks) {
#'     block_data <- filter(all_coefs, block == current_block)
#'
#'     if (nrow(block_data) == 0) {
#'       cat(sprintf("  %s: No data - skipping\n", current_block))
#'       next
#'     }
#'
#'     n_pred <- n_distinct(block_data$predictor)
#'     n_scores <- n_distinct(block_data$cvd_score)
#'     cat(sprintf("  %s: %d unique predictors across %d CVD scores\n",
#'                 current_block, n_pred, n_scores))
#'
#'     p <- plot_block_coefficients(block_data, current_block)
#'
#'     filename <- file.path(output_dir, sprintf("coef_block_%s.pdf", current_block))
#'     plot_width <- max(12, min(22, 6 + 0.18 * n_pred))
#'
#'     ggsave(filename, p, width = plot_width, height = 6)
#'   }
#'
#'   cat(sprintf("✓ Saved plots to %s/\n", output_dir))
#' }

#' # Re-run
#' create_block_plots(
#'   results_df = results_score_plus_blocks,
#'   output_dir = "coefficient_plots_by_block"
#' )
#'