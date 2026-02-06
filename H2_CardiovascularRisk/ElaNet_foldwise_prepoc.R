### Elastic Net Regression with Preprocessing Inside Each CV-Fold (No Data Leakage)
### Author: Luisa Delius
### Revised: Auto-incrementing output folders, kNN imputation via recipes, simplified code

# ============================================
# SETUP
# ============================================
library(glmnet)
library(tidyverse)
library(tidymodels)
library(qs2)
library(parallel)
library(flextable)
library(officer)
library(patchwork)
library(openxlsx)

set.seed(42)

# Set working directory
wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files/processed_data"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

# CREATE AUTO-INCREMENTING OUTPUT FOLDER
create_output_folder <- function(base_name = "elastic_net_results") {
  existing <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
  pattern <- paste0("^", base_name, "_run(\\d+)$")
  matches <- grep(pattern, existing, value = TRUE)
  
  if (length(matches) == 0) {
    run_number <- 1
  } else {
    run_numbers <- as.integer(gsub(pattern, "\\1", matches))
    run_number <- max(run_numbers) + 1
  }
  
  folder_name <- sprintf("%s_run%03d", base_name, run_number)
  dir.create(folder_name, showWarnings = FALSE)
  cat(sprintf("✓ Created output folder: %s\n", folder_name))
  return(folder_name)
}

OUTPUT_DIR <- create_output_folder()

save_output <- function(filename) {
  file.path(OUTPUT_DIR, filename)
}

# Parallel processing setup
CPUS <- parallel::detectCores() - 3 # change t0 -1 over night
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

# Predictor block datasets
df_body_comp_raw <- df_body_composition_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")

df_demographics_raw <- df_REDcap_demographics_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")

df_risk_factors_raw <- df_risk_factors_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")

df_lipids_raw <- df_lipidomics_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")

df_fatty_Acids_raw <- df_fatty_acids_ElaNet %>%
  select(-Statins, -Supplements) %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")

df_urine_nmr_raw <- df_urine_nmr_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")

datasets_list <- list(
  all_data = df_model2,
  # metabolites = df_model1,
  body_composition = df_body_comp_raw,
  sociodemographics_lifestyle = df_demographics_raw,
  clinical_risk_factors = df_risk_factors_raw,
  lipids = df_lipids_raw,
  fatty_acids = df_fatty_Acids_raw,
  urine_nmr = df_urine_nmr_raw
)

# Body Composition + Height/Weight
height_weight_vars <- df_QRISK3_input %>%
  select(Sample_ID, Height_cm, Weight_kg)

df_body_comp_extended <- df_body_composition_ElaNet %>%
  full_join(height_weight_vars, by = "Sample_ID") %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")

scores_without_ht_wt <- c("SCORE2_score", "ascvd_10y", "frs_10y")

# ============================================
# 2b. SANITIZE DATA TYPES
# ============================================
sanitize_df <- function(df) {
  for (col in names(df)) {
    if (col == "Sample_ID" || col %in% CVD_scores) next
    
    x <- df[[col]]
    
    if (is.logical(x)) {
      df[[col]] <- as.integer(x)
    } else if (inherits(x, c("Date", "POSIXct", "POSIXlt"))) {
      df[[col]] <- as.numeric(x)
    } else if (!is.numeric(x) && !is.factor(x)) {
      df[[col]] <- as.factor(x)
    }
  }
  df
}

cat("Sanitizing data types...\n")
datasets_list <- lapply(datasets_list, sanitize_df)
score_specific_datasets <- lapply(score_specific_datasets, sanitize_df)
df_body_comp_extended <- sanitize_df(df_body_comp_extended)
cat("✓ All datasets sanitized\n")

# ============================================
# 3. PREPROCESSING FUNCTIONS (FOLD-WISE with recipes kNN)
# ============================================

#' Filter columns by missingness threshold and rows by missingness threshold
filter_by_missingness <- function(df, exclude_cols, col_miss_thresh = 0.4, row_miss_thresh = 0.8) {
  
  pred_cols <- setdiff(names(df), exclude_cols)
  
  # Remove columns with >40% missing
  
  col_miss_rates <- sapply(df[pred_cols], function(x) mean(is.na(x)))
  cols_to_keep <- names(col_miss_rates)[col_miss_rates <= col_miss_thresh]
  cols_removed <- length(pred_cols) - length(cols_to_keep)
  
  if (cols_removed > 0) {
    cat(sprintf("    Removed %d columns with >%.0f%% missing\n", cols_removed, col_miss_thresh * 100))
  }
  
  # Remove rows with >80% missing
  if (length(cols_to_keep) > 0) {
    row_miss_rates <- rowMeans(is.na(df[cols_to_keep]))
    rows_to_keep <- which(row_miss_rates <= row_miss_thresh)
    rows_removed <- nrow(df) - length(rows_to_keep)
    
    if (rows_removed > 0) {
      cat(sprintf("    Removed %d rows with >%.0f%% missing\n", rows_removed, row_miss_thresh * 100))
    }
  } else {
    rows_to_keep <- 1:nrow(df)
  }
  
  df[rows_to_keep, c(intersect(exclude_cols, names(df)), cols_to_keep), drop = FALSE]
}


#' Preprocess training and test data using recipes (kNN imputation)
preprocess_train_test <- function(df_train, df_test, exclude_cols = NULL, k = 5) {
  
  # Identify predictor columns
  pred_cols <- setdiff(names(df_train), exclude_cols)
  
  # Create a recipe on training data
  rec <- recipe(~ ., data = df_train[, pred_cols, drop = FALSE]) %>%
    # kNN imputation for all predictors
    step_impute_knn(all_predictors(), neighbors = k) %>%
    # Convert character to factor
    step_string2factor(all_nominal_predictors()) %>%
    # Create dummy variables (one_hot = FALSE means k-1 dummies for k levels)
    step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
    # Center and scale numeric predictors
    step_normalize(all_numeric_predictors())
  
  # Prep on training data
  rec_prepped <- prep(rec, training = df_train[, pred_cols, drop = FALSE])
  
  # Bake (apply) to both train and test
  X_train <- bake(rec_prepped, new_data = df_train[, pred_cols, drop = FALSE]) %>%
    as.matrix()
  
  X_test <- bake(rec_prepped, new_data = df_test[, pred_cols, drop = FALSE]) %>%
    as.matrix()
  
  # Ensure same columns (handle any edge cases)
  common_cols <- intersect(colnames(X_train), colnames(X_test))
  X_train <- X_train[, common_cols, drop = FALSE]
  X_test <- X_test[, common_cols, drop = FALSE]
  
  list(
    X_train = X_train,
    X_test = X_test,
    recipe = rec_prepped,
    feature_names = common_cols
  )
}

# ============================================
# 4. PERMUTATION TEST FUNCTION
# ============================================

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
    
    if (is.null(preprocessed) || ncol(preprocessed$X_train) == 0) return(NA_real_)
    
    best_cv_fit <- NULL
    best_mse <- Inf
    
    for (a in alpha_grid) {
      cv_fit <- tryCatch({
        cv.glmnet(x = preprocessed$X_train, y = y_train, alpha = a,
                  nfolds = min(5, nfolds - 1),
                  type.measure = "mse", standardize = FALSE, family = "gaussian")
      }, error = function(e) NULL)
      
      if (!is.null(cv_fit) && min(cv_fit$cvm) < best_mse) {
        best_mse <- min(cv_fit$cvm)
        best_cv_fit <- cv_fit
      }
    }
    
    if (is.null(best_cv_fit)) return(NA_real_)
    
    y_oof[idx_test] <- as.numeric(
      predict(best_cv_fit, newx = preprocessed$X_test, s = "lambda.min")
    )
  }
  
  ss_res <- sum((y_shuffled - y_oof)^2)
  ss_tot <- sum((y_shuffled - mean(y_shuffled))^2)
  1 - ss_res / ss_tot
}


run_permutation_test <- function(df_raw, y, exclude_cols, observed_Q2,
                                 n_permutations = 100,
                                 alpha_grid = seq(0, 1, by = 0.1),
                                 nfolds = 10, seed = 42) {
  set.seed(seed)
  
  cat(sprintf("  Running permutation test (%d permutations) on %d cores...\n", 
              n_permutations, CPUS))
  
  shuffled_y_list <- lapply(1:n_permutations, function(p) sample(y))
  
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
  
  null_Q2_values <- unlist(null_Q2_values)
  null_Q2_values <- null_Q2_values[!is.na(null_Q2_values)]
  
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
# 6. CORE ELASTIC NET FUNCTION
# ============================================

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


run_elastic_net_model <- function(df_raw, outcome_col, cvd_score_name, dataset_name,
                                  exclude_cols = "Sample_ID",
                                  alpha_grid = seq(0, 1, by = 0.1),
                                  nfolds = 10, 
                                  n_permutations = 100,
                                  col_miss_thresh = 0.4,
                                  row_miss_thresh = 0.8,
                                  seed = 42) {
  set.seed(seed)
  
  # Remove rows with missing outcome
  df_complete <- df_raw %>% filter(!is.na(.data[[outcome_col]]))
  y <- df_complete[[outcome_col]]
  
  # Columns to exclude from predictors
  all_exclude <- unique(c(exclude_cols, outcome_col, CVD_scores))
  
  # Apply missingness filtering
  cat(sprintf("  Filtering by missingness (col >%.0f%%, row >%.0f%%)...\n", 
              col_miss_thresh * 100, row_miss_thresh * 100))
  
  df_filtered <- filter_by_missingness(df_complete, all_exclude, col_miss_thresh, row_miss_thresh)
  
  # Update y to match filtered rows
  y <- df_filtered[[outcome_col]]
  n_obs <- nrow(df_filtered)
  
  cat(sprintf("  Sample size after filtering: %d\n", n_obs))
  
  if (n_obs < nfolds * 2) {
    warning(sprintf("Too few observations (%d) for %d-fold CV", n_obs, nfolds))
    return(NULL)
  }
  
  # ---- PART A: Full-data model (for coefficients) ----
  # Use full data preprocessing to get feature names and final model
  preprocessed_full <- preprocess_train_test(df_filtered, df_filtered, all_exclude)
  X_full <- preprocessed_full$X_train
  n_pred <- ncol(X_full)
  predictor_names <- colnames(X_full)
  
  cat(sprintf("  Predictors after processing: %d\n", n_pred))
  
  if (n_pred == 0) {
    warning("No predictors remaining after preprocessing")
    return(NULL)
  }
  
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
  
  # ---- PART B: In-sample predictions (R²_Y) ----
  y_pred_in <- as.numeric(predict(best_cv_fit, newx = X_full, s = "lambda.min"))
  res_in <- y - y_pred_in
  rmse_in <- sqrt(mean(res_in^2))
  mae_in <- mean(abs(res_in))
  R2_Y <- 1 - sum(res_in^2) / sum((y - mean(y))^2)
  cor_in <- cor(y_pred_in, y)
  
  # ---- PART C: Out-of-fold predictions (Q²_Y) ----
  oof_results <- get_oof_predictions_clean(
    df_raw = df_filtered,
    y = y,
    exclude_cols = all_exclude,
    alpha_grid = alpha_grid,
    nfolds = nfolds,
    seed = seed + 100
  )
  
  y_pred_oof <- oof_results$predictions
  alpha_per_fold <- oof_results$alpha_per_fold
  Q2_per_fold <- oof_results$Q2_per_fold
  
  res_oof <- y - y_pred_oof
  rmse_oof <- sqrt(mean(res_oof^2))
  mae_oof <- mean(abs(res_oof))
  Q2_Y <- 1 - sum(res_oof^2) / sum((y - mean(y))^2)
  cor_oof <- cor(y_pred_oof, y)
  
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
  
  dev_explained <- best_cv_fit$glmnet.fit$dev.ratio[
    which.min(abs(best_cv_fit$glmnet.fit$lambda - lambda_min))
  ]
  
  # ---- PART E: Permutation test ----
  if (n_permutations > 0) {
    perm_results <- run_permutation_test(
      df_raw = df_filtered,
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
    cv_mse_min = cv_mse_min,
    cv_mse_1se = cv_mse_1se,
    R2_Y = R2_Y,
    in_sample_rmse = rmse_in,
    in_sample_mae = mae_in,
    in_sample_cor = cor_in,
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
    Q2_fold_mean = Q2_fold_mean,
    Q2_fold_sd = Q2_fold_sd,
    Q2_fold_summary = sprintf("%.3f ± %.3f", Q2_fold_mean, Q2_fold_sd),
    Q2_per_fold = list(Q2_per_fold),
    Q2_R2_ratio = Q2_R2_ratio,
    R2_Q2_gap = R2_Q2_gap,
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
    col_miss_threshold = col_miss_thresh,
    row_miss_threshold = row_miss_thresh,
    imputation_method = "kNN",
    seed = seed,
    date_run = as.character(Sys.time()),
    cv_fit_object = list(best_cv_fit)
  )
}

# ============================================
# 7. MODEL RUNNER FUNCTIONS
# ============================================

run_score_specific_models <- function(score_specific_datasets,
                                      alpha_grid = seq(0, 1, by = 0.1),
                                      nfolds = 10, 
                                      n_permutations = 100,
                                      seed = 42,
                                      output_file = "elastic_net_score_specific.qs2") {
  results_list <- list()
  n_models <- length(score_specific_datasets)
  
  cat("==========================================================\n")
  cat("ELASTIC NET: SCORE-SPECIFIC DATASETS\n")
  cat(sprintf("Total models: %d | Permutations: %d\n", n_models, n_permutations))
  cat("==========================================================\n")
  
  for (i in seq_along(score_specific_datasets)) {
    cvd_score <- names(score_specific_datasets)[i]
    df <- score_specific_datasets[[i]]
    dataset_name <- paste0(cvd_score, "_specific")
    
    cat(sprintf("\n[%d/%d] %s ~ %s\n", i, n_models, cvd_score, dataset_name))
    
    tryCatch({
      results_list[[i]] <- run_elastic_net_model(
        df_raw = df,
        outcome_col = cvd_score,
        cvd_score_name = cvd_score,
        dataset_name = dataset_name,
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
  qs_save(results_df, file = save_output(output_file))
  cat(sprintf("\n✓ Saved %d models to %s\n", nrow(results_df), output_file))
  results_df
}


run_common_datasets_models <- function(datasets_list, cvd_score_names,
                                       alpha_grid = seq(0, 1, by = 0.1),
                                       nfolds = 10, 
                                       n_permutations = 100,
                                       seed = 42,
                                       output_file = "elastic_net_common.qs2") {
  results_list <- list()
  counter <- 1
  n_models <- length(datasets_list) * length(cvd_score_names)
  
  cat("==========================================================\n")
  cat("ELASTIC NET: COMMON DATASETS × CVD SCORES\n")
  cat(sprintf("Total models: %d | Permutations: %d\n", n_models, n_permutations))
  cat("==========================================================\n")
  
  for (dataset_name in names(datasets_list)) {
    df <- datasets_list[[dataset_name]]
    
    for (cvd_score in cvd_score_names) {
      cat(sprintf("\n[%d/%d] %s ~ %s\n", counter, n_models, cvd_score, dataset_name))
      
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
  qs_save(results_df, file = save_output(output_file))
  cat(sprintf("\n✓ Saved %d models to %s\n", nrow(results_df), output_file))
  results_df
}


run_score_specific_plus_blocks_models <- function(score_specific_datasets, datasets_list,
                                                  alpha_grid = seq(0, 1, by = 0.1),
                                                  nfolds = 10, 
                                                  n_permutations = 100,
                                                  seed = 42,
                                                  output_file = "elastic_net_score_plus_blocks.qs2") {
  results_list <- list()
  counter <- 1
  n_models <- length(score_specific_datasets) * length(datasets_list)
  
  cat("==========================================================\n")
  cat("ELASTIC NET: SCORE-SPECIFIC + PREDICTOR BLOCK\n")
  cat(sprintf("Total models: %d | Permutations: %d\n", n_models, n_permutations))
  cat("==========================================================\n")
  
  for (cvd_score in names(score_specific_datasets)) {
    for (block_name in names(datasets_list)) {
      dataset_name <- paste0(cvd_score, "_specific+", block_name)
      
      cat(sprintf("\n[%d/%d] %s ~ %s\n", counter, n_models, cvd_score, dataset_name))
      
      df_combined <- build_score_plus_block(
        score = cvd_score, block_name = block_name,
        score_specific_datasets = score_specific_datasets,
        datasets_list = datasets_list
      )
      
      tryCatch({
        results_list[[counter]] <- run_elastic_net_model(
          df_raw = df_combined,
          outcome_col = cvd_score,
          cvd_score_name = cvd_score,
          dataset_name = dataset_name,
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
  qs_save(results_df, file = save_output(output_file))
  cat(sprintf("\n✓ Saved %d models to %s\n", nrow(results_df), output_file))
  results_df
}

# ============================================
# 8. RESULTS VIEWER
# ============================================
view_results_summary <- function(results_df) {
  cat("\n=== RESULTS SUMMARY ===\n\n")
  
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

results_score_specific <- run_score_specific_models(
  score_specific_datasets = score_specific_datasets,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 5, 
  n_permutations = 1000,
  seed = 42,
  output_file = "elastic_net_score_specific_5fold.qs2"
)

results_common <- run_common_datasets_models(
  datasets_list = datasets_list,
  cvd_score_names = CVD_scores,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 5, 
  n_permutations = 1000,
  seed = 42,
  output_file = "elastic_net_common_5fold.qs2"
)

results_score_plus_blocks <- run_score_specific_plus_blocks_models(
  score_specific_datasets = score_specific_datasets,
  datasets_list = datasets_list,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 5, 
  n_permutations = 1000,
  seed = 42,
  output_file = "elastic_net_score_plus_blocks_5fold.qs2"
)

# Body Composition + Height/Weight Models
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
      alpha_grid = seq(0, 1, by = 0.1),
      nfolds = 5,
      n_permutations = 1000,
      seed = 42
    )
    cat("  ✓ Success\n")
  }, error = function(e) {
    cat(sprintf("  ✗ ERROR: %s\n", e$message))
  })
}

results_body_comp_extended_df <- bind_rows(results_body_comp_extended)
qs_save(results_body_comp_extended_df, save_output("elastic_net_body_comp_extended_5fold.qs2"))

# Score-specific + body composition with height/weight
results_score_plus_body_extended <- list()

for (i in seq_along(scores_without_ht_wt)) {
  score <- scores_without_ht_wt[i]
  dataset_name <- paste0(score, "_specific+body_comp_with_ht_wt")
  
  cat(sprintf("\n[%d/3] %s ~ %s\n", i, score, dataset_name))
  
  base_df <- score_specific_datasets[[score]]
  body_comp_predictors <- df_body_comp_extended %>%
    select(-all_of(CVD_scores))
  
  df_combined <- base_df %>%
    left_join(body_comp_predictors, by = "Sample_ID") %>%
    filter(!is.na(.data[[score]]))
  
  tryCatch({
    results_score_plus_body_extended[[i]] <- run_elastic_net_model(
      df_raw = df_combined,
      outcome_col = score,
      cvd_score_name = score,
      dataset_name = dataset_name,
      alpha_grid = seq(0, 1, by = 0.1),
      nfolds = 5,
      n_permutations = 1000,
      seed = 42
    )
    cat("  ✓ Success\n")
  }, error = function(e) {
    cat(sprintf("  ✗ ERROR: %s\n", e$message))
  })
}

results_score_plus_body_extended_df <- bind_rows(results_score_plus_body_extended)
qs_save(results_score_plus_body_extended_df, save_output("elastic_net_score_plus_body_extended_5fold.qs2"))

# Combine all results
results_all <- bind_rows(
  results_score_specific, 
  results_common, 
  results_score_plus_blocks,
  results_body_comp_extended_df, 
  results_score_plus_body_extended_df
)

view_results_summary(results_all)

# ============================================
# 10. EXPORT TO EXCEL
# ============================================
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
    imputation_method, col_miss_threshold, row_miss_threshold,
    n_folds, seed, date_run
  )

# Top predictors sheet
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

# Q²_Y overview
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
writeData(wb, "Q2_Y_Overview", "Table A: Q²_Y (all models)", startRow = 1, startCol = 1)
writeData(wb, "Q2_Y_Overview", table_A_all_models, startRow = 2, startCol = 1)

start_row_B <- nrow(table_A_all_models) + 5
writeData(wb, "Q2_Y_Overview", "Table B: Δ Q²_Y = (score-specific + block) − (score-specific baseline)",
          startRow = start_row_B, startCol = 1)
writeData(wb, "Q2_Y_Overview", table_B_delta, startRow = start_row_B + 1, startCol = 1)

addStyle(wb, "Q2_Y_Overview", headerStyle,
         rows = c(2, start_row_B + 1),
         cols = 1:max(ncol(table_A_all_models), ncol(table_B_delta)),
         gridExpand = TRUE)
setColWidths(wb, "Q2_Y_Overview",
             cols = 1:max(ncol(table_A_all_models), ncol(table_B_delta)),
             widths = "auto")

saveWorkbook(wb, save_output("elastic_net_results.xlsx"), overwrite = TRUE)

cat("\n✓ Excel saved: elastic_net_results.xlsx\n")

# ============================================
# 11. VISUALIZE NON-ZERO COEFFICIENTS
# ============================================

#' Extract all non-zero coefficients from model results
extract_nonzero_coefs <- function(results_df) {
  coef_list <- list()
  
  for (i in 1:nrow(results_df)) {
    predictor_names <- results_df$all_predictor_names[[i]]
    coefficients <- results_df$all_coefficients[[i]]
    
    # Get non-zero coefficients
    nonzero_mask <- coefficients != 0
    
    if (sum(nonzero_mask) > 0) {
      coef_list[[i]] <- tibble(
        cvd_score = results_df$cvd_score[i],
        dataset_name = results_df$dataset_name[i],
        predictor = predictor_names[nonzero_mask],
        coefficient = coefficients[nonzero_mask],
        abs_coefficient = abs(coefficients[nonzero_mask]),
        direction = ifelse(coefficients[nonzero_mask] > 0, "Positive", "Negative")
      )
    }
  }
  
  bind_rows(coef_list)
}


#' Create bar plot for a single model's non-zero coefficients
plot_nonzero_coefs <- function(coef_data, cvd_score, dataset_name) {
  # Order predictors by absolute coefficient value
  coef_data <- coef_data %>%
    arrange(abs_coefficient) %>%
    mutate(predictor = factor(predictor, levels = predictor))
  
  # Create plot
  p <- ggplot(coef_data, aes(x = predictor, y = abs_coefficient, fill = direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = c("Positive" = "#FF8C00", "Negative" = "#4169E1"),
                      name = "Direction") +
    coord_flip() +
    labs(
      title = paste0("Non-Zero Coefficients: ", cvd_score),
      subtitle = paste0("Dataset: ", dataset_name, " (n = ", nrow(coef_data), " predictors)"),
      x = "Predictor",
      y = "Absolute Coefficient Value"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray40"),
      axis.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 8),
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}


#' Create coefficient plots for all score+block models
create_coefficient_plots <- function(results_df, output_dir = "coefficient_plots") {
  # Create output directory
  full_output_dir <- save_output(output_dir)
  if (!dir.exists(full_output_dir)) {
    dir.create(full_output_dir, recursive = TRUE)
  }
  
  # Extract all non-zero coefficients
  cat("\nExtracting non-zero coefficients...\n")
  all_coefs <- extract_nonzero_coefs(results_df)
  
  # Get unique model combinations
  models <- results_df %>%
    select(cvd_score, dataset_name) %>%
    distinct()
  
  cat(sprintf("Creating %d coefficient plots...\n", nrow(models)))
  
  plot_list <- list()
  
  for (i in 1:nrow(models)) {
    cvd <- models$cvd_score[i]
    dataset <- models$dataset_name[i]
    
    # Filter coefficients for this model
    model_coefs <- all_coefs %>%
      filter(cvd_score == cvd, dataset_name == dataset)
    
    if (nrow(model_coefs) == 0) {
      cat(sprintf("  [%d/%d] %s ~ %s: No non-zero coefficients\n", i, nrow(models), cvd, dataset))
      next
    }
    
    cat(sprintf("  [%d/%d] %s ~ %s: %d non-zero coefficients\n",
                i, nrow(models), cvd, dataset, nrow(model_coefs)))
    
    # Create plot
    p <- plot_nonzero_coefs(model_coefs, cvd, dataset)
    plot_list[[i]] <- p
    
    # Save plot
    filename <- file.path(full_output_dir, paste0("coef_", cvd, "_",
                                                  gsub("[^A-Za-z0-9_]", "_", dataset), ".pdf"))
    
    # Adjust plot height based on number of predictors
    plot_height <- max(6, min(20, 3 + 0.3 * nrow(model_coefs)))
    
    ggsave(filename, plot = p, width = 10, height = plot_height, device = "pdf")
  }
  
  cat(sprintf("\n✓ Saved %d plots to %s/\n", length(plot_list[!sapply(plot_list, is.null)]), output_dir))
  
  return(list(plots = plot_list, coefficients = all_coefs))
}


#' Smart truncation: keep start and end of long names
smart_truncate <- function(x, max_len = 30) {
  ifelse(
    nchar(x) > max_len,
    paste0(substr(x, 1, 12), "..", substr(x, nchar(x) - 12, nchar(x))),
    x
  )
}


#' Create faceted plot for one predictor block across all CVD scores
plot_block_coefficients <- function(coef_data, block_name) {
  
  # Adaptive text size for dense plots
  n_predictors <- n_distinct(coef_data$predictor)
  x_text_size <- case_when(
    block_name %in% c("all_data", "metabolites") ~ 5,
    n_predictors > 30 ~ 6,
    TRUE ~ 7
  )
  
  # Recode CVD score names for display
  coef_data <- coef_data %>%
    mutate(
      cvd_score_label = recode(cvd_score,
                               "ascvd_10y" = "ASCVD",
                               "frs_10y" = "Framingham",
                               "QRISK3_risk" = "QRISK3",
                               "SCORE2_score" = "SCORE2",
                               "mean_risk" = "Composite"
      )
    )
  
  # Prepare data with smart truncation
  coef_data <- coef_data %>%
    mutate(
      direction = ifelse(coefficient > 0, "Positive", "Negative"),
      abs_coef = abs(coefficient),
      predictor_short = smart_truncate(predictor, max_len = 30),
      predictor_id = paste(predictor_short, cvd_score_label, sep = "___")
    ) %>%
    arrange(cvd_score_label, desc(abs_coef))
  
  # Set factor levels in sorted order
  coef_data$predictor_id <- factor(coef_data$predictor_id, levels = unique(coef_data$predictor_id))
  
  # Order facets consistently
  coef_data$cvd_score_label <- factor(coef_data$cvd_score_label,
                                      levels = c("ASCVD", "Framingham", "QRISK3", "SCORE2", "Composite"))
  
  ggplot(coef_data, aes(x = predictor_id, y = abs_coef, fill = direction)) +
    geom_col(width = 0.7) +
    facet_wrap(~ cvd_score_label, scales = "free_x", nrow = 1) +
    scale_x_discrete(labels = function(x) gsub("___.*$", "", x)) +
    scale_fill_manual(
      values = c("Positive" = "#E07B39", "Negative" = "#3B7EA1"),
      name = "Effect Direction"
    ) +
    labs(
      title = sprintf("Non-Zero Coefficients by CVD Score: %s", block_name),
      x = NULL,
      y = "Absolute Coefficient"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = x_text_size, angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    )
}


#' Create separate plots for each predictor block
create_block_plots <- function(results_df,
                               blocks = c("all_data", "metabolites", "lipids", "fatty_acids", "urine_nmr", 
                                          "body_composition", "clinical_risk_factors", "sociodemographics_lifestyle"),
                               output_dir = "coefficient_plots_by_block") {
  
  full_output_dir <- save_output(output_dir)
  dir.create(full_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  all_coefs <- extract_nonzero_coefs(results_df) %>%
    mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
    filter(block %in% blocks) %>%
    filter(cvd_score != "mean_risk")  # Exclude composite score
  
  cat(sprintf("\nCreating plots for %d blocks (excluding mean_risk)...\n", length(blocks)))
  
  for (current_block in blocks) {
    block_data <- filter(all_coefs, block == current_block)
    
    if (nrow(block_data) == 0) {
      cat(sprintf("  %s: No data - skipping\n", current_block))
      next
    }
    
    n_pred <- n_distinct(block_data$predictor)
    n_scores <- n_distinct(block_data$cvd_score)
    cat(sprintf("  %s: %d unique predictors across %d CVD scores\n",
                current_block, n_pred, n_scores))
    
    p <- plot_block_coefficients(block_data, current_block)
    
    filename <- file.path(full_output_dir, sprintf("coef_block_%s.pdf", current_block))
    plot_width <- max(12, min(22, 6 + 0.18 * n_pred))
    
    ggsave(filename, p, width = plot_width, height = 6)
  }
  
  cat(sprintf("✓ Saved plots to %s/\n", output_dir))
}


# Run coefficient visualization
cat("\n==========================================================\n")
cat("VISUALIZING NON-ZERO COEFFICIENTS\n")
cat("==========================================================\n")

coef_plots_results <- create_coefficient_plots(
  results_df = results_score_plus_blocks,
  output_dir = "coefficient_plots_score_plus_blocks"
)

# Coefficient summary statistics
coef_summary <- coef_plots_results$coefficients %>%
  group_by(cvd_score, dataset_name) %>%
  summarise(
    n_nonzero = n(),
    n_positive = sum(direction == "Positive"),
    n_negative = sum(direction == "Negative"),
    mean_abs_coef = mean(abs_coefficient),
    max_abs_coef = max(abs_coefficient),
    .groups = "drop"
  ) %>%
  arrange(cvd_score, desc(n_nonzero))

print(coef_summary)

# Save coefficient summary
write.xlsx(coef_summary, save_output("coefficient_summary.xlsx"))

# Create block-level plots
create_block_plots(
  results_df = results_score_plus_blocks,
  output_dir = "coefficient_plots_by_block"
)


# ============================================
# 12. TABLE CREATION FUNCTIONS
# ============================================

cvd_score_order <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y", "mean_risk")

cvd_score_labels <- c(
  "QRISK3_risk" = "QRISK3",
  "SCORE2_score" = "SCORE2", 
  "frs_10y" = "Framingham",
  "ascvd_10y" = "ASCVD",
  "mean_risk" = "Composite"
)

dataset_labels_table <- c(
  "all_data" = "All predictors",
  "metabolites" = "All metabolites",
  "body_composition" = "Body composition",
  "sociodemographics_lifestyle" = "Sociodemographics & lifestyle",
  "clinical_risk_factors" = "Clinical risk factors",
  "lipids" = "Lipids",
  "fatty_acids" = "Fatty acids",
  "urine_nmr" = "Urine NMR",
  "body_comp_with_ht_wt" = "Body composition + Ht/Wt"
)

create_score_specific_table <- function(results_df, table_title) {
  
  table_data <- results_df %>%
    filter(cvd_score %in% cvd_score_order) %>%
    mutate(
      cvd_score = factor(cvd_score, levels = cvd_score_order),
      CVD_Score = cvd_score_labels[as.character(cvd_score)]
    ) %>%
    arrange(cvd_score)
  
  has_mae_sd <- "MAE_fold_sd" %in% names(table_data) && any(!is.na(table_data$MAE_fold_sd))
  has_rmse_sd <- "RMSE_fold_sd" %in% names(table_data) && any(!is.na(table_data$RMSE_fold_sd))
  
  formatted_table <- table_data %>%
    mutate(
      `R²` = sprintf("%.3f", R2_Y),
      `Q²` = ifelse(!is.na(Q2_fold_sd) & Q2_fold_sd > 0,
                    sprintf("%.3f ± %.3f", Q2_Y, Q2_fold_sd),
                    sprintf("%.3f", Q2_Y)),
      MAE = if (has_mae_sd) {
        ifelse(!is.na(MAE_fold_sd) & MAE_fold_sd > 0,
               sprintf("%.3f ± %.3f", Q2_Y_mae, MAE_fold_sd),
               sprintf("%.3f", Q2_Y_mae))
      } else {
        sprintf("%.3f", Q2_Y_mae)
      },
      RMSE = if (has_rmse_sd) {
        ifelse(!is.na(RMSE_fold_sd) & RMSE_fold_sd > 0,
               sprintf("%.3f ± %.3f", Q2_Y_rmse, RMSE_fold_sd),
               sprintf("%.3f", Q2_Y_rmse))
      } else {
        sprintf("%.3f", Q2_Y_rmse)
      },
      `Permutation p-value` = ifelse(is.na(permutation_p_value), 
                                     "—",
                                     ifelse(permutation_p_value < 0.001, 
                                            "<0.001",
                                            sprintf("%.3f", permutation_p_value)))
    ) %>%
    select(CVD_Score, `R²`, `Q²`, MAE, RMSE, `Permutation p-value`)
  
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
    width(j = 6, width = 1.5) %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1.5), part = "header") %>%
    padding(padding = 3, part = "all")
  
  return(ft)
}


create_multi_score_table <- function(results_df, table_title, 
                                     dataset_order = NULL,
                                     dataset_labels = dataset_labels_table) {
  
  table_data <- results_df %>%
    filter(cvd_score %in% cvd_score_order) %>%
    mutate(
      cvd_score = factor(cvd_score, levels = cvd_score_order),
      CVD_Score = cvd_score_labels[as.character(cvd_score)]
    )
  
  if (!is.null(dataset_order)) {
    table_data <- table_data %>%
      filter(dataset_name %in% dataset_order) %>%
      mutate(dataset_name = factor(dataset_name, levels = dataset_order))
  }
  
  if (nrow(table_data) == 0) {
    warning("No data remaining after filtering!")
    return(NULL)
  }
  
  table_data <- table_data %>% arrange(dataset_name, cvd_score)
  
  has_mae_sd <- "MAE_fold_sd" %in% names(table_data) && any(!is.na(table_data$MAE_fold_sd))
  has_rmse_sd <- "RMSE_fold_sd" %in% names(table_data) && any(!is.na(table_data$RMSE_fold_sd))
  
  formatted_table <- table_data %>%
    mutate(
      Dataset = ifelse(dataset_name %in% names(dataset_labels),
                       dataset_labels[as.character(dataset_name)],
                       as.character(dataset_name)),
      `R²` = sprintf("%.3f", R2_Y),
      `Q²` = ifelse(!is.na(Q2_fold_sd) & Q2_fold_sd > 0,
                    sprintf("%.3f ± %.3f", Q2_Y, Q2_fold_sd),
                    sprintf("%.3f", Q2_Y)),
      MAE = if (has_mae_sd) {
        ifelse(!is.na(MAE_fold_sd) & MAE_fold_sd > 0,
               sprintf("%.3f ± %.3f", Q2_Y_mae, MAE_fold_sd),
               sprintf("%.3f", Q2_Y_mae))
      } else {
        sprintf("%.3f", Q2_Y_mae)
      },
      RMSE = if (has_rmse_sd) {
        ifelse(!is.na(RMSE_fold_sd) & RMSE_fold_sd > 0,
               sprintf("%.3f ± %.3f", Q2_Y_rmse, RMSE_fold_sd),
               sprintf("%.3f", Q2_Y_rmse))
      } else {
        sprintf("%.3f", Q2_Y_rmse)
      },
      `Permutation p-value` = ifelse(is.na(permutation_p_value), 
                                     "—",
                                     ifelse(permutation_p_value < 0.001, 
                                            "<0.001",
                                            sprintf("%.3f", permutation_p_value)))
    ) %>%
    select(Dataset, CVD_Score, `R²`, `Q²`, MAE, RMSE, `Permutation p-value`)
  
  n_rows <- nrow(formatted_table)
  
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
    width(j = 7, width = 1.5) %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1.5), part = "header") %>%
    padding(padding = 3, part = "all")
  
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
# CREATE TABLES
# ============================================

cat("\n==========================================================\n")
cat("CREATING PUBLICATION TABLES\n")
cat("==========================================================\n")

# TABLE 1: Score-specific only
cat("\nTable 1: Score-specific models...\n")
table1 <- create_score_specific_table(
  results_df = results_score_specific,
  table_title = "Table 1: Elastic Net Performance – Score-Specific Input Variables"
)
cat("  ✓ Created\n")

# TABLE 2: Common datasets
cat("\nTable 2: Common datasets...\n")
datasets_table2 <- c("all_data", "metabolites", "lipids", "fatty_acids", "urine_nmr", 
                     "body_composition", "sociodemographics_lifestyle", "clinical_risk_factors")
datasets_table2 <- datasets_table2[datasets_table2 %in% unique(results_common$dataset_name)]

table2 <- create_multi_score_table(
  results_df = results_common %>% filter(dataset_name %in% datasets_table2),
  table_title = "Table 2: Elastic Net Performance – Predictor Blocks",
  dataset_order = datasets_table2
)
cat(sprintf("  ✓ Created (%d datasets)\n", length(datasets_table2)))

# TABLE 3: Score-specific + block
cat("\nTable 3: Score-specific + predictor blocks...\n")

blocks_table3 <- c("all_data", "metabolites", "lipids", "fatty_acids", "urine_nmr", 
                   "body_composition", "sociodemographics_lifestyle", "clinical_risk_factors")

dataset_order_table3 <- results_score_plus_blocks %>%
  mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
  filter(block %in% blocks_table3) %>%
  arrange(match(block, blocks_table3)) %>%
  pull(dataset_name) %>%
  unique()

dataset_labels_table3 <- dataset_labels_table
for (ds in dataset_order_table3) {
  block <- sub("^.*_specific\\+", "", ds)
  block_label <- ifelse(block %in% names(dataset_labels_table), 
                        dataset_labels_table[block], block)
  dataset_labels_table3[ds] <- paste0("Score inputs + ", block_label)
}

table3 <- create_multi_score_table(
  results_df = results_score_plus_blocks %>%
    mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
    filter(block %in% blocks_table3),
  table_title = "Table 3: Elastic Net Performance – Score-Specific + Predictor Blocks",
  dataset_order = dataset_order_table3,
  dataset_labels = dataset_labels_table3
)
cat(sprintf("  ✓ Created (%d dataset combinations)\n", length(dataset_order_table3)))

# Save tables
cat("\nSaving tables...\n")

save_as_docx(table1, path = save_output("Table1_score_specific.docx"))
save_as_docx(table2, path = save_output("Table2_predictor_blocks.docx"))
save_as_docx(table3, path = save_output("Table3_score_plus_blocks.docx"))

save_as_html(table1, path = save_output("Table1_score_specific.html"))
save_as_html(table2, path = save_output("Table2_predictor_blocks.html"))
save_as_html(table3, path = save_output("Table3_score_plus_blocks.html"))

# Combined document
doc <- read_docx() %>%
  body_add_par("Elastic Net Regression Results", style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(value = table1) %>%
  body_add_break() %>%
  body_add_flextable(value = table2) %>%
  body_add_break() %>%
  body_add_flextable(value = table3)

print(doc, target = save_output("Elastic_Net_All_Tables.docx"))

cat("\n✓ Tables saved\n")

# ============================================
# 13. VISUALIZATION
# ============================================

cat("\n==========================================================\n")
cat("CREATING VISUALIZATION PLOTS\n")
cat("==========================================================\n")

cvd_colors <- c(
  "QRISK3" = "#E69F00",
  "SCORE2" = "#56B4E9",
  "Framingham" = "#009E73",
  "ASCVD" = "#F0E442",
  "Composite" = "#D55E00"
)

# Prepare data for left column
score_specific_grouped <- results_score_specific %>%
  filter(cvd_score %in% cvd_score_order) %>%
  mutate(dataset_name = "score_specific", dataset_label = "Score-Specific")

datasets_for_left <- c("metabolites", "lipids", "fatty_acids", "urine_nmr", 
                       "body_composition", "clinical_risk_factors", 
                       "sociodemographics_lifestyle")
datasets_for_left <- datasets_for_left[datasets_for_left %in% unique(results_common$dataset_name)]

dataset_labels_left <- c("score_specific" = "Score-Specific", dataset_labels_table[datasets_for_left])
dataset_order_left <- rev(c("score_specific", datasets_for_left))

common_blocks_all_scores <- results_common %>%
  filter(dataset_name %in% datasets_for_left, cvd_score %in% cvd_score_order)

plot_data_left <- bind_rows(score_specific_grouped, common_blocks_all_scores) %>%
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

# Prepare data for right column
baseline_grouped <- results_score_specific %>%
  filter(cvd_score %in% cvd_score_order) %>%
  mutate(dataset_name = "score_specific", dataset_label = "Score-Specific", block = "score_specific")

blocks_for_plus <- datasets_for_left

score_plus_data <- results_score_plus_blocks %>%
  mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
  filter(block %in% blocks_for_plus, cvd_score %in% cvd_score_order) %>%
  mutate(
    dataset_name = block,
    dataset_label = ifelse(block %in% names(dataset_labels_table), dataset_labels_table[block], block)
  )

dataset_order_right <- rev(c("score_specific", blocks_for_plus))
dataset_labels_right <- c("score_specific" = "Score-Specific", dataset_labels_table[blocks_for_plus])

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

# Handle missing SD columns
if (!"MAE_fold_sd" %in% names(plot_data_left)) plot_data_left$MAE_fold_sd <- 0
if (!"MAE_fold_sd" %in% names(plot_data_right)) plot_data_right$MAE_fold_sd <- 0

# Create plots
p1_left_q2 <- plot_data_left %>%
  ggplot(aes(x = Q2_Y, y = dataset_label, fill = cvd_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(
    aes(xmin = pmax(0, Q2_Y - Q2_fold_sd), xmax = Q2_Y + Q2_fold_sd),
    position = position_dodge(width = 0.8), height = 0.3, linewidth = 0.4
  ) +
  geom_text(
    aes(x = Q2_Y + Q2_fold_sd, label = sprintf("%.2f", Q2_Y)),
    position = position_dodge(width = 0.8), hjust = -0.2, size = 2.2
  ) +
  scale_fill_manual(values = cvd_colors) +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Score-Specific & Predictor Blocks", x = expression(Q^2), y = NULL, fill = "CVD Score") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = element_text(size = 8), axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()
  )

p2_right_q2 <- plot_data_right %>%
  ggplot(aes(x = Q2_Y, y = dataset_label, fill = cvd_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(
    aes(xmin = pmax(0, Q2_Y - Q2_fold_sd), xmax = Q2_Y + Q2_fold_sd),
    position = position_dodge(width = 0.8), height = 0.3, linewidth = 0.4
  ) +
  geom_text(
    aes(x = Q2_Y + Q2_fold_sd, label = sprintf("%.2f", Q2_Y)),
    position = position_dodge(width = 0.8), hjust = -0.2, size = 2.2
  ) +
  scale_fill_manual(values = cvd_colors) +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Score-Specific + Predictor Blocks", x = expression(Q^2), y = NULL, fill = "CVD Score") +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = element_text(size = 8), axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")
  )

p3_left_mae <- plot_data_left %>%
  ggplot(aes(x = Q2_Y_mae, y = dataset_label, fill = cvd_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(
    aes(xmin = pmax(0, Q2_Y_mae - MAE_fold_sd), xmax = Q2_Y_mae + MAE_fold_sd),
    position = position_dodge(width = 0.8), height = 0.3, linewidth = 0.4
  ) +
  geom_text(
    aes(x = Q2_Y_mae + MAE_fold_sd, label = sprintf("%.2f", Q2_Y_mae)),
    position = position_dodge(width = 0.8), hjust = -0.2, size = 2.2
  ) +
  scale_fill_manual(values = cvd_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "", x = "MAE", y = NULL, fill = "CVD Score") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8), axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()
  )

p4_right_mae <- plot_data_right %>%
  ggplot(aes(x = Q2_Y_mae, y = dataset_label, fill = cvd_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(
    aes(xmin = pmax(0, Q2_Y_mae - MAE_fold_sd), xmax = Q2_Y_mae + MAE_fold_sd),
    position = position_dodge(width = 0.8), height = 0.3, linewidth = 0.4
  ) +
  geom_text(
    aes(x = Q2_Y_mae + MAE_fold_sd, label = sprintf("%.2f", Q2_Y_mae)),
    position = position_dodge(width = 0.8), hjust = -0.2, size = 2.2
  ) +
  scale_fill_manual(values = cvd_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "", x = "MAE", y = NULL, fill = "CVD Score") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 8), axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")
  )

combined_elastic_plot <- (p1_left_q2 | p2_right_q2) / (p3_left_mae | p4_right_mae) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Elastic Net Performance: Cross-Validated Predictive Ability",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

print(combined_elastic_plot)

ggsave(save_output("elastic_net_performance_comparison.png"),
       plot = combined_elastic_plot, width = 16, height = 12, dpi = 300, bg = "white")

ggsave(save_output("elastic_net_performance_comparison.pdf"),
       plot = combined_elastic_plot, width = 16, height = 12, device = "pdf")

cat("\n✓ Plots saved\n")

# Save combined results
qs_save(results_all, save_output("elastic_net_all_results.qs2"))

cat("\n==========================================================\n")
cat(sprintf("ALL OUTPUT FILES SAVED TO: %s\n", OUTPUT_DIR))
cat("==========================================================\n")

# ============================================
# 15. PERMUTATION FEATURE IMPORTANCE (PARALLELIZED)
# ============================================
# 
# Key insight: For elastic net, we ONLY calculate permutation importance
# for features with non-zero coefficients. Features with β = 0 are not
# used by the model, so permuting them has zero effect by definition.
#
# Parallelization: Uses mclapply to parallelize across features within each fold
# ============================================

#' Calculate permutation importance for a single feature
#' Helper function for parallel execution
#' @param feat Feature name to test
#' @param X_test Test feature matrix
#' @param y_test Test outcomes
#' @param y_pred_baseline Baseline predictions
#' @param mse_baseline Baseline MSE
#' @param cv_fit Fitted cv.glmnet object
#' @param n_permutations Number of permutations
#' @param seed Random seed for this feature
#' @return Tibble with feature importance results
calculate_single_feature_importance <- function(feat, X_test, y_test, 
                                                y_pred_baseline, mse_baseline,
                                                cv_fit, n_permutations, seed) {
  set.seed(seed)
  
  feat_idx <- which(colnames(X_test) == feat)
  
  if (length(feat_idx) == 0) {
    return(tibble(
      feature = feat,
      importance = NA_real_,
      importance_sd = NA_real_,
      p_value = NA_real_
    ))
  }
  
  # Permute and calculate MSE multiple times
  mse_permuted <- numeric(n_permutations)
  
  for (p in 1:n_permutations) {
    X_test_perm <- X_test
    X_test_perm[, feat_idx] <- sample(X_test_perm[, feat_idx])
    
    y_pred_perm <- as.numeric(predict(cv_fit, newx = X_test_perm, s = "lambda.min"))
    mse_permuted[p] <- mean((y_test - y_pred_perm)^2)
  }
  
  # Importance = mean increase in MSE when feature is permuted
  delta_mse <- mse_permuted - mse_baseline
  importance <- mean(delta_mse)
  importance_sd <- sd(delta_mse)
  
  # One-sided p-value: proportion of permutations where MSE didn't increase
  p_value <- mean(delta_mse <= 0)
  
  tibble(
    feature = feat,
    importance = importance,
    importance_sd = importance_sd,
    p_value = p_value
  )
}


#' Calculate permutation importance for a single fold (PARALLELIZED)
#' @param X_train Training feature matrix
#' @param X_test Test feature matrix  
#' @param y_train Training outcomes
#' @param y_test Test outcomes
#' @param cv_fit Fitted cv.glmnet object
#' @param features_to_test Character vector of feature names to test (only non-zero coefs)
#' @param n_permutations Number of permutations per feature
#' @param n_cores Number of CPU cores to use
#' @param seed Random seed
#' @return Tibble with feature, importance values, and p-values
calculate_fold_importance <- function(X_train, X_test, y_train, y_test,
                                      cv_fit, features_to_test,
                                      n_permutations = 100, 
                                      n_cores = 1,
                                      seed = 42) {
  set.seed(seed)
  
  # Baseline predictions and MSE
  y_pred_baseline <- as.numeric(predict(cv_fit, newx = X_test, s = "lambda.min"))
  mse_baseline <- mean((y_test - y_pred_baseline)^2)
  
  # Only test features that exist in X_test
  features_to_test <- intersect(features_to_test, colnames(X_test))
  
  if (length(features_to_test) == 0) {
    return(tibble(
      feature = character(),
      importance = numeric(),
      importance_sd = numeric(),
      p_value = numeric()
    ))
  }
  
  # Generate unique seeds for each feature (for reproducibility)
  feature_seeds <- seed + seq_along(features_to_test) * 1000
  
  # PARALLEL: Calculate importance for each feature
  importance_results <- parallel::mclapply(
    seq_along(features_to_test),
    function(i) {
      calculate_single_feature_importance(
        feat = features_to_test[i],
        X_test = X_test,
        y_test = y_test,
        y_pred_baseline = y_pred_baseline,
        mse_baseline = mse_baseline,
        cv_fit = cv_fit,
        n_permutations = n_permutations,
        seed = feature_seeds[i]
      )
    },
    mc.cores = n_cores
  )
  
  bind_rows(importance_results)
}


#' Get non-zero coefficient features from elastic net results
#' @param results_row Single row from results_df
#' @return Character vector of feature names with non-zero coefficients
get_nonzero_features <- function(results_row) {
  coefs <- results_row$all_coefficients[[1]]
  names <- results_row$all_predictor_names[[1]]
  
  if (is.null(coefs) || is.null(names)) return(character())
  
  names[coefs != 0]
}


#' Run permutation importance for a single model with nested CV (PARALLELIZED)
#' @param df_raw Raw data frame
#' @param outcome_col Name of outcome column
#' @param exclude_cols Columns to exclude from predictors
#' @param features_to_test Features with non-zero coefficients to test
#' @param alpha_grid Grid of alpha values for elastic net
#' @param nfolds Number of CV folds
#' @param n_permutations Number of permutations per feature per fold
#' @param n_cores Number of CPU cores for parallel processing
#' @param seed Random seed
#' @return List with per-fold and aggregated importance results
run_permutation_importance <- function(df_raw, outcome_col, exclude_cols,
                                       features_to_test,
                                       alpha_grid = seq(0, 1, by = 0.1),
                                       nfolds = 5,
                                       n_permutations = 100,
                                       n_cores = 1,
                                       col_miss_thresh = 0.4,
                                       row_miss_thresh = 0.8,
                                       seed = 42) {
  set.seed(seed)
  
  # Prepare data
  df_complete <- df_raw %>% filter(!is.na(.data[[outcome_col]]))
  y <- df_complete[[outcome_col]]
  
  all_exclude <- unique(c(exclude_cols, outcome_col, CVD_scores))
  
  # Apply missingness filtering (suppress output for cleaner logs)
  df_filtered <- suppressMessages(
    filter_by_missingness(df_complete, all_exclude, col_miss_thresh, row_miss_thresh)
  )
  y <- df_filtered[[outcome_col]]
  n <- length(y)
  
  if (n < nfolds * 2) {
    warning("Too few observations for CV")
    return(NULL)
  }
  
  # Create CV folds
  foldid <- sample(rep(1:nfolds, length.out = n))
  
  # Store results per fold
  fold_results <- list()
  
  for (k in 1:nfolds) {
    cat(sprintf("    Fold %d/%d: ", k, nfolds))
    
    idx_test <- which(foldid == k)
    idx_train <- setdiff(1:n, idx_test)
    
    df_train <- df_filtered[idx_train, , drop = FALSE]
    df_test <- df_filtered[idx_test, , drop = FALSE]
    y_train <- y[idx_train]
    y_test <- y[idx_test]
    
    # Preprocess (same as main elastic net)
    preprocessed <- tryCatch({
      preprocess_train_test(df_train, df_test, all_exclude)
    }, error = function(e) {
      cat(sprintf("Preprocessing failed: %s\n", e$message))
      return(NULL)
    })
    
    if (is.null(preprocessed) || ncol(preprocessed$X_train) == 0) {
      cat("No features after preprocessing\n")
      next
    }
    
    X_train <- preprocessed$X_train
    X_test <- preprocessed$X_test
    
    # Fit elastic net (same procedure as main analysis)
    best_cv_fit <- NULL
    best_mse <- Inf
    
    for (a in alpha_grid) {
      cv_fit <- tryCatch({
        cv.glmnet(x = X_train, y = y_train, alpha = a,
                  nfolds = min(5, nfolds - 1),
                  type.measure = "mse", standardize = FALSE, family = "gaussian")
      }, error = function(e) NULL)
      
      if (!is.null(cv_fit) && min(cv_fit$cvm) < best_mse) {
        best_mse <- min(cv_fit$cvm)
        best_cv_fit <- cv_fit
      }
    }
    
    if (is.null(best_cv_fit)) {
      cat("Model fitting failed\n")
      next
    }
    
    # Filter features_to_test to those available after preprocessing
    available_features <- colnames(X_test)
    features_this_fold <- features_to_test[features_to_test %in% available_features]
    
    # Also check for dummy-coded versions (feature_level format)
    for (f in features_to_test) {
      if (!(f %in% available_features)) {
        dummy_matches <- grep(paste0("^", f, "_"), available_features, value = TRUE)
        features_this_fold <- c(features_this_fold, dummy_matches)
      }
    }
    features_this_fold <- unique(features_this_fold)
    
    n_features <- length(features_this_fold)
    cat(sprintf("%d features (parallel on %d cores)... ", n_features, n_cores))
    
    if (n_features == 0) {
      cat("No matching features\n")
      next
    }
    
    # Calculate permutation importance for this fold (PARALLELIZED)
    fold_start <- Sys.time()
    
    fold_importance <- calculate_fold_importance(
      X_train = X_train,
      X_test = X_test,
      y_train = y_train,
      y_test = y_test,
      cv_fit = best_cv_fit,
      features_to_test = features_this_fold,
      n_permutations = n_permutations,
      n_cores = n_cores,
      seed = seed + k * 10000
    )
    
    fold_elapsed <- difftime(Sys.time(), fold_start, units = "secs")
    
    fold_importance$fold <- k
    fold_results[[k]] <- fold_importance
    
    cat(sprintf("done in %.1fs (top: %s, Δ MSE = %.4f)\n",
                as.numeric(fold_elapsed),
                fold_importance$feature[which.max(fold_importance$importance)],
                max(fold_importance$importance, na.rm = TRUE)))
  }
  
  # Combine fold results
  all_folds <- bind_rows(fold_results)
  
  if (nrow(all_folds) == 0) {
    warning("No results from any fold")
    return(NULL)
  }
  
  # Aggregate across folds
  aggregated <- all_folds %>%
    group_by(feature) %>%
    summarise(
      importance = mean(importance, na.rm = TRUE),
      importance_se = sd(importance, na.rm = TRUE) / sqrt(sum(!is.na(importance))),
      importance_sd_across_folds = sd(importance, na.rm = TRUE),
      importance_sd_within_fold = mean(importance_sd, na.rm = TRUE),
      # Combine p-values using Fisher's method
      p_value_fisher = tryCatch({
        pvals <- p_value[!is.na(p_value)]
        if (length(pvals) == 0) return(NA_real_)
        chi_stat <- -2 * sum(log(pmax(pvals, 1e-10)))
        pchisq(chi_stat, df = 2 * length(pvals), lower.tail = FALSE)
      }, error = function(e) NA_real_),
      p_value_mean = mean(p_value, na.rm = TRUE),
      n_folds = sum(!is.na(importance)),
      .groups = "drop"
    ) %>%
    mutate(
      p_value_bh = p.adjust(p_value_fisher, method = "BH"),
      significant_005 = p_value_bh < 0.05,
      significant_010 = p_value_bh < 0.10
    ) %>%
    arrange(desc(importance))
  
  list(
    per_fold = all_folds,
    aggregated = aggregated,
    n_folds_cv = nfolds,
    n_permutations = n_permutations,
    n_observations = n
  )
}


#' Build combined dataset for score-specific + block models
build_model_dataset <- function(cvd_score, dataset_name, 
                                score_specific_datasets, datasets_list) {
  
  if (grepl("_specific$", dataset_name)) {
    return(score_specific_datasets[[cvd_score]])
  }
  
  if (dataset_name %in% names(datasets_list)) {
    return(datasets_list[[dataset_name]])
  }
  
  block <- sub("^.*_specific\\+", "", dataset_name)
  
  if (!(block %in% names(datasets_list))) {
    warning(sprintf("Block '%s' not found in datasets_list", block))
    return(NULL)
  }
  
  base_df <- score_specific_datasets[[cvd_score]]
  extra_df <- datasets_list[[block]]
  extra_predictors <- prep_extra_predictors(extra_df)
  
  dplyr::left_join(base_df, extra_predictors, by = "Sample_ID")
}


#' Run permutation importance for multiple models (PARALLELIZED)
#' @param results_df Results from elastic net analysis
#' @param datasets_list List of predictor block datasets
#' @param score_specific_datasets List of score-specific datasets
#' @param n_permutations Number of permutations per feature per fold
#' @param nfolds Number of CV folds
#' @param n_cores Number of CPU cores for parallel processing
#' @param seed Random seed
#' @param models_to_run Optional tibble with cvd_score and dataset_name columns
#' @return List with all results
run_permutation_importance_batch <- function(results_df, 
                                             datasets_list,
                                             score_specific_datasets,
                                             n_permutations = 100,
                                             nfolds = 5,
                                             n_cores = 1,
                                             seed = 42,
                                             models_to_run = NULL) {
  
  if (is.null(models_to_run)) {
    models_to_run <- results_df %>%
      select(cvd_score, dataset_name) %>%
      distinct()
  }
  
  n_models <- nrow(models_to_run)
  
  cat("==========================================================\n")
  cat("PERMUTATION FEATURE IMPORTANCE (PARALLELIZED)\n")
  cat(sprintf("Models: %d | Permutations/feature/fold: %d | Folds: %d | Cores: %d\n", 
              n_models, n_permutations, nfolds, n_cores))
  cat("==========================================================\n")
  
  all_results <- list()
  
  for (i in 1:n_models) {
    cvd_score <- models_to_run$cvd_score[i]
    dataset_name <- models_to_run$dataset_name[i]
    
    cat(sprintf("\n[%d/%d] %s ~ %s\n", i, n_models, cvd_score, dataset_name))
    
    model_results <- results_df %>%
      filter(cvd_score == !!cvd_score, dataset_name == !!dataset_name)
    
    if (nrow(model_results) == 0) {
      cat("  ✗ No elastic net results found\n")
      next
    }
    
    features_to_test <- get_nonzero_features(model_results[1, ])
    n_nonzero <- length(features_to_test)
    
    cat(sprintf("  Non-zero coefficient features: %d\n", n_nonzero))
    
    if (n_nonzero == 0) {
      cat("  ✗ No features with non-zero coefficients (null model)\n")
      next
    }
    
    df_raw <- build_model_dataset(cvd_score, dataset_name, 
                                  score_specific_datasets, datasets_list)
    
    if (is.null(df_raw)) {
      cat("  ✗ Could not build dataset\n")
      next
    }
    
    model_start <- Sys.time()
    
    tryCatch({
      importance_results <- run_permutation_importance(
        df_raw = df_raw,
        outcome_col = cvd_score,
        exclude_cols = "Sample_ID",
        features_to_test = features_to_test,
        alpha_grid = seq(0, 1, by = 0.1),
        nfolds = nfolds,
        n_permutations = n_permutations,
        n_cores = n_cores,
        seed = seed + i * 100000
      )
      
      model_elapsed <- difftime(Sys.time(), model_start, units = "mins")
      
      if (!is.null(importance_results)) {
        importance_results$aggregated$cvd_score <- cvd_score
        importance_results$aggregated$dataset_name <- dataset_name
        importance_results$aggregated$n_observations <- importance_results$n_observations
        importance_results$aggregated$n_permutations <- importance_results$n_permutations
        importance_results$aggregated$n_folds_cv <- importance_results$n_folds_cv
        
        all_results[[i]] <- importance_results
        
        n_sig <- sum(importance_results$aggregated$significant_005, na.rm = TRUE)
        top_feat <- importance_results$aggregated$feature[1]
        top_imp <- importance_results$aggregated$importance[1]
        
        cat(sprintf("  ✓ Complete in %.1f min: %d significant (q < 0.05), top = %s (Δ MSE = %.4f)\n",
                    as.numeric(model_elapsed), n_sig, top_feat, top_imp))
      } else {
        cat("  ✗ No results returned\n")
      }
    }, error = function(e) {
      cat(sprintf("  ✗ ERROR: %s\n", e$message))
    })
  }
  
  aggregated_all <- bind_rows(lapply(all_results, function(x) {
    if (!is.null(x)) x$aggregated else NULL
  }))
  
  cat(sprintf("\n✓ Completed %d models\n", length(all_results[!sapply(all_results, is.null)])))
  
  list(
    all_results = all_results,
    aggregated = aggregated_all
  )
}


#' Smart truncation: keep start and end of long names
smart_truncate <- function(x, max_len = 30) {
  ifelse(
    nchar(x) > max_len,
    paste0(substr(x, 1, 12), "..", substr(x, nchar(x) - 12, nchar(x))),
    x
  )
}


#' Extract coefficient signs from elastic net results for direction labeling
extract_coefficient_directions <- function(results_df) {
  direction_list <- list()
  
  for (i in 1:nrow(results_df)) {
    predictor_names <- results_df$all_predictor_names[[i]]
    coefficients <- results_df$all_coefficients[[i]]
    
    if (is.null(predictor_names) || is.null(coefficients)) next
    
    direction_list[[i]] <- tibble(
      cvd_score = results_df$cvd_score[i],
      dataset_name = results_df$dataset_name[i],
      feature = predictor_names,
      coefficient = coefficients,
      direction = case_when(
        coefficients > 0 ~ "Risk-Increasing",
        coefficients < 0 ~ "Risk-Decreasing",
        TRUE ~ "Zero"
      )
    )
  }
  
  bind_rows(direction_list)
}


#' Create bar plot for a single model's important features
plot_importance_bars <- function(importance_data, coefficient_directions, 
                                 cvd_score, dataset_name) {
  
  plot_data <- importance_data %>%
    left_join(
      coefficient_directions %>%
        filter(cvd_score == !!cvd_score, dataset_name == !!dataset_name) %>%
        select(feature, coefficient, direction),
      by = "feature"
    ) %>%
    filter(!is.na(importance)) %>%
    arrange(desc(importance)) %>%
    mutate(
      feature = factor(feature, levels = rev(feature)),
      direction = ifelse(is.na(direction), "Unknown", direction),
      sig_label = case_when(
        p_value_bh < 0.001 ~ "***",
        p_value_bh < 0.01 ~ "**",
        p_value_bh < 0.05 ~ "*",
        TRUE ~ ""
      ),
      label_y = ifelse(importance > max(importance) * 0.08, 
                       importance * 0.5, 
                       importance + max(importance) * 0.02)
    )
  
  if (nrow(plot_data) == 0) {
    warning("No features to plot")
    return(NULL)
  }
  
  p <- ggplot(plot_data, aes(x = feature, y = importance, fill = direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sig_label, y = label_y), 
              hjust = 0.5, size = 3.5, color = "white", fontface = "bold") +
    scale_fill_manual(
      values = c("Risk-Increasing" = "#E07B39", 
                 "Risk-Decreasing" = "#3B7EA1",
                 "Unknown" = "grey50"),
      name = "Effect Direction",
      labels = c("Risk-Increasing" = "Risk-Increasing (β > 0)",
                 "Risk-Decreasing" = "Risk-Decreasing (β < 0)",
                 "Unknown" = "Unknown")
    ) +
    coord_flip() +
    labs(
      title = paste0("Permutation Feature Importance: ", cvd_score),
      subtitle = paste0("Dataset: ", dataset_name, " (", nrow(plot_data), 
                        " non-zero coefficient features | *: q < 0.05, **: q < 0.01, ***: q < 0.001)"),
      x = "Feature",
      y = "Importance (Δ MSE when permuted)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 7),
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}


#' Create importance plots for all models
create_importance_plots <- function(importance_results, coefficient_directions, 
                                    output_dir = "importance_plots") {
  full_output_dir <- save_output(output_dir)
  if (!dir.exists(full_output_dir)) {
    dir.create(full_output_dir, recursive = TRUE)
  }
  
  importance_df <- importance_results$aggregated
  
  models <- importance_df %>%
    select(cvd_score, dataset_name) %>%
    distinct()
  
  cat(sprintf("\nCreating %d importance plots...\n", nrow(models)))
  
  plot_list <- list()
  
  for (i in 1:nrow(models)) {
    cvd <- models$cvd_score[i]
    dataset <- models$dataset_name[i]
    
    model_importance <- importance_df %>%
      filter(cvd_score == cvd, dataset_name == dataset)
    
    if (nrow(model_importance) == 0) {
      cat(sprintf("  [%d/%d] %s ~ %s: No data\n", i, nrow(models), cvd, dataset))
      next
    }
    
    n_sig <- sum(model_importance$significant_005, na.rm = TRUE)
    cat(sprintf("  [%d/%d] %s ~ %s: %d features, %d significant (q < 0.05)\n",
                i, nrow(models), cvd, dataset, nrow(model_importance), n_sig))
    
    p <- plot_importance_bars(model_importance, coefficient_directions, cvd, dataset)
    
    if (is.null(p)) next
    
    plot_list[[i]] <- p
    
    filename <- file.path(full_output_dir, paste0("importance_", cvd, "_",
                                                  gsub("[^A-Za-z0-9_]", "_", dataset), ".pdf"))
    
    n_features <- nrow(model_importance)
    plot_height <- max(6, min(25, 3 + 0.3 * n_features))
    
    ggsave(filename, plot = p, width = 10, height = plot_height, device = "pdf")
  }
  
  cat(sprintf("\n✓ Saved %d plots to %s/\n", 
              length(plot_list[!sapply(plot_list, is.null)]), output_dir))
  
  return(plot_list)
}


#' Create faceted plot for one predictor block across all CVD scores
plot_block_importance <- function(importance_data, coefficient_directions, block_name,
                                  top_n_per_score = 15) {
  
  plot_data <- importance_data %>%
    left_join(
      coefficient_directions %>%
        select(cvd_score, dataset_name, feature, coefficient, direction),
      by = c("cvd_score", "dataset_name", "feature")
    ) %>%
    filter(!is.na(importance), importance > 0) %>%
    group_by(cvd_score) %>%
    arrange(desc(importance)) %>%
    slice_head(n = top_n_per_score) %>%
    ungroup()
  
  if (nrow(plot_data) == 0) {
    warning("No data to plot for block")
    return(NULL)
  }
  
  n_features <- n_distinct(plot_data$feature)
  x_text_size <- case_when(
    block_name %in% c("all_data", "metabolites") ~ 5,
    n_features > 30 ~ 6,
    TRUE ~ 7
  )
  
  plot_data <- plot_data %>%
    mutate(
      cvd_score_label = recode(cvd_score,
                               "ascvd_10y" = "ASCVD",
                               "frs_10y" = "Framingham",
                               "QRISK3_risk" = "QRISK3",
                               "SCORE2_score" = "SCORE2",
                               "mean_risk" = "Composite"
      ),
      direction = ifelse(is.na(direction), "Unknown", direction),
      sig_label = case_when(
        p_value_bh < 0.001 ~ "***",
        p_value_bh < 0.01 ~ "**",
        p_value_bh < 0.05 ~ "*",
        TRUE ~ ""
      ),
      feature_short = smart_truncate(feature, max_len = 30),
      feature_id = paste(feature_short, cvd_score_label, sep = "___")
    ) %>%
    arrange(cvd_score_label, desc(importance))
  
  plot_data$feature_id <- factor(plot_data$feature_id, levels = unique(plot_data$feature_id))
  plot_data$cvd_score_label <- factor(plot_data$cvd_score_label,
                                      levels = c("ASCVD", "Framingham", "QRISK3", "SCORE2", "Composite"))
  
  plot_data <- plot_data %>%
    group_by(cvd_score_label) %>%
    mutate(
      max_importance = max(importance),
      label_y = ifelse(importance > max_importance * 0.08,
                       importance * 0.5,
                       importance + max_importance * 0.02)
    ) %>%
    ungroup()
  
  ggplot(plot_data, aes(x = feature_id, y = importance, fill = direction)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sig_label, y = label_y), 
              size = 2.5, color = "white", fontface = "bold") +
    facet_wrap(~ cvd_score_label, scales = "free_x", nrow = 1) +
    scale_x_discrete(labels = function(x) gsub("___.*$", "", x)) +
    scale_fill_manual(
      values = c("Risk-Increasing" = "#E07B39", 
                 "Risk-Decreasing" = "#3B7EA1",
                 "Unknown" = "grey50"),
      name = "Effect Direction"
    ) +
    labs(
      title = sprintf("Permutation Feature Importance by CVD Score: %s", block_name),
      subtitle = sprintf("Top %d features per score | *: q < 0.05, **: q < 0.01, ***: q < 0.001", 
                         top_n_per_score),
      x = NULL,
      y = "Importance (Δ MSE)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "gray40", hjust = 0.5),
      axis.text.x = element_text(size = x_text_size, angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    )
}


#' Create separate plots for each predictor block
create_block_importance_plots <- function(importance_results, coefficient_directions,
                                          blocks = c("all_data", "metabolites", "lipids", "fatty_acids", "urine_nmr", 
                                                     "body_composition", "clinical_risk_factors", "sociodemographics_lifestyle"),
                                          output_dir = "importance_plots_by_block",
                                          top_n_per_score = 15) {
  
  full_output_dir <- save_output(output_dir)
  dir.create(full_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  importance_df <- importance_results$aggregated
  
  all_importance <- importance_df %>%
    mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
    filter(block %in% blocks) %>%
    filter(cvd_score != "mean_risk")
  
  cat(sprintf("\nCreating block importance plots for %d blocks (excluding mean_risk)...\n", length(blocks)))
  
  for (current_block in blocks) {
    block_data <- filter(all_importance, block == current_block)
    
    if (nrow(block_data) == 0) {
      cat(sprintf("  %s: No data - skipping\n", current_block))
      next
    }
    
    n_features <- n_distinct(block_data$feature)
    n_scores <- n_distinct(block_data$cvd_score)
    cat(sprintf("  %s: %d unique features across %d CVD scores\n",
                current_block, n_features, n_scores))
    
    p <- plot_block_importance(block_data, coefficient_directions, current_block, top_n_per_score)
    
    if (is.null(p)) next
    
    filename <- file.path(full_output_dir, sprintf("importance_block_%s.pdf", current_block))
    
    n_features_shown <- block_data %>%
      filter(importance > 0) %>%
      group_by(cvd_score) %>%
      slice_head(n = top_n_per_score) %>%
      ungroup() %>%
      n_distinct(feature)
    
    plot_width <- max(12, min(22, 6 + 0.18 * n_features_shown))
    
    ggsave(filename, p, width = plot_width, height = 7)
  }
  
  cat(sprintf("✓ Saved plots to %s/\n", output_dir))
}


# ============================================
# 16. EXECUTE PERMUTATION IMPORTANCE (PARALLELIZED)
# ============================================

cat("\n==========================================================\n")
cat("PERMUTATION IMPORTANCE: FULL RUN (PARALLELIZED)\n")
cat("Dataset: All score+block combinations | Permutations: 1000\n")
cat(sprintf("Using %d CPU cores\n", CPUS))
cat("==========================================================\n")

# Define all models from results_score_plus_blocks
models_for_full_run <- results_score_plus_blocks %>%
  select(cvd_score, dataset_name) %>%
  distinct()

cat(sprintf("Running %d models\n", nrow(models_for_full_run)))

# Time the full run
start_time <- Sys.time()

# Run permutation importance (full) - NOW PARALLELIZED
importance_results_full <- run_permutation_importance_batch(
  results_df = results_score_plus_blocks,
  datasets_list = datasets_list,
  score_specific_datasets = score_specific_datasets,
  n_permutations = 1000,
  nfolds = 5,
  n_cores = CPUS,  # USE ALL AVAILABLE CORES
  seed = 42,
  models_to_run = models_for_full_run
)

end_time <- Sys.time()
total_runtime <- difftime(end_time, start_time, units = "hours")
cat(sprintf("\nTotal runtime: %.2f hours\n", as.numeric(total_runtime)))

# Save results (multiple formats for safety)
cat("\nSaving results...\n")
qs_save(importance_results_full, save_output("permutation_importance_full_results.qs2"))
saveRDS(importance_results_full, save_output("permutation_importance_full_results.rds"))

# Export aggregated results to Excel
importance_full_excel <- importance_results_full$aggregated %>%
  select(
    cvd_score, dataset_name, feature,
    importance, importance_se, importance_sd_across_folds, importance_sd_within_fold,
    p_value_fisher, p_value_mean, p_value_bh,
    significant_005, significant_010,
    n_observations, n_permutations, n_folds_cv
  ) %>%
  arrange(cvd_score, dataset_name, desc(importance))

write.xlsx(importance_full_excel, save_output("permutation_importance_full_results.xlsx"))

# Extract coefficient directions for color-coding
cat("\nExtracting coefficient directions from elastic net results...\n")
coefficient_directions_full <- extract_coefficient_directions(results_score_plus_blocks)

# Create individual model plots
cat("\nCreating individual model plots...\n")
create_importance_plots(importance_results_full, coefficient_directions_full,
                        output_dir = "importance_plots_full")

# Create block-level plots
cat("\nCreating block-level plots...\n")
create_block_importance_plots(
  importance_results = importance_results_full,
  coefficient_directions = coefficient_directions_full,
  output_dir = "importance_plots_by_block_full",
  top_n_per_score = 15
)

# Summary statistics with direction breakdown
cat("\n==========================================================\n")
cat("FULL RUN SUMMARY\n")
cat("==========================================================\n")

importance_full_summary <- importance_results_full$aggregated %>%
  left_join(coefficient_directions_full, by = c("cvd_score", "dataset_name", "feature")) %>%
  group_by(cvd_score, dataset_name) %>%
  summarise(
    n_features_tested = n(),
    n_positive_importance = sum(importance > 0, na.rm = TRUE),
    n_significant_005 = sum(significant_005 & importance > 0, na.rm = TRUE),
    n_significant_010 = sum(significant_010 & importance > 0, na.rm = TRUE),
    n_risk_increasing = sum(importance > 0 & direction == "Risk-Increasing", na.rm = TRUE),
    n_risk_decreasing = sum(importance > 0 & direction == "Risk-Decreasing", na.rm = TRUE),
    pct_significant = round(100 * mean(significant_005 & importance > 0, na.rm = TRUE), 1),
    top_feature = feature[which.max(importance)],
    top_importance = round(max(importance, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  arrange(cvd_score, desc(n_significant_005))

print(importance_full_summary, n = Inf)

write.xlsx(importance_full_summary, save_output("permutation_importance_full_summary.xlsx"))

# Create combined summary table
importance_overview <- importance_full_summary %>%
  select(cvd_score, dataset_name, n_features_tested, n_significant_005, 
         n_risk_increasing, n_risk_decreasing, top_feature, top_importance) %>%
  mutate(
    block = sub("^.*_specific\\+", "", dataset_name),
    cvd_label = cvd_score_labels[cvd_score]
  )

write.xlsx(importance_overview, save_output("permutation_importance_overview.xlsx"))

cat("\n==========================================================\n")
cat("✓ PERMUTATION IMPORTANCE ANALYSIS COMPLETE\n")
cat(sprintf("  Results saved to: %s\n", OUTPUT_DIR))
cat(sprintf("  Total runtime: %.2f hours\n", as.numeric(total_runtime)))
cat("==========================================================\n")