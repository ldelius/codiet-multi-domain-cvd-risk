### Elastic Net Regression with Covariate Adjustment (Country, Statins, Supplements unpenalized)
### Preprocessing Inside Each CV-Fold (No Data Leakage)
### Author: Luisa Delius
### Based on original pipeline; adjusted for covariates, removed composite score & metabolites block

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
library(circlize)
library(cowplot)

set.seed(42)

wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files/processed_data"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

create_output_folder <- function(base_name = "elastic_net_covariate_adjusted") {
  existing <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
  pattern <- paste0("^", base_name, "_run(\\d+)$")
  matches <- grep(pattern, existing, value = TRUE)
  if (length(matches) == 0) { run_number <- 1
  } else { run_number <- max(as.integer(gsub(pattern, "\\1", matches))) + 1 }
  folder_name <- sprintf("%s_run%03d", base_name, run_number)
  dir.create(folder_name, showWarnings = FALSE)
  cat(sprintf("✓ Created output folder: %s\n", folder_name))
  return(folder_name)
}

OUTPUT_DIR <- create_output_folder()
save_output <- function(filename) file.path(OUTPUT_DIR, filename)

CPUS <- parallel::detectCores() - 2
cat(sprintf("Using %d CPU cores for parallel processing\n", CPUS))

# ============================================
# 1. DATA LOADING
# ============================================

df_fatty_acids_raw_full <- readRDS("df_fatty_acids_predictor_statin_suppl.rds")
df_risk_factors_raw_full <- readRDS("df_risk_factor_predictors.rds")

df_ascvd_frs_score2_input <- readRDS("ASCVD_SCORE2_Framingham_input.rds") %>%
  rename(Sample_ID = PatientID) %>%
  mutate(
    blood_pressure_treatment = factor(blood_pressure_treatment,
                                      levels = c(0, 1), labels = c("No", "Yes")))

df_QRISK3_input <- readRDS("QRISK3_calculation_input.rds") %>%
  rename(Sample_ID = PatientID)

# Source covariates BEFORE dropping them (df_QRISK3_input must be loaded first)
df_covariates <- df_fatty_acids_raw_full %>%
  select(Sample_ID, Statins, Supplements) %>%
  full_join(df_risk_factors_raw_full %>% select(Sample_ID, Country), by = "Sample_ID") %>%
  left_join(df_QRISK3_input %>% select(Sample_ID, Sex, Age), by = "Sample_ID") %>%
  distinct(Sample_ID, .keep_all = TRUE)

COVARIATE_COLS <- c("Country", "Statins", "Supplements")
COVARIATE_COLS_COMMON <- c("Country", "Statins", "Supplements", "Sex", "Age")
cat(sprintf("✓ Covariates table: %d samples\n", nrow(df_covariates)))

df_fatty_acids_ElaNet <- readRDS("df_fatty_acids_predictor_statin_suppl.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk, -Statins, -Supplements)

df_lipidomics_ElaNet <- readRDS("df_lipidomics_predictor_statin_suppl.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk, -Statins, -Supplements) %>%
  filter(!is.na(Sample_ID))

df_urine_nmr_ElaNet <- readRDS("df_urine_NMR_data.rds") %>%
  select(-starts_with("z_"))

df_body_composition_ElaNet <- readRDS("df_body_composition_metrics.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk) %>%
  full_join(df_risk_factors_raw_full %>% select(Body.fat, Fat.free.mass, Sample_ID), by = "Sample_ID") %>%
  filter(!is.na(Sample_ID))

df_REDcap_demographics_ElaNet <- readRDS("df_REDcap_demographics_ElaNet.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk, -recruitment_site, -any_of(COVARIATE_COLS_COMMON))

df_risk_factors_ElaNet <- readRDS("df_risk_factor_predictors.rds") %>%
  select(-starts_with("z_"), -QRISK3_risk,
         -Heart.Rate, -Age, -Total.Cholesterol.mg.dl, -HDL.mg.dl,
         -Body.Weight, -Height, -BMI, -Gender,
         -stress_resilience_status, -stress_index_status, -Age.Risk,
         -Country, -Body.fat, -Fat.free.mass)

df_all_cvd_risk_scores_ElaNet <- readRDS("df_all_risk_scores.rds") %>%
  rename(Sample_ID = PatientID) %>%
  select(-SCORE2_strat, -any_of("mean_risk"))

CVD_scores <- c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y")

# ============================================
# 2. CREATE DATASETS
# ============================================

# Helper: attach covariates to a NAMED LIST of data frames
attach_covariates <- function(datasets_list, covariate_cols = COVARIATE_COLS) {
  lapply(datasets_list, function(df) {
    covs_to_add <- df_covariates %>% select(Sample_ID, all_of(covariate_cols))
    df %>%
      select(-any_of(covariate_cols)) %>%
      left_join(covs_to_add, by = "Sample_ID")
  })
}

# Helper: attach covariates to a SINGLE data frame
attach_covariates_single <- function(df, covariate_cols = COVARIATE_COLS) {
  covs_to_add <- df_covariates %>% select(Sample_ID, all_of(covariate_cols))
  df %>%
    select(-any_of(covariate_cols)) %>%
    left_join(covs_to_add, by = "Sample_ID")
}

# ---- Build raw data frames (no covariates yet) ----

df_body_comp_raw <- df_body_composition_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")
df_demographics_raw <- df_REDcap_demographics_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")
df_risk_factors_raw <- df_risk_factors_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")
df_lipids_raw <- df_lipidomics_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")
df_fatty_acids_raw <- df_fatty_acids_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")
df_urine_nmr_raw <- df_urine_nmr_ElaNet %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID")

df_all_data <- df_fatty_acids_ElaNet %>%
  full_join(df_lipidomics_ElaNet, by = "Sample_ID") %>%
  full_join(df_urine_nmr_ElaNet, by = "Sample_ID") %>%
  full_join(df_body_composition_ElaNet, by = "Sample_ID") %>%
  full_join(df_risk_factors_ElaNet, by = "Sample_ID") %>%
  full_join(df_REDcap_demographics_ElaNet, by = "Sample_ID") %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID") %>%
  select(where(~!all(is.na(.))))

# ---- Assemble into list, then attach COMMON covariates (incl. Sex) ----

datasets_list <- list(
  all_data = df_all_data,
  body_composition = df_body_comp_raw,
  sociodemographics_lifestyle = df_demographics_raw,
  clinical_risk_factors = df_risk_factors_raw,
  lipids = df_lipids_raw,
  fatty_acids = df_fatty_acids_raw,
  urine_nmr = df_urine_nmr_raw
)
datasets_list <- attach_covariates(datasets_list, COVARIATE_COLS)

# ---- Score-specific input datasets (Sex stays as predictor, not covariate) ----

df_model0_QRISK3 <- df_QRISK3_input %>%
  mutate(townsend = replace_na(townsend, 0)) %>%
  full_join(df_all_cvd_risk_scores_ElaNet %>% select(Sample_ID, QRISK3_risk), by = "Sample_ID") %>%
  filter(!is.na(QRISK3_risk)) %>%
  attach_covariates_single(COVARIATE_COLS)

df_model0_ascvd <- df_ascvd_frs_score2_input %>%
  select(-Risk.region, -mean_LDL_mg_dl) %>%
  full_join(df_all_cvd_risk_scores_ElaNet %>% select(Sample_ID, ascvd_10y), by = "Sample_ID") %>%
  filter(!is.na(ascvd_10y)) %>%
  attach_covariates_single(COVARIATE_COLS)

df_model0_score2 <- df_ascvd_frs_score2_input %>%
  select(-mean_LDL_mg_dl, -blood_pressure_treatment, -race_ascvd) %>%
  full_join(df_all_cvd_risk_scores_ElaNet %>% select(Sample_ID, SCORE2_score), by = "Sample_ID") %>%
  filter(!is.na(SCORE2_score)) %>%
  attach_covariates_single(COVARIATE_COLS)

df_model0_frs <- df_ascvd_frs_score2_input %>%
  select(-Risk.region, -mean_LDL_mg_dl, -race_ascvd) %>%
  full_join(df_all_cvd_risk_scores_ElaNet %>% select(Sample_ID, frs_10y), by = "Sample_ID") %>%
  filter(!is.na(frs_10y)) %>%
  attach_covariates_single(COVARIATE_COLS)

score_specific_datasets <- list(
  QRISK3_risk = df_model0_QRISK3,
  SCORE2_score = df_model0_score2,
  frs_10y = df_model0_frs,
  ascvd_10y = df_model0_ascvd
)

# ---- Body Composition + Height/Weight ----

height_weight_vars <- df_QRISK3_input %>% select(Sample_ID, Height_cm, Weight_kg)
df_body_comp_extended <- df_body_composition_ElaNet %>%
  full_join(height_weight_vars, by = "Sample_ID") %>%
  full_join(df_all_cvd_risk_scores_ElaNet, by = "Sample_ID") %>%
  attach_covariates_single(COVARIATE_COLS_COMMON)

scores_without_ht_wt <- c("SCORE2_score", "ascvd_10y", "frs_10y")

cat("✓ All datasets created with covariates attached\n")

# ============================================
# 3. PREPROCESSING FUNCTIONS
# ============================================

filter_by_missingness <- function(df, exclude_cols, col_miss_thresh = 0.4, row_miss_thresh = 0.8) {
  pred_cols <- setdiff(names(df), exclude_cols)
  col_miss_rates <- sapply(df[pred_cols], function(x) mean(is.na(x)))
  cols_to_keep <- names(col_miss_rates)[col_miss_rates <= col_miss_thresh]
  cols_removed <- length(pred_cols) - length(cols_to_keep)
  if (cols_removed > 0) cat(sprintf("    Removed %d columns with >%.0f%% missing\n", cols_removed, col_miss_thresh * 100))
  if (length(cols_to_keep) > 0) {
    row_miss_rates <- rowMeans(is.na(df[cols_to_keep]))
    rows_to_keep <- which(row_miss_rates <= row_miss_thresh)
    rows_removed <- nrow(df) - length(rows_to_keep)
    if (rows_removed > 0) cat(sprintf("    Removed %d rows with >%.0f%% missing\n", rows_removed, row_miss_thresh * 100))
  } else { rows_to_keep <- 1:nrow(df) }
  df[rows_to_keep, c(intersect(exclude_cols, names(df)), cols_to_keep), drop = FALSE]
}

preprocess_train_test <- function(df_train, df_test, exclude_cols = NULL,
                                  covariate_cols = COVARIATE_COLS, k = 5) {
  pred_cols <- setdiff(names(df_train), exclude_cols)
  rec <- recipe(~ ., data = df_train[, pred_cols, drop = FALSE]) %>%
    step_impute_knn(all_predictors(), neighbors = k) %>%
    step_string2factor(all_nominal_predictors()) %>%
    step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
    step_normalize(all_numeric_predictors())
  rec_prepped <- prep(rec, training = df_train[, pred_cols, drop = FALSE])
  X_train <- bake(rec_prepped, new_data = df_train[, pred_cols, drop = FALSE]) %>% as.matrix()
  X_test <- bake(rec_prepped, new_data = df_test[, pred_cols, drop = FALSE]) %>% as.matrix()
  common_cols <- intersect(colnames(X_train), colnames(X_test))
  X_train <- X_train[, common_cols, drop = FALSE]
  X_test <- X_test[, common_cols, drop = FALSE]
  is_covariate <- grepl(paste0("^(", paste(covariate_cols, collapse = "|"), ")"), common_cols)
  list(X_train = X_train, X_test = X_test, recipe = rec_prepped,
       feature_names = common_cols, is_covariate = is_covariate)
}

# ============================================
# 3.5 HELPER: Alpha selection by minimum RMSE
# ============================================
# Selects the alpha with the lowest cross-validated RMSE at lambda.1se.
# Lambda selection uses the one-SE rule (lambda.1se) per Hastie et al.,
# while alpha is chosen purely by performance.

select_best_alpha <- function(cv_results, alpha_grid, nfolds) {
  # Extract RMSE at lambda.1se for each alpha
  cv_rmse_vals <- sapply(cv_results, function(x) {
    idx <- which(x$cv_fit$lambda == x$cv_fit$lambda.1se)
    sqrt(x$cv_fit$cvm[idx])
  })
  
  best_idx <- which.min(cv_rmse_vals)
  
  list(
    best_idx    = best_idx,
    best_alpha  = alpha_grid[best_idx],
    best_rmse   = cv_rmse_vals[best_idx],
    all_rmse    = cv_rmse_vals
  )
}


# ============================================
# 4. PERMUTATION TEST FUNCTION
# ============================================

compute_permutation_Q2 <- function(df_raw, y_shuffled, exclude_cols,
                                   alpha_grid = seq(0, 1, by = 0.1), nfolds = 10, seed = 42) {
  set.seed(seed)
  n <- length(y_shuffled)
  foldid <- sample(rep(1:nfolds, length.out = n))
  y_oof <- rep(NA_real_, n)
  for (k in 1:nfolds) {
    idx_test <- which(foldid == k); idx_train <- setdiff(1:n, idx_test)
    df_train <- df_raw[idx_train, , drop = FALSE]; df_test <- df_raw[idx_test, , drop = FALSE]
    y_train <- y_shuffled[idx_train]
    preprocessed <- tryCatch(preprocess_train_test(df_train, df_test, exclude_cols), error = function(e) NULL)
    if (is.null(preprocessed) || ncol(preprocessed$X_train) == 0) return(NA_real_)
    pen_factor <- ifelse(preprocessed$is_covariate, 0, 1)
    
    # Fit cv.glmnet for each alpha
    cv_results <- list()
    valid_count <- 0
    for (a_idx in seq_along(alpha_grid)) {
      cv_fit <- tryCatch(cv.glmnet(x = preprocessed$X_train, y = y_train, alpha = alpha_grid[a_idx],
                                   penalty.factor = pen_factor, nfolds = nfolds,
                                   type.measure = "mse", standardize = FALSE, family = "gaussian"), error = function(e) NULL)
      if (!is.null(cv_fit)) {
        valid_count <- valid_count + 1
        cv_results[[valid_count]] <- list(alpha = alpha_grid[a_idx], cv_fit = cv_fit)
      }
    }
    if (valid_count == 0) return(NA_real_)
    
    # Select alpha by minimum RMSE (with fallback)
    valid_alphas <- sapply(cv_results, function(x) x$alpha)
    alpha_selection <- tryCatch(
      select_best_alpha(cv_results, valid_alphas, nfolds),
      error = function(e) {
        rmse_vals <- sapply(cv_results, function(x) {
          idx <- which(x$cv_fit$lambda == x$cv_fit$lambda.1se)
          sqrt(x$cv_fit$cvm[idx])
        })
        list(best_idx = which.min(rmse_vals))
      }
    )
    
    best_cv_fit <- cv_results[[alpha_selection$best_idx]]$cv_fit
    y_oof[idx_test] <- as.numeric(predict(best_cv_fit, newx = preprocessed$X_test, s = "lambda.1se"))
  }
  ss_res <- sum((y_shuffled - y_oof)^2); ss_tot <- sum((y_shuffled - mean(y_shuffled))^2)
  1 - ss_res / ss_tot
}

run_permutation_test <- function(df_raw, y, exclude_cols, observed_Q2,
                                 n_permutations = 100, alpha_grid = seq(0, 1, by = 0.1),
                                 nfolds = 10, seed = 42) {
  set.seed(seed)
  cat(sprintf("  Running permutation test (%d permutations) on %d cores...\n", n_permutations, CPUS))
  shuffled_y_list <- lapply(1:n_permutations, function(p) sample(y))
  null_Q2_values <- parallel::mclapply(1:n_permutations, function(p) {
    compute_permutation_Q2(df_raw, shuffled_y_list[[p]], exclude_cols, alpha_grid, nfolds, seed + p)
  }, mc.cores = CPUS)
  null_Q2_values <- unlist(null_Q2_values); null_Q2_values <- null_Q2_values[!is.na(null_Q2_values)]
  p_value <- (sum(null_Q2_values >= observed_Q2) + 1) / (length(null_Q2_values) + 1)
  cat(sprintf("  Permutation test complete: p = %.4f\n", p_value))
  list(p_value = p_value, observed_Q2 = observed_Q2, null_distribution = null_Q2_values,
       n_permutations_completed = length(null_Q2_values))
}

# ============================================
# 5. HELPER FUNCTIONS
# ============================================
prep_extra_predictors <- function(extra_df) {
  extra_df[, setdiff(names(extra_df), CVD_scores), drop = FALSE]
}

build_score_plus_block <- function(score, block_name, score_specific_datasets, datasets_list) {
  base_df <- score_specific_datasets[[score]]
  extra_predictors <- prep_extra_predictors(datasets_list[[block_name]]) %>%
    select(-any_of(COVARIATE_COLS))
  dplyr::left_join(base_df, extra_predictors, by = "Sample_ID")
}

# ============================================
# 6. CORE ELASTIC NET FUNCTION (WITH COVARIATE ADJUSTMENT)
# ============================================

get_oof_predictions_clean <- function(df_raw, y, exclude_cols,
                                      covariate_cols = COVARIATE_COLS,
                                      alpha_grid = seq(0, 1, by = 0.1), nfolds = 10, seed = 42) {
  set.seed(seed)
  n <- length(y); y_mean_global <- mean(y)
  foldid_outer <- sample(rep(1:nfolds, length.out = n))
  y_oof <- rep(NA_real_, n)
  alpha_selected_per_fold <- numeric(nfolds)
  n_features_per_fold <- numeric(nfolds)
  Q2_per_fold <- numeric(nfolds)
  MAE_per_fold <- numeric(nfolds)
  RMSE_per_fold <- numeric(nfolds)
  R2_per_fold <- numeric(nfolds)
  
  for (k in 1:nfolds) {
    idx_test <- which(foldid_outer == k)
    idx_train <- setdiff(1:n, idx_test)
    
    df_train <- df_raw[idx_train, , drop = FALSE]
    df_test <- df_raw[idx_test, , drop = FALSE]
    y_train <- y[idx_train]
    y_test <- y[idx_test]
    
    preprocessed <- preprocess_train_test(df_train, df_test, exclude_cols, covariate_cols)
    X_train <- preprocessed$X_train
    X_test <- preprocessed$X_test
    pen_factor <- ifelse(preprocessed$is_covariate, 0, 1)
    
    n_features_per_fold[k] <- ncol(X_train)
    
    cv_results_fold <- lapply(alpha_grid, function(a) {
      cv_fit <- cv.glmnet(
        x = X_train, y = y_train, alpha = a,
        penalty.factor = pen_factor,
        nfolds = nfolds,
        type.measure = "mse", standardize = FALSE, family = "gaussian"
      )
      list(alpha = a, cv_fit = cv_fit)
    })
    
    # Select alpha by minimum RMSE at lambda.1se
    alpha_selection <- select_best_alpha(cv_results_fold, alpha_grid, nfolds)
    best_idx <- alpha_selection$best_idx
    alpha_selected_per_fold[k] <- alpha_grid[best_idx]
    best_cv_fit <- cv_results_fold[[best_idx]]$cv_fit
    
    y_pred_k <- as.numeric(predict(best_cv_fit, newx = X_test, s = "lambda.1se"))
    y_oof[idx_test] <- y_pred_k
    
    ss_res_k <- sum((y_test - y_pred_k)^2)
    ss_tot_k <- sum((y_test - y_mean_global)^2)
    Q2_per_fold[k] <- 1 - ss_res_k / ss_tot_k
    MAE_per_fold[k] <- mean(abs(y_test - y_pred_k))
    RMSE_per_fold[k] <- sqrt(mean((y_test - y_pred_k)^2))
    
    y_pred_train_k <- as.numeric(predict(best_cv_fit, newx = X_train, s = "lambda.1se"))
    ss_res_train <- sum((y_train - y_pred_train_k)^2)
    ss_tot_train <- sum((y_train - mean(y_train))^2)
    R2_per_fold[k] <- 1 - ss_res_train / ss_tot_train
  }
  
  list(
    predictions = y_oof,
    alpha_per_fold = alpha_selected_per_fold,
    n_features_per_fold = n_features_per_fold,
    foldid = foldid_outer,
    Q2_per_fold = Q2_per_fold,
    MAE_per_fold = MAE_per_fold,
    RMSE_per_fold = RMSE_per_fold,
    R2_per_fold = R2_per_fold
  )
}

run_elastic_net_model <- function(df_raw, outcome_col, cvd_score_name, dataset_name,
                                  exclude_cols = "Sample_ID",
                                  covariate_cols = COVARIATE_COLS,
                                  alpha_grid = seq(0, 1, by = 0.1),
                                  nfolds = 10, n_permutations = 100,
                                  col_miss_thresh = 0.4, row_miss_thresh = 0.8,
                                  seed = 42) {
  set.seed(seed)
  
  # Remove rows with missing outcome
  df_complete <- df_raw %>% filter(!is.na(.data[[outcome_col]]))
  y <- df_complete[[outcome_col]]
  all_exclude <- unique(c(exclude_cols, outcome_col, CVD_scores))
  
  # Missingness filtering — protect covariates from removal
  cat(sprintf("  Filtering by missingness (col >%.0f%%, row >%.0f%%)...\n",
              col_miss_thresh * 100, row_miss_thresh * 100))
  cols_to_protect <- unique(c(all_exclude, covariate_cols))
  df_filtered <- filter_by_missingness(df_complete, cols_to_protect, col_miss_thresh, row_miss_thresh)
  
  # Re-attach covariates if they were dropped
  for (cov in covariate_cols) {
    if (!cov %in% names(df_filtered) && cov %in% names(df_complete)) {
      df_filtered[[cov]] <- df_complete[[cov]][match(df_filtered$Sample_ID, df_complete$Sample_ID)]
    }
  }
  
  y <- df_filtered[[outcome_col]]
  n_obs <- nrow(df_filtered)
  cat(sprintf("  Sample size after filtering: %d\n", n_obs))
  
  if (n_obs < nfolds * 2) {
    warning(sprintf("Very few observations (%d) for %d-fold CV — proceeding with caution", n_obs, nfolds))
  }
  
  # ---- PART A: Full-data model (for coefficients) ----
  preprocessed_full <- preprocess_train_test(df_filtered, df_filtered, all_exclude, covariate_cols)
  X_full <- preprocessed_full$X_train
  n_pred <- ncol(X_full)
  predictor_names <- colnames(X_full)
  pen_factor_full <- ifelse(preprocessed_full$is_covariate, 0, 1)
  n_covariates_in_model <- sum(preprocessed_full$is_covariate)
  n_penalized <- sum(!preprocessed_full$is_covariate)
  
  cat(sprintf("  Predictors: %d (%d penalized, %d unpenalized covariates)\n",
              n_pred, n_penalized, n_covariates_in_model))
  
  if (n_pred == 0) {
    warning("No predictors remaining after preprocessing")
    return(NULL)
  }
  
  set.seed(seed)
  foldid_inner <- sample(rep(1:nfolds, length.out = n_obs))
  
  cv_results <- lapply(alpha_grid, function(a) {
    cv_fit <- cv.glmnet(
      x = X_full, y = y, alpha = a,
      penalty.factor = pen_factor_full,
      nfolds = nfolds, foldid = foldid_inner,
      type.measure = "mse", standardize = FALSE, family = "gaussian"
    )
    list(alpha = a, cv_fit = cv_fit)
  })
  
  # Select alpha by minimum RMSE at lambda.1se
  alpha_selection <- select_best_alpha(cv_results, alpha_grid, nfolds)
  best_alpha_idx <- alpha_selection$best_idx
  best_alpha_fulldata <- alpha_grid[best_alpha_idx]
  best_cv_fit <- cv_results[[best_alpha_idx]]$cv_fit
  
  cat(sprintf("  Alpha selection: min RMSE = %.4f at alpha = %.1f\n",
              alpha_selection$best_rmse, best_alpha_fulldata))
  
  lambda_min <- best_cv_fit$lambda.min
  lambda_1se <- best_cv_fit$lambda.1se
  cv_rmse_min <- sqrt(min(best_cv_fit$cvm))
  cv_rmse_1se <- sqrt(best_cv_fit$cvm[which(best_cv_fit$lambda == lambda_1se)])
  
  final_coef <- coef(best_cv_fit, s = "lambda.1se")
  intercept <- as.numeric(final_coef[1])
  coef_vector <- as.numeric(final_coef[-1])
  names(coef_vector) <- predictor_names
  n_nonzero <- sum(coef_vector != 0)
  
  # ---- PART B: In-sample predictions (R²_Y) ----
  y_pred_in <- as.numeric(predict(best_cv_fit, newx = X_full, s = "lambda.1se"))
  res_in <- y - y_pred_in
  rmse_in <- sqrt(mean(res_in^2))
  mae_in <- mean(abs(res_in))
  R2_Y <- 1 - sum(res_in^2) / sum((y - mean(y))^2)
  cor_in <- cor(y_pred_in, y)
  
  # ---- PART C: Out-of-fold predictions (Q²_Y) ----
  oof_results <- get_oof_predictions_clean(
    df_filtered, y, all_exclude, covariate_cols, alpha_grid, nfolds, seed + 100
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
  R2_per_fold <- oof_results$R2_per_fold
  R2_fold_mean <- mean(R2_per_fold)
  R2_fold_sd <- sd(R2_per_fold)
  
  # ---- PART D: Model quality metrics ----
  Q2_R2_ratio <- Q2_Y / R2_Y
  R2_Q2_gap <- R2_Y - Q2_Y
  dev_explained <- best_cv_fit$glmnet.fit$dev.ratio[
    which.min(abs(best_cv_fit$glmnet.fit$lambda - lambda_1se))
  ]
  
  # ---- PART E: Permutation test ----
  if (n_permutations > 0) {
    perm_results <- run_permutation_test(
      df_filtered, y, all_exclude, Q2_Y,
      n_permutations, alpha_grid, nfolds, seed + 200
    )
    perm_p_value <- perm_results$p_value
    perm_n_completed <- perm_results$n_permutations_completed
    perm_null_mean <- mean(perm_results$null_distribution)
  } else {
    perm_p_value <- NA_real_
    perm_n_completed <- 0
    perm_null_mean <- NA_real_
  }
  
  # ---- Separate covariate vs penalized coefficients ----
  nonzero_coefs <- coef_vector[coef_vector != 0]
  is_cov_nz <- grepl(paste0("^(", paste(covariate_cols, collapse = "|"), ")"), names(nonzero_coefs))
  penalized_nonzero <- nonzero_coefs[!is_cov_nz]
  covariate_nonzero <- nonzero_coefs[is_cov_nz]
  
  # All non-zero penalized predictors, sorted by absolute coefficient
  if (length(penalized_nonzero) > 0) {
    sorted_coefs <- sort(abs(penalized_nonzero), decreasing = TRUE)
    nonzero_predictor_names <- names(sorted_coefs)
    nonzero_predictor_coefs <- penalized_nonzero[nonzero_predictor_names]
  } else {
    nonzero_predictor_names <- NA
    nonzero_predictor_coefs <- NA
  }
  
  # ---- Return results tibble ----
  tibble(
    model_name = "elastic_net_cov_adjusted",
    cvd_score = cvd_score_name,
    dataset_name = dataset_name,
    n_observations = n_obs,
    n_predictors = n_pred,
    n_penalized_predictors = n_penalized,
    n_covariates = n_covariates_in_model,
    covariates_used = paste(covariate_cols, collapse = ", "),
    n_nonzero_coefs = n_nonzero,
    n_nonzero_penalized = length(penalized_nonzero),
    n_nonzero_covariates = length(covariate_nonzero),
    alpha_selection_method = "min_rmse",
    alpha_optimal_fulldata = best_alpha_fulldata,
    alpha_mean_nested = mean(alpha_per_fold),
    lambda_min = lambda_min,
    lambda_1se = lambda_1se,
    lambda_selection_method = "one_se",
    n_folds = nfolds,
    cv_rmse_min = cv_rmse_min,
    cv_rmse_1se = cv_rmse_1se,
    R2_Y = R2_Y,
    R2_fold_mean = R2_fold_mean,
    R2_fold_sd = R2_fold_sd,
    R2_per_fold = list(R2_per_fold),
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
    Q2_fold_summary = sprintf("%.3f +/- %.3f", Q2_fold_mean, Q2_fold_sd),
    Q2_per_fold = list(Q2_per_fold),
    Q2_R2_ratio = Q2_R2_ratio,
    R2_Q2_gap = R2_Q2_gap,
    permutation_p_value = perm_p_value,
    permutation_n = perm_n_completed,
    permutation_null_mean_Q2 = perm_null_mean,
    pct_deviance_explained = dev_explained * 100,
    intercept = intercept,
    nonzero_predictors = list(nonzero_predictor_names),
    nonzero_coefficients = list(nonzero_predictor_coefs),
    covariate_coefficients = list(covariate_nonzero),
    all_predictor_names = list(predictor_names),
    all_coefficients = list(coef_vector),
    alpha_per_fold_nested = list(alpha_per_fold),
    n_features_per_fold = list(oof_results$n_features_per_fold),
    col_miss_threshold = col_miss_thresh,
    row_miss_threshold = row_miss_thresh,
    imputation_method = "kNN",
    seed = seed,
    date_run = as.character(Sys.time()),
    cv_fit_object = list(best_cv_fit),
    all_cv_results = list(cv_results)
  )
}


# ============================================
# 7. MODEL RUNNER FUNCTIONS WITH INCREMENTAL SAVING
# ============================================

run_score_specific_models <- function(score_specific_datasets,
                                      alpha_grid = seq(0, 1, by = 0.1),
                                      nfolds = 10, n_permutations = 100,
                                      covariate_cols = COVARIATE_COLS,
                                      seed = 42,
                                      output_file = "elastic_net_score_specific.qs2") {
  results_list <- list()
  n_models <- length(score_specific_datasets)
  
  incremental_dir <- save_output("incremental_saves_specific")
  dir.create(incremental_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("==========================================================\n")
  cat("ELASTIC NET (COV-ADJUSTED): SCORE-SPECIFIC DATASETS\n")
  cat(sprintf("Models: %d | Perms: %d | Covariates: %s\n",
              n_models, n_permutations, paste(covariate_cols, collapse = ", ")))
  cat("==========================================================\n")
  
  for (i in seq_along(score_specific_datasets)) {
    cvd_score <- names(score_specific_datasets)[i]
    df <- score_specific_datasets[[i]]
    dataset_name <- paste0(cvd_score, "_specific")
    
    cat(sprintf("\n[%d/%d] %s ~ %s\n", i, n_models, cvd_score, dataset_name))
    
    incremental_file <- file.path(incremental_dir, sprintf("model_%02d_%s.rds", i, cvd_score))
    if (file.exists(incremental_file)) {
      cat("  Loading from incremental save\n")
      results_list[[i]] <- readRDS(incremental_file)
      next
    }
    
    tryCatch({
      result <- run_elastic_net_model(
        df_raw = df, outcome_col = cvd_score, cvd_score_name = cvd_score,
        dataset_name = dataset_name, covariate_cols = covariate_cols,
        alpha_grid = alpha_grid,
        nfolds = nfolds, n_permutations = n_permutations, seed = seed
      )
      results_list[[i]] <- result
      saveRDS(result, incremental_file)
      cat("  Success\n")
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      results_list[[i]] <- NULL
    })
  }
  
  results_df <- bind_rows(results_list)
  qs_save(results_df, file = save_output(output_file))
  saveRDS(results_df, file = save_output(gsub("\\.qs2$", ".rds", output_file)))
  cat(sprintf("\nSaved %d models to %s (+ .rds backup)\n", nrow(results_df), output_file))
  results_df
}


run_common_datasets_models <- function(datasets_list, cvd_score_names,
                                       alpha_grid = seq(0, 1, by = 0.1),
                                       nfolds = 10, n_permutations = 100,
                                       covariate_cols = COVARIATE_COLS,
                                       seed = 42,
                                       output_file = "elastic_net_common.qs2") {
  results_list <- list()
  counter <- 1
  n_models <- length(datasets_list) * length(cvd_score_names)
  
  incremental_dir <- save_output("incremental_saves_common")
  dir.create(incremental_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("==========================================================\n")
  cat("ELASTIC NET (COV-ADJUSTED): COMMON DATASETS x CVD SCORES\n")
  cat(sprintf("Models: %d | Perms: %d | Covariates: %s\n",
              n_models, n_permutations, paste(covariate_cols, collapse = ", ")))
  cat("==========================================================\n")
  
  for (dataset_name in names(datasets_list)) {
    df <- datasets_list[[dataset_name]]
    
    for (cvd_score in cvd_score_names) {
      cat(sprintf("\n[%d/%d] %s ~ %s\n", counter, n_models, cvd_score, dataset_name))
      
      incremental_file <- file.path(incremental_dir,
                                    sprintf("model_%03d_%s_%s.rds", counter, cvd_score, dataset_name))
      if (file.exists(incremental_file)) {
        cat("  Loading from incremental save\n")
        results_list[[counter]] <- readRDS(incremental_file)
        counter <- counter + 1
        next
      }
      
      if (!cvd_score %in% colnames(df)) {
        cat(sprintf("  Column '%s' not found\n", cvd_score))
        counter <- counter + 1
        next
      }
      
      tryCatch({
        result <- run_elastic_net_model(
          df_raw = df, outcome_col = cvd_score, cvd_score_name = cvd_score,
          dataset_name = dataset_name, covariate_cols = covariate_cols,
          alpha_grid = alpha_grid,
          nfolds = nfolds, n_permutations = n_permutations, seed = seed
        )
        results_list[[counter]] <- result
        saveRDS(result, incremental_file)
        cat("  Success\n")
      }, error = function(e) {
        cat(sprintf("  ERROR: %s\n", e$message))
        results_list[[counter]] <- NULL
      })
      counter <- counter + 1
    }
  }
  
  results_df <- bind_rows(results_list)
  qs_save(results_df, file = save_output(output_file))
  saveRDS(results_df, file = save_output(gsub("\\.qs2$", ".rds", output_file)))
  cat(sprintf("\nSaved %d models to %s (+ .rds backup)\n", nrow(results_df), output_file))
  results_df
}


run_score_specific_plus_blocks_models <- function(score_specific_datasets, datasets_list,
                                                  alpha_grid = seq(0, 1, by = 0.1),
                                                  nfolds = 10, n_permutations = 100,
                                                  covariate_cols = COVARIATE_COLS,
                                                  seed = 42,
                                                  output_file = "elastic_net_score_plus_blocks.qs2") {
  results_list <- list()
  counter <- 1
  n_models <- length(score_specific_datasets) * length(datasets_list)
  
  incremental_dir <- save_output("incremental_saves")
  dir.create(incremental_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("==========================================================\n")
  cat("ELASTIC NET (COV-ADJUSTED): SCORE-SPECIFIC + PREDICTOR BLOCK\n")
  cat(sprintf("Models: %d | Perms: %d | Covariates: %s\n",
              n_models, n_permutations, paste(covariate_cols, collapse = ", ")))
  cat("==========================================================\n")
  
  for (cvd_score in names(score_specific_datasets)) {
    for (block_name in names(datasets_list)) {
      dataset_name <- paste0(cvd_score, "_specific+", block_name)
      
      cat(sprintf("\n[%d/%d] %s ~ %s\n", counter, n_models, cvd_score, dataset_name))
      
      incremental_file <- file.path(incremental_dir,
                                    sprintf("model_%03d_%s_%s.rds", counter, cvd_score, block_name))
      if (file.exists(incremental_file)) {
        cat("  Loading from incremental save\n")
        results_list[[counter]] <- readRDS(incremental_file)
        counter <- counter + 1
        next
      }
      
      df_combined <- build_score_plus_block(
        cvd_score, block_name, score_specific_datasets, datasets_list
      )
      
      tryCatch({
        result <- run_elastic_net_model(
          df_raw = df_combined, outcome_col = cvd_score, cvd_score_name = cvd_score,
          dataset_name = dataset_name, covariate_cols = covariate_cols,
          alpha_grid = alpha_grid,
          nfolds = nfolds, n_permutations = n_permutations, seed = seed
        )
        results_list[[counter]] <- result
        saveRDS(result, incremental_file)
        cat("  Success\n")
      }, error = function(e) {
        cat(sprintf("  ERROR: %s\n", e$message))
        results_list[[counter]] <- NULL
      })
      counter <- counter + 1
    }
  }
  
  results_df <- bind_rows(results_list)
  qs_save(results_df, file = save_output(output_file))
  saveRDS(results_df, file = save_output(gsub("\\.qs2$", ".rds", output_file)))
  cat(sprintf("\nSaved %d models to %s (+ .rds backup)\n", nrow(results_df), output_file))
  results_df
}

# ============================================
# 8. RESULTS VIEWER
# ============================================

view_results_summary <- function(results_df) {
  cat("\n=== RESULTS SUMMARY (COVARIATE-ADJUSTED) ===\n\n")
  
  summary_table <- results_df %>%
    select(cvd_score, dataset_name, n_observations, n_predictors,
           n_penalized_predictors, n_covariates, n_nonzero_coefs,
           alpha_optimal_fulldata, R2_Y, Q2_Y, Q2_R2_ratio,
           permutation_p_value) %>%
    mutate(across(where(is.numeric), ~round(., 3)))
  print(summary_table, n = Inf)
  
  cat("\n=== NON-ZERO PREDICTORS PER MODEL ===\n\n")
  for (i in 1:nrow(results_df)) {
    cat(sprintf("%s ~ %s:\n", results_df$cvd_score[i], results_df$dataset_name[i]))
    
    # Covariate coefficients
    cov_coefs <- results_df$covariate_coefficients[[i]]
    if (length(cov_coefs) > 0) {
      cat("  [Covariates - unpenalized]\n")
      for (j in seq_along(cov_coefs))
        cat(sprintf("    %s (β = %.4f)\n", names(cov_coefs)[j], cov_coefs[j]))
    }
    
    # All non-zero penalized predictors
    nz_preds <- results_df$nonzero_predictors[[i]]
    nz_coefs <- results_df$nonzero_coefficients[[i]]
    if (!is.na(nz_preds[1])) {
      cat(sprintf("  [Penalized - %d non-zero]\n", length(nz_preds)))
      for (j in seq_along(nz_preds))
        cat(sprintf("    %d. %s (β = %.4f)\n", j, nz_preds[j], nz_coefs[j]))
    } else {
      cat("  No penalized predictors selected\n")
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
  nfolds = 5, n_permutations = 1000, seed = 42,
  output_file = "elastic_net_score_specific_5fold.qs2"
)

results_common <- run_common_datasets_models(
  datasets_list = datasets_list,
  cvd_score_names = CVD_scores,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 5, n_permutations = 1000, seed = 42,
  covariate_cols = COVARIATE_COLS,
  output_file = "elastic_net_common_5fold.qs2"
)

results_score_plus_blocks <- run_score_specific_plus_blocks_models(
  score_specific_datasets = score_specific_datasets,
  datasets_list = datasets_list,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 5, n_permutations = 0, seed = 42,
  output_file = "elastic_net_score_plus_blocks_5fold.qs2"
)

# Body Composition + Height/Weight
results_body_comp_extended <- list()
for (i in seq_along(scores_without_ht_wt)) {
  score <- scores_without_ht_wt[i]
  cat(sprintf("\n[%d/3] %s ~ body_comp_with_ht_wt\n", i, score))
  df_model <- df_body_comp_extended %>% filter(!is.na(.data[[score]]))
  tryCatch({
    results_body_comp_extended[[i]] <- run_elastic_net_model(
      df_raw = df_model, outcome_col = score, cvd_score_name = score,
      dataset_name = "body_comp_with_ht_wt", covariate_cols = COVARIATE_COLS_COMMON,
      alpha_grid = seq(0, 1, by = 0.1), nfolds = 5, n_permutations = 1000, seed = 42
    )
    cat("  ✓ Success\n")
  }, error = function(e) cat(sprintf("  ✗ ERROR: %s\n", e$message)))
}
results_body_comp_extended_df <- bind_rows(results_body_comp_extended)
qs_save(results_body_comp_extended_df, save_output("elastic_net_body_comp_extended_5fold.qs2"))

# Score-specific + body comp with height/weight
results_score_plus_body_extended <- list()
for (i in seq_along(scores_without_ht_wt)) {
  score <- scores_without_ht_wt[i]
  dataset_name <- paste0(score, "_specific+body_comp_with_ht_wt")
  cat(sprintf("\n[%d/3] %s ~ %s\n", i, score, dataset_name))
  
  base_df <- score_specific_datasets[[score]]
  body_comp_predictors <- df_body_comp_extended %>%
    select(-all_of(CVD_scores), -any_of(COVARIATE_COLS_COMMON))
  
  df_combined <- base_df %>%
    left_join(body_comp_predictors, by = "Sample_ID") %>%
    filter(!is.na(.data[[score]]))
  
  tryCatch({
    results_score_plus_body_extended[[i]] <- run_elastic_net_model(
      df_raw = df_combined, outcome_col = score, cvd_score_name = score,
      dataset_name = dataset_name,
      alpha_grid = seq(0, 1, by = 0.1), nfolds = 5, n_permutations = 1000, seed = 42
    )
    cat("  ✓ Success\n")
  }, error = function(e) cat(sprintf("  ✗ ERROR: %s\n", e$message)))
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
    n_observations, n_predictors, n_penalized_predictors, n_covariates, covariates_used,
    n_nonzero_coefs, n_nonzero_penalized, n_nonzero_covariates,
    alpha_selection_method, alpha_optimal_fulldata, alpha_mean_nested,
    lambda_min, lambda_1se, lambda_selection_method,
    R2_Y, R2_fold_mean, R2_fold_sd, in_sample_rmse, in_sample_mae, in_sample_cor,
    Q2_Y, Q2_fold_mean, Q2_fold_sd, Q2_fold_summary,
    Q2_Y_rmse, Q2_Y_mae, Q2_Y_cor,
    Q2_R2_ratio, R2_Q2_gap,
    permutation_p_value, permutation_n, permutation_null_mean_Q2,
    cv_rmse_min, cv_rmse_1se, pct_deviance_explained,
    imputation_method, col_miss_threshold, row_miss_threshold,
    n_folds, seed, date_run
  )

# Predictors sheet: all non-zero penalized + covariates
predictor_rows <- list()
for (i in 1:nrow(results_all)) {
  # Penalized non-zero predictors
  preds <- results_all$nonzero_predictors[[i]]
  coefs <- results_all$nonzero_coefficients[[i]]
  if (!is.null(preds) && length(preds) > 0 && !all(is.na(preds))) {
    for (j in seq_along(preds)) {
      if (!is.na(preds[j])) {
        predictor_rows[[length(predictor_rows) + 1]] <- tibble(
          cvd_score = results_all$cvd_score[i],
          dataset_name = results_all$dataset_name[i],
          rank = j,
          predictor_name = preds[j],
          coefficient = round(coefs[j], 4),
          type = "penalized"
        )
      }
    }
  }
  # Covariate coefficients
  cov_coefs <- results_all$covariate_coefficients[[i]]
  if (!is.null(cov_coefs) && length(cov_coefs) > 0) {
    for (j in seq_along(cov_coefs)) {
      predictor_rows[[length(predictor_rows) + 1]] <- tibble(
        cvd_score = results_all$cvd_score[i],
        dataset_name = results_all$dataset_name[i],
        rank = NA_integer_,
        predictor_name = names(cov_coefs)[j],
        coefficient = round(cov_coefs[j], 4),
        type = "covariate (unpenalized)"
      )
    }
  }
}
top_predictors_sheet <- bind_rows(predictor_rows)

# Q²_Y overview tables
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
  arrange(factor(dataset_group, levels = c("specific", "block_only", "specific+block")),
          dataset_name) %>%
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

# Build workbook
headerStyle <- createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
wb <- createWorkbook()

addWorksheet(wb, "All_Results")
writeData(wb, "All_Results", results_main)
addStyle(wb, "All_Results", headerStyle, rows = 1, cols = 1:ncol(results_main), gridExpand = TRUE)
freezePane(wb, "All_Results", firstRow = TRUE)
setColWidths(wb, "All_Results", cols = 1:ncol(results_main), widths = "auto")

addWorksheet(wb, "Nonzero_Predictors")
writeData(wb, "Nonzero_Predictors", top_predictors_sheet)
addStyle(wb, "Nonzero_Predictors", headerStyle, rows = 1, cols = 1:ncol(top_predictors_sheet), gridExpand = TRUE)
freezePane(wb, "Nonzero_Predictors", firstRow = TRUE)
setColWidths(wb, "Nonzero_Predictors", cols = 1:ncol(top_predictors_sheet), widths = "auto")

addWorksheet(wb, "Q2_Y_Overview")
writeData(wb, "Q2_Y_Overview", "Table A: Q²_Y (all models) — Covariate-Adjusted",
          startRow = 1, startCol = 1)
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

saveWorkbook(wb, save_output("elastic_net_results_cov_adjusted.xlsx"), overwrite = TRUE)
cat("\n✓ Excel saved: elastic_net_results_cov_adjusted.xlsx\n")

# ============================================
# 12. TABLE CREATION FUNCTIONS
# ============================================

cvd_score_order <- c("ascvd_10y", "frs_10y", "QRISK3_risk", "SCORE2_score")

cvd_score_labels <- c(
  "SCORE2_score" = "SCORE2",
  "QRISK3_risk" = "QRISK3",
  "frs_10y" = "Framingham",
  "ascvd_10y" = "ASCVD"
)

dataset_labels_table <- c(
  "all_data" = "Full predictor set",
  "body_composition" = "Body composition",
  "sociodemographics_lifestyle" = "Sociodemographics & lifestyle",
  "clinical_risk_factors" = "Clinical Measurements & Supplementary Biomarkers",
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
  
  formatted_table <- table_data %>%
    mutate(
      R2 = ifelse(!is.na(R2_fold_sd) & R2_fold_sd > 0,
                  sprintf("%.3f \u00B1 %.3f", R2_fold_mean, R2_fold_sd),
                  sprintf("%.3f", R2_Y)),
      Q2 = ifelse(!is.na(Q2_fold_sd) & Q2_fold_sd > 0,
                  sprintf("%.3f \u00B1 %.3f", Q2_Y, Q2_fold_sd),
                  sprintf("%.3f", Q2_Y)),
      MAE = ifelse(!is.na(MAE_fold_sd) & MAE_fold_sd > 0,
                   sprintf("%.3f \u00B1 %.3f", Q2_Y_mae, MAE_fold_sd),
                   sprintf("%.3f", Q2_Y_mae)),
      perm_p = ifelse(is.na(permutation_p_value), "\u2014",
                      ifelse(permutation_p_value < 0.001, "<0.001",
                             sprintf("%.3f", permutation_p_value)))
    ) %>%
    select(CVD_Score, R2, Q2, MAE, perm_p)
  
  formatted_table %>%
    flextable() %>%
    set_caption(caption = table_title) %>%
    set_header_labels(
      CVD_Score = "CVD Score",
      R2 = "R\u00B2 \u00B1 SD",
      Q2 = "Q\u00B2 \u00B1 SD",
      MAE = "MAE \u00B1 SD",
      perm_p = "Perm. p"
    ) %>%
    align(j = 1, align = "left", part = "all") %>%
    align(j = 2:5, align = "center", part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 9, part = "body") %>%
    flextable::fontsize(size = 10, part = "header") %>%
    bold(part = "header") %>%
    width(j = 1, width = 1.2) %>%
    width(j = 2:4, width = 1.3) %>%
    width(j = 5, width = 1.0) %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1.5), part = "header") %>%
    padding(padding = 3, part = "all")
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
  
  if (nrow(table_data) == 0) { warning("No data!"); return(NULL) }
  table_data <- table_data %>% arrange(dataset_name, cvd_score)
  
  formatted_table <- table_data %>%
    mutate(
      Dataset = ifelse(dataset_name %in% names(dataset_labels),
                       dataset_labels[as.character(dataset_name)],
                       as.character(dataset_name)),
      R2 = ifelse(!is.na(R2_fold_sd) & R2_fold_sd > 0,
                  sprintf("%.3f \u00B1 %.3f", R2_fold_mean, R2_fold_sd),
                  sprintf("%.3f", R2_Y)),
      Q2 = ifelse(!is.na(Q2_fold_sd) & Q2_fold_sd > 0,
                  sprintf("%.3f \u00B1 %.3f", Q2_Y, Q2_fold_sd),
                  sprintf("%.3f", Q2_Y)),
      MAE = ifelse(!is.na(MAE_fold_sd) & MAE_fold_sd > 0,
                   sprintf("%.3f \u00B1 %.3f", Q2_Y_mae, MAE_fold_sd),
                   sprintf("%.3f", Q2_Y_mae)),
      perm_p = ifelse(is.na(permutation_p_value), "\u2014",
                      ifelse(permutation_p_value < 0.001, "<0.001",
                             sprintf("%.3f", permutation_p_value)))
    ) %>%
    select(Dataset, CVD_Score, R2, Q2, MAE, perm_p)
  
  n_rows <- nrow(formatted_table)
  
  ft <- formatted_table %>%
    flextable() %>%
    merge_v(j = "Dataset") %>%
    set_caption(caption = table_title) %>%
    set_header_labels(
      CVD_Score = "CVD Score",
      R2 = "R\u00B2 \u00B1 SD",
      Q2 = "Q\u00B2 \u00B1 SD",
      MAE = "MAE \u00B1 SD",
      perm_p = "Perm. p"
    ) %>%
    align(j = 1:2, align = "left", part = "all") %>%
    align(j = 3:6, align = "center", part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 9, part = "body") %>%
    flextable::fontsize(size = 10, part = "header") %>%
    bold(part = "header") %>%
    bold(j = 1, part = "body") %>%
    width(j = 1, width = 1.8) %>%
    width(j = 2, width = 1.0) %>%
    width(j = 3:5, width = 1.3) %>%
    width(j = 6, width = 1.0) %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1.5), part = "header") %>%
    padding(padding = 3, part = "all")
  
  if (n_rows > 4) {
    hline_rows <- seq(4, n_rows - 1, by = 4)
    if (length(hline_rows) > 0) {
      ft <- ft %>%
        hline(i = hline_rows,
              border = fp_border(width = 0.5, color = "#CCCCCC"), part = "body")
    }
  }
  ft
}


# ============================================
# CREATE AND SAVE TABLES
# ============================================

cat("\n==========================================================\n")
cat("CREATING PUBLICATION TABLES\n")
cat("==========================================================\n")

table1 <- create_score_specific_table(
  results_score_specific,
  "Table 1: Elastic Net – Score-Specific (Covariate-Adjusted)"
)

datasets_table2 <- c("all_data", "lipids", "fatty_acids", "urine_nmr",
                     "body_composition", "sociodemographics_lifestyle", "clinical_risk_factors")
datasets_table2 <- datasets_table2[datasets_table2 %in% unique(results_common$dataset_name)]

table2 <- create_multi_score_table(
  results_common %>% filter(dataset_name %in% datasets_table2),
  "Table 2: Elastic Net – Predictor Blocks (Covariate-Adjusted)",
  dataset_order = datasets_table2
)

blocks_table3 <- datasets_table2
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
  results_score_plus_blocks %>%
    mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
    filter(block %in% blocks_table3),
  "Table 3: Elastic Net – Score-Specific + Blocks (Covariate-Adjusted)",
  dataset_order = dataset_order_table3,
  dataset_labels = dataset_labels_table3
)

# Save tables
save_as_docx(table1, path = save_output("Table1_score_specific.docx"))
save_as_docx(table2, path = save_output("Table2_predictor_blocks.docx"))
save_as_docx(table3, path = save_output("Table3_score_plus_blocks.docx"))

save_as_html(table1, path = save_output("Table1_score_specific.html"))
save_as_html(table2, path = save_output("Table2_predictor_blocks.html"))
save_as_html(table3, path = save_output("Table3_score_plus_blocks.html"))

# Combined document
doc <- read_docx() %>%
  body_add_par("Elastic Net Results — Covariate-Adjusted", style = "heading 1") %>%
  body_add_par(sprintf("Covariates (unpenalized): %s", paste(COVARIATE_COLS, collapse = ", "))) %>%
  body_add_par("") %>%
  body_add_flextable(value = table1) %>%
  body_add_break() %>%
  body_add_flextable(value = table2) %>%
  body_add_break() %>%
  body_add_flextable(value = table3)

print(doc, target = save_output("Elastic_Net_All_Tables_CovAdj.docx"))
cat("\n✓ All tables saved\n")





# ============================================
# 14. GROUPED-BY-SCORE PERFORMANCE PLOT COMPARING SCORE-SPECIFIC WITH SCORE-SPECIFIC + DATA SETS
# ============================================

cat("\n==========================================================\n")
cat("CREATING GROUPED-BY-SCORE PERFORMANCE PLOT\n")
cat("==========================================================\n")

# --- Define predictor set order and labels ---
block_order <- c("lipids", "fatty_acids", "urine_nmr",
                 "body_composition", "clinical_risk_factors",
                 "sociodemographics_lifestyle")
block_order <- block_order[block_order %in% unique(results_common$dataset_name)]

model_type_order <- c(
  "score_specific",
  "score_specific+all_data",
  paste0("score_specific+", block_order)
)

model_type_labels <- c(
  "score_specific" = "Score-specific only",
  "score_specific+all_data" = "Score inputs + Full set",
  "score_specific+lipids" = "Score inputs + Lipids",
  "score_specific+fatty_acids" = "Score inputs + Fatty acids",
  "score_specific+urine_nmr" = "Score inputs + Urine NMR",
  "score_specific+body_composition" = "Score inputs + Body composition",
  "score_specific+clinical_risk_factors" = "Score inputs + Clinical Measurements & Suppl. Biomarkers",
  "score_specific+sociodemographics_lifestyle" = "Score inputs + Sociodemographics & Lifestyle"
)

# --- Predictor set colours (avoiding CVD score colours) ---
predictor_colors <- c(
  "Score-specific only" = "#888888",
  "Score inputs + Full set" = "#2D2D2D",
  "Score inputs + Lipids" = "#E63946",
  "Score inputs + Fatty acids" = "#F4A261",
  "Score inputs + Urine NMR" = "#E9C46A",
  "Score inputs + Body composition" = "#457B9D",
  "Score inputs + Clinical Measurements & Suppl. Biomarkers" = "#2A9D8F",
  "Score inputs + Sociodemographics & Lifestyle" = "#7B2D8E"
)

# --- Assemble data ---
d_score_specific <- results_score_specific %>%
  filter(cvd_score %in% cvd_score_order) %>%
  mutate(model_type = "score_specific")

d_score_plus <- results_score_plus_blocks %>%
  filter(cvd_score %in% cvd_score_order) %>%
  mutate(
    block = sub("^.*_specific\\+", "", dataset_name),
    model_type = paste0("score_specific+", block)
  )

plot_data_grouped <- bind_rows(d_score_specific, d_score_plus) %>%
  filter(model_type %in% model_type_order) %>%
  mutate(
    cvd_label = cvd_score_labels[as.character(cvd_score)],
    cvd_label = factor(cvd_label, levels = rev(cvd_score_labels)),
    model_type = factor(model_type, levels = model_type_order),
    model_label = model_type_labels[as.character(model_type)],
    model_label = factor(model_label,
                         levels = model_type_labels[model_type_order])
  )

if (!"Q2_fold_sd" %in% names(plot_data_grouped)) plot_data_grouped$Q2_fold_sd <- 0
if (!"MAE_fold_sd" %in% names(plot_data_grouped)) plot_data_grouped$MAE_fold_sd <- 0

q2_ymin <- floor(min(plot_data_grouped$Q2_Y - plot_data_grouped$Q2_fold_sd, na.rm = TRUE) * 10) / 10

# --- Plot A: Q² ± SD (vertical) ---
p_q2_vertical <- plot_data_grouped %>%
  ggplot(aes(x = cvd_label, y = Q2_Y, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_errorbar(
    aes(ymin = Q2_Y - Q2_fold_sd, ymax = Q2_Y + Q2_fold_sd),
    position = position_dodge(width = 0.85), width = 0.3, linewidth = 0.3
  ) +
  geom_text(
    aes(y = pmax(Q2_Y, Q2_Y + Q2_fold_sd), label = sprintf("%.2f", Q2_Y)),
    position = position_dodge(width = 0.85), vjust = -0.4, size = 2.1
  ) +
  scale_fill_manual(values = predictor_colors, breaks = unname(model_type_labels[model_type_order])) +
  coord_cartesian(ylim = c(q2_ymin, 1.08)) +
  scale_y_continuous(breaks = seq(q2_ymin, 1, by = 0.1)) +
  labs(title = expression(bold(paste("A) ", Q^2, " ± SD"))),
       y = expression(bold(Q^2)), x = NULL, fill = "Model") +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0, size = 15),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2, override.aes = list(size = 3)))

# --- Plot B: MAE ± SD (vertical) ---
p_mae_vertical <- plot_data_grouped %>%
  ggplot(aes(x = cvd_label, y = Q2_Y_mae, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_errorbar(
    aes(ymin = pmax(0, Q2_Y_mae - MAE_fold_sd), ymax = Q2_Y_mae + MAE_fold_sd),
    position = position_dodge(width = 0.85), width = 0.3, linewidth = 0.3
  ) +
  geom_text(
    aes(y = Q2_Y_mae + MAE_fold_sd, label = sprintf("%.2f", Q2_Y_mae)),
    position = position_dodge(width = 0.85), vjust = -0.4, size = 2.1
  ) +
  scale_fill_manual(values = predictor_colors, breaks = unname(model_type_labels[model_type_order])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = expression(bold(paste("B) MAE ± SD"))),
       y = expression(bold("MAE")), x = NULL, fill = "Model") +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0, size = 15),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2, override.aes = list(size = 3)))

# --- Combine with shared legend ---
combined_vertical_plot <- (p_q2_vertical | p_mae_vertical) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.justification = "center")

print(combined_vertical_plot)

ggsave(save_output("elastic_net_grouped_by_score_vertical.png"),
       plot = combined_vertical_plot, width = 14, height = 7, dpi = 300, bg = "white")
ggsave(save_output("elastic_net_grouped_by_score_vertical.pdf"),
       plot = combined_vertical_plot, width = 14, height = 7, device = "pdf")

cat("\n✓ Vertical grouped-by-score plots saved\n")






# ============================================
# TUNING DIAGNOSTIC PLOTS — FROM ACTUAL MODEL FITS
# ============================================

library(tidyverse)
library(glmnet)
library(patchwork)

# ============================================
# 1. EXTRACT TUNING DATA FROM STORED RESULTS
# ============================================

extract_tuning_from_results <- function(results_row, dataset_display_name = NULL) {
  
  cv_results <- results_row$all_cv_results[[1]]
  best_cv_fit <- results_row$cv_fit_object[[1]]
  best_alpha <- results_row$alpha_optimal_fulldata
  n_pred <- results_row$n_predictors
  n_obs <- results_row$n_observations
  nfolds <- results_row$n_folds
  
  if (is.null(dataset_display_name)) {
    dataset_display_name <- results_row$dataset_name
  }
  
  # Alpha panel data
  alpha_grid <- sapply(cv_results, function(x) x$alpha)
  alpha_rmse_df <- tibble(
    alpha = alpha_grid,
    cv_rmse = sapply(cv_results, function(x) {
      idx <- which(x$cv_fit$lambda == x$cv_fit$lambda.1se)
      sqrt(x$cv_fit$cvm[idx])
    }),
    cv_rmse_se = sapply(cv_results, function(x) {
      idx <- which(x$cv_fit$lambda == x$cv_fit$lambda.1se)
      mse_sd <- x$cv_fit$cvsd[idx]
      rmse_val <- sqrt(x$cv_fit$cvm[idx])
      (mse_sd / sqrt(nfolds)) / (2 * rmse_val)
    })
  )
  
  best_alpha_idx <- which.min(alpha_rmse_df$cv_rmse)
  
  # Lambda panel data
  lambda_df <- tibble(
    log_lambda = log(best_cv_fit$lambda),
    cv_rmse = sqrt(best_cv_fit$cvm),
    cv_rmse_se = sapply(seq_along(best_cv_fit$lambda), function(i) {
      mse_sd <- best_cv_fit$cvsd[i]
      rmse_val <- sqrt(best_cv_fit$cvm[i])
      (mse_sd / sqrt(nfolds)) / (2 * rmse_val)
    }),
    nzero = best_cv_fit$nzero
  )
  
  nzero_at_min <- best_cv_fit$nzero[which(best_cv_fit$lambda == best_cv_fit$lambda.min)]
  nzero_at_1se <- best_cv_fit$nzero[which(best_cv_fit$lambda == best_cv_fit$lambda.1se)]
  
  list(
    alpha_rmse_df    = alpha_rmse_df,
    best_alpha_idx   = best_alpha_idx,
    best_alpha       = best_alpha,
    best_cv_fit      = best_cv_fit,
    lambda_df        = lambda_df,
    nzero_at_min     = nzero_at_min,
    nzero_at_1se     = nzero_at_1se,
    n_raw_predictors = n_pred,
    n_obs            = n_obs,
    dataset_name     = dataset_display_name
  )
}

# ============================================
# 2. PANEL BUILDERS (unchanged)
# ============================================

build_alpha_panel <- function(res, show_y_label = TRUE) {
  df <- res$alpha_rmse_df
  best_idx <- res$best_alpha_idx
  selected <- df[best_idx, ]
  
  p <- ggplot(df, aes(x = alpha, y = cv_rmse)) +
    geom_line(colour = "grey60", linewidth = 0.5) +
    geom_point(colour = "grey40", size = 2, alpha = 0.7) +
    geom_errorbar(aes(ymin = cv_rmse - cv_rmse_se, ymax = cv_rmse + cv_rmse_se),
                  colour = "grey70", alpha = 0.5, width = 0) +
    geom_point(data = selected, colour = "#E41A1C", size = 3.5, shape = 18) +
    geom_vline(xintercept = selected$alpha,
               colour = "#E41A1C", linetype = "dashed", alpha = 0.5) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    labs(
      title = res$dataset_name,
      subtitle = sprintf("p = %d, n = %d", res$n_raw_predictors, res$n_obs),
      x = expression(alpha ~ "(mixing parameter)"),
      y = if (show_y_label) "CV RMSE" else NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 11, colour = "grey40", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      panel.grid.minor = element_blank()
    )
  p
}

build_lambda_panel <- function(res, show_y_label = TRUE) {
  df <- res$lambda_df
  cv_fit <- res$best_cv_fit
  
  lambda_min_val <- log(cv_fit$lambda.min)
  lambda_1se_val <- log(cv_fit$lambda.1se)
  
  sel_idx <- which.min(abs(df$log_lambda - lambda_1se_val))
  selected <- df[sel_idx, ]
  
  n_breaks <- min(6, nrow(df))
  break_idx <- seq(1, nrow(df), length.out = n_breaks) %>% round()
  
  p <- ggplot(df, aes(x = log_lambda, y = cv_rmse)) +
    geom_line(colour = "grey60", linewidth = 0.5) +
    geom_errorbar(aes(ymin = cv_rmse - cv_rmse_se, ymax = cv_rmse + cv_rmse_se),
                  colour = "grey70", alpha = 0.3, width = 0) +
    geom_vline(xintercept = lambda_min_val, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = lambda_1se_val, linetype = "dashed", colour = "#E41A1C", alpha = 0.5, linewidth = 0.4) +
    geom_point(data = selected, colour = "#E41A1C", size = 3.5, shape = 18) +
    annotate("text", x = lambda_min_val, y = max(df$cv_rmse) * 0.98,
             label = sprintf("p==%d", res$nzero_at_min),
             parse = TRUE, hjust = 1.1, size = 2.5, colour = "grey40") +
    annotate("text", x = lambda_1se_val, y = max(df$cv_rmse) * 0.98,
             label = sprintf("p==%d", res$nzero_at_1se),
             parse = TRUE, hjust = -0.1, size = 2.5, colour = "#E41A1C") +
    scale_x_continuous(
      sec.axis = sec_axis(
        transform = ~ .,
        breaks = df$log_lambda[break_idx],
        labels = df$nzero[break_idx]
      )
    ) +
    labs(
      x = expression(log(lambda)),
      y = if (show_y_label) "CV RMSE" else NULL,
      subtitle = sprintf("\u03b1 = %.1f  |  \u03bb.1se: %d coefs", res$best_alpha, res$nzero_at_1se)
    ) +
    theme_minimal(base_size = 11) +
    theme(
      text = element_text(family = "Arial"),
      plot.subtitle = element_text(size = 11, colour = "grey40", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      axis.title.x.top = element_blank(),
      axis.text.x.top = element_text(size = 7),
      panel.grid.minor = element_blank()
    )
  p
}

# ============================================
# 3. ASSEMBLE FIGURE
# ============================================

assemble_tuning_figure <- function(tuning_results, add_title = TRUE) {
  n <- length(tuning_results)
  
  alpha_panels <- lapply(seq_len(n), function(i) {
    build_alpha_panel(tuning_results[[i]], show_y_label = (i == 1))
  })
  lambda_panels <- lapply(seq_len(n), function(i) {
    build_lambda_panel(tuning_results[[i]], show_y_label = (i == 1))
  })
  
  alpha_panels[[1]] <- alpha_panels[[1]] +
    labs(tag = "A") +
    theme(plot.tag = element_text(face = "bold", family = "Arial", size = 14))
  lambda_panels[[1]] <- lambda_panels[[1]] +
    labs(tag = "B") +
    theme(plot.tag = element_text(face = "bold", family = "Arial", size = 14))
  
  alpha_row <- wrap_plots(alpha_panels, nrow = 1)
  lambda_row <- wrap_plots(lambda_panels, nrow = 1)
  combined <- alpha_row / lambda_row
  
  if (add_title) {
    combined <- combined +
      plot_annotation(
        title = "Elastic Net Hyperparameter Tuning: QRISK3",
        subtitle = "A: \u03b1 selected by minimum CV RMSE | B: \u03bb selected by one-SE rule at chosen \u03b1",
        theme = theme(
          plot.title = element_text(face = "bold", family = "Arial", size = 14),
          plot.subtitle = element_text(family = "Arial", size = 10, colour = "grey40")
        )
      )
  }
  combined
}

# ============================================
# 4. EXTRACT FROM ACTUAL RESULTS
# ============================================

# Dataset display names (same order as before)
dataset_display_names <- list(
  "QRISK3_risk_specific" = "Score-Specific",
  "all_data" = "Full Predictor Set",
  "body_composition" = "Body Composition",
  "sociodemographics_lifestyle" = "Sociodemographics",
  "clinical_risk_factors" = "Clinical Risk Factors",
  "lipids" = "Lipids",
  "fatty_acids" = "Fatty Acids",
  "urine_nmr" = "Urine NMR"
)

# Filter QRISK3 results from score-specific and common
qrisk3_score_specific <- results_score_specific %>%
  filter(cvd_score == "QRISK3_risk")

qrisk3_common <- results_common %>%
  filter(cvd_score == "QRISK3_risk")

tuning_configs <- list(
  list(row = qrisk3_score_specific, name = "Score-Specific"),
  list(row = qrisk3_common %>% filter(dataset_name == "all_data"), name = "Full Predictor Set"),
  list(row = qrisk3_common %>% filter(dataset_name == "body_composition"), name = "Body Composition"),
  list(row = qrisk3_common %>% filter(dataset_name == "sociodemographics_lifestyle"), name = "Sociodemographics"),
  list(row = qrisk3_common %>% filter(dataset_name == "clinical_risk_factors"), name = "Clinical Risk Factors"),
  list(row = qrisk3_common %>% filter(dataset_name == "lipids"), name = "Lipids"),
  list(row = qrisk3_common %>% filter(dataset_name == "fatty_acids"), name = "Fatty Acids"),
  list(row = qrisk3_common %>% filter(dataset_name == "urine_nmr"), name = "Urine NMR")
)

tuning_results_all <- lapply(tuning_configs, function(cfg) {
  cat(sprintf("Extracting tuning data for QRISK3 ~ %s...\n", cfg$name))
  extract_tuning_from_results(cfg$row, dataset_display_name = cfg$name)
})

# ============================================
# 5. PLOT 1: ALL 8 DATASETS (exploratory)
# ============================================

plot_all <- assemble_tuning_figure(tuning_results_all, add_title = TRUE)

ggsave(save_output("tuning_diagnostics_QRISK3_all.png"),
       plot = plot_all, width = 24, height = 9, dpi = 300)
ggsave(save_output("tuning_diagnostics_QRISK3_all.pdf"),
       plot = plot_all, width = 24, height = 9, dpi = 300, device = cairo_pdf)
cat("Saved: tuning_diagnostics_QRISK3_all (.png + .pdf)\n")

# ============================================
# 6. PLOT 2: 5 SELECTED DATASETS (thesis)
# ============================================

thesis_indices <- c(1, 2, 3, 4, 6)
tuning_results_thesis <- tuning_results_all[thesis_indices]

plot_thesis <- assemble_tuning_figure(tuning_results_thesis, add_title = FALSE)

ggsave(save_output("tuning_diagnostics_QRISK3_thesis.png"),
       plot = plot_thesis, width = 16, height = 7, dpi = 300)
ggsave(save_output("tuning_diagnostics_QRISK3_thesis.pdf"),
       plot = plot_thesis, width = 16, height = 7, dpi = 300, device = cairo_pdf)
cat("Saved: tuning_diagnostics_QRISK3_thesis (.png + .pdf)\n")








# ============================================
# 11. VISUALIZE NON-ZERO COEFFICIENTS
# ============================================

extract_nonzero_coefs <- function(results_df) {
  coef_list <- list()
  for (i in 1:nrow(results_df)) {
    pnames <- results_df$all_predictor_names[[i]]
    coeffs <- results_df$all_coefficients[[i]]
    nz <- coeffs != 0
    if (sum(nz) > 0) {
      is_cov <- grepl(paste0("^(", paste(COVARIATE_COLS, collapse = "|"), ")"), pnames[nz])
      coef_list[[i]] <- tibble(
        cvd_score = results_df$cvd_score[i],
        dataset_name = results_df$dataset_name[i],
        predictor = pnames[nz],
        coefficient = coeffs[nz],
        abs_coefficient = abs(coeffs[nz]),
        direction = ifelse(coeffs[nz] > 0, "Positive", "Negative"),
        predictor_type = ifelse(is_cov, "covariate", "penalized")
      )
    }
  }
  bind_rows(coef_list)
}


plot_nonzero_coefs <- function(coef_data, cvd_score, dataset_name) {
  # Only plot penalized predictors (exclude covariates)
  coef_data <- coef_data %>%
    filter(predictor_type == "penalized") %>%
    arrange(abs_coefficient) %>%
    mutate(predictor = factor(predictor, levels = predictor))
  
  if (nrow(coef_data) == 0) return(NULL)
  
  ggplot(coef_data, aes(x = predictor, y = abs_coefficient, fill = direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = c("Positive" = "#FF8C00", "Negative" = "#4169E1"),
                      name = "Direction") +
    coord_flip() +
    labs(
      title = paste0("Non-Zero Penalized Coefficients: ", cvd_score),
      subtitle = paste0("Dataset: ", dataset_name,
                        " (adjusted for ", paste(COVARIATE_COLS, collapse = ", "), ")"),
      x = "Predictor", y = "Absolute Coefficient Value"
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
}


create_coefficient_plots <- function(results_df, output_dir = "coefficient_plots") {
  full_output_dir <- save_output(output_dir)
  dir.create(full_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\nExtracting non-zero coefficients...\n")
  all_coefs <- extract_nonzero_coefs(results_df)
  
  models <- results_df %>% select(cvd_score, dataset_name) %>% distinct()
  cat(sprintf("Creating %d coefficient plots...\n", nrow(models)))
  
  plot_list <- list()
  for (i in 1:nrow(models)) {
    cvd <- models$cvd_score[i]
    dataset <- models$dataset_name[i]
    model_coefs <- all_coefs %>%
      filter(cvd_score == cvd, dataset_name == dataset, predictor_type == "penalized")
    
    if (nrow(model_coefs) == 0) {
      cat(sprintf("  [%d/%d] %s ~ %s: No non-zero coefficients\n",
                  i, nrow(models), cvd, dataset))
      next
    }
    
    cat(sprintf("  [%d/%d] %s ~ %s: %d non-zero coefficients\n",
                i, nrow(models), cvd, dataset, nrow(model_coefs)))
    
    p <- plot_nonzero_coefs(model_coefs, cvd, dataset)
    plot_list[[i]] <- p
    
    filename <- file.path(full_output_dir,
                          paste0("coef_", cvd, "_", gsub("[^A-Za-z0-9_]", "_", dataset), ".pdf"))
    plot_height <- max(6, min(20, 3 + 0.3 * nrow(model_coefs)))
    ggsave(filename, plot = p, width = 10, height = plot_height, device = "pdf")
  }
  
  cat(sprintf("\n✓ Saved %d plots to %s/\n",
              length(plot_list[!sapply(plot_list, is.null)]), output_dir))
  list(plots = plot_list, coefficients = all_coefs)
}


smart_truncate <- function(x, max_len = 30) {
  ifelse(nchar(x) > max_len,
         paste0(substr(x, 1, 12), "..", substr(x, nchar(x) - 12, nchar(x))),
         x)
}


plot_block_coefficients <- function(coef_data, block_name) {
  n_predictors <- n_distinct(coef_data$predictor)
  x_text_size <- case_when(
    block_name == "all_data" ~ 5,
    n_predictors > 30 ~ 6,
    TRUE ~ 7
  )
  
  coef_data <- coef_data %>%
    mutate(
      cvd_score_label = recode(cvd_score,
                               "ascvd_10y" = "ASCVD",
                               "frs_10y" = "Framingham",
                               "QRISK3_risk" = "QRISK3",
                               "SCORE2_score" = "SCORE2"),
      direction = ifelse(coefficient > 0, "Positive", "Negative"),
      abs_coef = abs(coefficient),
      predictor_short = smart_truncate(predictor, max_len = 30),
      predictor_id = paste(predictor_short, cvd_score_label, sep = "___")
    ) %>%
    arrange(cvd_score_label, desc(abs_coef))
  
  coef_data$predictor_id <- factor(coef_data$predictor_id, levels = unique(coef_data$predictor_id))
  coef_data$cvd_score_label <- factor(coef_data$cvd_score_label,
                                      levels = c("ASCVD", "Framingham", "QRISK3", "SCORE2"))
  
  ggplot(coef_data, aes(x = predictor_id, y = abs_coef, fill = direction)) +
    geom_col(width = 0.7) +
    facet_wrap(~ cvd_score_label, scales = "free_x", nrow = 1) +
    scale_x_discrete(labels = function(x) gsub("___.*$", "", x)) +
    scale_fill_manual(values = c("Positive" = "#E07B39", "Negative" = "#3B7EA1"),
                      name = "Effect Direction") +
    labs(
      title = sprintf("Non-Zero Penalized Coefficients: %s (Adjusted for %s)",
                      block_name, paste(COVARIATE_COLS, collapse = ", ")),
      x = NULL, y = "Absolute Coefficient"
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


create_block_plots <- function(results_df,
                               blocks = c("all_data", "lipids", "fatty_acids", "urine_nmr",
                                          "body_composition", "clinical_risk_factors",
                                          "sociodemographics_lifestyle"),
                               output_dir = "coefficient_plots_by_block") {
  full_output_dir <- save_output(output_dir)
  dir.create(full_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  all_coefs <- extract_nonzero_coefs(results_df) %>%
    filter(predictor_type == "penalized") %>%
    mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
    filter(block %in% blocks)
  
  cat(sprintf("\nCreating plots for %d blocks...\n", length(blocks)))
  
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
  
  cat(sprintf("✓ Saved block plots to %s/\n", output_dir))
}


# ---- Run coefficient visualization ----
cat("\n==========================================================\n")
cat("VISUALIZING NON-ZERO COEFFICIENTS\n")
cat("==========================================================\n")

coef_plots_results <- create_coefficient_plots(
  results_df = results_score_plus_blocks,
  output_dir = "coefficient_plots_score_plus_blocks"
)

coef_summary <- coef_plots_results$coefficients %>%
  group_by(cvd_score, dataset_name) %>%
  summarise(
    n_nonzero = n(),
    n_positive = sum(direction == "Positive"),
    n_negative = sum(direction == "Negative"),
    n_covariates = sum(predictor_type == "covariate"),
    n_penalized = sum(predictor_type == "penalized"),
    mean_abs_coef = mean(abs_coefficient),
    max_abs_coef = max(abs_coefficient),
    .groups = "drop"
  ) %>%
  arrange(cvd_score, desc(n_nonzero))

print(coef_summary)
write.xlsx(coef_summary, save_output("coefficient_summary.xlsx"))

create_block_plots(
  results_df = results_score_plus_blocks,
  output_dir = "coefficient_plots_by_block"
)

# ---- Block-only plots (results_common) ----
cat("\n--- Block-only coefficient plots ---\n")

coef_plots_blocks_only <- create_coefficient_plots(
  results_df = results_common,
  output_dir = "coefficient_plots_blocks_only"
)

# Faceted block-only plots: dataset_name IS the block name directly
create_block_plots_common <- function(results_df,
                                      blocks = c("all_data", "lipids", "fatty_acids", "urine_nmr",
                                                 "body_composition", "clinical_risk_factors",
                                                 "sociodemographics_lifestyle"),
                                      output_dir = "coefficient_plots_blocks_only_faceted") {
  full_output_dir <- save_output(output_dir)
  dir.create(full_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  all_coefs <- extract_nonzero_coefs(results_df) %>%
    filter(predictor_type == "penalized") %>%
    filter(dataset_name %in% blocks)
  
  cat(sprintf("\nCreating faceted plots for %d blocks (block-only models)...\n", length(blocks)))
  
  for (current_block in blocks) {
    block_data <- filter(all_coefs, dataset_name == current_block)
    if (nrow(block_data) == 0) {
      cat(sprintf("  %s: No data - skipping\n", current_block))
      next
    }
    
    n_pred <- n_distinct(block_data$predictor)
    n_scores <- n_distinct(block_data$cvd_score)
    cat(sprintf("  %s: %d unique predictors across %d CVD scores\n",
                current_block, n_pred, n_scores))
    
    p <- plot_block_coefficients(block_data, current_block)
    filename <- file.path(full_output_dir, sprintf("coef_block_only_%s.pdf", current_block))
    plot_width <- max(12, min(22, 6 + 0.18 * n_pred))
    ggsave(filename, p, width = plot_width, height = 6)
  }
  
  cat(sprintf("✓ Saved block-only faceted plots to %s/\n", output_dir))
}

create_block_plots_common(
  results_df = results_common,
  output_dir = "coefficient_plots_blocks_only_faceted"
)

# ============================================
# 13. PERFORMANCE COMPARISON PLOTS
# ============================================

cat("\n==========================================================\n")
cat("CREATING PERFORMANCE COMPARISON PLOTS\n")
cat("==========================================================\n")

cvd_colors <- c(
  "QRISK3" = "#F8766D",
  "SCORE2" = "#00BA38",
  "ASCVD" = "#619CFF",
  "Framingham" = "#B79F00"
)

# --- Left column: Score-specific & predictor blocks ---
score_specific_grouped <- results_score_specific %>%
  filter(cvd_score %in% cvd_score_order) %>%
  mutate(dataset_name = "score_specific", dataset_label = "Score-Specific")

datasets_for_left <- c("lipids", "fatty_acids", "urine_nmr",
                       "body_composition", "clinical_risk_factors",
                       "sociodemographics_lifestyle")
datasets_for_left <- datasets_for_left[datasets_for_left %in% unique(results_common$dataset_name)]

dataset_labels_left <- c("score_specific" = "Score-Specific",
                         dataset_labels_table[datasets_for_left])
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
    dataset_label = factor(dataset_label,
                           levels = unique(dataset_labels_left[dataset_order_left]))
  )

# --- Right column: Score-specific + predictor blocks ---
baseline_grouped <- results_score_specific %>%
  filter(cvd_score %in% cvd_score_order) %>%
  mutate(dataset_name = "score_specific", dataset_label = "Score-Specific",
         block = "score_specific")

blocks_for_plus <- datasets_for_left

score_plus_data <- results_score_plus_blocks %>%
  mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
  filter(block %in% blocks_for_plus, cvd_score %in% cvd_score_order) %>%
  mutate(
    dataset_name = block,
    dataset_label = ifelse(block %in% names(dataset_labels_table),
                           dataset_labels_table[block], block)
  )

dataset_order_right <- rev(c("score_specific", blocks_for_plus))
dataset_labels_right <- c("score_specific" = "Score-Specific",
                          dataset_labels_table[blocks_for_plus])

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
    dataset_label = factor(dataset_label,
                           levels = unique(dataset_labels_right[dataset_order_right]))
  )

# Handle missing SD columns
if (!"MAE_fold_sd" %in% names(plot_data_left)) plot_data_left$MAE_fold_sd <- 0
if (!"MAE_fold_sd" %in% names(plot_data_right)) plot_data_right$MAE_fold_sd <- 0

# --- Q² plots ---
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
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Score-Specific & Predictor Blocks",
       x = expression(Q^2), y = NULL, fill = "CVD Score") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
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
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Score-Specific + Predictor Blocks",
       x = expression(Q^2), y = NULL, fill = "CVD Score") +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold")
  )

# --- MAE plots ---
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
    axis.text = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
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
    axis.text = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold")
  )

# --- Combine ---
combined_elastic_plot <- (p1_left_q2 | p2_right_q2) / (p3_left_mae | p4_right_mae) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Elastic Net Performance: Cross-Validated Predictive Ability (Covariate-Adjusted)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

print(combined_elastic_plot)

ggsave(save_output("elastic_net_performance_comparison.png"),
       plot = combined_elastic_plot, width = 16, height = 12, dpi = 300, bg = "white")
ggsave(save_output("elastic_net_performance_comparison.pdf"),
       plot = combined_elastic_plot, width = 16, height = 12, device = "pdf")

cat("\n✓ Performance plots saved\n")

# ============================================
# FINAL SAVE
# ============================================

qs_save(results_all, save_output("elastic_net_all_results.qs2"))
saveRDS(results_all, save_output("elastic_net_all_results.rds"))

cat("\n==========================================================\n")
cat(sprintf("ALL OUTPUT FILES SAVED TO: %s\n", OUTPUT_DIR))
cat(sprintf("Covariates (unpenalized): %s\n", paste(COVARIATE_COLS, collapse = ", ")))
cat("==========================================================\n")


# ============================================
# 14. PERMUTATION FEATURE IMPORTANCE (PARALLELIZED, COV-ADJUSTED)
# ============================================
# Only tests features with non-zero PENALIZED coefficients.
# Covariates (Country, Statins, Supplements) are unpenalized and not tested.
# ============================================

#' Calculate permutation importance for a single feature
calculate_single_feature_importance <- function(feat, X_test, y_test,
                                                y_pred_baseline, mae_baseline,
                                                cv_fit, n_permutations, seed) {
  set.seed(seed)
  feat_idx <- which(colnames(X_test) == feat)
  if (length(feat_idx) == 0) {
    return(tibble(feature = feat, importance = NA_real_,
                  importance_sd = NA_real_, p_value = NA_real_))
  }
  
  mae_permuted <- numeric(n_permutations)
  for (p in 1:n_permutations) {
    X_test_perm <- X_test
    X_test_perm[, feat_idx] <- sample(X_test_perm[, feat_idx])
    y_pred_perm <- as.numeric(predict(cv_fit, newx = X_test_perm, s = "lambda.1se"))
    mae_permuted[p] <- mean(abs(y_test - y_pred_perm))
  }
  
  delta_mae <- mae_permuted - mae_baseline
  tibble(
    feature = feat,
    importance = mean(delta_mae),
    importance_sd = sd(delta_mae),
    p_value = mean(delta_mae <= 0)
  )
}


#' Calculate permutation importance for a single fold (PARALLELIZED across features)
calculate_fold_importance <- function(X_train, X_test, y_train, y_test,
                                      cv_fit, features_to_test,
                                      n_permutations = 1000,
                                      n_cores = 1, seed = 42) {
  set.seed(seed)
  
  y_pred_baseline <- as.numeric(predict(cv_fit, newx = X_test, s = "lambda.1se"))
  mae_baseline <- mean(abs(y_test - y_pred_baseline))
  
  features_to_test <- intersect(features_to_test, colnames(X_test))
  if (length(features_to_test) == 0) {
    return(tibble(feature = character(), importance = numeric(),
                  importance_sd = numeric(), p_value = numeric()))
  }
  
  feature_seeds <- seed + seq_along(features_to_test) * 1000
  
  importance_results <- parallel::mclapply(
    seq_along(features_to_test),
    function(i) {
      calculate_single_feature_importance(
        features_to_test[i], X_test, y_test,
        y_pred_baseline, mae_baseline,
        cv_fit, n_permutations, feature_seeds[i]
      )
    },
    mc.cores = n_cores
  )
  
  bind_rows(importance_results)
}


#' Get non-zero PENALIZED features (excluding covariates)
get_nonzero_features <- function(results_row, covariate_cols = COVARIATE_COLS) {
  coefs <- results_row$all_coefficients[[1]]
  feat_names <- results_row$all_predictor_names[[1]]
  if (is.null(coefs) || is.null(feat_names)) return(character())
  
  nonzero_mask <- coefs != 0
  is_covariate <- grepl(paste0("^(", paste(covariate_cols, collapse = "|"), ")"), feat_names)
  
  feat_names[nonzero_mask & !is_covariate]
}


#' Build dataset for a given model (handles score-specific, block-only, and score+block)
build_model_dataset <- function(cvd_score, dataset_name,
                                score_specific_datasets, datasets_list) {
  if (grepl("_specific$", dataset_name)) {
    return(score_specific_datasets[[cvd_score]])
  }
  if (dataset_name %in% names(datasets_list)) {
    return(datasets_list[[dataset_name]])
  }
  # Score-specific + block
  block <- sub("^.*_specific\\+", "", dataset_name)
  if (!(block %in% names(datasets_list))) {
    warning(sprintf("Block '%s' not found", block))
    return(NULL)
  }
  build_score_plus_block(cvd_score, block, score_specific_datasets, datasets_list)
}


#' Run permutation importance for a single model with nested CV
run_permutation_importance <- function(df_raw, outcome_col, exclude_cols,
                                       features_to_test,
                                       covariate_cols = COVARIATE_COLS,
                                       alpha_grid = seq(0, 1, by = 0.1),
                                       nfolds = 5,
                                       n_permutations = 1000,
                                       n_cores = 1,
                                       col_miss_thresh = 0.4,
                                       row_miss_thresh = 0.8,
                                       seed = 42) {
  set.seed(seed)
  
  df_complete <- df_raw %>% filter(!is.na(.data[[outcome_col]]))
  y <- df_complete[[outcome_col]]
  all_exclude <- unique(c(exclude_cols, outcome_col, CVD_scores))
  
  # Protect covariates from missingness filter
  cols_to_protect <- unique(c(all_exclude, covariate_cols))
  df_filtered <- suppressMessages(
    filter_by_missingness(df_complete, cols_to_protect, col_miss_thresh, row_miss_thresh)
  )
  
  # Re-attach covariates if dropped
  for (cov in covariate_cols) {
    if (!cov %in% names(df_filtered) && cov %in% names(df_complete)) {
      df_filtered[[cov]] <- df_complete[[cov]][match(df_filtered$Sample_ID, df_complete$Sample_ID)]
    }
  }
  
  y <- df_filtered[[outcome_col]]
  n <- length(y)
  
  if (n < nfolds * 2) {
    warning(sprintf("⚠ Very few observations (%d) for %d-fold CV", n, nfolds))
  }
  
  foldid <- sample(rep(1:nfolds, length.out = n))
  fold_results <- list()
  
  for (k in 1:nfolds) {
    cat(sprintf("    Fold %d/%d: ", k, nfolds))
    
    idx_test <- which(foldid == k)
    idx_train <- setdiff(1:n, idx_test)
    
    df_train <- df_filtered[idx_train, , drop = FALSE]
    df_test <- df_filtered[idx_test, , drop = FALSE]
    y_train <- y[idx_train]
    y_test <- y[idx_test]
    
    preprocessed <- tryCatch(
      preprocess_train_test(df_train, df_test, all_exclude, covariate_cols),
      error = function(e) { cat(sprintf("Preprocessing failed: %s\n", e$message)); NULL }
    )
    
    if (is.null(preprocessed) || ncol(preprocessed$X_train) == 0) {
      cat("No features after preprocessing\n")
      next
    }
    
    X_train <- preprocessed$X_train
    X_test <- preprocessed$X_test
    pen_factor <- ifelse(preprocessed$is_covariate, 0, 1)
    
    # Fit elastic net with covariate adjustment
    best_cv_fit <- NULL
    best_mse <- Inf
    for (a in alpha_grid) {
      cv_fit <- tryCatch(
        cv.glmnet(x = X_train, y = y_train, alpha = a,
                  penalty.factor = pen_factor,
                  nfolds = nfolds,
                  type.measure = "mse", standardize = FALSE, family = "gaussian"),
        error = function(e) NULL
      )
      if (!is.null(cv_fit) && min(cv_fit$cvm) < best_mse) {
        best_mse <- min(cv_fit$cvm)
        best_cv_fit <- cv_fit
      }
    }
    
    if (is.null(best_cv_fit)) { cat("Model fitting failed\n"); next }
    
    # Match features to test (including dummy-coded versions)
    available_features <- colnames(X_test)
    features_this_fold <- features_to_test[features_to_test %in% available_features]
    for (f in features_to_test) {
      if (!(f %in% available_features)) {
        dummy_matches <- grep(paste0("^", f, "_"), available_features, value = TRUE)
        features_this_fold <- c(features_this_fold, dummy_matches)
      }
    }
    features_this_fold <- unique(features_this_fold)
    
    n_features <- length(features_this_fold)
    cat(sprintf("%d features (parallel on %d cores)... ", n_features, n_cores))
    if (n_features == 0) { cat("No matching features\n"); next }
    
    fold_start <- Sys.time()
    fold_importance <- calculate_fold_importance(
      X_train, X_test, y_train, y_test,
      best_cv_fit, features_this_fold,
      n_permutations, n_cores, seed + k * 10000
    )
    fold_elapsed <- difftime(Sys.time(), fold_start, units = "secs")
    
    fold_importance$fold <- k
    fold_results[[k]] <- fold_importance
    
    cat(sprintf("done in %.1fs (top: %s, Δ MAE = %.4f)\n",
                as.numeric(fold_elapsed),
                fold_importance$feature[which.max(fold_importance$importance)],
                max(fold_importance$importance, na.rm = TRUE)))
  }
  
  all_folds <- bind_rows(fold_results)
  if (nrow(all_folds) == 0) { warning("No results from any fold"); return(NULL) }
  
  # Aggregate across folds
  aggregated <- all_folds %>%
    group_by(feature) %>%
    summarise(
      importance = mean(importance, na.rm = TRUE),
      importance_se = sd(importance, na.rm = TRUE) / sqrt(sum(!is.na(importance))),
      importance_sd_across_folds = sd(importance, na.rm = TRUE),
      importance_sd_within_fold = mean(importance_sd, na.rm = TRUE),
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
    per_fold = all_folds, aggregated = aggregated,
    n_folds_cv = nfolds, n_permutations = n_permutations, n_observations = n
  )
}


#' Run permutation importance batch WITH INCREMENTAL SAVING
run_permutation_importance_batch <- function(results_df,
                                             datasets_list,
                                             score_specific_datasets,
                                             n_permutations = 1000,
                                             nfolds = 5,
                                             n_cores = 1,
                                             covariate_cols = COVARIATE_COLS,
                                             seed = 42,
                                             models_to_run = NULL) {
  
  if (is.null(models_to_run)) {
    models_to_run <- results_df %>%
      select(cvd_score, dataset_name) %>% distinct()
  }
  
  n_models <- nrow(models_to_run)
  
  incremental_dir <- save_output("incremental_saves_importance")
  dir.create(incremental_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("==========================================================\n")
  cat("PERMUTATION FEATURE IMPORTANCE (COV-ADJUSTED, PARALLELIZED)\n")
  cat(sprintf("Models: %d | Perms/feature/fold: %d | Folds: %d | Cores: %d\n",
              n_models, n_permutations, nfolds, n_cores))
  cat(sprintf("Covariates (unpenalized, not tested): %s\n", paste(covariate_cols, collapse = ", ")))
  cat(sprintf("Incremental saves to: %s/\n", incremental_dir))
  cat("==========================================================\n")
  
  all_results <- list()
  
  for (i in 1:n_models) {
    cvd_score <- models_to_run$cvd_score[i]
    dataset_name <- models_to_run$dataset_name[i]
    
    cat(sprintf("\n[%d/%d] %s ~ %s\n", i, n_models, cvd_score, dataset_name))
    
    # Check for existing incremental save
    safe_dataset_name <- gsub("[^A-Za-z0-9_]", "_", dataset_name)
    incremental_file <- file.path(incremental_dir,
                                  sprintf("importance_%03d_%s_%s.rds", i, cvd_score, safe_dataset_name))
    
    if (file.exists(incremental_file)) {
      cat("  ↪ Loading from incremental save\n")
      all_results[[i]] <- readRDS(incremental_file)
      next
    }
    
    model_results <- results_df %>%
      filter(cvd_score == !!cvd_score, dataset_name == !!dataset_name)
    
    if (nrow(model_results) == 0) { cat("  ✗ No elastic net results found\n"); next }
    
    features_to_test <- get_nonzero_features(model_results[1, ], covariate_cols)
    n_nonzero <- length(features_to_test)
    cat(sprintf("  Non-zero penalized features: %d\n", n_nonzero))
    
    if (n_nonzero == 0) { cat("  ✗ No penalized features with non-zero coefficients\n"); next }
    
    df_raw <- build_model_dataset(cvd_score, dataset_name,
                                  score_specific_datasets, datasets_list)
    if (is.null(df_raw)) { cat("  ✗ Could not build dataset\n"); next }
    
    model_start <- Sys.time()
    
    tryCatch({
      importance_results <- run_permutation_importance(
        df_raw = df_raw, outcome_col = cvd_score, exclude_cols = "Sample_ID",
        features_to_test = features_to_test, covariate_cols = covariate_cols,
        alpha_grid = seq(0, 1, by = 0.1),
        nfolds = nfolds, n_permutations = n_permutations,
        n_cores = n_cores, seed = seed + i * 100000
      )
      
      model_elapsed <- difftime(Sys.time(), model_start, units = "mins")
      
      if (!is.null(importance_results)) {
        importance_results$aggregated$cvd_score <- cvd_score
        importance_results$aggregated$dataset_name <- dataset_name
        importance_results$aggregated$n_observations <- importance_results$n_observations
        importance_results$aggregated$n_permutations <- importance_results$n_permutations
        importance_results$aggregated$n_folds_cv <- importance_results$n_folds_cv
        
        all_results[[i]] <- importance_results
        saveRDS(importance_results, incremental_file)
        
        n_sig <- sum(importance_results$aggregated$significant_005, na.rm = TRUE)
        top_feat <- importance_results$aggregated$feature[1]
        top_imp <- importance_results$aggregated$importance[1]
        
        cat(sprintf("  ✓ Complete in %.1f min: %d significant (q < 0.05), top = %s (Δ MAE = %.4f)\n",
                    as.numeric(model_elapsed), n_sig, top_feat, top_imp))
        cat(sprintf("  ✓ Saved to %s\n", basename(incremental_file)))
      } else { cat("  ✗ No results returned\n") }
    }, error = function(e) cat(sprintf("  ✗ ERROR: %s\n", e$message)))
  }
  
  aggregated_all <- bind_rows(lapply(all_results, function(x) {
    if (!is.null(x)) x$aggregated else NULL
  }))
  
  cat(sprintf("\n✓ Completed %d models\n", length(all_results[!sapply(all_results, is.null)])))
  
  list(all_results = all_results, aggregated = aggregated_all)
}


# ============================================
# IMPORTANCE PLOTTING FUNCTIONS
# ============================================

#' Extract coefficient directions for labeling importance plots
extract_coefficient_directions <- function(results_df) {
  direction_list <- list()
  for (i in 1:nrow(results_df)) {
    pnames <- results_df$all_predictor_names[[i]]
    coeffs <- results_df$all_coefficients[[i]]
    if (is.null(pnames) || is.null(coeffs)) next
    direction_list[[i]] <- tibble(
      cvd_score = results_df$cvd_score[i],
      dataset_name = results_df$dataset_name[i],
      feature = pnames, coefficient = coeffs,
      direction = case_when(
        coeffs > 0 ~ "Risk-Increasing", coeffs < 0 ~ "Risk-Decreasing", TRUE ~ "Zero"
      )
    )
  }
  bind_rows(direction_list)
}


#' Bar plot for a single model's important features
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
        p_value_bh < 0.001 ~ "***", p_value_bh < 0.01 ~ "**",
        p_value_bh < 0.05 ~ "*", TRUE ~ ""
      ),
      label_y = ifelse(importance > max(importance) * 0.08,
                       importance * 0.5, importance + max(importance) * 0.02)
    )
  
  if (nrow(plot_data) == 0) return(NULL)
  
  ggplot(plot_data, aes(x = feature, y = importance, fill = direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sig_label, y = label_y),
              hjust = 0.5, size = 3.5, color = "white", fontface = "bold") +
    scale_fill_manual(
      values = c("Risk-Increasing" = "#E07B39", "Risk-Decreasing" = "#3B7EA1", "Unknown" = "grey50"),
      name = "Effect Direction",
      labels = c("Risk-Increasing" = "β > 0", "Risk-Decreasing" = "β < 0", "Unknown" = "Unknown")
    ) +
    coord_flip() +
    labs(
      title = paste0("Permutation Feature Importance: ", cvd_score),
      subtitle = paste0("Dataset: ", dataset_name,
                        " | Adjusted for ", paste(COVARIATE_COLS, collapse = ", "),
                        " | *q<0.05 **q<0.01 ***q<0.001"),
      x = "Feature", y = "Importance (Δ MAE when permuted)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.text.y = element_text(size = 7),
      legend.position = "top",
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()
    )
}


#' Create importance plots for all models
create_importance_plots <- function(importance_results, coefficient_directions,
                                    output_dir = "importance_plots") {
  full_output_dir <- save_output(output_dir)
  dir.create(full_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  importance_df <- importance_results$aggregated
  models <- importance_df %>% select(cvd_score, dataset_name) %>% distinct()
  cat(sprintf("\nCreating %d importance plots...\n", nrow(models)))
  
  plot_list <- list()
  for (i in 1:nrow(models)) {
    cvd <- models$cvd_score[i]; dataset <- models$dataset_name[i]
    model_importance <- importance_df %>% filter(cvd_score == cvd, dataset_name == dataset)
    if (nrow(model_importance) == 0) next
    
    n_sig <- sum(model_importance$significant_005, na.rm = TRUE)
    cat(sprintf("  [%d/%d] %s ~ %s: %d features, %d significant\n",
                i, nrow(models), cvd, dataset, nrow(model_importance), n_sig))
    
    p <- plot_importance_bars(model_importance, coefficient_directions, cvd, dataset)
    if (is.null(p)) next
    plot_list[[i]] <- p
    
    filename <- file.path(full_output_dir, paste0("importance_", cvd, "_",
                                                  gsub("[^A-Za-z0-9_]", "_", dataset), ".pdf"))
    ggsave(filename, p, width = 10, height = max(6, min(25, 3 + 0.3 * nrow(model_importance))),
           device = "pdf")
  }
  
  cat(sprintf("✓ Saved %d plots\n", length(plot_list[!sapply(plot_list, is.null)])))
  plot_list
}


#' Faceted block importance plot — all features with positive importance
plot_block_importance <- function(importance_data, coefficient_directions,
                                  block_name) {
  plot_data <- importance_data %>%
    left_join(
      coefficient_directions %>%
        select(cvd_score, dataset_name, feature, coefficient, direction),
      by = c("cvd_score", "dataset_name", "feature")
    ) %>%
    filter(!is.na(importance), importance > 0)
  
  if (nrow(plot_data) == 0) return(NULL)
  
  n_features <- n_distinct(plot_data$feature)
  x_text_size <- case_when(block_name == "all_data" ~ 5, n_features > 30 ~ 6, TRUE ~ 7)
  
  plot_data <- plot_data %>%
    mutate(
      cvd_score_label = recode(cvd_score, "ascvd_10y" = "ASCVD", "frs_10y" = "Framingham",
                               "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2"),
      direction = ifelse(is.na(direction), "Unknown", direction),
      sig_label = case_when(p_value_bh < 0.001 ~ "***", p_value_bh < 0.01 ~ "**",
                            p_value_bh < 0.05 ~ "*", TRUE ~ ""),
      feature_short = smart_truncate(feature, 30),
      feature_id = paste(feature_short, cvd_score_label, sep = "___")
    ) %>%
    arrange(cvd_score_label, desc(importance))
  
  plot_data$feature_id <- factor(plot_data$feature_id, levels = unique(plot_data$feature_id))
  plot_data$cvd_score_label <- factor(plot_data$cvd_score_label,
                                      levels = c("ASCVD", "Framingham", "QRISK3", "SCORE2"))
  
  plot_data <- plot_data %>%
    group_by(cvd_score_label) %>%
    mutate(max_imp = max(importance),
           label_y = ifelse(importance > max_imp * 0.08,
                            importance * 0.5, importance + max_imp * 0.02)) %>%
    ungroup()
  
  ggplot(plot_data, aes(x = feature_id, y = importance, fill = direction)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sig_label, y = label_y), size = 2.5, color = "white", fontface = "bold") +
    facet_wrap(~ cvd_score_label, scales = "free", nrow = 1) +
    scale_x_discrete(labels = function(x) gsub("___.*$", "", x)) +
    scale_fill_manual(
      values = c("Risk-Increasing" = "#E07B39", "Risk-Decreasing" = "#3B7EA1", "Unknown" = "grey50"),
      name = "Effect Direction"
    ) +
    labs(
      title = sprintf("Permutation Importance: %s (Adjusted for %s)",
                      block_name, paste(COVARIATE_COLS, collapse = ", ")),
      subtitle = sprintf("All features with positive importance | *q<0.05 **q<0.01 ***q<0.001"),
      x = NULL, y = "Importance (Δ MAE)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "gray40", hjust = 0.5),
      axis.text.x = element_text(size = x_text_size, angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom", panel.grid.major.x = element_blank()
    )
}


#' Create block-level importance plots
create_block_importance_plots <- function(importance_results, coefficient_directions,
                                          blocks = c("all_data", "lipids", "fatty_acids", "urine_nmr",
                                                     "body_composition", "clinical_risk_factors",
                                                     "sociodemographics_lifestyle"),
                                          output_dir = "importance_plots_by_block") {
  full_output_dir <- save_output(output_dir)
  dir.create(full_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  importance_df <- importance_results$aggregated %>%
    mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
    filter(block %in% blocks)
  
  cat(sprintf("\nCreating block importance plots for %d blocks...\n", length(blocks)))
  
  for (current_block in blocks) {
    block_data <- filter(importance_df, block == current_block)
    if (nrow(block_data) == 0) { cat(sprintf("  %s: No data\n", current_block)); next }
    
    n_features <- n_distinct(block_data$feature)
    cat(sprintf("  %s: %d features across %d scores\n",
                current_block, n_features, n_distinct(block_data$cvd_score)))
    
    p <- plot_block_importance(block_data, coefficient_directions, current_block)
    if (is.null(p)) next
    
    n_shown <- block_data %>% filter(importance > 0) %>%
      pull(feature) %>% n_distinct()
    
    ggsave(file.path(full_output_dir, sprintf("importance_block_%s.pdf", current_block)),
           p, width = max(12, min(30, 6 + 0.18 * n_shown)), height = 7)
  }
  cat(sprintf("✓ Saved block importance plots\n"))
}

# ============================================
# 15. EXECUTE PERMUTATION IMPORTANCE --> does only score-specifc + block? move pn to 16. for block only
# ============================================

cat("\n==========================================================\n")
cat("PERMUTATION IMPORTANCE: SCORE+BLOCK MODELS (PARALLELIZED)\n")
cat(sprintf("Permutations: 1000 | Cores: %d\n", CPUS))
cat("==========================================================\n")

models_for_full_run <- results_score_plus_blocks %>%
  select(cvd_score, dataset_name) %>% distinct()
cat(sprintf("Running %d models\n", nrow(models_for_full_run)))

start_time <- Sys.time()

importance_results_full <- run_permutation_importance_batch(
  results_df = results_score_plus_blocks,
  datasets_list = datasets_list,
  score_specific_datasets = score_specific_datasets,
  n_permutations = 1000, nfolds = 5, n_cores = CPUS,
  covariate_cols = COVARIATE_COLS, seed = 42,
  models_to_run = models_for_full_run
)

end_time <- Sys.time()
total_runtime_spb <- difftime(end_time, start_time, units = "hours")
cat(sprintf("\nScore+block runtime: %.2f hours\n", as.numeric(total_runtime_spb)))

# Save
qs_save(importance_results_full, save_output("permutation_importance_full_results.qs2"))
saveRDS(importance_results_full, save_output("permutation_importance_full_results.rds"))

importance_full_excel <- importance_results_full$aggregated %>%
  select(cvd_score, dataset_name, feature,
         importance, importance_se, importance_sd_across_folds, importance_sd_within_fold,
         p_value_fisher, p_value_mean, p_value_bh,
         significant_005, significant_010,
         n_observations, n_permutations, n_folds_cv) %>%
  arrange(cvd_score, dataset_name, desc(importance))
write.xlsx(importance_full_excel, save_output("permutation_importance_full_results.xlsx"))

# Plots for score+block
coefficient_directions_full <- extract_coefficient_directions(results_score_plus_blocks)

create_importance_plots(importance_results_full, coefficient_directions_full,
                        output_dir = "importance_plots_full")

create_block_importance_plots(
  importance_results = importance_results_full,
  coefficient_directions = coefficient_directions_full,
  output_dir = "importance_plots_by_block_full"
)

# Summary
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

importance_overview <- importance_full_summary %>%
  select(cvd_score, dataset_name, n_features_tested, n_significant_005,
         n_risk_increasing, n_risk_decreasing, top_feature, top_importance) %>%
  mutate(
    block = sub("^.*_specific\\+", "", dataset_name),
    cvd_label = cvd_score_labels[cvd_score]
  )
write.xlsx(importance_overview, save_output("permutation_importance_overview.xlsx"))

# ============================================
# 16. PERMUTATION IMPORTANCE: COMMON (BLOCK-ONLY) MODELS
# ============================================

cat("\n==========================================================\n")
cat("PERMUTATION IMPORTANCE: COMMON (BLOCK-ONLY) MODELS\n")
cat(sprintf("Permutations: 1000 | Cores: %d\n", CPUS))
cat("==========================================================\n")

models_common_run <- results_common %>%
  select(cvd_score, dataset_name) %>% distinct()
cat(sprintf("Running %d models\n", nrow(models_common_run)))

start_time_common <- Sys.time()

importance_results_common <- run_permutation_importance_batch(
  results_df = results_common,
  datasets_list = datasets_list,
  score_specific_datasets = score_specific_datasets,
  n_permutations = 1000, nfolds = 5, n_cores = CPUS,
  covariate_cols = COVARIATE_COLS, seed = 42,
  models_to_run = models_common_run
)

end_time_common <- Sys.time()
total_runtime_common <- difftime(end_time_common, start_time_common, units = "hours")
cat(sprintf("\nCommon runtime: %.2f hours\n", as.numeric(total_runtime_common)))

# Save
qs_save(importance_results_common, save_output("permutation_importance_common_results.qs2"))
saveRDS(importance_results_common, save_output("permutation_importance_common_results.rds"))

importance_common_excel <- importance_results_common$aggregated %>%
  select(cvd_score, dataset_name, feature,
         importance, importance_se, importance_sd_across_folds, importance_sd_within_fold,
         p_value_fisher, p_value_mean, p_value_bh,
         significant_005, significant_010,
         n_observations, n_permutations, n_folds_cv) %>%
  arrange(cvd_score, dataset_name, desc(importance))
write.xlsx(importance_common_excel, save_output("permutation_importance_common_results.xlsx"))

# Plots for common
coefficient_directions_common <- extract_coefficient_directions(results_common)

create_importance_plots(importance_results_common, coefficient_directions_common,
                        output_dir = "importance_plots_common")

# Block-level plots for common (dataset_name IS the block)
create_block_importance_plots_common <- function(importance_results, coefficient_directions,
                                                 blocks = c("all_data", "lipids", "fatty_acids", "urine_nmr",
                                                            "body_composition", "clinical_risk_factors",
                                                            "sociodemographics_lifestyle"),
                                                 output_dir = "importance_plots_common_by_block") {
  full_output_dir <- save_output(output_dir)
  dir.create(full_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  importance_df <- importance_results$aggregated %>%
    filter(dataset_name %in% blocks)
  
  cat(sprintf("\nCreating block importance plots for %d blocks (block-only)...\n", length(blocks)))
  
  for (current_block in blocks) {
    block_data <- filter(importance_df, dataset_name == current_block)
    if (nrow(block_data) == 0) { cat(sprintf("  %s: No data\n", current_block)); next }
    
    n_features <- n_distinct(block_data$feature)
    cat(sprintf("  %s: %d features across %d scores\n",
                current_block, n_features, n_distinct(block_data$cvd_score)))
    
    p <- plot_block_importance(block_data, coefficient_directions, current_block)
    if (is.null(p)) next
    
    n_shown <- block_data %>% filter(importance > 0) %>%
      pull(feature) %>% n_distinct()
    
    ggsave(file.path(full_output_dir, sprintf("importance_block_only_%s.pdf", current_block)),
           p, width = max(12, min(30, 6 + 0.18 * n_shown)), height = 7)
  }
  cat("✓ Saved common block importance plots\n")
}

create_block_importance_plots_common(
  importance_results = importance_results_common,
  coefficient_directions = coefficient_directions_common,
  output_dir = "importance_plots_common_by_block"
)

cat("\n==========================================================\n")
cat("✓ ALL PERMUTATION IMPORTANCE ANALYSES COMPLETE\n")
cat(sprintf("  Score+block runtime: %.2f hours\n", as.numeric(total_runtime_spb)))
cat(sprintf("  Common runtime: %.2f hours\n", as.numeric(total_runtime_common)))
cat(sprintf("  Results saved to: %s\n", OUTPUT_DIR))
cat("==========================================================\n")

# ============================================
# 17. IMPORTANCE TABLE FUNCTIONS
# ============================================

#' Create a combined importance table for one dataset across all CVD scores
create_dataset_importance_table <- function(importance_df,
                                            coefficient_df,
                                            dataset_pattern,
                                            min_importance = 0,
                                            result_type = "score_plus_block") {
  if (result_type == "score_plus_block") {
    data_filtered <- importance_df %>%
      filter(grepl(paste0("_specific\\+", dataset_pattern, "$"), dataset_name))
  } else {
    data_filtered <- importance_df %>%
      filter(dataset_name == dataset_pattern)
  }
  
  if (nrow(data_filtered) == 0) {
    warning(sprintf("No data found for pattern: %s", dataset_pattern))
    return(NULL)
  }
  
  data_joined <- data_filtered %>%
    left_join(
      coefficient_df %>%
        select(cvd_score, dataset_name, feature, coefficient, direction),
      by = c("cvd_score", "dataset_name", "feature")
    )
  
  data_positive <- data_joined %>% filter(importance > min_importance)
  if (nrow(data_positive) == 0) {
    warning(sprintf("No features with importance > %s for: %s", min_importance, dataset_pattern))
    return(NULL)
  }
  
  data_positive <- data_positive %>%
    mutate(
      sig_stars = case_when(
        p_value_bh < 0.001 ~ "***", p_value_bh < 0.01 ~ "**",
        p_value_bh < 0.05 ~ "*", TRUE ~ ""
      ),
      cvd_label = recode(cvd_score,
                         "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                         "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD")
    )
  
  feature_order <- data_positive %>%
    group_by(feature) %>%
    summarise(avg_importance = mean(importance, na.rm = TRUE),
              n_scores = n(), .groups = "drop") %>%
    arrange(desc(avg_importance))
  
  score_order <- c("QRISK3", "SCORE2", "Framingham", "ASCVD")
  
  values_wide <- data_positive %>%
    select(feature, cvd_label, sig_stars) %>%
    pivot_wider(names_from = cvd_label, values_from = sig_stars, values_fill = NA_character_)
  
  directions_wide <- data_positive %>%
    select(feature, cvd_label, direction) %>%
    pivot_wider(names_from = cvd_label, values_from = direction, values_fill = NA_character_)
  
  table_final <- feature_order %>%
    select(feature, avg_importance) %>%
    left_join(values_wide, by = "feature") %>%
    arrange(desc(avg_importance))
  
  directions_final <- feature_order %>%
    select(feature) %>%
    left_join(directions_wide, by = "feature")
  
  score_cols_present <- intersect(score_order, names(table_final))
  table_final <- table_final %>% select(feature, all_of(score_cols_present))
  directions_final <- directions_final %>% select(feature, any_of(score_cols_present))
  
  list(
    table = table_final, directions = directions_final,
    n_features = nrow(table_final), dataset = dataset_pattern,
    score_cols = score_cols_present
  )
}


#' Render importance table as a formatted flextable
render_importance_flextable <- function(table_data, title = "Permutation Feature Importance") {
  tbl <- table_data$table
  directions <- table_data$directions
  score_cols <- table_data$score_cols
  
  ft <- flextable(tbl) %>%
    set_caption(caption = title) %>%
    set_header_labels(feature = "Feature") %>%
    align(j = 1, align = "left", part = "all") %>%
    align(j = score_cols, align = "center", part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 9, part = "body") %>%
    flextable::fontsize(size = 10, part = "header") %>%
    bold(part = "header") %>%
    width(j = "feature", width = 2.5) %>%
    width(j = score_cols, width = 1.0) %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1.5), part = "header") %>%
    padding(padding = 3, part = "all")
  
  # Color cells by direction
  for (col in score_cols) {
    if (col %in% names(directions)) {
      for (row_idx in 1:nrow(directions)) {
        dir_val <- directions[[col]][row_idx]
        if (!is.na(dir_val)) {
          bg_color <- case_when(
            dir_val == "Risk-Increasing" ~ "#FDEBD0",
            dir_val == "Risk-Decreasing" ~ "#D6EAF8",
            TRUE ~ "white"
          )
          ft <- bg(ft, i = row_idx, j = col, bg = bg_color)
        }
      }
    }
  }
  
  # Add footnote
  ft <- add_footer_lines(ft, values = c(
    "* q < 0.05, ** q < 0.01, *** q < 0.001 (BH-corrected, Fisher's method across folds)",
    "Orange = Risk-Increasing (β > 0), Blue = Risk-Decreasing (β < 0)",
    "Features sorted by average importance across scores"
  ))
  ft <- flextable::fontsize(ft, size = 7, part = "footer")
  ft <- italic(ft, part = "footer")
  
  ft
}


#' Create and save all dataset tables
create_all_importance_tables <- function(importance_results,
                                         coefficient_directions,
                                         datasets = c("all_data", "lipids", "fatty_acids",
                                                      "urine_nmr", "body_composition",
                                                      "clinical_risk_factors",
                                                      "sociodemographics_lifestyle"),
                                         output_dir = "importance_tables",
                                         result_type = "score_plus_block",
                                         dataset_labels = NULL) {
  full_output_dir <- save_output(output_dir)
  dir.create(full_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  importance_df <- importance_results$aggregated
  
  if (is.null(dataset_labels)) {
    dataset_labels <- c(
      "all_data" = "Full predictor set",
      "lipids" = "Lipids",
      "fatty_acids" = "Fatty Acids",
      "urine_nmr" = "Urine NMR Metabolites",
      "body_composition" = "Body Composition",
      "clinical_risk_factors" = "Clinical Measurements & Supplementary Biomarkers",
      "sociodemographics_lifestyle" = "Sociodemographics & Lifestyle"
    )
  }
  
  type_label <- if (result_type == "score_plus_block") "Score-Specific + " else ""
  
  cat(sprintf("\nCreating importance tables for %d datasets (%s)...\n",
              length(datasets),
              if (result_type == "score_plus_block") "score+block" else "blocks only"))
  
  all_tables <- list()
  
  for (ds in datasets) {
    cat(sprintf("  %s: ", ds))
    
    table_data <- create_dataset_importance_table(
      importance_df = importance_df,
      coefficient_df = coefficient_directions,
      dataset_pattern = ds, min_importance = 0,
      result_type = result_type
    )
    
    if (is.null(table_data)) { cat("No data\n"); next }
    cat(sprintf("%d features\n", table_data$n_features))
    
    title_text <- paste0("Permutation Feature Importance: ", type_label,
                         ifelse(ds %in% names(dataset_labels), dataset_labels[ds], ds))
    
    ft <- render_importance_flextable(table_data, title = title_text)
    
    save_as_docx(ft, path = file.path(full_output_dir, paste0("importance_table_", ds, ".docx")))
    save_as_html(ft, path = file.path(full_output_dir, paste0("importance_table_", ds, ".html")))
    
    all_tables[[ds]] <- list(data = table_data, flextable = ft)
  }
  
  # Combined document
  cat("\nCreating combined document...\n")
  combined_doc <- read_docx()
  combined_doc <- body_add_par(combined_doc, "Permutation Feature Importance Tables", style = "heading 1")
  combined_doc <- body_add_par(combined_doc,
                               paste0(type_label, "Predictor Blocks"),
                               style = "heading 2")
  combined_doc <- body_add_par(combined_doc,
                               sprintf("Adjusted for: %s", paste(COVARIATE_COLS, collapse = ", ")))
  combined_doc <- body_add_par(combined_doc, "")
  
  for (ds in names(all_tables)) {
    title_text <- ifelse(ds %in% names(dataset_labels), dataset_labels[ds], ds)
    combined_doc <- body_add_par(combined_doc, title_text, style = "heading 3")
    combined_doc <- body_add_flextable(combined_doc, value = all_tables[[ds]]$flextable)
    combined_doc <- body_add_break(combined_doc)
  }
  
  print(combined_doc, target = file.path(full_output_dir, "All_Importance_Tables_Combined.docx"))
  
  cat(sprintf("✓ Tables saved to %s/\n", output_dir))
  all_tables
}


# ============================================
# 18. CREATE TABLES FOR BOTH TYPES
# ============================================

cat("\n=== Creating tables for SCORE+BLOCKS results ===\n")
importance_tables_score_plus_blocks <- create_all_importance_tables(
  importance_results = importance_results_full,
  coefficient_directions = coefficient_directions_full,
  output_dir = "importance_tables_score_plus_blocks",
  result_type = "score_plus_block"
)

cat("\n=== Creating tables for COMMON (blocks only) results ===\n")
importance_tables_common <- create_all_importance_tables(
  importance_results = importance_results_common,
  coefficient_directions = coefficient_directions_common,
  output_dir = "importance_tables_common",
  result_type = "common"
)

cat("\n✓ All importance tables created\n")
# ============================================
# 19. SIGNIFICANT-ONLY IMPORTANCE TABLES
# ============================================

#' Create importance table with ONLY features significant in at least one score
create_significant_only_table <- function(importance_df,
                                          coefficient_df,
                                          dataset_pattern,
                                          min_importance = 0,
                                          result_type = "score_plus_block") {
  if (result_type == "score_plus_block") {
    data_filtered <- importance_df %>%
      filter(grepl(paste0("_specific\\+", dataset_pattern, "$"), dataset_name))
  } else {
    data_filtered <- importance_df %>%
      filter(dataset_name == dataset_pattern)
  }
  
  if (nrow(data_filtered) == 0) {
    warning(sprintf("No data found for pattern: %s", dataset_pattern))
    return(NULL)
  }
  
  data_joined <- data_filtered %>%
    left_join(
      coefficient_df %>%
        select(cvd_score, dataset_name, feature, coefficient, direction),
      by = c("cvd_score", "dataset_name", "feature")
    )
  
  data_positive <- data_joined %>% filter(importance > min_importance)
  if (nrow(data_positive) == 0) {
    warning(sprintf("No features with importance > %s for: %s", min_importance, dataset_pattern))
    return(NULL)
  }
  
  # Features significant in at least one score
  significant_features <- data_positive %>%
    filter(p_value_bh < 0.05) %>% pull(feature) %>% unique()
  
  if (length(significant_features) == 0) {
    warning(sprintf("No significant features for: %s", dataset_pattern))
    return(NULL)
  }
  
  data_significant <- data_positive %>%
    filter(feature %in% significant_features) %>%
    mutate(
      sig_stars = case_when(
        p_value_bh < 0.001 ~ "***", p_value_bh < 0.01 ~ "**",
        p_value_bh < 0.05 ~ "*", TRUE ~ ""
      ),
      cvd_label = recode(cvd_score,
                         "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                         "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD")
    )
  
  feature_order <- data_significant %>%
    group_by(feature) %>%
    summarise(avg_importance = mean(importance, na.rm = TRUE),
              n_significant = sum(p_value_bh < 0.05, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(n_significant), desc(avg_importance))
  
  score_order <- c("QRISK3", "SCORE2", "Framingham", "ASCVD")
  
  values_wide <- data_significant %>%
    select(feature, cvd_label, sig_stars) %>%
    pivot_wider(names_from = cvd_label, values_from = sig_stars, values_fill = NA_character_)
  
  directions_wide <- data_significant %>%
    select(feature, cvd_label, direction) %>%
    pivot_wider(names_from = cvd_label, values_from = direction, values_fill = NA_character_)
  
  table_final <- feature_order %>%
    select(feature, avg_importance, n_significant) %>%
    left_join(values_wide, by = "feature") %>%
    arrange(desc(n_significant), desc(avg_importance))
  
  directions_final <- feature_order %>%
    select(feature) %>%
    left_join(directions_wide, by = "feature")
  
  score_cols_present <- intersect(score_order, names(table_final))
  table_final <- table_final %>% select(feature, all_of(score_cols_present))
  directions_final <- directions_final %>% select(feature, any_of(score_cols_present))
  
  list(
    table = table_final, directions = directions_final,
    n_features = nrow(table_final),
    n_total_positive = n_distinct(data_positive$feature),
    dataset = dataset_pattern, score_cols = score_cols_present
  )
}


#' Render significant-only table as flextable
render_significant_flextable <- function(table_data, title = NULL) {
  if (is.null(table_data)) return(NULL)
  
  tbl <- table_data$table
  dirs <- table_data$directions
  score_cols <- table_data$score_cols
  
  for (col in score_cols) tbl[[col]] <- ifelse(is.na(tbl[[col]]), "", tbl[[col]])
  
  ft <- flextable(tbl) %>%
    set_header_labels(feature = "Feature") %>%
    align(j = 1, align = "left", part = "all") %>%
    align(j = score_cols, align = "center", part = "all") %>%
    bold(part = "header") %>%
    flextable::fontsize(size = 9, part = "body") %>%
    flextable::fontsize(size = 10, part = "header") %>%
    font(fontname = "Arial", part = "all") %>%
    width(j = "feature", width = 2.5) %>%
    width(j = score_cols, width = 1.0) %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1.5), part = "header") %>%
    padding(padding = 4, part = "all")
  
  for (col in score_cols) {
    for (i in 1:nrow(tbl)) {
      direction <- dirs[[col]][i]
      if (!is.na(direction) && direction %in% c("Risk-Increasing", "Risk-Decreasing")) {
        bg_color <- if (direction == "Risk-Increasing") "#FDDBC7" else "#D1E5F0"
        text_color <- if (direction == "Risk-Increasing") "#B35806" else "#2166AC"
        ft <- ft %>%
          bg(i = i, j = col, bg = bg_color, part = "body") %>%
          color(i = i, j = col, color = text_color, part = "body") %>%
          bold(i = i, j = col, part = "body")
      }
    }
  }
  
  if (!is.null(title)) ft <- set_caption(ft, caption = title)
  
  footnote_text <- sprintf(
    "Showing %d features significant (q < 0.05) in ≥1 score (of %d with positive importance) | * q<0.05, ** q<0.01, *** q<0.001 | Orange = β > 0, Blue = β < 0 | Empty = not selected or importance ≤ 0",
    table_data$n_features, table_data$n_total_positive
  )
  ft <- ft %>%
    add_footer_lines(values = footnote_text) %>%
    flextable::fontsize(size = 8, part = "footer") %>%
    italic(part = "footer")
  ft
}


#' Create and save significant-only tables for all datasets
create_significant_tables <- function(importance_results, coefficient_directions,
                                      datasets = c("all_data", "lipids", "fatty_acids",
                                                   "urine_nmr", "body_composition",
                                                   "clinical_risk_factors",
                                                   "sociodemographics_lifestyle"),
                                      output_dir = "importance_tables_significant",
                                      result_type = "score_plus_block",
                                      dataset_labels = NULL) {
  full_output_dir <- save_output(output_dir)
  dir.create(full_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  importance_df <- importance_results$aggregated
  if (is.null(dataset_labels)) {
    dataset_labels <- c(
      "all_data" = "Full predictor set", "lipids" = "Lipids",
      "fatty_acids" = "Fatty Acids", "urine_nmr" = "Urine NMR Metabolites",
      "body_composition" = "Body Composition",
      "clinical_risk_factors" = "Clinical Measurements & Supplementary Biomarkers",
      "sociodemographics_lifestyle" = "Sociodemographics & Lifestyle"
    )
  }
  
  type_label <- if (result_type == "score_plus_block") "Score-Specific + " else ""
  cat(sprintf("\nCreating SIGNIFICANT-ONLY tables for %d datasets (%s)...\n",
              length(datasets),
              if (result_type == "score_plus_block") "score+block" else "blocks only"))
  
  all_tables <- list()
  for (ds in datasets) {
    cat(sprintf("  %s: ", ds))
    table_data <- create_significant_only_table(
      importance_df, coefficient_directions, ds, 0, result_type
    )
    if (is.null(table_data)) { cat("No significant features\n"); next }
    cat(sprintf("%d significant features (of %d positive)\n",
                table_data$n_features, table_data$n_total_positive))
    
    title_text <- paste0("Significant Features: ", type_label,
                         ifelse(ds %in% names(dataset_labels), dataset_labels[ds], ds))
    ft <- render_significant_flextable(table_data, title = title_text)
    
    save_as_docx(ft, path = file.path(full_output_dir, paste0("significant_table_", ds, ".docx")))
    save_as_html(ft, path = file.path(full_output_dir, paste0("significant_table_", ds, ".html")))
    all_tables[[ds]] <- list(data = table_data, flextable = ft)
  }
  
  if (length(all_tables) > 0) {
    cat("\nCreating combined document...\n")
    combined_doc <- read_docx()
    combined_doc <- body_add_par(combined_doc, "Significant Features Only", style = "heading 1")
    combined_doc <- body_add_par(combined_doc,
                                 paste0(type_label, "Predictor Blocks"),
                                 style = "heading 2")
    combined_doc <- body_add_par(combined_doc,
                                 sprintf("Adjusted for: %s | Features with q < 0.05 in ≥1 CVD score",
                                         paste(COVARIATE_COLS, collapse = ", ")))
    combined_doc <- body_add_par(combined_doc, "")
    for (ds in names(all_tables)) {
      title_text <- ifelse(ds %in% names(dataset_labels), dataset_labels[ds], ds)
      combined_doc <- body_add_par(combined_doc, title_text, style = "heading 3")
      combined_doc <- body_add_flextable(combined_doc, value = all_tables[[ds]]$flextable)
      combined_doc <- body_add_break(combined_doc)
    }
    print(combined_doc, target = file.path(full_output_dir, "Significant_Features_Combined.docx"))
  }
  
  cat(sprintf("✓ Tables saved to %s/\n", output_dir))
  all_tables
}


# ---- Run significant tables ----
cat("\n=== SIGNIFICANT-ONLY tables for SCORE+BLOCKS ===\n")
significant_tables_score_plus_blocks <- create_significant_tables(
  importance_results_full, coefficient_directions_full,
  output_dir = "significant_tables_score_plus_blocks",
  result_type = "score_plus_block"
)

cat("\n=== SIGNIFICANT-ONLY tables for COMMON (blocks only) ===\n")
significant_tables_common <- create_significant_tables(
  importance_results_common, coefficient_directions_common,
  output_dir = "significant_tables_common",
  result_type = "common"
)



# ============================================
# 20. UNIFIED BAR IMPORTANCE PLOTS
# ============================================
# Four filter modes:
#   a) all_selected:           non-zero elastic net coefficient
#   b) significant:            + significant (q < 0.05)
#   c) positive_delta:         + positive mean ΔMAE
#   d) positive_delta_3folds:  + positive ΔMAE in ≥3 folds
#
# Each mode × each result type (common, score+block) = separate output folder.
# ============================================

library(ggplot2)

#' Extract per-fold importance data from batch results
extract_per_fold_data <- function(importance_results) {
  all_fold_data <- list()
  for (i in seq_along(importance_results$all_results)) {
    res <- importance_results$all_results[[i]]
    if (is.null(res)) next
    fold_df <- res$per_fold
    if (is.null(fold_df) || nrow(fold_df) == 0) next
    fold_df$cvd_score <- res$aggregated$cvd_score[1]
    fold_df$dataset_name <- res$aggregated$dataset_name[1]
    all_fold_data[[i]] <- fold_df
  }
  bind_rows(all_fold_data)
}


#' Reference level legend for sociodemographics dummy variables
get_sociodemographic_references <- function() {
  c(
    "Living Status" = "ref: Living alone",
    "Employment Status" = "ref: Not working/Other",
    "Education Level" = "ref: Level 0 (lowest)",
    "Marital status" = "ref: Married/In partnership",
    "Working time" = "ref: Full time",
    "Annual net salary" = "ref: Category 0 (lowest)",
    "Naps during day" = "ref: No",
    "Servings nuts per week" = "ref: < 3",
    "Wine per week" = "ref: < 7 glasses",
    "Olive oil given day" = "ref: < 4 tbsp",
    "Olive oil as main culinary fat" = "ref: No"
  )
}


#' Create direction-aware bar plot for one dataset block
#'
#' @param filter_mode One of: "all_selected", "significant", "positive_delta", "positive_delta_3folds"
create_importance_barplot <- function(fold_data, coefficient_df, aggregated_df,
                                      block_name, block_label = NULL,
                                      cvd_colors = NULL,
                                      result_type = "common",
                                      filter_mode = "positive_delta_3folds") {
  
  if (is.null(cvd_colors)) {
    cvd_colors <- c(
      "QRISK3" = "#F8766D", "SCORE2" = "#00BA38",
      "Framingham" = "#B79F00", "ASCVD" = "#619CFF"
    )
  }
  cvd_order <- c("QRISK3", "SCORE2", "Framingham", "ASCVD")
  if (is.null(block_label)) block_label <- block_name
  
  min_nonzero_folds <- if (filter_mode == "positive_delta_3folds") 3 else 1
  
  # ---- Filter fold data to this block ----
  if (result_type == "score_plus_block") {
    plot_data <- fold_data %>%
      filter(grepl(paste0("_specific\\+", block_name, "$"), dataset_name))
  } else {
    plot_data <- fold_data %>%
      filter(dataset_name == block_name)
  }
  if (nrow(plot_data) == 0) return(NULL)
  
  # ---- Join with coefficient directions (step 1: non-zero coef) ----
  plot_data <- plot_data %>%
    left_join(
      coefficient_df %>% select(cvd_score, dataset_name, feature, direction),
      by = c("cvd_score", "dataset_name", "feature")
    ) %>%
    filter(!is.na(direction), direction != "Zero")
  if (nrow(plot_data) == 0) return(NULL)
  
  # ---- Labels ----
  plot_data <- plot_data %>%
    mutate(
      cvd_label = recode(cvd_score,
                         "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                         "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD"),
      cvd_label = factor(cvd_label, levels = cvd_order)
    )
  
  # ---- Apply filter_mode ----
  n_excluded_combos <- 0
  filter_description <- ""
  
  if (filter_mode == "all_selected") {
    # Keep everything (all folds, including importance == 0)
    filter_description <- "All features with non-zero elastic net coefficients"
    
  } else if (filter_mode == "significant") {
    # Keep only features significant in aggregated results
    if (result_type == "score_plus_block") {
      agg_sig <- aggregated_df %>%
        filter(grepl(paste0("_specific\\+", block_name, "$"), dataset_name), p_value_bh < 0.05)
    } else {
      agg_sig <- aggregated_df %>%
        filter(dataset_name == block_name, p_value_bh < 0.05)
    }
    sig_combos <- agg_sig %>%
      mutate(cvd_label = recode(cvd_score,
                                "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                                "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD")) %>%
      select(feature, cvd_label)
    
    plot_data <- plot_data %>% semi_join(sig_combos, by = c("feature", "cvd_label"))
    filter_description <- "Significant features only (q < 0.05)"
    
  } else if (filter_mode == "positive_delta") {
    # Keep only folds with importance > 0
    plot_data <- plot_data %>% filter(importance > 0)
    filter_description <- "Features with positive \u0394MAE"
    
  } else if (filter_mode == "positive_delta_3folds") {
    # Keep only folds with importance > 0, then require ≥3 such folds
    plot_data <- plot_data %>% filter(importance > 0)
    
    fold_counts <- plot_data %>%
      group_by(feature, cvd_label) %>%
      summarise(n_nonzero = n(), .groups = "drop")
    valid_combos <- fold_counts %>% filter(n_nonzero >= min_nonzero_folds)
    n_excluded_combos <- nrow(fold_counts) - nrow(valid_combos)
    
    plot_data <- plot_data %>%
      semi_join(valid_combos, by = c("feature", "cvd_label"))
    filter_description <- "Features with positive \u0394MAE in \u22653 folds"
  }
  
  if (nrow(plot_data) == 0) return(NULL)
  
  # ---- Summarise: mean ± SE per feature per score ----
  summary_data <- plot_data %>%
    group_by(feature, cvd_label, direction) %>%
    summarise(
      mean_imp = mean(importance, na.rm = TRUE),
      se_imp = sd(importance, na.rm = TRUE) / sqrt(sum(!is.na(importance))),
      n_folds = n(),
      .groups = "drop"
    ) %>%
    mutate(
      signed_mean = ifelse(direction == "Risk-Increasing", mean_imp, -mean_imp),
      ymin = ifelse(direction == "Risk-Increasing",
                    mean_imp - se_imp, -(mean_imp + se_imp)),
      ymax = ifelse(direction == "Risk-Increasing",
                    mean_imp + se_imp, -(mean_imp - se_imp))
    )
  
  # Sign fold-level points
  plot_data <- plot_data %>%
    mutate(signed_importance = ifelse(direction == "Risk-Increasing", importance, -importance))
  
  # ---- Join significance from aggregated data ----
  if (result_type == "score_plus_block") {
    agg_filtered <- aggregated_df %>%
      filter(grepl(paste0("_specific\\+", block_name, "$"), dataset_name))
  } else {
    agg_filtered <- aggregated_df %>%
      filter(dataset_name == block_name)
  }
  
  sig_data <- agg_filtered %>%
    mutate(
      cvd_label = recode(cvd_score,
                         "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                         "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD"),
      sig_star = case_when(
        p_value_bh < 0.001 ~ "***", p_value_bh < 0.01 ~ "**",
        p_value_bh < 0.05 ~ "*", TRUE ~ ""
      )
    ) %>%
    select(feature, cvd_label, sig_star)
  
  summary_data <- summary_data %>%
    left_join(sig_data, by = c("feature", "cvd_label")) %>%
    mutate(sig_star = ifelse(is.na(sig_star), "", sig_star))
  
  # ---- Feature ordering ----
  feat_inc <- summary_data %>%
    filter(direction == "Risk-Increasing") %>%
    group_by(feature) %>% summarise(max_imp = max(mean_imp), .groups = "drop") %>%
    arrange(desc(max_imp)) %>% pull(feature)
  
  feat_dec <- summary_data %>%
    filter(direction == "Risk-Decreasing") %>%
    group_by(feature) %>% summarise(max_imp = max(mean_imp), .groups = "drop") %>%
    arrange(desc(max_imp)) %>% pull(feature)
  
  feat_dec <- setdiff(feat_dec, feat_inc)
  feature_levels <- c(feat_inc, feat_dec)
  
  clean_name <- function(x) gsub("_", " ", gsub("\\.", " ", x))
  clean_levels <- clean_name(feature_levels)
  
  summary_data <- summary_data %>%
    mutate(feature_clean = factor(clean_name(feature), levels = clean_levels)) %>%
    filter(feature %in% feature_levels)
  
  plot_data <- plot_data %>%
    mutate(feature_clean = factor(clean_name(feature), levels = clean_levels)) %>%
    filter(feature %in% feature_levels)
  
  n_features <- length(feature_levels)
  
  # ---- No-gap positioning ----
  bar_width <- 0.18
  gap_between_features <- 0.3
  
  position_rows <- list()
  current_x <- 1
  feature_centers <- numeric(n_features)
  
  for (f_idx in seq_along(feature_levels)) {
    feat <- feature_levels[f_idx]
    scores_present <- summary_data %>%
      filter(feature == feat) %>% arrange(cvd_label) %>%
      pull(cvd_label) %>% as.character()
    
    n_bars <- length(scores_present)
    if (n_bars == 0) next
    
    group_width <- n_bars * bar_width + (n_bars - 1) * 0.02
    start_x <- current_x - group_width / 2 + bar_width / 2
    feature_centers[f_idx] <- current_x
    
    for (b_idx in seq_along(scores_present)) {
      x_pos <- start_x + (b_idx - 1) * (bar_width + 0.02)
      position_rows[[length(position_rows) + 1]] <- tibble(
        feature = feat, cvd_label = scores_present[b_idx], x_pos = x_pos
      )
    }
    current_x <- current_x + gap_between_features + group_width / 2 + bar_width
  }
  
  positions <- bind_rows(position_rows) %>%
    mutate(cvd_label = factor(cvd_label, levels = cvd_order))
  
  feature_tick_positions <- tibble(
    feature = feature_levels,
    feature_clean = factor(clean_levels, levels = clean_levels),
    x_tick = feature_centers
  )
  
  summary_data <- summary_data %>%
    left_join(positions, by = c("feature", "cvd_label"))
  plot_data <- plot_data %>%
    left_join(positions, by = c("feature", "cvd_label"))
  
  # Star positions
  summary_data <- summary_data %>%
    mutate(
      star_y = ifelse(direction == "Risk-Increasing",
                      ymax + max(abs(summary_data$signed_mean), na.rm = TRUE) * 0.05,
                      ymin - max(abs(summary_data$signed_mean), na.rm = TRUE) * 0.05)
    )
  
  # ---- Dynamic sizing ----
  plot_width <- if (n_features > 50) {
    max(20, 4 + n_features * 0.4)
  } else {
    max(10, 4 + n_features * 0.8)
  }
  
  x_text_size <- if (n_features > 80) 5 else if (n_features > 50) 6 else if (n_features > 30) 9 else 12
  
  type_prefix <- if (result_type == "score_plus_block") "Score-Specific + " else ""
  covariate_label <- paste(COVARIATE_COLS, collapse = ", ")
  
  # ---- Build plot ----
  p <- ggplot() +
    geom_col(data = summary_data,
             aes(x = x_pos, y = signed_mean, fill = cvd_label),
             width = bar_width, colour = "white", linewidth = 0.2) +
    geom_errorbar(data = summary_data %>% filter(!is.na(se_imp) & se_imp > 0),
                  aes(x = x_pos, ymin = ymin, ymax = ymax),
                  width = bar_width * 0.5, linewidth = 0.4, colour = "grey30") +
    geom_point(data = plot_data,
               aes(x = x_pos, y = signed_importance, fill = cvd_label),
               shape = 21, size = 1.2, alpha = 0.6, colour = "grey30",
               position = position_jitter(width = bar_width * 0.12, height = 0, seed = 42)) +
    geom_text(data = summary_data %>% filter(sig_star != ""),
              aes(x = x_pos, y = star_y, label = sig_star),
              size = 3.5, fontface = "bold", colour = "black") +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey30", linewidth = 0.5) +
    scale_fill_manual(values = cvd_colors, name = "CVD Score", drop = FALSE) +
    scale_x_continuous(
      breaks = feature_tick_positions$x_tick,
      labels = feature_tick_positions$feature_clean,
      expand = expansion(add = 0.8)
    ) +
    labs(
      title = sprintf("Permutation Feature Importance: %s%s", type_prefix, block_label),
      subtitle = paste0(filter_description,
                        " | Mean \u00b1 SE | 1000 permutations/fold | Adjusted for ",
                        covariate_label),
      x = NULL,
      y = expression(Delta ~ "MAE (signed by effect direction)")
      ) +
    annotate("text", x = min(feature_tick_positions$x_tick) - 0.5, y = Inf,
             label = "Risk-Increasing  \u2191",
             hjust = 0, vjust = 1.5, size = 3.5, color = "grey40", fontface = "italic") +
    annotate("text", x = min(feature_tick_positions$x_tick) - 0.5, y = -Inf,
             label = "Risk-Decreasing  \u2193",
             hjust = 0, vjust = -0.5, size = 3.5, color = "grey40", fontface = "italic") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(size = 9, color = "gray40", hjust = 0.5),
      axis.text.x = element_text(size = x_text_size, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 11),
      axis.title.y = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  # ---- Sociodemographic reference level legend ----
  if (block_name == "sociodemographics_lifestyle") {
    refs <- get_sociodemographic_references()
    features_in_plot <- feature_levels
    refs_to_show <- character()
    for (ref_name in names(refs)) {
      pattern <- gsub(" ", "[_ .]", tolower(ref_name))
      if (any(grepl(pattern, tolower(features_in_plot)))) {
        refs_to_show <- c(refs_to_show, paste0(ref_name, ": ", refs[ref_name]))
      }
    }
    if (length(refs_to_show) > 0) {
      ref_text <- paste("Reference levels |", paste(refs_to_show, collapse = " | "))
      p <- p + labs(caption = ref_text) +
        theme(plot.caption = element_text(size = 7.5, color = "grey40", hjust = 0,
                                          face = "italic", margin = margin(t = 8)))
    }
  }
  
  # Footnote about excluded combos
  if (n_excluded_combos > 0) {
    exclude_note <- sprintf("%d feature-score combinations excluded (appeared in <%d folds)",
                            n_excluded_combos, min_nonzero_folds)
    existing_caption <- p$labels$caption
    new_caption <- if (!is.null(existing_caption)) {
      paste0(existing_caption, "\n", exclude_note)
    } else exclude_note
    
    p <- p + labs(caption = new_caption) +
      theme(plot.caption = element_text(size = 7.5, color = "grey40", hjust = 0,
                                        face = "italic", margin = margin(t = 8)))
  }
  
  list(plot = p, n_features = n_features, plot_width = plot_width,
       n_excluded_combos = n_excluded_combos)
}


#' Create and save bar plots for all dataset blocks
create_all_importance_barplots <- function(importance_results,
                                           coefficient_directions,
                                           blocks = c("all_data", "lipids", "fatty_acids",
                                                      "urine_nmr", "body_composition",
                                                      "clinical_risk_factors",
                                                      "sociodemographics_lifestyle"),
                                           output_dir = "importance_barplots",
                                           result_type = "common",
                                           filter_mode = "positive_delta_3folds",
                                           dataset_labels = NULL) {
  
  full_output_dir <- save_output(output_dir)
  dir.create(full_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (is.null(dataset_labels)) {
    dataset_labels <- c(
      "all_data" = "Full predictor set", "lipids" = "Lipids",
      "fatty_acids" = "Fatty Acids", "urine_nmr" = "Urine NMR Metabolites",
      "body_composition" = "Body Composition",
      "clinical_risk_factors" = "Clinical Measurements & Supplementary Biomarkers",
      "sociodemographics_lifestyle" = "Sociodemographics & Lifestyle"
    )
  }
  
  type_label <- if (result_type == "score_plus_block") "score+block" else "blocks only"
  cat(sprintf("\n=== Extracting per-fold importance data (%s) ===\n", type_label))
  fold_data <- extract_per_fold_data(importance_results)
  cat(sprintf("  Total fold-level records: %d\n", nrow(fold_data)))
  
  aggregated_df <- importance_results$aggregated
  
  cat(sprintf("\nCreating bar plots for %d blocks (%s, filter: %s)...\n",
              length(blocks), type_label, filter_mode))
  all_plots <- list()
  
  for (block in blocks) {
    label <- ifelse(block %in% names(dataset_labels), dataset_labels[block], block)
    cat(sprintf("  %s: ", label))
    
    result <- create_importance_barplot(
      fold_data = fold_data,
      coefficient_df = coefficient_directions,
      aggregated_df = aggregated_df,
      block_name = block,
      block_label = label,
      result_type = result_type,
      filter_mode = filter_mode
    )
    
    if (is.null(result)) {
      cat("No features passing filter\n")
      next
    }
    
    cat(sprintf("%d features (%d combos excluded)\n",
                result$n_features, result$n_excluded_combos))
    
    ggsave(
      file.path(full_output_dir, sprintf("importance_barplot_%s.pdf", block)),
      result$plot,
      width = min(50, result$plot_width),
      height = 7,
      limitsize = FALSE
    )
    all_plots[[block]] <- result
  }
  
  cat(sprintf("\n\u2713 Bar plots saved to %s/\n", output_dir))
  all_plots
}


# ============================================
# 21. RUN ALL BAR PLOTS
# ============================================

filter_modes <- c("all_selected", "significant", "positive_delta", "positive_delta_3folds")

for (mode in filter_modes) {
  cat(sprintf("\n########## Filter mode: %s ##########\n", mode))
  
  cat(sprintf("\n--- COMMON (blocks only) | %s ---\n", mode))
  create_all_importance_barplots(
    importance_results = importance_results_common,
    coefficient_directions = coefficient_directions_common,
    output_dir = sprintf("barplots_common_%s", mode),
    result_type = "common",
    filter_mode = mode
  )
  
  cat(sprintf("\n--- SCORE+BLOCK | %s ---\n", mode))
  create_all_importance_barplots(
    importance_results = importance_results_full,
    coefficient_directions = coefficient_directions_full,
    output_dir = sprintf("barplots_score_plus_block_%s", mode),
    result_type = "score_plus_block",
    filter_mode = mode
  )
}

cat("\n\u2713 All bar plots created (4 modes \u00d7 2 types = 8 output folders)\n")


# ============================================
# 21b. COMBINED FIGURE — COMMON positive_delta_3folds
# ============================================

library(patchwork)

create_combined_importance_figure <- function(importance_results,
                                              coefficient_directions,
                                              blocks = NULL,
                                              result_type = "common",
                                              filter_mode = "positive_delta_3folds",
                                              output_file = "combined_importance_common_positive_delta_3folds.pdf") {
  
  fold_data <- extract_per_fold_data(importance_results)
  aggregated_df <- importance_results$aggregated
  
  dataset_labels <- c(
    "all_data" = "Full predictor set", "lipids" = "Lipids",
    "fatty_acids" = "Fatty Acids", "urine_nmr" = "Urine NMR Metabolites",
    "body_composition" = "Body Composition",
    "clinical_risk_factors" = "Clinical Measurements & Supplementary Biomarkers",
    "sociodemographics_lifestyle" = "Sociodemographics & Lifestyle"
  )
  
  if (is.null(blocks)) {
    blocks <- c("all_data", "lipids", "fatty_acids", "urine_nmr",
                "body_composition", "clinical_risk_factors",
                "sociodemographics_lifestyle")
  }
  
  # Collect plots that have results
  plot_list <- list()
  feature_counts <- c()
  
  for (block in blocks) {
    label <- ifelse(block %in% names(dataset_labels), dataset_labels[block], block)
    
    result <- create_importance_barplot(
      fold_data = fold_data,
      coefficient_df = coefficient_directions,
      aggregated_df = aggregated_df,
      block_name = block,
      block_label = label,
      result_type = result_type,
      filter_mode = filter_mode
    )
    
    if (!is.null(result)) {
      # Adjust individual plots for combined figure:
      # - smaller title, remove subtitle, shrink axis text, remove legend
      p_adjusted <- result$plot +
        labs(subtitle = NULL) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          axis.text.x = element_text(
            size = if (result$n_features > 30) 5 else if (result$n_features > 15) 7 else 9,
            angle = 45, hjust = 1
          ),
          axis.text.y = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          legend.position = "none",
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
        )
      
      plot_list[[block]] <- p_adjusted
      feature_counts <- c(feature_counts, result$n_features)
      names(feature_counts)[length(feature_counts)] <- block
    }
  }
  
  if (length(plot_list) == 0) {
    cat("No plots to combine\n")
    return(NULL)
  }
  
  cat(sprintf("Combining %d panels: %s\n",
              length(plot_list),
              paste(sprintf("%s (%d feat)", names(feature_counts), feature_counts), collapse = ", ")))
  
  # ---- Layout based on relative feature counts ----
  n_panels <- length(plot_list)
  
  # Use relative widths proportional to feature count (with minimum)
  rel_widths <- pmax(feature_counts / max(feature_counts), 0.3)
  
  # Arrange: put wider panels on top row, narrower on bottom
  order_idx <- order(feature_counts, decreasing = TRUE)
  ordered_plots <- plot_list[order_idx]
  ordered_widths <- rel_widths[order_idx]
  
  if (n_panels <= 2) {
    # Side by side
    combined <- wrap_plots(ordered_plots, nrow = 1, widths = ordered_widths)
  } else if (n_panels <= 4) {
    # 2 rows: larger panels on top, smaller on bottom
    n_top <- ceiling(n_panels / 2)
    n_bot <- n_panels - n_top
    
    top_plots <- ordered_plots[1:n_top]
    bot_plots <- ordered_plots[(n_top + 1):n_panels]
    top_widths <- ordered_widths[1:n_top]
    bot_widths <- ordered_widths[(n_top + 1):n_panels]
    
    top_row <- wrap_plots(top_plots, nrow = 1, widths = top_widths)
    bot_row <- wrap_plots(bot_plots, nrow = 1, widths = bot_widths)
    
    combined <- top_row / bot_row
  } else {
    # 3 rows for 5+ panels
    n_per_row <- ceiling(n_panels / 3)
    rows <- list()
    for (r in 1:3) {
      start <- (r - 1) * n_per_row + 1
      end <- min(r * n_per_row, n_panels)
      if (start > n_panels) break
      rows[[r]] <- wrap_plots(ordered_plots[start:end], nrow = 1,
                              widths = ordered_widths[start:end])
    }
    combined <- Reduce("/", rows)
  }
  
  # Add shared legend from first plot
  legend_plot <- plot_list[[1]] +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(size = 9))
  shared_legend <- cowplot::get_legend(legend_plot)
  
  # Add overall title, subtitle, and shared legend
  covariate_label <- paste(COVARIATE_COLS, collapse = ", ")
  type_prefix <- if (result_type == "score_plus_block") "Score-Specific + " else ""
  
  final <- combined /
    wrap_elements(shared_legend) +
    plot_layout(heights = c(rep(1, if (n_panels <= 2) 1 else if (n_panels <= 4) 2 else 3), 0.08)) +
    plot_annotation(
      title = expression(bold("Permutation Feature Importance (positive" ~ Delta * "MAE in" >= "3 folds)")),
      subtitle = sprintf("Adjusted for %s | Mean +/- SE | 1000 permutations/fold",
                         covariate_label),
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5)
      )
    )
  
  # ---- Save ----
  total_features <- sum(feature_counts)
  fig_width <- min(50, max(14, total_features * 0.25))
  fig_height <- if (n_panels <= 2) 8 else if (n_panels <= 4) 14 else 20
  
  output_path <- save_output(output_file)
  ggsave(output_path, final, width = fig_width, height = fig_height, limitsize = FALSE)
  
  cat(sprintf("\u2713 Combined figure saved to %s (%.0f x %.0f inches)\n",
              output_file, fig_width, fig_height))
  
  invisible(final)
}


# ---- Run ----
cat("\n=== Creating combined importance figure (common, positive_delta_3folds) ===\n")

combined_fig <- create_combined_importance_figure(
  importance_results = importance_results_common,
  coefficient_directions = coefficient_directions_common,
  blocks = c("all_data", "urine_nmr", "body_composition", "sociodemographics_lifestyle"),
  result_type = "common",
  filter_mode = "positive_delta_3folds",
  output_file = "combined_importance_common_positive_delta_3folds.pdf"
)




# ============================================
# FIGURE 1: COMMON PREDICTORS — positive_delta_3folds
# FIGURE 2: SCORE-SPECIFIC + BLOCK — positive_delta_3folds
# ============================================

library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)


# ================================================================
# SHARED PANEL BUILDER
# ================================================================

build_importance_panel <- function(fold_data, coefficient_directions, aggregated_df,
                                   block_name, result_type,
                                   score_colours, cvd_order,
                                   remove_score_specific = FALSE,
                                   score_specific_regex = NULL) {
  
  # Filter fold data to this block
  if (result_type == "score_plus_block") {
    plot_data <- fold_data %>%
      filter(grepl(paste0("_specific\\+", block_name, "$"), dataset_name))
  } else {
    plot_data <- fold_data %>% filter(dataset_name == block_name)
  }
  if (nrow(plot_data) == 0) return(NULL)
  
  # Remove score-specific features for Row B of figure 2
  if (remove_score_specific && !is.null(score_specific_regex)) {
    n_before <- length(unique(plot_data$feature))
    plot_data <- plot_data %>% filter(!grepl(score_specific_regex, feature))
    n_after <- length(unique(plot_data$feature))
    cat(sprintf("    Removed %d score-specific features, %d remaining\n",
                n_before - n_after, n_after))
    if (nrow(plot_data) == 0) return(NULL)
  }
  
  # Join coefficient directions
  plot_data <- plot_data %>%
    left_join(
      coefficient_directions %>% select(cvd_score, dataset_name, feature, direction),
      by = c("cvd_score", "dataset_name", "feature")
    ) %>%
    filter(!is.na(direction), direction != "Zero")
  if (nrow(plot_data) == 0) return(NULL)
  
  # Recode score labels
  plot_data <- plot_data %>%
    mutate(
      cvd_label = recode(cvd_score,
                         "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                         "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD"),
      cvd_label = factor(cvd_label, levels = cvd_order)
    )
  
  # Filter: positive delta in >=3 folds
  plot_data <- plot_data %>% filter(importance > 0)
  fold_counts <- plot_data %>%
    group_by(feature, cvd_label) %>%
    summarise(n_nonzero = n(), .groups = "drop")
  valid_combos <- fold_counts %>% filter(n_nonzero >= 3)
  n_excluded <- nrow(fold_counts) - nrow(valid_combos)
  plot_data <- plot_data %>% semi_join(valid_combos, by = c("feature", "cvd_label"))
  if (nrow(plot_data) == 0) return(NULL)
  
  # Summarise
  summary_data <- plot_data %>%
    group_by(feature, cvd_label, direction) %>%
    summarise(
      mean_imp = mean(importance, na.rm = TRUE),
      se_imp = sd(importance, na.rm = TRUE) / sqrt(sum(!is.na(importance))),
      n_folds = n(), .groups = "drop"
    ) %>%
    mutate(
      signed_mean = ifelse(direction == "Risk-Increasing", mean_imp, -mean_imp),
      ymin = ifelse(direction == "Risk-Increasing", mean_imp - se_imp, -(mean_imp + se_imp)),
      ymax = ifelse(direction == "Risk-Increasing", mean_imp + se_imp, -(mean_imp - se_imp))
    )
  
  plot_data <- plot_data %>%
    mutate(signed_importance = ifelse(direction == "Risk-Increasing", importance, -importance))
  
  # Significance stars
  if (result_type == "score_plus_block") {
    agg_filtered <- aggregated_df %>%
      filter(grepl(paste0("_specific\\+", block_name, "$"), dataset_name))
  } else {
    agg_filtered <- aggregated_df %>% filter(dataset_name == block_name)
  }
  sig_data <- agg_filtered %>%
    mutate(
      cvd_label = recode(cvd_score,
                         "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                         "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD"),
      sig_star = case_when(
        p_value_bh < 0.001 ~ "***", p_value_bh < 0.01 ~ "**",
        p_value_bh < 0.05 ~ "*", TRUE ~ "")
    ) %>% select(feature, cvd_label, sig_star)
  
  summary_data <- summary_data %>%
    left_join(sig_data, by = c("feature", "cvd_label")) %>%
    mutate(sig_star = ifelse(is.na(sig_star), "", sig_star))
  
  # Feature ordering
  feat_inc <- summary_data %>% filter(direction == "Risk-Increasing") %>%
    group_by(feature) %>% summarise(max_imp = max(mean_imp), .groups = "drop") %>%
    arrange(desc(max_imp)) %>% pull(feature)
  feat_dec <- summary_data %>% filter(direction == "Risk-Decreasing") %>%
    group_by(feature) %>% summarise(max_imp = max(mean_imp), .groups = "drop") %>%
    arrange(desc(max_imp)) %>% pull(feature)
  feat_dec <- setdiff(feat_dec, feat_inc)
  feature_levels <- c(feat_inc, feat_dec)
  
  # ---- Feature display name lookup ----
  # Maps raw feature names (with underscores) to clean display names.
  # Features not listed here fall back to default cleaning (underscores → spaces, first letter capitalised).
  feature_display_names <- c(
    # Clinical / score-specific
    "systolic"                            = "Systolic BP",
    "blood_pressure_treatment_Yes"        = "Blood Pressure Treatment",
    "blood_pressure_treatment"            = "Blood Pressure Treatment",
    "SmokingStatusQRISK3"                 = "Smoking",
    "mean_Total_Cholesterol_mg_dl"        = "Total Cholesterol",
    "ratio_chol_hdl"                      = "Cholesterol/HDL",
    "creatinine"                          = "Creatinine",
    "mean_LDL_mg_dl"                      = "LDL",
    "LDL.mg.dl"                           = "LDL",
    "TAG.mg.dl"                           = "TAG",
    "mean_HDL_mg_dl"                      = "HDL",
    "mean_hrt"                            = "Mean Heart Rate",
    
    # Lipids
    "lpe_20_4"                            = "LPE(20:4)",
    "pe_o_22_2_20_4"                      = "PE-O(22:2/20:4)",
    
    # Fatty acids
    "rbc_linoleic_18_2n6"                 = "RBC Linoleic Acid (18:2n-6)",
    "d_nervonic_24_1n9"                   = "d Nervonic Acid (24:1n-9)",
    "d_dgla_20_3n6"                       = "d DGLA (20:3n-6)",
    "rbc_dpa_22_5n3"                      = "RBC Docosapentaenoic Acid (22:5n-3)",
    
    # Urine NMR
    "allantoin"                           = "Allantoin",
    "lactic_acid"                         = "Lactic Acid",
    "creatine"                            = "Creatine",
    "proline_betaine"                     = "Proline Betaine",
    
    # Body composition
    "ecw_tbw"                                              = "ECW/TBW",
    "bmr__basal_metabolic_rate_"                           = "Basal Metabolic Rate",
    "svr_skeletal_muscle_mass_visceral_fat_area_ratio_"    = "Skeletal Muscle Mass/Visceral Fat Area",
    "smi__skeletal_muscle_index_"                          = "Skeletal Muscle Index",
    "vfa__visceral_fat_area_"                              = "Visceral Fat Area",
    "bfm_of_trunk"                                         = "Body Fat Mass of Trunk",
    "Fat.free.mass"                                        = "Fat-Free Mass",
    "Body.fat"                                             = "Body Fat",
    "measured_circumference_of_right_thigh"                = "Right Thigh Circumference",
    "measured_circumference_of_left_thigh"                 = "Left Thigh Circumference",
    "sleep_quality"                                        = "Sleep Quality",
    
    # Sociodemographics
    "naps_during_day_yes"                                  = "Naps During Day",
    "Education_Level_4"                                    = "Doctoral",
    "Education_Level_3"                                    = "University (Bachelor or Equivalent)",
    "annual_net_salary_4"                                  = "Annual Net Salary >34,400\u201343,000",
    "Living_Status_Living.with.a.partner"                  = "Living With a Partner",
    "Employment_Status_Retired"                            = "Retired"
  )
  
  # Capitalise first letter of each word as fallback
  capitalise_words <- function(s) {
    gsub("(^|\\s)(\\w)", "\\1\\U\\2", s, perl = TRUE)
  }
  
  clean_name <- function(x) {
    # Build a lowercase lookup for case-insensitive matching
    lookup_lower <- setNames(feature_display_names, tolower(names(feature_display_names)))
    sapply(x, function(feat) {
      if (feat %in% names(feature_display_names)) {
        feature_display_names[feat]
      } else if (tolower(feat) %in% names(lookup_lower)) {
        lookup_lower[tolower(feat)]
      } else {
        capitalise_words(gsub("_", " ", gsub("\\.", " ", feat)))
      }
    }, USE.NAMES = FALSE)
  }
  
  clean_levels <- clean_name(feature_levels)
  
  summary_data <- summary_data %>%
    mutate(feature_clean = factor(clean_name(feature), levels = clean_levels)) %>%
    filter(feature %in% feature_levels)
  plot_data <- plot_data %>%
    mutate(feature_clean = factor(clean_name(feature), levels = clean_levels)) %>%
    filter(feature %in% feature_levels)
  
  n_features <- length(feature_levels)
  
  # No-gap positioning
  bar_width <- 0.18
  gap_between_features <- 0.3
  position_rows <- list()
  current_x <- 1
  feature_centers <- numeric(n_features)
  
  for (f_idx in seq_along(feature_levels)) {
    feat <- feature_levels[f_idx]
    scores_present <- summary_data %>%
      filter(feature == feat) %>% arrange(cvd_label) %>%
      pull(cvd_label) %>% as.character()
    n_bars <- length(scores_present)
    if (n_bars == 0) next
    group_width <- n_bars * bar_width + (n_bars - 1) * 0.02
    start_x <- current_x - group_width / 2 + bar_width / 2
    feature_centers[f_idx] <- current_x
    for (b_idx in seq_along(scores_present)) {
      x_pos <- start_x + (b_idx - 1) * (bar_width + 0.02)
      position_rows[[length(position_rows) + 1]] <- tibble(
        feature = feat, cvd_label = scores_present[b_idx], x_pos = x_pos)
    }
    current_x <- current_x + gap_between_features + group_width / 2 + bar_width
  }
  
  positions <- bind_rows(position_rows) %>%
    mutate(cvd_label = factor(cvd_label, levels = cvd_order))
  feature_tick_positions <- tibble(
    feature = feature_levels,
    feature_clean = factor(clean_levels, levels = clean_levels),
    x_tick = feature_centers)
  
  summary_data <- summary_data %>% left_join(positions, by = c("feature", "cvd_label"))
  plot_data <- plot_data %>% left_join(positions, by = c("feature", "cvd_label"))
  
  summary_data <- summary_data %>%
    mutate(star_y = ifelse(direction == "Risk-Increasing",
                           ymax + max(abs(signed_mean), na.rm = TRUE) * 0.05,
                           ymin - max(abs(signed_mean), na.rm = TRUE) * 0.05))
  
  # Bigger feature text, adaptive
  x_text_size <- if (n_features > 80) 10 else if (n_features > 50) 11 else if (n_features > 30) 12 else if (n_features > 15) 13 else 14
  
  # Build plot — no title, no subtitle (will be in figure legend)
  p <- ggplot() +
    geom_col(data = summary_data, aes(x = x_pos, y = signed_mean, fill = cvd_label),
             width = bar_width, colour = "white", linewidth = 0.2) +
    geom_errorbar(data = summary_data %>% filter(!is.na(se_imp) & se_imp > 0),
                  aes(x = x_pos, ymin = ymin, ymax = ymax),
                  width = bar_width * 0.5, linewidth = 0.4, colour = "grey30") +
    geom_point(data = plot_data, aes(x = x_pos, y = signed_importance, fill = cvd_label),
               shape = 21, size = 1.0, alpha = 0.5, colour = "grey30",
               position = position_jitter(width = bar_width * 0.12, height = 0, seed = 42)) +
    geom_text(data = summary_data %>% filter(sig_star != ""),
              aes(x = x_pos, y = star_y, label = sig_star),
              size = 4, fontface = "bold", colour = "black") +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey30", linewidth = 0.5) +
    scale_fill_manual(values = score_colours, name = "CVD Risk Score",
                      limits = cvd_order, drop = FALSE) +
    scale_x_continuous(breaks = feature_tick_positions$x_tick,
                       labels = feature_tick_positions$feature_clean,
                       expand = expansion(add = 0.5)) +
    labs(x = NULL, y = expression(Delta ~ "MAE (signed by effect direction)")) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.02))) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(size = x_text_size, angle = 30, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 5, r = 8, b = 2, l = 8))
  
  y_range <- range(c(summary_data$ymin, summary_data$ymax,
                     plot_data$signed_importance), na.rm = TRUE)
  list(plot = p, n_features = n_features, y_range = y_range, n_excluded = n_excluded)
}


# ================================================================
# SHARED ASSEMBLY HELPER
# ================================================================

assemble_combined_figure <- function(panel_a, panels_b,
                                     ylim_a, ylim_b,
                                     score_colours, cvd_order,
                                     output_file, fig_height = 12) {
  
  # Apply y-limits with panel tags A, B, C, D...
  plot_a <- panel_a$plot + coord_cartesian(ylim = ylim_a, clip = "off") +
    labs(tag = "A") +
    theme(plot.tag = element_text(face = "bold", size = 16))
  
  b_list <- list()
  panel_letters <- LETTERS[2:(1 + length(panels_b))]  # B, C, D, ...
  for (i in seq_along(panels_b)) {
    p_b <- panels_b[[i]]$plot +
      coord_cartesian(ylim = ylim_b, clip = "off") +
      labs(tag = panel_letters[i]) +
      theme(plot.tag = element_text(face = "bold", size = 16))
    # Only keep y-axis title on first B panel
    if (i > 1) {
      p_b <- p_b + labs(y = NULL)
    }
    b_list[[i]] <- p_b
  }
  
  # Row B widths proportional to feature count
  b_counts <- sapply(panels_b, function(x) x$n_features)
  b_widths <- pmax(b_counts / max(b_counts), 0.25)
  
  row_a <- plot_a
  row_b <- wrap_plots(b_list, nrow = 1, widths = b_widths)
  
  # Legend (no direction labels — will be in figure caption)
  dummy_data <- tibble(x = 1:4, y = 1:4, cvd_label = factor(cvd_order, levels = cvd_order))
  legend_plot <- ggplot(dummy_data, aes(x, y, fill = cvd_label)) + geom_col() +
    scale_fill_manual(values = score_colours, limits = cvd_order,
                      name = "CVD Risk Score", drop = FALSE) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 11),
          legend.key.size = unit(0.5, "cm"))
  legend_row <- cowplot::get_legend(legend_plot)
  
  # Height ratio
  height_ratio_b <- diff(ylim_b) / diff(ylim_a)
  
  combined <- (row_a / row_b / wrap_elements(legend_row)) +
    plot_layout(heights = c(1, height_ratio_b, 0.06))
  
  # Save
  total_b <- sum(b_counts)
  fig_width <- min(50, max(16, max(panel_a$n_features, total_b) * 0.35))
  
  output_path <- save_output(output_file)
  ggsave(output_path, combined, width = fig_width, height = fig_height, limitsize = FALSE)
  cat(sprintf("\n\u2713 Figure saved to %s (%.0f x %.0f in)\n", output_file, fig_width, fig_height))
  invisible(combined)
}


# ================================================================
# FIGURE 1: COMMON PREDICTORS
# ================================================================

create_combined_figure_common <- function(importance_results,
                                          coefficient_directions,
                                          output_file = "combined_importance_common_v2.pdf") {
  
  score_colours <- c("ASCVD" = "#F8766D", "Framingham" = "#7CAE00",
                     "QRISK3" = "#C77CFF", "SCORE2" = "#00BFC4")
  cvd_order <- c("QRISK3", "SCORE2", "Framingham", "ASCVD")
  
  fold_data <- extract_per_fold_data(importance_results)
  aggregated_df <- importance_results$aggregated
  
  cat("\n--- Figure 1: Building Row A ---\n")
  panel_a <- build_importance_panel(fold_data, coefficient_directions, aggregated_df,
                                    "all_data", "common", score_colours, cvd_order)
  if (is.null(panel_a)) { cat("Row A: no results\n"); return(NULL) }
  cat(sprintf("  all_data: %d features\n", panel_a$n_features))
  
  cat("\n--- Figure 1: Building Row B ---\n")
  row_b_blocks <- c("body_composition", "sociodemographics_lifestyle", "urine_nmr")
  panels_b <- list()
  for (block in row_b_blocks) {
    res <- build_importance_panel(fold_data, coefficient_directions, aggregated_df,
                                  block, "common", score_colours, cvd_order)
    if (!is.null(res)) {
      panels_b[[block]] <- res
      cat(sprintf("  %s: %d features\n", block, res$n_features))
    }
  }
  if (length(panels_b) == 0) { cat("No Row B panels\n"); return(NULL) }
  
  assemble_combined_figure(
    panel_a, panels_b,
    ylim_a = c(-0.4, 0.4),
    ylim_b = c(-0.6, 0.4),
    score_colours, cvd_order,
    output_file, fig_height = 12
  )
}


# ================================================================
# FIGURE 2: SCORE-SPECIFIC + BLOCK
# ================================================================

create_combined_figure_score_plus_block <- function(importance_results,
                                                    coefficient_directions,
                                                    output_file = "combined_importance_score_plus_block_v2.pdf") {
  
  score_colours <- c("ASCVD" = "#F8766D", "Framingham" = "#7CAE00",
                     "QRISK3" = "#C77CFF", "SCORE2" = "#00BFC4")
  cvd_order <- c("QRISK3", "SCORE2", "Framingham", "ASCVD")
  
  # Score-specific feature stems
  score_specific_stems <- c(
    "Sex", "Age",
    "EthnicityCodeQRISK3", "SmokingStatusQRISK3",
    "diabetes2", "Diabetes",
    "Weight_kg", "Height_cm",
    "blood_pressure_treatment",
    "Severe_mental_illness",
    "ratio_chol_hdl", "systolic", "townsend",
    "race_ascvd",
    "mean_Total_Cholesterol_mg_dl", "mean_HDL_mg_dl", "mean_LDL_mg_dl",
    "Risk.region"
  )
  escaped_stems <- gsub("\\.", "\\\\.", score_specific_stems)
  score_specific_regex <- paste0("^(", paste(escaped_stems, collapse = "|"), ")")
  
  fold_data <- extract_per_fold_data(importance_results)
  aggregated_df <- importance_results$aggregated
  
  cat("\n--- Figure 2: Building Row A (all features) ---\n")
  panel_a <- build_importance_panel(fold_data, coefficient_directions, aggregated_df,
                                    "all_data", "score_plus_block", score_colours, cvd_order,
                                    remove_score_specific = FALSE)
  if (is.null(panel_a)) { cat("Row A: no results\n"); return(NULL) }
  cat(sprintf("  all_data: %d features\n", panel_a$n_features))
  
  cat("\n--- Figure 2: Building Row B (score-specific excluded) ---\n")
  all_b_blocks <- c("lipids", "fatty_acids", "urine_nmr", "body_composition",
                    "clinical_risk_factors", "sociodemographics_lifestyle")
  panels_b <- list()
  for (block in all_b_blocks) {
    cat(sprintf("  %s:\n", block))
    res <- build_importance_panel(fold_data, coefficient_directions, aggregated_df,
                                  block, "score_plus_block", score_colours, cvd_order,
                                  remove_score_specific = TRUE,
                                  score_specific_regex = score_specific_regex)
    if (!is.null(res)) {
      panels_b[[block]] <- res
      cat(sprintf("    -> %d block-specific features in plot\n", res$n_features))
    } else {
      cat("    -> No block-specific features passing filter (skipped)\n")
    }
  }
  if (length(panels_b) == 0) { cat("No Row B panels\n"); return(NULL) }
  
  # Row A y-axis: include negative range so labels sit near bars
  pad_a <- diff(panel_a$y_range) * 0.05
  ylim_a <- c(-0.5,
              ceiling((panel_a$y_range[2] + pad_a) * 10) / 10)
  
  cat(sprintf("\n  Y-axis: A = (%.2f, %.2f), B = (-0.25, 0.25)\n", ylim_a[1], ylim_a[2]))
  
  assemble_combined_figure(
    panel_a, panels_b,
    ylim_a = ylim_a,
    ylim_b = c(-0.25, 0.25),
    score_colours, cvd_order,
    output_file, fig_height = 12
  )
}


# ================================================================
# RUN BOTH
# ================================================================

cat("\n=== Figure 1: Common predictors ===\n")
fig_common <- create_combined_figure_common(
  importance_results = importance_results_common,
  coefficient_directions = coefficient_directions_common,
  output_file = "combined_importance_common_v2.pdf"
)

cat("\n=== Figure 2: Score-specific + block ===\n")
fig_spb <- create_combined_figure_score_plus_block(
  importance_results = importance_results_full,
  coefficient_directions = coefficient_directions_full,
  output_file = "combined_importance_score_plus_block_v2.pdf"
)