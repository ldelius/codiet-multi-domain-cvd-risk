################################################################################
### Elastic Net Regression with Covariate Adjustment (Country, Statins, Supplements unpenalised)
### Description: Fits covariate-adjusted Elastic Net models predicting CVD risk scores from multi-domain predictor sets.
### Model evaluation; permutation test and permutation-based feature importance.
### tables and figures.
### Based on code in ElaNet_foldwise_preproc.R.
### Author: Luisa Delius
################################################################################


# ============================================
# Set everything up
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
library(cowplot)

set.seed(42)

wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files/processed_data"
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

CPUS <- parallel::detectCores() - 3

# ============================================
# 1. Load data
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

# Source covariates before dropping them (df_QRISK3_input must be loaded first)
df_covariates <- df_fatty_acids_raw_full %>%
  select(Sample_ID, Statins, Supplements) %>%
  full_join(df_risk_factors_raw_full %>% select(Sample_ID, Country), by = "Sample_ID") %>%
  left_join(df_QRISK3_input %>% select(Sample_ID, Sex, Age), by = "Sample_ID") %>%
  distinct(Sample_ID, .keep_all = TRUE)

COVARIATE_COLS <- c("Country", "Statins", "Supplements")

df_fatty_acids_ElaNet <- df_fatty_acids_raw_full %>%
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
  select(-starts_with("z_"), -QRISK3_risk, -recruitment_site, -any_of(COVARIATE_COLS))

df_risk_factors_ElaNet <- df_risk_factors_raw_full %>%
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
# 2. Create predictor sets
# ============================================

attach_covariates <- function(datasets_list, covariate_cols = COVARIATE_COLS) {
  lapply(datasets_list, function(df) {
    covs_to_add <- df_covariates %>% select(Sample_ID, all_of(covariate_cols))
    df %>%
      select(-any_of(covariate_cols)) %>%
      left_join(covs_to_add, by = "Sample_ID")
  })
}

attach_covariates_single <- function(df, covariate_cols = COVARIATE_COLS) {
  covs_to_add <- df_covariates %>% select(Sample_ID, all_of(covariate_cols))
  df %>%
    select(-any_of(covariate_cols)) %>%
    left_join(covs_to_add, by = "Sample_ID")
}

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

# ============================================
# 3. Preprocessing and helper functions
# ============================================
# missingness exclusion
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

# preprocessing recipe
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

# Alpha selection by minimum RMSE
select_best_alpha <- function(cv_results, alpha_grid) {
  cv_rmse_vals <- sapply(cv_results, function(x) {
    idx <- which(x$cv_fit$lambda == x$cv_fit$lambda.1se)
    sqrt(x$cv_fit$cvm[idx])
  })
  best_idx <- which.min(cv_rmse_vals)
  list(best_idx = best_idx, best_alpha = alpha_grid[best_idx],
       best_rmse = cv_rmse_vals[best_idx], all_rmse = cv_rmse_vals)
}

# Stratified fold assignment for continuous outcomes
create_stratified_folds <- function(y, nfolds, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(y)
  # Create quantile-based strata (nfolds bins)
  strata <- cut(y, breaks = quantile(y, probs = seq(0, 1, length.out = nfolds + 1)),
                include.lowest = TRUE, labels = FALSE)
  # Handle any NA from ties at boundaries
  strata[is.na(strata)] <- sample(1:nfolds, sum(is.na(strata)), replace = TRUE)
  foldid <- integer(n)
  for (s in unique(strata)) {
    idx <- which(strata == s)
    foldid[idx] <- sample(rep(1:nfolds, length.out = length(idx)))
  }
  foldid
}

# ============================================
# 4. Permutation test function
# ============================================

compute_permutation_Q2 <- function(df_raw, y_shuffled, exclude_cols,
                                   alpha_grid = seq(0, 1, by = 0.1), nfolds = 10, seed = 42) {
  set.seed(seed)
  n <- length(y_shuffled)
  foldid <- create_stratified_folds(y_shuffled, nfolds, seed = seed)
  y_oof <- rep(NA_real_, n)
  for (k in 1:nfolds) {
    idx_test <- which(foldid == k); idx_train <- setdiff(1:n, idx_test)
    df_train <- df_raw[idx_train, , drop = FALSE]; df_test <- df_raw[idx_test, , drop = FALSE]
    y_train <- y_shuffled[idx_train]
    preprocessed <- tryCatch(preprocess_train_test(df_train, df_test, exclude_cols), error = function(e) NULL)
    if (is.null(preprocessed) || ncol(preprocessed$X_train) == 0) return(NA_real_)
    pen_factor <- ifelse(preprocessed$is_covariate, 0, 1)
    cv_results <- list(); valid_count <- 0
    for (a_idx in seq_along(alpha_grid)) {
      cv_fit <- tryCatch(cv.glmnet(x = preprocessed$X_train, y = y_train, alpha = alpha_grid[a_idx],
                                   penalty.factor = pen_factor, nfolds = nfolds,
                                   type.measure = "mse", standardize = FALSE, family = "gaussian"), error = function(e) NULL)
      if (!is.null(cv_fit)) { valid_count <- valid_count + 1; cv_results[[valid_count]] <- list(alpha = alpha_grid[a_idx], cv_fit = cv_fit) }
    }
    if (valid_count == 0) return(NA_real_)
    valid_alphas <- sapply(cv_results, function(x) x$alpha)
    alpha_selection <- select_best_alpha(cv_results, valid_alphas)
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
# 5. Helper functions
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
# 6. Core Elastic Net function
# ============================================

get_oof_predictions_clean <- function(df_raw, y, exclude_cols,
                                      covariate_cols = COVARIATE_COLS,
                                      alpha_grid = seq(0, 1, by = 0.1), nfolds = 10, seed = 42) {
  set.seed(seed)
  n <- length(y); y_mean_global <- mean(y)
  foldid_outer <- create_stratified_folds(y, nfolds, seed = seed)
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
    y_train <- y[idx_train]; y_test <- y[idx_test]
    
    preprocessed <- preprocess_train_test(df_train, df_test, exclude_cols, covariate_cols)
    X_train <- preprocessed$X_train; X_test <- preprocessed$X_test
    pen_factor <- ifelse(preprocessed$is_covariate, 0, 1)
    n_features_per_fold[k] <- ncol(X_train)
    
    cv_results_fold <- lapply(alpha_grid, function(a) {
      cv_fit <- cv.glmnet(x = X_train, y = y_train, alpha = a, penalty.factor = pen_factor,
                          nfolds = nfolds, type.measure = "mse", standardize = FALSE, family = "gaussian")
      list(alpha = a, cv_fit = cv_fit)
    })
    
    alpha_selection <- select_best_alpha(cv_results_fold, alpha_grid)
    alpha_selected_per_fold[k] <- alpha_grid[alpha_selection$best_idx]
    best_cv_fit <- cv_results_fold[[alpha_selection$best_idx]]$cv_fit
    
    y_pred_k <- as.numeric(predict(best_cv_fit, newx = X_test, s = "lambda.1se"))
    y_oof[idx_test] <- y_pred_k
    
    ss_res_k <- sum((y_test - y_pred_k)^2)
    ss_tot_k <- sum((y_test - y_mean_global)^2)
    Q2_per_fold[k] <- 1 - ss_res_k / ss_tot_k
    MAE_per_fold[k] <- mean(abs(y_test - y_pred_k))
    RMSE_per_fold[k] <- sqrt(mean((y_test - y_pred_k)^2))
    
    y_pred_train_k <- as.numeric(predict(best_cv_fit, newx = X_train, s = "lambda.1se"))
    R2_per_fold[k] <- 1 - sum((y_train - y_pred_train_k)^2) / sum((y_train - mean(y_train))^2)
  }
  
  list(predictions = y_oof, alpha_per_fold = alpha_selected_per_fold,
       n_features_per_fold = n_features_per_fold, foldid = foldid_outer,
       Q2_per_fold = Q2_per_fold, MAE_per_fold = MAE_per_fold,
       RMSE_per_fold = RMSE_per_fold, R2_per_fold = R2_per_fold)
}

run_elastic_net_model <- function(df_raw, outcome_col, cvd_score_name, dataset_name,
                                  exclude_cols = "Sample_ID",
                                  covariate_cols = COVARIATE_COLS,
                                  alpha_grid = seq(0, 1, by = 0.1),
                                  nfolds = 10, n_permutations = 100,
                                  col_miss_thresh = 0.4, row_miss_thresh = 0.8,
                                  seed = 42) {
  set.seed(seed)
  df_complete <- df_raw %>% filter(!is.na(.data[[outcome_col]]))
  y <- df_complete[[outcome_col]]
  all_exclude <- unique(c(exclude_cols, outcome_col, CVD_scores))
  
  cat(sprintf("  Filtering by missingness (col >%.0f%%, row >%.0f%%)...\n",
              col_miss_thresh * 100, row_miss_thresh * 100))
  cols_to_protect <- unique(c(all_exclude, covariate_cols))
  df_filtered <- filter_by_missingness(df_complete, cols_to_protect, col_miss_thresh, row_miss_thresh)
  
  for (cov in covariate_cols) {
    if (!cov %in% names(df_filtered) && cov %in% names(df_complete)) {
      df_filtered[[cov]] <- df_complete[[cov]][match(df_filtered$Sample_ID, df_complete$Sample_ID)]
    }
  }
  
  y <- df_filtered[[outcome_col]]; n_obs <- nrow(df_filtered)
  cat(sprintf("  Sample size after filtering: %d\n", n_obs))
  if (n_obs < nfolds * 2) warning(sprintf("Very few observations (%d) for %d-fold CV", n_obs, nfolds))
  
  # Full-data model (for coefficients)
  preprocessed_full <- preprocess_train_test(df_filtered, df_filtered, all_exclude, covariate_cols)
  X_full <- preprocessed_full$X_train; n_pred <- ncol(X_full)
  predictor_names <- colnames(X_full)
  pen_factor_full <- ifelse(preprocessed_full$is_covariate, 0, 1)
  n_covariates_in_model <- sum(preprocessed_full$is_covariate)
  n_penalized <- sum(!preprocessed_full$is_covariate)
  
  cat(sprintf("  Predictors: %d (%d penalized, %d unpenalized covariates)\n",
              n_pred, n_penalized, n_covariates_in_model))
  if (n_pred == 0) { warning("No predictors remaining after preprocessing"); return(NULL) }
  
  set.seed(seed); foldid_inner <- create_stratified_folds(y, nfolds, seed = seed)
  cv_results <- lapply(alpha_grid, function(a) {
    cv_fit <- cv.glmnet(x = X_full, y = y, alpha = a, penalty.factor = pen_factor_full,
                        nfolds = nfolds, foldid = foldid_inner,
                        type.measure = "mse", standardize = FALSE, family = "gaussian")
    list(alpha = a, cv_fit = cv_fit)
  })
  
  alpha_selection <- select_best_alpha(cv_results, alpha_grid)
  best_alpha_fulldata <- alpha_grid[alpha_selection$best_idx]
  best_cv_fit <- cv_results[[alpha_selection$best_idx]]$cv_fit
  cat(sprintf("  Alpha selection: min RMSE = %.4f at alpha = %.1f\n",
              alpha_selection$best_rmse, best_alpha_fulldata))
  
  lambda_min <- best_cv_fit$lambda.min; lambda_1se <- best_cv_fit$lambda.1se
  cv_rmse_min <- sqrt(min(best_cv_fit$cvm))
  cv_rmse_1se <- sqrt(best_cv_fit$cvm[which(best_cv_fit$lambda == lambda_1se)])
  
  final_coef <- coef(best_cv_fit, s = "lambda.1se")
  intercept <- as.numeric(final_coef[1])
  coef_vector <- as.numeric(final_coef[-1]); names(coef_vector) <- predictor_names
  n_nonzero <- sum(coef_vector != 0)
  
  # In-sample predictions (R²_Y)
  y_pred_in <- as.numeric(predict(best_cv_fit, newx = X_full, s = "lambda.1se"))
  res_in <- y - y_pred_in
  rmse_in <- sqrt(mean(res_in^2)); mae_in <- mean(abs(res_in))
  R2_Y <- 1 - sum(res_in^2) / sum((y - mean(y))^2); cor_in <- cor(y_pred_in, y)
  
  # Out-of-fold predictions (Q²_Y)
  oof_results <- get_oof_predictions_clean(df_filtered, y, all_exclude, covariate_cols, alpha_grid, nfolds, seed + 100)
  y_pred_oof <- oof_results$predictions; alpha_per_fold <- oof_results$alpha_per_fold
  Q2_per_fold <- oof_results$Q2_per_fold
  res_oof <- y - y_pred_oof
  rmse_oof <- sqrt(mean(res_oof^2)); mae_oof <- mean(abs(res_oof))
  Q2_Y <- 1 - sum(res_oof^2) / sum((y - mean(y))^2); cor_oof <- cor(y_pred_oof, y)
  Q2_fold_mean <- mean(Q2_per_fold); Q2_fold_sd <- sd(Q2_per_fold)
  MAE_per_fold <- oof_results$MAE_per_fold; RMSE_per_fold <- oof_results$RMSE_per_fold
  MAE_fold_mean <- mean(MAE_per_fold); MAE_fold_sd <- sd(MAE_per_fold)
  RMSE_fold_mean <- mean(RMSE_per_fold); RMSE_fold_sd <- sd(RMSE_per_fold)
  R2_per_fold <- oof_results$R2_per_fold
  R2_fold_mean <- mean(R2_per_fold); R2_fold_sd <- sd(R2_per_fold)
  
  # Model quality metrics
  Q2_R2_ratio <- Q2_Y / R2_Y; R2_Q2_gap <- R2_Y - Q2_Y
  dev_explained <- best_cv_fit$glmnet.fit$dev.ratio[which.min(abs(best_cv_fit$glmnet.fit$lambda - lambda_1se))]
  
  # Permutation test
  if (n_permutations > 0) {
    perm_results <- run_permutation_test(df_filtered, y, all_exclude, Q2_Y, n_permutations, alpha_grid, nfolds, seed + 200)
    perm_p_value <- perm_results$p_value; perm_n_completed <- perm_results$n_permutations_completed
    perm_null_mean <- mean(perm_results$null_distribution)
  } else { perm_p_value <- NA_real_; perm_n_completed <- 0; perm_null_mean <- NA_real_ }
  
  # Separate covariate vs penalized coefficients
  nonzero_coefs <- coef_vector[coef_vector != 0]
  is_cov_nz <- grepl(paste0("^(", paste(covariate_cols, collapse = "|"), ")"), names(nonzero_coefs))
  penalized_nonzero <- nonzero_coefs[!is_cov_nz]; covariate_nonzero <- nonzero_coefs[is_cov_nz]
  if (length(penalized_nonzero) > 0) {
    sorted_coefs <- sort(abs(penalized_nonzero), decreasing = TRUE)
    nonzero_predictor_names <- names(sorted_coefs)
    nonzero_predictor_coefs <- penalized_nonzero[nonzero_predictor_names]
  } else { nonzero_predictor_names <- NA; nonzero_predictor_coefs <- NA }
  
  tibble(
    model_name = "elastic_net_cov_adjusted", cvd_score = cvd_score_name,
    dataset_name = dataset_name, n_observations = n_obs, n_predictors = n_pred,
    n_penalized_predictors = n_penalized, n_covariates = n_covariates_in_model,
    covariates_used = paste(covariate_cols, collapse = ", "),
    n_nonzero_coefs = n_nonzero, n_nonzero_penalized = length(penalized_nonzero),
    n_nonzero_covariates = length(covariate_nonzero),
    alpha_selection_method = "min_rmse", alpha_optimal_fulldata = best_alpha_fulldata,
    alpha_mean_nested = mean(alpha_per_fold),
    lambda_min = lambda_min, lambda_1se = lambda_1se, lambda_selection_method = "one_se",
    n_folds = nfolds, cv_rmse_min = cv_rmse_min, cv_rmse_1se = cv_rmse_1se,
    R2_Y = R2_Y, R2_fold_mean = R2_fold_mean, R2_fold_sd = R2_fold_sd,
    R2_per_fold = list(R2_per_fold),
    in_sample_rmse = rmse_in, in_sample_mae = mae_in, in_sample_cor = cor_in,
    Q2_Y = Q2_Y, Q2_Y_rmse = rmse_oof, Q2_Y_mae = mae_oof,
    MAE_fold_mean = MAE_fold_mean, MAE_fold_sd = MAE_fold_sd,
    RMSE_fold_mean = RMSE_fold_mean, RMSE_fold_sd = RMSE_fold_sd,
    MAE_per_fold = list(MAE_per_fold), RMSE_per_fold = list(RMSE_per_fold),
    Q2_Y_cor = cor_oof, Q2_fold_mean = Q2_fold_mean, Q2_fold_sd = Q2_fold_sd,
    Q2_per_fold = list(Q2_per_fold), Q2_R2_ratio = Q2_R2_ratio, R2_Q2_gap = R2_Q2_gap,
    permutation_p_value = perm_p_value, permutation_n = perm_n_completed,
    permutation_null_mean_Q2 = perm_null_mean, pct_deviance_explained = dev_explained * 100,
    intercept = intercept, nonzero_predictors = list(nonzero_predictor_names),
    nonzero_coefficients = list(nonzero_predictor_coefs),
    covariate_coefficients = list(covariate_nonzero),
    all_predictor_names = list(predictor_names), all_coefficients = list(coef_vector),
    alpha_per_fold_nested = list(alpha_per_fold),
    n_features_per_fold = list(oof_results$n_features_per_fold),
    col_miss_threshold = col_miss_thresh, row_miss_threshold = row_miss_thresh,
    imputation_method = "kNN", seed = seed, date_run = as.character(Sys.time()),
    cv_fit_object = list(best_cv_fit), all_cv_results = list(cv_results)
  )
}

# ============================================
# 7. Model runner functions with incremental saving
# ============================================

run_score_specific_models <- function(score_specific_datasets, alpha_grid = seq(0, 1, by = 0.1),
                                      nfolds = 10, n_permutations = 100, covariate_cols = COVARIATE_COLS,
                                      seed = 42, output_file = "elastic_net_score_specific.qs2") {
  results_list <- list(); n_models <- length(score_specific_datasets)
  incremental_dir <- save_output("incremental_saves_specific")
  dir.create(incremental_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Elastic Net (cov-adjusted): Score-specific datasets\n")
  cat(sprintf("Models: %d | Perms: %d | Covariates: %s\n", n_models, n_permutations, paste(covariate_cols, collapse = ", ")))
  for (i in seq_along(score_specific_datasets)) {
    cvd_score <- names(score_specific_datasets)[i]; df <- score_specific_datasets[[i]]
    dataset_name <- paste0(cvd_score, "_specific")
    cat(sprintf("\n[%d/%d] %s ~ %s\n", i, n_models, cvd_score, dataset_name))
    incremental_file <- file.path(incremental_dir, sprintf("model_%02d_%s.rds", i, cvd_score))
    if (file.exists(incremental_file)) { cat("  Loading from incremental save\n"); results_list[[i]] <- readRDS(incremental_file); next }
    tryCatch({
      result <- run_elastic_net_model(df_raw = df, outcome_col = cvd_score, cvd_score_name = cvd_score,
                                      dataset_name = dataset_name, covariate_cols = covariate_cols,
                                      alpha_grid = alpha_grid, nfolds = nfolds, n_permutations = n_permutations, seed = seed)
      results_list[[i]] <- result; saveRDS(result, incremental_file); cat("  Success\n")
    }, error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); results_list[[i]] <- NULL })
  }
  results_df <- bind_rows(results_list)
  qs_save(results_df, file = save_output(output_file))
  saveRDS(results_df, file = save_output(gsub("\\.qs2$", ".rds", output_file)))
  cat(sprintf("\nSaved %d models to %s (+ .rds backup)\n", nrow(results_df), output_file)); results_df
}

run_common_datasets_models <- function(datasets_list, cvd_score_names, alpha_grid = seq(0, 1, by = 0.1),
                                       nfolds = 10, n_permutations = 100, covariate_cols = COVARIATE_COLS,
                                       seed = 42, output_file = "elastic_net_common.qs2") {
  results_list <- list(); counter <- 1; n_models <- length(datasets_list) * length(cvd_score_names)
  incremental_dir <- save_output("incremental_saves_common")
  dir.create(incremental_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Elastic Net (cov-adjusted): Common datasets x CVD scores\n")
  cat(sprintf("Models: %d | Perms: %d | Covariates: %s\n", n_models, n_permutations, paste(covariate_cols, collapse = ", ")))
  for (dataset_name in names(datasets_list)) {
    df <- datasets_list[[dataset_name]]
    for (cvd_score in cvd_score_names) {
      cat(sprintf("\n[%d/%d] %s ~ %s\n", counter, n_models, cvd_score, dataset_name))
      incremental_file <- file.path(incremental_dir, sprintf("model_%03d_%s_%s.rds", counter, cvd_score, dataset_name))
      if (file.exists(incremental_file)) { cat("  Loading from incremental save\n"); results_list[[counter]] <- readRDS(incremental_file); counter <- counter + 1; next }
      if (!cvd_score %in% colnames(df)) { cat(sprintf("  Column '%s' not found\n", cvd_score)); counter <- counter + 1; next }
      tryCatch({
        result <- run_elastic_net_model(df_raw = df, outcome_col = cvd_score, cvd_score_name = cvd_score,
                                        dataset_name = dataset_name, covariate_cols = covariate_cols,
                                        alpha_grid = alpha_grid, nfolds = nfolds, n_permutations = n_permutations, seed = seed)
        results_list[[counter]] <- result; saveRDS(result, incremental_file); cat("  Success\n")
      }, error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); results_list[[counter]] <- NULL })
      counter <- counter + 1
    }
  }
  results_df <- bind_rows(results_list)
  qs_save(results_df, file = save_output(output_file))
  saveRDS(results_df, file = save_output(gsub("\\.qs2$", ".rds", output_file)))
  cat(sprintf("\nSaved %d models to %s (+ .rds backup)\n", nrow(results_df), output_file)); results_df
}

run_score_specific_plus_blocks_models <- function(score_specific_datasets, datasets_list,
                                                  alpha_grid = seq(0, 1, by = 0.1),
                                                  nfolds = 10, n_permutations = 100,
                                                  covariate_cols = COVARIATE_COLS,
                                                  seed = 42, output_file = "elastic_net_score_plus_blocks.qs2") {
  results_list <- list(); counter <- 1; n_models <- length(score_specific_datasets) * length(datasets_list)
  incremental_dir <- save_output("incremental_saves_score_plus_blocks")
  dir.create(incremental_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Elastic Net (cov-adjusted): Score-specific + predictor block\n")
  cat(sprintf("Models: %d | Perms: %d | Covariates: %s\n", n_models, n_permutations, paste(covariate_cols, collapse = ", ")))
  for (cvd_score in names(score_specific_datasets)) {
    for (block_name in names(datasets_list)) {
      dataset_name <- paste0(cvd_score, "_specific+", block_name)
      cat(sprintf("\n[%d/%d] %s ~ %s\n", counter, n_models, cvd_score, dataset_name))
      incremental_file <- file.path(incremental_dir, sprintf("model_%03d_%s_%s.rds", counter, cvd_score, block_name))
      if (file.exists(incremental_file)) { cat("  Loading from incremental save\n"); results_list[[counter]] <- readRDS(incremental_file); counter <- counter + 1; next }
      df_combined <- build_score_plus_block(cvd_score, block_name, score_specific_datasets, datasets_list)
      tryCatch({
        result <- run_elastic_net_model(df_raw = df_combined, outcome_col = cvd_score, cvd_score_name = cvd_score,
                                        dataset_name = dataset_name, covariate_cols = covariate_cols,
                                        alpha_grid = alpha_grid, nfolds = nfolds, n_permutations = n_permutations, seed = seed)
        results_list[[counter]] <- result; saveRDS(result, incremental_file); cat("  Success\n")
      }, error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); results_list[[counter]] <- NULL })
      counter <- counter + 1
    }
  }
  results_df <- bind_rows(results_list)
  qs_save(results_df, file = save_output(output_file))
  saveRDS(results_df, file = save_output(gsub("\\.qs2$", ".rds", output_file)))
  cat(sprintf("\nSaved %d models to %s (+ .rds backup)\n", nrow(results_df), output_file)); results_df
}

# ============================================
# 8. Execute pipeline
# ============================================

results_score_specific <- run_score_specific_models(
  score_specific_datasets = score_specific_datasets,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 5, n_permutations = 0, seed = 42, # set to 1000 for thesis run
  output_file = "elastic_net_score_specific_5fold.qs2"
)

results_common <- run_common_datasets_models(
  datasets_list = datasets_list,
  cvd_score_names = CVD_scores,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 5, n_permutations = 0, seed = 42, # set to 1000 for thesis run
  covariate_cols = COVARIATE_COLS,
  output_file = "elastic_net_common_5fold.qs2"
)

results_score_plus_blocks <- run_score_specific_plus_blocks_models(
  score_specific_datasets = score_specific_datasets,
  datasets_list = datasets_list,
  alpha_grid = seq(0, 1, by = 0.1),
  nfolds = 5, n_permutations = 0, seed = 42, # set to 1000 for thesis run
  output_file = "elastic_net_score_plus_blocks_5fold.qs2"
)

# Combine all results
results_all <- bind_rows(
  results_score_specific, results_common, results_score_plus_blocks
)

# ============================================
# 9. Export to excel
# ============================================

results_main <- results_all %>%
  select(
    cvd_score, dataset_name, model_name,
    n_observations, n_predictors, n_penalized_predictors, n_covariates, covariates_used,
    n_nonzero_coefs, n_nonzero_penalized, n_nonzero_covariates,
    alpha_selection_method, alpha_optimal_fulldata, alpha_mean_nested,
    lambda_min, lambda_1se, lambda_selection_method,
    R2_Y, R2_fold_mean, R2_fold_sd, in_sample_rmse, in_sample_mae, in_sample_cor,
    Q2_Y, Q2_fold_mean, Q2_fold_sd,
    Q2_Y_rmse, Q2_Y_mae, Q2_Y_cor,
    Q2_R2_ratio, R2_Q2_gap,
    permutation_p_value, permutation_n, permutation_null_mean_Q2,
    cv_rmse_min, cv_rmse_1se, pct_deviance_explained,
    imputation_method, col_miss_threshold, row_miss_threshold,
    n_folds, seed, date_run
  )

predictor_rows <- list()
for (i in 1:nrow(results_all)) {
  preds <- results_all$nonzero_predictors[[i]]; coefs <- results_all$nonzero_coefficients[[i]]
  if (!is.null(preds) && length(preds) > 0 && !all(is.na(preds))) {
    for (j in seq_along(preds)) {
      if (!is.na(preds[j])) {
        predictor_rows[[length(predictor_rows) + 1]] <- tibble(
          cvd_score = results_all$cvd_score[i], dataset_name = results_all$dataset_name[i],
          rank = j, predictor_name = preds[j], coefficient = round(coefs[j], 4), type = "penalized")
      }
    }
  }
  cov_coefs <- results_all$covariate_coefficients[[i]]
  if (!is.null(cov_coefs) && length(cov_coefs) > 0) {
    for (j in seq_along(cov_coefs)) {
      predictor_rows[[length(predictor_rows) + 1]] <- tibble(
        cvd_score = results_all$cvd_score[i], dataset_name = results_all$dataset_name[i],
        rank = NA_integer_, predictor_name = names(cov_coefs)[j],
        coefficient = round(cov_coefs[j], 4), type = "covariate (unpenalized)")
    }
  }
}
top_predictors_sheet <- bind_rows(predictor_rows)

table_A_all_models <- results_all %>%
  select(dataset_name, cvd_score, Q2_Y) %>% mutate(Q2_Y = round(Q2_Y, 3)) %>%
  distinct(dataset_name, cvd_score, .keep_all = TRUE) %>%
  pivot_wider(names_from = cvd_score, values_from = Q2_Y) %>%
  mutate(dataset_group = case_when(
    grepl("_specific\\+", dataset_name) ~ "specific+block",
    grepl("_specific$", dataset_name) ~ "specific", TRUE ~ "block_only")) %>%
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

headerStyle <- createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
wb <- createWorkbook()
addWorksheet(wb, "All_Results"); writeData(wb, "All_Results", results_main)
addStyle(wb, "All_Results", headerStyle, rows = 1, cols = 1:ncol(results_main), gridExpand = TRUE)
freezePane(wb, "All_Results", firstRow = TRUE)
setColWidths(wb, "All_Results", cols = 1:ncol(results_main), widths = "auto")
addWorksheet(wb, "Nonzero_Predictors"); writeData(wb, "Nonzero_Predictors", top_predictors_sheet)
addStyle(wb, "Nonzero_Predictors", headerStyle, rows = 1, cols = 1:ncol(top_predictors_sheet), gridExpand = TRUE)
freezePane(wb, "Nonzero_Predictors", firstRow = TRUE)
setColWidths(wb, "Nonzero_Predictors", cols = 1:ncol(top_predictors_sheet), widths = "auto")
addWorksheet(wb, "Q2_Y_Overview")
writeData(wb, "Q2_Y_Overview", "Table A: Q²_Y (all models) — Covariate-Adjusted", startRow = 1, startCol = 1)
writeData(wb, "Q2_Y_Overview", table_A_all_models, startRow = 2, startCol = 1)
start_row_B <- nrow(table_A_all_models) + 5
writeData(wb, "Q2_Y_Overview", "Table B: Δ Q²_Y = (score-specific + block) − (score-specific baseline)", startRow = start_row_B, startCol = 1)
writeData(wb, "Q2_Y_Overview", table_B_delta, startRow = start_row_B + 1, startCol = 1)
addStyle(wb, "Q2_Y_Overview", headerStyle, rows = c(2, start_row_B + 1),
         cols = 1:max(ncol(table_A_all_models), ncol(table_B_delta)), gridExpand = TRUE)
saveWorkbook(wb, save_output("elastic_net_results_cov_adjusted.xlsx"), overwrite = TRUE)
cat("\n✓ Excel saved: elastic_net_results_cov_adjusted.xlsx\n")

# ============================================
# 10. Creates Tables
# ============================================

cvd_score_order <- c("ascvd_10y", "frs_10y", "QRISK3_risk", "SCORE2_score")
cvd_score_labels <- c("SCORE2_score" = "SCORE2", "QRISK3_risk" = "QRISK3",
                      "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD")
dataset_labels_table <- c(
  "all_data" = "Full predictor set", "body_composition" = "Body composition",
  "sociodemographics_lifestyle" = "Sociodemographics & lifestyle",
  "clinical_risk_factors" = "Clinical Measurements & Supplementary Biomarkers",
  "lipids" = "Lipids", "fatty_acids" = "Fatty acids", "urine_nmr" = "Urine NMR")

create_score_specific_table <- function(results_df, table_title) {
  table_data <- results_df %>%
    filter(cvd_score %in% cvd_score_order) %>%
    mutate(cvd_score = factor(cvd_score, levels = cvd_score_order),
           CVD_Score = cvd_score_labels[as.character(cvd_score)]) %>%
    arrange(cvd_score)
  formatted_table <- table_data %>%
    mutate(
      R2 = ifelse(!is.na(R2_fold_sd) & R2_fold_sd > 0, sprintf("%.3f \u00B1 %.3f", R2_fold_mean, R2_fold_sd), sprintf("%.3f", R2_Y)),
      Q2 = ifelse(!is.na(Q2_fold_sd) & Q2_fold_sd > 0, sprintf("%.3f \u00B1 %.3f", Q2_Y, Q2_fold_sd), sprintf("%.3f", Q2_Y)),
      MAE = ifelse(!is.na(MAE_fold_sd) & MAE_fold_sd > 0, sprintf("%.3f \u00B1 %.3f", Q2_Y_mae, MAE_fold_sd), sprintf("%.3f", Q2_Y_mae)),
      perm_p = ifelse(is.na(permutation_p_value), "\u2014", ifelse(permutation_p_value < 0.001, "<0.001", sprintf("%.3f", permutation_p_value)))
    ) %>% select(CVD_Score, R2, Q2, MAE, perm_p)
  formatted_table %>% flextable() %>%
    set_caption(caption = table_title) %>%
    set_header_labels(CVD_Score = "CVD Score", R2 = "R\u00B2 \u00B1 SD", Q2 = "Q\u00B2 \u00B1 SD", MAE = "MAE \u00B1 SD", perm_p = "Perm. p") %>%
    align(j = 1, align = "left", part = "all") %>% align(j = 2:5, align = "center", part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 9, part = "body") %>% flextable::fontsize(size = 10, part = "header") %>%
    bold(part = "header") %>% width(j = 1, width = 1.2) %>% width(j = 2:4, width = 1.3) %>% width(j = 5, width = 1.0) %>%
    border_remove() %>% hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1.5), part = "header") %>% padding(padding = 3, part = "all")
}

create_multi_score_table <- function(results_df, table_title, dataset_order = NULL, dataset_labels = dataset_labels_table) {
  table_data <- results_df %>%
    filter(cvd_score %in% cvd_score_order) %>%
    mutate(cvd_score = factor(cvd_score, levels = cvd_score_order),
           CVD_Score = cvd_score_labels[as.character(cvd_score)])
  if (!is.null(dataset_order)) {
    table_data <- table_data %>% filter(dataset_name %in% dataset_order) %>%
      mutate(dataset_name = factor(dataset_name, levels = dataset_order))
  }
  if (nrow(table_data) == 0) { warning("No data!"); return(NULL) }
  table_data <- table_data %>% arrange(dataset_name, cvd_score)
  formatted_table <- table_data %>%
    mutate(
      Dataset = ifelse(dataset_name %in% names(dataset_labels), dataset_labels[as.character(dataset_name)], as.character(dataset_name)),
      R2 = ifelse(!is.na(R2_fold_sd) & R2_fold_sd > 0, sprintf("%.3f \u00B1 %.3f", R2_fold_mean, R2_fold_sd), sprintf("%.3f", R2_Y)),
      Q2 = ifelse(!is.na(Q2_fold_sd) & Q2_fold_sd > 0, sprintf("%.3f \u00B1 %.3f", Q2_Y, Q2_fold_sd), sprintf("%.3f", Q2_Y)),
      MAE = ifelse(!is.na(MAE_fold_sd) & MAE_fold_sd > 0, sprintf("%.3f \u00B1 %.3f", Q2_Y_mae, MAE_fold_sd), sprintf("%.3f", Q2_Y_mae)),
      perm_p = ifelse(is.na(permutation_p_value), "\u2014", ifelse(permutation_p_value < 0.001, "<0.001", sprintf("%.3f", permutation_p_value)))
    ) %>% select(Dataset, CVD_Score, R2, Q2, MAE, perm_p)
  n_rows <- nrow(formatted_table)
  ft <- formatted_table %>% flextable() %>% merge_v(j = "Dataset") %>%
    set_caption(caption = table_title) %>%
    set_header_labels(CVD_Score = "CVD Score", R2 = "R\u00B2 \u00B1 SD", Q2 = "Q\u00B2 \u00B1 SD", MAE = "MAE \u00B1 SD", perm_p = "Perm. p") %>%
    align(j = 1:2, align = "left", part = "all") %>% align(j = 3:6, align = "center", part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 9, part = "body") %>% flextable::fontsize(size = 10, part = "header") %>%
    bold(part = "header") %>% bold(j = 1, part = "body") %>%
    width(j = 1, width = 1.8) %>% width(j = 2, width = 1.0) %>% width(j = 3:5, width = 1.3) %>% width(j = 6, width = 1.0) %>%
    border_remove() %>% hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1.5), part = "header") %>% padding(padding = 3, part = "all")
  if (n_rows > 4) { hline_rows <- seq(4, n_rows - 1, by = 4)
  if (length(hline_rows) > 0) ft <- ft %>% hline(i = hline_rows, border = fp_border(width = 0.5, color = "#CCCCCC"), part = "body") }
  ft
}

# Create and save tables
table1 <- create_score_specific_table(results_score_specific, "Table 1: Elastic Net – Score-Specific (Covariate-Adjusted)")
datasets_table2 <- c("all_data", "lipids", "fatty_acids", "urine_nmr", "body_composition", "sociodemographics_lifestyle", "clinical_risk_factors")
datasets_table2 <- datasets_table2[datasets_table2 %in% unique(results_common$dataset_name)]
table2 <- create_multi_score_table(results_common %>% filter(dataset_name %in% datasets_table2),
                                   "Table 2: Elastic Net – Predictor Blocks (Covariate-Adjusted)", dataset_order = datasets_table2)
blocks_table3 <- datasets_table2
dataset_order_table3 <- results_score_plus_blocks %>%
  mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>%
  filter(block %in% blocks_table3) %>% arrange(match(block, blocks_table3)) %>% pull(dataset_name) %>% unique()
dataset_labels_table3 <- dataset_labels_table
for (ds in dataset_order_table3) {
  block <- sub("^.*_specific\\+", "", ds)
  block_label <- ifelse(block %in% names(dataset_labels_table), dataset_labels_table[block], block)
  dataset_labels_table3[ds] <- paste0("Score inputs + ", block_label)
}
table3 <- create_multi_score_table(
  results_score_plus_blocks %>% mutate(block = sub("^.*_specific\\+", "", dataset_name)) %>% filter(block %in% blocks_table3),
  "Table 3: Elastic Net – Score-Specific + Blocks (Covariate-Adjusted)",
  dataset_order = dataset_order_table3, dataset_labels = dataset_labels_table3)

save_as_docx(table1, path = save_output("Table1_score_specific.docx"))
save_as_docx(table2, path = save_output("Table2_predictor_blocks.docx"))
save_as_docx(table3, path = save_output("Table3_score_plus_blocks.docx"))
save_as_html(table1, path = save_output("Table1_score_specific.html"))
save_as_html(table2, path = save_output("Table2_predictor_blocks.html"))
save_as_html(table3, path = save_output("Table3_score_plus_blocks.html"))

doc <- read_docx() %>%
  body_add_par("Elastic Net Results — Covariate-Adjusted", style = "heading 1") %>%
  body_add_par(sprintf("Covariates (unpenalized): %s", paste(COVARIATE_COLS, collapse = ", "))) %>%
  body_add_par("") %>% body_add_flextable(value = table1) %>% body_add_break() %>%
  body_add_flextable(value = table2) %>% body_add_break() %>% body_add_flextable(value = table3)
print(doc, target = save_output("Elastic_Net_All_Tables_CovAdj.docx"))

# ============================================
# 11. Plot for model performance comparison
# ============================================

block_order <- c("lipids", "fatty_acids", "urine_nmr", "body_composition", "clinical_risk_factors", "sociodemographics_lifestyle")
block_order <- block_order[block_order %in% unique(results_common$dataset_name)]
model_type_order <- c("score_specific", "score_specific+all_data", paste0("score_specific+", block_order))
model_type_labels <- c(
  "score_specific" = "Traditional inputs only", "score_specific+all_data" = "Traditional inputs + Full set",
  "score_specific+lipids" = "Traditional inputs + Lipids", "score_specific+fatty_acids" = "Traditional inputs + Fatty acids",
  "score_specific+urine_nmr" = "Traditional inputs + Urine NMR", "score_specific+body_composition" = "Traditional inputs + Body composition",
  "score_specific+clinical_risk_factors" = "Traditional inputs + Clinical Measurements & Suppl. Biomarkers",
  "score_specific+sociodemographics_lifestyle" = "Traditional inputs + Sociodemographics & Lifestyle")
predictor_colors <- c(
  "Traditional inputs only" = "#888888", "Traditional inputs + Full set" = "#2D2D2D",
  "Traditional inputs + Lipids" = "#E63946", "Traditional inputs + Fatty acids" = "#F4A261",
  "Traditional inputs + Urine NMR" = "#E9C46A", "Traditional inputs + Body composition" = "#457B9D",
  "Traditional inputs + Clinical Measurements & Suppl. Biomarkers" = "#2A9D8F",
  "Traditional inputs + Sociodemographics & Lifestyle" = "#7B2D8E")

d_score_specific <- results_score_specific %>% filter(cvd_score %in% cvd_score_order) %>% mutate(model_type = "score_specific")
d_score_plus <- results_score_plus_blocks %>% filter(cvd_score %in% cvd_score_order) %>%
  mutate(block = sub("^.*_specific\\+", "", dataset_name), model_type = paste0("score_specific+", block))
plot_data_grouped <- bind_rows(d_score_specific, d_score_plus) %>%
  filter(model_type %in% model_type_order) %>%
  mutate(cvd_label = cvd_score_labels[as.character(cvd_score)],
         cvd_label = factor(cvd_label, levels = rev(cvd_score_labels)),
         model_type = factor(model_type, levels = model_type_order),
         model_label = model_type_labels[as.character(model_type)],
         model_label = factor(model_label, levels = model_type_labels[model_type_order]))
if (!"Q2_fold_sd" %in% names(plot_data_grouped)) plot_data_grouped$Q2_fold_sd <- 0
if (!"MAE_fold_sd" %in% names(plot_data_grouped)) plot_data_grouped$MAE_fold_sd <- 0
q2_ymin <- floor(min(plot_data_grouped$Q2_Y - plot_data_grouped$Q2_fold_sd, na.rm = TRUE) * 10) / 10

p_q2_vertical <- plot_data_grouped %>%
  ggplot(aes(x = cvd_label, y = Q2_Y, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_errorbar(aes(ymin = Q2_Y - Q2_fold_sd, ymax = Q2_Y + Q2_fold_sd),
                position = position_dodge(width = 0.85), width = 0.3, linewidth = 0.3) +
  scale_fill_manual(values = predictor_colors, breaks = unname(model_type_labels[model_type_order])) +
  coord_cartesian(ylim = c(q2_ymin, 1.02)) + scale_y_continuous(breaks = seq(q2_ymin, 1, by = 0.1)) +
  labs(title = expression(bold(paste("A) ", Q^2, " ± SD"))), y = expression(bold(Q^2)), x = NULL, fill = "Model") +
  theme_minimal() + theme(text = element_text(family = "Arial"), legend.position = "bottom",
                          legend.text = element_text(size = 11), legend.title = element_text(size = 13, face = "bold"),
                          plot.title = element_text(hjust = 0, size = 15), axis.text.x = element_text(size = 12, face = "bold"),
                          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14, face = "bold"),
                          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), plot.title.position = "plot") +
  guides(fill = guide_legend(nrow = 2, override.aes = list(size = 3)))

p_mae_vertical <- plot_data_grouped %>%
  ggplot(aes(x = cvd_label, y = Q2_Y_mae, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_errorbar(aes(ymin = pmax(0, Q2_Y_mae - MAE_fold_sd), ymax = Q2_Y_mae + MAE_fold_sd),
                position = position_dodge(width = 0.85), width = 0.3, linewidth = 0.3) +
  scale_fill_manual(values = predictor_colors, breaks = unname(model_type_labels[model_type_order])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = expression(bold(paste("B) MAE ± SD"))), y = expression(bold("MAE")), x = NULL, fill = "Model") +
  theme_minimal() + theme(text = element_text(family = "Arial"), legend.position = "bottom",
                          legend.text = element_text(size = 11), legend.title = element_text(size = 13, face = "bold"),
                          plot.title = element_text(hjust = 0, size = 15), axis.text.x = element_text(size = 12, face = "bold"),
                          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14, face = "bold"),
                          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), plot.title.position = "plot") +
  guides(fill = guide_legend(nrow = 2, override.aes = list(size = 3)))

combined_vertical_plot <- (p_q2_vertical | p_mae_vertical) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.justification = "center")

ggsave(save_output("elastic_net_grouped_by_score_vertical.png"), plot = combined_vertical_plot, width = 14, height = 7, dpi = 300, bg = "white")
ggsave(save_output("elastic_net_grouped_by_score_vertical.pdf"), plot = combined_vertical_plot, width = 14, height = 7, device = "pdf")
cat("\n✓ Vertical grouped-by-score plots saved\n")

# ============================================
# 12. Tuning diagnostics plots
# ============================================

extract_tuning_from_results <- function(results_row, dataset_display_name = NULL) {
  cv_results <- results_row$all_cv_results[[1]]; best_cv_fit <- results_row$cv_fit_object[[1]]
  best_alpha <- results_row$alpha_optimal_fulldata; n_pred <- results_row$n_predictors
  n_obs <- results_row$n_observations; nfolds <- results_row$n_folds
  if (is.null(dataset_display_name)) dataset_display_name <- results_row$dataset_name
  alpha_grid <- sapply(cv_results, function(x) x$alpha)
  alpha_rmse_df <- tibble(
    alpha = alpha_grid,
    cv_rmse = sapply(cv_results, function(x) { idx <- which(x$cv_fit$lambda == x$cv_fit$lambda.1se); sqrt(x$cv_fit$cvm[idx]) }),
    cv_rmse_se = sapply(cv_results, function(x) { idx <- which(x$cv_fit$lambda == x$cv_fit$lambda.1se); mse_sd <- x$cv_fit$cvsd[idx]; rmse_val <- sqrt(x$cv_fit$cvm[idx]); (mse_sd / sqrt(nfolds)) / (2 * rmse_val) }))
  best_alpha_idx <- which.min(alpha_rmse_df$cv_rmse)
  lambda_df <- tibble(
    log_lambda = log(best_cv_fit$lambda), cv_rmse = sqrt(best_cv_fit$cvm),
    cv_rmse_se = sapply(seq_along(best_cv_fit$lambda), function(i) { mse_sd <- best_cv_fit$cvsd[i]; rmse_val <- sqrt(best_cv_fit$cvm[i]); (mse_sd / sqrt(nfolds)) / (2 * rmse_val) }),
    nzero = best_cv_fit$nzero)
  nzero_at_min <- best_cv_fit$nzero[which(best_cv_fit$lambda == best_cv_fit$lambda.min)]
  nzero_at_1se <- best_cv_fit$nzero[which(best_cv_fit$lambda == best_cv_fit$lambda.1se)]
  list(alpha_rmse_df = alpha_rmse_df, best_alpha_idx = best_alpha_idx, best_alpha = best_alpha,
       best_cv_fit = best_cv_fit, lambda_df = lambda_df, nzero_at_min = nzero_at_min,
       nzero_at_1se = nzero_at_1se, n_raw_predictors = n_pred, n_obs = n_obs, dataset_name = dataset_display_name)
}

build_alpha_panel <- function(res, show_y_label = TRUE) {
  df <- res$alpha_rmse_df; selected <- df[res$best_alpha_idx, ]
  ggplot(df, aes(x = alpha, y = cv_rmse)) +
    geom_line(colour = "grey60", linewidth = 0.5) + geom_point(colour = "grey40", size = 2, alpha = 0.7) +
    geom_errorbar(aes(ymin = cv_rmse - cv_rmse_se, ymax = cv_rmse + cv_rmse_se), colour = "grey70", alpha = 0.5, width = 0) +
    geom_point(data = selected, colour = "#E41A1C", size = 3.5, shape = 18) +
    geom_vline(xintercept = selected$alpha, colour = "#E41A1C", linetype = "dashed", alpha = 0.5) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    labs(title = res$dataset_name, subtitle = sprintf("p = %d, n = %d", res$n_raw_predictors, res$n_obs),
         x = expression(alpha), y = if (show_y_label) "CV RMSE" else NULL) +
    theme_minimal(base_size = 11) + theme(text = element_text(family = "Arial"),
                                          plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
                                          plot.subtitle = element_text(size = 13, colour = "grey40", hjust = 0.5),
                                          axis.title = element_text(size = 12), axis.text = element_text(size = 11), panel.grid.minor = element_blank())
}

build_lambda_panel <- function(res, show_y_label = TRUE) {
  df <- res$lambda_df; cv_fit <- res$best_cv_fit
  lambda_min_val <- log(cv_fit$lambda.min); lambda_1se_val <- log(cv_fit$lambda.1se)
  sel_idx <- which.min(abs(df$log_lambda - lambda_1se_val)); selected <- df[sel_idx, ]
  n_breaks <- min(6, nrow(df)); break_idx <- seq(1, nrow(df), length.out = n_breaks) %>% round()
  ggplot(df, aes(x = log_lambda, y = cv_rmse)) +
    geom_line(colour = "grey60", linewidth = 0.5) +
    geom_errorbar(aes(ymin = cv_rmse - cv_rmse_se, ymax = cv_rmse + cv_rmse_se), colour = "grey70", alpha = 0.3, width = 0) +
    geom_vline(xintercept = lambda_min_val, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = lambda_1se_val, linetype = "dashed", colour = "#E41A1C", alpha = 0.5, linewidth = 0.4) +
    geom_point(data = selected, colour = "#E41A1C", size = 3.5, shape = 18) +
    annotate("text", x = lambda_min_val, y = max(df$cv_rmse) * 0.98, label = sprintf("p==%d", res$nzero_at_min), parse = TRUE, hjust = 1.1, size = 3.5, colour = "grey40") +
    annotate("text", x = lambda_1se_val, y = max(df$cv_rmse) * 0.98, label = sprintf("p==%d", res$nzero_at_1se), parse = TRUE, hjust = -0.1, size = 3.5, colour = "#E41A1C") +
    scale_x_continuous(sec.axis = sec_axis(transform = ~ ., breaks = df$log_lambda[break_idx], labels = df$nzero[break_idx])) +
    labs(x = expression(log(lambda)), y = if (show_y_label) "CV RMSE" else NULL,
         subtitle = sprintf("\u03b1 = %.1f  |  \u03bb.1se: %d coefs", res$best_alpha, res$nzero_at_1se)) +
    theme_minimal(base_size = 11) + theme(text = element_text(family = "Arial"),
                                          plot.subtitle = element_text(size = 13, colour = "grey40", hjust = 0.5),
                                          axis.title = element_text(size = 12), axis.text = element_text(size = 11),
                                          axis.title.x.top = element_blank(), axis.text.x.top = element_text(size = 9), panel.grid.minor = element_blank())
}

assemble_tuning_figure <- function(tuning_results, add_title = TRUE) {
  n <- length(tuning_results)
  alpha_panels <- lapply(seq_len(n), function(i) build_alpha_panel(tuning_results[[i]], show_y_label = (i == 1)))
  lambda_panels <- lapply(seq_len(n), function(i) build_lambda_panel(tuning_results[[i]], show_y_label = (i == 1)))
  alpha_panels[[1]] <- alpha_panels[[1]] + labs(tag = "A") + theme(plot.tag = element_text(face = "bold", family = "Arial", size = 14))
  lambda_panels[[1]] <- lambda_panels[[1]] + labs(tag = "B") + theme(plot.tag = element_text(face = "bold", family = "Arial", size = 14))
  combined <- wrap_plots(alpha_panels, nrow = 1) / wrap_plots(lambda_panels, nrow = 1)
  if (add_title) combined <- combined + plot_annotation(
    title = "Elastic Net Hyperparameter Tuning",
    theme = theme(plot.title = element_text(face = "bold", family = "Arial", size = 14),
                  plot.subtitle = element_text(family = "Arial", size = 10, colour = "grey40")))
  combined
}


# create plots for supplementary
for (score in c("QRISK3_risk", "SCORE2_score", "frs_10y", "ascvd_10y")) {
  
  score_label <- cvd_score_labels[score]
  
  ss <- results_score_specific %>% filter(cvd_score == score)
  cm <- results_common %>% filter(cvd_score == score)
  
  tuning_configs <- list(
    list(row = ss, name = "Traditional Inputs"),
    list(row = cm %>% filter(dataset_name == "all_data"), name = "Multi-Domain Predictor Set"),
    list(row = cm %>% filter(dataset_name == "body_composition"), name = "Body Composition"),
    list(row = cm %>% filter(dataset_name == "sociodemographics_lifestyle"), name = "Sociodemographics"),
    list(row = cm %>% filter(dataset_name == "clinical_risk_factors"), name = "Clinical Risk Factors"),
    list(row = cm %>% filter(dataset_name == "lipids"), name = "Lipids"),
    list(row = cm %>% filter(dataset_name == "fatty_acids"), name = "Fatty Acids"),
    list(row = cm %>% filter(dataset_name == "urine_nmr"), name = "Urine NMR"))
  
  tuning_results <- lapply(tuning_configs, function(cfg) {
    extract_tuning_from_results(cfg$row, dataset_display_name = cfg$name)
  })
  
  p <- assemble_tuning_figure(tuning_results, add_title = TRUE)
  p <- p + plot_annotation(
    title = paste("Elastic Net Hyperparameter Tuning:", score_label),
    theme = theme(plot.title = element_text(face = "bold", family = "Arial", size = 14),
                  plot.subtitle = element_text(family = "Arial", size = 10, colour = "grey40")))
  
  ggsave(save_output(paste0("tuning_diagnostics_", score, "_all.png")),
         plot = p, width = 24, height = 9, dpi = 300, bg = "white")
  ggsave(save_output(paste0("tuning_diagnostics_", score, "_all.pdf")),
         plot = p, width = 24, height = 9, dpi = 300, device = cairo_pdf)
  
  cat("Saved tuning plot for", score_label, "\n")
}


### tuning plot table for supplementary
tuning_table <- results_all %>%
  mutate(
    score_label = cvd_score_labels[as.character(cvd_score)],
    dataset_label = case_when(
      grepl("_specific$", dataset_name) ~ "Traditional Inputs",
      dataset_name == "all_data" ~ "Multi-Domain Predictor Set",
      dataset_name == "body_composition" ~ "Body Composition",
      dataset_name == "sociodemographics_lifestyle" ~ "Sociodemographics & Lifestyle",
      dataset_name == "clinical_risk_factors" ~ "Clinical Measurements & Suppl. Biomarkers",
      dataset_name == "lipids" ~ "Lipids",
      dataset_name == "fatty_acids" ~ "Fatty Acids",
      dataset_name == "urine_nmr" ~ "Urine NMR",
      TRUE ~ dataset_name
    )
  ) %>%
  transmute(
    Dataset = dataset_label,
    `CVD Score` = score_label,
    p = n_predictors,
    `α` = sprintf("%.1f", alpha_optimal_fulldata),
    `λ.1se` = sprintf("%.4f", lambda_1se),
    `Coefs retained` = paste0(n_nonzero_penalized, " / ", n_penalized_predictors)
  ) %>%
  arrange(factor(Dataset, levels = c(
    "Multi-Domain Predictor Set", "Lipids", "Fatty Acids", "Urine NMR",
    "Body Composition", "Sociodemographics & Lifestyle",
    "Clinical Measurements & Suppl. Biomarkers", "Traditional Inputs")),
    factor(`CVD Score`, levels = c("ASCVD", "Framingham", "QRISK3", "SCORE2")))

ft <- flextable(tuning_table) %>%
  merge_v(j = "Dataset") %>%
  style(part = "body", pr_t = fp_text(font.size = 9, font.family = "Arial")) %>%
  style(part = "header", pr_t = fp_text(font.size = 9, font.family = "Arial", bold = TRUE)) %>%
  bold(j = "Dataset", part = "body") %>%
  align(j = 3:6, align = "center", part = "all") %>%
  align(j = 1:2, align = "left", part = "all") %>%
  width(j = 1, width = 2.5) %>%
  width(j = 2, width = 1.2) %>%
  width(j = 3, width = 0.5) %>%
  width(j = 4, width = 0.5) %>%
  width(j = 5, width = 0.8) %>%
  width(j = 6, width = 1.2) %>%
  border_remove() %>%
  hline_top(border = fp_border(width = 2), part = "all") %>%
  hline_bottom(border = fp_border(width = 2), part = "body") %>%
  hline(i = 1, border = fp_border(width = 1), part = "header")

doc <- read_docx() %>% body_add_flextable(value = ft)
print(doc, target = save_output("Supplementary_Table_Tuning_Parameters.docx"))


# ============================================
# 13. Permutation feature importance
# ============================================

calculate_single_feature_importance <- function(feat, X_test, y_test,
                                                y_pred_baseline, mae_baseline,
                                                cv_fit, n_permutations, seed) {
  set.seed(seed)
  feat_idx <- which(colnames(X_test) == feat)
  if (length(feat_idx) == 0) return(tibble(feature = feat, importance = NA_real_, importance_sd = NA_real_, p_value = NA_real_))
  mae_permuted <- numeric(n_permutations)
  for (p in 1:n_permutations) {
    X_test_perm <- X_test; X_test_perm[, feat_idx] <- sample(X_test_perm[, feat_idx])
    y_pred_perm <- as.numeric(predict(cv_fit, newx = X_test_perm, s = "lambda.1se"))
    mae_permuted[p] <- mean(abs(y_test - y_pred_perm))
  }
  delta_mae <- mae_permuted - mae_baseline
  tibble(feature = feat, importance = mean(delta_mae), importance_sd = sd(delta_mae), p_value = mean(delta_mae <= 0))
}

calculate_fold_importance <- function(X_train, X_test, y_train, y_test,
                                      cv_fit, features_to_test,
                                      n_permutations = 1000, n_cores = 1, seed = 42) {
  set.seed(seed)
  y_pred_baseline <- as.numeric(predict(cv_fit, newx = X_test, s = "lambda.1se"))
  mae_baseline <- mean(abs(y_test - y_pred_baseline))
  features_to_test <- intersect(features_to_test, colnames(X_test))
  if (length(features_to_test) == 0) return(tibble(feature = character(), importance = numeric(), importance_sd = numeric(), p_value = numeric()))
  feature_seeds <- seed + seq_along(features_to_test) * 1000
  importance_results <- parallel::mclapply(seq_along(features_to_test), function(i) {
    calculate_single_feature_importance(features_to_test[i], X_test, y_test, y_pred_baseline, mae_baseline, cv_fit, n_permutations, feature_seeds[i])
  }, mc.cores = n_cores)
  bind_rows(importance_results)
}

get_nonzero_features <- function(results_row, covariate_cols = COVARIATE_COLS) {
  coefs <- results_row$all_coefficients[[1]]; feat_names <- results_row$all_predictor_names[[1]]
  if (is.null(coefs) || is.null(feat_names)) return(character())
  nonzero_mask <- coefs != 0
  is_covariate <- grepl(paste0("^(", paste(covariate_cols, collapse = "|"), ")"), feat_names)
  feat_names[nonzero_mask & !is_covariate]
}

build_model_dataset <- function(cvd_score, dataset_name, score_specific_datasets, datasets_list) {
  if (grepl("_specific$", dataset_name)) return(score_specific_datasets[[cvd_score]])
  if (dataset_name %in% names(datasets_list)) return(datasets_list[[dataset_name]])
  block <- sub("^.*_specific\\+", "", dataset_name)
  if (!(block %in% names(datasets_list))) { warning(sprintf("Block '%s' not found", block)); return(NULL) }
  build_score_plus_block(cvd_score, block, score_specific_datasets, datasets_list)
}

run_permutation_importance <- function(df_raw, outcome_col, exclude_cols, features_to_test,
                                       covariate_cols = COVARIATE_COLS, alpha_grid = seq(0, 1, by = 0.1),
                                       nfolds = 5, n_permutations = 1000, n_cores = 1,
                                       col_miss_thresh = 0.4, row_miss_thresh = 0.8, seed = 42) {
  set.seed(seed)
  df_complete <- df_raw %>% filter(!is.na(.data[[outcome_col]]))
  y <- df_complete[[outcome_col]]
  all_exclude <- unique(c(exclude_cols, outcome_col, CVD_scores))
  cols_to_protect <- unique(c(all_exclude, covariate_cols))
  df_filtered <- suppressMessages(filter_by_missingness(df_complete, cols_to_protect, col_miss_thresh, row_miss_thresh))
  for (cov in covariate_cols) {
    if (!cov %in% names(df_filtered) && cov %in% names(df_complete))
      df_filtered[[cov]] <- df_complete[[cov]][match(df_filtered$Sample_ID, df_complete$Sample_ID)]
  }
  y <- df_filtered[[outcome_col]]; n <- length(y)
  if (n < nfolds * 2) warning(sprintf("Very few observations (%d) for %d-fold CV", n, nfolds))
  foldid <- create_stratified_folds(y, nfolds, seed = seed); fold_results <- list()
  for (k in 1:nfolds) {
    cat(sprintf("    Fold %d/%d: ", k, nfolds))
    idx_test <- which(foldid == k); idx_train <- setdiff(1:n, idx_test)
    df_train <- df_filtered[idx_train, , drop = FALSE]; df_test <- df_filtered[idx_test, , drop = FALSE]
    y_train <- y[idx_train]; y_test <- y[idx_test]
    preprocessed <- tryCatch(preprocess_train_test(df_train, df_test, all_exclude, covariate_cols),
                             error = function(e) { cat(sprintf("Preprocessing failed: %s\n", e$message)); NULL })
    if (is.null(preprocessed) || ncol(preprocessed$X_train) == 0) { cat("No features after preprocessing\n"); next }
    X_train <- preprocessed$X_train; X_test <- preprocessed$X_test
    pen_factor <- ifelse(preprocessed$is_covariate, 0, 1)
    cv_results_imp <- list(); valid_count <- 0
    for (a_idx in seq_along(alpha_grid)) {
      cv_fit <- tryCatch(cv.glmnet(x = X_train, y = y_train, alpha = alpha_grid[a_idx],
                                   penalty.factor = pen_factor, nfolds = nfolds,
                                   type.measure = "mse", standardize = FALSE, family = "gaussian"), error = function(e) NULL)
      if (!is.null(cv_fit)) { valid_count <- valid_count + 1; cv_results_imp[[valid_count]] <- list(alpha = alpha_grid[a_idx], cv_fit = cv_fit) }
    }
    if (valid_count == 0) { cat("Model fitting failed\n"); next }
    valid_alphas <- sapply(cv_results_imp, function(x) x$alpha)
    alpha_sel <- select_best_alpha(cv_results_imp, valid_alphas)
    best_cv_fit <- cv_results_imp[[alpha_sel$best_idx]]$cv_fit
    available_features <- colnames(X_test)
    features_this_fold <- features_to_test[features_to_test %in% available_features]
    for (f in features_to_test) {
      if (!(f %in% available_features)) features_this_fold <- c(features_this_fold, grep(paste0("^", f, "_"), available_features, value = TRUE))
    }
    features_this_fold <- unique(features_this_fold); n_features <- length(features_this_fold)
    cat(sprintf("%d features (parallel on %d cores)... ", n_features, n_cores))
    if (n_features == 0) { cat("No matching features\n"); next }
    fold_start <- Sys.time()
    fold_importance <- calculate_fold_importance(X_train, X_test, y_train, y_test, best_cv_fit, features_this_fold, n_permutations, n_cores, seed + k * 10000)
    fold_elapsed <- difftime(Sys.time(), fold_start, units = "secs")
    fold_importance$fold <- k; fold_results[[k]] <- fold_importance
    cat(sprintf("done in %.1fs\n", as.numeric(fold_elapsed)))
  }
  all_folds <- bind_rows(fold_results)
  if (nrow(all_folds) == 0) { warning("No results from any fold"); return(NULL) }
  aggregated <- all_folds %>% group_by(feature) %>%
    summarise(importance = mean(importance, na.rm = TRUE),
              importance_se = sd(importance, na.rm = TRUE) / sqrt(sum(!is.na(importance))),
              importance_sd_across_folds = sd(importance, na.rm = TRUE),
              importance_sd_within_fold = mean(importance_sd, na.rm = TRUE),
              p_value_fisher = tryCatch({ pvals <- p_value[!is.na(p_value)]; if (length(pvals) == 0) return(NA_real_)
              chi_stat <- -2 * sum(log(pmax(pvals, 1e-10))); pchisq(chi_stat, df = 2 * length(pvals), lower.tail = FALSE) }, error = function(e) NA_real_),
              p_value_mean = mean(p_value, na.rm = TRUE), n_folds = sum(!is.na(importance)), .groups = "drop") %>%
    mutate(p_value_bh = p.adjust(p_value_fisher, method = "BH"),
           significant_005 = p_value_bh < 0.05, significant_010 = p_value_bh < 0.10) %>%
    arrange(desc(importance))
  list(per_fold = all_folds, aggregated = aggregated, n_folds_cv = nfolds, n_permutations = n_permutations, n_observations = n)
}

run_permutation_importance_batch <- function(results_df, datasets_list, score_specific_datasets,
                                             n_permutations = 1000, nfolds = 5, n_cores = 1,
                                             covariate_cols = COVARIATE_COLS, seed = 42, models_to_run = NULL) {
  if (is.null(models_to_run)) models_to_run <- results_df %>% select(cvd_score, dataset_name) %>% distinct()
  n_models <- nrow(models_to_run)
  incremental_dir <- save_output("incremental_saves_importance")
  dir.create(incremental_dir, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("Models: %d | Perms/feature/fold: %d | Folds: %d | Cores: %d\n", n_models, n_permutations, nfolds, n_cores))
  all_results <- list()
  for (i in 1:n_models) {
    cvd_score <- models_to_run$cvd_score[i]; dataset_name <- models_to_run$dataset_name[i]
    cat(sprintf("\n[%d/%d] %s ~ %s\n", i, n_models, cvd_score, dataset_name))
    safe_dataset_name <- gsub("[^A-Za-z0-9_]", "_", dataset_name)
    incremental_file <- file.path(incremental_dir, sprintf("importance_%03d_%s_%s.rds", i, cvd_score, safe_dataset_name))
    if (file.exists(incremental_file)) { cat("  Loading from incremental save\n"); all_results[[i]] <- readRDS(incremental_file); next }
    model_results <- results_df %>% filter(cvd_score == !!cvd_score, dataset_name == !!dataset_name)
    if (nrow(model_results) == 0) { cat("  No elastic net results found\n"); next }
    features_to_test <- get_nonzero_features(model_results[1, ], covariate_cols)
    cat(sprintf("  Non-zero penalized features: %d\n", length(features_to_test)))
    if (length(features_to_test) == 0) { cat("  No penalized features with non-zero coefficients\n"); next }
    df_raw <- build_model_dataset(cvd_score, dataset_name, score_specific_datasets, datasets_list)
    if (is.null(df_raw)) { cat("  Could not build dataset\n"); next }
    tryCatch({
      importance_results <- run_permutation_importance(
        df_raw = df_raw, outcome_col = cvd_score, exclude_cols = "Sample_ID",
        features_to_test = features_to_test, covariate_cols = covariate_cols,
        alpha_grid = seq(0, 1, by = 0.1), nfolds = nfolds, n_permutations = n_permutations,
        n_cores = n_cores, seed = seed + i * 100000)
      if (!is.null(importance_results)) {
        importance_results$aggregated$cvd_score <- cvd_score
        importance_results$aggregated$dataset_name <- dataset_name
        importance_results$aggregated$n_observations <- importance_results$n_observations
        importance_results$aggregated$n_permutations <- importance_results$n_permutations
        importance_results$aggregated$n_folds_cv <- importance_results$n_folds_cv
        all_results[[i]] <- importance_results; saveRDS(importance_results, incremental_file)
        cat(sprintf("  ✓ Saved to %s\n", basename(incremental_file)))
      }
    }, error = function(e) cat(sprintf("  ERROR: %s\n", e$message)))
  }
  aggregated_all <- bind_rows(lapply(all_results, function(x) if (!is.null(x)) x$aggregated else NULL))
  cat(sprintf("\n✓ Completed %d models\n", length(all_results[!sapply(all_results, is.null)])))
  list(all_results = all_results, aggregated = aggregated_all)
}

extract_coefficient_directions <- function(results_df) {
  direction_list <- list()
  for (i in 1:nrow(results_df)) {
    pnames <- results_df$all_predictor_names[[i]]; coeffs <- results_df$all_coefficients[[i]]
    if (is.null(pnames) || is.null(coeffs)) next
    direction_list[[i]] <- tibble(cvd_score = results_df$cvd_score[i], dataset_name = results_df$dataset_name[i],
                                  feature = pnames, coefficient = coeffs,
                                  direction = case_when(coeffs > 0 ~ "Risk-Increasing", coeffs < 0 ~ "Risk-Decreasing", TRUE ~ "Zero"))
  }
  bind_rows(direction_list)
}

# ============================================
# 14. Execute permutation importance: Score-specific+ predictor block models
# ============================================

models_for_full_run <- results_score_plus_blocks %>% select(cvd_score, dataset_name) %>% distinct()
cat(sprintf("\nRunning %d score+block models\n", nrow(models_for_full_run)))
start_time <- Sys.time()
importance_results_full <- run_permutation_importance_batch(
  results_df = results_score_plus_blocks, datasets_list = datasets_list,
  score_specific_datasets = score_specific_datasets,
  n_permutations = 1000, nfolds = 5, n_cores = CPUS, covariate_cols = COVARIATE_COLS, seed = 42,
  models_to_run = models_for_full_run)
total_runtime_spb <- difftime(Sys.time(), start_time, units = "hours")
cat(sprintf("Score+block runtime: %.2f hours\n", as.numeric(total_runtime_spb)))

qs_save(importance_results_full, save_output("permutation_importance_full_results.qs2"))
saveRDS(importance_results_full, save_output("permutation_importance_full_results.rds"))
write.xlsx(importance_results_full$aggregated %>%
             select(cvd_score, dataset_name, feature, importance, importance_se, importance_sd_across_folds,
                    importance_sd_within_fold, p_value_fisher, p_value_mean, p_value_bh,
                    significant_005, significant_010, n_observations, n_permutations, n_folds_cv) %>%
             arrange(cvd_score, dataset_name, desc(importance)),
           save_output("permutation_importance_full_results.xlsx"))
coefficient_directions_full <- extract_coefficient_directions(results_score_plus_blocks)

# ============================================
# 15. Execute permutation importance: common (predictor blocks only) models
# ============================================

models_common_run <- results_common %>% select(cvd_score, dataset_name) %>% distinct()
cat(sprintf("\nRunning %d common models\n", nrow(models_common_run)))
start_time_common <- Sys.time()
importance_results_common <- run_permutation_importance_batch(
  results_df = results_common, datasets_list = datasets_list,
  score_specific_datasets = score_specific_datasets,
  n_permutations = 1000, nfolds = 5, n_cores = CPUS, covariate_cols = COVARIATE_COLS, seed = 42,
  models_to_run = models_common_run)
total_runtime_common <- difftime(Sys.time(), start_time_common, units = "hours")
cat(sprintf("Common runtime: %.2f hours\n", as.numeric(total_runtime_common)))

qs_save(importance_results_common, save_output("permutation_importance_common_results.qs2"))
saveRDS(importance_results_common, save_output("permutation_importance_common_results.rds"))
write.xlsx(importance_results_common$aggregated %>%
             select(cvd_score, dataset_name, feature, importance, importance_se, importance_sd_across_folds,
                    importance_sd_within_fold, p_value_fisher, p_value_mean, p_value_bh,
                    significant_005, significant_010, n_observations, n_permutations, n_folds_cv) %>%
             arrange(cvd_score, dataset_name, desc(importance)),
           save_output("permutation_importance_common_results.xlsx"))
coefficient_directions_common <- extract_coefficient_directions(results_common)

# ============================================
# 16. Feature importance thesis figures
# ============================================
extract_per_fold_data <- function(importance_results) {
  all_fold_data <- list()
  for (i in seq_along(importance_results$all_results)) {
    res <- importance_results$all_results[[i]]
    if (is.null(res)) next; fold_df <- res$per_fold
    if (is.null(fold_df) || nrow(fold_df) == 0) next
    fold_df$cvd_score <- res$aggregated$cvd_score[1]
    fold_df$dataset_name <- res$aggregated$dataset_name[1]
    all_fold_data[[i]] <- fold_df
  }
  bind_rows(all_fold_data)
}

feature_display_names <- c(
  "systolic" = "Systolic BP", "blood_pressure_treatment_Yes" = "Blood Pressure Treatment",
  "blood_pressure_treatment" = "Blood Pressure Treatment", "SmokingStatusQRISK3" = "Smoking",
  "mean_Total_Cholesterol_mg_dl" = "Total Cholesterol", "ratio_chol_hdl" = "Cholesterol/HDL",
  "creatinine" = "Creatinine", "mean_LDL_mg_dl" = "LDL cholesterol", "LDL.mg.dl" = "LDL", "TAG.mg.dl" = "Triacyglycerol",
  "mean_HDL_mg_dl" = "HDL cholesterol", "mean_hrt" = "Mean Heart Rate",
  "lpe_20_4" = "LPE (20:4)", "pe_o_22_2_20_4" = "PE-O (22:2/20:4)",
  "rbc_linoleic_18_2n6" = "RBC Linoleic Acid (18:2n-6)", "d_nervonic_24_1n9" = "DBS Nervonic Acid (24:1n-9)",
  "d_dgla_20_3n6" = "DBS DGLA (20:3n-6)", "rbc_dpa_22_5n3" = "RBC Docosapentaenoic Acid (22:5n-3)",
  "allantoin" = "Allantoin", "lactic_acid" = "Lactic Acid", "creatine" = "Creatine",
  "proline_betaine" = "Proline Betaine", "ecw_tbw" = "ECW/TBW",
  "bmr__basal_metabolic_rate_" = "Basal Metabolic Rate",
  "svr_skeletal_muscle_mass_visceral_fat_area_ratio_" = "Skeletal Muscle Mass/Visceral Fat Area",
  "smi__skeletal_muscle_index_" = "Skeletal Muscle Index", "vfa__visceral_fat_area_" = "Visceral Fat Area",
  "bfm_of_trunk" = "Body Fat Mass of Trunk", "Fat.free.mass" = "Fat-Free Mass", "Body.fat" = "Body Fat",
  "measured_circumference_of_right_thigh" = "Right Thigh Circumference",
  "measured_circumference_of_left_thigh" = "Left Thigh Circumference", "sleep_quality" = "Sleep Quality",
  "naps_during_day_yes" = "Naps During Day", "Education_Level_4" = "Doctoral",
  "Education_Level_3" = "University (Bachelor or Equivalent)",
  "annual_net_salary_4" = "Annual Net Salary >34,400\u201343,000",
  "Living_Status_Living.with.a.partner" = "Living With a Partner", "Employment_Status_Retired" = "Retired",
  "pe_o_19_1_20_5" = "PE-O (19:1/20:5)", "lps_18_0" = "LPS (18:0)",
  "rbc_n6_n3" = "RBC N6/N3", "rbc_margaric_17_0" = "RBC Margaric Acid (17:0)",
  "Fasted.glucose.mmol.l" = "Fasted Glucose", "CRP.mg.dl" = "CRP",
  "ALT.unit.L" = "ALT", "AST.unit.L" = "AST",
  "rbc_12_0" = "RBC Lauric Acid (12:0)")

capitalise_words <- function(s) gsub("(^|\\s)(\\w)", "\\1\\U\\2", s, perl = TRUE)
clean_feature_name <- function(x) {
  lookup_lower <- setNames(feature_display_names, tolower(names(feature_display_names)))
  sapply(x, function(feat) {
    if (feat %in% names(feature_display_names)) feature_display_names[feat]
    else if (tolower(feat) %in% names(lookup_lower)) lookup_lower[tolower(feat)]
    else capitalise_words(gsub("_", " ", gsub("\\.", " ", feat)))
  }, USE.NAMES = FALSE)
}

build_importance_panel <- function(fold_data, coefficient_directions, aggregated_df,
                                   block_name, result_type, score_colours, cvd_order,
                                   remove_score_specific = FALSE, score_specific_regex = NULL) {
  if (result_type == "score_plus_block") {
    plot_data <- fold_data %>% filter(grepl(paste0("_specific\\+", block_name, "$"), dataset_name))
  } else { plot_data <- fold_data %>% filter(dataset_name == block_name) }
  if (nrow(plot_data) == 0) return(NULL)
  if (remove_score_specific && !is.null(score_specific_regex)) {
    plot_data <- plot_data %>% filter(!grepl(score_specific_regex, feature))
    if (nrow(plot_data) == 0) return(NULL)
  }
  plot_data <- plot_data %>%
    left_join(coefficient_directions %>% select(cvd_score, dataset_name, feature, direction),
              by = c("cvd_score", "dataset_name", "feature")) %>%
    filter(!is.na(direction), direction != "Zero")
  if (nrow(plot_data) == 0) return(NULL)
  plot_data <- plot_data %>% mutate(
    cvd_label = recode(cvd_score, "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                       "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD"),
    cvd_label = factor(cvd_label, levels = cvd_order))
  plot_data <- plot_data %>% filter(importance > 0)
  fold_counts <- plot_data %>% group_by(feature, cvd_label) %>% summarise(n_nonzero = n(), .groups = "drop")
  valid_combos <- fold_counts %>% filter(n_nonzero >= 3)
  n_excluded <- nrow(fold_counts) - nrow(valid_combos)
  plot_data <- plot_data %>% semi_join(valid_combos, by = c("feature", "cvd_label"))
  if (nrow(plot_data) == 0) return(NULL)
  
  summary_data <- plot_data %>% group_by(feature, cvd_label, direction) %>%
    summarise(mean_imp = mean(importance, na.rm = TRUE),
              se_imp = sd(importance, na.rm = TRUE) / sqrt(sum(!is.na(importance))),
              n_folds = n(), .groups = "drop") %>%
    mutate(signed_mean = ifelse(direction == "Risk-Increasing", mean_imp, -mean_imp),
           ymin = ifelse(direction == "Risk-Increasing", mean_imp - se_imp, -(mean_imp + se_imp)),
           ymax = ifelse(direction == "Risk-Increasing", mean_imp + se_imp, -(mean_imp - se_imp)))
  plot_data <- plot_data %>% mutate(signed_importance = ifelse(direction == "Risk-Increasing", importance, -importance))
  
  if (result_type == "score_plus_block") {
    agg_filtered <- aggregated_df %>% filter(grepl(paste0("_specific\\+", block_name, "$"), dataset_name))
  } else { agg_filtered <- aggregated_df %>% filter(dataset_name == block_name) }
  sig_data <- agg_filtered %>% mutate(
    cvd_label = recode(cvd_score, "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                       "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD"),
    sig_star = case_when(p_value_bh < 0.001 ~ "***", p_value_bh < 0.01 ~ "**", p_value_bh < 0.05 ~ "*", TRUE ~ "")) %>%
    select(feature, cvd_label, sig_star)
  summary_data <- summary_data %>% left_join(sig_data, by = c("feature", "cvd_label")) %>%
    mutate(sig_star = ifelse(is.na(sig_star), "", sig_star))
  
  feat_inc <- summary_data %>% filter(direction == "Risk-Increasing") %>% group_by(feature) %>%
    summarise(max_imp = max(mean_imp), .groups = "drop") %>% arrange(desc(max_imp)) %>% pull(feature)
  feat_dec <- summary_data %>% filter(direction == "Risk-Decreasing") %>% group_by(feature) %>%
    summarise(max_imp = max(mean_imp), .groups = "drop") %>% arrange(desc(max_imp)) %>% pull(feature)
  feature_levels <- c(feat_inc, setdiff(feat_dec, feat_inc))
  clean_levels <- clean_feature_name(feature_levels)
  summary_data <- summary_data %>% mutate(feature_clean = factor(clean_feature_name(feature), levels = clean_levels)) %>% filter(feature %in% feature_levels)
  plot_data <- plot_data %>% mutate(feature_clean = factor(clean_feature_name(feature), levels = clean_levels)) %>% filter(feature %in% feature_levels)
  n_features <- length(feature_levels)
  
  # No-gap positioning
  bar_width <- 0.18; gap_between_features <- 0.3; position_rows <- list(); current_x <- 1; feature_centers <- numeric(n_features)
  for (f_idx in seq_along(feature_levels)) {
    feat <- feature_levels[f_idx]
    scores_present <- summary_data %>% filter(feature == feat) %>% arrange(cvd_label) %>% pull(cvd_label) %>% as.character()
    n_bars <- length(scores_present); if (n_bars == 0) next
    group_width <- n_bars * bar_width + (n_bars - 1) * 0.02
    start_x <- current_x - group_width / 2 + bar_width / 2; feature_centers[f_idx] <- current_x
    for (b_idx in seq_along(scores_present)) {
      x_pos <- start_x + (b_idx - 1) * (bar_width + 0.02)
      position_rows[[length(position_rows) + 1]] <- tibble(feature = feat, cvd_label = scores_present[b_idx], x_pos = x_pos)
    }
    current_x <- current_x + gap_between_features + group_width / 2 + bar_width
  }
  positions <- bind_rows(position_rows) %>% mutate(cvd_label = factor(cvd_label, levels = cvd_order))
  feature_tick_positions <- tibble(feature = feature_levels, feature_clean = factor(clean_levels, levels = clean_levels), x_tick = feature_centers)
  score_specific_stems <- c("Sex", "Age", "EthnicityCodeQRISK3", "SmokingStatusQRISK3", "diabetes2", "Diabetes",
                            "Weight_kg", "Height_cm", "blood_pressure_treatment", "Severe_mental_illness",
                            "ratio_chol_hdl", "systolic", "townsend", "race_ascvd",
                            "mean_Total_Cholesterol_mg_dl", "mean_HDL_mg_dl", "mean_LDL_mg_dl", "Risk.region")
  escaped_stems <- gsub("\\.", "\\\\.", score_specific_stems)
  ss_regex <- paste0("^(", paste(escaped_stems, collapse = "|"), ")")
  feature_tick_positions <- feature_tick_positions %>%
    mutate(is_score_specific = grepl(ss_regex, feature),
           label_face = ifelse(is_score_specific, "bold", "plain"))
  summary_data <- summary_data %>% left_join(positions, by = c("feature", "cvd_label"))
  plot_data <- plot_data %>% left_join(positions, by = c("feature", "cvd_label"))
  summary_data <- summary_data %>% mutate(star_y = ifelse(direction == "Risk-Increasing",
                                                          ymax + max(abs(signed_mean), na.rm = TRUE) * 0.05, ymin - max(abs(signed_mean), na.rm = TRUE) * 0.05))
  x_text_size <- if (n_features > 80) 10 else if (n_features > 50) 11 else if (n_features > 30) 12 else if (n_features > 15) 13 else 14
  
  p <- ggplot() +
    geom_col(data = summary_data, aes(x = x_pos, y = signed_mean, fill = cvd_label), width = bar_width, colour = "white", linewidth = 0.2) +
    geom_errorbar(data = summary_data %>% filter(!is.na(se_imp) & se_imp > 0), aes(x = x_pos, ymin = ymin, ymax = ymax), width = bar_width * 0.5, linewidth = 0.4, colour = "grey30") +
    geom_point(data = plot_data, aes(x = x_pos, y = signed_importance, fill = cvd_label), shape = 21, size = 1.0, alpha = 0.5, colour = "grey30", position = position_jitter(width = bar_width * 0.12, height = 0, seed = 42)) +
    geom_text(data = summary_data %>% filter(sig_star != ""), aes(x = x_pos, y = star_y, label = sig_star), size = 4, fontface = "bold", colour = "black") +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey30", linewidth = 0.5) +
    scale_fill_manual(values = score_colours, name = "CVD Risk Score", limits = cvd_order, drop = FALSE) +
    scale_x_continuous(breaks = feature_tick_positions$x_tick, labels = feature_tick_positions$feature_clean, expand = expansion(add = 0.5)) +
    labs(x = NULL, y = expression(Delta ~ "MAE (signed by effect direction)")) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.02))) +
    theme_minimal(base_size = 12) + theme(
      axis.text.x = element_text(size = x_text_size, angle = 30, hjust = 1, vjust = 1,
                                 face = feature_tick_positions$label_face),
      axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 13),
      legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = margin(t = 5, r = 8, b = 2, l = 8))
  y_range <- range(c(summary_data$ymin, summary_data$ymax, plot_data$signed_importance), na.rm = TRUE)
  list(plot = p, n_features = n_features, y_range = y_range, n_excluded = n_excluded)
}

assemble_combined_figure <- function(panel_a, panels_b, ylim_a, ylim_b, score_colours, cvd_order, output_file, fig_height = 12) {
  plot_a <- panel_a$plot + coord_cartesian(ylim = ylim_a, clip = "off") + labs(tag = "A") + theme(plot.tag = element_text(face = "bold", size = 16))
  b_list <- list(); panel_letters <- LETTERS[2:(1 + length(panels_b))]
  for (i in seq_along(panels_b)) {
    p_b <- panels_b[[i]]$plot + coord_cartesian(ylim = ylim_b, clip = "off") + labs(tag = panel_letters[i]) + theme(plot.tag = element_text(face = "bold", size = 16))
    if (i > 1) p_b <- p_b + labs(y = NULL)
    b_list[[i]] <- p_b
  }
  b_counts <- sapply(panels_b, function(x) x$n_features); b_widths <- pmax(b_counts / max(b_counts), 0.25)
  row_b <- wrap_plots(b_list, nrow = 1, widths = b_widths)
  dummy_data <- tibble(x = 1:4, y = 1:4, cvd_label = factor(cvd_order, levels = cvd_order))
  legend_plot <- ggplot(dummy_data, aes(x, y, fill = cvd_label)) + geom_col() +
    scale_fill_manual(values = score_colours, limits = cvd_order, name = "CVD Risk Score", drop = FALSE) +
    guides(fill = guide_legend(nrow = 1)) + theme(legend.position = "bottom", legend.title = element_text(face = "bold", size = 14), legend.text = element_text(size = 13))
  shared_legend <- cowplot::get_legend(legend_plot)
  height_ratio_b <- diff(ylim_b) / diff(ylim_a)
  combined <- (plot_a / row_b / wrap_elements(shared_legend)) + plot_layout(heights = c(1, height_ratio_b, 0.06))
  total_b <- sum(b_counts); fig_width <- min(50, max(16, max(panel_a$n_features, total_b) * 0.35))
  ggsave(save_output(output_file), combined, width = fig_width, height = fig_height, limitsize = FALSE)
  cat(sprintf("\n\u2713 Figure saved to %s (%.0f x %.0f in)\n", output_file, fig_width, fig_height))
  invisible(combined)
}


### Create figure 1: predictor blocks only
create_combined_figure_common <- function(importance_results, coefficient_directions,
                                          output_file = "combined_importance_common.pdf") {
  score_colours <- c("ASCVD" = "#F8766D", "Framingham" = "#7CAE00", "QRISK3" = "#00BFC4", "SCORE2" = "#C77CFF")
  cvd_order <- c("QRISK3", "SCORE2", "Framingham", "ASCVD")
  fold_data <- extract_per_fold_data(importance_results); aggregated_df <- importance_results$aggregated
  cat("\n--- Figure 1: Building Row A ---\n")
  panel_a <- build_importance_panel(fold_data, coefficient_directions, aggregated_df, "all_data", "common", score_colours, cvd_order)
  if (is.null(panel_a)) { cat("Row A: no results\n"); return(NULL) }
  cat(sprintf("  all_data: %d features\n", panel_a$n_features))
  cat("\n--- Figure 1: Building Row B ---\n")
  row_b_blocks <- c("body_composition", "sociodemographics_lifestyle", "urine_nmr")
  panels_b <- list()
  for (block in row_b_blocks) {
    res <- build_importance_panel(fold_data, coefficient_directions, aggregated_df, block, "common", score_colours, cvd_order)
    if (!is.null(res)) { panels_b[[block]] <- res; cat(sprintf("  %s: %d features\n", block, res$n_features)) }
  }
  if (length(panels_b) == 0) { cat("No Row B panels\n"); return(NULL) }
  assemble_combined_figure(panel_a, panels_b, ylim_a = c(-0.4, 0.4), ylim_b = c(-0.6, 0.4), score_colours, cvd_order, output_file, fig_height = 12)
}


### Create figure 2: score-specific + all predictor blocks
create_combined_figure_score_plus_block <- function(importance_results, coefficient_directions,
                                                    output_file = "combined_importance_score_plus_block.pdf") {
  score_colours <- c("ASCVD" = "#F8766D", "Framingham" = "#7CAE00", "QRISK3" = "#00BFC4", "SCORE2" = "#C77CFF")
  cvd_order <- c("QRISK3", "SCORE2", "Framingham", "ASCVD")
  score_specific_stems <- c("Sex", "Age", "EthnicityCodeQRISK3", "SmokingStatusQRISK3", "diabetes2", "Diabetes",
                            "Weight_kg", "Height_cm", "blood_pressure_treatment", "Severe_mental_illness",
                            "ratio_chol_hdl", "systolic", "townsend", "race_ascvd",
                            "mean_Total_Cholesterol_mg_dl", "mean_HDL_mg_dl", "mean_LDL_mg_dl", "Risk.region")
  escaped_stems <- gsub("\\.", "\\\\.", score_specific_stems)
  score_specific_regex <- paste0("^(", paste(escaped_stems, collapse = "|"), ")")
  fold_data <- extract_per_fold_data(importance_results); aggregated_df <- importance_results$aggregated
  cat("\n--- Figure 2: Building Row A (all features) ---\n")
  panel_a <- build_importance_panel(fold_data, coefficient_directions, aggregated_df, "all_data", "score_plus_block", score_colours, cvd_order, remove_score_specific = FALSE)
  if (is.null(panel_a)) { cat("Row A: no results\n"); return(NULL) }
  cat(sprintf("  all_data: %d features\n", panel_a$n_features))
  cat("\n--- Figure 2: Building Row B (score-specific excluded) ---\n")
  all_b_blocks <- c("lipids", "fatty_acids", "urine_nmr", "body_composition", "clinical_risk_factors", "sociodemographics_lifestyle")
  panels_b <- list()
  for (block in all_b_blocks) {
    cat(sprintf("  %s:\n", block))
    res <- build_importance_panel(fold_data, coefficient_directions, aggregated_df, block, "score_plus_block", score_colours, cvd_order, remove_score_specific = TRUE, score_specific_regex = score_specific_regex)
    if (!is.null(res)) { panels_b[[block]] <- res; cat(sprintf("    -> %d block-specific features\n", res$n_features)) }
    else cat("    -> No features passing filter (skipped)\n")
  }
  if (length(panels_b) == 0) { cat("No Row B panels\n"); return(NULL) }
  pad_a <- diff(panel_a$y_range) * 0.05
  ylim_a <- c(-0.5, ceiling((panel_a$y_range[2] + pad_a) * 10) / 10)
  assemble_combined_figure(panel_a, panels_b, ylim_a = ylim_a, ylim_b = c(-0.25, 0.25), score_colours, cvd_order, output_file, fig_height = 12)
}


### Run both figures
fig_common <- create_combined_figure_common(
  importance_results = importance_results_common,
  coefficient_directions = coefficient_directions_common,
  output_file = "combined_importance_common.png")

fig_spb <- create_combined_figure_score_plus_block(
  importance_results = importance_results_full,
  coefficient_directions = coefficient_directions_full,
  output_file = "combined_importance_score_plus_block.png")


### final save
qs_save(results_all, save_output("elastic_net_all_results.qs2"))
saveRDS(results_all, save_output("elastic_net_all_results.rds"))

# ============================================
# poster figure
# ============================================

# Make sure we have everything
fold_data <- extract_per_fold_data(importance_results_full)
aggregated_df <- importance_results_full$aggregated

score_colours <- c("ASCVD" = "#F8766D", "Framingham" = "#8AAD40", "QRISK3" = "#4DB8B0", "SCORE2" = "#C77CFF")
cvd_order <- c("ASCVD", "Framingham", "QRISK3", "SCORE2")

score_specific_stems <- c("Sex", "Age", "EthnicityCodeQRISK3", "SmokingStatusQRISK3", "diabetes2", "Diabetes",
                          "Weight_kg", "Height_cm", "blood_pressure_treatment", "Severe_mental_illness",
                          "ratio_chol_hdl", "systolic", "townsend", "race_ascvd",
                          "mean_Total_Cholesterol_mg_dl", "mean_HDL_mg_dl", "mean_LDL_mg_dl", "Risk.region")
escaped_stems <- gsub("\\.", "\\\\.", score_specific_stems)
score_specific_regex <- paste0("^(", paste(escaped_stems, collapse = "|"), ")")

# Row 1: Full predictor set (all features including score-specific)
panel_all <- build_importance_panel(fold_data, coefficient_directions_full, aggregated_df,
                                    "all_data", "score_plus_block", score_colours, cvd_order,
                                    remove_score_specific = FALSE)

# Row 2: Clinical measurements + Body composition (block-specific features only)
panel_clinical <- build_importance_panel(fold_data, coefficient_directions_full, aggregated_df,
                                         "clinical_risk_factors", "score_plus_block", score_colours, cvd_order,
                                         remove_score_specific = TRUE, score_specific_regex = score_specific_regex)

panel_body <- build_importance_panel(fold_data, coefficient_directions_full, aggregated_df,
                                     "body_composition", "score_plus_block", score_colours, cvd_order,
                                     remove_score_specific = TRUE, score_specific_regex = score_specific_regex)

# Row 3: Sociodemographic + Urine NMR
panel_socio <- build_importance_panel(fold_data, coefficient_directions_full, aggregated_df,
                                      "sociodemographics_lifestyle", "score_plus_block", score_colours, cvd_order,
                                      remove_score_specific = TRUE, score_specific_regex = score_specific_regex)

panel_urine <- build_importance_panel(fold_data, coefficient_directions_full, aggregated_df,
                                      "urine_nmr", "score_plus_block", score_colours, cvd_order,
                                      remove_score_specific = TRUE, score_specific_regex = score_specific_regex)

# Row 4: Fatty acids + Lipids
panel_fa <- build_importance_panel(fold_data, coefficient_directions_full, aggregated_df,
                                   "fatty_acids", "score_plus_block", score_colours, cvd_order,
                                   remove_score_specific = TRUE, score_specific_regex = score_specific_regex)

panel_lipids <- build_importance_panel(fold_data, coefficient_directions_full, aggregated_df,
                                       "lipids", "score_plus_block", score_colours, cvd_order,
                                       remove_score_specific = TRUE, score_specific_regex = score_specific_regex)

# Shared y-limits
ylim_all <- c(-0.5, ceiling((panel_all$y_range[2] + 0.05) * 10) / 10)
ylim_blocks <- c(-0.25, 0.25)

# Build plots with coord_cartesian and tags
p_all <- panel_all$plot + coord_cartesian(ylim = ylim_all, clip = "off") +
  labs(tag = "A") + theme(plot.tag = element_text(face = "bold", size = 14))

p_clinical <- panel_clinical$plot + coord_cartesian(ylim = ylim_blocks, clip = "off") +
  labs(tag = "B") + theme(plot.tag = element_text(face = "bold", size = 14))
p_body <- panel_body$plot + coord_cartesian(ylim = ylim_blocks, clip = "off") +
  labs(tag = "C", y = NULL) + theme(plot.tag = element_text(face = "bold", size = 14))

p_socio <- panel_socio$plot + coord_cartesian(ylim = ylim_blocks, clip = "off") +
  labs(tag = "D") + theme(plot.tag = element_text(face = "bold", size = 14))
p_urine <- panel_urine$plot + coord_cartesian(ylim = ylim_blocks, clip = "off") +
  labs(tag = "E", y = NULL) + theme(plot.tag = element_text(face = "bold", size = 14))

p_fa <- panel_fa$plot + coord_cartesian(ylim = ylim_blocks, clip = "off") +
  labs(tag = "F") + theme(plot.tag = element_text(face = "bold", size = 14))
p_lipids <- panel_lipids$plot + coord_cartesian(ylim = ylim_blocks, clip = "off") +
  labs(tag = "G", y = NULL) + theme(plot.tag = element_text(face = "bold", size = 14))

# Shared legend
dummy_data <- tibble(x = 1:4, y = 1:4, cvd_label = factor(cvd_order, levels = cvd_order))
legend_plot <- ggplot(dummy_data, aes(x, y, fill = cvd_label)) + geom_col() +
  scale_fill_manual(values = score_colours, limits = cvd_order, name = "CVD Risk Score", drop = FALSE) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom", 
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.6, "cm"))
shared_legend <- cowplot::get_legend(legend_plot)

# Width ratios for paired rows (proportional to feature count)
w_r2 <- c(panel_clinical$n_features, panel_body$n_features)
w_r3 <- c(panel_socio$n_features, panel_urine$n_features)
w_r4 <- c(panel_fa$n_features, panel_lipids$n_features)

# Assemble: tall narrow layout
poster_importance <- (p_all) /
  (p_clinical | p_body) /
  (p_socio | p_urine) /
  (p_fa | p_lipids) /
  wrap_elements(shared_legend) +
  plot_layout(
    heights = c(1.2, 0.6, 0.6, 0.6, 0.08),
    widths = list(NULL, w_r2, w_r3, w_r4, NULL)
  )

row2 <- (p_clinical | p_body) + plot_layout(widths = c(panel_clinical$n_features, panel_body$n_features))
row3 <- (p_socio | p_urine) + plot_layout(widths = c(panel_socio$n_features, panel_urine$n_features))
row4 <- (p_fa | p_lipids) + plot_layout(widths = c(panel_fa$n_features, panel_lipids$n_features))

poster_importance <- p_all / row2 / row3 / row4 / wrap_elements(shared_legend) +
  plot_layout(heights = c(1.2, 0.6, 0.6, 0.6, 0.12))

ggsave(
  filename = file.path(figures_path, "importance_poster.png"),
  plot = poster_importance,
  width = 10,
  height = 16,
  dpi = 600,
  bg = "white"
)


# ============================================
# Poster: Dot plot predictive performance by predictor block
# ============================================
library(ggnewscale)
dataset_order_dotplot <- c("all_data", "lipids", "fatty_acids", "urine_nmr",
                           "body_composition", "sociodemographics_lifestyle", "clinical_risk_factors")

dataset_labels_dotplot <- c(
  "all_data" = "Multi-Domain\nPredictor Set",
  "lipids" = "Lipid Profiles",
  "fatty_acids" = "Fatty Acids",
  "urine_nmr" = "Urinary\nMetabolites",
  "body_composition" = "Body Composition",
  "sociodemographics_lifestyle" = "Sociodemographic\n& Lifestyle",
  "clinical_risk_factors" = "Clinical Measurements"
)

cvd_score_order <- c("ascvd_10y", "frs_10y", "QRISK3_risk", "SCORE2_score")

cvd_score_labels_dotplot <- c("SCORE2_score" = "SCORE2", "QRISK3_risk" = "QRISK3",
                              "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD")

cvd_shapes <- c("ASCVD" = 21, "Framingham" = 22, "QRISK3" = 24, "SCORE2" = 23)

dotplot_data <- results_common %>%
  filter(dataset_name %in% dataset_order_dotplot,
         cvd_score %in% cvd_score_order) %>%
  mutate(
    dataset_label = factor(dataset_labels_dotplot[dataset_name],
                           levels = rev(dataset_labels_dotplot[dataset_order_dotplot])),
    cvd_label = cvd_score_labels_dotplot[cvd_score],
    cvd_label = factor(cvd_label, levels = c("ASCVD", "Framingham", "QRISK3", "SCORE2")),
    significant = !is.na(permutation_p_value) & permutation_p_value < 0.05,
    Q2_fold_sd = replace_na(Q2_fold_sd, 0),
    fill_group = ifelse(significant, as.character(cvd_label), "ns")
  )

# Override significance (permutations not run in this batch)
dotplot_data <- dotplot_data %>%
  mutate(
    significant = case_when(
      dataset_name == "body_composition" ~ TRUE,
      dataset_name == "all_data" & cvd_score == "QRISK3_risk" ~ TRUE,
      TRUE ~ FALSE
    ),
    fill_group = ifelse(significant, as.character(cvd_label), "ns")
  )

p_dotplot <- ggplot(dotplot_data, aes(x = Q2_Y, y = dataset_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.6) +
  geom_hline(yintercept = seq(1.5, length(dataset_order_dotplot) - 0.5, by = 1),
             colour = "grey90", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = Q2_Y - Q2_fold_sd, xmax = Q2_Y + Q2_fold_sd,
                     group = cvd_label),
                 height = 0, linewidth = 0.4, colour = "grey60",
                 position = position_dodge(width = 0.6)) +
  geom_point(aes(shape = cvd_label,
                 fill = fill_group,
                 colour = cvd_label,
                 alpha = ifelse(significant, "p < 0.05", "p ≥ 0.05")),
             size = 3.5, stroke = 0.8,
             position = position_dodge(width = 0.6)) +
  scale_colour_manual(
    values = c("ASCVD" = "#F8766D", "Framingham" = "#8AAD40",
               "QRISK3" = "#4DB8B0", "SCORE2" = "#C77CFF"),
    name = "CVD Risk Score"
  ) +
  scale_fill_manual(
    values = c("ASCVD" = "#F8766D", "Framingham" = "#8AAD40",
               "QRISK3" = "#4DB8B0", "SCORE2" = "#C77CFF",
               "ns" = "white"),
    guide = "none"
  ) +
  scale_shape_manual(values = cvd_shapes, name = "CVD Risk Score") +
  scale_alpha_manual(
    values = c("p < 0.05" = 1, "p ≥ 0.05" = 1),
    name = "Significance",
    guide = guide_legend(
      override.aes = list(
        shape = 24,
        fill = c("grey30", "white"),
        colour = "grey30",
        size = 4
      )
    )
  ) +
  scale_x_continuous(breaks = seq(-0.20, 0.20, by = 0.10)) +
  labs(x = expression(Q^2 ~ "(cross-validated " * R^2 * ")"), y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.justification = "left",
    legend.margin = margin(0, 0, 0, -110),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.spacing.y = unit(4, "pt"),
    legend.box.spacing = unit(2, "pt"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 5)
  ) +
  guides(
    colour = guide_legend(nrow = 1, override.aes = list(size = 4)),
    shape = guide_legend(nrow = 1, override.aes = list(size = 4))
  ) +
  coord_cartesian(xlim = c(-0.30, 0.25), clip = "off")

ggsave(save_output("poster_dotplot_predictive_performance.png"),
       plot = p_dotplot, width = 5, height = 6.5, dpi = 600, bg = "white")

cat("\n✓ Poster dot plot saved\n")




# ============================================
# Poster: Body composition feature importance (standalone)
# ============================================

score_colours_poster <- c("ASCVD" = "#F8766D", "Framingham" = "#8AAD40", "QRISK3" = "#4DB8B0", "SCORE2" = "#C77CFF")
cvd_order_poster <- c("ASCVD", "Framingham", "QRISK3", "SCORE2")

# Extract fold data and coefficient directions
fold_data_body <- extract_per_fold_data(importance_results_common)
aggregated_df_body <- importance_results_common$aggregated
coef_directions_body <- extract_coefficient_directions(results_common)

# Filter to body composition only
plot_data_body <- fold_data_body %>%
  filter(dataset_name == "body_composition") %>%
  left_join(coef_directions_body %>%
              filter(dataset_name == "body_composition") %>%
              select(cvd_score, feature, direction),
            by = c("cvd_score", "feature")) %>%
  filter(!is.na(direction), direction != "Zero") %>%
  mutate(
    cvd_label = recode(cvd_score, "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                       "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD"),
    cvd_label = factor(cvd_label, levels = cvd_order_poster)
  ) %>%
  filter(importance > 0)

# Require feature present in >= 3 folds
fold_counts <- plot_data_body %>%
  group_by(feature, cvd_label) %>%
  summarise(n_nonzero = n(), .groups = "drop")
valid_combos <- fold_counts %>% filter(n_nonzero >= 3)
plot_data_body <- plot_data_body %>% semi_join(valid_combos, by = c("feature", "cvd_label"))

# Summarise across folds
summary_body <- plot_data_body %>%
  group_by(feature, cvd_label, direction) %>%
  summarise(
    mean_imp = mean(importance, na.rm = TRUE),
    se_imp = sd(importance, na.rm = TRUE) / sqrt(sum(!is.na(importance))),
    .groups = "drop"
  ) %>%
  mutate(
    signed_mean = ifelse(direction == "Risk-Increasing", mean_imp, -mean_imp),
    ymin = ifelse(direction == "Risk-Increasing", mean_imp - se_imp, -(mean_imp + se_imp)),
    ymax = ifelse(direction == "Risk-Increasing", mean_imp + se_imp, -(mean_imp - se_imp))
  )

# Significance from aggregated results
sig_body <- aggregated_df_body %>%
  filter(dataset_name == "body_composition") %>%
  mutate(
    cvd_label = recode(cvd_score, "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                       "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD"),
    sig_star = case_when(
      p_value_bh < 0.001 ~ "***",
      p_value_bh < 0.01 ~ "**",
      p_value_bh < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(feature, cvd_label, sig_star)

summary_body <- summary_body %>%
  left_join(sig_body, by = c("feature", "cvd_label")) %>%
  mutate(sig_star = ifelse(is.na(sig_star), "", sig_star))

# Order features: risk-increasing first (by max importance), then risk-decreasing
feat_inc <- summary_body %>%
  filter(direction == "Risk-Increasing") %>%
  group_by(feature) %>%
  summarise(max_imp = max(mean_imp), .groups = "drop") %>%
  arrange(desc(max_imp)) %>%
  pull(feature)
feat_dec <- summary_body %>%
  filter(direction == "Risk-Decreasing") %>%
  group_by(feature) %>%
  summarise(max_imp = max(mean_imp), .groups = "drop") %>%
  arrange(desc(max_imp)) %>%
  pull(feature)
feature_levels <- c(feat_inc, setdiff(feat_dec, feat_inc))
clean_levels <- clean_feature_name(feature_levels)

summary_body <- summary_body %>%
  filter(feature %in% feature_levels) %>%
  mutate(feature_clean = factor(clean_feature_name(feature), levels = clean_levels))

n_features <- length(feature_levels)

# Build x positions (no-gap grouped bars)
bar_width <- 0.18
gap_between_features <- 0.3
position_rows <- list()
current_x <- 1
feature_centers <- numeric(n_features)

for (f_idx in seq_along(feature_levels)) {
  feat <- feature_levels[f_idx]
  scores_present <- summary_body %>%
    filter(feature == feat) %>%
    arrange(cvd_label) %>%
    pull(cvd_label) %>%
    as.character()
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
  mutate(cvd_label = factor(cvd_label, levels = cvd_order_poster))

feature_tick_positions <- tibble(
  feature = feature_levels,
  feature_clean = factor(clean_levels, levels = clean_levels),
  x_tick = feature_centers
)

# Join positions
summary_body <- summary_body %>% left_join(positions, by = c("feature", "cvd_label"))

# Star y positions
summary_body <- summary_body %>%
  mutate(star_y = ifelse(
    direction == "Risk-Increasing",
    ymax + 0.04,
    ymin - 0.04
  ))

# Build plot
p_body_poster <- ggplot() +
  geom_col(data = summary_body,
           aes(x = x_pos, y = signed_mean, fill = cvd_label),
           width = bar_width, colour = "white", linewidth = 0.2) +
  geom_errorbar(data = summary_body %>% filter(!is.na(se_imp) & se_imp > 0),
                aes(x = x_pos, ymin = ymin, ymax = ymax),
                width = bar_width * 0.5, linewidth = 0.4, colour = "grey30") +
  geom_text(data = summary_body %>% filter(sig_star != ""),
            aes(x = x_pos, y = star_y, label = sig_star),
            size = 6, fontface = "bold", colour = "black") +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey30", linewidth = 0.5) +
  scale_fill_manual(values = score_colours_poster, name = "CVD Risk Score",
                    limits = cvd_order_poster, drop = FALSE) +
  scale_x_continuous(breaks = feature_tick_positions$x_tick,
                     labels = feature_tick_positions$feature_clean,
                     expand = expansion(add = 0.5)) +
  labs(x = NULL, y = expression(Delta ~ "MAE (signed by effect direction)")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 15, angle = 30, hjust = 1, vjust = 1),  # was 12
    axis.text.y = element_text(size = 14),    # was 11
    axis.title.y = element_text(size = 16),   # was 13
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 13),  # was 10
    legend.text = element_text(size = 12),    # was 9
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 25, r = 10, b = 5, l = 10)
  ) +
  guides(fill = guide_legend(nrow = 1))

ggsave(save_output("poster_body_composition_importance.png"),
       plot = p_body_poster, width = 7.5, height = 5, dpi = 600, bg = "white")

cat("\n✓ Body composition importance plot saved\n")






# ============================================
# Poster: Feature importance (A on top, combined B-G below)
# ============================================

feature_display_names["systolic"] <- "Systolic Blood Pressure"

score_colours_fi <- c("ASCVD" = "#F8766D", "Framingham" = "#8AAD40", "QRISK3" = "#4DB8B0", "SCORE2" = "#C77CFF")
cvd_order_fi <- c("ASCVD", "Framingham", "QRISK3", "SCORE2")

fold_data_fi <- extract_per_fold_data(importance_results_full)
aggregated_df_fi <- importance_results_full$aggregated
coef_directions_fi <- extract_coefficient_directions(results_score_plus_blocks)

score_specific_stems <- c("Sex", "Age", "EthnicityCodeQRISK3", "SmokingStatusQRISK3", "diabetes2", "Diabetes",
                          "Weight_kg", "Height_cm", "blood_pressure_treatment", "Severe_mental_illness",
                          "ratio_chol_hdl", "systolic", "townsend", "race_ascvd",
                          "mean_Total_Cholesterol_mg_dl", "mean_HDL_mg_dl", "mean_LDL_mg_dl", "Risk.region")
escaped_stems <- gsub("\\.", "\\\\.", score_specific_stems)
score_specific_regex <- paste0("^(", paste(escaped_stems, collapse = "|"), ")")

domain_labels_fi <- c(
  "body_composition" = "Body\nComposition",
  "clinical_risk_factors" = "Clinical\nMeasurements",
  "sociodemographics_lifestyle" = "Sociodemographic\n& Lifestyle",
  "fatty_acids" = "Fatty\nAcids",
  "urine_nmr" = "Urinary\nMetabolites",
  "lipids" = "Lipid\nProfiles"
)

# ---- Helper: process one block ----
process_block_fi <- function(block_name, fold_data, coef_directions, aggregated_df,
                             cvd_order, remove_ss = TRUE, ss_regex = NULL) {
  plot_data <- fold_data %>%
    filter(grepl(paste0("_specific\\+", block_name, "$"), dataset_name))
  if (nrow(plot_data) == 0) return(NULL)
  if (remove_ss && !is.null(ss_regex)) {
    plot_data <- plot_data %>% filter(!grepl(ss_regex, feature))
    if (nrow(plot_data) == 0) return(NULL)
  }
  plot_data <- plot_data %>%
    left_join(coef_directions %>% select(cvd_score, dataset_name, feature, direction),
              by = c("cvd_score", "dataset_name", "feature")) %>%
    filter(!is.na(direction), direction != "Zero") %>%
    mutate(cvd_label = recode(cvd_score, "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                              "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD"),
           cvd_label = factor(cvd_label, levels = cvd_order)) %>%
    filter(importance > 0)
  fold_counts <- plot_data %>% group_by(feature, cvd_label) %>%
    summarise(n_nonzero = n(), .groups = "drop")
  valid_combos <- fold_counts %>% filter(n_nonzero >= 3)
  plot_data <- plot_data %>% semi_join(valid_combos, by = c("feature", "cvd_label"))
  if (nrow(plot_data) == 0) return(NULL)
  
  summary <- plot_data %>%
    group_by(feature, cvd_label, direction) %>%
    summarise(mean_imp = mean(importance, na.rm = TRUE),
              se_imp = sd(importance, na.rm = TRUE) / sqrt(sum(!is.na(importance))),
              .groups = "drop") %>%
    mutate(signed_mean = ifelse(direction == "Risk-Increasing", mean_imp, -mean_imp),
           ymin = ifelse(direction == "Risk-Increasing", mean_imp - se_imp, -(mean_imp + se_imp)),
           ymax = ifelse(direction == "Risk-Increasing", mean_imp + se_imp, -(mean_imp - se_imp)))
  
  agg_filtered <- aggregated_df %>%
    filter(grepl(paste0("_specific\\+", block_name, "$"), dataset_name))
  sig_data <- agg_filtered %>%
    mutate(cvd_label = recode(cvd_score, "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                              "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD"),
           sig_star = case_when(p_value_bh < 0.001 ~ "***", p_value_bh < 0.01 ~ "**",
                                p_value_bh < 0.05 ~ "*", TRUE ~ "")) %>%
    select(feature, cvd_label, sig_star)
  
  summary <- summary %>%
    left_join(sig_data, by = c("feature", "cvd_label")) %>%
    mutate(sig_star = ifelse(is.na(sig_star), "", sig_star),
           domain = block_name)
  summary
}

# ---- Panel A: full predictor set ----
panel_a_data <- fold_data_fi %>%
  filter(grepl("_specific\\+all_data$", dataset_name)) %>%
  left_join(coef_directions_fi %>% select(cvd_score, dataset_name, feature, direction),
            by = c("cvd_score", "dataset_name", "feature")) %>%
  filter(!is.na(direction), direction != "Zero") %>%
  mutate(cvd_label = recode(cvd_score, "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                            "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD"),
         cvd_label = factor(cvd_label, levels = cvd_order_fi)) %>%
  filter(importance > 0)

fold_counts_a <- panel_a_data %>% group_by(feature, cvd_label) %>%
  summarise(n_nonzero = n(), .groups = "drop")
valid_a <- fold_counts_a %>% filter(n_nonzero >= 3)
panel_a_data <- panel_a_data %>% semi_join(valid_a, by = c("feature", "cvd_label"))

summary_a <- panel_a_data %>%
  group_by(feature, cvd_label, direction) %>%
  summarise(mean_imp = mean(importance, na.rm = TRUE),
            se_imp = sd(importance, na.rm = TRUE) / sqrt(sum(!is.na(importance))),
            .groups = "drop") %>%
  mutate(signed_mean = ifelse(direction == "Risk-Increasing", mean_imp, -mean_imp),
         ymin = ifelse(direction == "Risk-Increasing", mean_imp - se_imp, -(mean_imp + se_imp)),
         ymax = ifelse(direction == "Risk-Increasing", mean_imp + se_imp, -(mean_imp - se_imp)),
         feature_clean = clean_feature_name(feature))

sig_a <- aggregated_df_fi %>%
  filter(grepl("_specific\\+all_data$", dataset_name)) %>%
  mutate(cvd_label = recode(cvd_score, "QRISK3_risk" = "QRISK3", "SCORE2_score" = "SCORE2",
                            "frs_10y" = "Framingham", "ascvd_10y" = "ASCVD"),
         sig_star = case_when(p_value_bh < 0.001 ~ "***", p_value_bh < 0.01 ~ "**",
                              p_value_bh < 0.05 ~ "*", TRUE ~ "")) %>%
  select(feature, cvd_label, sig_star)

summary_a <- summary_a %>%
  left_join(sig_a, by = c("feature", "cvd_label")) %>%
  mutate(sig_star = ifelse(is.na(sig_star), "", sig_star),
         is_score_specific = grepl(score_specific_regex, feature))

# Order Panel A features
feat_inc_a <- summary_a %>% filter(direction == "Risk-Increasing") %>%
  group_by(feature) %>% summarise(max_imp = max(mean_imp), .groups = "drop") %>%
  arrange(desc(max_imp)) %>% pull(feature)
feat_dec_a <- summary_a %>% filter(direction == "Risk-Decreasing") %>%
  group_by(feature) %>% summarise(max_imp = max(mean_imp), .groups = "drop") %>%
  arrange(desc(max_imp)) %>% pull(feature)
feature_levels_a <- c(feat_inc_a, setdiff(feat_dec_a, feat_inc_a))
clean_levels_a <- clean_feature_name(feature_levels_a)

summary_a <- summary_a %>%
  filter(feature %in% feature_levels_a) %>%
  mutate(feature_clean = factor(clean_feature_name(feature), levels = clean_levels_a))

n_feat_a <- length(feature_levels_a)

# Panel A x-positions
bar_width <- 0.18
gap_between <- 0.3
pos_rows_a <- list(); current_x <- 1; centers_a <- numeric(n_feat_a)
for (f_idx in seq_along(feature_levels_a)) {
  feat <- feature_levels_a[f_idx]
  scores_present <- summary_a %>% filter(feature == feat) %>% arrange(cvd_label) %>%
    pull(cvd_label) %>% as.character()
  n_bars <- length(scores_present); if (n_bars == 0) next
  group_width <- n_bars * bar_width + (n_bars - 1) * 0.02
  start_x <- current_x - group_width / 2 + bar_width / 2
  centers_a[f_idx] <- current_x
  for (b_idx in seq_along(scores_present)) {
    x_pos <- start_x + (b_idx - 1) * (bar_width + 0.02)
    pos_rows_a[[length(pos_rows_a) + 1]] <- tibble(feature = feat, cvd_label = scores_present[b_idx], x_pos = x_pos)
  }
  current_x <- current_x + gap_between + group_width / 2 + bar_width
}
positions_a <- bind_rows(pos_rows_a) %>% mutate(cvd_label = factor(cvd_label, levels = cvd_order_fi))
ticks_a <- tibble(feature = feature_levels_a, feature_clean = factor(clean_levels_a, levels = clean_levels_a), x_tick = centers_a)
ticks_a <- ticks_a %>% mutate(is_ss = grepl(score_specific_regex, feature),
                              label_face = ifelse(is_ss, "bold", "plain"))

summary_a <- summary_a %>% left_join(positions_a, by = c("feature", "cvd_label"))
summary_a <- summary_a %>% mutate(star_y = ifelse(direction == "Risk-Increasing", ymax + 0.06, ymin - 0.06))

p_a <- ggplot() +
  geom_col(data = summary_a, aes(x = x_pos, y = signed_mean, fill = cvd_label),
           width = bar_width, colour = "white", linewidth = 0.2) +
  geom_errorbar(data = summary_a %>% filter(!is.na(se_imp) & se_imp > 0),
                aes(x = x_pos, ymin = ymin, ymax = ymax),
                width = bar_width * 0.5, linewidth = 0.4, colour = "grey30") +
  geom_text(data = summary_a %>% filter(sig_star != ""),
            aes(x = x_pos, y = star_y, label = sig_star),
            size = 4, fontface = "bold", colour = "black") +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey30", linewidth = 0.5) +
  scale_fill_manual(values = score_colours_fi, name = "CVD Risk Score",
                    limits = cvd_order_fi, drop = FALSE) +
  scale_x_continuous(breaks = ticks_a$x_tick, labels = ticks_a$feature_clean,
                     expand = expansion(add = 0.5)) +
  labs(x = NULL, y = expression(Delta ~ "MAE"), title = "Multi-Domain Predictors") +
  theme_minimal(base_size = 11) +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 12, angle = 30, hjust = 1, vjust = 1,
                                   face = ticks_a$label_face),  # was 9
        axis.text.y = element_text(size = 12),    # was 10
        axis.title.y = element_text(size = 14),   # was 11
        plot.title = element_text(size = 11, face = "plain", colour = "grey30"),  # was 9
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 5, r = 8, b = 0, l = 8))

# ---- Combined B-G panel ----
blocks <- c("body_composition", "clinical_risk_factors", "sociodemographics_lifestyle",
            "fatty_acids", "urine_nmr", "lipids")

all_blocks <- bind_rows(lapply(blocks, function(b) {
  process_block_fi(b, fold_data_fi, coef_directions_fi, aggregated_df_fi,
                   cvd_order_fi, remove_ss = TRUE, ss_regex = score_specific_regex)
}))

all_blocks <- all_blocks %>%
  mutate(feature_clean = clean_feature_name(feature),
         domain_label = domain_labels_fi[domain],
         domain_label = factor(domain_label, levels = unname(domain_labels_fi)))

# Order features within domains
block_feature_order <- all_blocks %>%
  group_by(domain, feature) %>%
  summarise(max_imp = max(abs(signed_mean)), .groups = "drop") %>%
  arrange(domain, desc(max_imp))

# Build x-positions with domain gaps
pos_rows_bg <- list(); current_x <- 1; centers_bg <- numeric(0)
feature_list_bg <- character(0); domain_boundaries <- list()

for (d in blocks) {
  feats_in_domain <- block_feature_order %>% filter(domain == d) %>% pull(feature)
  if (length(feats_in_domain) == 0) next
  domain_start <- current_x
  
  for (feat in feats_in_domain) {
    scores_present <- all_blocks %>% filter(feature == feat, domain == d) %>%
      arrange(cvd_label) %>% pull(cvd_label) %>% as.character()
    n_bars <- length(scores_present); if (n_bars == 0) next
    group_width <- n_bars * bar_width + (n_bars - 1) * 0.02
    start_x <- current_x - group_width / 2 + bar_width / 2
    centers_bg <- c(centers_bg, current_x)
    feature_list_bg <- c(feature_list_bg, feat)
    for (b_idx in seq_along(scores_present)) {
      x_pos <- start_x + (b_idx - 1) * (bar_width + 0.02)
      pos_rows_bg[[length(pos_rows_bg) + 1]] <- tibble(
        feature = feat, cvd_label = scores_present[b_idx], x_pos = x_pos, domain = d)
    }
    current_x <- current_x + gap_between + group_width / 2 + bar_width
  }
  domain_end <- current_x - gap_between
  domain_boundaries[[d]] <- c(start = domain_start, end = domain_end,
                              mid = (domain_start + domain_end) / 2)
  current_x <- current_x + 0.3  # extra gap between domains
}

positions_bg <- bind_rows(pos_rows_bg) %>% mutate(cvd_label = factor(cvd_label, levels = cvd_order_fi))
clean_list_bg <- clean_feature_name(feature_list_bg)
ticks_bg <- tibble(feature = feature_list_bg,
                   feature_clean = clean_list_bg,
                   x_tick = centers_bg)

all_blocks <- all_blocks %>% left_join(positions_bg, by = c("feature", "cvd_label", "domain"))
all_blocks <- all_blocks %>% mutate(star_y = ifelse(direction == "Risk-Increasing", ymax + 0.03, ymin - 0.03))

# Domain label positions
domain_annot <- bind_rows(lapply(names(domain_boundaries), function(d) {
  tibble(domain = d,
         x_mid = domain_boundaries[[d]]["mid"],
         x_start = domain_boundaries[[d]]["start"] - 0.3,
         x_end = domain_boundaries[[d]]["end"] + 0.1,
         label = domain_labels_fi[d])
}))

# Y position for domain labels (above the plot)
y_top_bg <- max(c(all_blocks$ymax, all_blocks$signed_mean), na.rm = TRUE) + 0.08

p_bg <- ggplot() +
  geom_col(data = all_blocks, aes(x = x_pos, y = signed_mean, fill = cvd_label),
           width = bar_width, colour = "white", linewidth = 0.2) +
  geom_errorbar(data = all_blocks %>% filter(!is.na(se_imp) & se_imp > 0),
                aes(x = x_pos, ymin = ymin, ymax = ymax),
                width = bar_width * 0.5, linewidth = 0.4, colour = "grey30") +
  geom_text(data = all_blocks %>% filter(sig_star != ""),
            aes(x = x_pos, y = star_y, label = sig_star),
            size = 3.5, fontface = "bold", colour = "black") +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey30", linewidth = 0.5) +
  # Domain labels on top
  geom_text(data = domain_annot, aes(x = x_mid, y = y_top_bg, label = label),
            size = 3.5, fontface = "plain", colour = "grey30", lineheight = 0.85) +  # was 3
  # Domain separator lines
  geom_segment(data = domain_annot %>% filter(row_number() > 1),
               aes(x = x_start - 0.1, xend = x_start - 0.1,
                   y = -Inf, yend = Inf),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = score_colours_fi, name = "CVD Risk Score",
                    limits = cvd_order_fi, drop = FALSE) +
  scale_x_continuous(breaks = ticks_bg$x_tick, labels = ticks_bg$feature_clean,
                     expand = expansion(add = c(1.5, 0.5))) +
  labs(x = NULL, y = expression(Delta ~ "MAE")) +
  theme_minimal(base_size = 11) +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 11, angle = 30, hjust = 1, vjust = 1),  # was 8
        axis.text.y = element_text(size = 12),    # was 10
        axis.title.y = element_text(size = 14),   # was 11
        plot.title = element_text(size = 11, face = "bold"),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 20, r = 8, b = 2, l = 8)) +
  coord_cartesian(clip = "off")

# ---- Shared legend ----
dummy_legend_fi <- ggplot(data.frame(x = 1:4, y = 1:4,
                                     cvd = factor(cvd_order_fi, levels = cvd_order_fi)),
                          aes(x, y, fill = cvd)) +
  geom_col() +
  scale_fill_manual(values = score_colours_fi, name = "CVD Risk Score",
                    limits = cvd_order_fi, drop = FALSE) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 13),  # was 9
        legend.text = element_text(size = 12))  # was 8
shared_legend_fi <- cowplot::get_legend(dummy_legend_fi)

# ---- Combine ----
poster_fi_combined <- p_a / p_bg / wrap_elements(shared_legend_fi) +
  plot_layout(heights = c(1, 0.5, 0.08))

ggsave(save_output("poster_feature_importance_combined.png"),
       plot = poster_fi_combined, width = 12, height = 7, dpi = 600, bg = "white")

cat("\n✓ Poster feature importance combined figure saved\n")



cat("\n✓ Poster dot plot saved\n")