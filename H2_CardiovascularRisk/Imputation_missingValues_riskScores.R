### Imputation of missing values in risk score calculation plus risk scores themselves (missForest)
### Author: Luisa Delius

# Load packages
library(tidyverse)
library(missForest)
library(parallel)
library(doParallel)

# Set working directory
wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

# Upload the data
df_all_risk_scores <- readRDS(file.path(wkdir, "processed_data", "df_all_risk_scores.rds"))
df_qrisk3_input <- readRDS(file.path(wkdir, "processed_data", "QRISK3_calculation_input.rds"))
df_ASCVD_Framingham_SCORE2_input <- readRDS(file.path(wkdir, "processed_data", "ASCVD_SCORE2_Framingham_input.rds"))

# Join the dfs and make sure, age etc are in only once
df_joined_for_imputation <- df_qrisk3_input %>%
  full_join(df_ASCVD_Framingham_SCORE2_input %>%
      select(
        PatientID,
        race_ascvd,
        mean_Total_Cholesterol_mg_dl,
        mean_HDL_mg_dl,
        Risk.region,
        SmokingStatusQRISK3
      ) %>%
      rename(SmokerASCVDFraminghamSCORE2 = SmokingStatusQRISK3),   # rename BEFORE join
    by = "PatientID"
  ) %>%
  full_join(
    df_all_risk_scores,
    by = "PatientID"
  )

# exclude PatientID (as not a biological predictor) and keep only columns with ≤50% missing
# Keep ID separately so we can add it back later if needed
patient_ids <- df_joined_for_imputation$PatientID

df_for_mf <- df_joined_for_imputation %>%
  select(-PatientID) %>%   # remove PatientID
  select(where(~ mean(is.na(.x)) <= 0.40)) %>%   # keep only columns with ≤ 40% missingness
  mutate(across(where(is.character), as.factor)) %>%
  as.data.frame()

# set up parallel processing
n_cores <- detectCores() - 1 # Detect number of available CPU cores (minus 1 to keep system responsive)

cl <- makeCluster(n_cores) # Create a cluster of workers, missForest will distribute computations across these cores

registerDoParallel(cl) # Register this cluster so foreach/doParallel can use it

# Run missForest imputation
set.seed(1234) # Ensures reproducibility → random forests are stochastic

mf_res <- missForest(
  df_for_mf,    # the cleaned dataset for imputation
  maxiter     = 10,           # number of imputation iterations → default to start
  ntree       = 100,          # number of trees per forest → default to start
  parallelize = "variables"   # tells missForest to parallelise across variables
)

stopCluster(cl) # Stop the parallel cluster after use, otherwise it will keep running outside od missForest

mf_res$OOBerror

# join the real data with the imputed values
df_imputed <- df_joined_for_imputation %>%
  mutate(across(
    all_of(colnames(df_for_mf)),
    ~ mf_res$ximp[, cur_column()]
  )) %>%
  mutate(PatientID = patient_ids)

# plot the distribution of the imputed Risk Scores
## create function for plotting
plot_risk_distribution <- function(data, var, title = NULL, n_bins = 20) {
  
  v <- data[[var]]   # vector of the variable (like SCORE2$SCORE2_score)

  # stats for subtitle & line
  var_mean  <- mean(v, na.rm = TRUE)
  range_diff <- diff(range(v, na.rm = TRUE))
  n         <- nrow(data)
  
  ggplot(data, aes(x = .data[[var]])) +
    geom_histogram(
      bins  = n_bins,
      fill  = "lightblue",
      color = "white"
    ) +
    labs(
      title = ifelse(
        is.null(title),
        paste("Distribution of", var, "10-year Risk"),
        title
      ),
      subtitle = paste(
        "Mean =", round(var_mean, 2), "%",
        "| N =", n
      ),
      x = "10-year cardiovascular risk (%)",
      y = "Number of participants"
    ) +
    geom_density(
      aes(
        y = after_stat(density) *
          n * range_diff / n_bins
      ),
      color = "darkblue",
      size  = 1.2
    ) +
    geom_vline(
      xintercept = var_mean,
      color      = "red",
      linetype   = "dashed",
      size       = 1
    )
}

# plot the several risk scores
plot_risk_distribution(df_imputed, "QRISK3_2017", 
                       title = "Distribution of QRISK3 (After Imputation)")

plot_risk_distribution(df_imputed, "SCORE2_score", 
                       title = "Distribution of SCORE2 (After Imputation)")

plot_risk_distribution(df_imputed, "ascvd_10y",
                       title = "Distribution of ASCVD 10-year Risk (After Imputation)")

plot_risk_distribution(df_imputed, "frs_10y",
                       title = "Distribution of Framingham 10-year Risk (After Imputation)")

plot_risk_distribution(df_imputed, "mean_risk",
  title = "Distribution of Composite Mean Risk (After Imputation)")


