### QRISK3 and lipidomics data
### Author: Luisa Delius

# Load packages
library(tidyverse)
library(broom)
library(patchwork)
library(QRISK3)
library(readxl)

# Set working directory
wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

# 2 data sets
## File 1: "CoDiet_lipidomics_2025_04_09.csv" contains intact lipid species measured by LC-MS/MS, each defined by two acyl or alkyl chains (e.g. 16:0/20:3), representing complex lipid composition.
## file 2: "CoDiet_lipidomics_dbs_rbc_2025_09_24.csv" contains individual fatty acids measured by GC-MS in red blood cells and dried blood spots, representing fatty acid composition.

# 1. Upload the data/excel sheet
df_intact_lipid_species <- read_csv("CoDiet_lipidomics_2025_04_09.csv")
df_individual_fatty_acids <- read_csv("CoDiet_lipidomics_dbs_rbc_2025_09_24.csv")
QRISK3_sample_ID <- readRDS("QRISK3_sample_ID.RDS") ### this has to be exchanged if i change anything of the QRISK input data of courese!!!
df_statins_supplements <- read_excel("Statins_Supplements.xlsx")

# 2. Prepare both data sets for Regression Modelling
df_lipidomics_mean <- df_intact_lipid_species %>%
  rename(Sample_ID = patient) %>%   # rename ID column
  mutate(Sample_ID = str_replace_all(Sample_ID, "-", "_")) %>%   # replace '-' with '_' in IDs
  rename_with(~ str_replace_all(.x, "[-/:; ]", "_")) %>% # cleans the column names for the later formula
  group_by(Sample_ID) %>%   # average across duplicate Sample_IDs
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")  # <- drops grouping
  
# remove outliers 1at and 99th percentile for all numeric columns
df_lipidomics_trimmed <- df_lipidomics_mean %>%
  mutate(
    across(
      where(is.numeric),
      \(x) replace(
        x,
        x < quantile(x, 0.01, na.rm = TRUE) |
          x > quantile(x, 0.99, na.rm = TRUE),
        NA_real_
      )
    )
  )

removed_outliers_lipidomics <- df_lipidomics_mean %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to  = "variable",
    values_to = "original_value"
  ) %>%
  group_by(variable) %>%  # thresholds are per lipid
  filter(
    original_value < quantile(original_value, 0.01, na.rm = TRUE) |
      original_value > quantile(original_value, 0.99, na.rm = TRUE)
  ) %>%
  ungroup()

df_lipidomics_scaled_and_QRISK3 <- df_lipidomics_trimmed %>%
  mutate(across(where(is.double), ~ scale(.x), .names = "z_{col}")) %>%   # create z-scaled versions of numeric predictors (keep QRISK unscaled)
  inner_join(QRISK3_sample_ID, by = "Sample_ID")

######################## Fatty acids ###############################################

# 2. Prepare fatty acids data for Regression Modelling ---------------------

# 2.1 Clean & average duplicates
df_fatty_acids_mean <- df_individual_fatty_acids %>%
  rename(Sample_ID = patient) %>%   # rename ID column
  mutate(Sample_ID = str_replace_all(Sample_ID, "-", "_")) %>%   # replace '-' with '_' in IDs
  rename_with(~ str_replace_all(.x, "[-/:; ]", "_")) %>%         # clean column names
  group_by(Sample_ID) %>%                                        # average across duplicate Sample_IDs
  summarise(across(where(is.numeric), mean, na.rm = TRUE),
            .groups = "drop")                                    # drop grouping

# 2.2 Remove outliers at 1st and 99th percentile for all numeric columns
df_fatty_acids_trimmed <- df_fatty_acids_mean %>%
  mutate(
    across(
      where(is.numeric),
      \(x) replace(
        x,
        x < quantile(x, 0.01, na.rm = TRUE) |
          x > quantile(x, 0.99, na.rm = TRUE),
        NA_real_
      )
    )
  )

# 2.3 Save removed outliers for inspection
removed_outliers_fatty_acids <- df_fatty_acids_mean %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to  = "variable",
    values_to = "original_value"
  ) %>%
  group_by(variable) %>%  # thresholds are per fatty acid
  filter(
    original_value < quantile(original_value, 0.01, na.rm = TRUE) |
      original_value > quantile(original_value, 0.99, na.rm = TRUE)
  ) %>%
  ungroup()

# 2.4 Z-scale AFTER trimming, then join QRISK3
df_fatty_acids_scaled_and_QRISK3 <- df_fatty_acids_trimmed %>%
  mutate(across(where(is.numeric), ~ scale(.x), .names = "z_{col}")) %>%   # z-scaled versions
  inner_join(QRISK3_sample_ID, by = "Sample_ID")


# 3.Running a Generalized Linear Model each ---------------- #####
# GLM with Gamma distribution and log link
# For GLMs don’t test normality, as the residuals are not supposed to be normal

## Create predictor df
# Use all z-scored features as predictors for both dfs
lipid_predictors <- df_lipidomics_scaled_and_QRISK3 %>%
  select(starts_with("z_")) %>%
  names()
lipid_predictors

fatty_acid_predictors <- df_fatty_acids_scaled_and_QRISK3 %>%
  select(starts_with("z_")) %>%
  names()
fatty_acid_predictors


## Running the glm
# working with lipid predictors
glm_output_lipids <- map_dfr(lipid_predictors, function(var) { # loops over all the predictors, fits a seperate model each and bins together
  df_pair <- df_lipidomics_scaled_and_QRISK3 %>%
    select(QRISK3_risk, all_of(var)) %>%
    drop_na()
  
  model <- glm(
    as.formula(paste("QRISK3_risk ~", var)),
    data   = df_pair,
    family = Gamma(link = "log") # chose Gamma distribution given the positive, continuous, right-skewed data.
  ) # why link = log? Gamma distribution only makes sense for positive means, so by log link the predicted means are all positive. Each 1-unit increase in X increases the expected QRISK3 by e^b1.
  
  # Return tidy estimates; exponentiate to get multiplicative effects (ratios)
  tidy(model, conf.int = TRUE, exponentiate = TRUE) %>% # exponentiates b1 (logarithm of the multiplicative effect) --> easier to read and interpret
    filter(term != "(Intercept)") %>%
    mutate(
      predictor = var,
      n = nrow(df_pair)
    )
}) %>%
  # p.value here is the Wald test p-value for the GLM coefficient
  mutate(p.adj = p.adjust(p.value, method = "BH"),
         significant = if_else(conf.low > 1 | conf.high < 1, "Significant", "Not significant")
  ) %>%
  arrange(p.adj)

glm_output_lipids


# Running GLM with fatty acids
glm_output_fatty_acids <- map_dfr(fatty_acid_predictors, function(var) { # loops over all the predictors, fits a seperate model each and bins together
  df_pair <- df_fatty_acids_scaled_and_QRISK3 %>%
    select(QRISK3_risk, all_of(var)) %>%
    drop_na()
  
  model <- glm(
    as.formula(paste("QRISK3_risk ~", var)),
    data   = df_pair,
    family = Gamma(link = "log") # chose Gamma distribution given the positive, continuous, right-skewed data.
  ) # why link = log? Gamma distribution only makes sense for positive means, so by log link the predicted means are all positive. Each 1-unit increase in X increases the expected QRISK3 by e^b1.
  
  # Return tidy estimates; exponentiate to get multiplicative effects (ratios)
  tidy(model, conf.int = TRUE, exponentiate = TRUE) %>% # exponentiates b1 (logarithm of the multiplicative effect) --> easier to read and interpret
    filter(term != "(Intercept)") %>%
    mutate(
      predictor = var,
      n = nrow(df_pair)
    )
}) %>%
  # p.value here is the Wald test p-value for the GLM coefficient
  mutate(p.adj = p.adjust(p.value, method = "BH"),
         significant = if_else(conf.low > 1 | conf.high < 1, "Significant", "Not significant")
  ) %>%
  arrange(p.adj)

glm_output_fatty_acids


# 4. Plotting the GLM results ---------------- #####
## 4.1 plotting lipids
ord_lipids <- glm_output_lipids %>%
  arrange(estimate) %>%
  pull(predictor)

ggplot(glm_output_lipids,
       aes(x = estimate, y = factor(predictor, levels = ord_lipids), colour = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed") +     # no-effect line is 1 (ratio scale)
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point() +
  scale_colour_manual(values = c("Significant" = "firebrick3",
                                 "Not significant" = "black")) +
  labs(
    x = "Ratio of expected QRISK3 (exp(beta))",
    y = "Predictor",
    title = "GLM (95% CIs) with lipids"
  ) +
  theme(axis.text.y = element_text(size = 6))


## 4.2 plotting fatty acids
ord_fatty_acids <- glm_output_fatty_acids %>%
  arrange(estimate) %>%
  pull(predictor)

ggplot(glm_output_fatty_acids,
       aes(x = estimate, y = factor(predictor, levels = ord_fatty_acids), colour = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed") +     # no-effect line is 1 (ratio scale)
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point() +
  scale_colour_manual(values = c("Significant" = "firebrick3",
                                 "Not significant" = "black")) +
  labs(
    x = "Ratio of expected QRISK3 (exp(beta))",
    y = "Predictor",
    title = "Individual fatty acids"
  ) +
  theme(axis.text.y = element_text(size = 6))

## 4.3 individuals plots for fatty acids
make_plot_glm <- function(var) {
  d <- df_fatty_acids_scaled_and_QRISK3 %>%
    select(QRISK3_risk, all_of(var)) %>%
    drop_na()
  
  gg <- ggplot(d, aes(x = .data[[var]], y = QRISK3_risk)) +
    labs(x = var, y = "QRISK3_risk")
  
  if (is.factor(d[[var]])) {
    gg + geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.5, size = 1.4)
  } else {
    gg + geom_point(size = 1.0, alpha = 0.7) +
      geom_smooth(method = "glm",
                  method.args = list(family = Gamma(link = "log")),
                  se = TRUE)
  }
}

# build all plots
plots_glm <- purrr::map(fatty_acid_predictors, make_plot_glm)

# show 6 per page and print all pages:
lapply(
  split(plots_glm, ceiling(seq_along(plots_glm) / 12)),
  \(p) print(wrap_plots(p, ncol = 4, nrow = 3))
)

# 4.4 individual plots for complete lipids 
make_plot_glm <- function(var) {
  d <- df_lipidomics_scaled_and_QRISK3 %>%
    select(QRISK3_risk, all_of(var)) %>%
    drop_na()
  
  gg <- ggplot(d, aes(x = .data[[var]], y = QRISK3_risk)) +
    labs(x = var, y = "QRISK3_risk")
  
  if (is.factor(d[[var]])) {
    gg + geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.5, size = 1.4)
  } else {
    gg + geom_point(size = 1.0, alpha = 0.7) +
      geom_smooth(method = "glm",
                  method.args = list(family = Gamma(link = "log")),
                  se = TRUE)
  }
}

# build all plots
plots_glm <- purrr::map(lipid_predictors, make_plot_glm)

# show 6 per page and print all pages:
lapply(
  split(plots_glm, ceiling(seq_along(plots_glm) / 12)),
  \(p) print(wrap_plots(p, ncol = 4, nrow = 3))
)

# 5. GLM with fixed effects Statins and Supplements --------#################
## prepare the Statin and Supplemnt data
df_keep_statins_supplements <- df_statins_supplements %>%
  select(Sample_ID, Statins, Supplements) %>%
  mutate(Statins = factor(Statins),
        Supplements = factor(Supplements),
        Sample_ID = str_replace_all(Sample_ID, "-", "_")
  )

## join the columns to the main data
df_fatty_acids_scaled_QRISK3_random <- df_fatty_acids_scaled_and_QRISK3 %>%
  left_join(df_keep_statins_supplements, by = "Sample_ID")

df_lipidomics_scaled_QRISK3_random <- df_lipidomics_scaled_and_QRISK3 %>%
  left_join(df_keep_statins_supplements, by = "Sample_ID")

## save predictors + statins + supplements as RDS
saveRDS(df_fatty_acids_scaled_QRISK3_random, 
        file = file.path(wkdir, "processed_data", "df_fatty_acids_predictor_statin_suppl.rds"))

saveRDS(df_lipidomics_scaled_QRISK3_random,
        file = file.path(wkdir, "processed_data", "df_lipidomics_predictor_statin_suppl.rds"))


## run the glm with fixed effect on FATTY ACIDS
glm_output_fatty_acids <- map_dfr(fatty_acid_predictors, function(var) { # loops over all the predictors, fits a seperate model each and bins together
  df_pair <- df_fatty_acids_scaled_QRISK3_random %>%
    select(QRISK3_risk, all_of(var), Statins, Supplements) %>%
    drop_na()
  
model <- glm(
    as.formula(paste("QRISK3_risk ~", var, "+ Statins + Supplements")),
    data   = df_pair,
    family = Gamma(link = "log") # chose Gamma distribution given the positive, continuous, right-skewed data.
  ) # why link = log? Gamma distribution only makes sense for positive means, so by log link the predicted means are all positive. Each 1-unit increase in X increases the expected QRISK3 by e^b1.
  
  # Return tidy estimates; exponentiate to get multiplicative effects (ratios)
  tidy(model, conf.int = TRUE, exponentiate = TRUE) %>% # exponentiates b1 (logarithm of the multiplicative effect) --> easier to read and interpret
    filter(term == var) %>%
    mutate(
      predictor = var,
      n = nrow(df_pair)
    )
}) %>%
  # p.value here is the Wald test p-value for the GLM coefficient
  mutate(p.adj = p.adjust(p.value, method = "BH"),
         significant = if_else(conf.low > 1 | conf.high < 1, "Significant", "Not significant")
  ) %>%
  arrange(p.adj)

glm_output_fatty_acids

## forest plot FATTY ACIDS with fixed effect
ord_fatty_acids <- glm_output_fatty_acids %>%
  arrange(estimate) %>%
  pull(predictor)

ggplot(glm_output_fatty_acids,
       aes(x = estimate, y = factor(predictor, levels = ord_fatty_acids), colour = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed") +     # no-effect line is 1 (ratio scale)
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point() +
  scale_colour_manual(values = c("Significant" = "firebrick3",
                                 "Not significant" = "black")) +
  labs(
    x = "Ratio of expected QRISK3 (exp(beta))",
    y = "Predictor",
    title = "Fatty acids (adjusted for Statins + Supplements)",
    caption = "Model: QRISK3_risk ~ predictor + Statins + Supplements (Gamma, log link)
              Red = BH-adjusted p < 0.05; Black = BH-adjusted p ≥ 0.05
              Ratios of expected QRISK3 (exp(beta)) with 95% CIs
              Outliers excluded (1st/99th percentile per predictor)
              Sample size (n) varies by predictor due to drop_na() and outlier trimming"
  ) +
  theme(axis.text.y = element_text(size = 6),
        plot.caption.position = "plot",        # ensures caption is placed at bottom right
        plot.caption = element_text(size = 6, hjust = 1)                          # smaller caption tex 
        )

## run the glm with fixed effect on LIPIDS
glm_output_lipids <- map_dfr(lipid_predictors, function(var) { # loops over all the predictors, fits a seperate model each and bins together
  df_pair <- df_lipidomics_scaled_QRISK3_random %>%
    select(QRISK3_risk, all_of(var), Statins, Supplements) %>%
    drop_na()
  
model <- glm(
    as.formula(paste("QRISK3_risk ~", var, "+ Statins + Supplements")),
    data   = df_pair,
    family = Gamma(link = "log") # chose Gamma distribution given the positive, continuous, right-skewed data.
  ) # why link = log? Gamma distribution only makes sense for positive means, so by log link the predicted means are all positive. Each 1-unit increase in X increases the expected QRISK3 by e^b1.
  
  # Return tidy estimates; exponentiate to get multiplicative effects (ratios)
  tidy(model, conf.int = TRUE, exponentiate = TRUE) %>% # exponentiates b1 (logarithm of the multiplicative effect) --> easier to read and interpret
    filter(term == var) %>%
    mutate(
      predictor = var,
      n = nrow(df_pair)
    )
}) %>%
  # p.value here is the Wald test p-value for the GLM coefficient
  mutate(p.adj = p.adjust(p.value, method = "BH"),
         significant = if_else(conf.low > 1 | conf.high < 1, "Significant", "Not significant")
  ) %>%
  arrange(p.adj)

glm_output_lipids

## forest plot LIPIDS with fixed effect
ord_lipids <- glm_output_lipids %>%
  arrange(estimate) %>%
  pull(predictor)

ggplot(glm_output_lipids,
       aes(x = estimate, y = factor(predictor, levels = ord_lipids), colour = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed") +     # no-effect line is 1 (ratio scale)
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point() +
  scale_colour_manual(values = c("Significant" = "firebrick3",
                                 "Not significant" = "black")) +
  labs(
    x = "Ratio of expected QRISK3 (exp(beta))",
    y = "Predictor",
    title = "LIPIDS (adjusted for Statins + Supplements)",
    caption = "Model: QRISK3_risk ~ predictor + Statins + Supplements (Gamma, log link)
              Red = BH-adjusted p < 0.05; Black = BH-adjusted p ≥ 0.05
              Ratios of expected QRISK3 (exp(beta)) with 95% CIs
              Outliers excluded (1st/99th percentile per predictor)
              Sample size (n) varies by predictor due to drop_na() and outlier trimming"
  ) +
  theme(axis.text.y = element_text(size = 6),
        plot.caption.position = "plot",        # ensures caption is placed at bottom right
        plot.caption = element_text(size = 6, hjust = 1)                          # smaller caption tex 
  )