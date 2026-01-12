### Cardiovascular Risk Scores and Predictor Associations GLMM with 
### fixed effects: statins, supplements, sex and 
### random effect: country

### Author: Luisa Delius

# Load packages
library(tidyverse)
library(broom)
library(patchwork)
library(readxl)
library(lme4)
library(broom.mixed)

set.seed(42)

# Set working directory
wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files/processed_data"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

# 1. Upload the data/excel sheet
df_fatty_acids_predictors <- readRDS("df_fatty_acids_predictor_statin_suppl.rds") %>% arrange(Sample_ID)
df_lipidomics_predictors <- readRDS("df_lipidomics_predictor_statin_suppl.rds") %>% arrange(Sample_ID)
df_REDcap_demographics_predictors <- readRDS("df_REDcap_demographics_predictor.rds") %>% arrange(Sample_ID)
df_body_composition_metrics <- readRDS("df_body_composition_metrics.rds") %>% arrange(Sample_ID)
df_urine_nmr_data <- readRDS("df_urine_NMR_data.rds") %>% arrange(Sample_ID)

df_risk_factors_predictors <- readRDS("df_risk_factor_predictors.rds") %>% 
  select(-Age.Risk, -stress_resilience_status, -stress_index_status) %>% # remove the factors we dont want to include
  arrange(Sample_ID)

df_all_cvd_risk_scores <- readRDS("df_all_risk_scores.rds") %>% 
  rename(Sample_ID = PatientID) %>%
  select(-SCORE2_strat) %>%
  arrange(Sample_ID)

# 2. Preparation for running GLM
df_sex_country_statins_supplements <- df_risk_factors_predictors %>% # joint df with the fixed and random effects I will use
  select(Sample_ID, Gender, Country) %>%
  full_join(df_lipidomics_predictors %>% select(Sample_ID, Statins, Supplements), by = "Sample_ID")

df_cvd_scores_and_fatty_acids <- df_all_cvd_risk_scores %>%
  full_join(df_fatty_acids_predictors %>% select(-QRISK3_risk), by = "Sample_ID") %>%
  full_join(df_sex_country_statins_supplements %>% select(Sample_ID, Gender, Country), by = "Sample_ID")

df_cvd_scores_and_lipidomics <- df_all_cvd_risk_scores %>%
  full_join(df_lipidomics_predictors %>% select(-QRISK3_risk), by = "Sample_ID") %>%
  full_join(df_sex_country_statins_supplements %>% select(Sample_ID, Gender, Country), by = "Sample_ID")

df_cvd_scores_and_REDcap <- df_all_cvd_risk_scores %>%
  full_join(df_REDcap_demographics_predictors %>% select(-QRISK3_risk, -recruitment_site), by = "Sample_ID") %>%
  full_join(df_sex_country_statins_supplements, by = "Sample_ID")

df_cvd_scores_and_risk_factors <- df_all_cvd_risk_scores %>%
  full_join(df_risk_factors_predictors %>% select(-QRISK3_risk), by = "Sample_ID") %>%
  full_join(df_sex_country_statins_supplements %>% select(Sample_ID, Statins, Supplements), by = "Sample_ID")

df_cvd_scores_and_body_composition <- df_all_cvd_risk_scores %>%
  full_join(df_body_composition_metrics %>% select(-QRISK3_risk), by = "Sample_ID") %>%
  full_join(df_sex_country_statins_supplements, by = "Sample_ID")

df_cvd_scores_and_urine_nmr <- df_all_cvd_risk_scores %>%
  full_join(df_urine_nmr_data, by = "Sample_ID") %>%
  full_join(df_sex_country_statins_supplements, by = "Sample_ID")


# 3. Creating the predictors for the model
lipid_predictors <- df_lipidomics_predictors %>%
  select(starts_with("z_")) %>%
  names()
lipid_predictors

fatty_acid_predictors <- df_fatty_acids_predictors %>%
  select(starts_with("z_")) %>%
  names()
fatty_acid_predictors

REDcap_numeric_predictors <- df_REDcap_demographics_predictors %>%
  select(
    starts_with("z_")) %>%
  names()
REDcap_numeric_predictors

REDcap_factor_predictors <- df_REDcap_demographics_predictors %>%
  select(where(is.factor)) %>%
  select(-recruitment_site) %>%
  names()
REDcap_factor_predictors

risk_factor_predictors <- df_risk_factors_predictors %>%
  select(starts_with("z_")) %>%
  names()
risk_factor_predictors

body_composition_predictors <- df_body_composition_metrics %>%
  select(starts_with("z_")) %>%
  names()
body_composition_predictors

urine_nmr_predictors <- df_urine_nmr_data %>%
  select(starts_with("z_")) %>%
  names()
urine_nmr_predictors

# 4. define outcomes for the model
outcomes <- c(
  "QRISK3_risk",
  "SCORE2_score",
  "ascvd_10y",
  "frs_10y",
  "mean_risk"      # thats the composite one
)

# 5. Model function (GLM with fixed effect Country, Statins, supplements, sex) for all predictor-outcome combinations
# 5.1 GLM for numeric predictor
run_glm_num <- function(data, outcomes, num_predictors, fixed_effects) {
  
  fixed_part <- paste(fixed_effects, collapse = " + ")
  
  map_dfr(num_predictors, function(var) {
    map_dfr(outcomes, function(outcome) {
      
      df_pair <- data %>%
        select(all_of(outcome), all_of(var), all_of(fixed_effects)) %>%
        drop_na()
      
      # Skip if no data
      if (nrow(df_pair) == 0) return(tibble())
      
      # Skip if any fixed effect has < 2 levels
      for (fe in fixed_effects) {
        if (n_distinct(df_pair[[fe]]) < 2) {
          return(tibble())
        }
      }
      
      model <- glm(
        as.formula(paste(outcome, "~", var, "+", fixed_part)),
        data = df_pair,
        family = Gamma(link = "log")
      )
      
      tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(term == var) %>%
        mutate(
          predictor = var,
          level = NA_character_,
          outcome = outcome,
          n = nrow(df_pair)
        )
    })
  })
}

# 5.2 GLM for factor predictors
run_glm_fac <- function(data, outcomes, factor_predictors, fixed_effects) {
  
  fixed_part <- paste(fixed_effects, collapse = " + ")
  
  map_dfr(factor_predictors, function(var) {
    map_dfr(outcomes, function(outcome) {
      
      df_pair <- data %>%
        select(all_of(outcome), all_of(var), all_of(fixed_effects)) %>%
        drop_na()
      
      levs <- levels(df_pair[[var]])
      if (length(levs) == 0L) return(tibble())
      
      map_dfr(levs, function(lev) {
        
        df_tmp <- df_pair %>%
          mutate(dummy = if_else(.data[[var]] == lev, 1, 0))
        
        # Skip if no data
        if (nrow(df_tmp) == 0) return(tibble())
        
        model <- glm(
          as.formula(paste(outcome, "~ dummy +", fixed_part)),
          data = df_tmp,
          family = Gamma(link = "log")
        )
        
        tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
          filter(term == "dummy") %>%
          mutate(
            predictor = var,
            level = lev,
            term = paste0(var, " = ", lev),
            outcome = outcome,
            n = nrow(df_tmp)
          )
      })
    })
  })
}

# 5.3 Run the GLM
# Define fixed effects once
fixed_effects <- c("Statins", "Supplements", "Gender", "Country")

# Run GLM models - Numeric predictors
glm_lipid_num     <- run_glm_num(df_cvd_scores_and_lipidomics, outcomes, lipid_predictors, fixed_effects)
glm_fatty_num     <- run_glm_num(df_cvd_scores_and_fatty_acids, outcomes, fatty_acid_predictors, fixed_effects)
glm_urine_num     <- run_glm_num(df_cvd_scores_and_urine_nmr, outcomes, urine_nmr_predictors, fixed_effects)
glm_REDcap_num    <- run_glm_num(df_cvd_scores_and_REDcap, outcomes, REDcap_numeric_predictors, fixed_effects)
glm_risk_num      <- run_glm_num(df_cvd_scores_and_risk_factors, outcomes, risk_factor_predictors, fixed_effects)
glm_body_comp_num <- run_glm_num(df_cvd_scores_and_body_composition, outcomes, body_composition_predictors, fixed_effects)

# Factor predictors
glm_REDcap_fac <- run_glm_fac(df_cvd_scores_and_REDcap, outcomes, REDcap_factor_predictors, fixed_effects)

# Combined table
all_results_glm <- bind_rows(
  glm_lipid_num     %>% mutate(predictor_set = "Lipids"),
  glm_fatty_num     %>% mutate(predictor_set = "Fatty acids"),
  glm_urine_num     %>% mutate(predictor_set = "Urine NMR"),
  glm_REDcap_num    %>% mutate(predictor_set = "REDCap numeric"),
  glm_REDcap_fac    %>% mutate(predictor_set = "REDCap factors"),
  glm_risk_num      %>% mutate(predictor_set = "Risk factors"),
  glm_body_comp_num %>% mutate(predictor_set = "Body composition")
) %>%
  mutate(
    term_plot = if_else(is.na(level), predictor, term),
    p.adjusted = p.adjust(p.value, method = "BH"),
    significant = if_else(p.adjusted < 0.05, "Significant", "Not significant")
  )

# Summary of the analysis
print("Total associations tested:")
print(nrow(all_results_glm))

print("Significant (p.adjusted < 0.05):")
print(sum(all_results_glm$significant == "Significant", na.rm = TRUE))


# Function to create clean supplementary table
create_supp_table <- function(data, predictor_set_name) {
  data %>%
    filter(predictor_set == predictor_set_name) %>%
    select(term_plot, outcome, estimate, conf.low, conf.high, 
           p.value, p.adjusted, n, significant) %>%
    arrange(p.adjusted, outcome)
}

# Create tables for each predictor set
supp_table_lipids <- create_supp_table(all_results_glm, "Lipids")
supp_table_fatty_acids <- create_supp_table(all_results_glm, "Fatty acids")
supp_table_urine_nmr <- create_supp_table(all_results_glm, "Urine NMR")
supp_table_redcap_num <- create_supp_table(all_results_glm, "REDCap numeric")
supp_table_redcap_fac <- create_supp_table(all_results_glm, "REDCap factors")
supp_table_risk_factors <- create_supp_table(all_results_glm, "Risk factors")
supp_table_body_comp <- create_supp_table(all_results_glm, "Body composition")

View(supp_table_lipids)

# 6.1 plot FDR-corrected significant results
# Prepare data for heatmap - REMOVE NAs
heatmap_data <- all_results_glm %>%
  filter(!is.na(estimate)) %>%  # REMOVE NAs
  mutate(
    log_effect = log(estimate),
    is_significant = p.adjusted < 0.05
  )

# Filter to only predictors with at least one significant result
predictors_with_sig <- heatmap_data %>%
  filter(is_significant) %>%
  pull(term_plot) %>%
  unique()

heatmap_data_filtered <- heatmap_data %>%
  filter(term_plot %in% predictors_with_sig)

# Order predictors by significance
term_order <- heatmap_data_filtered %>%
  filter(is_significant) %>%
  group_by(term_plot) %>%
  summarise(
    n_sig = n(),
    mean_effect = mean(abs(log_effect))
  ) %>%
  arrange(desc(n_sig), desc(mean_effect)) %>%
  pull(term_plot)

# Create plot
ggplot(heatmap_data_filtered,
       aes(x = outcome, 
           y = factor(term_plot, levels = rev(term_order)), 
           fill = log_effect)) +
  geom_tile(colour = "white", size = 0.5) +
  geom_text(
    data = heatmap_data_filtered %>% filter(is_significant),
    aes(label = "*"),
    size = 4,
    vjust = 0.75,
    colour = "black"
  ) +
  facet_grid(predictor_set ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#0072B2",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    name = "log(Rate Ratio)",
    breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
    labels = c("-0.8", "-0.4", "0", "+0.4", "+0.8")
  ) +
  labs(
    x = "Cardiovascular risk score",
    y = "Predictor",
    title = "GLM: Associations between biomarkers and CVD risk scores",
    subtitle = "* indicates FDR-adjusted p < 0.05",
    caption = "Model: outcome ~ predictor + Statins + Supplemets + Sex + Country (Gamma GLM, log link)\nMultiple testing correction: Benjamini-Hochberg"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 8, hjust = 0),
    strip.text.y = element_text(size = 9, face = "bold", angle = 0),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave("heatmap_significant_predictors.png", width = 12, height = 16, dpi = 300)

# 6.2 Plot predictors significant at least 3x for p.value < 0.05 (& p.adjusted > 0.05)
# Get FDR-significant predictors to exclude
fdr_significant_predictors <- all_results_glm %>%
  filter(p.adjusted < 0.05) %>%
  pull(predictor) %>%
  unique()

# Filter: nominal significant, NOT FDR-significant, replicated ≥3 scores
exploratory_replicated <- all_results_glm %>%
  filter(p.value < 0.05 & p.adjusted > 0.05,
         !predictor %in% fdr_significant_predictors) %>%
  group_by(term_plot, predictor, predictor_set) %>%
  filter(n() >= 3) %>%  # Keep only terms significant in ≥3 scores
  ungroup()

# Prepare heatmap data
heatmap_data <- exploratory_replicated %>%
  filter(!is.na(estimate)) %>%
  mutate(log_effect = log(estimate))

# Order by mean absolute effect size
term_order <- heatmap_data %>%
  group_by(term_plot, predictor_set) %>%
  summarise(mean_abs_effect = mean(abs(log_effect)), .groups = "drop") %>%
  arrange(predictor_set, desc(mean_abs_effect)) %>%
  pull(term_plot)

# Create heatmap
ggplot(heatmap_data,
       aes(x = outcome, 
           y = factor(term_plot, levels = rev(term_order)), 
           fill = log_effect)) +
  geom_tile(colour = "white", size = 0.5) +
  facet_grid(predictor_set ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#0072B2",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    name = "log(Rate Ratio)",
    breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
    labels = c("-0.8", "-0.4", "0", "+0.4", "+0.8")
  ) +
  labs(
    x = "Cardiovascular risk score",
    y = "Predictor",
    title = "GLM: Exploratory Associations",
    subtitle = "Predictors with nominal significance (p < 0.05, FDR p > 0.05) in ≥3 CVD risk scores",
    caption = "Model: outcome ~ predictor + Statins + Supplements + Sex + Country (Gamma GLM, log link)\nExcludes FDR-significant predictors | Blue = protective, Orange = risk factor"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 8, hjust = 0),
    strip.text.y = element_text(size = 9, face = "bold", angle = 0),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave("exploratory_heatmap_3plus_scores.png", width = 10, height = 8, dpi = 300)