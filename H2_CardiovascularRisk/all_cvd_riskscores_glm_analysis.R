### Cardiovascular Risk Scores and Predictor Associations GLM (+random/fixed effects)
### Author: Luisa Delius

# Load packages
library(tidyverse)
library(broom)
library(patchwork)

# Set working directory
wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files/processed_data"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

# Upload the data/excel sheet
df_all_cvd_risk_scores <- readRDS("df_all_risk_scores.rds")
df_fatty_acids_predictors <- readRDS("df_fatty_acids_predictor_statin_suppl.rds")
df_lipidomics_predictors <- readRDS("df_lipidomics_predictor_statin_suppl.rds")
df_risk_factors_predictors <- readRDS("df_risk_factor_predictors.rds")
df_REDcap_demographics_predictors <- readRDS("df_REDcap_demographics_predictor.rds")
df_body_composition_metrics <- readRDS("df_body_composition_metrics.rds")

# Preparation for running GLM
# 1. do a joint df with everything I want as predictors and all the risk scores. also include statins and supplements for lipids and fatty acids.
df_cvd_scores_and_fatty_acids <- df_all_cvd_risk_scores %>%
  rename(Sample_ID = PatientID) %>% select(-SCORE2_strat) %>%
  full_join(df_fatty_acids_predictors %>% select(-QRISK3_2017), by = "Sample_ID")
 
df_cvd_scores_and_lipidomics <- df_all_cvd_risk_scores %>%
  rename(Sample_ID = PatientID) %>% select(-SCORE2_strat) %>%
  full_join(df_lipidomics_predictors %>% select(-QRISK3_2017), by = "Sample_ID") 

df_cvd_scores_and_REDcap <- df_all_cvd_risk_scores %>%
  rename(Sample_ID = PatientID) %>% select(-SCORE2_strat) %>%
  full_join(df_REDcap_demographics_predictors %>% select(-QRISK3_2017), by = "Sample_ID") 

df_cvd_scores_and_risk_factors <- df_all_cvd_risk_scores %>%
  rename(Sample_ID = PatientID) %>% select(-SCORE2_strat) %>%
  full_join(df_risk_factors_predictors %>% select(-QRISK3_2017), by = "Sample_ID") 

df_cvd_scores_and_body_composition <- df_all_cvd_risk_scores %>%
  rename(Sample_ID = PatientID) %>% select(-SCORE2_strat) %>%
  full_join(df_body_composition_metrics %>% select(-QRISK3_2017), by = "Sample_ID")

# 2. Creating the predictors for the model
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
  names()
REDcap_factor_predictors

risk_factor_num_predictors <- df_risk_factors_predictors %>%
  select(starts_with("z_")) %>%
  names()
risk_factor_num_predictors

risk_factor_factor_predictors <- df_risk_factors_predictors %>% # once i do glmm, no need to keep that
  select(where(is.factor)) %>%
  names()
risk_factor_factor_predictors

body_composition_predictors <- df_body_composition_metrics %>%
  select(starts_with("z_")) %>%
  names()
body_composition_predictors

# 3. define outcomes for the model
outcomes <- c(
  "QRISK3_2017",
  "SCORE2_score",
  "ascvd_10y",
  "frs_10y",
  "mean_risk"      # thats the composite one
)

# 4. GLM functions
## 4.1 GLM function for numeric predictors
run_glm_num <- function(data, outcomes, num_predictors) {
  
  map_dfr(num_predictors, function(var) { # looping over each numeric predictor
    map_dfr(outcomes, function(outcome) { # looping over each outcome
      
      df_pair <- data %>%
        select(all_of(outcome), all_of(var)) %>%
        drop_na()
      
      model <- glm(
        as.formula(paste(outcome, "~", var)),
        data   = df_pair,
        family = Gamma(link = "log")
      )
      
      tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(term != "(Intercept)") %>%
        mutate(
          predictor = var,
          level     = NA_character_,
          outcome   = outcome,
          n         = nrow(df_pair),
          significant = if_else(conf.low > 1 | conf.high < 1, "Significant", "Not significant")
        )
    })
  })
}

# running the GLM over numeric predictors
glm_fatty_num  <- run_glm_num(df_cvd_scores_and_fatty_acids, outcomes, fatty_acid_predictors)
glm_lipid_num  <- run_glm_num(df_cvd_scores_and_lipidomics,  outcomes, lipid_predictors)
glm_REDcap_num <- run_glm_num(df_cvd_scores_and_REDcap,      outcomes, REDcap_numeric_predictors)
glm_risk_num   <- run_glm_num(df_cvd_scores_and_risk_factors, outcomes, risk_factor_num_predictors)
glm_body_comp_num <- run_glm_num(df_cvd_scores_and_body_composition, outcomes, body_composition_predictors)

## 4.2 GLM function for factor predictors
run_glm_fac <- function(data, outcomes, factor_predictors) {
  
  map_dfr(factor_predictors, function(var) {
    map_dfr(outcomes, function(outcome) {
      
      df_pair <- data %>%
        select(all_of(outcome), all_of(var)) %>%
        drop_na()
      
      levs <- levels(df_pair[[var]])
      if (length(levs) == 0L) return(tibble())
      
      map_dfr(levs, function(lev) {
        
        df_tmp <- df_pair %>%
          mutate(dummy = if_else(.data[[var]] == lev, 1, 0))
        
        model <- glm(
          as.formula(paste(outcome, "~ dummy")),
          data   = df_tmp,
          family = Gamma(link = "log")
        )
        
        tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
          filter(term == "dummy") %>%
          mutate(
            predictor = var,
            level     = lev,
            term      = paste0(var, " = ", lev),
            outcome   = outcome,
            n         = nrow(df_tmp),
            significant = if_else(conf.low > 1 | conf.high < 1, "Significant", "Not significant")
          )
      })
    })
  })
}

# running the GLM over factor predictors
glm_REDcap_fac <- run_glm_fac(df_cvd_scores_and_REDcap,      outcomes, REDcap_factor_predictors)
glm_risk_fac   <- run_glm_fac(df_cvd_scores_and_risk_factors, outcomes, risk_factor_factor_predictors)


# 5. Forest plot of the GLM results
## 5.1 Preparing a Forest plot function
make_forest_one_outcome <- function(df, outcome_name, title_prefix) {
  
  df_out <- df %>%
    filter(outcome == outcome_name) %>%
    mutate(term_plot = if_else(is.na(level), predictor, term))
  
  ord <- df_out %>% # order by estimate
    arrange(estimate) %>%
    pull(term_plot) %>%
    unique()
  
  p <- ggplot(df_out,
    aes(
      x      = estimate,
      y      = factor(term_plot, levels = ord),
      colour = significant
    )
  ) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_errorbar(
      aes(xmin = conf.low, xmax = conf.high),
      height      = 0.2
    ) +
    geom_point() +
    scale_colour_manual(
      values = c(
        "Significant"    = "firebrick3",
        "Not significant" = "black"
      )
    ) +
    labs(
      x     = "Ratio of expected risk score (exp(beta))",
      y     = "Predictor / level",
      title = paste0(title_prefix, " – ", outcome_name)
    ) +
    theme(axis.text.y = element_text(size = 6))
}

make_all_forest_plots <- function(df, title_prefix, outcomes) {
  setNames(
    map(
      outcomes,
      ~ make_forest_one_outcome(df, .x, title_prefix)
    ),
    outcomes
  )
}

## 5.2 Running the Forest Plot functions
lipid_plots <- make_all_forest_plots(
  df           = glm_lipid_num,
  title_prefix = "Lipids",
  outcomes     = outcomes)

fatty_plots <- make_all_forest_plots(
  df           = glm_fatty_num,
  title_prefix = "Fatty acids",
  outcomes     = outcomes)

REDcap_all <- bind_rows(glm_REDcap_num, glm_REDcap_fac)
REDcap_plots <- make_all_forest_plots(
  df           = REDcap_all,
  title_prefix = "REDCap demographics",
  outcomes     = outcomes)

risk_fact_all <- bind_rows(glm_risk_num, 
                            glm_risk_fac %>% filter(term != "Age.Risk = III"))
risk_fact_plots <- make_all_forest_plots(
  df           = risk_fact_all,
  title_prefix = "Clinical risk factors",
  outcomes     = outcomes)

body_comp_plots <- make_all_forest_plots(
  df           = glm_body_comp_num,
  title_prefix = "Body composition metrics",
  outcomes     = outcomes)

## 5.3 Plotting
lipid_plots[["QRISK3_2017"]]
lipid_plots[["SCORE2_score"]]
lipid_plots[["ascvd_10y"]]
lipid_plots[["frs_10y"]]
lipid_plots[["mean_risk"]]

fatty_plots[["QRISK3_2017"]]
fatty_plots[["SCORE2_score"]]

REDcap_plots[["QRISK3_2017"]]
REDcap_plots[["SCORE2_score"]]

risk_fact_plots[["QRISK3_2017"]]
risk_fact_plots[["SCORE2_score"]]

body_comp_plots[["SCORE2_score"]]

# 6. Combined table about significant ones
all_results <- bind_rows( # combine all results
  glm_lipid_num  %>% mutate(predictor_set = "Lipids"),
  glm_fatty_num  %>% mutate(predictor_set = "Fatty acids"),
  glm_REDcap_num %>% mutate(predictor_set = "REDCap"),
  glm_REDcap_fac %>% mutate(predictor_set = "REDCap"),
  glm_risk_num   %>% mutate(predictor_set = "Risk factors"),
  glm_risk_fac   %>% mutate(predictor_set = "Risk factors"),
  glm_body_comp_num %>% mutate(predictor_set = "Body composition")
) %>%
  mutate(
    term_plot = if_else(is.na(level), predictor, term)  
  )

# I only want to have significant predictors for the following
sig_only <- all_results %>% 
  filter(significant == "Significant") %>%
  mutate(
    direction = case_when(
      estimate < 1 ~ "Lower risk",
      estimate > 1 ~ "Higher risk",
      TRUE         ~ NA_character_
    )
  )

# define y-axis order: grouped by predictor set and internaly by amount of significance
term_levels <- sig_only %>%
  count(predictor_set, term_plot) %>%
  arrange(predictor_set, desc(n), term_plot) %>%
  pull(term_plot)

# plot
ggplot(sig_only,
  aes(
    x    = outcome,
    y    = factor(term_plot, levels = term_levels),
    fill = direction
  )
) +
  geom_tile(
    colour = "white"
  ) +
  scale_fill_manual(values = c(
    "Lower risk"  = "#3cb371",  
    "Higher risk" = "#d73027" 
  )) +
  labs(
    x = "Cardiovascular risk score",
    y = "Predictor",
    fill = "",
    title = "Significant predictors across cardiovascular risk scores"
  ) +
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y      = element_text(size = 6),
    panel.background = element_rect(fill = "white")
    )

# 7. Kendall correlation of the significant result of the risk scores
# Direction of the predictors
wide_dir <- all_results %>%
  mutate( # encode direction per predictor oucome
    dir_val = case_when(
      significant == "Significant" & estimate > 1 ~  1L,
      significant == "Significant" & estimate < 1 ~ -1L,
      TRUE                                       ~  0L # 0 for not significant
    )
  ) %>%
  select(term_plot, outcome, dir_val) %>%
  pivot_wider(
    names_from  = outcome, # column = risk score
    values_from = dir_val # value = direction of the predictor
  )

# Matrix of risk scores (columns) across predictors (rows)
risk_mat <- wide_dir %>%
  select(all_of(outcomes))

# Kendall correlation matrix
cor_kendall <- cor(risk_mat, method = "kendall", use = "pairwise.complete.obs")

# Plot Heatmpa style
as.data.frame(cor_kendall) %>% # turn correlation matrix → data frame
  rownames_to_column("var1") %>% # rownames → column
  pivot_longer(-var1, names_to = "var2", values_to = "tau") %>% # take all columns except var1
  ggplot(aes(var1, var2, fill = tau)) +
  geom_tile() +
  geom_text(aes(label = round(tau, 2))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1)) +
  coord_equal() + # makes x and y scales the same
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Kendall correlation between risk scores\n(based on significant predictor directions)",
    x = "",
    y = ""
  )

# Plot Network graph style
edges <- as.data.frame(cor_kendall) %>% # turn Kendall Matrix into an edge list
  rownames_to_column("from") %>%
  pivot_longer(-from, names_to = "to", values_to = "tau")

# put nodes on a circle
nodes <- data.frame(
  name = scores,
  angle = seq(0, 2*pi, length.out = n+1)[-1]
) %>%
  mutate(
    x = cos(angle),
    y = sin(angle)
  )

# join coordinates for drawing lines
edges <- edges %>%
  left_join(nodes, by = c("from" = "name")) %>%
  rename(x1 = x, y1 = y) %>%
  left_join(nodes, by = c("to" = "name")) %>%
  rename(x2 = x, y2 = y)

ggplot() +
  geom_segment(
    data = edges,
    aes(x = x1, y = y1, xend = x2, yend = y2,
        colour = tau, size = abs(tau))
  ) +
  geom_point(
    data = nodes,
    aes(x = x, y = y),
    size = 5
  ) +
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = name),
    vjust = -1
  ) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                         limits = c(-1, 1)) +
  scale_size(range = c(0.5, 3)) +
  coord_equal() +
  theme_void() +
  ggtitle("Network-style plot of Kendall correlations between risk scores")
