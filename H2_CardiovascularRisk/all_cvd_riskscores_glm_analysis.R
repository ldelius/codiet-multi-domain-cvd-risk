### Cardiovascular Risk Scores and Predictor Associations GLM (+random/fixed effects)
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

# Upload the data/excel sheet
df_all_cvd_risk_scores <- readRDS("df_all_risk_scores.rds") %>% arrange(PatientID)
# df_all_cvd_risk_scores_old <- readRDS("df_all_risk_scores_old.rds")  %>% arrange(PatientID)
df_fatty_acids_predictors <- readRDS("df_fatty_acids_predictor_statin_suppl.rds") %>% arrange(Sample_ID)
df_lipidomics_predictors <- readRDS("df_lipidomics_predictor_statin_suppl.rds") %>% arrange(Sample_ID)
df_risk_factors_predictors <- readRDS("df_risk_factor_predictors.rds") %>% arrange(Sample_ID)
df_REDcap_demographics_predictors <- readRDS("df_REDcap_demographics_predictor.rds") %>% arrange(Sample_ID)
df_body_composition_metrics <- readRDS("df_body_composition_metrics.rds") %>% arrange(Sample_ID)

# Preparation for running GLM
df_cvd_scores_and_fatty_acids <- df_all_cvd_risk_scores %>%
  rename(Sample_ID = PatientID) %>% select(-SCORE2_strat) %>%
  full_join(df_fatty_acids_predictors %>% select(-QRISK3_risk), by = "Sample_ID")
 
df_cvd_scores_and_lipidomics <- df_all_cvd_risk_scores %>%
  rename(Sample_ID = PatientID) %>% select(-SCORE2_strat) %>%
  full_join(df_lipidomics_predictors %>% select(-QRISK3_risk), by = "Sample_ID") 

df_cvd_scores_and_REDcap <- df_all_cvd_risk_scores %>%
  rename(Sample_ID = PatientID) %>% select(-SCORE2_strat) %>%
  full_join(df_REDcap_demographics_predictors %>% select(-QRISK3_risk), by = "Sample_ID") 

df_cvd_scores_and_risk_factors <- df_all_cvd_risk_scores %>%
  rename(Sample_ID = PatientID) %>% select(-SCORE2_strat) %>%
  full_join(df_risk_factors_predictors %>% select(-QRISK3_risk), by = "Sample_ID") 

df_cvd_scores_and_body_composition <- df_all_cvd_risk_scores %>%
  rename(Sample_ID = PatientID) %>% select(-SCORE2_strat) %>%
  full_join(df_body_composition_metrics %>% select(-QRISK3_risk), by = "Sample_ID")

df_sex_country <- df_risk_factors_predictors %>%
  select(Sample_ID, Gender, Country)

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
  select(-recruitment_site) %>%
  names()
REDcap_factor_predictors

risk_factor_num_predictors <- df_risk_factors_predictors %>%
  select(starts_with("z_")) %>%
  names()
risk_factor_num_predictors

risk_factor_factor_predictors <- df_risk_factors_predictors %>% # once i do glmm, no need to keep that
  select(where(is.factor)) %>%
  select(-Gender, -Country) %>% # only for the glmm, delete this for glm
  names()
risk_factor_factor_predictors

body_composition_predictors <- df_body_composition_metrics %>%
  select(starts_with("z_")) %>%
  names()
body_composition_predictors

# 3. define outcomes for the model
outcomes <- c(
  "QRISK3_risk",
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
lipid_plots[["QRISK3_risk"]]
lipid_plots[["SCORE2_score"]]
lipid_plots[["ascvd_10y"]]
lipid_plots[["frs_10y"]]
lipid_plots[["mean_risk"]]

fatty_plots[["QRISK3_risk"]]

REDcap_plots[["QRISK3_risk"]]

risk_fact_plots[["QRISK3_risk"]]

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


# Order: direction first, then count of significance
term_levels <- sig_only %>%
  mutate(direction_order = if_else(estimate > 1, 1, 2)) %>%
  group_by(predictor_set, term_plot, direction_order) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(predictor_set, direction_order, desc(n), term_plot) %>%
  pull(term_plot)

# Plot
ggplot(sig_only,
       aes(x = outcome, y = factor(term_plot, levels = term_levels), fill = direction)
) +
  geom_tile(colour = "white") +
  facet_grid(predictor_set ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Lower risk" = "#3cb371", "Higher risk" = "#d73027")) +
  labs(
    x = "Cardiovascular risk score",
    y = "Predictor",
    fill = "",
    title = "GLM: Significant predictors across cardiovascular risk scores",
    caption = "Model: outcome ~ predictor (Gamma, log link)\nOutliers excluded (1st/99th percentile) for Lipids & Fatty acids only"
  ) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    panel.background = element_rect(fill = "white"),
    plot.caption = element_text(size = 6, hjust = 1),
    strip.text.y = element_text(size = 8, face = "bold", angle = 0)
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



######## --------------- ########## -------------------- ############
# 8. Run GLMM (mixed effect with Gender and Country) for risk_factor, sociodemographic, body composition
# and mixed effectwith statins and supplements for lipidomics and fatty acids

# 8.1 Add Sex and Country to the df | or statins & supplements
df_cvd_scores_and_REDcap <- df_cvd_scores_and_REDcap %>%
  left_join(df_sex_country, by = "Sample_ID") %>%
  select(-recruitment_site)

df_cvd_scores_and_body_composition <- df_cvd_scores_and_body_composition %>%
  left_join(df_sex_country, by = "Sample_ID")


# 8.2 glmm function
# GLMM for numeric predictors
run_glmm_num <- function(data, outcomes, num_predictors, random_effects) {
  
  random_part <- paste0("(1|", random_effects, ")", collapse = " + ")
  
  map_dfr(num_predictors, function(var) {
    map_dfr(outcomes, function(outcome) {
      
      df_pair <- data %>%
        select(all_of(outcome), all_of(var), all_of(random_effects)) %>%
        drop_na()
      
      # Skip if any random effect has < 2 levels
      for (re in random_effects) {
        if (n_distinct(df_pair[[re]]) < 2) {
          message("Skipping ", var, " ~ ", outcome, ": ", re, " has < 2 levels")
          return(tibble())
        }
      }     
      
      model <- lme4::glmer(
        as.formula(paste(outcome, "~", var, "+", random_part)),
        data    = df_pair,
        family  = Gamma(link = "log"),
        control = glmerControl(optimizer = "bobyqa")
      )
      
      broom.mixed::tidy(model, conf.int = TRUE, exponentiate = TRUE, effects = "fixed") %>%
        filter(term != "(Intercept)") %>%
        mutate(
          predictor   = var,
          level       = NA_character_,
          outcome     = outcome,
          n           = nrow(df_pair),
          significant = if_else(conf.low > 1 | conf.high < 1, "Significant", "Not significant")
        )
    })
  })
}

# GLMM for factor predictors (with dummies)
run_glmm_fac <- function(data, outcomes, factor_predictors, random_effects) {
  
  random_part <- paste0("(1|", random_effects, ")", collapse = " + ")
  
  map_dfr(factor_predictors, function(var) {
    map_dfr(outcomes, function(outcome) {
      
      df_pair <- data %>%
        select(all_of(outcome), all_of(var), all_of(random_effects)) %>%
        drop_na()
      
      levs <- levels(df_pair[[var]])
      if (length(levs) == 0L) return(tibble())
      
      map_dfr(levs, function(lev) {
        
        df_tmp <- df_pair %>%
          mutate(dummy = if_else(.data[[var]] == lev, 1, 0))
        
        model <- lme4::glmer(
          as.formula(paste(outcome, "~ dummy +", random_part)),
          data    = df_tmp,
          family  = Gamma(link = "log"),
          control = glmerControl(optimizer = "bobyqa")
        )
        
        broom.mixed::tidy(model, conf.int = TRUE, exponentiate = TRUE, effects = "fixed") %>%
          filter(term == "dummy") %>%
          mutate(
            predictor   = var,
            level       = lev,
            term        = paste0(var, " = ", lev),
            outcome     = outcome,
            n           = nrow(df_tmp),
            significant = if_else(conf.low > 1 | conf.high < 1, "Significant", "Not significant")
          )
      })
    })
  })
}

# 8.3 Run the GLMM model
# Numeric
glmm_REDcap_num    <- run_glmm_num(df_cvd_scores_and_REDcap, outcomes, REDcap_numeric_predictors, c("Gender", "Country"))
glmm_risk_num      <- run_glmm_num(df_cvd_scores_and_risk_factors, outcomes, risk_factor_num_predictors, c("Gender", "Country"))
glmm_body_comp_num <- run_glmm_num(df_cvd_scores_and_body_composition, outcomes, body_composition_predictors, c("Gender", "Country"))
glmm_lipid_num     <- run_glmm_num(df_cvd_scores_and_lipidomics, outcomes, lipid_predictors, c("Statins", "Supplements"))
glmm_fatty_num     <- run_glmm_num(df_cvd_scores_and_fatty_acids, outcomes, fatty_acid_predictors, c("Statins", "Supplements"))

# Factor
glmm_REDcap_fac <- run_glmm_fac(df_cvd_scores_and_REDcap, outcomes, REDcap_factor_predictors, c("Gender", "Country"))
glmm_risk_fac   <- run_glmm_fac(df_cvd_scores_and_risk_factors, outcomes, risk_factor_factor_predictors, c("Gender", "Country"))

# 8.4 Combined plot of significant results
# Combined table for GLMM results
all_results_glmm <- bind_rows(
  glmm_lipid_num     %>% mutate(predictor_set = "Lipids"),
  glmm_fatty_num     %>% mutate(predictor_set = "Fatty acids"),
  glmm_REDcap_num    %>% mutate(predictor_set = "REDCap"),
  glmm_REDcap_fac    %>% mutate(predictor_set = "REDCap"),
  glmm_risk_num      %>% mutate(predictor_set = "Risk factors"),
  glmm_risk_fac      %>% mutate(predictor_set = "Risk factors"),
  glmm_body_comp_num %>% mutate(predictor_set = "Body composition")
) %>%
  mutate(
    term_plot = if_else(is.na(level), predictor, term),
    p.adjusted = p.adjust(p.value, method = "BH"),
    significant_adjusted = if_else(p.adjusted < 0.05, "Significant", "Not significant")
  )

#### Predictors that are significant based on p-adj (BH) ####
# Keep only results significant after BH adjustment
all_results_glmm_sig <- all_results_glmm %>%
  filter(p.adjusted < 0.05)

# Create factor levels for current significant results
term_levels_padj <- sig_padj_glmm %>%
  arrange(predictor_set, term_plot) %>%
  pull(term_plot) %>%
  unique()

# Plot
ggplot(sig_padj_glmm,
       aes(x = outcome, y = factor(term_plot, levels = term_levels_padj), fill = direction)
) +
  geom_tile(colour = "white") +
  facet_grid(predictor_set ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Lower risk" = "#3cb371", "Higher risk" = "#d73027")) +
  labs(
    x = "Cardiovascular risk score",
    y = "Predictor",
    fill = "",
    title = "GLMM: Predictors with BH-adjusted p < 0.05",
    caption = "Model: outcome ~ predictor + (1|random effects) (Gamma, log link)\nRandom effects: Gender & Country (REDCap, Risk factors, Body composition)\nRandom effects: Statins & Supplements (Lipids, Fatty acids)\nOutliers excluded (1st/99th percentile) for Lipids & Fatty acids only"
  ) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    panel.background = element_rect(fill = "white"),
    plot.caption = element_text(size = 6, hjust = 1),
    strip.text.y = element_text(size = 8, face = "bold", angle = 0)
  )

# Significant predictors only
# sig_only_glmm <- all_results_glmm %>% 
#   filter(significant == "Significant") %>%
#   mutate(
#     direction = case_when(
#       estimate < 1 ~ "Lower risk",
#       estimate > 1 ~ "Higher risk",
#       TRUE         ~ NA_character_
#     )
#   )

# Order: direction first, then count of significance
term_levels_facet <- sig_only_glmm %>%
  mutate(direction_order = if_else(estimate > 1, 1, 2)) %>%
  group_by(predictor_set, term_plot, direction_order) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(predictor_set, direction_order, desc(n), term_plot) %>%
  pull(term_plot)

# Plot
ggplot(sig_only_glmm,
       aes(x = outcome, y = factor(term_plot, levels = term_levels_facet), fill = direction)
) +
  geom_tile(colour = "white") +
  facet_grid(predictor_set ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Lower risk" = "#3cb371", "Higher risk" = "#d73027")) +
  labs(
    x = "Cardiovascular risk score",
    y = "Predictor",
    fill = "",
    title = "GLMM: Significant predictors across cardiovascular risk scores",
    caption = "Model: outcome ~ predictor + (1|random effects) (Gamma, log link)\nRandom effects: Gender & Country (REDCap, Risk factors, Body composition)\nRandom effects: Statins & Supplements (Lipids, Fatty acids)\nOutliers excluded (1st/99th percentile) for Lipids & Fatty acids only"
  ) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    panel.background = element_rect(fill = "white"),
    plot.caption = element_text(size = 6, hjust = 1),
    strip.text.y = element_text(size = 8, face = "bold", angle = 0)
)

# 8.4.2 ONly plotting what is significant at least twice
# Filter for predictors significant at least twice
sig_twice_glmm <- sig_only_glmm %>%
  group_by(term_plot) %>%
  filter(n() >= 2) %>%
  ungroup()

# Order
term_levels_twice <- sig_twice_glmm %>%
  mutate(direction_order = if_else(estimate > 1, 1, 2)) %>%
  group_by(predictor_set, term_plot, direction_order) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(predictor_set, direction_order, desc(n), term_plot) %>%
  pull(term_plot)

# Plot
ggplot(sig_twice_glmm,
       aes(x = outcome, y = factor(term_plot, levels = term_levels_twice), fill = direction)
) +
  geom_tile(colour = "white") +
  facet_grid(predictor_set ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Lower risk" = "#3cb371", "Higher risk" = "#d73027")) +
  labs(
    x = "Cardiovascular risk score",
    y = "Predictor",
    fill = "",
    title = "GLMM: Predictors significant in ≥2 risk scores",
    caption = "Model: outcome ~ predictor + (1|random effects) (Gamma, log link)\nRandom effects: Gender & Country (REDCap, Risk factors, Body composition)\nRandom effects: Statins & Supplements (Lipids, Fatty acids)\nOutliers excluded (1st/99th percentile) for Lipids & Fatty acids only"
  ) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    panel.background = element_rect(fill = "white"),
    plot.caption = element_text(size = 6, hjust = 1),
    strip.text.y = element_text(size = 8, face = "bold", angle = 0)
  )

# 8.4.3 plot at least 2 times significant and with colour coding of association strength
# Step 1: Add log effect size to your data
sig_twice_glmm <- sig_twice_glmm %>%
  mutate(
    log_effect = log(estimate),  # Log of exp(β) - negative = protective, positive = adverse
    direction = if_else(estimate < 1, "Lower risk", "Higher risk")
  )

# Step 3: Plot with adjusted color scale
ggplot(sig_twice_glmm,
       aes(x = outcome, 
           y = factor(term_plot, levels = term_levels_twice), 
           fill = log_effect)) +
  geom_tile(colour = "white") +
  facet_grid(predictor_set ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#0072B2",           # Blue for protective
    mid = "white",             # White for neutral
    high = "#D55E00",          # Orange for adverse
    midpoint = 0,
    limits = c(-1.0, 1.1),     # Symmetric, slightly beyond your range
    name = "log(exp(β))",
    breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
    labels = c("-0.8", "-0.4", "0", "+0.4", "+0.8")
  ) +
  labs(
    x = "Cardiovascular risk score",
    y = "Predictor",
    title = "GLMM: Predictors significant in ≥2 risk scores",
    caption = "Model: outcome ~ predictor + (1|random effects) (Gamma, log link)\nRandom effects: Gender & Country (REDCap, Risk factors, Body composition)\nRandom effects: Statins & Supplements (Lipids, Fatty acids)\nOutliers excluded (1st/99th percentile) for Lipids & Fatty acids only"
  ) +
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 6),
    panel.background = element_rect(fill = "white"),
    plot.caption = element_text(size = 6, hjust = 1),
    strip.text.y = element_text(size = 7, face = "bold", angle = 0),
    legend.position = "right"
  )


# 8.5 GLMM forest plots with QRISK3 only
body_comp_glmm_plot <- make_forest_one_outcome(
  df           = glmm_body_comp_num,
  outcome_name = "QRISK3_risk",
  title_prefix = "GLMM: Body composition metrics"
) +
  labs(caption = "Model: QRISK3 ~ predictor + (1|Gender) + (1|Country) (Gamma, log link)") +
  theme(plot.caption = element_text(size = 6, hjust = 1))
body_comp_glmm_plot

# Lipids
lipid_glmm_plot <- make_forest_one_outcome(
  df           = glmm_lipid_num,
  outcome_name = "QRISK3_risk",
  title_prefix = "GLMM: Lipids"
) +
  labs(caption = "Model: QRISK3 ~ predictor + (1|Statins) + (1|Supplements) (Gamma, log link)\nOutliers excluded (1st/99th percentile)") +
  theme(plot.caption = element_text(size = 6, hjust = 1))
lipid_glmm_plot

# Fatty acids
fatty_glmm_plot <- make_forest_one_outcome(
  df           = glmm_fatty_num,
  outcome_name = "QRISK3_risk",
  title_prefix = "GLMM: Fatty acids"
) +
  labs(caption = "Model: QRISK3 ~ predictor + (1|Statins) + (1|Supplements) (Gamma, log link)\nOutliers excluded (1st/99th percentile)") +
  theme(plot.caption = element_text(size = 6, hjust = 1))
fatty_glmm_plot

# REDcap
REDcap_glmm_all <- bind_rows(glmm_REDcap_num, glmm_REDcap_fac)
REDcap_glmm_plot <- make_forest_one_outcome(
  df           = REDcap_glmm_all,
  outcome_name = "QRISK3_risk",
  title_prefix = "GLMM: REDCap demographics"
) +
  labs(caption = "Model: QRISK3 ~ predictor + (1|Gender) + (1|Country) (Gamma, log link)") +
  theme(plot.caption = element_text(size = 6, hjust = 1))
REDcap_glmm_plot

# Risk factors
risk_glmm_all <- bind_rows(glmm_risk_num, 
                           glmm_risk_fac %>% filter(term != "Age.Risk = III"))
risk_glmm_plot <- make_forest_one_outcome(
  df           = risk_glmm_all,
  outcome_name = "QRISK3_risk",
  title_prefix = "GLMM: Clinical risk factors"
) +
  labs(caption = "Model: QRISK3 ~ predictor + (1|Gender) + (1|Country) (Gamma, log link)") +
  theme(plot.caption = element_text(size = 6, hjust = 1))
risk_glmm_plot