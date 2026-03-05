################################################################################
### Cardiovascular Risk Scores and Predictor Associations: 
### Generalised Linear Models with Fixed Effects (statin, supplements, study site + Age)
### Author: Luisa Delius
################################################################################

# ============================================
# Set everything up
# ============================================
library(tidyverse)
library(broom)
library(patchwork)
library(flextable)
library(DHARMa)

set.seed(42)

wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files/processed_data"
setwd(wkdir)

# ============================================
# 1. Load data
# ============================================
df_fatty_acids_predictors <- readRDS("df_fatty_acids_predictor_statin_suppl.rds") %>% arrange(Sample_ID)
df_lipidomics_predictors  <- readRDS("df_lipidomics_predictor_statin_suppl.rds") %>% arrange(Sample_ID)
df_REDcap_demographics    <- readRDS("df_REDcap_demographics_predictor.rds") %>% arrange(Sample_ID)
df_urine_nmr_data         <- readRDS("df_urine_NMR_data.rds") %>% arrange(Sample_ID)

df_risk_factors <- readRDS("df_risk_factor_predictors.rds") %>%
  select(-Age.Risk, -stress_resilience_status, -stress_index_status, -Heart.Rate,
         -z_Heart.Rate, -z_Age, -z_Total.Cholesterol.mg.dl, -z_HDL.mg.dl,
         -z_Body.Weight, -z_Height, -z_BMI) %>%
  arrange(Sample_ID)

df_all_cvd_risk_scores <- readRDS("df_all_risk_scores.rds") %>%
  rename(Sample_ID = PatientID) %>%
  select(-SCORE2_strat) %>%
  arrange(Sample_ID)

df_body_composition <- readRDS("df_body_composition_metrics.rds") %>%
  arrange(Sample_ID) %>%
  full_join(df_risk_factors %>% select(z_Body.fat, z_Fat.free.mass, Sample_ID), by = "Sample_ID")

# Load the original risk factors to get Age before it was dropped
df_age_bmi <- readRDS("df_risk_factor_predictors.rds") %>%
  select(Sample_ID, Age) %>%
  mutate(z_Age = as.numeric(scale(Age))) %>%
  arrange(Sample_ID)

# Prepare Merged Data frames for GLM 
df_fixed_effects <- df_risk_factors %>%
  select(Sample_ID, Gender, Country) %>%
  full_join(df_lipidomics_predictors %>% select(Sample_ID, Statins, Supplements), by = "Sample_ID") %>%
  full_join(df_age_bmi %>% select(Sample_ID, z_Age), by = "Sample_ID")

df_cvd_and_fatty_acids <- df_all_cvd_risk_scores %>%
  full_join(df_fatty_acids_predictors %>% select(-QRISK3_risk), by = "Sample_ID") %>%
  full_join(df_fixed_effects %>% select(Sample_ID, Gender, Country, z_Age), by = "Sample_ID")

df_cvd_and_lipidomics <- df_all_cvd_risk_scores %>%
  full_join(df_lipidomics_predictors %>% select(-QRISK3_risk), by = "Sample_ID") %>%
  full_join(df_fixed_effects %>% select(Sample_ID, Gender, Country, z_Age), by = "Sample_ID")

df_cvd_and_REDcap <- df_all_cvd_risk_scores %>%
  full_join(df_REDcap_demographics %>% select(-QRISK3_risk, -recruitment_site), by = "Sample_ID") %>%
  full_join(df_fixed_effects, by = "Sample_ID")

df_cvd_and_risk_factors <- df_all_cvd_risk_scores %>%
  full_join(df_risk_factors %>% select(-QRISK3_risk), by = "Sample_ID") %>%
  full_join(df_fixed_effects %>% select(Sample_ID, Statins, Supplements, z_Age), by = "Sample_ID")

df_cvd_and_body_comp <- df_all_cvd_risk_scores %>%
  full_join(df_body_composition %>% select(-QRISK3_risk), by = "Sample_ID") %>%
  full_join(df_fixed_effects, by = "Sample_ID")

df_cvd_and_urine_nmr <- df_all_cvd_risk_scores %>%
  full_join(df_urine_nmr_data, by = "Sample_ID") %>%
  full_join(df_fixed_effects, by = "Sample_ID")

# ============================================
# 3. Define Predictors and Outcomes
# ============================================
lipid_predictors        <- names(select(df_lipidomics_predictors, starts_with("z_")))
fatty_acid_predictors   <- names(select(df_fatty_acids_predictors, starts_with("z_")))
REDcap_numeric_preds    <- names(select(df_REDcap_demographics, starts_with("z_")))
REDcap_factor_preds     <- names(select(df_REDcap_demographics, where(is.factor), -recruitment_site))
risk_factor_predictors  <- names(select(df_risk_factors, starts_with("z_"), -z_Body.fat, -z_Fat.free.mass))
body_comp_predictors    <- names(select(df_body_composition, starts_with("z_")))
urine_nmr_predictors    <- names(select(df_urine_nmr_data, starts_with("z_")))

# Exclude z_Age (covariate, not a predictor)
lipid_predictors       <- setdiff(lipid_predictors, c("z_Age"))
fatty_acid_predictors  <- setdiff(fatty_acid_predictors, c("z_Age"))
REDcap_numeric_preds   <- setdiff(REDcap_numeric_preds, c("z_Age"))
risk_factor_predictors <- setdiff(risk_factor_predictors, c("z_Age"))
body_comp_predictors   <- setdiff(body_comp_predictors, c("z_Age"))
urine_nmr_predictors   <- setdiff(urine_nmr_predictors, c("z_Age"))

fixed_effects <- c("Statins", "Supplements", "Gender", "Country", "z_Age")

outcomes <- c("QRISK3_risk", "SCORE2_score", "ascvd_10y", "frs_10y")

# ============================================
# 4. GLM Functions
# ============================================
# 4.1 GLM for numeric predictors
run_glm_num <- function(data, outcomes, num_predictors, fixed_effects) {
  fixed_part <- paste(fixed_effects, collapse = " + ")
  models_list <- list()
  
  results <- map_dfr(num_predictors, function(var) {
    map_dfr(outcomes, function(outcome) {
      df_pair <- data %>%
        select(all_of(c(outcome, var, fixed_effects))) %>%
        drop_na()
      
      if (nrow(df_pair) == 0) return(tibble())
      for (fe in fixed_effects) {
        if (is.factor(df_pair[[fe]]) && n_distinct(df_pair[[fe]]) < 2) return(tibble())
      }
      
      model <- glm(
        as.formula(paste(outcome, "~", var, "+", fixed_part)),
        data = df_pair, family = Gamma(link = "log")
      )
      
      models_list[[paste(var, outcome, sep = "___")]] <<- model
      
      tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(term == var) %>%
        mutate(predictor = var, level = NA_character_, outcome = outcome, n = nrow(df_pair))
    })
  })
  
  list(results = results, models = models_list)
}

# 4.2 GLM for factor predictors
run_glm_fac <- function(data, outcomes, factor_predictors, fixed_effects) {
  fixed_part <- paste(fixed_effects, collapse = " + ")
  
  map_dfr(factor_predictors, function(var) {
    map_dfr(outcomes, function(outcome) {
      df_pair <- data %>%
        select(all_of(c(outcome, var, fixed_effects))) %>%
        drop_na()
      
      if (length(levels(df_pair[[var]])) < 2 || nrow(df_pair) == 0) return(tibble())
      
      model <- glm(
        as.formula(paste(outcome, "~", var, "+", fixed_part)),
        data = df_pair, family = Gamma(link = "log")
      )
      
      ref_level <- levels(df_pair[[var]])[1]
      
      tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(str_starts(term, var)) %>%
        mutate(
          predictor = var,
          level = str_remove(term, paste0("^", var)),
          comparison = paste0(level, " vs ", ref_level),
          outcome = outcome, n = nrow(df_pair)
        )
    })
  })
}

# ============================================
# 5. Run GLMs
# ============================================
glm_lipid     <- run_glm_num(df_cvd_and_lipidomics, outcomes, lipid_predictors, fixed_effects)
glm_fatty     <- run_glm_num(df_cvd_and_fatty_acids, outcomes, fatty_acid_predictors, fixed_effects)
glm_urine     <- run_glm_num(df_cvd_and_urine_nmr, outcomes, urine_nmr_predictors, fixed_effects)
glm_REDcap    <- run_glm_num(df_cvd_and_REDcap, outcomes, REDcap_numeric_preds, fixed_effects)
glm_risk      <- run_glm_num(df_cvd_and_risk_factors, outcomes, risk_factor_predictors, fixed_effects)
glm_body_comp <- run_glm_num(df_cvd_and_body_comp, outcomes, body_comp_predictors, fixed_effects)
glm_REDcap_fac <- run_glm_fac(df_cvd_and_REDcap, outcomes, REDcap_factor_preds, fixed_effects)

glm_lipid_num     <- glm_lipid$results
glm_fatty_num     <- glm_fatty$results
glm_urine_num     <- glm_urine$results
glm_REDcap_num    <- glm_REDcap$results
glm_risk_num      <- glm_risk$results
glm_body_comp_num <- glm_body_comp$results

all_models <- c(glm_lipid$models, glm_fatty$models, glm_urine$models,
                glm_REDcap$models, glm_risk$models, glm_body_comp$models)

incremental_dir <- file.path(wkdir, "glm_incremental_results_age_adjusted")
dir.create(incremental_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================
# 6. Combine Results and Apply Multiple Testing Correction
# ============================================
all_results_glm <- bind_rows(
  glm_lipid_num     %>% mutate(predictor_set = "Lipids", comparison = NA_character_),
  glm_fatty_num     %>% mutate(predictor_set = "Fatty acids", comparison = NA_character_),
  glm_urine_num     %>% mutate(predictor_set = "Urine NMR", comparison = NA_character_),
  glm_REDcap_num    %>% mutate(predictor_set = "REDCap numeric", comparison = NA_character_),
  glm_REDcap_fac    %>% mutate(predictor_set = "REDCap factors"),
  glm_risk_num      %>% mutate(predictor_set = "Risk factors", comparison = NA_character_),
  glm_body_comp_num %>% mutate(predictor_set = "Body composition", comparison = NA_character_)
) %>%
  mutate(
    term_plot = if_else(is.na(level), predictor, term),
    p.adjusted = p.adjust(p.value, method = "BH"),
    significant = if_else(p.adjusted < 0.05, "Significant", "Not significant")
  )

# Save Combined Results 
saveRDS(all_results_glm, file.path(incremental_dir, "all_results_glm_combined_age_adj.rds"))

# ============================================
# 7. Plotting Helpers
# ============================================
rename_predictors <- function(term_plot) {
  case_when(
    term_plot == "z_whtr_waist_height_ratio_"      ~ "Waist-to-Height Ratio",
    term_plot == "z_LDL.mg.dl"                     ~ "LDL cholesterol",
    term_plot == "z_Hba1C"                         ~ "HbA1c",
    term_plot == "z_TAG.mg.dl"                    ~ "Triacylglycerol",
    term_plot == "z_d_stearic_18_0"                ~ "Stearic acid (18:0)",
    term_plot == "z_fahfa_18_1_18_0"              ~ "FAHFA 18:1/18:0",
    term_plot == "z_tartaric_acid"                ~ "Tartaric Acid",
    TRUE ~ term_plot
  )
}

rename_outcomes <- function(outcome) {
  case_when(
    outcome == "ascvd_10y"     ~ "ASCVD",
    outcome == "frs_10y"       ~ "Framingham",
    outcome == "QRISK3_risk"   ~ "QRISK3",
    outcome == "SCORE2_score"  ~ "SCORE2",
    TRUE ~ outcome
  )
}

rename_predictor_sets <- function(predictor_set) {
  case_when(
    predictor_set == "Body composition"  ~ "Body Composition",
    predictor_set == "Fatty acids"       ~ "Fatty Acids",
    predictor_set == "REDCap factors"    ~ "Lifestyle &\nSociodemographic Factors",
    predictor_set == "REDCap numeric"    ~ "Lifestyle &\nSociodemographic Factors",
    predictor_set == "Risk factors"      ~ "Clinical Measurements &\nSupplementary Biomarkers",
    predictor_set == "Lipids"            ~ "Lipid Species",
    predictor_set == "Urine NMR"        ~ "Urinary Metabolites",
    TRUE ~ predictor_set
  )
}

# ============================================
# 8 Combined Figure: Age-Adjusted GLM Heatmap
# ============================================
# Data Preparation
heatmap_data <- all_results_glm %>%
  filter(!is.na(estimate)) %>%
  mutate(
    log_effect      = log(estimate),
    is_fdr_sig      = p.adjusted < 0.05,
    is_nominal_sig  = p.value < 0.05 & p.adjusted >= 0.05,
    outcome         = rename_outcomes(outcome),
    predictor_set   = rename_predictor_sets(predictor_set),
    term_plot_clean = rename_predictors(term_plot),
    is_categorical  = grepl("Lifestyle", predictor_set)
  )

fdr_sig_terms <- heatmap_data %>%
  filter(is_fdr_sig) %>%
  pull(term_plot_clean) %>%
  unique()

fdr_sig_predictors <- all_results_glm %>%
  filter(p.adjusted < 0.05) %>%
  pull(predictor) %>%
  unique()

heatmap_exploratory <- heatmap_data %>%
  filter(!predictor %in% fdr_sig_predictors)

preds_3plus <- heatmap_exploratory %>%
  filter(is_nominal_sig) %>%
  group_by(term_plot_clean) %>%
  filter(n() >= 3) %>%
  pull(term_plot_clean) %>%
  unique()

# Panel data (continuous only for age-adjusted)
panel_a_cont <- heatmap_data %>% filter(term_plot_clean %in% fdr_sig_terms, !is_categorical)
panel_b_cont <- heatmap_exploratory %>% filter(term_plot_clean %in% preds_3plus, !is_categorical)

# Ordering
order_a_cont <- panel_a_cont %>%
  filter(is_fdr_sig) %>%
  group_by(term_plot_clean, predictor_set) %>%
  summarise(n_sig = n(), mean_effect = mean(abs(log_effect)), .groups = "drop") %>%
  arrange(predictor_set, desc(n_sig), desc(mean_effect)) %>%
  pull(term_plot_clean)

order_b_cont <- panel_b_cont %>%
  group_by(term_plot_clean, predictor_set) %>%
  summarise(mean_effect = mean(abs(log_effect)), .groups = "drop") %>%
  arrange(predictor_set, desc(mean_effect)) %>%
  pull(term_plot_clean)

# Use the same max_abs as the main figure for direct comparability
max_abs_cont <- 0.78

# Shared Theme
legend_theme <- theme(
  legend.position    = "right",
  legend.title       = element_text(size = 10),
  legend.text        = element_text(size = 9),
  legend.key.height  = unit(0.6, "cm"),
  legend.key.width   = unit(0.4, "cm"),
  legend.margin      = margin(0, 0, 0, 0),
  legend.box.margin  = margin(0, 0, 0, 0)
)

theme_heatmap <- theme_minimal(base_size = 14) +
  theme(
    axis.text.y      = element_text(size = 13),
    axis.title.y     = element_text(size = 11, face = "bold"),
    strip.text.y     = element_text(size = 12, face = "bold", angle = 0, lineheight = 0.85),
    panel.grid       = element_blank(),
    panel.spacing    = unit(0.5, "lines")
  ) +
  legend_theme

fill_continuous <- scale_fill_gradient2(
  low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
  name = "log(Rate Ratio)",
  limits = c(-max_abs_cont, max_abs_cont),
  oob = scales::squish,
  breaks = seq(-max_abs_cont, max_abs_cont, length.out = 5),
  labels = function(x) sprintf("%.2f", x)
)

# Title labels
title_a <- wrap_elements(
  full = grid::textGrob(
    label = expression(bold("A) Significant Associations After Multiple Testing Correction (") * bold(italic(p)[adj]) * bold(" < 0.05, *)")),
    x = 0, hjust = 0,
    gp = grid::gpar(fontsize = 14)
  )
)

title_b <- wrap_elements(
  full = grid::textGrob(
    label = expression(bold("B) Exploratory Associations (nominal ") * bold(italic(p)) * bold(" < 0.05, ● ; ≥3 CVD scores)")),
    x = 0, hjust = 0,
    gp = grid::gpar(fontsize = 14)
  )
)

# Panel A: Continuous (FDR-significant)
p_a_cont <- ggplot(panel_a_cont,
                   aes(x = outcome,
                       y = factor(term_plot_clean, levels = rev(order_a_cont)),
                       fill = log_effect)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(data = panel_a_cont %>% filter(is_fdr_sig),
            aes(label = "*"), size = 8, vjust = 0.75, colour = "black") +
  facet_grid(predictor_set ~ ., scales = "free_y", space = "free_y") +
  fill_continuous +
  labs(x = NULL, y = "Continuous Predictors") +
  theme_heatmap +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(2, 5, 2, 5)
  )

# Panel B: Continuous (exploratory)
p_b_cont <- ggplot(panel_b_cont,
                   aes(x = outcome,
                       y = factor(term_plot_clean, levels = rev(order_b_cont)),
                       fill = log_effect)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_point(data = panel_b_cont %>% filter(is_nominal_sig),
             shape = 21, size = 3, fill = "black", colour = "white", stroke = 0.5) +
  facet_grid(predictor_set ~ ., scales = "free_y", space = "free_y") +
  fill_continuous +
  labs(x = NULL, y = "Continuous Predictors") +
  theme_heatmap +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 13, face = "bold"),
    plot.margin  = margin(2, 5, 5, 5)
  )

# Stack everything 
n_a_cont <- length(unique(panel_a_cont$term_plot_clean))
n_b_cont <- length(unique(panel_b_cont$term_plot_clean))

combined_final <- title_a / p_a_cont / title_b / p_b_cont +
  plot_layout(
    heights = c(1, n_a_cont, 1, n_b_cont)
  )

# Save
total_rows <- n_a_cont + n_b_cont
plot_height <- total_rows * 0.45 + 5

ggsave("combined_glm_heatmap_age_adj.png", combined_final,
       width = 13, height = plot_height, dpi = 300)


# Save individual plots
individual_fdr <- title_a / p_a_cont +
  plot_layout(heights = c(1, n_a_cont)) +
  plot_annotation(
    caption = paste0(
      "Model: outcome ~ predictor + Statins + Supplements + Sex + Country + Age (Gamma GLM, log link)\n",
      "Multiple testing correction: Benjamini–Hochberg"
    ),
    theme = theme(plot.caption = element_text(size = 11, hjust = 0))
  )

ggsave("heatmap_significant_predictors_split_age_adj.png", individual_fdr,
       width = 12, height = n_a_cont * 0.4 + 3, dpi = 300)

individual_exp <- title_b / p_b_cont +
  plot_layout(heights = c(1, n_b_cont)) +
  plot_annotation(
    caption = paste0(
      "Model: outcome ~ predictor + Statins + Supplements + Sex + Country + Age (Gamma GLM, log link)\n",
      "Excludes FDR-significant predictors"
    ),
    theme = theme(plot.caption = element_text(size = 11, hjust = 0))
  )

ggsave("exploratory_heatmap_3plus_scores_split_age_adj.png", individual_exp,
       width = 12, height = n_b_cont * 0.4 + 3, dpi = 300)


# ============================================
# 9. Supplementary Tables
# ============================================
create_wide_supp_table <- function(data, predictor_set_name) {
  df_wide <- data %>%
    filter(predictor_set == predictor_set_name) %>%
    mutate(
      term_clean = str_remove(term_plot, "^z_"),
      RR_CI  = sprintf("%.2f(%.2f-%.2f)", estimate, conf.low, conf.high),
      p_val  = ifelse(p.value < 0.001, sprintf("%.2e", p.value), sprintf("%.3f", p.value)),
      p_adj  = ifelse(p.adjusted < 0.001, sprintf("%.2e", p.adjusted), sprintf("%.3f", p.adjusted))
    ) %>%
    select(term_clean, outcome, RR_CI, p_val, p_adj) %>%
    pivot_wider(names_from = outcome, values_from = c(RR_CI, p_val, p_adj), names_sep = "___")
  
  df_wide <- df_wide %>%
    select(
      term_clean,
      starts_with("RR_CI___QRISK3"), starts_with("p_val___QRISK3"), starts_with("p_adj___QRISK3"),
      starts_with("RR_CI___SCORE2"), starts_with("p_val___SCORE2"), starts_with("p_adj___SCORE2"),
      starts_with("RR_CI___ascvd"),  starts_with("p_val___ascvd"),  starts_with("p_adj___ascvd"),
      starts_with("RR_CI___frs"),    starts_with("p_val___frs"),    starts_with("p_adj___frs")
    )
  
  ft <- df_wide %>%
    flextable() %>%
    set_header_labels(
      term_clean = "Predictor",
      RR_CI___QRISK3_risk = "RR(95%CI)", p_val___QRISK3_risk = "p", p_adj___QRISK3_risk = "p.adj",
      RR_CI___SCORE2_score = "RR(95%CI)", p_val___SCORE2_score = "p", p_adj___SCORE2_score = "p.adj",
      RR_CI___ascvd_10y = "RR(95%CI)", p_val___ascvd_10y = "p", p_adj___ascvd_10y = "p.adj",
      RR_CI___frs_10y = "RR(95%CI)", p_val___frs_10y = "p", p_adj___frs_10y = "p.adj"
    ) %>%
    add_header_row(
      values = c("", "QRISK3", "SCORE2", "ASCVD", "Framingham"),
      colwidths = c(1, 3, 3, 3, 3)
    ) %>%
    theme_booktabs() %>%
    width(j = 1, width = 1.5) %>%
    width(j = c(2, 5, 8, 11), width = 1.0) %>%
    width(j = c(3, 6, 9, 12), width = 0.5) %>%
    width(j = c(4, 7, 10, 13), width = 0.5) %>%
    align(j = 1, align = "left", part = "all") %>%
    align(j = 2:ncol(df_wide), align = "center", part = "all") %>%
    bold(part = "header") %>%
    merge_h(part = "header") %>%
    vline(j = c(1, 4, 7, 10), border = fp_border(width = 1.5), part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    hrule(rule = "exact")
  
  return(ft)
}

ft_lipids       <- create_wide_supp_table(all_results_glm, "Lipids")
ft_fatty_acids  <- create_wide_supp_table(all_results_glm, "Fatty acids")
ft_urine_nmr    <- create_wide_supp_table(all_results_glm, "Urine NMR")
ft_redcap_num   <- create_wide_supp_table(all_results_glm, "REDCap numeric")
ft_redcap_fac   <- create_wide_supp_table(all_results_glm, "REDCap factors")
ft_risk_factors <- create_wide_supp_table(all_results_glm, "Risk factors")
ft_body_comp    <- create_wide_supp_table(all_results_glm, "Body composition")

save_as_docx(
  "Table S2: Lipid Species (Age-Adj)"    = ft_lipids,
  "Table S3: Fatty Acids (Age-Adj)"      = ft_fatty_acids,
  "Table S4: Urine NMR (Age-Adj)"        = ft_urine_nmr,
  "Table S5: REDCap Numeric (Age-Adj)"   = ft_redcap_num,
  "Table S6: REDCap Factors (Age-Adj)"   = ft_redcap_fac,
  "Table S7: Risk Factors (Age-Adj)"     = ft_risk_factors,
  "Table S8: Body Composition (Age-Adj)" = ft_body_comp,
  path = "supplementary_tables_GLM_all_age_adj.docx",
  pr_section = prop_section(page_size = page_size(orient = "landscape"))
)

##################################################################################
# 10. GLM Evaluation by deviance stats and DHARMa
##################################################################################

# 10.1. Deviance Statistics
deviance_all <- map_dfr(names(all_models), function(key) {
  model <- all_models[[key]]
  parts <- str_split(key, "___", simplify = TRUE)
  
  tibble(
    predictor      = parts[1],
    outcome        = parts[2],
    n              = nrow(model$data),
    deviance       = model$deviance,
    df_residual    = model$df.residual,
    deviance_ratio = model$deviance / model$df.residual,
    AIC            = AIC(model)
  )
})

deviance_summary <- deviance_all %>%
  group_by(outcome) %>%
  summarise(
    n_models         = n(),
    mean_dev_ratio   = mean(deviance_ratio),
    median_dev_ratio = median(deviance_ratio),
    sd_dev_ratio     = sd(deviance_ratio)
  ) %>%
  arrange(desc(mean_dev_ratio))

print("Deviance Summary by Outcome (Age-Adjusted):")
print(deviance_summary)

# 10.2. DHARMa Diagnostics
# Run on one representative predictor per outcome
pred_name <- risk_factor_predictors[1]

dharma_results <- map_dfr(outcomes, function(outcome) {
  key <- paste(pred_name, outcome, sep = "___")
  model <- all_models[[key]]
  
  sim_resid <- simulateResiduals(model, n = 1000, plot = FALSE)
  
  tibble(
    outcome        = outcome,
    predictor      = pred_name,
    n              = nrow(model$data),
    uniformity_p   = testUniformity(sim_resid, plot = FALSE)$p.value,
    dispersion_p   = testDispersion(sim_resid, plot = FALSE)$p.value,
    dispersion_stat = testDispersion(sim_resid, plot = FALSE)$statistic,
    outlier_p      = testOutliers(sim_resid, plot = FALSE)$p.value
  )
})

print("DHARMa Diagnostics by Outcome (Age-Adjusted):")
print(dharma_results)

# 10.3. Combined Summary Table
diagnostics_summary <- deviance_summary %>%
  left_join(
    dharma_results %>% select(outcome, uniformity_p, dispersion_stat, dispersion_p, outlier_p),
    by = "outcome"
  ) %>%
  mutate(
    outcome_label = case_when(
      outcome == "QRISK3_risk"  ~ "QRISK3",
      outcome == "SCORE2_score" ~ "SCORE2",
      outcome == "ascvd_10y"    ~ "ASCVD",
      outcome == "frs_10y"      ~ "Framingham"
    )
  ) %>%
  select(outcome_label, n_models, median_dev_ratio,
         uniformity_p, dispersion_stat, dispersion_p, outlier_p)

print("Combined Model Diagnostics (Age-Adjusted):")
print(diagnostics_summary)

write_csv(deviance_summary, "deviance_summary_age_adj.csv")
write_csv(dharma_results, "dharma_diagnostics_age_adj.csv")
write_csv(diagnostics_summary, "combined_diagnostics_summary_age_adj.csv")

saveRDS(deviance_all,        file.path(incremental_dir, "deviance_all_age_adj.rds"))
saveRDS(dharma_results,      file.path(incremental_dir, "dharma_results_age_adj.rds"))
saveRDS(diagnostics_summary, file.path(incremental_dir, "diagnostics_summary_age_adj.rds"))