### Cardiovascular Risk Scores and Predictor Associations: 
### Generalised Linear Models with Fixed Effects (Age-Adjusted)
### Author: Luisa Delius

# ─── Packages ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(broom)
library(patchwork)
library(flextable)
library(DHARMa)

set.seed(42)

# ─── 1. Load Data ─────────────────────────────────────────────────────────────
wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files/processed_data"
setwd(wkdir)

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

# ─── 1b. Load Age and z-scale it ─────────────────────────────────────────────
# Load the original risk factors to get Age before it was dropped
df_age_bmi <- readRDS("df_risk_factor_predictors.rds") %>%
  select(Sample_ID, Age) %>%
  mutate(z_Age = as.numeric(scale(Age))) %>%
  arrange(Sample_ID)

# df_bmi <- readRDS("df_risk_factor_predictors.rds") %>%
#   select(Sample_ID, z_BMI) %>%
#   arrange(Sample_ID)

# ─── 2. Prepare Merged DataFrames for GLM ────────────────────────────────────
df_fixed_effects <- df_risk_factors %>%
  select(Sample_ID, Gender, Country) %>%
  full_join(df_lipidomics_predictors %>% select(Sample_ID, Statins, Supplements), by = "Sample_ID") %>%
  full_join(df_age_bmi %>% select(Sample_ID, z_Age), by = "Sample_ID")
  # full_join(df_bmi %>% select(Sample_ID, z_BMI), by = "Sample_ID")

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

# ─── 3. Define Predictors and Outcomes ────────────────────────────────────────
lipid_predictors        <- names(select(df_lipidomics_predictors, starts_with("z_")))
fatty_acid_predictors   <- names(select(df_fatty_acids_predictors, starts_with("z_")))
REDcap_numeric_preds    <- names(select(df_REDcap_demographics, starts_with("z_")))
REDcap_factor_preds     <- names(select(df_REDcap_demographics, where(is.factor), -recruitment_site))
risk_factor_predictors  <- names(select(df_risk_factors, starts_with("z_"), -z_Body.fat, -z_Fat.free.mass))
body_comp_predictors    <- names(select(df_body_composition, starts_with("z_")))
urine_nmr_predictors    <- names(select(df_urine_nmr_data, starts_with("z_")))

# Safety: remove z_Age from any predictor set in case it slipped in
# Add this line alongside the existing setdiff calls:
lipid_predictors       <- setdiff(lipid_predictors, c("z_Age"))
fatty_acid_predictors  <- setdiff(fatty_acid_predictors, c("z_Age"))
REDcap_numeric_preds   <- setdiff(REDcap_numeric_preds, c("z_Age"))
risk_factor_predictors <- setdiff(risk_factor_predictors, c("z_Age"))
body_comp_predictors   <- setdiff(body_comp_predictors, c("z_Age"))
urine_nmr_predictors   <- setdiff(urine_nmr_predictors, c("z_Age"))

fixed_effects <- c("Statins", "Supplements", "Gender", "Country", "z_Age")

outcomes <- c("QRISK3_risk", "SCORE2_score", "ascvd_10y", "frs_10y")

# ─── 4. GLM Functions ────────────────────────────────────────────────────────
# 4.1 GLM for numeric (continuous) predictors
run_glm_num <- function(data, outcomes, num_predictors, fixed_effects) {
  fixed_part <- paste(fixed_effects, collapse = " + ")
  
  map_dfr(num_predictors, function(var) {
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
      
      tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(term == var) %>%
        mutate(predictor = var, level = NA_character_, outcome = outcome, n = nrow(df_pair))
    })
  })
}

# 4.2 GLM for factor (categorical) predictors
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

# ─── 5. Run GLMs ─────────────────────────────────────────────────────────────
glm_lipid_num     <- run_glm_num(df_cvd_and_lipidomics, outcomes, lipid_predictors, fixed_effects)
glm_fatty_num     <- run_glm_num(df_cvd_and_fatty_acids, outcomes, fatty_acid_predictors, fixed_effects)
glm_urine_num     <- run_glm_num(df_cvd_and_urine_nmr, outcomes, urine_nmr_predictors, fixed_effects)
glm_REDcap_num    <- run_glm_num(df_cvd_and_REDcap, outcomes, REDcap_numeric_preds, fixed_effects)
glm_risk_num      <- run_glm_num(df_cvd_and_risk_factors, outcomes, risk_factor_predictors, fixed_effects)
glm_body_comp_num <- run_glm_num(df_cvd_and_body_comp, outcomes, body_comp_predictors, fixed_effects)
glm_REDcap_fac    <- run_glm_fac(df_cvd_and_REDcap, outcomes, REDcap_factor_preds, fixed_effects)

# ─── 5b. Incremental Saves ───────────────────────────────────────────────────
incremental_dir <- file.path(wkdir, "glm_incremental_results_age_adjusted")
dir.create(incremental_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(glm_lipid_num,     file.path(incremental_dir, "glm_lipid_num_age_adj.rds"))
saveRDS(glm_fatty_num,     file.path(incremental_dir, "glm_fatty_num_age_adj.rds"))
saveRDS(glm_urine_num,     file.path(incremental_dir, "glm_urine_num_age_adj.rds"))
saveRDS(glm_REDcap_num,    file.path(incremental_dir, "glm_REDcap_num_age_adj.rds"))
saveRDS(glm_risk_num,      file.path(incremental_dir, "glm_risk_num_age_adj.rds"))
saveRDS(glm_body_comp_num, file.path(incremental_dir, "glm_body_comp_num_age_adj.rds"))
saveRDS(glm_REDcap_fac,    file.path(incremental_dir, "glm_REDcap_fac_age_adj.rds"))
cat("Incremental GLM results saved to:", incremental_dir, "\n")

# ─── 6. Combine Results and Apply Multiple Testing Correction ────────────────
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

cat("Total associations tested:", nrow(all_results_glm), "\n")
cat("Significant (BH-adjusted p < 0.05):", sum(all_results_glm$significant == "Significant", na.rm = TRUE), "\n")

# ─── 6b. Save Combined Results ───────────────────────────────────────────────
saveRDS(all_results_glm, file.path(incremental_dir, "all_results_glm_combined_age_adj.rds"))
cat("Combined results saved.\n")

# ─── 7. Plotting Helpers ─────────────────────────────────────────────────────
rename_predictors <- function(term_plot) {
  case_when(
    term_plot == "z_whtr_waist_height_ratio_"      ~ "Waist / Height",
    term_plot == "z_LDL.mg.dl"                     ~ "LDL",
    term_plot == "z_Hba1C"                         ~ "HbA1c",
    term_plot == "z_TAG.mg.dl"                    ~ "TAG",
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

# ─── 8. Figure 1: FDR-Significant Heatmap ────────────────────────────────────
heatmap_data <- all_results_glm %>%
  filter(!is.na(estimate)) %>%
  mutate(
    log_effect     = log(estimate),
    is_significant = p.adjusted < 0.05,
    outcome        = rename_outcomes(outcome),
    predictor_set  = rename_predictor_sets(predictor_set),
    term_plot_clean = rename_predictors(term_plot),
    is_categorical = grepl("Lifestyle", predictor_set)
  )

predictors_with_sig <- heatmap_data %>%
  filter(is_significant) %>%
  pull(term_plot_clean) %>%
  unique()

# ── Guard: skip FDR heatmap if nothing is significant ────────────────────────
if (length(predictors_with_sig) > 0) {
  
  heatmap_data_filtered <- heatmap_data %>% filter(term_plot_clean %in% predictors_with_sig)
  heatmap_continuous    <- heatmap_data_filtered %>% filter(!is_categorical)
  heatmap_categorical   <- heatmap_data_filtered %>% filter(is_categorical)
  
  has_continuous  <- nrow(heatmap_continuous) > 0
  has_categorical <- nrow(heatmap_categorical) > 0
  
  if (has_continuous) {
    term_order_continuous <- heatmap_continuous %>%
      filter(is_significant) %>%
      group_by(term_plot_clean, predictor_set) %>%
      summarise(n_sig = n(), mean_effect = mean(abs(log_effect)), .groups = "drop") %>%
      arrange(predictor_set, desc(n_sig), desc(mean_effect)) %>%
      pull(term_plot_clean)
    
    max_abs_cont <- max(abs(heatmap_continuous$log_effect))
    
    p_continuous <- ggplot(heatmap_continuous,
                           aes(x = outcome,
                               y = factor(term_plot_clean, levels = rev(term_order_continuous)),
                               fill = log_effect)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_text(data = heatmap_continuous %>% filter(is_significant),
                aes(label = "*"), size = 8, vjust = 0.75, colour = "black") +
      facet_grid(predictor_set ~ ., scales = "free_y", space = "free_y") +
      scale_fill_gradient2(
        low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
        name = "log(Rate Ratio)",
        limits = c(-max_abs_cont, max_abs_cont),
        breaks = seq(-max_abs_cont, max_abs_cont, length.out = 5),
        labels = function(x) sprintf("%.2f", x)
      ) +
      labs(x = if (has_categorical) NULL else "Cardiovascular risk score",
           y = "Continuous Predictors") +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = if (has_categorical) element_blank() else element_text(size = 14, angle = 45, hjust = 1),
        axis.ticks.x = if (has_categorical) element_blank() else element_line(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 11, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold", angle = 0, lineheight = 0.85),
        legend.title = element_text(size = 13), legend.text = element_text(size = 12),
        legend.position = "right", panel.grid = element_blank(),
        plot.margin = margin(5, 5, 0, 5), panel.spacing = unit(0.5, "lines")
      )
  }
  
  if (has_categorical) {
    term_order_categorical <- heatmap_categorical %>%
      filter(is_significant) %>%
      group_by(term_plot_clean) %>%
      summarise(n_sig = n(), mean_effect = mean(abs(log_effect)), .groups = "drop") %>%
      arrange(desc(n_sig), desc(mean_effect)) %>%
      pull(term_plot_clean)
    
    max_abs_cat <- max(abs(heatmap_categorical$log_effect))
    
    p_categorical <- ggplot(heatmap_categorical,
                            aes(x = outcome,
                                y = factor(term_plot_clean, levels = rev(term_order_categorical)),
                                fill = log_effect)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_text(data = heatmap_categorical %>% filter(is_significant),
                aes(label = "*"), size = 8, vjust = 0.75, colour = "black") +
      facet_grid(predictor_set ~ ., scales = "free_y", space = "free_y") +
      scale_fill_gradient2(
        low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
        name = "log(Rate Ratio)",
        limits = c(-max_abs_cat, max_abs_cat),
        breaks = seq(-max_abs_cat, max_abs_cat, length.out = 5),
        labels = function(x) sprintf("%.2f", x)
      ) +
      labs(x = "Cardiovascular risk score", y = "Categorical Predictors") +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13, face = "bold"),
        axis.title.y = element_text(size = 11, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold", angle = 0, lineheight = 0.85),
        legend.title = element_text(size = 13), legend.text = element_text(size = 12),
        legend.position = "right", panel.grid = element_blank(),
        plot.margin = margin(0, 5, 5, 5), panel.spacing = unit(0.5, "lines")
      )
  }
  
  # Combine panels
  if (has_continuous && has_categorical) {
    n_cont <- length(unique(heatmap_continuous$term_plot_clean))
    n_cat  <- length(unique(heatmap_categorical$term_plot_clean))
    combined_fdr <- p_continuous / p_categorical +
      plot_layout(heights = c(n_cont, n_cat))
    plot_height <- (n_cont + n_cat) * 0.4 + 3
  } else if (has_continuous) {
    combined_fdr <- p_continuous
    n_cont <- length(unique(heatmap_continuous$term_plot_clean))
    plot_height <- n_cont * 0.4 + 3
  } else {
    combined_fdr <- p_categorical
    n_cat <- length(unique(heatmap_categorical$term_plot_clean))
    plot_height <- n_cat * 0.4 + 3
  }
  
  combined_fdr <- combined_fdr +
    plot_annotation(
      title    = "GLM: Associations between biomarkers and CVD risk scores",
      subtitle = "* indicates FDR-adjusted p < 0.05",
      caption  = paste0(
        "Model: outcome ~ predictor + Statins + Supplements + Sex + Country + Age (Gamma GLM, log link)\n",
        "Multiple testing correction: Benjamini-Hochberg"
      ),
      theme = theme(
        plot.title = element_text(size = 17, face = "bold"),
        plot.subtitle = element_text(size = 14),
        plot.caption = element_text(size = 11, hjust = 0)
      )
    )
  
  ggsave("heatmap_significant_predictors_split_age_adj.png", combined_fdr,
         width = 12, height = plot_height, dpi = 300)
  
} else {
  cat("No FDR-significant associations found. Skipping FDR heatmap.\n")
}

# ─── 9. Figure 2: Exploratory Heatmap (≥3 Nominal Significant) ──────────────
fdr_sig_predictors <- all_results_glm %>%
  filter(p.adjusted < 0.05) %>%
  pull(predictor) %>%
  unique()

heatmap_exploratory <- all_results_glm %>%
  filter(!predictor %in% fdr_sig_predictors, !is.na(estimate)) %>%
  mutate(
    log_effect      = log(estimate),
    is_nominal_sig  = p.value < 0.05,
    outcome         = rename_outcomes(outcome),
    predictor_set   = rename_predictor_sets(predictor_set),
    term_plot_clean = rename_predictors(term_plot),
    is_categorical  = grepl("Lifestyle", predictor_set)
  )

# Keep predictors with ≥3 nominally significant associations
preds_3plus <- heatmap_exploratory %>%
  filter(is_nominal_sig) %>%
  group_by(term_plot_clean) %>%
  filter(n() >= 3) %>%
  pull(term_plot_clean) %>%
  unique()

# ── Guard: skip exploratory heatmap if no predictors qualify ─────────────────
if (length(preds_3plus) > 0) {
  
  heatmap_exp_filtered <- heatmap_exploratory %>% filter(term_plot_clean %in% preds_3plus)
  heatmap_cont_exp <- heatmap_exp_filtered %>% filter(!is_categorical)
  heatmap_cat_exp  <- heatmap_exp_filtered %>% filter(is_categorical)
  
  has_cont_exp <- nrow(heatmap_cont_exp) > 0
  has_cat_exp  <- nrow(heatmap_cat_exp) > 0
  
  if (has_cont_exp) {
    term_order_cont_exp <- heatmap_cont_exp %>%
      group_by(term_plot_clean, predictor_set) %>%
      summarise(mean_effect = mean(abs(log_effect)), .groups = "drop") %>%
      arrange(predictor_set, desc(mean_effect)) %>%
      pull(term_plot_clean)
    
    max_abs_cont_exp <- max(abs(heatmap_cont_exp$log_effect))
    
    p_cont_exp <- ggplot(heatmap_cont_exp,
                         aes(x = outcome,
                             y = factor(term_plot_clean, levels = rev(term_order_cont_exp)),
                             fill = log_effect)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_point(data = heatmap_cont_exp %>% filter(is_nominal_sig),
                 shape = 21, size = 3, fill = "black", colour = "white", stroke = 0.5) +
      facet_grid(predictor_set ~ ., scales = "free_y", space = "free_y") +
      scale_fill_gradient2(
        low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
        name = "log(Rate Ratio)",
        limits = c(-max_abs_cont_exp, max_abs_cont_exp),
        breaks = seq(-max_abs_cont_exp, max_abs_cont_exp, length.out = 5),
        labels = function(x) sprintf("%.2f", x)
      ) +
      labs(x = if (has_cat_exp) NULL else "Cardiovascular risk score",
           y = "Continuous Predictors") +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = if (has_cat_exp) element_blank() else element_text(size = 14, angle = 45, hjust = 1),
        axis.ticks.x = if (has_cat_exp) element_blank() else element_line(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 11, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold", angle = 0, lineheight = 0.85),
        legend.title = element_text(size = 13), legend.text = element_text(size = 12),
        legend.position = "right", panel.grid = element_blank(),
        plot.margin = margin(5, 5, 0, 5), panel.spacing = unit(0.5, "lines")
      )
  }
  
  if (has_cat_exp) {
    term_order_cat_exp <- heatmap_cat_exp %>%
      group_by(term_plot_clean) %>%
      summarise(mean_effect = mean(abs(log_effect)), .groups = "drop") %>%
      arrange(desc(mean_effect)) %>%
      pull(term_plot_clean)
    
    max_abs_cat_exp <- max(abs(heatmap_cat_exp$log_effect))
    
    p_cat_exp <- ggplot(heatmap_cat_exp,
                        aes(x = outcome,
                            y = factor(term_plot_clean, levels = rev(term_order_cat_exp)),
                            fill = log_effect)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_point(data = heatmap_cat_exp %>% filter(is_nominal_sig),
                 shape = 21, size = 3, fill = "black", colour = "white", stroke = 0.5) +
      facet_grid(predictor_set ~ ., scales = "free_y", space = "free_y") +
      scale_fill_gradient2(
        low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
        name = "log(Rate Ratio)",
        limits = c(-max_abs_cat_exp, max_abs_cat_exp),
        breaks = seq(-max_abs_cat_exp, max_abs_cat_exp, length.out = 5),
        labels = function(x) sprintf("%.2f", x)
      ) +
      labs(x = "Cardiovascular risk score", y = "Categorical Predictors") +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13, face = "bold"),
        axis.title.y = element_text(size = 11, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold", angle = 0, lineheight = 0.85),
        legend.title = element_text(size = 13), legend.text = element_text(size = 12),
        legend.position = "right", panel.grid = element_blank(),
        plot.margin = margin(0, 5, 5, 5), panel.spacing = unit(0.5, "lines")
      )
  }
  
  # Combine panels
  if (has_cont_exp && has_cat_exp) {
    n_cont_exp <- length(unique(heatmap_cont_exp$term_plot_clean))
    n_cat_exp  <- length(unique(heatmap_cat_exp$term_plot_clean))
    combined_exploratory <- p_cont_exp / p_cat_exp +
      plot_layout(heights = c(n_cont_exp, n_cat_exp))
    plot_height_exp <- (n_cont_exp + n_cat_exp) * 0.4 + 3
  } else if (has_cont_exp) {
    combined_exploratory <- p_cont_exp
    n_cont_exp <- length(unique(heatmap_cont_exp$term_plot_clean))
    plot_height_exp <- n_cont_exp * 0.4 + 3
  } else {
    combined_exploratory <- p_cat_exp
    n_cat_exp <- length(unique(heatmap_cat_exp$term_plot_clean))
    plot_height_exp <- n_cat_exp * 0.4 + 3
  }
  
  combined_exploratory <- combined_exploratory +
    plot_annotation(
      title    = "GLM: Exploratory Associations",
      subtitle = "● indicates predictors with ≥3 nominal significant associations (p < 0.05, FDR p > 0.05)",
      caption  = paste0(
        "Model: outcome ~ predictor + Statins + Supplements + Sex + Country + Age (Gamma GLM, log link)\n",
        "Excludes FDR-significant predictors"
      ),
      theme = theme(
        plot.title = element_text(size = 17, face = "bold"),
        plot.subtitle = element_text(size = 14),
        plot.caption = element_text(size = 11, hjust = 0)
      )
    )
  
  ggsave("exploratory_heatmap_3plus_scores_split_age_adj.png", combined_exploratory,
         width = 12, height = plot_height_exp, dpi = 300)
  
} else {
  cat("No predictors with ≥3 nominal significant associations. Skipping exploratory heatmap.\n")
}

# ─── 10. Supplementary Tables ────────────────────────────────────────────────
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
# GLM Evaluation by deviance stats and DHARMa
##################################################################################

# ─── 1. Deviance Statistics ──────────────────────────────────────────────────
calc_deviance_stats <- function(model, predictor, outcome, n) {
  tibble(
    predictor = predictor,
    outcome = outcome,
    n = n,
    deviance = model$deviance,
    df_residual = model$df.residual,
    deviance_ratio = model$deviance / model$df.residual,
    AIC = AIC(model)
  )
}

run_deviance_for_set <- function(data, predictors, predictor_set_name) {
  map_dfr(predictors, function(var) {
    map_dfr(outcomes, function(outcome) {
      df_pair <- data %>%
        select(all_of(c(outcome, var, fixed_effects))) %>%
        drop_na()
      
      if (nrow(df_pair) == 0) return(tibble())
      for (fe in fixed_effects) {
        if (is.factor(df_pair[[fe]]) && n_distinct(df_pair[[fe]]) < 2) return(tibble())
      }
      
      model <- glm(as.formula(paste(outcome, "~", var, "+", paste(fixed_effects, collapse = " + "))),
                   data = df_pair, family = Gamma(link = "log"))
      calc_deviance_stats(model, var, outcome, nrow(df_pair))
    })
  }) %>% mutate(predictor_set = predictor_set_name)
}

deviance_all <- bind_rows(
  run_deviance_for_set(df_cvd_and_lipidomics, lipid_predictors, "Lipids"),
  run_deviance_for_set(df_cvd_and_fatty_acids, fatty_acid_predictors, "Fatty acids"),
  run_deviance_for_set(df_cvd_and_risk_factors, risk_factor_predictors, "Risk factors"),
  run_deviance_for_set(df_cvd_and_body_comp, body_comp_predictors, "Body composition"),
  run_deviance_for_set(df_cvd_and_urine_nmr, urine_nmr_predictors, "Urine NMR"),
  run_deviance_for_set(df_cvd_and_REDcap, REDcap_numeric_preds, "REDCap numeric")
)

deviance_summary <- deviance_all %>%
  group_by(outcome) %>%
  summarise(
    n_models = n(),
    mean_dev_ratio = mean(deviance_ratio),
    median_dev_ratio = median(deviance_ratio),
    sd_dev_ratio = sd(deviance_ratio)
  ) %>%
  arrange(desc(mean_dev_ratio))

print("Deviance Summary by Outcome (Age-Adjusted):")
print(deviance_summary)

# ─── 2. DHARMa Diagnostics ───────────────────────────────────────────────────
test_dharma_per_outcome <- function() {
  pred_name <- risk_factor_predictors[1]
  
  map_dfr(outcomes, function(outcome) {
    df_pair <- df_cvd_and_risk_factors %>%
      select(all_of(c(outcome, pred_name)), Statins, Supplements, Gender, Country, z_Age) %>%
      drop_na()
    
    model <- glm(
      as.formula(paste(outcome, "~", pred_name, "+ Statins + Supplements + Gender + Country + z_Age")),
      data = df_pair, 
      family = Gamma(link = "log")
    )
    
    sim_resid <- simulateResiduals(model, n = 1000, plot = FALSE)
    
    tibble(
      outcome = outcome,
      predictor = pred_name,
      n = nrow(df_pair),
      uniformity_p = testUniformity(sim_resid, plot = FALSE)$p.value,
      dispersion_p = testDispersion(sim_resid, plot = FALSE)$p.value,
      dispersion_stat = testDispersion(sim_resid, plot = FALSE)$statistic,
      outlier_p = testOutliers(sim_resid, plot = FALSE)$p.value
    )
  })
}

dharma_results <- test_dharma_per_outcome()

print("DHARMa Diagnostics by Outcome (Age-Adjusted):")
print(dharma_results)

# ─── 3. Combined Summary Table ───────────────────────────────────────────────
diagnostics_summary <- deviance_summary %>%
  left_join(
    dharma_results %>% select(outcome, uniformity_p, dispersion_stat, dispersion_p, outlier_p),
    by = "outcome"
  ) %>%
  mutate(
    outcome_label = case_when(
      outcome == "QRISK3_risk" ~ "QRISK3",
      outcome == "SCORE2_score" ~ "SCORE2",
      outcome == "ascvd_10y" ~ "ASCVD",
      outcome == "frs_10y" ~ "Framingham"
    )
  ) %>%
  select(outcome_label, n_models, median_dev_ratio, 
         uniformity_p, dispersion_stat, dispersion_p, outlier_p)

print("Combined Model Diagnostics (Age-Adjusted):")
print(diagnostics_summary)

# Save results
write_csv(deviance_summary, "deviance_summary_age_adj.csv")
write_csv(dharma_results, "dharma_diagnostics_age_adj.csv")
write_csv(diagnostics_summary, "combined_diagnostics_summary_age_adj.csv")

# ─── Save Diagnostics ────────────────────────────────────────────────────────
saveRDS(deviance_all,         file.path(incremental_dir, "deviance_all_age_adj.rds"))
saveRDS(dharma_results,       file.path(incremental_dir, "dharma_results_age_adj.rds"))
saveRDS(diagnostics_summary,  file.path(incremental_dir, "diagnostics_summary_age_adj.rds"))
cat("Diagnostics saved.\n")