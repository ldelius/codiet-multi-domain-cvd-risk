### Cardiovascular Risk Scores and Predictor Associations: 
### Generalised Linear Models with Fixed Effects
### Author: Luisa Delius

# ─── Packages ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(broom)
library(patchwork)
library(flextable)
library(DHARMa)
library(readxl)

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

df_hypertension_treatment <- read_excel("../QRISK3_data.xlsx", sheet = "bp_medication") %>%
  mutate(Sample_ID = str_replace_all(Sample_ID, "-", "_"))

df_qrisk3_input <- readRDS("QRISK3_calculation_input.rds") %>%
  rename(Sample_ID = PatientID) %>%
  arrange(Sample_ID) %>%
  mutate(z_systolic = as.numeric(scale(systolic)))

# ─── 2. Prepare Merged DataFrames for GLM ────────────────────────────────────
df_fixed_effects <- df_risk_factors %>%
  select(Sample_ID, Gender, Country) %>%
  full_join(df_lipidomics_predictors %>% select(Sample_ID, Statins, Supplements), by = "Sample_ID")

df_cvd_and_fatty_acids <- df_all_cvd_risk_scores %>%
  full_join(df_fatty_acids_predictors %>% select(-QRISK3_risk), by = "Sample_ID") %>%
  full_join(df_fixed_effects %>% select(Sample_ID, Gender, Country), by = "Sample_ID")

df_cvd_and_lipidomics <- df_all_cvd_risk_scores %>%
  full_join(df_lipidomics_predictors %>% select(-QRISK3_risk), by = "Sample_ID") %>%
  full_join(df_fixed_effects %>% select(Sample_ID, Gender, Country), by = "Sample_ID")

df_cvd_and_REDcap <- df_all_cvd_risk_scores %>%
  full_join(df_REDcap_demographics %>% select(-QRISK3_risk, -recruitment_site), by = "Sample_ID") %>%
  full_join(df_fixed_effects, by = "Sample_ID")

df_cvd_and_risk_factors <- df_all_cvd_risk_scores %>%
  full_join(df_risk_factors %>% select(-QRISK3_risk), by = "Sample_ID") %>%
  full_join(df_fixed_effects %>% select(Sample_ID, Statins, Supplements), by = "Sample_ID")

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

outcomes <- c("QRISK3_risk", "SCORE2_score", "ascvd_10y", "frs_10y")

fixed_effects <- c("Statins", "Supplements", "Gender", "Country")

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
        if (n_distinct(df_pair[[fe]]) < 2) return(tibble())
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

# Save incremental results
incremental_dir <- file.path(wkdir, "glm_incremental_results")
dir.create(incremental_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(glm_lipid_num,     file.path(incremental_dir, "glm_lipid_num.rds"))
saveRDS(glm_fatty_num,     file.path(incremental_dir, "glm_fatty_num.rds"))
saveRDS(glm_urine_num,     file.path(incremental_dir, "glm_urine_num.rds"))
saveRDS(glm_REDcap_num,    file.path(incremental_dir, "glm_REDcap_num.rds"))
saveRDS(glm_risk_num,      file.path(incremental_dir, "glm_risk_num.rds"))
saveRDS(glm_body_comp_num, file.path(incremental_dir, "glm_body_comp_num.rds"))
saveRDS(glm_REDcap_fac,    file.path(incremental_dir, "glm_REDcap_fac.rds"))

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

saveRDS(all_results_glm, file.path(incremental_dir, "all_results_glm_combined.rds"))

# ─── 7. Display Helpers ──────────────────────────────────────────────────────
rename_predictors <- function(term_plot) {
  case_when(
    term_plot == "z_ecw_tbw"                      ~ "Extracellular water / Total body water",
    term_plot == "z_mean_hrt"                      ~ "Mean heart rate",
    term_plot == "z_ecm_bcm"                       ~ "Extracellular mass / Body cell mass",
    term_plot == "z_rbc_22_4n_6"                   ~ "RBC adrenic acid (22:4n-6)",
    term_plot == "z_rbc_dpa_22_5n3"                ~ "RBC docosapentaenoic acid (22:5n-3)",
    term_plot == "z_rbc_eicosadienoic_20_2n6"      ~ "RBC eicosadienoic acid (20:2n-6)",
    term_plot == "z_rbc_epa_20_5n3"                ~ "RBC eicosapentaenoic acid (20:5n-3)",
    term_plot == "z_pe_o_19_1_20_5"                ~ "PE-O (19:1/20:5)",
    term_plot == "naps_during_dayyes"              ~ "Daytime naps",
    term_plot == "Living_StatusLiving with a partner" ~ "Living with partner",
    term_plot == "Employment_StatusRetired"         ~ "Retired",
    term_plot == "z_AGE.reader"                    ~ "AGE Reader",
    term_plot == "z_Trigonelline"                   ~ "Trigonelline",
    term_plot == "z_Hba1C"                         ~ "HbA1c",
    term_plot == "z_ALT.unit.L"                    ~ "Alanine aminotransferase (ALT)",
    term_plot == "z_LDL.mg.dl"                     ~ "LDL",
    term_plot == "z_3_hydroxybutyric_acid"         ~ "3-Hydroxybutyric acid",
    term_plot == "z_hippuric_acid"                 ~ "Hippuric acid",
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

heatmap_data_filtered <- heatmap_data %>% filter(term_plot_clean %in% predictors_with_sig)
heatmap_continuous    <- heatmap_data_filtered %>% filter(!is_categorical)
heatmap_categorical   <- heatmap_data_filtered %>% filter(is_categorical)

# Order by number of significant associations, then effect size
term_order_continuous <- heatmap_continuous %>%
  filter(is_significant) %>%
  group_by(term_plot_clean, predictor_set) %>%
  summarise(n_sig = n(), mean_effect = mean(abs(log_effect)), .groups = "drop") %>%
  arrange(predictor_set, desc(n_sig), desc(mean_effect)) %>%
  pull(term_plot_clean)

term_order_categorical <- heatmap_categorical %>%
  filter(is_significant) %>%
  group_by(term_plot_clean) %>%
  summarise(n_sig = n(), mean_effect = mean(abs(log_effect)), .groups = "drop") %>%
  arrange(desc(n_sig), desc(mean_effect)) %>%
  pull(term_plot_clean)

# Continuous panel
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
  labs(x = NULL, y = "Continuous Predictors") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 11, face = "bold"),
    strip.text.y = element_text(size = 12, face = "bold", angle = 0, lineheight = 0.85),
    legend.title = element_text(size = 13), legend.text = element_text(size = 12),
    legend.position = "right", panel.grid = element_blank(),
    plot.margin = margin(5, 5, 0, 5), panel.spacing = unit(0.5, "lines")
  )

# Categorical panel
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

n_cont <- length(unique(heatmap_continuous$term_plot_clean))
n_cat  <- length(unique(heatmap_categorical$term_plot_clean))

combined_fdr <- p_continuous / p_categorical +
  plot_layout(heights = c(n_cont, n_cat)) +
  plot_annotation(
    title    = "GLM: Associations between biomarkers and CVD risk scores",
    subtitle = "* indicates FDR-adjusted p < 0.05",
    caption  = paste0(
      "Model: outcome ~ predictor + Statins + Supplements + Sex + Country (Gamma GLM, log link)\n",
      "Multiple testing correction: Benjamini-Hochberg | Separate colour scales for continuous and categorical predictors\n",
      "Reference levels: Daytime naps (vs no naps), Living with partner/family (vs living alone), Retired/Non-working (vs employed)"
    ),
    theme = theme(
      plot.title = element_text(size = 17, face = "bold"),
      plot.subtitle = element_text(size = 14),
      plot.caption = element_text(size = 11, hjust = 0)
    )
  )

plot_height <- (n_cont + n_cat) * 0.4 + 3
ggsave("heatmap_significant_predictors_split.png", combined_fdr,
       width = 12, height = plot_height, dpi = 300)

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

preds_3plus <- heatmap_exploratory %>%
  filter(is_nominal_sig) %>%
  group_by(term_plot_clean) %>%
  filter(n() >= 3) %>%
  pull(term_plot_clean) %>%
  unique()

heatmap_exp_filtered <- heatmap_exploratory %>% filter(term_plot_clean %in% preds_3plus)
heatmap_cont_exp <- heatmap_exp_filtered %>% filter(!is_categorical)
heatmap_cat_exp  <- heatmap_exp_filtered %>% filter(is_categorical)

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
  labs(x = if (nrow(heatmap_cat_exp) > 0) NULL else "Cardiovascular risk score",
       y = "Continuous Predictors") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = if (nrow(heatmap_cat_exp) > 0) element_blank() else element_text(size = 14, angle = 45, hjust = 1),
    axis.ticks.x = if (nrow(heatmap_cat_exp) > 0) element_blank() else element_line(),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 11, face = "bold"),
    strip.text.y = element_text(size = 12, face = "bold", angle = 0, lineheight = 0.85),
    legend.title = element_text(size = 13), legend.text = element_text(size = 12),
    legend.position = "right", panel.grid = element_blank(),
    plot.margin = margin(5, 5, 0, 5), panel.spacing = unit(0.5, "lines")
  )

if (nrow(heatmap_cat_exp) > 0) {
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
  
  n_cont_exp <- length(unique(heatmap_cont_exp$term_plot_clean))
  n_cat_exp  <- length(unique(heatmap_cat_exp$term_plot_clean))
  
  combined_exploratory <- p_cont_exp / p_cat_exp +
    plot_layout(heights = c(n_cont_exp, n_cat_exp))
  
  plot_height_exp <- (n_cont_exp + n_cat_exp) * 0.4 + 3
} else {
  combined_exploratory <- p_cont_exp
  n_cont_exp <- length(unique(heatmap_cont_exp$term_plot_clean))
  plot_height_exp <- n_cont_exp * 0.4 + 3
}

combined_exploratory <- combined_exploratory +
  plot_annotation(
    title    = "GLM: Exploratory Associations",
    subtitle = "● indicates predictors with ≥3 nominal significant associations (p < 0.05, FDR p > 0.05)",
    caption  = paste0(
      "Model: outcome ~ predictor + Statins + Supplements + Sex + Country (Gamma GLM, log link)\n",
      "Excludes FDR-significant predictors"
    ),
    theme = theme(
      plot.title = element_text(size = 17, face = "bold"),
      plot.subtitle = element_text(size = 14),
      plot.caption = element_text(size = 11, hjust = 0)
    )
  )

ggsave("exploratory_heatmap_3plus_scores_split.png", combined_exploratory,
       width = 12, height = plot_height_exp, dpi = 300)

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
      starts_with("RR_CI___frs"),    starts_with("p_val___frs"),    starts_with("p_adj___frs"),
    )
  
  df_wide %>%
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
}

save_as_docx(
  "Table S2: Lipid Species"    = create_wide_supp_table(all_results_glm, "Lipids"),
  "Table S3: Fatty Acids"      = create_wide_supp_table(all_results_glm, "Fatty acids"),
  "Table S4: Urine NMR"        = create_wide_supp_table(all_results_glm, "Urine NMR"),
  "Table S5: REDCap Numeric"   = create_wide_supp_table(all_results_glm, "REDCap numeric"),
  "Table S6: REDCap Factors"   = create_wide_supp_table(all_results_glm, "REDCap factors"),
  "Table S7: Risk Factors"     = create_wide_supp_table(all_results_glm, "Risk factors"),
  "Table S8: Body Composition" = create_wide_supp_table(all_results_glm, "Body composition"),
  path = "supplementary_tables_GLM_all.docx",
  pr_section = prop_section(page_size = page_size(orient = "landscape"))
)


################################################################################
# MODEL DIAGNOSTICS: Deviance and DHARMa Tests
################################################################################

# ─── 11. Deviance Statistics ─────────────────────────────────────────────────
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

# Generic function to compute deviance for a predictor set
run_deviance <- function(data, predictors, fixed_effects, outcomes, predictor_set_label) {
  fixed_part <- paste(fixed_effects, collapse = " + ")
  
  map_dfr(predictors, function(var) {
    map_dfr(outcomes, function(outcome) {
      df_pair <- data %>%
        select(all_of(c(outcome, var, fixed_effects))) %>%
        drop_na()
      
      if (nrow(df_pair) == 0) return(tibble())
      for (fe in fixed_effects) {
        if (n_distinct(df_pair[[fe]]) < 2) return(tibble())
      }
      
      model <- glm(as.formula(paste(outcome, "~", var, "+", fixed_part)),
                   data = df_pair, family = Gamma(link = "log"))
      calc_deviance_stats(model, var, outcome, nrow(df_pair))
    })
  }) %>%
    mutate(predictor_set = predictor_set_label)
}

deviance_all <- bind_rows(
  run_deviance(df_cvd_and_lipidomics,   lipid_predictors,       fixed_effects, outcomes, "Lipids"),
  run_deviance(df_cvd_and_fatty_acids,  fatty_acid_predictors,  fixed_effects, outcomes, "Fatty acids"),
  run_deviance(df_cvd_and_risk_factors, risk_factor_predictors, fixed_effects, outcomes, "Risk factors"),
  run_deviance(df_cvd_and_body_comp,    body_comp_predictors,   fixed_effects, outcomes, "Body composition"),
  run_deviance(df_cvd_and_urine_nmr,    urine_nmr_predictors,  fixed_effects, outcomes, "Urine NMR"),
  run_deviance(df_cvd_and_REDcap,       REDcap_numeric_preds,  fixed_effects, outcomes, "REDCap numeric")
)

deviance_summary <- deviance_all %>%
  group_by(outcome) %>%
  summarise(
    n_models = n(),
    mean_dev_ratio = mean(deviance_ratio),
    median_dev_ratio = median(deviance_ratio),
    sd_dev_ratio = sd(deviance_ratio),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_dev_ratio))

cat("\nDeviance Summary by Outcome:\n")
print(deviance_summary)

# ─── 12. DHARMa Diagnostics ─────────────────────────────────────────────────
dharma_results <- map_dfr(outcomes, function(outcome) {
  pred_name <- risk_factor_predictors[1]
  
  df_pair <- df_cvd_and_risk_factors %>%
    select(all_of(c(outcome, pred_name, fixed_effects))) %>%
    drop_na()
  
  model <- glm(
    as.formula(paste(outcome, "~", pred_name, "+", paste(fixed_effects, collapse = " + "))),
    data = df_pair, family = Gamma(link = "log")
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

cat("\nDHARMa Diagnostics by Outcome:\n")
print(dharma_results)

# ─── 13. Combined Diagnostics Summary ────────────────────────────────────────
diagnostics_summary <- deviance_summary %>%
  left_join(
    dharma_results %>% select(outcome, uniformity_p, dispersion_stat, dispersion_p, outlier_p),
    by = "outcome"
  ) %>%
  mutate(
    outcome_label = rename_outcomes(outcome)
  ) %>%
  select(outcome_label, n_models, median_dev_ratio,
         uniformity_p, dispersion_stat, dispersion_p, outlier_p)

cat("\nCombined Model Diagnostics:\n")
print(diagnostics_summary)

write_csv(deviance_summary, "deviance_summary.csv")
write_csv(dharma_results, "dharma_diagnostics.csv")
write_csv(diagnostics_summary, "combined_diagnostics_summary.csv")

saveRDS(deviance_all,         file.path(incremental_dir, "deviance_all.rds"))
saveRDS(dharma_results,       file.path(incremental_dir, "dharma_results.rds"))
saveRDS(diagnostics_summary,  file.path(incremental_dir, "diagnostics_summary.rds"))


################################################################################
# SENSITIVITY ANALYSES
################################################################################

# ─── 14. Build Master Dataframe for Sensitivity Analyses ─────────────────────
# Reload df_risk_factors without the column exclusions (need z_Age etc.)
df_risk_factors_full <- readRDS("df_risk_factor_predictors.rds") %>% arrange(Sample_ID)

df_master <- df_all_cvd_risk_scores %>%
  full_join(df_risk_factors_full, by = "Sample_ID") %>%
  full_join(df_lipidomics_predictors %>% select(Sample_ID, Statins, Supplements, starts_with("z_")),
            by = "Sample_ID") %>%
  full_join(df_REDcap_demographics %>% select(Sample_ID, starts_with("z_"), where(is.factor), -recruitment_site),
            by = "Sample_ID") %>%
  full_join(df_body_composition, by = "Sample_ID") %>%
  full_join(df_urine_nmr_data, by = "Sample_ID") %>%
  full_join(df_hypertension_treatment %>% select(Sample_ID, blood_pressure_treatment),
            by = "Sample_ID") %>%
  full_join(df_qrisk3_input %>% select(Sample_ID, systolic, z_systolic),
            by = "Sample_ID")

base_fixed <- c("Statins", "Supplements", "Gender", "Country")

# ─── 15. Sensitivity GLM Function ────────────────────────────────────────────
run_sensitivity_glm <- function(data, outcomes, base_fixed,
                                predictor_spec, extra_fixed,
                                analysis_label) {
  fixed_base_part  <- paste(base_fixed, collapse = " + ")
  fixed_extra_part <- paste(c(base_fixed, extra_fixed), collapse = " + ")
  
  var       <- predictor_spec$var
  var_type  <- predictor_spec$type
  term_lab  <- predictor_spec$term_label
  
  map_dfr(outcomes, function(outcome) {
    cols_needed <- unique(c(outcome, var, base_fixed, extra_fixed))
    df_pair <- data %>% select(all_of(cols_needed)) %>% drop_na()
    
    if (nrow(df_pair) < 10) return(tibble())
    for (fe in c(base_fixed, extra_fixed)) {
      if (is.factor(df_pair[[fe]]) || is.character(df_pair[[fe]])) {
        if (n_distinct(df_pair[[fe]]) < 2) return(tibble())
      }
    }
    
    m_base <- glm(as.formula(paste(outcome, "~", var, "+", fixed_base_part)),
                  data = df_pair, family = Gamma(link = "log"))
    m_adj  <- glm(as.formula(paste(outcome, "~", var, "+", fixed_extra_part)),
                  data = df_pair, family = Gamma(link = "log"))
    
    if (var_type == "numeric") {
      t_base <- tidy(m_base, conf.int = TRUE, exponentiate = TRUE) %>% filter(term == var)
      t_adj  <- tidy(m_adj,  conf.int = TRUE, exponentiate = TRUE) %>% filter(term == var)
    } else {
      t_base <- tidy(m_base, conf.int = TRUE, exponentiate = TRUE) %>% filter(str_starts(term, var))
      t_adj  <- tidy(m_adj,  conf.int = TRUE, exponentiate = TRUE) %>% filter(str_starts(term, var))
    }
    
    if (nrow(t_base) == 0 || nrow(t_adj) == 0) return(tibble())
    
    bind_cols(
      t_base %>%
        transmute(term, predictor_label = term_lab, outcome = outcome, n = nrow(df_pair),
                  RR_base = estimate, CI_low_base = conf.low, CI_high_base = conf.high, p_base = p.value),
      t_adj %>%
        transmute(RR_adj = estimate, CI_low_adj = conf.low, CI_high_adj = conf.high, p_adj = p.value)
    ) %>%
      mutate(
        delta_RR = RR_adj - RR_base,
        pct_change_RR = ((RR_adj - RR_base) / (RR_base - 1)) * 100,
        direction_base = if_else(RR_base > 1, "↑", "↓"),
        direction_adj  = if_else(RR_adj > 1, "↑", "↓"),
        direction_changed = direction_base != direction_adj,
        analysis = analysis_label
      )
  })
}

format_sensitivity_table <- function(results) {
  results %>%
    mutate(
      outcome_label = rename_outcomes(outcome),
      RR_CI_base = sprintf("%.3f (%.3f–%.3f)", RR_base, CI_low_base, CI_high_base),
      RR_CI_adj  = sprintf("%.3f (%.3f–%.3f)", RR_adj, CI_low_adj, CI_high_adj),
      p_base_fmt = ifelse(p_base < 0.001, sprintf("%.2e", p_base), sprintf("%.4f", p_base)),
      p_adj_fmt  = ifelse(p_adj < 0.001, sprintf("%.2e", p_adj), sprintf("%.4f", p_adj)),
      delta_RR_fmt = sprintf("%+.4f", delta_RR),
      direction_change = if_else(direction_changed, "YES", "")
    ) %>%
    select(analysis, predictor_label, outcome_label, n,
           RR_CI_base, p_base_fmt, RR_CI_adj, p_adj_fmt,
           delta_RR_fmt, direction_change)
}

# ─── 16. Run Sensitivity Analyses ────────────────────────────────────────────

# A1: Add Age
sens_hr_age <- run_sensitivity_glm(
  df_master, outcomes, base_fixed,
  list(var = "z_mean_hrt", type = "numeric", term_label = "Mean Heart Rate"),
  "z_Age", "+ Age"
)

sens_naps_age <- run_sensitivity_glm(
  df_master, outcomes, base_fixed,
  list(var = "naps_during_day", type = "factor", term_label = "Daytime Naps (yes vs no)"),
  "z_Age", "+ Age"
)

sens_retired_age <- run_sensitivity_glm(
  df_master, outcomes, base_fixed,
  list(var = "Employment_Status", type = "factor", term_label = "Employment Status"),
  "z_Age", "+ Age"
)

results_age <- bind_rows(sens_hr_age, sens_naps_age, sens_retired_age)

cat("\n═══ A1: Effect of adding AGE as fixed effect ═══\n")
print(format_sensitivity_table(results_age), n = Inf, width = Inf)

# A2: Add BP Medication (Heart Rate only)
sens_hr_bpmed <- run_sensitivity_glm(
  df_master, outcomes, base_fixed,
  list(var = "z_mean_hrt", type = "numeric", term_label = "Mean Heart Rate"),
  "blood_pressure_treatment", "+ BP Medication"
)

cat("\n═══ A2: Effect of adding BP MEDICATION ═══\n")
print(format_sensitivity_table(sens_hr_bpmed), n = Inf, width = Inf)

# A3: Body composition / BP adjustments (Heart Rate)
body_comp_adjustments <- list(
  list(extra = "z_smi__skeletal_muscle_index_",                       label = "+ Muscle Index"),
  list(extra = "z_svr_skeletal_muscle_mass_visceral_fat_area_ratio_", label = "+ Muscle/Fat Ratio"),
  list(extra = "z_vfa__visceral_fat_area_",                           label = "+ Visceral Fat"),
  list(extra = "z_whtr_waist_height_ratio_",                          label = "+ Waist-to-Height Ratio"),
  list(extra = "z_systolic",                                          label = "+ Systolic BP")
)

results_bodycomp_hr <- map_dfr(body_comp_adjustments, function(adj) {
  run_sensitivity_glm(
    df_master, outcomes, base_fixed,
    list(var = "z_mean_hrt", type = "numeric", term_label = "Mean Heart Rate"),
    adj$extra, adj$label
  )
})

cat("\n═══ A3: Body composition / BP adjustments (Heart Rate) ═══\n")
print(format_sensitivity_table(results_bodycomp_hr), n = Inf, width = Inf)

# A4: Visceral Fat — Heart Rate AND Trigonelline
sens_trigonelline_vfat <- run_sensitivity_glm(
  df_master, outcomes, base_fixed,
  list(var = "z_Trigonelline", type = "numeric", term_label = "Trigonelline"),
  "z_vfa__visceral_fat_area_", "+ Visceral Fat"
)

cat("\n═══ A4: Visceral Fat adjustment — HR and Trigonelline ═══\n")
print(format_sensitivity_table(
  bind_rows(results_bodycomp_hr %>% filter(analysis == "+ Visceral Fat"), sens_trigonelline_vfat)
), n = Inf, width = Inf)

# Save all sensitivity results
all_sensitivity <- bind_rows(results_age, sens_hr_bpmed, results_bodycomp_hr, sens_trigonelline_vfat)
write_csv(format_sensitivity_table(all_sensitivity), "sensitivity_analysis_results.csv")
cat("\nSaved: sensitivity_analysis_results.csv\n")


################################################################################
# ASSOCIATION PLOTS
################################################################################

# ─── 17. Plot Helpers ────────────────────────────────────────────────────────
theme_assoc <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )

plot_scatter <- function(data, x_var, y_var, x_label, y_label, title) {
  df_plot <- data %>% select(all_of(c(x_var, y_var))) %>% drop_na()
  names(df_plot) <- c("x", "y")
  
  cor_test <- cor.test(df_plot$x, df_plot$y, method = "pearson")
  r <- cor_test$estimate
  p <- cor_test$p.value
  p_label <- if (p < 0.001) sprintf("r = %.2f, p < 0.001", r) else sprintf("r = %.2f, p = %.3f", r, p)
  
  ggplot(df_plot, aes(x = x, y = y)) +
    geom_point(alpha = 0.4, size = 1.5, colour = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, colour = "#D55E00", linewidth = 0.8) +
    annotate("text", x = -Inf, y = Inf, label = p_label,
             hjust = -0.05, vjust = 1.5, size = 3.5, fontface = "italic") +
    labs(x = x_label, y = y_label, title = title) +
    theme_assoc
}

plot_boxplot <- function(data, x_var, y_var, x_label, y_label, title) {
  df_plot <- data %>% select(all_of(c(x_var, y_var))) %>% drop_na()
  names(df_plot) <- c("x", "y")
  df_plot$x <- factor(df_plot$x)
  
  n_groups <- n_distinct(df_plot$x)
  
  if (n_groups == 2) {
    test <- t.test(y ~ x, data = df_plot, var.equal = TRUE)
    test_label <- sprintf("Student's t-test p = %s",
                          if (test$p.value < 0.001) sprintf("%.2e", test$p.value) else sprintf("%.3f", test$p.value))
  } else {
    test <- summary(aov(y ~ x, data = df_plot))
    p_val <- test[[1]][["Pr(>F)"]][1]
    test_label <- sprintf("One-way ANOVA p = %s",
                          if (p_val < 0.001) sprintf("%.2e", p_val) else sprintf("%.3f", p_val))
  }
  
  n_per_group <- df_plot %>% count(x)
  label_map <- setNames(paste0(n_per_group$x, "\n(n=", n_per_group$n, ")"), n_per_group$x)
  df_plot$x_label <- factor(label_map[as.character(df_plot$x)], levels = label_map)
  
  ggplot(df_plot, aes(x = x_label, y = y)) +
    geom_boxplot(fill = "steelblue", alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, colour = "red") +
    annotate("text", x = -Inf, y = Inf, label = test_label,
             hjust = -0.05, vjust = 1.5, size = 3.5, fontface = "italic") +
    labs(x = x_label, y = y_label, title = title) +
    theme_assoc
}

# ─── 18. Heart Rate Association Panel (2×4) ──────────────────────────────────
combined_hr <- (
  plot_boxplot(df_master, "blood_pressure_treatment", "mean_hrt",
               "BP Medication", "Mean Heart Rate (bpm)", "Heart Rate vs BP Medication") |
    plot_boxplot(df_master, "Gender", "mean_hrt",
                 "Sex", "Mean Heart Rate (bpm)", "Heart Rate vs Sex") |
    plot_scatter(df_master, "Age", "mean_hrt",
                 "Age (years)", "Mean Heart Rate (bpm)", "Heart Rate vs Age") |
    plot_scatter(df_master, "systolic", "mean_hrt",
                 "Systolic BP (mmHg)", "Mean Heart Rate (bpm)", "Heart Rate vs Systolic BP")
) / (
  plot_scatter(df_master, "smi__skeletal_muscle_index_", "mean_hrt",
               "Skeletal Muscle Index", "Mean Heart Rate (bpm)", "Heart Rate vs Muscle Index") |
    plot_scatter(df_master, "svr_skeletal_muscle_mass_visceral_fat_area_ratio_", "mean_hrt",
                 "Muscle/Visceral Fat Ratio", "Mean Heart Rate (bpm)", "Heart Rate vs Muscle/Fat Ratio") |
    plot_scatter(df_master, "vfa__visceral_fat_area_", "mean_hrt",
                 "Visceral Fat Area", "Mean Heart Rate (bpm)", "Heart Rate vs Visceral Fat") |
    plot_scatter(df_master, "whtr_waist_height_ratio_", "mean_hrt",
                 "Waist-to-Height Ratio", "Mean Heart Rate (bpm)", "Heart Rate vs Waist-to-Height Ratio")
) +
  plot_annotation(
    title = "Associations with Mean Heart Rate",
    subtitle = "Boxplots: Student's t-test (red diamond = mean) | Scatter: Pearson correlation with linear fit",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

ggsave("heart_rate_associations_panel.png", combined_hr, width = 18, height = 10, dpi = 300)

# ─── 19. Age vs Naps and Employment ──────────────────────────────────────────
combined_age <- (
  plot_boxplot(df_master, "naps_during_day", "Age",
               "Daytime Naps", "Age (years)", "Age vs Daytime Naps") |
    plot_boxplot(df_master, "Employment_Status", "Age",
                 "Employment Status", "Age (years)", "Age vs Employment Status")
) +
  plot_annotation(
    title = "Age by Lifestyle / Sociodemographic Factors",
    subtitle = "Student's t-test (2 groups) or One-way ANOVA (>2 groups) | Red diamond = mean",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

ggsave("age_naps_employment_panel.png", combined_age, width = 14, height = 6, dpi = 300)

cat("\nAll outputs saved.\n")