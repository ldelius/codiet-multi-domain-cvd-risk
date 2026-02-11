### Cardiovascular Risk Scores and Predictor Associations: 
### Generalised Linear Models with Fixed Effects
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

outcomes <- c("QRISK3_risk", "SCORE2_score", "ascvd_10y", "frs_10y", "mean_risk")

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

# ─── 7. Plotting Helpers ─────────────────────────────────────────────────────
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
    outcome == "mean_risk"     ~ "Composite Score",
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

# Combine
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

# Keep predictors with ≥3 nominally significant associations
preds_3plus <- heatmap_exploratory %>%
  filter(is_nominal_sig) %>%
  group_by(term_plot_clean) %>%
  filter(n() >= 3) %>%
  pull(term_plot_clean) %>%
  unique()

heatmap_exp_filtered <- heatmap_exploratory %>% filter(term_plot_clean %in% preds_3plus)
heatmap_cont_exp <- heatmap_exp_filtered %>% filter(!is_categorical)
heatmap_cat_exp  <- heatmap_exp_filtered %>% filter(is_categorical)

# Order
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

# Combine (with or without categorical panel)
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
      "Excludes FDR-significant predictors | Separate scales for continuous and categorical"
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
# Wide-format table for each predictor set
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
  
  # Reorder columns by outcome
  df_wide <- df_wide %>%
    select(
      term_clean,
      starts_with("RR_CI___QRISK3"), starts_with("p_val___QRISK3"), starts_with("p_adj___QRISK3"),
      starts_with("RR_CI___SCORE2"), starts_with("p_val___SCORE2"), starts_with("p_adj___SCORE2"),
      starts_with("RR_CI___ascvd"),  starts_with("p_val___ascvd"),  starts_with("p_adj___ascvd"),
      starts_with("RR_CI___frs"),    starts_with("p_val___frs"),    starts_with("p_adj___frs"),
      starts_with("RR_CI___mean"),   starts_with("p_val___mean"),   starts_with("p_adj___mean")
    )
  
  ft <- df_wide %>%
    flextable() %>%
    set_header_labels(
      term_clean = "Predictor",
      RR_CI___QRISK3_risk = "RR(95%CI)", p_val___QRISK3_risk = "p", p_adj___QRISK3_risk = "p.adj",
      RR_CI___SCORE2_score = "RR(95%CI)", p_val___SCORE2_score = "p", p_adj___SCORE2_score = "p.adj",
      RR_CI___ascvd_10y = "RR(95%CI)", p_val___ascvd_10y = "p", p_adj___ascvd_10y = "p.adj",
      RR_CI___frs_10y = "RR(95%CI)", p_val___frs_10y = "p", p_adj___frs_10y = "p.adj",
      RR_CI___mean_risk = "RR(95%CI)", p_val___mean_risk = "p", p_adj___mean_risk = "p.adj"
    ) %>%
    add_header_row(
      values = c("", "QRISK3", "SCORE2", "ASCVD", "Framingham", "Composite"),
      colwidths = c(1, 3, 3, 3, 3, 3)
    ) %>%
    theme_booktabs() %>%
    width(j = 1, width = 1.5) %>%
    width(j = c(2, 5, 8, 11, 14), width = 1.0) %>%
    width(j = c(3, 6, 9, 12, 15), width = 0.5) %>%
    width(j = c(4, 7, 10, 13, 16), width = 0.5) %>%
    align(j = 1, align = "left", part = "all") %>%
    align(j = 2:ncol(df_wide), align = "center", part = "all") %>%
    bold(part = "header") %>%
    merge_h(part = "header") %>%
    vline(j = c(1, 4, 7, 10, 13), border = fp_border(width = 1.5), part = "all") %>%
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
  "Table S2: Lipid Species"    = ft_lipids,
  "Table S3: Fatty Acids"      = ft_fatty_acids,
  "Table S4: Urine NMR"        = ft_urine_nmr,
  "Table S5: REDCap Numeric"   = ft_redcap_num,
  "Table S6: REDCap Factors"   = ft_redcap_fac,
  "Table S7: Risk Factors"     = ft_risk_factors,
  "Table S8: Body Composition" = ft_body_comp,
  path = "supplementary_tables_GLM_all.docx",
  pr_section = prop_section(page_size = page_size(orient = "landscape"))
)

##################################################################################
# GLM Evaluation by devaince stas and DHARMa
##################################################################################
# ═══════════════════════════════════════════════════════════════════════════
# MODEL DIAGNOSTICS: Deviance and DHARMa Tests
# ═══════════════════════════════════════════════════════════════════════════

# ─── 1. Deviance Statistics ──────────────────────────────────────────────────
# Function to calculate deviance statistics
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

# Run for all predictor sets
deviance_all <- bind_rows(
  map_dfr(lipid_predictors, function(var) {
    map_dfr(outcomes, function(outcome) {
      df_pair <- df_cvd_and_lipidomics %>%
        select(all_of(c(outcome, var, fixed_effects))) %>%
        drop_na()
      
      if (nrow(df_pair) == 0) return(tibble())
      
      # Check each fixed effect has at least 2 levels
      for (fe in fixed_effects) {
        if (n_distinct(df_pair[[fe]]) < 2) return(tibble())
      }
      
      model <- glm(as.formula(paste(outcome, "~", var, "+", paste(fixed_effects, collapse = " + "))),
                   data = df_pair, family = Gamma(link = "log"))
      calc_deviance_stats(model, var, outcome, nrow(df_pair))
    })
  }) %>% mutate(predictor_set = "Lipids"),
  
  map_dfr(fatty_acid_predictors, function(var) {
    map_dfr(outcomes, function(outcome) {
      df_pair <- df_cvd_and_fatty_acids %>%
        select(all_of(c(outcome, var, fixed_effects))) %>%
        drop_na()
      
      if (nrow(df_pair) == 0) return(tibble())
      
      # Check each fixed effect has at least 2 levels
      for (fe in fixed_effects) {
        if (n_distinct(df_pair[[fe]]) < 2) return(tibble())
      }
      
      model <- glm(as.formula(paste(outcome, "~", var, "+", paste(fixed_effects, collapse = " + "))),
                   data = df_pair, family = Gamma(link = "log"))
      calc_deviance_stats(model, var, outcome, nrow(df_pair))
    })
  }) %>% mutate(predictor_set = "Fatty acids"),
  
  map_dfr(risk_factor_predictors, function(var) {
    map_dfr(outcomes, function(outcome) {
      df_pair <- df_cvd_and_risk_factors %>%
        select(all_of(c(outcome, var, fixed_effects))) %>%
        drop_na()
      
      if (nrow(df_pair) == 0) return(tibble())
      
      # Check each fixed effect has at least 2 levels
      for (fe in fixed_effects) {
        if (n_distinct(df_pair[[fe]]) < 2) return(tibble())
      }
      
      model <- glm(as.formula(paste(outcome, "~", var, "+", paste(fixed_effects, collapse = " + "))),
                   data = df_pair, family = Gamma(link = "log"))
      calc_deviance_stats(model, var, outcome, nrow(df_pair))
    })
  }) %>% mutate(predictor_set = "Risk factors"),
  
  map_dfr(body_comp_predictors, function(var) {
    map_dfr(outcomes, function(outcome) {
      df_pair <- df_cvd_and_body_comp %>%
        select(all_of(c(outcome, var, fixed_effects))) %>%
        drop_na()
      
      if (nrow(df_pair) == 0) return(tibble())
      
      for (fe in fixed_effects) {
        if (n_distinct(df_pair[[fe]]) < 2) return(tibble())
      }
      
      model <- glm(as.formula(paste(outcome, "~", var, "+", paste(fixed_effects, collapse = " + "))),
                   data = df_pair, family = Gamma(link = "log"))
      calc_deviance_stats(model, var, outcome, nrow(df_pair))
    })
  }) %>% mutate(predictor_set = "Body composition"),
  
  map_dfr(urine_nmr_predictors, function(var) {
    map_dfr(outcomes, function(outcome) {
      df_pair <- df_cvd_and_urine_nmr %>%
        select(all_of(c(outcome, var, fixed_effects))) %>%
        drop_na()
      
      if (nrow(df_pair) == 0) return(tibble())
      
      for (fe in fixed_effects) {
        if (n_distinct(df_pair[[fe]]) < 2) return(tibble())
      }
      
      model <- glm(as.formula(paste(outcome, "~", var, "+", paste(fixed_effects, collapse = " + "))),
                   data = df_pair, family = Gamma(link = "log"))
      calc_deviance_stats(model, var, outcome, nrow(df_pair))
    })
  }) %>% mutate(predictor_set = "Urine NMR"),
  
  map_dfr(REDcap_numeric_preds, function(var) {
    map_dfr(outcomes, function(outcome) {
      df_pair <- df_cvd_and_REDcap %>%
        select(all_of(c(outcome, var, fixed_effects))) %>%
        drop_na()
      
      if (nrow(df_pair) == 0) return(tibble())
      
      for (fe in fixed_effects) {
        if (n_distinct(df_pair[[fe]]) < 2) return(tibble())
      }
      
      model <- glm(as.formula(paste(outcome, "~", var, "+", paste(fixed_effects, collapse = " + "))),
                   data = df_pair, family = Gamma(link = "log"))
      calc_deviance_stats(model, var, outcome, nrow(df_pair))
    })
  }) %>% mutate(predictor_set = "REDCap numeric")
)

# Summary by outcome
deviance_summary <- deviance_all %>%
  group_by(outcome) %>%
  summarise(
    n_models = n(),
    mean_dev_ratio = mean(deviance_ratio),
    median_dev_ratio = median(deviance_ratio),
    sd_dev_ratio = sd(deviance_ratio)
  ) %>%
  arrange(desc(mean_dev_ratio))

print("Deviance Summary by Outcome:")
print(deviance_summary)

# ─── 2. DHARMa Diagnostics ───────────────────────────────────────────────────

# Test one model per outcome (using first risk factor predictor)
test_dharma_per_outcome <- function() {
  pred_name <- risk_factor_predictors[1]  # z_AGE.reader
  
  map_dfr(outcomes, function(outcome) {
    df_pair <- df_cvd_and_risk_factors %>%
      select(all_of(c(outcome, pred_name)), Statins, Supplements, Gender, Country) %>%
      drop_na()
    
    model <- glm(
      as.formula(paste(outcome, "~", pred_name, "+ Statins + Supplements + Gender + Country")),
      data = df_pair, 
      family = Gamma(link = "log")
    )
    
    # DHARMa tests
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

print("DHARMa Diagnostics by Outcome:")
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
      outcome == "frs_10y" ~ "Framingham",
      outcome == "mean_risk" ~ "Composite"
    )
  ) %>%
  select(outcome_label, n_models, median_dev_ratio, 
         uniformity_p, dispersion_stat, dispersion_p, outlier_p)

print("Combined Model Diagnostics:")
print(diagnostics_summary)

# Save results
write_csv(deviance_summary, "deviance_summary.csv")
write_csv(dharma_results, "dharma_diagnostics.csv")
write_csv(diagnostics_summary, "combined_diagnostics_summary.csv")



##################################################################################
##################################################################################
# 7. Investigate why a high Hear Rate is associated with low CVD risk (thats wrong)
# 7.1 is bp treatment responsible 
df_hypertension_treatment <- read_excel("../QRISK3_data.xlsx",sheet = "bp_medication") %>%
  mutate(Sample_ID = str_replace_all(Sample_ID, "-", "_"))

# For z_Heart.Rate from df_cvd_scores_and_risk_factors
map_dfr(outcomes, function(outcome) {
  df_test <- df_cvd_scores_and_risk_factors %>%
    left_join(df_hypertension_treatment %>% select(Sample_ID, blood_pressure_treatment), 
              by = "Sample_ID") %>%
    select(all_of(outcome), z_Heart.Rate, blood_pressure_treatment, 
           Statins, Supplements, Gender, Country) %>%
    drop_na()
  
  # Model WITHOUT medication adjustment
  model_no_med <- glm(as.formula(paste(outcome, "~ z_Heart.Rate + Statins + Supplements + Gender + Country")),
                      data = df_test, family = Gamma(link = "log"))
  
  # Model WITH medication adjustment
  model_with_med <- glm(as.formula(paste(outcome, "~ z_Heart.Rate + blood_pressure_treatment + Statins + Supplements + Gender + Country")),
                        data = df_test, family = Gamma(link = "log"))
  
  tibble(
    outcome = outcome,
    HR_coef_no_med = coef(model_no_med)["z_Heart.Rate"],
    HR_coef_with_med = coef(model_with_med)["z_Heart.Rate"],
    HR_p_no_med = summary(model_no_med)$coefficients["z_Heart.Rate", "Pr(>|t|)"], # i previously removed Heart Rate so code will stop running here..
    HR_p_with_med = summary(model_with_med)$coefficients["z_Heart.Rate", "Pr(>|t|)"]
  )
})


# plotting bp treatment vs hear rate
plot_data <- df_cvd_scores_and_risk_factors %>%
  left_join(df_hypertension_treatment %>% select(Sample_ID, blood_pressure_treatment), 
            by = "Sample_ID") %>%
  drop_na(blood_pressure_treatment, Heart.Rate) %>%
  mutate(blood_pressure_treatment = factor(blood_pressure_treatment, levels = c("No", "Yes")))

# Get stats
stats <- plot_data %>%
  group_by(blood_pressure_treatment) %>%
  summarise(n = n(), mean = mean(Heart.Rate))

# Student's t-test (equal variances assumed)
tt <- t.test(Heart.Rate ~ blood_pressure_treatment, data = plot_data, var.equal = TRUE)

ggplot(plot_data, aes(x = blood_pressure_treatment, y = Heart.Rate)) +
  geom_boxplot(fill = "steelblue", alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.2) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  annotate("text", x = 1.5, y = max(plot_data$Heart.Rate) * 1.05,
           label = paste0("Student's t-test: p = ", signif(tt$p.value, 3))) +
  scale_x_discrete(labels = paste0(c("No", "Yes"), "\n(n = ", stats$n, ")")) +
  labs(x = "Blood pressure treatment", y = "Mean heart rate") +
  theme_minimal()

# 7.2 is fittness responsible?