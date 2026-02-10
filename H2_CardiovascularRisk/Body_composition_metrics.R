### Preparation of Body Composition metrics for CVD Risk Score Associations
### Author: Luisa Delius

# Load packages
library(tidyverse)
library(broom)
library(patchwork)
library(readxl)

# Set working directory
wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

# upload data
df_body_composition_metrics <- read_excel("Body_composition_matrix_biosensorsMicrcaya.xlsx")
QRISK3_sample_ID <- readRDS("QRISK3_sample_ID.RDS")

# Prepare the data for later CVD Risk Score Calculations
df_body_composition_metrics <- df_body_composition_metrics %>%
  rename_with(~ str_replace_all(.x, "[^a-zA-Z0-9_]", "_")) %>%
  rename(Sample_ID = volunteer_id) %>%
  mutate(Sample_ID = toupper(Sample_ID)) %>%
  mutate(Sample_ID = str_replace_all(Sample_ID, "\\s+", "")) %>%  # FIX SPACING HERE!
  group_by(Sample_ID) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.x)), .names = "z_{col}")) %>%
  left_join(QRISK3_sample_ID, by = "Sample_ID")

# save as RDS
saveRDS(df_body_composition_metrics, file.path(wkdir, "processed_data", "df_body_composition_metrics.rds"))

# create predictor
body_composition_predictors <- df_body_composition_metrics %>%
  select(starts_with("z_")) %>%
  names()
body_composition_predictors

# run glm
glm_output_body_composition <- map_dfr(body_composition_predictors, function(var) {
  df_pair <- df_body_composition_metrics %>%
    select(QRISK3_risk, all_of(var)) %>%
    drop_na()
  
  model <- glm(
    as.formula(paste("QRISK3_risk ~", var)),
    data   = df_pair,
    family = Gamma(link = "log")
  )
  
  tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      predictor = var,
      n = nrow(df_pair)
    )
}) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    significant = if_else(conf.low > 1 | conf.high < 1, "Significant", "Not significant")
  ) %>%
  arrange(p.adj)

glm_output_body_composition

# plot glm resultord_body_comp <- glm_output_body_composition %>%
ord_body_comp <- glm_output_body_composition %>%
  arrange(estimate) %>%
  pull(predictor)

ggplot(glm_output_body_composition,
       aes(x = estimate, y = factor(predictor, levels = ord_body_comp), colour = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point() +
  scale_colour_manual(values = c("Significant" = "firebrick3",
                                 "Not significant" = "black")) +
  labs(
    x = "Ratio of expected QRISK3 (exp(beta))",
    y = "Predictor",
    title = "GLM (95% CIs) with body composition metrics"
  ) +
  theme(axis.text.y = element_text(size = 6))
