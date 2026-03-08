### Cardiovascular Risk Score Calculations
### Author: Luisa Delius

# this file is to calculate several cardiovascular risk scores
# 1. Framingham score
# 2. ASCVD
# 3. SCORE2
# 4. Combine the results

# Load packages
library(tidyverse)
library(broom)
library(patchwork)
library(QRISK3)
library(readxl)
library(CVrisk)
library(RiskScorescvd)
library(flextable)

# Set working directory
wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

# checking the various risk factor functions
?ascvd_10y_frs
?ascvd_10y_accaha
?SCORE2_scores
?QRISK3_2017

# Upload the data/excel sheet
excel_sheets("cardiovascular_riskfactor_calculations.xlsx")
df_maindata <- read_excel("cardiovascular_riskfactor_calculations.xlsx", sheet = "Sheet1")
df_Cholesterol_HDL  <- read_excel("cardiovascular_riskfactor_calculations.xlsx", sheet = "Cholesterol_HDL")
df_Blood_Pressure     <- read_excel("cardiovascular_riskfactor_calculations.xlsx", sheet = "Blood_Pressure")
df_age_risk_factor_sheet <- read_excel("cardiovascular_riskfactor_calculations.xlsx", sheet = "Age_risk_factor_sheet")
df_risk_region <- read_excel("cardiovascular_riskfactor_calculations.xlsx", sheet = "Risk_region")
QRISK3_sample_ID <- readRDS("QRISK3_sample_ID.RDS") ### this has to be exchanged if i change anything of the QRISK input data of courese!!!
QRISK3_calculation_input <- readRDS("processed_data/QRISK3_calculation_input.rds")
df_statins_supplements <- read_excel("Statins_Supplements.xlsx")
df_risk_factor_sheet <-readRDS("processed_data/df_risk_factor_predictors.rds")
  
# cleaning the data, so that it can be used for risk score calculations
## working on blood_pressure: change to SampleID, - to _, ...
blood_pressure_keep <- df_Blood_Pressure %>%
  rename(PatientID = patient) %>%   # make the ID name consistent
  mutate(
    PatientID = str_replace_all(PatientID, "-", "_")  # CD-003 → CD_003
  ) %>%
  select(
    PatientID, systolic
  )

## Working on cholesterol and HDL
lipids_keep <- df_Cholesterol_HDL %>%
  mutate(
    # Create patient-level ID (remove visit suffix _V1/_V2/_V3)
    PatientID = str_remove(SampleID, "_V[1-3]$"),
    Total_Cholesterol_mg_dl = parse_number(Total_Cholesterol_mg_dl), # parse_number() turns character columns in numeric columns
    HDL_mg_dl               = parse_number(HDL_mg_dl),
    LDL_mg_dl               = parse_number(LDL_mg_dl)
  ) %>%
  group_by(PatientID) %>%
  summarise( # calculates the mean each and the ratio of the means
    mean_Total_Cholesterol_mg_dl = mean(Total_Cholesterol_mg_dl, na.rm = TRUE),
    mean_HDL_mg_dl = mean(HDL_mg_dl, na.rm = TRUE),
    mean_LDL_mg_dl = mean(LDL_mg_dl, na.rm = TRUE),
        .groups = "drop" # removes all grouping (the data got grouped by (group_by() to calculate each patients mean))
  )

## working with patient age
patient_age <- df_age_risk_factor_sheet %>%
  mutate(PatientID = str_remove(SampleID, "_V[1-3]+$")) %>%
  group_by(PatientID) %>%
  summarise(Age = mean(Age, na.rm = TRUE)) %>%
  ungroup() # not needed, but i will keep it after every group to not forget about it.

## working with the main data sheet
maindata_keep <- df_maindata %>%
  rename(PatientID = volunteer_id) %>%   # make the ID name consistent
  mutate(
    PatientID = str_replace_all(PatientID, "-", "_"),  # CD-003 → CD_003
    blood_pressure_treatment = ifelse(blood_pressure_treatment == "Yes", 1, 0),
    SmokingStatusQRISK3 = ifelse(SmokingStatusQRISK3 %in% c(1, 2), 0, 1),
    Diabetes = ifelse(Diabetes == "Yes", 1, 0),
    race_ascvd = case_when(
      EthnicityCodeQRISK3 == 1 ~ "white",          # White or not stated
      EthnicityCodeQRISK3 %in% c(6, 7) ~ "aa",     # Black Caribbean + Black African
      TRUE ~ "other")
  ) %>%
  select(
    PatientID, Sex, race_ascvd, SmokingStatusQRISK3, Diabetes, blood_pressure_treatment 
  )

risk_region_keep <- df_risk_region %>%
  rename(PatientID = volunteer_id) %>%   # make the ID name consistent
  mutate(
    PatientID = str_replace_all(PatientID, "-", "_"),
    Risk.region = case_when(
      recruitment_site %in% c("CBG", "ICL", "UVEG") ~ "low",
      recruitment_site %in% c("CORK", "AUTH") ~ "moderate",
      TRUE ~ NA_character_
    )
  ) %>%
  select(PatientID, Risk.region)    # keep only the relevant columns

## combining the cleaned data
risk_factor_input <- maindata_keep %>%
  full_join(lipids_keep, by = "PatientID") %>%
  full_join(blood_pressure_keep, by = "PatientID") %>%
  full_join(patient_age, by = "PatientID") %>%
  full_join(risk_region_keep, by = "PatientID")

saveRDS(risk_factor_input, file.path(wkdir, "processed_data", "ASCVD_SCORE2_Framingham_input.rds"))

###################### Risk score calulations ######################
## 1. Framingham Score
Framingham <- data.frame(
  PatientID = risk_factor_input$PatientID,
  gender = tolower(risk_factor_input$Sex),
  age = risk_factor_input$Age,
  hdl = risk_factor_input$mean_HDL_mg_dl,
  totchol = risk_factor_input$mean_Total_Cholesterol_mg_dl,
  sbp = risk_factor_input$systolic,
  bp_med = risk_factor_input$blood_pressure_treatment,
  smoker = risk_factor_input$SmokingStatusQRISK3,
  diabetes = risk_factor_input$Diabetes
) %>%
  filter(!if_any(everything(), is.na)) %>% # cannot have NA values
  filter(age >= 30 & age <= 74) %>%  # Framingham age range
  mutate(
    frs_10y = ascvd_10y_frs(
      gender, age, hdl, totchol, sbp, bp_med, smoker, diabetes
    )
  ) # must be between 30 and 74 for Framingham


### checking distribution and mean of Framingham score
# Histogram + density for Framingham scores
plot_framigham <- ggplot(Framingham, aes(x = frs_10y)) +
  geom_histogram( # create the normal histogram
    bins = 20,
    fill = "#8AAD40",
    color = "white"
  ) + 
  labs( # create the titles 
    title = NULL,
    subtitle = NULL,
    x = "10-year cardiovascular risk (%)",
    y = "Number of participants"
  ) +
  geom_density( # distribution
    aes(y = ..density.. * nrow(Framingham) * 
        diff(range(Framingham$frs_10y, na.rm = TRUE)) / 20),
    color = "darkblue",
    size = 1.2
  ) +
  geom_vline(
    aes(xintercept = mean(frs_10y, na.rm = TRUE)), # mean Framingham value
    color = "red", linetype = "dashed", size = 1
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
plot_framigham

# Statistics for the Framingham distribution
shapiro.test(Framingham$frs_10y)  # normality test
summary(Framingham$frs_10y)

##---------------------------------------------------------------------##
## 2. ASCVD Score
# at the moment adjusted to exclude LDL otuside 70-189 mg/dl
ASCVD <- data.frame(
  PatientID = risk_factor_input$PatientID,
  race    = risk_factor_input$race_ascvd,   # "white", "aa", "other"
  gender  = tolower(risk_factor_input$Sex),
  age     = risk_factor_input$Age,
  hdl     = risk_factor_input$mean_HDL_mg_dl,
  totchol = risk_factor_input$mean_Total_Cholesterol_mg_dl,
  sbp     = risk_factor_input$systolic,
  bp_med  = risk_factor_input$blood_pressure_treatment,
  smoker  = risk_factor_input$SmokingStatusQRISK3,
  diabetes = risk_factor_input$Diabetes,
  ldl     = risk_factor_input$mean_LDL_mg_dl   # <- included, but NOT used in score
) %>%
  filter(!if_any(everything(), is.na)) %>%    # cannot have NA values
  filter(
    age >= 40 & age <= 79,  # age range according to website
    ldl >= 70 & ldl <= 189 # LDL according to website
    ) %>%          
  mutate(
    ascvd_10y = ascvd_10y_accaha(
      race, gender, age, totchol, hdl, sbp, bp_med, smoker, diabetes
    )
  )

# just checking the people that got excluded because of their values outside specific ranges
#ASCVD %>%
#  filter(is.na(ascvd_10y)) %>%
#  select(PatientID, age, totchol, hdl, sbp)

# checking how many people we loose because of LDL exclusion
#excluded <- risk_factor_input %>%
#  filter(is.na(mean_LDL_mg_dl) |
#           mean_LDL_mg_dl < 70 |
#           mean_LDL_mg_dl > 189) %>%
#  select(PatientID, mean_LDL_mg_dl)
#excluded

### checking distribution and mean of ASCVD score
# Histogram + density for ASCVD scores
plot_ASCVD <- ggplot(ASCVD, aes(x = ascvd_10y)) +
  geom_histogram( # create the normal histogram
    bins = 20,
    fill = "#F8766D",
    color = "white"
  ) + 
  labs( # create the titles 
    title = NULL,
    subtitle = NULL,
    x = "10-year cardiovascular risk (%)",
    y = "Number of participants"
  ) +
  geom_density( # distribution
    aes(y = ..density.. * nrow(ASCVD) * 
          diff(range(ASCVD$ascvd_10y, na.rm = TRUE)) / 20),
    color = "darkblue",
    size = 1.2
  ) +
  geom_vline(
    aes(xintercept = mean(ascvd_10y, na.rm = TRUE)), # mean ASCVD value
    color = "red", linetype = "dashed", size = 1
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
plot_ASCVD

# Statistics for the ASCVD distribution
shapiro.test(ASCVD$ascvd_10y)  # normality test
summary(ASCVD$ascvd_10y)
##---------------------------------------------------------------------##
## 3. SCORE2
SCORE2_input <- data.frame(
  PatientID   = risk_factor_input$PatientID,
  Gender  = tolower(risk_factor_input$Sex),
  Age     = risk_factor_input$Age,
  smoker  = risk_factor_input$SmokingStatusQRISK3,
  systolic.bp = risk_factor_input$systolic,
  diabetes = risk_factor_input$Diabetes,
  total.chol = risk_factor_input$mean_Total_Cholesterol_mg_dl * 0.02586,
  total.hdl	= risk_factor_input$mean_HDL_mg_dl * 0.02586,
  region  = risk_factor_input$Risk.region
) %>%
  filter(!if_any(everything(), is.na)) %>%
  filter(Age >= 40 & Age <= 69)

# Calculate for LOW region
df_low <- SCORE2_input %>% filter(region == "low")

scores_low <- SCORE2_scores(
  data        = df_low %>% select(Age, Gender, smoker,
                                  systolic.bp, diabetes,
                                  total.chol, total.hdl),
  Risk.region = "Low",
  classify    = TRUE
)

df_low_res <- bind_cols(df_low, scores_low)

# Calculate for MODERATE region
df_mod <- SCORE2_input %>% filter(region == "moderate")

scores_mod <- SCORE2_scores(
  data        = df_mod %>% select(Age, Gender, smoker,
                                  systolic.bp, diabetes,
                                  total.chol, total.hdl),
  Risk.region = "Moderate",  # MUST be capitalized
  classify    = TRUE
)

df_mod_res <- bind_cols(df_mod, scores_mod)

# combine the scores and PatientID
SCORE2 <- bind_rows(df_low_res, df_mod_res) %>%
  select(PatientID, starts_with("SCORE2"))

### checking distribution and mean of SCORE2 
# Histogram + density for SCORE2
plot_SCORE2 <- ggplot(SCORE2, aes(x = SCORE2_score)) +
  geom_histogram(
    bins = 20,
    fill = "#C77CFF",
    color = "white"
  ) + 
  labs(
    title = NULL,
    subtitle = NULL,
    x = "10-year cardiovascular risk (%)",
    y = "Number of participants"
  ) +
  geom_density(
    aes(y = ..density.. * nrow(SCORE2) *
          diff(range(SCORE2$SCORE2_score, na.rm = TRUE)) / 20),
    color = "darkblue",
    size = 1.2
  ) +
  geom_vline(
    aes(xintercept = mean(SCORE2_score, na.rm = TRUE)),
    color = "red", linetype = "dashed", size = 1
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
plot_SCORE2

# Statistics for the SCORE2 distribution
shapiro.test(SCORE2$SCORE2_score)   # normality test
summary(SCORE2$SCORE2_score)        # summary stats
##---------------------------------------------------------------------##
# Plot QRISK3 dsitribtuion
# QRISK3
plot_qrisk3 <- ggplot(QRISK3_sample_ID, aes(x = QRISK3_risk)) +
  geom_histogram(
    bins = 20,
    fill = "#4DB8B0",
    color = "white"
  ) + 
  labs(
    title = NULL,
    subtitle = NULL,
    x = "10-year cardiovascular risk (%)",
    y = "Number of participants"
  ) +
  geom_density(
    aes(y = ..density.. * nrow(QRISK3_sample_ID) * 
          diff(range(QRISK3_sample_ID$QRISK3_risk)) / 20),
    color = "darkblue",
    size = 1.2
  ) +
  geom_vline(
    aes(xintercept = mean(QRISK3_risk, na.rm = TRUE)),
    color = "red", linetype = "dashed", size = 1
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
plot_qrisk3
##---------------------------------------------------------------------##
# 4. Combine the results
## create a joint data frame
df_all_risk_scores <- SCORE2 %>%
  full_join(ASCVD %>% select(PatientID, ascvd_10y), by = "PatientID") %>%
  full_join(Framingham %>% select(PatientID, frs_10y),by = "PatientID") %>%
  full_join(QRISK3_sample_ID %>% rename(PatientID = Sample_ID), by = "PatientID") %>%
  mutate(mean_risk = rowMeans( #calculate an average value for each patient --> composite score
    select(., SCORE2_score, ascvd_10y, frs_10y, QRISK3_risk),
    na.rm = TRUE
  ))

saveRDS(df_all_risk_scores, file = file.path(wkdir, "processed_data", "df_all_risk_scores.rds"))

## histogram of the composite risk score
plot_composite <- ggplot(df_all_risk_scores, aes(x = mean_risk)) +
  geom_histogram(
    bins = 20,
    fill = "lightblue",
    color = "white"
  ) + 
  labs(
    title = "Distribution of the Composite Cardiovascular Risk Score",
    subtitle = paste(
      "Mean =", round(mean(df_all_risk_scores$mean_risk, na.rm = TRUE), 2), "%",
      "| N =", sum(!is.na(df_all_risk_scores$mean_risk))
    ),
    x = "10-year cardiovascular risk (%)",
    y = "Number of participants"
  ) +
    geom_density(
      aes(y = ..density.. * nrow(df_all_risk_scores) *
            diff(range(df_all_risk_scores$mean_risk, na.rm = TRUE)) / 20),
      color = "darkblue",
      size = 1.2
  ) +
  geom_vline(
    aes(xintercept = mean(mean_risk, na.rm = TRUE)),
    color = "red", linetype = "dashed", size = 1
  )
plot_composite

###########################################################################
## Plotting the CVD score results in combined figures
###########################################################################

figures_path <- "/Users/luisadelius/Documents/Code/project_one/Figures"

### plotting the distribution plots in a combined figure
combined_plot <- (plot_ASCVD | plot_framigham | plot_qrisk3 | plot_SCORE2) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.margin = margin(5, 10, 5, 10)
  )
combined_plot <- combined_plot & ylim(0, 30)
combined_plot

ggsave(
  filename = file.path(figures_path, "CVD_risk_distributions_wo_composite.png"),
  plot = combined_plot,
  width = 16,
  height = 4.5,
  dpi = 300
)

###########################################################################
########### Violine plot + Kendall's Tau correlation matrix ##############
###########################################################################
## plotting the risk scores in a combined figure
df_long <- df_all_risk_scores %>%
  pivot_longer(
    cols = c(SCORE2_score, ascvd_10y, frs_10y, QRISK3_risk),
    names_to = "score_type",
    values_to = "risk"
  ) %>%
  select(PatientID, score_type, risk)


## create Kendall's Tau correlation matrix
# subset relevant columns
risk_mat <- df_all_risk_scores %>%
  select(SCORE2_score, ascvd_10y, frs_10y, QRISK3_risk)

# Kendall correlation matrix
cor_kendall <- cor(risk_mat, method = "kendall", use = "pairwise.complete.obs") 
#with pairwise complete each correlation uses the maximum available data. If i do "complete" instead, only participatns with no NA will be used -> loosing some data.

cor_kendall # call the matrix


plot_A <- ggplot(df_long, aes(x = score_type, y = risk)) +
  geom_violin(aes(fill = score_type), alpha = 0.4, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  geom_line(aes(group = PatientID), colour = "black", alpha = 0.25) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_fill_manual(
    values = c(
      "ascvd_10y"    = "#F8766D",
      "frs_10y"      = "#8AAD40",
      "QRISK3_risk"  = "#4DB8B0",
      "SCORE2_score" = "#C77CFF"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "SCORE2_score" = "SCORE2",
      "QRISK3_risk"  = "QRISK3",
      "ascvd_10y"    = "ASCVD",
      "frs_10y"      = "Framingham"
    )
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "10-year Risk (%)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

plot_B <- as.data.frame(cor_kendall) %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "tau") %>%
  ggplot(aes(var1, var2, fill = tau)) +
  geom_tile() +
  geom_text(aes(label = round(tau, 2)), size = 5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1)) +
  scale_x_discrete(
    labels = c(
      "SCORE2_score" = "SCORE2",
      "QRISK3_risk"  = "QRISK3",
      "ascvd_10y"    = "ASCVD",
      "frs_10y"      = "Framingham"
    )
  ) +
  scale_y_discrete(
    labels = c(
      "SCORE2_score" = "SCORE2",
      "QRISK3_risk"  = "QRISK3",
      "ascvd_10y"    = "ASCVD",
      "frs_10y"      = "Framingham"
    )
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  labs(x = "", y = "")

# Combine: A is ~65% width, B is ~35%
combined_fig <- (plot_A | plot_B) +
  plot_layout(widths = c(3, 1.5)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 18, face = "bold")
  )

combined_fig

ggsave(
  filename = file.path(figures_path, "CVD_risk_violin_correlation.png"),
  plot = combined_fig,
  width = 14,
  height = 6,
  dpi = 300
)

###########################################################################
# Violin plot only for poster presentation
###########################################################################
plot_A_poster <- ggplot(df_long, aes(x = score_type, y = risk)) +
  geom_violin(aes(fill = score_type), alpha = 0.4, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  geom_line(aes(group = PatientID), colour = "black", alpha = 0.25) +
  geom_point(size = 2, alpha = 0.7) +
  scale_fill_manual(
    values = c(
      "ascvd_10y"    = "#F8766D",
      "frs_10y"      = "#8AAD40",
      "QRISK3_risk"  = "#4DB8B0",
      "SCORE2_score" = "#C77CFF"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "SCORE2_score" = "SCORE2",
      "QRISK3_risk"  = "QRISK3",
      "ascvd_10y"    = "ASCVD",
      "frs_10y"      = "Framingham"
    )
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "10-year Risk (%)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(angle = 25, hjust = 1, size = 20),
    axis.text.y  = element_text(size = 20),
    axis.title.y = element_text(size = 24)
  )

ggsave(
  filename = file.path(figures_path, "CVD_risk_violin_poster.tiff"),
  plot     = plot_A_poster,
  width    = 8,
  height   = 7,
  dpi      = 600,
  compression = "lzw"
)

###########################################################################
#. Supplementary Figure with CVD score input and outcomes
###########################################################################
CVD_input_plus_outcome <- df_all_risk_scores %>%
  select(-SCORE2_strat)%>%
  full_join(risk_factor_input, by = "PatientID") %>%
  select(-Sex, -race_ascvd, -SmokingStatusQRISK3, -Diabetes, 
         -blood_pressure_treatment, -mean_LDL_mg_dl, -systolic, -Age) %>%
  full_join(QRISK3_calculation_input, by = "PatientID") 

# # Prepare the data
# CVD_table <- CVD_input_plus_outcome %>%
#   filter(!is.na(mean_risk)) %>%  # Exclude participants with NA in mean_risk
#   select(
#     PatientID, Age, Sex, systolic, blood_pressure_treatment, 
#     SmokingStatusQRISK3, Weight_kg, Height_cm, diabetes2, 
#     Severe_mental_illness, EthnicityCodeQRISK3, townsend, 
#     Risk.region, ratio_chol_hdl, mean_Total_Cholesterol_mg_dl, 
#     mean_HDL_mg_dl, SCORE2_score, ascvd_10y, frs_10y, 
#     QRISK3_risk, mean_risk
#   ) %>%
#   arrange(PatientID) %>%  # Sort by PatientID
#   mutate(across(where(is.numeric), ~round(., 2))) %>%  # Round to 2 decimal places
#   rename(
#     Composite = mean_risk,
#     SCORE2 = SCORE2_score,
#     ASCVD = ascvd_10y,
#     Framingham = frs_10y,
#     QRISK3 = QRISK3_risk
#   )

# # Create flextable with smaller size and minimal spacing
# ft_cvd <- flextable(CVD_table) %>%
#   theme_booktabs() %>%
#   colformat_double(digits = 2, na_str = "NA") %>%
#   height_all(height = 0.15) %>%  # Minimal row height
#   hrule(rule = "exact") %>%  # Enforce exact row heights
#   set_table_properties(width = 1, layout = "autofit") %>%
#   autofit()
# 
# ft_cvd
# 
# # Save the table in landscape orientation
# save_as_docx(
#   ft_cvd, 
#   path = file.path(figures_path, "CVD_input_output_table.docx"),
#   pr_section = prop_section(
#     page_size = page_size(orient = "landscape"),
#     page_margins = page_mar(bottom = 0.5, top = 0.5, right = 0.5, left = 0.5)
#   )
# )




###########################################################################
# 6. QRISK3 calculation with statin exclusion
###########################################################################
# Option 1: Fix the Sample_ID in df_statins_supplements before joining
qrisk_input <- QRISK3_calculation_input %>%
  rename(Sample_ID = PatientID) %>%
  left_join(
    df_statins_supplements %>% 
      mutate(Sample_ID = str_replace_all(Sample_ID, "-", "_")) %>%
      select(Sample_ID, Statins),
    by = "Sample_ID"
  )

qrisk_input_wo_statins <- qrisk_input %>%
  filter(Statins == "No" | is.na(Statins))  

# Create dataframe to calculate QRISK3 (every condition gets their values assigned)
# Using statin-free participants only
QRISK3_wo_statins <- data.frame(
  gender = if_else(qrisk_input_wo_statins$Sex == "Female", 1,
                   if_else(qrisk_input_wo_statins$Sex == "Male", 0, NA_integer_)),
  age = qrisk_input_wo_statins$Age,
  b_AF = 0, b_atypicalantipsy = 0, b_corticosteroids = 0,
  b_impotence2 = 0, b_migraine = 0, b_ra = 0, b_renal = 0,
  b_semi = if_else(qrisk_input_wo_statins$Severe_mental_illness == "Yes", 1, 0),
  b_sle = 0,
  b_treatedhyp = if_else(qrisk_input_wo_statins$blood_pressure_treatment == "Yes", 1, 0),
  b_type1 = 0,
  b_type2 = if_else(qrisk_input_wo_statins$diabetes2 == "Yes", 1, 0),
  weight = qrisk_input_wo_statins$Weight_kg,
  height = qrisk_input_wo_statins$Height_cm,
  ethrisk = qrisk_input_wo_statins$EthnicityCodeQRISK3,
  fh_cvd = 0,
  rati = qrisk_input_wo_statins$ratio_chol_hdl,
  sbp  = qrisk_input_wo_statins$systolic,
  sbps5 = 10,
  smoke_cat = qrisk_input_wo_statins$SmokingStatusQRISK3,
  town = ifelse(is.na(qrisk_input_wo_statins$townsend), 0, qrisk_input_wo_statins$townsend),
  Sample_ID = qrisk_input_wo_statins$Sample_ID) %>%  # Note: using Sample_ID not PatientID
  mutate(ID = rownames(.)) %>%
  filter(!if_any(everything(), is.na)) %>% # cannot have NA values
  filter(age >= 25 & age <= 84) # must be between 25 and 84

# Calculate the QRISK3
QRISK3_scores_wo_statins <- QRISK3_2017(
  data = QRISK3_wo_statins, 
  patid = "ID", 
  gender = "gender", 
  age = "age", 
  atrial_fibrillation = "b_AF", 
  atypical_antipsy = "b_atypicalantipsy", 
  regular_steroid_tablets = "b_corticosteroids", 
  erectile_disfunction = "b_impotence2", 
  migraine = "b_migraine", 
  rheumatoid_arthritis = "b_ra", 
  chronic_kidney_disease = "b_renal", 
  severe_mental_illness = "b_semi", 
  systemic_lupus_erythematosis = "b_sle", 
  blood_pressure_treatment = "b_treatedhyp", 
  diabetes1 = "b_type1", 
  diabetes2 = "b_type2", 
  weight = "weight", 
  height = "height", 
  ethiniciy = "ethrisk", 
  heart_attack_relative = "fh_cvd", 
  cholesterol_HDL_ratio = "rati", 
  systolic_blood_pressure = "sbp", 
  std_systolic_blood_pressure = "sbps5", 
  smoke = "smoke_cat", 
  townsend = "town"
)


QRISK3_sample_ID_wo_statins <- QRISK3_wo_statins %>%
  inner_join(QRISK3_scores_wo_statins, by = "ID") %>%
  select(Sample_ID, QRISK3_2017) %>%
  rename(QRISK3_wo_statins = QRISK3_2017)

# Save your dataframe with patient IDs and QRISK3 scores
saveRDS(QRISK3_sample_ID_wo_statins, "QRISK3_sample_ID_wo_statins.RDS")

ggplot(QRISK3_sample_ID_wo_statins, aes(x = QRISK3_wo_statins)) +
  geom_histogram(
    bins = 20,
    fill = "lightblue",
    color = "white"
  ) + 
  labs(
    title = "Distribution of QRISK3 Scores w/o Statins",
    subtitle = paste("Mean =", round(mean(QRISK3_sample_ID_wo_statins$QRISK3_wo_statins, na.rm = TRUE), 2), "%",
                     "| N =", nrow(QRISK3_sample_ID_wo_statins)),
    x = "10-year cardiovascular risk (%)",
    y = "Number of participants"
  ) +
  geom_density(
    aes(y = ..density.. * nrow(QRISK3_sample_ID_wo_statins) * 
          diff(range(QRISK3_sample_ID_wo_statins$QRISK3_wo_statins)) / 20),
    color = "darkblue",
    size = 1.2
  ) +
  geom_vline(
    aes(xintercept = mean(QRISK3_wo_statins, na.rm = TRUE)),
    color = "red", linetype = "dashed", size = 1
  )
