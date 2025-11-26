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
    HDL_mg_dl               = parse_number(HDL_mg_dl)
  ) %>%
  group_by(PatientID) %>%
  summarise( # calculates the mean each and the ratio of the means
    mean_Total_Cholesterol_mg_dl = mean(Total_Cholesterol_mg_dl, na.rm = TRUE),
    mean_HDL_mg_dl = mean(HDL_mg_dl, na.rm = TRUE),
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
      recruitment_site %in% c("CBG", "ICL") ~ "low",
      recruitment_site %in% c("CORK", "AUTH", "UVEG") ~ "moderate",
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
ggplot(Framingham, aes(x = frs_10y)) +
  geom_histogram( # create the normal histogram
    bins = 20,
    fill = "lightblue",
    color = "white"
  ) + 
  labs( # create the titles 
    title = "Distribution of Framingham Scores",
    subtitle = paste("Mean =", round(mean(Framingham$frs_10y, na.rm = TRUE), 2), "%",
    "| N =", nrow(Framingham)),
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
  )

# Statistics for the Framingham distribution
shapiro.test(Framingham$frs_10y)  # normality test
summary(Framingham$frs_10y)

##---------------------------------------------------------------------##
## 2. ASCVD Score
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
  diabetes = risk_factor_input$Diabetes
) %>%
  filter(!if_any(everything(), is.na)) %>%    # cannot have NA values
  filter(age >= 40 & age <= 79) %>%           # age range according to website
  mutate(
    ascvd_10y = ascvd_10y_accaha(
      race, gender, age, totchol, hdl, sbp, bp_med, smoker, diabetes
    )
  )

# just checking the people that got excluded because they have values outside specific ranges
ASCVD %>%
  filter(is.na(ascvd_10y)) %>%
  select(age, totchol, hdl, sbp)

### checking distribution and mean of ASCVD score
# Histogram + density for ASCVD scores
ggplot(ASCVD, aes(x = ascvd_10y)) +
  geom_histogram( # create the normal histogram
    bins = 20,
    fill = "lightblue",
    color = "white"
  ) + 
  labs( # create the titles 
    title = "Distribution of ASCVD Scores",
    subtitle = paste("Mean =", round(mean(ASCVD$ascvd_10y, na.rm = TRUE), 2), "%",
    "| N =", nrow(ASCVD)),
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
  )

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
ggplot(SCORE2, aes(x = SCORE2_score)) +
  geom_histogram(
    bins = 20,
    fill = "lightblue",
    color = "white"
  ) + 
  labs(
    title = "Distribution of SCORE2 10-year Risk",
    subtitle = paste(
      "Mean =", round(mean(SCORE2$SCORE2_score, na.rm = TRUE), 2), "%",
      "| N =", nrow(SCORE2)
    ),
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
  )


# Statistics for the SCORE2 distribution
shapiro.test(SCORE2$SCORE2_score)   # normality test
summary(SCORE2$SCORE2_score)        # summary stats
##---------------------------------------------------------------------##

# 4. Combine the results
## create a joint data frame
df_all_risk_scores <- SCORE2 %>%
  full_join(ASCVD %>% select(PatientID, ascvd_10y), by = "PatientID") %>%
  full_join(Framingham %>% select(PatientID, frs_10y),by = "PatientID") %>%
  full_join(QRISK3_sample_ID %>% rename(PatientID = Sample_ID), by = "PatientID") %>%
  mutate(mean_risk = rowMeans( #calculate an average value for each patient --> composite score
    select(., SCORE2_score, ascvd_10y, frs_10y, QRISK3_2017),
    na.rm = TRUE
  ))

## histogram of the composite risk score
ggplot(df_all_risk_scores, aes(x = mean_risk)) +
  geom_histogram(
    bins = 20,
    fill = "lightblue",
    color = "white"
  ) + 
  labs(
    title = "Distribution of Composite Cardiovascular Risk Score",
    subtitle = paste(
      "Mean =", round(mean(df_all_risk_scores$mean_risk, na.rm = TRUE), 2), "%",
      "| N =", sum(!is.na(df_all_risk_scores$mean_risk))
    ),
    x = "Composite 10-year cardiovascular risk (%)",
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


## plotting the risk scores in a combined figure
df_long <- df_all_risk_scores %>%
  pivot_longer(
    cols = c(SCORE2_score, ascvd_10y, frs_10y, QRISK3_2017),
    names_to = "score_type",
    values_to = "risk"
  ) %>%
  select(PatientID, score_type, risk)


ggplot(df_long, aes(x = score_type, y = risk)) +
  # 1. violins (overall distribution)
  geom_violin(aes(fill = score_type), alpha = 0.4, colour = NA) +
  
  # 2. boxplot (median + IQR)
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  
  # 3. patient lines (connects a patient's scores)
  geom_line(aes(group = PatientID), colour = "black", alpha = 0.25) +
  
  # 4. dots (individual values)
  geom_point(size = 1.5, alpha = 0.7) +
  
  # change legend
  scale_fill_discrete(
    labels = c(
      "SCORE2_score" = "SCORE2",
      "QRISK3_2017"  = "QRISK3",
      "ascvd_10y"    = "ASCVD",
      "frs_10y"      = "Framingham"
    )
  ) +
  
  labs(
    title = "Cardiovascular Risk Scores per Patient",
    x = "Score Type",
    y = "10-year Risk (%)"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )



