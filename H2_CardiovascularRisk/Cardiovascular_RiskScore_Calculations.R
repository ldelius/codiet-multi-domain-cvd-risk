### Cardiovascular Risk Score Calculations
### Author: Luisa Delius

# this file is to calculate several cardiovascular risk scores
# 1. Framingham score
# 2. ASCVD
# 3. SCORE2

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

## combining the cleaned data
risk_factor_input <- maindata_keep %>%
  full_join(lipids_keep, by = "PatientID") %>%
  full_join(blood_pressure_keep, by = "PatientID") %>%
  full_join(patient_age, by = "PatientID")

###################### Risk score calulations ######################
## 1. Framingham Score
Framingham <- data.frame(
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

unique(tolower(risk_factor_input$Sex))

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
    subtitle = paste("Mean =", round(mean(Framingham$frs_10y, na.rm = TRUE), 2), "%"),
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
