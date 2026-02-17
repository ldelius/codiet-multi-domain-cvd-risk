### QRISK3 score Calculation
### Author: Luisa Delius

# Load packages
library(tidyverse)
library(broom)
library(patchwork)
library(QRISK3)
library(readxl)
library(QDiabetes)

# Set working directory
wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

# Upload the data/excel sheet
excel_sheets("QRISK3_data.xlsx")
df_maindata <- read_excel("QRISK3_data.xlsx", sheet = "Sheet1")
df_Cholesterol_HDL  <- read_excel("QRISK3_data.xlsx", sheet = "Cholesterol_HDL")
df_Blood_Pressure     <- read_excel("QRISK3_data.xlsx", sheet = "Blood_Pressure")
df_age_risk_factor_sheet <- read_excel("QRISK3_data.xlsx", sheet = "Age_risk_factor_sheet")
df_UK_postcodes <- read_excel(path = "cardiovascular_riskfactor_calculations.xlsx", sheet = "Risk_region")


names(df_Cholesterol_HDL) # gives you column titles 
glimpse(df_Cholesterol_HDL)
summary(df_Cholesterol_HDL[, c("Total_Cholesterol_mg_dl","HDL_mg_dl")])


# Step 1.1: Calculate ratio cholesterol/HDL for each person only once (each patient has 2-3 meassuremnts in the list, so calculating the mean of them)
lipids_by_patient <- df_Cholesterol_HDL %>%
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
    ratio_chol_hdl = mean(Total_Cholesterol_mg_dl, na.rm = TRUE) /
      mean(HDL_mg_dl, na.rm = TRUE),
    .groups = "drop" # removes all grouping (the data got grouped by (group_by() to calculate each patients mean))
  )

head(lipids_by_patient)
str(lipids_by_patient)

# Step 1.2.: create a df with one age for each person (currently every patient has an age per visit)
patient_age <- df_age_risk_factor_sheet %>%
  mutate(PatientID = str_remove(SampleID, "_V[1-3]+$")) %>%
  group_by(PatientID) %>%
  summarise(Age = mean(Age, na.rm = TRUE)) %>%
  ungroup() # not needed, but i will keep it after every group to not forget about it.
  

# Step 1.3 Calculate townsend score for ICL participants
put_TDS_or_NA <- function(postcode) {
  tryCatch(getTDS(postcode), error = function(e) NA_real_)
}

df_townsend <- df_UK_postcodes %>%
  rename(PatientID = volunteer_id) %>%   # unify ID naming
  mutate(
    PatientID = str_replace_all(PatientID, "-", "_"),     
    postcodes = str_trim(as.character(postcodes)),        # clean postcode
    townsend = sapply(postcodes, put_TDS_or_NA)           # Townsend or NA
  )


# Step 2: Prepare everything for joining to one dataframe with all the information
## create same key for joining & keep only the coulmns wanted for joining
maindata_keep <- df_maindata %>%
  rename(PatientID = volunteer_id) %>%   # make the ID name consistent
  mutate(
    PatientID = str_replace_all(PatientID, "-", "_")  # CD-003 → CD_003
  ) %>%
  select(
    PatientID, Sex, EthnicityCodeQRISK3, SmokingStatusQRISK3, diabetes2, Weight_kg,
    Height_cm, blood_pressure_treatment, Severe_mental_illness 
  )


blood_pressure_keep <- df_Blood_Pressure %>%
  rename(PatientID = patient) %>%   # make the ID name consistent
  mutate(
    PatientID = str_replace_all(PatientID, "-", "_")  # CD-003 → CD_003
  ) %>%
  select(
    PatientID, systolic
  )

lipids_by_patient_keep <- lipids_by_patient %>%
  select(
    PatientID, ratio_chol_hdl
    )


# Step 3: Combining the columns needed for the QRISK3 Calculation
qrisk_input <- maindata_keep %>%
  full_join(lipids_by_patient_keep, by = "PatientID") %>%
  full_join(blood_pressure_keep, by = "PatientID") %>%
  full_join(patient_age, by = "PatientID") %>%
  full_join(df_townsend %>% select(PatientID, townsend), by = "PatientID")

saveRDS(qrisk_input, file.path(wkdir, "processed_data", "QRISK3_calculation_input.rds"))

# Create dataframe to calculate QRISK3 (every condition gets their values assigned)
QRISK3 <- data.frame(
                    gender = if_else(qrisk_input$Sex == "Female", 1,
                      if_else(qrisk_input$Sex == "Male", 0, NA_integer_)),
                     age = qrisk_input$Age,
                     b_AF = 0, b_atypicalantipsy = 0, b_corticosteroids = 0,
                     b_impotence2 = 0, b_migraine = 0, b_ra = 0, b_renal = 0,
                     b_semi = if_else(qrisk_input$Severe_mental_illness == "Yes", 1, 0),
                     b_sle = 0,
                     b_treatedhyp = if_else(qrisk_input$blood_pressure_treatment == "Yes", 1, 0),
                     b_type1 = 0,
                     b_type2 = if_else(qrisk_input$diabetes2 == "Yes", 1, 0),
                     weight = qrisk_input$Weight_kg,
                     height = qrisk_input$Height_cm,
                     ethrisk = qrisk_input$EthnicityCodeQRISK3,
                     fh_cvd = 0,
                     rati = qrisk_input$ratio_chol_hdl,
                     sbp  = qrisk_input$systolic,
                     sbps5 = 10,
                     smoke_cat = qrisk_input$SmokingStatusQRISK3,
                    town = ifelse(is.na(qrisk_input$townsend), 0, qrisk_input$townsend),
                    Sample_ID = qrisk_input$PatientID) %>%
  mutate(ID = rownames(.)) %>% ## in der QRISK df starten die IDs bei 2, weil ich in der nächsten Zeile die erste Zeile wegen einem NA raugeschmissen habe.
  filter(!if_any(everything(), is.na)) %>% # cannot have NA values
  filter(age >= 25 & age <= 84) # must be between 25 and 84

# Calculate the QRISK3
QRISK3_test <- QRISK3_2017(data=QRISK3, patid="ID", gender="gender", age="age", 
                           atrial_fibrillation="b_AF", atypical_antipsy="b_atypicalantipsy", 
                           regular_steroid_tablets="b_corticosteroids", 
                           erectile_disfunction="b_impotence2", migraine="b_migraine", 
                           rheumatoid_arthritis="b_ra", chronic_kidney_disease="b_renal", 
                           severe_mental_illness="b_semi", systemic_lupus_erythematosis="b_sle", 
                           blood_pressure_treatment="b_treatedhyp", diabetes1="b_type1", 
                           diabetes2="b_type2", weight="weight", height="height", 
                           ethiniciy="ethrisk", heart_attack_relative="fh_cvd", 
                           cholesterol_HDL_ratio="rati", systolic_blood_pressure="sbp", 
                           std_systolic_blood_pressure="sbps5", smoke="smoke_cat", townsend="town")

# link the QRISK back to the partient!
QRISK3_sample_ID <- QRISK3 %>%
  inner_join(QRISK3_test, by = "ID") %>%
  select(Sample_ID, QRISK3_2017) %>%
  rename(QRISK3_risk = QRISK3_2017)
## prints out the The calculated 10-year cardiovascular risk (%) (to compare, mine is 0.1%)

# Save your dataframe with patient IDs and QRISK3 scores
saveRDS(QRISK3_sample_ID, "QRISK3_sample_ID.RDS")


# Histogram of the QRISK values
ggplot(QRISK3_sample_ID, aes(x = QRISK3_risk)) +
  geom_histogram( #create the normal histogram
    bins = 20,
    fill = "lightblue",
    color = "white"
  ) + 
  labs( # create the titles 
    title = "Distribution of QRISK3 Scores",
    subtitle = paste("Mean =", round(mean(QRISK3_sample_ID$QRISK3_risk, na.rm = TRUE), 2), "%",
    "| N =", nrow(QRISK3_sample_ID)),
    x = "10-year cardiovascular risk (%)",
    y = "Number of participants"   )+
  geom_density( # distribution
    aes(y = ..density.. * nrow(QRISK3_sample_ID) * diff(range(QRISK3_sample_ID$QRISK3_risk)) / 20),
    color = "darkblue",
    size = 1.2
  ) +
  geom_vline(aes(xintercept = mean(QRISK3_risk, na.rm = TRUE)), # mean QRISK value
             color = "red", linetype = "dashed", size = 1
  )

# Statistics to the QRISK3 distribution
shapiro.test(QRISK3_sample_ID$QRISK3_risk) # does the data follow a normal distribution? H0 = data floows normal distribution --> <0.05 means reject H0 --> data is not normal.
summary(QRISK3_sample_ID$QRISK3_risk)


###############################################################################
###############################################################################

### --------- Preparing the risk_factor dataset -----------###
# Load the data
df_risk_factors <- read_csv("risk_factors.csv")

# Step 1: Prepare the data, so that I can summarize the visits of each patient in Step 2.
## I have several visits per patient, but QRISK3 only calculated per patient once. Therefore, I have to change the risk_factor so that i have one value per patient as well.
df_patient_risk <-  df_risk_factors %>%
  mutate(
    PatientID = str_remove(SampleID, "_V[1-3]+$"), # new column PatientID without visit
  ## change categorials to factors to be able to calculate with them
    Country = factor(Country),
    Gender = factor(Gender),
    Age.Risk = factor(Age.Risk, levels = c("Normal", "I", "II", "III")), ## how to treat NA??
    stress_resilience_status = factor(stress_resilience_status,
                                      levels = c("poor", "Bad", "Normal", "Good", "Excellent"),
                                      ordered = TRUE),
    stress_index_status = factor(stress_index_status,
                                 levels = c("poor", "Bad", "Normal", "Good", "Excellent"),
                                 ordered = TRUE)
  )
  
# This shows patients where the same categorical column had different values across visits
inconsistent_categoricals <- df_patient_risk %>%
  group_by(PatientID) %>%
  summarise(across(where(is.factor), n_distinct), .groups = "drop") %>% # how many different values appear for each patient?
  filter(if_any(everything(), ~ .x > 1))
# I changed the age and gender accordingly after running those lines of code

# sort(unique(df_patient_risk$stress_index_status)) # i had used that line of code to find the level of each column out


# Step 2: Create a df with the summarized per patient data, ready to use for further analysis.
df_patient_risk_for_lm <- df_patient_risk %>%
  group_by(PatientID) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), # numeric columns → mean across visits
    
    # categorical columns: when possible: mode; otherwise the middle one
    across(where(is.factor), function(.x) { # for each factor apply this function for all the patient values
            x <- na.omit(.x) # removes missing values
            
            if (length(x) == 0) {
              return(factor(NA, levels = levels(.x), ordered = is.ordered(.x)))
            } # important to keep factor and not go back to chr
             
            freq_tab <- table(x)
            winners  <- names(freq_tab)[freq_tab == max(freq_tab)]
            factor(
              if (length(winners) == 1) {
                winners
              } else {
                levs    <- levels(.x)              # preserve original level order
                present <- levs[levs %in% winners]
                present[ceiling(length(present) / 2)]
              },
              levels  = levels(.x),
              ordered = is.ordered(.x)
            )
    }),
    
    .groups = "drop"
  )

# Step 3: "Ensure the Urine tPred scores are going in the right direction.
# Phenylacetylglutamine score (animal protein) and O-Acetylcarnitine score (red meat) should be negative;
# everything else should be positive" --> therefore i am flipping those two around.
df_patient_risk_for_lm <- df_patient_risk_for_lm %>%
mutate(
  Phenylacetylglutamine = -Phenylacetylglutamine,
  `O-Acetylcarnitine`   = -`O-Acetylcarnitine`
)



# add the QRISK3_risk results to the df with risk_factors
df_patient_risk_for_lm <- df_patient_risk_for_lm %>%
  rename(Sample_ID = PatientID) %>%
  left_join(QRISK3_sample_ID, by = "Sample_ID")

df_patient_risk_for_lm %>%
  summarise(
    mean_age = mean(Age, na.rm = TRUE),
    sd_age   = sd(Age, na.rm = TRUE)
  )

# Standardise the numeric data (zscale) for the lm
## I z-scaled all the numeric columns, however not the QRISK3 score. I did that, so that the predictors are comparable. Didn’t do the outcome because I want the original, interpretable units/results.
df_scaled_patient_risk_for_lm <- df_patient_risk_for_lm %>% 
  mutate(across(where(is.double) & !matches("QRISK3_risk"), ~ scale(.x), .names = "z_{col}")) %>%
  rename_with(function(x) gsub("[- ]", "_", x))

saveRDS(df_scaled_patient_risk_for_lm,
  file = file.path(wkdir, "processed_data", "df_risk_factor_predictors.rds"))


### --------- Running and Plotting normal Linear Regression Model -----------###
# Run several LMs with QRISK and z_scores
predictors <- df_scaled_patient_risk_for_lm %>%
  select(
    starts_with("z_"),
    where(is.factor)
    ) %>% # selects all columns whose names begin with z_ & factors. Those variables shall be modeled against my selected group.
  names() #creates a simple list of column names
predictors

# one lm per predictor, dropping NAs only for that predictor + outcome
lm_output <- map_dfr(predictors, function(var) { # runs the function over each variable and collects them in a df.
  df_pair <- df_scaled_patient_risk_for_lm %>%
    select(QRISK3_risk, all_of(var)) %>% # selects always one variable and the QRISK score
    drop_na() # and drops the NAs specific for those columns. 
  
  model <- lm(as.formula(paste("QRISK3_risk ~", var)), data = df_pair)
  sw <- shapiro.test(residuals(model))  # this is to test the disitribution of the residuals of the lm in order to see, whether I used the fitting linear model.
  
  tidy(model, conf.int = TRUE) %>% # returns a clean tibble
    filter(term != "(Intercept)") %>% # keep predictor terms only
    mutate(predictor = var,
           n = nrow(df_pair),   # add name and row (context)
           shapiro_W = unname(sw$statistic),   # <- added distribtion test result to the table
           shapiro_p = sw$p.value              
           )
}) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>% # multiple testing adjustment 
  arrange(p.adj) # arranges by adjusted pvalue

lm_output

summary(lm_output$shapiro_p)

###### --------------- Running a Generalized Linear Model ---------------- #####
# GLM with Gamma distribution and log link
# For GLMs don’t test normality, as the residuals are not supposed to be normal

# 1. Create predictors 
# Numeric predictors: all z_ columns
num_predictors <- df_scaled_patient_risk_for_lm %>%
  select(starts_with("z_")) %>%
  names()

# Factor predictors
factor_predictors <- df_scaled_patient_risk_for_lm %>%
  select(where(is.factor)) %>%
  names()
factor_predictors

# Combined for plotting
plot_predictors <- c(num_predictors, factor_predictors)

# 2. GLM for numeric predictors
glm_output_num <- map_dfr(num_predictors, function(var) {
  
  df_pair <- df_scaled_patient_risk_for_lm %>%
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
      level     = NA_character_,
      term      = term,
      n         = nrow(df_pair)
    )
})


# 3. GLM for factor predictors (dummy style, one GLM per level)
glm_output_fac <- map_dfr(factor_predictors, function(var) {
  
  df_pair <- df_scaled_patient_risk_for_lm %>%
    select(QRISK3_risk, all_of(var)) %>%
    drop_na()
  
  levs <- levels(df_pair[[var]])
  
  map_dfr(levs, function(lev) {
    
    df_tmp <- df_pair %>%
      mutate(dummy = if_else(.data[[var]] == lev, 1, 0))
    
    model <- glm(
      QRISK3_risk ~ dummy,
      data   = df_tmp,
      family = Gamma(link = "log")
    )
    
    tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term == "dummy") %>%
      mutate(
        predictor = var,
        level     = lev,
        term      = paste0(var, " = ", lev),
        n         = nrow(df_tmp)
      )
  })
})

# 4. Combine, adjust p-values, significance, ordering
glm_output <- bind_rows(glm_output_num, glm_output_fac) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    significant = if_else(conf.low > 1 | conf.high < 1,
                          "Significant", "Not significant")
  ) %>%
  arrange(p.adj)


###### -- Running GLMM (Generalised Linear Mixed Model) with Random Effects (Gender & Country)-- #####
predictors <- df_scaled_patient_risk_for_lm %>%
  select(
    starts_with("z_"),
    where(is.factor)
  ) %>% # selects all columns whose names begin with z_ & factors. Those variables shall be modeled against my selected group.
  select(-Gender, -Country) %>%     # do not want to model Gender as fixed effect here
  names() #creates a simple list of column names
predictors

glmm_gender_country <- map_dfr(predictors, function(var) {
  
  df_pair <- df_scaled_patient_risk_for_lm %>%
    select(QRISK3_risk, Gender, Country, all_of(var)) %>%
    drop_na()
  
  if (nrow(df_pair) == 0 || n_distinct(df_pair$Country) < 2) {
    message("Skipping ", var, ": only one Country level in subset after drop_na().")
    return(tibble())
  }
  
  form <- as.formula(paste("QRISK3_risk ~", var, "+ (1 | Gender) + (1 | Country)")) # fixed effect = predictor, random effect = Gender (&Country)
  
  model <- lme4::glmer(
    form,
    data   = df_pair,
    family = Gamma(link = "log")
  )
  
  coef_df <- as.data.frame(summary(model)$coefficients)
  coef_df$term <- rownames(coef_df)
  
  coef_df %>%
    filter(term != "(Intercept)") %>%
    mutate(
      conf.low  = Estimate - 1.96 * `Std. Error`,
      conf.high = Estimate + 1.96 * `Std. Error`,
      predictor = var,
      n         = nrow(df_pair)
    )
  
}) %>%
  mutate(p.adj = p.adjust(`Pr(>|z|)`, method = "BH")) %>%
  arrange(p.adj)


###### ---------------- Plotting my LM models ---------------- ######
# one plot per predictor (handles NA per predictor)
make_plot <- function(var) {
  d <- df_scaled_patient_risk_for_lm %>%
    select(QRISK3_risk, all_of(var))%>%
    drop_na()
  
  gg <- ggplot(d, aes(x = .data[[var]], y = QRISK3_risk)) +
    labs(x = var, y = "QRISK3_risk")
  
  if (is.factor(d[[var]])) {
    gg + geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.5, size = 1.4)
  } else {
    gg + geom_point(size = 1.6, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE)
  }
}

# build all plots
plots <- map(predictors, make_plot) # for each element in predictors, run make_plot() and collect the ggplots in a list "plots".

# show 6 per page and print all pages:
lapply(split(plots, ceiling(seq_along(plots) / 6)),
       \(p) print(wrap_plots(p, ncol = 3)))



## create a confidence plot with lm_output (lm_output already created above and arranged by p.adj)
ord <- rev(lm_output$term)  # keep plotted order with the most significant one at top 

ggplot(lm_output,
       aes(x = estimate, y = factor(term, levels = ord))) +
  geom_vline(xintercept = 0, linetype = "dashed") + # this adds the "no effect" reference line at 0.
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point() +
  labs(
    x = "Effect on QRISK3_risk",
    y = "Predictor",
    title = "LM with 95% CIs (BH-adjusted ordering)"
  ) +
  theme(axis.text.y = element_text(size = 6)) # I added this to make the yaxis text smaller, so that it doesnt overlap.


###### ---------------- Plotting my GLM models (Gamma, log link) ---------------- ######
# order: most according to estimate
ord <- glm_output %>%
  arrange(estimate) %>%
  pull(term)

# creating single-predictor plots
make_plot_glm <- function(var) {
  d <- df_scaled_patient_risk_for_lm %>%
    select(QRISK3_risk, all_of(var)) %>%
    drop_na()
  
  gg <- ggplot(d, aes(x = .data[[var]], y = QRISK3_risk)) +
    labs(x = var, y = "QRISK3_risk")
  
  if (is.factor(d[[var]])) {
    gg + geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.5, size = 1.4)
  } else {
    gg + geom_point(size = 1.0, alpha = 0.7) +
      geom_smooth(method = "glm",
                  method.args = list(family = Gamma(link = "log")),
                  se = TRUE)
  }
}

# build all plots
plots_glm <- purrr::map(plot_predictors, make_plot_glm)

# show 12 per page and print all pages:
lapply(
  split(plots_glm, ceiling(seq_along(plots_glm) / 12)),
  \(p) print(wrap_plots(p, ncol = 4, nrow = 3)))

# forest plot
ggplot(glm_output %>%
         filter(!(predictor == "Age.Risk" & level == "III")),   # drop only Age.Risk III
       aes(x = estimate, y = factor(term, levels = ord), colour = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed") +     # no-effect line is 1 (ratio scale)
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point() +
  scale_colour_manual(values = c("Significant" = "firebrick3",
                                 "Not significant" = "black")) +
  labs(
    x = "Ratio of expected QRISK3 (exp(beta))",
    y = "Predictor",
    title = "GLM with 95% CIs (Estimate ordering)"
  ) +
  theme(axis.text.y = element_text(size = 6))

##### plot only significant and country
ggplot(
  glm_output %>%
    filter(significant == "Significant" | str_detect(term, "Country")),
  aes(x = estimate, y = factor(term, levels = ord), colour = significant)
) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point() +
  scale_colour_manual(values = c("Significant" = "firebrick3")) +
  labs(
    x = "Ratio of expected QRISK3 (exp(beta))",
    y = "Predictor",
    title = "GLM — Significant & Country predictors"
  ) +
  theme(axis.text.y = element_text(size = 7))

###### ---------------- Plotting GLM Random Effect Model (Gender & Country) ---------------- ######
glmm_plot <- glmm_gender_country %>%
  mutate(
    ratio       = exp(Estimate),
    ci.low      = exp(conf.low),
    ci.high     = exp(conf.high),
    significant = if_else(p.adj < 0.05, "Significant", "Not significant")
  )

ord <- glmm_plot$term[order(glmm_plot$ratio)]

ggplot(glmm_plot %>%
         filter(term != "Age.RiskIII"),
       aes(x = ratio,
           y = factor(term, levels = ord),
           colour = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbarh(aes(xmin = ci.low, xmax = ci.high), height = 0.2) +
  geom_point() +
  scale_colour_manual(values = c("Significant" = "firebrick3",
                                 "Not significant" = "black")) +
  labs(
    x = "Ratio of expected QRISK3 (exp(beta))",
    y = "Predictor",
    title   = "GLMM with random intercept for Gender & Country",
    caption = "Model: QRISK3_risk ~ predictor + (1 | Gender) + (1 | Country)\nBlack = p ≥ 0.05, Red = p < 0.05 (BH correction)\n95% CIs on ratio scale\nn varies by predictor\nSkipping z_Fasted.NEFA.mmol.L and z_Fasted.Insulin.pmol.l as only one Country level in subset after drop_na()."
  ) +
  theme(
    axis.text.y = element_text(size = 6),
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 1, vjust = 1, size = 8)
  )


### ----- check whether weight and fat free mass correlate ----- ###
# 1. check distribution to decide for a normality test (using histogram)
df_correlate_fatfreemass_weight <- df_patient_risk_for_lm %>%
  rename(PatientID = Sample_ID) %>%   # rename so the join key matches
  select(PatientID, Fat.free.mass) %>%
  inner_join(qrisk_input %>% select(PatientID, Weight_kg),
             by = "PatientID")

hist(df_correlate_fatfreemass_weight$Fat.free.mass,
     main = "Distribution of Fat-free Mass",
     xlab = "Fat-free mass (kg)")

hist(df_correlate_fatfreemass_weight$Weight_kg,
     main = "Distribution of Weight",
     xlab = "Weight (kg)")

# 2. Pearson Correlation
# data is similar enough to normal distribution --> I will use Pearson correlation
cor.test(
  df_correlate_fatfreemass_weight$Fat.free.mass,
  df_correlate_fatfreemass_weight$Weight_kg,
  method = "pearson",
  use = "complete.obs"
)

plot(df_correlate_fatfreemass_weight$Fat.free.mass,
     df_correlate_fatfreemass_weight$Weight_kg,
     xlab = "Fat-free mass (kg)", ylab = "Weight (kg)",
     main = "Fat-free Mass vs Weight")
abline(lm(Weight_kg ~ Fat.free.mass, data = df_correlate_fatfreemass_weight), col = "red")



###############################################################################
## BMI Calculation (concidering RACE)
###############################################################################
df_bmi <- df_patient_risk_for_lm %>%
  select(PatientID, BMI) %>%
  full_join(qrisk_input %>% select(PatientID, EthnicityCodeQRISK3), by = "PatientID") %>%
  mutate(
    BMI = as.numeric(BMI),
    EthnicityCodeQRISK3 = as.numeric(EthnicityCodeQRISK3),
    bmi_category = case_when(
      is.na(BMI) | is.na(EthnicityCodeQRISK3) ~ NA_character_,
      # Standard thresholds (White/Not stated: codes 1 and 9)
      EthnicityCodeQRISK3 %in% c(1, 9) & BMI >= 40   ~ "Obesity Class 3",
      EthnicityCodeQRISK3 %in% c(1, 9) & BMI >= 35   ~ "Obesity Class 2",
      EthnicityCodeQRISK3 %in% c(1, 9) & BMI >= 30   ~ "Obesity Class 1",
      EthnicityCodeQRISK3 %in% c(1, 9) & BMI >= 25   ~ "Overweight",
      EthnicityCodeQRISK3 %in% c(1, 9) & BMI >= 18.5 ~ "Normal weight",
      EthnicityCodeQRISK3 %in% c(1, 9)               ~ "Underweight",
      # Lower thresholds (South Asian, Black, other minority ethnic groups: codes 2-8)
      EthnicityCodeQRISK3 %in% 2:8 & BMI >= 37.5 ~ "Obesity Class 3",
      EthnicityCodeQRISK3 %in% 2:8 & BMI >= 32.5 ~ "Obesity Class 2",
      EthnicityCodeQRISK3 %in% 2:8 & BMI >= 27.5 ~ "Obesity Class 1",
      EthnicityCodeQRISK3 %in% 2:8 & BMI >= 23   ~ "Overweight",
      EthnicityCodeQRISK3 %in% 2:8 & BMI >= 18.5 ~ "Normal weight",
      EthnicityCodeQRISK3 %in% 2:8                ~ "Underweight",
      TRUE ~ NA_character_
    ),
    bmi_category = factor(bmi_category,
                          levels = c("Underweight", "Normal weight", "Overweight",
                                     "Obesity Class 1", "Obesity Class 2", "Obesity Class 3"))
  )

# Summary table
bmi_summary <- df_bmi %>%
  filter(!is.na(bmi_category)) %>%
  count(bmi_category, .drop = FALSE) %>%
  mutate(pct = round(n / sum(n) * 100, 1))
print(bmi_summary)

# Breakdown by ethnicity group
bmi_by_ethnicity <- df_bmi %>%
  filter(!is.na(bmi_category)) %>%
  mutate(eth_group = ifelse(EthnicityCodeQRISK3 %in% c(1, 9), "Standard (1,9)", "Lower threshold (2-8)")) %>%
  count(eth_group, bmi_category, .drop = FALSE) %>%
  group_by(eth_group) %>%
  mutate(pct = round(n / sum(n) * 100, 1))
print(bmi_by_ethnicity)
