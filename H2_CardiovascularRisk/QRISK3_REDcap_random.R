### QRISK3 and REDcap data set columns
### Author: Luisa Delius

# Load packages
library(tidyverse)
library(broom)
library(patchwork)
library(QRISK3)
library(readxl)

# Set working directory
wkdir <- "/Users/luisadelius/Documents/Code/project_one/Teams_Files"
knitr::opts_knit$set(root.dir = wkdir)
setwd(wkdir)

# upload data
df_other_compare <- read_excel("QRISK3_data.xlsx", sheet = "other_compare")
QRISK3_sample_ID <- readRDS("QRISK3_sample_ID.RDS") ### this has to be exchanged if i change anything of the QRISK input data of courese!!!


# Prepare/clean the data for lm and join with QRISK3 score
df_other_cleaned_and_QRISK <- df_other_compare %>%
  rename(Sample_ID = volunteer_id) %>%   # rename ID column
  mutate(
    Sample_ID = str_replace_all(Sample_ID, "-", "_"),   # replace '-' with '_' in IDs
 
 #  Some general cleaning  
    Marital_status = factor(Marital_status),
    Living_Status = str_trim(Living_Status),
    Employment_Status = str_trim(Employment_Status),
    Working_time = str_trim(Working_time),
    naps_during_day = str_to_lower(naps_during_day),

 #  Wine per Week    
    wine_per_week = str_replace_all(wine_per_week, ">", "≥"),
    wine_per_week = factor(wine_per_week, levels = c("<7 glasses", "≥7 glasses")),
  
 #  Olive Oil per Day  
    olive_oil_given_day = str_replace_all(
      olive_oil_given_day,
      c( 
        "tablespoons" = "tbsp",   # mutate here the various inputs like <4 tablespoons <4 tbsp
        "Tablespoons" = "tbsp",
        " " = "",      # remove extra spaces
        ">" = "≥"      # replace > with ≥
      )
    ),
    olive_oil_given_day = factor(olive_oil_given_day,
                                levels = c("<4tbsp", "≥4tbsp")),
 
 #  Sleep hours on workdays   
    hours_sleep_workdays = str_replace_all(hours_sleep_workdays, ",", "."),
    hours_sleep_workdays = map_chr(hours_sleep_workdays, function(x) {
      # Handle missing values first
      if (is.na(x) || x == "") {
        return(NA_character_)
      }
      # Remove extra spaces
      x <- str_trim(x)
      # If it's a range, take the mean
      if (str_detect(x, "[-–]")) {
        parts <- as.numeric(unlist(str_split(x, "[-–]")))
        mean(parts, na.rm = TRUE)
      } else {
        x
      }
    }),
    hours_sleep_workdays = as.numeric(hours_sleep_workdays),  # Convert everything to numeric
   
 #  anualy net salary
    annual_net_salary = str_replace_all(
      annual_net_salary,
      c("," = "", "£" = "", "\\s+" = "") # remove , and GBP, and space with jsut nothing --> so everything should be in the same shape
    ),
    annual_net_salary = factor(annual_net_salary,
                               levels = c(
                                 "0-10320",
                                 ">10320-17200",
                                 ">17200-25800",
                                 ">25800-34400",
                                 ">34400-43000",
                                 ">43000-51600",
                                 ">51600"
                               ),
                               ordered = TRUE),

 #  Education level
 Education_Level = factor(Education_Level,
                          levels = c(
                            "Illiterate or primary education",
                            "Secondary",
                            "Short-cycle tertiary education",
                            "University (bachelor or equivalent)",
                            "Doctoral"),
                          ordered = TRUE),
 #  Nuts per week   
 servings_nuts_per_week = factor(servings_nuts_per_week,
                                 levels = c("<3", "≥3"))
  ) %>%
  
 # create factors that do not need specific level  
  mutate(
    across(
      c(recruitment_site, Marital_status, Living_Status, Employment_Status,
        Working_time, naps_during_day, olive_oil_as_main_culinary_fat,
        servings_nuts_per_week),
      as.factor)) %>%
    
# z-scaling
  mutate(across(where(is.numeric), ~ as.numeric(scale(.x)), .names = "z_{col}")) %>%

 # join with the QRISK3 score 
  inner_join(QRISK3_sample_ID, by = "Sample_ID")


# Remove the following
df_other_cleaned_and_QRISK <- df_other_cleaned_and_QRISK %>%
  filter(
    Employment_Status != "Student",
    Working_time != "Shift worker",
    Marital_status != "Widowed"
  )

saveRDS(df_other_cleaned_and_QRISK,
  file = file.path(wkdir, "processed_data", "df_REDcap_demographics_predictor.rds"))

#####Preperaing and running the glm #######
## i will run one GLM per numeric predictor
## in order to not loose a level in my forest plots, I will run a GLM per factor level 

# 1. Identify predictors
numeric_predictors <- df_other_cleaned_and_QRISK %>%
  select(
    starts_with("z_")) %>%
    names()

factor_predictors <- df_other_cleaned_and_QRISK %>%
  select(where(is.factor)) %>%
  names() #creates a simple list of column names

# 2. GLMs for numeric predictors
glm_output_num <- map_dfr(numeric_predictors, function(var) { # loops over all the predictors, fits a seperate model each and bins together
  df_pair <- df_other_cleaned_and_QRISK %>%
    select(QRISK3_risk, all_of(var)) %>%
    drop_na()
  
  model <- glm(
    as.formula(paste("QRISK3_risk ~", var)),
    data   = df_pair,
    family = Gamma(link = "log") # chose Gamma distribution given the positive, continuous, right-skewed data.
  ) # why link = log? Gamma distribution only makes sense for positive means, so by log link the predicted means are all positive. Each 1-unit increase in X increases the expected QRISK3 by e^b1.
  
  # Return tidy estimates; exponentiate to get multiplicative effects (ratios)
  tidy(model, conf.int = TRUE, exponentiate = TRUE) %>% # exponentiates b1 (logarithm of the multiplicative effect) --> easier to read and interpret
    filter(term != "(Intercept)") %>%
    mutate(
      predictor = var,
      n = nrow(df_pair)
    )
}) 

# 3. GLM factor predictors
glm_output_fac <- map_dfr(factor_predictors, function(var) {
  df_pair <- df_other_cleaned_and_QRISK %>%
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
        term      = paste0(var, " = ", lev),  # <— THIS is the important line
        n         = nrow(df_tmp)
      )
  })
})

# 4. Combine numeric + factor results
glm_output <- bind_rows(glm_output_num, glm_output_fac) %>%
  mutate(
    term_pretty = if_else(
      is.na(level),
      predictor,
      paste0(predictor, " = ", level)
    )
  ) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    significant = if_else(conf.low > 1 | conf.high < 1,
                          "Significant", "Not significant")
  ) %>%
  arrange(p.adj)


# ----- Plotting my GLM models (Gamma, log link) -----
# 1. do single plots
plot_predictors <- c(numeric_predictors, factor_predictors)

make_plot_glm <- function(var) {
  d <- df_other_cleaned_and_QRISK %>%
    select(QRISK3_risk, all_of(var)) %>%
    drop_na()
  
  gg <- ggplot(d, aes(x = .data[[var]], y = QRISK3_risk)) +
    labs(x = var, y = "QRISK3_risk")
  
  if (is.factor(d[[var]])) {
    gg + geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.5, size = 1.4) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    gg + geom_point(size = 1.6, alpha = 0.7) +
      geom_smooth(method = "glm",
                  method.args = list(family = Gamma(link = "log")),
                  se = TRUE)
  }
}

# build all plots
plots_glm <- purrr::map(plot_predictors, make_plot_glm)

# show 6 per page and print all pages:
lapply(split(plots_glm, ceiling(seq_along(plots_glm) / 6)),
       \(p) print(wrap_plots(p, ncol = 3)))

## 2. Forest plot
ggplot(
  glm_output %>% 
    filter(is.finite(estimate), is.finite(conf.low), is.finite(conf.high)) %>% 
    arrange(estimate) %>% 
    mutate(term = factor(term, levels = term)),
  aes(x = estimate, y = term, colour = significant)
) +
  geom_vline(xintercept = 1, linetype = "dashed") +
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
