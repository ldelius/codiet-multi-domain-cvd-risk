### Preparation of urine NMR data for CVD Risk Score Associations
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
df_urine_NMR_data <- read_excel("unified-nmr-targeted-urine_v2.xlsx")

# Prepare the unrine NMR data for the associations
df_urine_NMR_data <- df_urine_NMR_data %>%
  rename_with(~ gsub("[-,]", "_", .x)) %>%
  rename(Sample_ID = patient) %>%
  mutate(Sample_ID = gsub("-", "_", Sample_ID)) %>%
  mutate(Sample_ID = gsub("^(CD)(\\d)", "\\1_\\2", Sample_ID)) %>%
  group_by(Sample_ID) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.x)), .names = "z_{.col}"))

# save as RDS
saveRDS(df_urine_NMR_data, file.path(wkdir, "processed_data", "df_urine_NMR_data.rds"))
