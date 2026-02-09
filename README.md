# Associations with Cardiovascular Risk in the CoDiet Cohort
### An MRes Thesis Project (Biomedical Research – Data Science, Imperial College London)

## Abstract

WIP: Cardiovascular disease (CVD) remains a leading cause of morbidity and mortality worldwide (Add Source here). ... 
This thesis investigates features associated with established CVD risk scores using multi-domain data from the CoDiet cohort, including lipidomics, fatty acids, urine NMR, body composition, clinical risk factors and socio-demographic variables.
We applied generalised liner model (GLM) and elastic net regression approaches to assess associations and predictive performance and identify contributors to CVD risk.

<p align="center">
  <img src="figures/graphical_abstract.png" width="700">
</p>
**Figure 1.** Graphical abstract illustrating the data sources and workflow used in this thesis.
<br><br>
This repository contains the full analysis pipeline used for data preprocessing, model training, evaluation, and visualisation.

## Repository Structure / Mapping Thesis Sections to Code
#### Chapter 3.1 – Cardiovascular Risk Scores

The calculation and aggregation of cardiovascular risk scores reported in Chapter 3.1 of the thesis are implemented as follows:

- **QRISK3 Score**
  - Calculation script: [`QRISK3_calculation.R`](H2_CardiovascularRisk/QRISK3_calculation.R)
  - Input data: [`QRISK3_calculation_input.rds`](Teams_Files/processed_data/QRISK3_calculation_input.rds)
  - *Note:* QRISK3 was calculated separately from the other risk scores, as it was the first score implemented in the project.

- **ASCVD, SCORE2, Framingham, Composite score**
  - Calculation script: [`Cardiovascular_RiskScore_Calculation.R`](H2_CardiovascularRisk/Cardiovascular_RiskScore_Calculation.R)
  - Input data: [`ASCVD_SCORE2_Framingham_input.rds`](Teams_Files/processed_data/ASCVD_SCORE2_Framingham_input.rds)

- **Combined output**
  - [`df_all_risk_scores.rds`](Teams_Files/processed_data/df_all_risk_scores.rds)
  - Contains the final results for all calculated cardiovascular risk scores used in downstream analyses.

- **Figures**
  - All figures presented in Chapter 3.1 (*Cardiovascular Risk Scores*) were generated within  
    `Cardiovascular_RiskScore_Calculation.R`.


## Reproducibility

- R version: X.X.X 
- Key packages: glmnet, tidymodels, XXX
- Random seed: 42  
