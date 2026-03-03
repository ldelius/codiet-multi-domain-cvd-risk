# Multi-Domain Analysis of Predictors for Cardiovascular Risk Estimation in an Overweight and Obese Population
### An MRes Thesis Project (Biomedical Research – Data Science, Imperial College London)

## Abstract

Cardiovascular disease (CVD) remains the leading cause of mortality worldwide. Overweight and obesity are major and increasingly prevalent risk factors, yet established CVD risk scores do not adequately capture their risk and may underestimate CVD risk in this population. This work examines a multi-domain dataset comprising metabolomics (urinary metabolites, lipid profiles, fatty acids), clinical measurements, body composition metrics, and sociodemographic factors for their association with CVD risk and predictive utility, as assessed by four established 10-year CVD risk scores (ASCVD, Framingham, QRISK3, SCORE2). The overall study design is summarised in Figure 1.

<p align="center">
  <img src="figures/graphical_abstract.png" width="700">
</p>
**Figure 1.** Overview of the study design.
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
  - All figures presented in Chapter 3.1 (*Cardiovascular Risk Scores*) were generated within [`Cardiovascular_RiskScore_Calculation.R`](H2_CardiovascularRisk/Cardiovascular_RiskScore_Calculation.R).
  - The supplementary figure showing CVD risk score input variables and CVD Score results was generated within the same script.

#### Predictor Preparation

Predictor datasets used for GLM, Regression Model Comparison, and Elastic Net Analysis were preprocessed separately from outcome generation to ensure consistency across modelling approaches.  
All predictor preparation steps are documented in dedicated preprocessing scripts, with analysis-ready datasets stored as `.rds` files.

- **Sociodemographic Predictors**
  - Preprocessing script: [`QRISK3_REDcap_random.R`](H2_CardiovascularRisk/QRISK3_REDcap_random.R)
Two analysis-ready datasets were generated for different modelling frameworks:
  - GLM: [`df_REDcap_demographics_predictor.rds`](Teams_Files/processed_data/df_REDcap_demographics_predictor.rds) (Contains z-scaled sociodemographic predictors with variables and observations excluded based on small sample sizes.)
  - Elastic Net Regression: [`df_REDcap_demographics_ElaNet.rds`](Teams_Files/processed_data/df_REDcap_demographics_ElaNet.rds) (Contains cleaned but non–z-scaled sociodemographic predictors without exclusion.)

- **Risk Factor Dataset**
  - This dataset comprises sociodemographic variables, body composition metrics, clinical measurements, and selected faecal and urinary metabolites.
  - Preprocessing script: [`QRISK3_calculation.R`](H2_CardiovascularRisk/QRISK3_calculation.R)
  - The resulting analysis-ready dataset is stored as: [`df_risk_factor_predictors.rds`](Teams_Files/processed_data/df_risk_factor_predictors.rds) (Contains z-scaled continuous variables, processed factor variables, and the transformed urine tpred score.)
  - Predictors from this dataset that were conceptually better aligned with other data domains were joint those datasets for the GLM and elastic net regression analyses.

- **Fatty Acid and Lipidomics Predictors**
  - Preprocessing script: [`QRISK3_lipidomics.R`](H2_CardiovascularRisk/QRISK3_lipidomics.R)
  - The resulting analysis-ready datasets are stored as (Outliers were removed using the 1st/99th percentile thresholds and continuous variables were z-scaled.):
    - Fatty acids: [`df_fatty_acids_predictor_statin_suppl.rds`](Teams_Files/processed_data/df_fatty_acids_predictor_statin_suppl.rds)  
    - Lipidomics: [`df_lipidomics_predictor_statin_suppl.rds`](Teams_Files/processed_data/df_lipidomics_predictor_statin_suppl.rds)

- **Urine NMR Predictors**
  - Preprocessing script: [`urine_nmr_data.R`](H2_CardiovascularRisk/urine_nmr_data.R)
  - The resulting z-scaled, analysis-ready dataset is stored as: [`df_urine_NMR_data.rds`](Teams_Files/processed_data/df_urine_NMR_data.rds)

- **Body Composition Predictors**
  - Preprocessing script: [`Body_composition_metrics.R`](H2_CardiovascularRisk/Body_composition_metrics.R)
  - The resulting z-scaled, analysis-ready dataset is stored as: [`df_body_composition_metrics.rds`](Teams_Files/processed_data/df_body_composition_metrics.rds)

#### Chapter 3.2 – GLM Analysis

Associations between predictors and CVD scores reported in Chapter 3.2 of the thesis were analysed using generalised linear models.
[`FixedEffectAssociations.R`](H2_CardiovascularRisk/FixedEffectAssociations.R) includes:
- Data preparation, model fitting, and output generation for the GLM analyses
- Figure generation for figures presented in Chapter 3.2
- Produces supplementary tables summarising covariate effects and full GLM results
- Model Evaluation and respective table

The same analysis, additionally adjusting for Age as a fixed effect, was conducted in: [`FixedEffectAge.R`](H2_CardiovascularRisk/FixedEffectAge.R).


#### Chapter 3.3 – ML Model comparison analysis

ML comparison that was used in thesis is stored in [`ML_comp_parameter_tuning.R`](H2_CardiovascularRisk/ML_comp_parameter_tuning.R). Previous comparison with 20-hould out participants and all datasets (not shown in thesis) are stored in [`tidymodels_comparison.R`](H2_CardiovascularRisk/tidymodels_comparison.R).
- Key pipline parameter: XXX


#### Chapter 3.4 – Elastic Net Regression Analysis
The File XY includes:
- Data preparation, model fitting, and output generation
- Figure generation for figures in Elastic Net Regression Analysis Chapter

## Reproducibility

- R version: X.X.X 
- Key packages: glmnet, tidymodels, XXX
- Random seed: 42  
