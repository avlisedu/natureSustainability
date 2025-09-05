########## Nature Sustainability ########## 
################################################################################
# Title: Statistical modeling of Dengue cases in the state of Pernambuco, Brazil: 
#         an approach using panel data and Zero-Inflated Negative Binomial Regression (ZINB)
# Description: This script performs data processing, exploratory analysis, and
#              estimation of Zero-Inflated Negative Binomial (ZINB) models 
#              with fixed and random effects to investigate the determinants 
#              of dengue incidence across municipalities in Pernambuco State.
#              The analysis includes variance inflation factor (VIF) testing,
#              equidispersion testing (Cameron & Trivedi), and regional model 
#              estimation.
#
# Authors: Eduardo Silva, Erlandson Ferreira Saraiva and Maisa Mendonça Silva 
# Affiliation: Graduate Program in Production Engineering (PPGEP/UFPE)
#              Research Group on Information and Decision Systems (GPSID)
#
# Date: 2025-08-12
# Language: R (Version ≥ 4.2)
# Required packages: See 'packages' section below
#
# Notes:
# - All variables are standardized according to the dataset documentation.
# - Results are reported in the Supplementary Materials of the manuscript
#   submitted to Nature Sustainability.
################################################################################

# Required Packages ----------------------------------------------------------------
packages <- c("glmmTMB", "stargazer","readr", "readxl", "lme4", "MASS", 
              "plm", "nortest", "car", "GGally", "ggplot2", "dplyr", 
              "pscl", "performance", "lmtest", "jtools", "sjPlot", "overdisp", "ggcorrplot", "DHARMa", "AER")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Data Import and Processing ----------------------------------------------------------------
dataset <- read_excel("dataset_dengue.xlsx", sheet = "panel", col_names = TRUE)

filtered_panel_data <- dataset %>%
  dplyr::filter(!is.na(dengue_cases) & !is.na(solid_waste_per_1000) & !is.na(temp_max) & !is.na(temp_min) & 
                  !is.na(precipitation) & !is.na(basic_sanitation_plan) & !is.na(humidity) & 
                  !is.na(water_coverage))

final_panel_data <- pdata.frame(filtered_panel_data, index = c("municipality", "yyyymm"))

head(final_panel_data)
summary(final_panel_data)

# Creating lag for precipitation only
final_panel_data <- final_panel_data %>%  
  group_by(municipality) %>%  
  mutate(precipitation_lagOne = lag(precipitation, n = 1)) %>% 
  ungroup() 

# Variance Inflation Factor (VIF) analysis
mod_vif <- lm(dengue_cases ~ temp_max + temp_min + humidity + precipitation +
                solid_waste_per_1000 + waste_collection + sewage_coverage +
                water_coverage + basic_sanitation_plan + community_health_agent_coverage,
              data = final_panel_data)

vif_df <- data.frame(
  Variable = names(vif(mod_vif)),
  VIF = round(vif(mod_vif), 2),
  row.names = NULL
)
print(vif_df)

# Variance Inflation Factor (VIF) analysis indicated no significant multicollinearity among the independent variables (all VIFs < 1.51), ensuring the stability and interpretability of the model coefficients.

# Equidispersion Test ----------------------------------------------------------------
overdisp(
  x = final_panel_data, 
  dependent.position = 10, 
  predictor.position = c(12, 13, 14, 15, 16, 17, 18, 19, 20, 11)
)
# The overdispersion test was performed using the Cameron & Trivedi (1990) 
# approach on the final panel dataset. The test yielded a 
# Lambda t-statistic of 6.50 (p-value < 0.001), 
# providing strong evidence against the null hypothesis of equidispersion.

# Zero-Inflated Negative Binomial ----------------------------------------------------------------
# Complete Fixed Effects Model
model_01 <- glmmTMB(
  dengue_cases ~ temp_max + temp_min + humidity + precipitation_lagOne,
  ziformula = ~ waste_collection + sewage_coverage + water_coverage + 
    as.factor(basic_sanitation_plan) + community_health_agent_coverage + solid_waste_per_1000,
  family = nbinom2,
  offset = log(population),
  data = final_panel_data
)

summary(model_01)
confint(model_01)
model_01_residuals <- simulateResiduals(fittedModel = model_01)
plotQQunif(model_01_residuals)

# Complete Random Effects Model
model_02 <- glmmTMB(
  dengue_cases ~ temp_max + temp_min + humidity + precipitation_lagOne +
    (1 | municipality),
  ziformula = ~ waste_collection + sewage_coverage + water_coverage +
    as.factor(basic_sanitation_plan) + community_health_agent_coverage + solid_waste_per_1000,
  family = nbinom2,
  offset = log(population),
  data = final_panel_data
)

summary(model_02)
confint(model_02)
model_02_residuals <- simulateResiduals(fittedModel = model_02)
plotQQunif(model_02_residuals)
# AIC | BIC: 12263.9  12355.9 < Best model

# Random Effects – Only significant variables
model_03 <- glmmTMB(dengue_cases ~ temp_max + precipitation_lagOne +
                      (1 | municipality/yyyymm),
                    ziformula = ~ water_coverage + solid_waste_per_1000,
                    family = nbinom2,
                    offset = log(population),
                    data = final_panel_data)
summary(model_03)
confint(model_03)
model_03_residuals <- simulateResiduals(fittedModel = model_03)
plot(model_03_residuals)
# AIC | BIC: 12288.2  12343.4

# Random Effects – Only significant variables with structural change
model_zinbase <- glmmTMB(Dengue ~ Max + Prec_lag1 + 
                           (1 | Municipality/yyyymm),
                         ziformula = ~ Water + Solid_Waste_per_1000 + Max + Prec_lag1,
                         family = nbinom2,
                         offset = log(Population),
                         data = panel_data)
summary(model_zinbase)
confint(model_zinbase)
# AIC | BIC: 12251.8  12319.3
model_zinbase <- simulateResiduals(fittedModel = model_zinbase)
plot(model_zinbase)

# Zero-inflation probability
options(scipen = 999)
prob_data <- data.frame(
  Municipality = model_zinbase$frame$Municipality,  
  yyyymm = model_zinbase$frame$yyyymm,  
  Prob_ZINBASE = round(predict(model_zinbase, type = "zprob"), 4)  
)
write.csv(prob_data, "probability.csv", row.names = FALSE, fileEncoding = "UTF-8")

# Regional Analysis - using the zinbase model ----------------------------------------------------------------
# Split dataset by region
unique_regions <- unique(final_panel_data$region)
for (reg in unique_regions) {
  if (!is.na(reg)) {
    filtered_df <- final_panel_data %>% filter(region == reg)
    assign(paste0("panel_data_", reg), filtered_df)
    print(paste("Dataset created:", paste0("panel_data_", reg)))
  }
}

# ============= Region 0 - Noronha Island
# Not available

# ============= Region 1 - Hinterland of Pernambuco
model_zinbase1 <- glmmTMB(Dengue ~ Max + Prec_lag1 + 
                            (1 | Municipality/yyyymm),
                          ziformula = ~ Water + Solid_Waste_per_1000 + Max + Prec_lag1,
                          family = nbinom2,
                          offset = log(Population),
                          data = panel_data_1)
summary(model_zinbase1)
confint(model_zinbase1)
# AIC | BIC: 3017.7   3071.2
model_zinbase <- simulateResiduals(fittedModel = model_zinbase)
plot(model_zinbase)

# ============= Region 2 - São Francisco Region of Pernambuco
model_zinbase2 <- glmmTMB(Dengue ~ Max + Prec_lag1 + 
                            (1 | Municipality/yyyymm),
                          ziformula = ~ Water + Solid_Waste_per_1000 + Max + Prec_lag1,
                          family = nbinom2,
                          offset = log(Population),
                          data = panel_data_2)
summary(model_zinbase2)
confint(model_zinbase2)
# AIC | BIC: 1465.8   1510.5
model_zinbase <- simulateResiduals(fittedModel = model_zinbase)
plot(model_zinbase)

# ============= Region 3 - Agreste of Pernambuco
model_zinbase3 <- glmmTMB(Dengue ~ Max + Prec_lag1 + 
                            (1 | Municipality/yyyymm),
                          ziformula = ~ Water + Solid_Waste_per_1000 + Max + Prec_lag1,
                          family = nbinom2,
                          offset = log(Population),
                          data = panel_data_3)
summary(model_zinbase3)
confint(model_zinbase3)
# AIC | BIC: 1465.8   1510.5
model_zinbase <- simulateResiduals(fittedModel = model_zinbase)
plot(model_zinbase)

# ============= Region 4 - Forest Zone of Pernambuco
model_zinbase4 <- glmmTMB(Dengue ~ Max + Prec_lag1 + 
                            (1 | Municipality/yyyymm),
                          ziformula = ~ Water + Solid_Waste_per_1000 + Max + Prec_lag1,
                          family = nbinom2,
                          offset = log(Population),
                          data = panel_data_4)
summary(model_zinbase4)
confint(model_zinbase4)
# AIC | BIC: Did not converge
model_zinbase <- simulateResiduals(fittedModel = model_zinbase)
plot(model_zinbase)

# ============= Region 5 - Recife Metropolitan Area
model_zinbase5 <- glmmTMB(Dengue ~ Max + Prec_lag1 + 
                            (1 | Municipality/yyyymm),
                          ziformula = ~ Water + Solid_Waste_per_1000 + Max + Prec_lag1,
                          family = nbinom2,
                          offset = log(Population),
                          data = panel_data_5)
summary(model_zinbase5)
confint(model_zinbase5)
model_zinbase <- simulateResiduals(fittedModel = model_zinbase)
plot(model_zinbase)

