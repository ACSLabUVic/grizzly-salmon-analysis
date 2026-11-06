#-------------------------------------------------------------------------------
# Grizzly bear salmon consumption - Biology Letters
# 
# This script fits the beta GLMMs, compares the candidate models (AICc and likelihood-ratio
# tests), checks model fit with DHARMa, and reports the results of the top model (m1).
# 
# 
#-------------------------------------------------------------------------------

library(here)
library(tidyr)
library(glmmTMB)
library(bbmle) 
library(MuMIn)
library(dplyr)
library(DHARMa)
library(ggeffects)
library(insight)
library(performance)

# ----------------------

##  Salmon.med       = response variable
##  Fixed effects    = sex , biomass density, log biomass, biomass, diversity, mean EVI, relative age
##                          all have been centered and scaled with a mean of 0
##  Random effects   = individual, year, watershed

#covariote names: 
# sex
# LogBM_z
# EVI_Mean_z
# diversity_z
# relative_age_z

out_dir <- here("output")
dir.create(out_dir, showWarnings = FALSE)

# -----------------------
# read in df for modelling
df <- read.csv(here("data", "processed", "covariotes_no_locations_2026-11-04.csv"))

df$sex <- as.factor(df$sex)

## data checks: 
nrow(df)
n_distinct(df$individual)
n_distinct(df$WSHDFID)
range(df$Salmon.med)
summary(df$hill_diversity_z)

# ----------------------
### MODEL 1: NULL MODEL  - random + sex + biomass
null <- glmmTMB(Salmon.med ~  + (1 | WSHDFID ) + (1 | year) + (1 | individual), data = df, family = beta_family(link = "logit"))

# H1a: individual and landscape drivers
# salmon consumption is driven by biomass availability, sex, and relative age
m1 <- glmmTMB(Salmon.med ~ sex + LogBM_z + relative_age_z + 
                (1 | WSHDFID ) + (1 | year) + (1 | individual), 
              data = df, family = beta_family(link = "logit"))

# H1b Salmon community complexity: 
# diversity extends seasonal access, driving consumption more than biomass alone (Service et al. 2019)
m2_shannon <- glmmTMB(Salmon.med ~ LogBM_z + diversity_z + sex + relative_age_z +
               (1 | WSHDFID) + (1 | year) + (1 | individual),
              data = df, family = beta_family(link = "logit"))

m2_hill <- glmmTMB(Salmon.med ~ LogBM_z + hill_diversity_z + sex + relative_age_z +
                (1 | WSHDFID) + (1 | year) + (1 | individual),
              data = df, family = beta_family(link = "logit"))

# H2a: Foraging experience
# the effect of biomass availability is modulated by age/experience
m3 <- glmmTMB(Salmon.med ~ sex + LogBM_z + relative_age_z + LogBM_z:relative_age_z + 
                (1 | WSHDFID ) + (1 | year) + (1 | individual), 
              data = df, family = beta_family(link = "logit"))

# H2b: Social dominance   
# age-related increase in consumption differs between sexes
m4 <- glmmTMB(Salmon.med ~ sex + LogBM_z + relative_age_z + relative_age_z:sex + 
                (1 | WSHDFID ) + (1 | year) + (1 | individual), 
              data = df, family = beta_family(link = "logit"))

# H3: Terrestrial habitat quality
# terrestrial vegetation productivity explains additional variation
m5 <- glmmTMB(Salmon.med ~ sex + LogBM_z + EVI_Mean_z + relative_age_z +  
                (1 | WSHDFID ) + (1 | year) + (1 | individual), 
              data = df, family = beta_family(link = "logit"))

# Shannon versus Hill (same number of parameters, so compare by AIC)
bbmle::AICtab(m2_shannon, m2_hill, weights = TRUE, delta = TRUE)
confint(m2_hill, parm = "beta_")

summary(m2_hill)
summary(m2_shannon)

MuMIn::AICc(m2_shannon, m2_hill)                 # exact AICc for both models
confint(m2_shannon, parm = "beta_")["diversity_z", ]
confint(m2_hill, parm = "beta_")["hill_diversity_z", ]

models <- list(null = null, m1 = m1, m2_hill = m2_hill, m3 = m3, m4 = m4, m5 = m5)

# ----------------------
fit_status <- function(m) {
  data.frame(
    convergence = m$fit$convergence,   # 0 = converged
    pdHess      = m$sdr$pdHess,        # TRUE = positive-definite Hessian
    singular    = check_singularity(m),
    nobs        = nobs(m)              # must be the same for every model
  )
}
status <- do.call(rbind, lapply(models, fit_status))
status

# Random-effect variances
lapply(models, VarCorr)

## Model Selection: AICc Table:
bbmle::AICctab(
  null, m1, m2_hill, m3, m4, m5,
  weights = TRUE,
  delta   = TRUE,
  sort    = TRUE
)

# the same table as a data frame (df, logLik, AIC, AICc, delta AICc, weight), saved for the SI
model_selection <- data.frame(
  model  = names(models),
  df     = sapply(models, function(m) attr(logLik(m), "df")),
  logLik = sapply(models, function(m) as.numeric(logLik(m))),
  AIC    = sapply(models, AIC),
  AICc   = sapply(models, MuMIn::AICc),
  row.names = NULL
) %>%
  mutate(dAICc  = AICc - min(AICc),
         weight = exp(-0.5 * dAICc) / sum(exp(-0.5 * dAICc))) %>%
  arrange(dAICc)
model_selection
write.csv(model_selection, file.path(out_dir, "table_postmed_model_selection_NOV4.csv"), row.names = FALSE)

## ---  Likelihood ratio tests: ---
anova(m1, m2_hill)  # m2 = m1 + hill_diversity_z
anova(m1, m3)  # m3 = m1 + biomass-age interaction
anova(m1, m4)  # m4 = m1 + sex by age interatcion
anova(m1, m5)  # m5 = m1 + EVI

## drop 1 term from m1 at a time
m1_no_sex <- update(m1, . ~ . - sex)
m1_no_biomass <- update(m1, . ~ . - LogBM_z)
m1_no_age <- update(m1, . ~ . - relative_age_z)

anova(m1_no_sex, m1)      # LRT for sex
anova(m1_no_biomass, m1)  # LRT for LogBM_z
anova(m1_no_age, m1)      # LRT for relative_age_z

## --- DHARMA model diagnostics: m1 ---
set.seed(123)

## Response range 
range(model.frame(m1)$Salmon.med)

## simulate scaled residuals: 
sim_m1 <- simulateResiduals(fittedModel = m1, n = 1000)

## QQ plot (left) and residuals vs predicted (right)
plot(sim_m1)

## formal test
testUniformity(sim_m1)
testDispersion(sim_m1)
testOutliers(sim_m1)
check_collinearity(m1)

## residuals against each predictor
mf <- model.frame(m1)               # keeps rows aligned with the model
plotResiduals(sim_m1, mf$LogBM_z)
plotResiduals(sim_m1, mf$relative_age_z)
plotResiduals(sim_m1, mf$sex)
# random
plotResiduals(sim_m1, mf$WSHDFID)
plotResiduals(sim_m1, mf$year)

#save for SI:
png(file.path(out_dir, "FigS_DHARMa_m1.png"), width = 2400, height = 1200, res = 300)
plot(sim_m1)
dev.off()

png(file.path(out_dir, "FigS_DHARMa_predictors.png"), width = 3000, height = 1000, res = 300)
par(mfrow = c(1, 3))
plotResiduals(sim_m1, mf$LogBM_z)
plotResiduals(sim_m1, mf$relative_age_z)
plotResiduals(sim_m1, mf$sex)
dev.off()
## 
lapply(list(m1 = m1, m2_hill = m2_hill, m5 = m5), check_collinearity)

## --- m1 results: ---
summary(m1)

ci_m1 <- confint(m1)
ci_m1
write.csv(data.frame(term = rownames(ci_m1), ci_m1, row.names = NULL),
          file.path(out_dir, "tableS3_m1_coefficients.csv"), row.names = FALSE)

r2_nakagawa(m1)

# Partial (drop-one) marginal R2
marg <- function(m) r2_nakagawa(m)$R2_marginal
data.frame(
  term       = c("sex", "LogBM_z", "relative_age_z"),
  R2_without = c(marg(m1_no_sex), marg(m1_no_biomass), marg(m1_no_age)),
  R2_full    = marg(m1)
) %>% mutate(partial_R2 = R2_full - R2_without)


## --- Generate predictions from m1 on the response scale --- ##
## used as model check only
pred_check <- df %>% mutate(predicted = predict(m1, type = "response"))

# Median predicted proportion of salmon in diet by sex
pred_check %>%
  summarise(
    median_consumption = median(predicted, na.rm = TRUE),
    min_consump = min(predicted, na.rm = TRUE),
    max_consump = max(predicted, na.rm = TRUE)
  )

# Predicted consumption for an average bear?"
ggpredict(m1, terms = "sex")


# the average predicted consumption across all bears?"
pred_check %>%
  group_by(sex) %>%
  summarise(mean_predicted = mean(predicted, na.rm = TRUE))

#-------------------------------------------------------------------------------
# Check how much watershed random effect matters
# by comparing model with and without it
m1_no_watershed <- glmmTMB(
  Salmon.med ~ sex + LogBM_z + relative_age_z +
    (1 | year) + (1 | individual),
  data   = df,
  family = beta_family(link = "logit")
)

# Compare
anova(m1_no_watershed, m1)

# Check if coefficients change meaningfully
fixef(m1)$cond
fixef(m1_no_watershed)$cond

## 
sum(df$Salmon.med == 0)
sum(df$Salmon.med == 1)
range(df$Salmon.med)


