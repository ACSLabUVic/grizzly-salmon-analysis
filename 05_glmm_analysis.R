#-------------------------------------------------------------------------------
# Grizzly bear salmon consumption
# Does salmon consumption increase with bear age
#-------------------------------------------------------------------------------

library(tidyr)
library(lattice)
library(glmmTMB)
library(bbmle) 
library(MuMIn)
library(ggeffects)
library(sjPlot)
library(fitdistrplus)
library(dplyr)
library(corrplot)
library(coefplot)
library(DHARMa)

# ----------------------

##Salmon.med = response variable
##Fixed effects = sex , biomass density, log biomass, biomass, diversity, mean EVI, relative age
## all have been centered and scaled with a mean of 0
##Random effects = individual, year, watershed

#covariote names: 
# sex
# LogBM_z
# EVI_Mean_z
# diversity_z
# relative_age_z

## six models created from hypotheses:

# -----------------------
# read in df for modelling
df <- read.csv(here("data", "processed", "all_cov_scaled_no_locations.csv"))

df$sex <- as.factor(df$sex)

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
m2 <- glmmTMB(Salmon.med ~ LogBM_z + diversity_z + sex + relative_age_z +
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

# ----------------------
## AIC Table:
bbmle::AICtab(
  null, m1, m2, m3, m4, m5,
  weights = TRUE,
  delta   = TRUE,
  sort    = TRUE
)

r.squaredGLMM(m1)
summary(m1)
confint(m1)

# ----------------------
# Plots

# coefficient plot from m1: 
library(dotwhisker)

summary(m1)

coef_data_m1 <- data.frame(
  term = c(
    "Sex (male)",
    "Log salmon biomass",
    "Relative age"
  ),
  estimate = c(
    0.51604,   # sexmale
    0.05827,   # LogBM_z
    0.29132    # relative_age_z
  ),
  std.error = c(
    0.09463,   # sexmale
    0.04584,   # LogBM_z
    0.04835    # relative_age_z
  ),
  model = "m1"
)

#  Create the plot 
coef_plot <- dwplot(
  coef_data_m1,
  dot_args     = list(size = 3),
  whisker_args = list(linewidth = 0.8)
) +
  # Reference line at zero — effects crossing this are non-significant
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  # Colour significant vs non-significant predictors differently
  scale_colour_manual(values = c("m1" = "black")) +
  scale_x_continuous(limits = c(-0.1, 0.7), breaks = seq(-0.1, 0.7, 0.1)) +
  labs(
    x = "Standardised coefficient estimate (± 95% CI)",
    y = ""
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position   = "none",
    axis.text.y       = element_text(size = 13),
    axis.text.x       = element_text(size = 13),
    axis.title.x      = element_text(size = 13),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  )

coef_plot

## --- Generate predictions from m1 on the response scale --- ##
df$predicted <- predict(m1, type = "response")

# Median predicted proportion of salmon in diet by sex
df %>%
  summarise(
    median_consumption = median(predicted, na.rm = TRUE),
    min_consump = min(predicted, na.rm = TRUE),
    max_consump = max(predicted, na.rm = TRUE)
  )
#^ across the whole dataset, this is the range of predicted salmon consumption
# "What is the predicted consumption for an average bear?"
ggpredict(m1, terms = "sex")


#  mean across all observations
#    "What is the average predicted consumption across all bears?"
df %>%
  group_by(sex) %>%
  summarise(mean_predicted = mean(predicted, na.rm = TRUE))

## -----------------------------------------------------------------------------
## --- CREATE FIG 2 --- ##

library(ggeffects)
library(ggplot2)
library(gridExtra)

# ---- Recreate scaling constants ---- #
bm_mean  <- mean(df$LogBM, na.rm = TRUE)
bm_sd    <- sd(df$LogBM,   na.rm = TRUE)
age_mean <- mean(df$relative_age, na.rm = TRUE)
age_sd   <- sd(df$relative_age,   na.rm = TRUE)


# ---- Get marginal predictions for both sexes ---- #
ggpred_biomass <- ggpredict(m1, terms = c("LogBM_z [all]", "sex"))
ggpred_age     <- ggpredict(m1, terms = c("relative_age_z [all]", "sex"))
ggpred_sex     <- ggpredict(m1, terms = "sex")

# ---- Back-transform x values ---- #
ggpred_biomass$x_real <- ggpred_biomass$x * bm_sd  + bm_mean
ggpred_age$x_real     <- ggpred_age$x     * age_sd + age_mean

# ---- Rename group column for clean legend labels ---- #
ggpred_biomass$group <- factor(ggpred_biomass$group, 
                               levels = c("female", "male"), 
                               labels = c("Female", "Male"))
ggpred_age$group     <- factor(ggpred_age$group,     
                               levels = c("female", "male"), 
                               labels = c("Female", "Male"))
ggpred_sex$x         <- factor(ggpred_sex$x,         
                               levels = c("female", "male"), 
                               labels = c("Female", "Male"))

# ---- Shared colour scale ---- #
sex_colours <- c("Female" = "#FF83FA", "Male" = "#1C86EE")

# ---- Plot 1: Log Salmon Biomass ---- #
p1 <- ggplot(ggpred_biomass, aes(x = x_real, y = predicted, 
                                 colour = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = sex_colours) +
  scale_fill_manual(values   = sex_colours) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.4, 1.0, 0.2)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(
    title  = "(a)",
    x      = "Log biomass density",
    y      = "Proportion of salmon in annual diet",
    colour = "Sex", fill = "Sex"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text    = element_text(size = 14),
    axis.title   = element_text(size = 14),
    plot.title   = element_text(size = 16),
    legend.position = "none"  
  )

# ---- Plot 2: Sex ---- #
p2 <- ggplot(ggpred_sex, aes(x = x, y = predicted, colour = x)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 1) +
  scale_colour_manual(values = sex_colours) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.4, 1.0, 0.2)) +
  labs(
    title  = "(b)",
    x      = "Sex",
    y      = NULL,
    colour = "Sex"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text    = element_text(size = 14),
    axis.title   = element_text(size = 14),
    plot.title   = element_text(size = 16),
    legend.position = "none"
  )

# ---- Plot 3: Relative Age ---- #
p3 <- ggplot(ggpred_age, aes(x = x_real, y = predicted, 
                             colour = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = sex_colours) +
  scale_fill_manual(values   = sex_colours) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.4, 1.0, 0.2)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(
    title  = "(c)",
    x      = "Relative age",
    y      = NULL,
    colour = "Sex", fill = "Sex"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text    = element_text(size = 14),
    axis.title   = element_text(size = 14),
    plot.title   = element_text(size = 16),
    legend.position = "right"  # show legend on last panel only
  )

library(patchwork)
library(egg)

final_plot <- ggarrange(p1, p2, p3,
                        nrow   = 1)

## -----------------------------------------------------------------------------

# Get range from marginal predictions
ggpred_range <- ggpredict(m1, 
                          terms = c("LogBM_z [all]", 
                                    "relative_age_z [all]", 
                                    "sex"),
                          bias_correction = TRUE)

range(ggpred_range$predicted)

# Get marginal range holding biomass at mean
ggpred_age <- ggpredict(m1, 
                        terms = c("relative_age_z [all]", "sex"),
                        bias_correction = TRUE)
range(ggpred_age$predicted)
median(ggpred_age$predicted)

## --------------------------------------------------
## --- Calculating the absolute percentage points change: FOR RESULTS **---
library(plotrix)

dff=subset(df,df$sex=="female")
dfm=subset(df,df$sex=="male")

par(mfrow=c(1,2))
minage_z=min(df$relative_age_z)  #minimum age of the bears sampled

# Assign individual coefficients
# Pull intercept and coefficients
coefs <- fixef(m1)$cond  
intercept <- fixef(m1)$cond["(Intercept)"]
age_coef <- summary(m1)$coefficients$cond["relative_age_z", "Estimate"] 
sex_coef <- summary(m1)$coefficients$cond["sexmale", "Estimate"]  #fixed_effects["sexmale"] 


# Back-transform relative_age range to real units
age_min <- min(df$relative_age)   
age_max_f <- max(dff$relative_age) 
age_max_m <- max(dfm$relative_age) 

age_min_z  <- (age_min    - mean(df$relative_age)) / sd(df$relative_age)
age_maxf_z <- (age_max_f  - mean(df$relative_age)) / sd(df$relative_age)
age_maxm_z <- (age_max_m  - mean(df$relative_age)) / sd(df$relative_age)

# Predicted proportions on response scale
female_start <- plogis(intercept + age_coef * age_min_z)
female_end   <- plogis(intercept + age_coef * age_maxf_z)

male_start   <- plogis(intercept + age_coef * age_min_z   + sex_coef)
male_end     <- plogis(intercept + age_coef * age_maxm_z  + sex_coef)

# Print results
cat("=== ABSOLUTE PERCENTAGE POINT CHANGE ===\n\n")

cat("Females:\n")
cat("  Start (relative age", age_min, "):", round(female_start * 100, 1), "%\n")
cat("  End   (relative age", age_max_f, "):", round(female_end   * 100, 1), "%\n")
cat("  Change:", round((female_end - female_start) * 100, 1), "percentage points\n\n")

cat("Males:\n")
cat("  Start (relative age", age_min, "):", round(male_start * 100, 1), "%\n")
cat("  End   (relative age", age_max_m, "):", round(male_end   * 100, 1), "%\n")
cat("  Change:", round((male_end - male_start) * 100, 1), "percentage points\n")

## --------- AT WHAT RELATIVE AGE DOES CONSUMPTION INCREASE BY 15 %  --------
## absolute percentage change ##
# Create a fine sequence of ages
age_seq <- seq(age_min, max(df$relative_age), by = 0.1)
age_seq_z <- (age_seq - mean(df$relative_age)) / sd(df$relative_age)

# Predicted proportions across the full age range
pred_female <- plogis(intercept + age_coef * age_seq_z)
pred_male   <- plogis(intercept + age_coef * age_seq_z + sex_coef)

# Find the age where the increase first hits 15 percentage points
target_pp <- 0.15

age_15pp_female <- age_seq[which(pred_female >= (female_start + target_pp))[1]]
age_15pp_male   <- age_seq[which(pred_male   >= (male_start   + target_pp))[1]]

cat("=== AGE AT WHICH CONSUMPTION INCREASES BY 15 PERCENTAGE POINTS ===\n\n")
cat("Females:\n")
cat("  Start:", round(female_start * 100, 1), "%\n")
cat("  Reaches", round((female_start + target_pp) * 100, 1), "% at relative age:", age_15pp_female, "\n\n")
cat("Males:\n")
cat("  Start:", round(male_start * 100, 1), "%\n")
cat("  Reaches", round((male_start + target_pp) * 100, 1), "% at relative age:", age_15pp_male, "\n")

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

