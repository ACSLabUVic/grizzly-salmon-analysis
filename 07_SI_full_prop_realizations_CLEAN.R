###
# Full error propagation (Supplementary Information sensitivity analysis)
#
# Each realization draws, jointly:
#    1. a whole MixSIAR posterior iteration of salmon proportion for every bear-year;
#    2. an age at first detection for every bear (scaled beta distribution whose mean and
#       variance decline over the study, as in the original realization code);
#  then refits the top beta GLMM with absolute age in place of relative age. From each refit
#  we simulate fixed effects (GLMM estimation uncertainty) and, for each simulation, draw
#    3. von Bertalanffy mass-at-age parameters (Kingsley et al. 1983, Table 1: SEs and
#       error correlations), and
#    4. intake-equation parameters (% salmon -> kg salmon per kg bear).
#  All draws are pooled; intervals are the 2.5th and 97.5th percentiles of the pooled draws.
#  Sources are added one at a time so the contribution of each can be seen.
#
# Inputs : data/processed/all_cov_scaled_no_locations.csv      (written by script 04)
#          data/processed/gb_all_posteriors_no_location_2026-10-20.csv  (MixSIAR draws, no coordinates)
#          data/raw/Deacy_Robbins_intake_data.csv              (intake data; check permission to share)
#          output/main_text_draws.rds                          (written by script 06; run 06 first)
# Outputs: output/tableS_*.csv, output/FigS20*, FigS21*, FigS30*, output/SI_summary.rds
#
# No location data are used in this script.
#
# Quick test:  Sys.setenv(N_REAL = 50)  before running (the default is 2000 refits)
###

library(here)
library(glmmTMB)
library(MASS)       # mvrnorm
library(parallel)
library(tidyr)
library(ggplot2)
library(patchwork)
library(dplyr)      # load dplyr after MASS so that dplyr::select and dplyr::filter are not masked


## ---- Settings ---- ##
set.seed(20261028)
n_real         <- as.integer(Sys.getenv("N_REAL", 2000))  # realizations (one GLMM refit each)
n_sims_per_fit <- 5                                       # fixed-effect draws per refit
n_cores        <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 2)  # mclapply needs 1 core on Windows
probs          <- c(0.025, 0.5, 0.975)
ages_pred      <- 3:25        # absolute ages for prediction

out_dir <- here("output")     # same folder as 06
dir.create(out_dir, showWarnings = FALSE)

sex_cols <- c(female = "#0072B2", male = "#E69F00")
theme_ms <- theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.tag = element_text(face = "plain"),
        plot.tag.location = "panel",
        plot.tag.position = c(0.065, 1.025))


## ---- 1. DATA ---- ##
# Covariates (550 bear-years), written by script 04; Salmon.med is the MixSIAR posterior median
df <- read.csv(here("data", "processed", "all_cov_scaled_no_locations.csv"))
df$sex <- factor(df$sex, levels = c("female", "male"))

# Full MixSIAR posteriors for salmon proportion (638 bear-years x 3000 iterations), no coordinates
post <- data.table::fread(here("data", "processed", "gb_all_posteriors_no_location_2026-10-20.csv"),
                          data.table = FALSE)
draw_cols <- grep("^iteration_", names(post))
stopifnot(length(draw_cols) > 0,
          !anyDuplicated(post[, c("individual", "year")]))

# Align posterior rows to the 550 bear-years used in the GLMM
idx <- match(paste(df$individual, df$year), paste(post$individual, post$year))
stopifnot(!anyNA(idx))
Y <- as.matrix(post[idx, draw_cols])   # bear-years x iterations, rows aligned with df
rm(post)

# Checks: draws reproduce the posterior medians used in the main analysis,
# and every draw is strictly inside (0, 1) as the beta likelihood requires
stopifnot(max(abs(apply(Y, 1, median) - df$Salmon.med)) < 0.02)
stopifnot(all(Y > 0 & Y < 1))

# Scaling constants used to build the *_z predictors in script 04
bm_center <- mean(df$LogBM); bm_scale <- sd(df$LogBM)
stopifnot(max(abs((df$LogBM - bm_center) / bm_scale - df$LogBM_z)) < 1e-6)


## ---- 2. AGE AT FIRST DETECTION: REALIZATIONS ---- ##
## Age at first detection is a whole number of years between 3 and 15, drawn as
## 3 + 12 * Beta(a, b), rounded. The mean (7 yr) and variance (10 yr^2) apply to bears first
## detected 2009-2011, then decline linearly to a mean of 4 yr and variance of 3 yr^2 for bears
## first detected in 2019 (later first detections are more likely to be young bears).
first_year <- ave(df$year, df$individual, FUN = min)
stopifnot(all(df$relative_age == df$year - first_year + 1))
years_since_first <- df$relative_age - 1

min_age_first <- 3; max_age_first <- 15
beta_shape <- function(fy) {
  t  <- pmax(0, fy - 2011)
  mu <- 7 + (4 - 7) / 8 * t                     # mean, years
  v  <- 10 + (3 - 10) / 8 * t                   # variance, years^2
  m  <- (mu - min_age_first) / (max_age_first - min_age_first)
  vs <- v / (max_age_first - min_age_first)^2
  a  <- m^2 * (1 - m) / vs - m
  cbind(a = a, b = a * (1 - m) / m)
}
ind_first <- df %>% distinct(individual) %>%
  mutate(first_year = first_year[match(individual, df$individual)])
shapes <- beta_shape(ind_first$first_year)
stopifnot(all(shapes > 0))
draw_age_first <- function() {
  round(min_age_first + rbeta(nrow(ind_first), shapes[, "a"], shapes[, "b"]) *
          (max_age_first - min_age_first))
}

# All realizations are drawn here, in the main process, so results are reproducible
age_first_mat <- replicate(n_real, draw_age_first())               # individuals x realizations
row_ind       <- match(df$individual, ind_first$individual)
age_mat       <- age_first_mat[row_ind, , drop = FALSE] + years_since_first   # bear-years x realizations
draw_idx      <- sample(ncol(Y), n_real, replace = n_real > ncol(Y))

# One fixed scaling for age across all realizations, so coefficients are comparable
age_center <- mean(age_mat); age_scale <- sd(age_mat)

# Check the realized age distributions
age_first_check <- data.frame(first_year = rep(ind_first$first_year, n_real),
                              age_first  = as.vector(age_first_mat)) %>%
  group_by(first_year) %>%
  summarise(n_bears = n() / n_real, mean = mean(age_first), var = var(age_first),
            p_age3 = mean(age_first == 3), p_ge10 = mean(age_first >= 10), .groups = "drop")
print(age_first_check, digits = 3)
cat("Realized bear-year ages: median", median(age_mat), "; 97.5th pct", quantile(age_mat, 0.975),
    "; share > 20 yr", round(mean(age_mat > 20), 4), "; max", max(age_mat), "\n")
write.csv(age_first_check, file.path(out_dir, "tableS_age_at_first_detection_realized.csv"), row.names = FALSE)


## ---- 3. GLMM REFITS (cached) ---- ##
fit_glmm_age <- function(r) {
  d <- df
  d$y     <- Y[, draw_idx[r]]
  d$age_z <- (age_mat[, r] - age_center) / age_scale
  m <- glmmTMB(y ~ sex + LogBM_z + age_z + (1 | WSHDFID) + (1 | year) + (1 | individual),
               data = d, family = beta_family(link = "logit"),
               control = glmmTMBControl(parallel = 1))
  list(beta = fixef(m)$cond, vcov = vcov(m)$cond,
       ok = isTRUE(m$sdr$pdHess) && all(is.finite(vcov(m)$cond)))
}

# the cache file name carries a hash of everything the refits depend on, so changed data
# can never be silently read from an old cache
data_hash <- substr(rlang::hash(list(df[, c("individual", "year", "sex", "WSHDFID", "LogBM_z", "relative_age")],
                                     n_real, draw_idx, age_center, age_scale)), 1, 8)
fits_file <- file.path(out_dir, paste0("fits_SI_age_mixsiar_", n_real, "_", data_hash, ".rds"))
if (file.exists(fits_file)) {
  fits <- readRDS(fits_file)
  stopifnot(length(fits) == n_real)
} else {
  fits <- mclapply(seq_len(n_real), fit_glmm_age, mc.cores = n_cores)
  saveRDS(fits, fits_file)
}
ok <- vapply(fits, function(f) is.list(f) && isTRUE(f$ok), logical(1))
cat("Refits converged:", sum(ok), "of", length(ok), "\n")
if (mean(ok) < 0.95) warning("Fewer than 95% of the refits converged; check the models before using the results")
fits <- fits[ok]

B_full <- do.call(rbind, lapply(fits, function(f) mvrnorm(n_sims_per_fit, f$beta, f$vcov)))
n_sims <- nrow(B_full)
# pred_prop() below builds its design matrix in this column order
stopifnot(identical(colnames(B_full), c("(Intercept)", "sexmale", "LogBM_z", "age_z")))

coef_SI <- data.frame(term = colnames(B_full),
                      est = apply(B_full, 2, median),
                      lwr = apply(B_full, 2, quantile, 0.025),
                      upr = apply(B_full, 2, quantile, 0.975), row.names = NULL)
age_row <- coef_SI$term == "age_z"
coef_SI$per_year     <- ifelse(age_row, coef_SI$est / age_scale, NA)
coef_SI$per_year_lwr <- ifelse(age_row, coef_SI$lwr / age_scale, NA)
coef_SI$per_year_upr <- ifelse(age_row, coef_SI$upr / age_scale, NA)
print(coef_SI, digits = 3)
write.csv(coef_SI, file.path(out_dir, "tableS_glmm_coefficients_full.csv"), row.names = FALSE)


## ---- 4. GROWTH (MASS-AT-AGE) PARAMETER DRAWS ---- ##
## W(age) = W_inf * (1 - exp(-(k * age + x)))^3, spring parameters, Kingsley et al. (1983) Table 1.
## Parameter uncertainty uses the reported standard errors (not SE * sqrt(n)) and error correlations.
vbgm <- function(age, W_inf, k, x) W_inf * (1 - exp(-k * age - x))^3
growth_mean <- list(female = c(W_inf = 105.0, k = 0.409, x = 0.506),
                    male   = c(W_inf = 190.6, k = 0.283, x = 0.407))
growth_se   <- list(female = c(5.7, 0.049, 0.0374), male = c(11.1, 0.026, 0.0282))
growth_cor  <- list(female = c(-0.684, -0.043, -0.421), male = c(-0.702, -0.054, -0.516))  # (W,k), (W,x), (k,x)
make_vcov <- function(se, r) {
  R <- diag(3); R[1, 2] <- R[2, 1] <- r[1]; R[1, 3] <- R[3, 1] <- r[2]; R[2, 3] <- R[3, 2] <- r[3]
  diag(se) %*% R %*% diag(se)
}
G <- lapply(c(female = "female", male = "male"), function(s)
  mvrnorm(n_sims, growth_mean[[s]], make_vcov(growth_se[[s]], growth_cor[[s]])))
G_point <- lapply(growth_mean, function(p) matrix(p, n_sims, 3, byrow = TRUE))
stopifnot(all(sapply(G, function(g) all(g > 0))))    # no negative asymptotes or growth rates

mass_mat <- function(g, ages) {                       # sims x ages
  g[, 1] * (1 - exp(-outer(g[, 2], ages) - g[, 3]))^3
}


## ---- 5. INTAKE-EQUATION PARAMETER DRAWS ---- ##
## %salmon = a + b * ln(kg salmon per kg bear). Central values are the published equations used in
## the main text. Their sampling covariance comes from refitting the same form to the intake data.
## Male: the 79 "VD.et al" males reproduce the published male equation exactly, so we use that fit.
## Female: refitting the 36 females gives a = 19.0, b = 23.2 rather than the published 14.61, 26.73;
## we keep the published values (as in the main text) and borrow the refit covariance.
intake_file <- here("data", "raw", "Deacy_Robbins_intake_data.csv")
intake_published <- list(female = c(a = 14.61, b = 26.73), male = c(a = 43.26, b = 16.50))
intake_data <- read.csv(intake_file, fileEncoding = "UTF-8-BOM")[, 1:4]
names(intake_data) <- c("sex", "pct_salmon", "kg_per_kg", "source")
intake_data <- filter(intake_data, sex %in% c("female", "male"), !is.na(pct_salmon), !is.na(kg_per_kg))
fit_intake <- function(d) nls(pct_salmon ~ a + b * log(kg_per_kg), data = d, start = list(a = 20, b = 20))  # non-linear least squares
intake_fit <- list(
  female   = fit_intake(filter(intake_data, sex == "female")),
  male     = fit_intake(filter(intake_data, sex == "male", source == "VD.et al")),
  male_all = fit_intake(filter(intake_data, sex == "male"))       # all males, kept for comparison
)
print(lapply(intake_fit, coef))
stopifnot(max(abs(coef(intake_fit$male) - intake_published$male)) < 0.01)

H <- list(female = mvrnorm(n_sims, intake_published$female, vcov(intake_fit$female)),
          male   = mvrnorm(n_sims, intake_published$male,   vcov(intake_fit$male)))
H_point <- lapply(intake_published, function(p) matrix(p, n_sims, 2, byrow = TRUE))


## ---- 6. PREDICTIONS ---- ##
pred_prop <- function(B, s, ages, LogBM_z = 0) {       # sims x ages
  X <- cbind(1, as.numeric(s == "male"), LogBM_z, (ages - age_center) / age_scale)
  plogis(B %*% t(X))
}
P <- lapply(c(female = "female", male = "male"), function(s) pred_prop(B_full, s, ages_pred))
prop_SI <- bind_rows(lapply(c("female", "male"), function(s) {
  q <- apply(P[[s]][, ages_pred %in% c(3, 5, 10, 13, 15, 20)], 2, quantile, probs)
  at_mean <- quantile(plogis(B_full %*% c(1, s == "male", 0, 0)), probs)
  data.frame(sex = s, age = c(3, 5, 10, 13, 15, 20, NA), lwr = c(q[1, ], at_mean[1]),
             med = c(q[2, ], at_mean[2]), upr = c(q[3, ], at_mean[3]))
}))
print(prop_SI, digits = 3)   # age NA = at mean age and mean biomass
write.csv(prop_SI, file.path(out_dir, "tableS_diet_proportion_full.csv"), row.names = FALSE)

kg_draws <- function(s, growth, intake) {
  I <- exp((100 * P[[s]] - intake[[s]][, 1]) / intake[[s]][, 2])   # kg salmon per kg bear
  M <- mass_mat(growth[[s]], ages_pred)                              # kg bear
  list(I = I, M = M, KG = I * M)
}
stages <- list(
  "MixSIAR + GLMM + age"                 = list(growth = G_point, intake = H_point),
  "+ mass-at-age parameters"             = list(growth = G,       intake = H_point),
  "+ intake-equation parameters (full)"  = list(growth = G,       intake = H)
)
KGs <- lapply(stages, function(st) lapply(c(female = "female", male = "male"),
                                          function(s) kg_draws(s, st$growth, st$intake)))

col_of <- function(a) which(ages_pred == a)
# each row is one comparison window (rows 4-5 = full sampled-age range):
spans <- data.frame(age_start = c(3, 5, 10, 3, 3),
                    age_end   = c(13, 15, 20, 20, 25))

kg_by_age <- bind_rows(lapply(names(KGs), function(st) bind_rows(lapply(c("female", "male"), function(s) {
  q <- apply(KGs[[st]][[s]]$KG, 2, quantile, probs)
  data.frame(stage = st, sex = s, age = ages_pred, lwr = q[1, ], med = q[2, ], upr = q[3, ])
}))))
write.csv(kg_by_age, file.path(out_dir, "tableS_kg_by_age_full.csv"), row.names = FALSE)

fold_SI <- bind_rows(lapply(names(KGs), function(st) bind_rows(lapply(c("female", "male"), function(s) {
  bind_rows(lapply(seq_len(nrow(spans)), function(i) {
    K  <- KGs[[st]][[s]]$KG
    k0 <- K[, col_of(spans$age_start[i])]; k1 <- K[, col_of(spans$age_end[i])]
    data.frame(stage = st, sex = s, span = paste0(spans$age_start[i], "-", spans$age_end[i]),
               kg_start = t(quantile(k0, probs)), kg_end = t(quantile(k1, probs)),
               fold = t(quantile(k1 / k0, probs)))
  }))
}))))
print(fold_SI, digits = 3)


## ---- 7. COMPARE WITH THE MAIN-TEXT RESULT (script 06) ---- ##
# Main-text intervals (MixSIAR + GLMM, fixed age scenarios)
main <- readRDS(file.path(out_dir, "main_text_draws.rds"))   # written by 06 (see its saveRDS block)
stopifnot(nrow(main$KG) == nrow(main$kg_grid))
spans_main <- spans[1:3, ]
fold_main <- bind_rows(lapply(c("female", "male"), function(s) bind_rows(lapply(seq_len(nrow(spans_main)), function(i) {
  g  <- main$kg_grid
  i0 <- which(g$sex == s & g$bio_age_start == spans_main$age_start[i] & g$relative_age == 1)
  i1 <- which(g$sex == s & g$bio_age_start == spans_main$age_start[i] & g$relative_age == 11)
  stopifnot(length(i0) == 1, length(i1) == 1)
  data.frame(stage = "Main text: MixSIAR + GLMM (fixed age scenario)", sex = s,
             span = paste0(spans_main$age_start[i], "-", spans_main$age_end[i]),
             kg_start = t(quantile(main$KG[i0, ], probs)), kg_end = t(quantile(main$KG[i1, ], probs)),
             fold = t(quantile(main$KG[i1, ] / main$KG[i0, ], probs)))
}))))
fold_all <- bind_rows(fold_main, fold_SI) %>%
  mutate(years = as.numeric(sub(".*-", "", span)) - as.numeric(sub("-.*", "", span)),
         pct_per_year = (fold.50.^(1 / years) - 1) * 100)

print(fold_all)
write.csv(fold_all, file.path(out_dir, "tableS_fold_change_main_vs_SI.csv"), row.names = FALSE)


## ---- 8. FIGURES ---- ##
# (Within/between analysis and its figures are in 06; this script is the full-propagation sensitivity analysis only.)

# Fig. S20: realized ages
figS20a <- data.frame(first_year = factor(rep(ind_first$first_year, n_real)),
                      age_first = as.vector(age_first_mat)) %>%
  count(first_year, age_first) %>% group_by(first_year) %>% mutate(p = n / sum(n)) %>%
  ggplot(aes(age_first, p)) + geom_col(fill = "grey45", width = 0.8) +
  facet_wrap(~ first_year, nrow = 2) + scale_x_continuous(breaks = c(3, 9, 15)) +
  labs(tag = "a", x = "Realized age at first detection (years)", y = "Proportion of draws") + theme_ms
figS20b <- data.frame(age = as.vector(age_mat), sex = rep(df$sex, n_real)) %>%
  count(sex, age) %>% group_by(sex) %>% mutate(p = n / sum(n)) %>%
  ggplot(aes(age, p, fill = sex)) + geom_col(position = position_dodge(0.8), width = 0.8) +
  geom_vline(xintercept = 20.5, linetype = "22", colour = "grey50") +
  scale_fill_manual(values = sex_cols, labels = c("Female", "Male"), name = NULL) +
  labs(tag = "b", x = "Realized age of sampled bear-years (years)", y = "Proportion of draws") +
  theme_ms + theme(legend.position = "inside", legend.position.inside = c(0.85, 0.8))
figS20 <- figS20a / figS20b + plot_layout(heights = c(1.3, 1))
print(figS20)
ggsave(file.path(out_dir, "FigS20_age_realizations.png"), figS20, width = 8, height = 7, dpi = 300)
ggsave(file.path(out_dir, "FigS20_age_realizations.pdf"), figS20, width = 8, height = 7)

# Fig. S21: diet proportion vs biomass, sex and absolute age (MixSIAR + age + GLMM)
# = main-text Fig. 2 recreated with uncertainty from all sources
LogBM_seq <- seq(min(df$LogBM), max(df$LogBM), length.out = 60)
bm_dat <- bind_rows(lapply(c("female", "male"), function(s) {
  Pm <- sapply(LogBM_seq, function(b) plogis(B_full %*% c(1, s == "male", (b - bm_center) / bm_scale, 0)))
  q <- apply(Pm, 2, quantile, probs)
  data.frame(sex = s, LogBM = LogBM_seq, lwr = q[1, ], med = q[2, ], upr = q[3, ])
}))
sex_dat <- bind_rows(lapply(c("female", "male"), function(s) {
  q <- quantile(plogis(B_full %*% c(1, s == "male", 0, 0)), probs)
  data.frame(sex = s, lwr = q[1], med = q[2], upr = q[3])
}))
age_dat <- bind_rows(lapply(c("female", "male"), function(s) {
  q <- apply(P[[s]], 2, quantile, probs)
  data.frame(sex = s, age = ages_pred, lwr = q[1, ], med = q[2, ], upr = q[3, ])
}))
figS21a <- ggplot(bm_dat, aes(LogBM, med, colour = sex, fill = sex)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) + geom_line(linewidth = 0.8) +
  scale_colour_manual(values = sex_cols, labels = c("Female", "Male"), name = NULL) +
  scale_fill_manual(values = sex_cols, labels = c("Female", "Male"), name = NULL) +
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(tag = "a", x = expression(log~salmon~biomass~density~(kg~km^-2)),
       y = "Proportion of salmon in annual diet") +
  theme_ms + theme(legend.position = "inside", legend.position.inside = c(0.8, 0.12),
                   legend.background = element_blank())
figS21b <- ggplot(sex_dat, aes(sex, med, colour = sex)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.12, linewidth = 0.7) + geom_point(size = 2.8) +
  scale_colour_manual(values = sex_cols, guide = "none") + scale_x_discrete(labels = c("Female", "Male")) +
  coord_cartesian(ylim = c(0.4, 1)) + labs(tag = "b", x = NULL, y = NULL) + theme_ms
figS21c <- ggplot(age_dat, aes(age, med, colour = sex, fill = sex)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) + geom_line(linewidth = 0.8) +
  scale_colour_manual(values = sex_cols, guide = "none") + scale_fill_manual(values = sex_cols, guide = "none") +
  coord_cartesian(ylim = c(0.4, 1)) + labs(tag = "c", x = "Age (years)", y = NULL) + theme_ms
figS21 <- figS21a + figS21b + figS21c + plot_layout(widths = c(1.2, 0.6, 1.2))
print(figS21)
ggsave(file.path(out_dir, "FigS21_m1_full_propagation.png"), figS21, width = 10, height = 3.8, dpi = 300)

# Fig. S30. Panel a: kg salmon per bear by age; the grey bands show how much wider the 95% interval gets
# when growth-curve and intake-equation uncertainty are added to MixSIAR + GLMM + age.
# Panel b: fold change in each window, by source of uncertainty, next to the main-text result.
figS30a_dat <- kg_by_age %>%
  filter(stage %in% c("MixSIAR + GLMM + age", "+ intake-equation parameters (full)")) %>%
  mutate(sex_lab = ifelse(sex == "female", "Female", "Male"),
         stage = factor(stage, levels = c("+ intake-equation parameters (full)", "MixSIAR + GLMM + age"),
                        labels = c("All sources", "MixSIAR + GLMM + age only")))
figS30a <- ggplot(figS30a_dat, aes(age, med)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = stage), colour = NA) +
  geom_line(data = filter(figS30a_dat, stage == "All sources"), aes(colour = sex), linewidth = 0.9) +
  facet_wrap(~ sex_lab, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("All sources" = "grey88", "MixSIAR + GLMM + age only" = "grey68"), name = "95% interval") +
  scale_colour_manual(values = sex_cols, guide = "none") +
  scale_x_continuous(breaks = seq(3, 25, 2)) + scale_y_continuous(labels = scales::comma) +
  labs(tag = "a", x = "Age (years)", y = expression(Salmon~consumption~(kg~bear^-1~yr^-1))) +
  theme_ms + theme(legend.position = "top")
stage_levels <- c("Main text: MixSIAR + GLMM (fixed age scenario)", names(stages))
figS30b_dat <- fold_all %>%
  mutate(stage = factor(stage, levels = rev(stage_levels)),
         sex_lab = ifelse(sex == "female", "Female", "Male"),
         span = factor(paste("Age", span), levels = paste("Age", c("3-13", "5-15", "10-20", "3-20", "3-25"))))
figS30b <- ggplot(figS30b_dat, aes(fold.50., stage, colour = sex_lab)) +
  geom_errorbar(aes(xmin = fold.2.5., xmax = fold.97.5.), orientation = "y", width = 0.2,
                linewidth = 0.6, position = position_dodge(0.5)) +
  geom_point(size = 2, position = position_dodge(0.5)) +
  facet_wrap(~ span, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = setNames(sex_cols, c("Female", "Male")), name = NULL) +
  labs(tag = "b", x = "Fold change in salmon consumed over the age window", y = NULL) +
  theme_ms + theme(legend.position = "top")
figS30 <- free(figS30a) / figS30b + plot_layout(heights = c(1, 1.1))
print(figS30)
ggsave(file.path(out_dir, "FigS30_kg_full_propagation.png"), figS30, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_dir, "FigS30_kg_full_propagation.pdf"), figS30, width = 10, height = 8)


## ---- 9. SAVE SUMMARY AND PACKAGE VERSIONS ---- ##
saveRDS(list(coef_SI = coef_SI, fold_all = fold_all,
             age_center = age_center, age_scale = age_scale, n_fits = length(fits),
             intake_fit = lapply(intake_fit, function(f) list(coef = coef(f), vcov = vcov(f)))),
        file.path(out_dir, "SI_summary.rds"))

sessionInfo()
