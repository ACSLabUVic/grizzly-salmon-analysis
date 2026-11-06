####
# Propagates stable-isotope mixing-model (MixSIAR) uncertainty into the beta GLMM
#     by refitting the top model to whole posterior iterations of the mixing model
#     (multiple-imputation style), and adds GLMM estimation uncertainty by simulating
#     fixed effects from each refit. All draws are pooled; intervals are the 2.5th and
#     97.5th percentiles of the pooled draws.
#  2. Converts predicted diet proportions to kg salmon per bear per year under the
#     three age-at-first-detection scenarios (3, 5, 10 yr), carrying every draw through
#     so fold changes and differences get their own intervals.
#  3. Decomposes the age-related increase in kg salmon into the part due to body
#     growth and the part due to a higher proportion of salmon in the diet
#     (symmetric / Shapley decomposition, so there is no unallocated interaction).
#
#  Mass-at-age (Kingsley et al. 1983) and the intake equations are held at their point
#  estimates here; their uncertainty, and uncertainty in age at first detection,
#  are propagated in the SI script 07_SI_full_error_propagation.R.
#
#  Updated: October 28, 2026
####

## ---- Libraries ----
library(here)
library(glmmTMB)
library(MASS)       # mvrnorm
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(MuMIn)

select <- dplyr::select   # MASS masks dplyr::select

## ---- Settings ----
set.seed(20261027)
n_draws         <- as.integer(Sys.getenv("N_DRAWS", 1000))  # MixSIAR posterior iterations used (of 3000); one GLMM refit each
n_sims_per_draw <- 10     # fixed-effect draws per refit -> n_draws * n_sims_per_draw pooled draws
n_cores         <- max(1, detectCores() - 2)
probs           <- c(0.025, 0.5, 0.975)

out_dir <- here("output")
dir.create(out_dir, showWarnings = FALSE)


comp_cols   <- c("Body growth" = "#B07AA1", "Diet composition" = "#59A14F")


## ---- DATA ---- ##
# Covariates (550 bear-years); Salmon.med is the MixSIAR posterior median
df <- read.csv(here("data", "processed", "covariotes_no_locations_2026-11-04.csv"))
df$sex <- factor(df$sex, levels = c("female", "male"))

# Full MixSIAR posteriors for salmon proportion (638 bear-years x 3000 iterations) 
post <- data.table::fread(here("data", "processed", "gb_all_posteriors_no_location_2026-10-20.csv"), data.table = FALSE)
draw_cols <- grep("^iteration_", names(post))

# Align posterior rows to the 550 bear-years used in the GLMM
idx <- match(paste(df$individual, df$year), paste(post$individual, post$year))
stopifnot(!anyNA(idx))
Y <- as.matrix(post[idx, draw_cols])   # 550 x 3000, rows aligned with df
rm(post)

# Checks: draws reproduce the posterior medians used in the original analysis,
# and every draw is strictly inside (0, 1) as the beta likelihood requires
stopifnot(max(abs(apply(Y, 1, median) - df$Salmon.med)) < 0.02)
stopifnot(all(Y > 0 & Y < 1))

df$Salmon.2.5  <- apply(Y, 1, quantile, 0.025)
df$Salmon.97.5 <- apply(Y, 1, quantile, 0.975)

# Columns of Y are whole MCMC iterations. Using a whole column for all bears keeps the
# posterior correlation among bears (shared source and discrimination uncertainty).
draw_ids <- round(seq(1, ncol(Y), length.out = n_draws))

# Scaling constants used to build the *_z predictors
age_center <- mean(df$relative_age); age_scale <- sd(df$relative_age)
bm_center  <- mean(df$LogBM);        bm_scale  <- sd(df$LogBM)
stopifnot(max(abs((df$relative_age - age_center) / age_scale - df$relative_age_z)) < 1e-6)
stopifnot(max(abs((df$LogBM - bm_center) / bm_scale - df$LogBM_z)) < 1e-6)

# Observed range of relative age by sex (females are not observed beyond relative age 7)
max_rel_age <- tapply(df$relative_age, df$sex, max)
print(max_rel_age)


## ---- GLMM: point-estimate fit (original analysis, Salmon.med) ---- ##
glmm_formula <- ~ sex + LogBM_z + relative_age_z +
  (1 | WSHDFID) + (1 | year) + (1 | individual)

fit_glmm <- function(y) {
  d <- df
  d$y <- y
  m <- glmmTMB(update(glmm_formula, y ~ .), data = d,
               family = beta_family(link = "logit"),
               control = glmmTMBControl(parallel = 1))
  list(beta = fixef(m)$cond, vcov = vcov(m)$cond,
       ok = isTRUE(m$sdr$pdHess) && all(is.finite(vcov(m)$cond)))
}

m1 <- glmmTMB(update(glmm_formula, Salmon.med ~ .), data = df,
              family = beta_family(link = "logit"))
summary(m1)


## ---- GLMM: refit to MixSIAR posterior iterations (cached) ---- ##
fits_file <- file.path(out_dir, paste0("fits_mixsiar_draws_", n_draws, ".rds"))
if (file.exists(fits_file)) {
  fits <- readRDS(fits_file)
  stopifnot(length(fits) == n_draws)
} else {
  fits <- mclapply(draw_ids, function(j) fit_glmm(Y[, j]), mc.cores = n_cores)
  saveRDS(fits, fits_file)
}
ok <- vapply(fits, function(f) is.list(f) && isTRUE(f$ok), logical(1))
cat("Refits converged (positive-definite Hessian):", sum(ok), "of", length(ok), "\n")
fits <- fits[ok]

# Pooled fixed-effect draws: MixSIAR uncertainty (between refits) +
# GLMM estimation uncertainty (within refit)
B_prop <- do.call(rbind, lapply(fits, function(f)
  mvrnorm(n_sims_per_draw, f$beta, f$vcov)))
# For comparison: GLMM uncertainty only (posterior medians treated as known)
B_med <- mvrnorm(nrow(B_prop), fixef(m1)$cond, vcov(m1)$cond)


## ---- Coefficient table ---- ##
beta_hat <- t(sapply(fits, `[[`, "beta"))
W_within <- colMeans(t(sapply(fits, function(f) diag(f$vcov))))
B_betw   <- apply(beta_hat, 2, var)
m_fit    <- length(fits)
T_total  <- W_within + (1 + 1 / m_fit) * B_betw

coef_table <- data.frame(
  term            = colnames(B_prop),
  est_median_only = fixef(m1)$cond,
  lwr_median_only = confint(m1, parm = "beta_")[, 1],
  upr_median_only = confint(m1, parm = "beta_")[, 2],
  est_propagated  = apply(B_prop, 2, median),
  lwr_propagated  = apply(B_prop, 2, quantile, 0.025),
  upr_propagated  = apply(B_prop, 2, quantile, 0.975),
  # share of the total variance that comes from mixing-model uncertainty
  frac_var_mixsiar = (1 + 1 / m_fit) * B_betw / T_total,
  row.names = NULL
)
# Per-year slope of relative age on the logit scale (comparable with the SI)
coef_table$per_year_note <- ifelse(coef_table$term == "relative_age_z",
                                   paste0("per year: ", round(coef_table$est_propagated / age_scale, 3)), "")
print(coef_table, digits = 3)

coef_table$width_ratio <- (coef_table$upr_propagated - coef_table$lwr_propagated) /
  (coef_table$upr_median_only - coef_table$lwr_median_only)
coef_table$pct_wider   <- 100 * (coef_table$width_ratio - 1)
print(coef_table[, c("term", "pct_wider", "frac_var_mixsiar")], digits = 3)
write.csv(coef_table, file.path(out_dir, "table_glmm_coefficients_propagated_2026-11-04.csv"), row.names = FALSE)


## ---- Prediction helper ---- ##
# Returns a matrix (rows = newdata rows, columns = draws) of predicted salmon proportions
predict_draws <- function(B, newdata) {
  X <- model.matrix(~ sex + LogBM_z + relative_age_z,
                    data = transform(newdata, sex = factor(sex, levels = c("female", "male"))))
  plogis(X %*% t(B[, colnames(X)]))
}

summarise_draws <- function(P, prefix = "") {
  q <- t(apply(P, 1, quantile, probs, na.rm = TRUE))
  setNames(as.data.frame(q), paste0(prefix, c("lwr", "med", "upr")))
}


## ---- FIGURE 2: diet proportion vs biomass, sex and relative age ---- ##
sex_cols <- c(female = "#0072B2", male = "#E69F00")

# (b) biomass, at mean relative age
grid_bm <- expand.grid(sex = c("female", "male"),
                       LogBM = seq(min(df$LogBM), max(df$LogBM), length.out = 60),
                       relative_age_z = 0) %>%
  mutate(LogBM_z = (LogBM - bm_center) / bm_scale)
grid_bm <- bind_cols(grid_bm, summarise_draws(predict_draws(B_prop, grid_bm)))

# (c) sex, at mean biomass and mean relative age
grid_sex <- data.frame(sex = c("female", "male"), LogBM_z = 0, relative_age_z = 0)
grid_sex <- bind_cols(grid_sex,
                      summarise_draws(predict_draws(B_prop, grid_sex)),
                      summarise_draws(predict_draws(B_med, grid_sex), "medonly_"))

# (d) relative age, at mean biomass
grid_age <- expand.grid(sex = c("female", "male"),
                        relative_age = seq(1, 11, by = 0.25), LogBM_z = 0) %>%
  mutate(relative_age_z = (relative_age - age_center) / age_scale,
         observed = relative_age <= max_rel_age[as.character(sex)])
grid_age <- bind_cols(grid_age,
                      summarise_draws(predict_draws(B_prop, grid_age)),
                      summarise_draws(predict_draws(B_med, grid_age), "medonly_"))
#observed range: 
prop_table <- bind_rows(
  grid_sex %>% transmute(sex = as.character(sex), relative_age = NA, lwr, med, upr),
  grid_age %>% filter(relative_age %in% c(1, 7, 11)) %>%
    transmute(sex = as.character(sex), relative_age, lwr, med, upr))
print(prop_table, digits = 3)   # relative_age NA = at mean relative age and mean biomass
write.csv(prop_table, file.path(out_dir, "table_diet_proportion_propagated_2026-11-04.csv"), row.names = FALSE)

range(prop_table$med) # reported in results
# 0.6280787 0.9046010
median(prop_table$med)
#  0.778155

# Bear-year data with their MixSIAR 95% credible intervals
pts <- df %>% mutate(sex = factor(sex, levels = c("female", "male")),
                     relative_age_jit = relative_age + runif(n(), -0.3, 0.3))



theme_ms <- theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.tag = element_text(face = "plain"),
        plot.tag.location = "panel",
        plot.tag.position = c(0.065, 1.025))

p2b <- ggplot(grid_bm, aes(LogBM, med, colour = sex, fill = sex)) +
  geom_linerange(data = pts, aes(x = LogBM, ymin = Salmon.2.5, ymax = Salmon.97.5),
                 inherit.aes = FALSE, colour = "grey88", linewidth = 0.15) +
  geom_point(data = pts, aes(x = LogBM, y = Salmon.med, colour = sex),
             inherit.aes = FALSE, size = 0.7, alpha = 0.5) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = sex_cols, labels = c("Female", "Male"), name = NULL) +
  scale_fill_manual(values = sex_cols, labels = c("Female", "Male"), name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(tag = "b)", x = expression(Log~salmon~biomass~density~(kg~km^-2)),
       y = "Proportion of salmon in annual diet") +
  theme_ms + theme(legend.position = "top")

p2c <- ggplot(grid_sex, aes(sex, med, colour = sex)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.12, linewidth = 0.7) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = sex_cols, guide = "none") +
  scale_x_discrete(labels = c("Female", "Male")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(tag = "c)", x = NULL, y = NULL) +
  theme_ms

p2d <- ggplot(grid_age, aes(relative_age, med, colour = sex, fill = sex)) +
  geom_linerange(data = pts, aes(x = relative_age_jit, ymin = Salmon.2.5, ymax = Salmon.97.5),
                 inherit.aes = FALSE, colour = "grey88", linewidth = 0.15) +
  geom_point(data = pts, aes(x = relative_age_jit, y = Salmon.med, colour = sex),
             inherit.aes = FALSE, size = 0.7, alpha = 0.5) +
  geom_ribbon(data = filter(grid_age, observed), aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) +
  geom_line(data = filter(grid_age, observed), linewidth = 0.8) +
  geom_line(data = filter(grid_age, !observed), linewidth = 0.6, linetype = "22") +
  scale_colour_manual(values = sex_cols, guide = "none") +
  scale_fill_manual(values = sex_cols, guide = "none") +
  scale_x_continuous(breaks = 1:11) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(tag = "d)", x = "Relative age\n(years since first detection)", y = NULL) +
  theme_ms


## Coefficient panel built from the PROPAGATED estimates
coef_dat <- coef_table %>%
  filter(term != "(Intercept)") %>%
  mutate(term = factor(term,
                       levels = c("sexmale", "LogBM_z", "relative_age_z"),
                       labels = c("Sex (male)", "Log biomass\ndensity", "Relative age")))

p2coef <- ggplot(coef_dat, aes(x = est_propagated, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lwr_propagated, xmax = upr_propagated),
                orientation = "y", width = 0.15, linewidth = 0.7) +
  geom_point(size = 2.5) +
  labs(tag = "a)", x = "Standardized coefficient\n(logit scale, 95% CI)", y = NULL) +
  theme_ms


fig2 <- p2coef + p2b + p2c + p2d + plot_layout(widths = c(0.8, 1, 0.7, 1), guides = "collect") &
  theme(legend.position = "top")

print(fig2)
ggsave(file.path(out_dir, "Fig2_diet_propagated_2026-11-04.png"), fig2, width = 10, height = 3.8, dpi = 300)
ggsave(file.path(out_dir, "Fig2_diet_propagated_2026-11-04.pdf"), fig2, width = 10, height = 3.8)

## BIOENERGETIC CALCULATIONS
## ---- KINGSLEY MASS FUNCTION ---- ##
# W(age) = W_inf * (1 - exp(-(k * age + x)))^3 ; spring parameters, Kingsley et al. (1983) Table 1
vbgm_king83 <- function(age, W_inf, k, x) W_inf * (1 - exp(-k * age - x))^3

kingsley83_param <- list(
  female = list(W_inf = 105.0, k = 0.409, x = 0.506),
  male   = list(W_inf = 190.6, k = 0.283, x = 0.407)
)
mass_at_age <- function(age, sex) {
  p <- kingsley83_param[[sex]]
  vbgm_king83(age, p$W_inf, p$k, p$x)
}

## ---- CONVERSION FUNCTIONS (% salmon -> kg salmon per kg bear) ---- ##
# female from Deacy et al. (2018); male from C. Robbins (unpublished)
# inverted from: %Salmon = a + b * ln(X)
intake_param <- list(female = c(a = 14.61, b = 26.73), male = c(a = 43.26, b = 16.50))
convert_intake <- function(pct_salmon, sex) {
  p <- intake_param[[sex]]
  exp((pct_salmon - p["a"]) / p["b"])
}


## ---- KG SALMON PER BEAR, ALL DRAWS ---- ##
anchor_windows <- data.frame(
  window_label  = c("Age 3-13", "Age 5-15", "Age 10-20"),
  bio_age_start = c(3, 5, 10)
)

kg_grid <- expand.grid(sex = c("female", "male"), relative_age = 1:11,
                       window_label = anchor_windows$window_label,
                       stringsAsFactors = FALSE) %>%
  left_join(anchor_windows, by = "window_label") %>%
  mutate(biological_age = bio_age_start + relative_age - 1,
         relative_age_z = (relative_age - age_center) / age_scale,
         LogBM_z = 0,
         mass_kg = mapply(mass_at_age, biological_age, sex))

P_kg   <- predict_draws(B_prop, kg_grid)                 # proportion, rows = kg_grid, cols = draws
I_kg   <- P_kg                                           # kg salmon per kg bear
for (s in c("female", "male")) {
  r <- kg_grid$sex == s
  I_kg[r, ] <- convert_intake(100 * P_kg[r, ], s)
}
KG <- I_kg * kg_grid$mass_kg                              # kg salmon per bear per year

kg_summary <- bind_cols(kg_grid,
                        summarise_draws(P_kg, "prop_"),
                        summarise_draws(I_kg, "intake_"),
                        summarise_draws(KG, "kg_")) %>%
  mutate(observed = relative_age <= max_rel_age[sex])

saveRDS(list(B_prop = B_prop, B_med = B_med, kg_grid = kg_grid, KG = KG,
             age_center = age_center, age_scale = age_scale),
        file.path(out_dir, "main_text_draws.rds"))

write.csv(kg_summary, file.path(out_dir, "table_kg_by_age_propagated_2026-11-04.csv"), row.names = FALSE)


## ---- FOLD CHANGE AND DIFFERENCE OVER EACH 10-YEAR WINDOW ---- ##
# Computed draw by draw, so start and end share the same draw (their correlation is kept)
row_of <- function(win, s, ra) which(kg_grid$window_label == win & kg_grid$sex == s & kg_grid$relative_age == ra)

window_contrast <- function(ra_end) {
  bind_rows(lapply(anchor_windows$window_label, function(win) {
    bind_rows(lapply(c("female", "male"), function(s) {
      i0 <- row_of(win, s, 1); i1 <- row_of(win, s, ra_end)
      k0 <- KG[i0, ]; k1 <- KG[i1, ]
      q <- function(v) quantile(v, probs)
      data.frame(window_label = win, sex = s,
                 bio_age_start = kg_grid$biological_age[i0], bio_age_end = kg_grid$biological_age[i1],
                 extrapolated = ra_end > max_rel_age[s],
                 prop_start = t(q(P_kg[i0, ])), prop_end = t(q(P_kg[i1, ])),
                 kg_start = t(q(k0)), kg_end = t(q(k1)),
                 diff_kg = t(q(k1 - k0)), fold = t(q(k1 / k0)),
                 N_start_kg = t(q(k0 * 0.0281)), P_start_kg = t(q(k0 * 0.0035)),
                 N_end_kg   = t(q(k1 * 0.0281)), P_end_kg   = t(q(k1 * 0.0035)))
    }))
  }))
}
ratio_table <- window_contrast(11)
# Females are observed only to relative age 7: also report the within-data contrast
ratio_table_female_obs <- window_contrast(7) %>% filter(sex == "female")
print(ratio_table, digits = 3)
print(ratio_table_female_obs, digits = 3)
write.csv(ratio_table, file.path(out_dir, "table_fold_change_propagated_2026-11-04.csv"), row.names = FALSE)
write.csv(ratio_table_female_obs, file.path(out_dir, "table_fold_change_females_observed_range.csv"), row.names = FALSE)

# Same contrasts with GLMM-only uncertainty, to show what propagating MixSIAR adds
P_med  <- predict_draws(B_med, kg_grid)
KG_med <- P_med
for (s in c("female", "male")) {
  r <- kg_grid$sex == s
  KG_med[r, ] <- convert_intake(100 * P_med[r, ], s) * kg_grid$mass_kg[r]
}
width_compare <- bind_rows(lapply(anchor_windows$window_label, function(win) {
  bind_rows(lapply(c("female", "male"), function(s) {
    i0 <- row_of(win, s, 1); i1 <- row_of(win, s, 11)
    data.frame(window_label = win, sex = s,
               fold_medonly = t(quantile(KG_med[i1, ] / KG_med[i0, ], probs)),
               fold_propagated = t(quantile(KG[i1, ] / KG[i0, ], probs)),
               kg_end_medonly = t(quantile(KG_med[i1, ], probs)),
               kg_end_propagated = t(quantile(KG[i1, ], probs)))
  }))
}))
print(width_compare, digits = 3)
write.csv(width_compare, file.path(out_dir, "table_interval_width_medonly_vs_propagated.csv"), row.names = FALSE)


## ---- FIGURE 3: kg salmon per bear by relative age and scenario ---- ##
window_cols <- c("Age 3-13" = "#1baf7a", "Age 5-15" = "#eda100", "Age 10-20" = "#6250d6")
fig3_dat <- kg_summary %>%
  mutate(window_label = factor(window_label, levels = anchor_windows$window_label),
         sex_lab = ifelse(sex == "female", "Female", "Male"))

fig3 <- ggplot(fig3_dat, aes(relative_age, kg_med, colour = window_label, fill = window_label)) +
  geom_ribbon(data = filter(fig3_dat, observed), aes(ymin = kg_lwr, ymax = kg_upr),
              alpha = 0.15, colour = NA) +
  geom_line(data = filter(fig3_dat, observed), linewidth = 0.9) +
  geom_line(data = filter(fig3_dat, relative_age >= max_rel_age[sex]),
            linewidth = 0.6, linetype = "22") +
  geom_ribbon(data = filter(fig3_dat, relative_age >= max_rel_age[sex]),
              aes(ymin = kg_lwr, ymax = kg_upr), alpha = 0.05, colour = NA) +
  facet_wrap(~ sex_lab, ncol = 2) +
  scale_colour_manual(values = window_cols, name = "Biological age\nscenario") +
  scale_fill_manual(values = window_cols, name = "Biological age\nscenario") +
  scale_x_continuous(breaks = 1:11) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Relative age (years since first detection)",
       y = expression(Predicted~salmon~consumption~(kg~bear^-1~yr^-1))) +
  theme_ms + theme(legend.position = "right")
print(fig3)
ggsave(file.path(out_dir, "Fig3_kg_by_age_propagated_2026-11-04.png"), fig3, width = 8.5, height = 3.8, dpi = 300)
ggsave(file.path(out_dir, "Fig3_kg_by_age_propagated_2026-11-04.pdf"), fig3, width = 8.5, height = 3.8)


## ---- DECOMPOSITION: body growth vs diet composition ---- ##
## kg = M (body mass) x I (kg salmon per kg bear). Between the start of a window (M0, I0)
## and a later age (M1, I1) the increase is split symmetrically (Shapley / midpoint):
##   growth part = (M1 - M0) * (I0 + I1) / 2
##   diet part   = (I1 - I0) * (M0 + M1) / 2
## The two parts sum exactly to M1*I1 - M0*I0, so there is no leftover interaction term;
## the growth x diet interaction is shared equally. The multiplicative version,
## log(T1/T0) = log(M1/M0) + log(I1/I0), is reported alongside.

decomp_draws <- function(win, s, ra) {
  i0 <- row_of(win, s, 1); i1 <- row_of(win, s, ra)
  M0 <- kg_grid$mass_kg[i0]; M1 <- kg_grid$mass_kg[i1]
  I0 <- I_kg[i0, ];          I1 <- I_kg[i1, ]
  growth <- (M1 - M0) * (I0 + I1) / 2
  diet   <- (I1 - I0) * (M0 + M1) / 2
  list(growth = growth, diet = diet, total = growth + diet,
       share_growth = growth / (growth + diet),
       log_share_growth = log(M1 / M0) / log((M1 * I1) / (M0 * I0)),
       fold_growth = rep(M1 / M0, length(I0)), fold_diet = I1 / I0)
}

decomp_path <- bind_rows(lapply(anchor_windows$window_label, function(win) {
  bind_rows(lapply(c("female", "male"), function(s) {
    bind_rows(lapply(1:11, function(ra) {   # ra = 1 is the reference age (both parts 0)
      d <- decomp_draws(win, s, ra)
      i1 <- row_of(win, s, ra)
      data.frame(window_label = win, sex = s, relative_age = ra,
                 biological_age = kg_grid$biological_age[i1],
                 component = c("Body growth", "Diet composition"),
                 lwr = c(quantile(d$growth, 0.025), quantile(d$diet, 0.025)),
                 med = c(median(d$growth), median(d$diet)),
                 upr = c(quantile(d$growth, 0.975), quantile(d$diet, 0.975)))
    }))
  }))
}))

decomp_table <- bind_rows(lapply(anchor_windows$window_label, function(win) {
  bind_rows(lapply(c("female", "male"), function(s) {
    d <- decomp_draws(win, s, 11)
    q <- function(v) quantile(v, probs)
    data.frame(window_label = win, sex = s, extrapolated = 11 > max_rel_age[s],
               increase_kg = t(q(d$total)),
               growth_kg = t(q(d$growth)), diet_kg = t(q(d$diet)),
               pct_growth = t(q(100 * d$share_growth)),
               pct_diet = t(q(100 * (1 - d$share_growth))),
               pct_growth_log = t(q(100 * d$log_share_growth)),
               fold_growth = d$fold_growth[1], fold_diet = t(q(d$fold_diet)))
  }))
}))
print(decomp_table, digits = 3)
write.csv(decomp_table, file.path(out_dir, "table_decomposition_propagated_2026-11-04.csv"), row.names = FALSE)


## ---- FIGURE 4: decomposition ---- ##
fig4a_dat <- decomp_path %>%
  filter(window_label == "Age 3-13") %>%
  mutate(sex_lab = ifelse(sex == "female", "Female", "Male"),
         observed = relative_age <= max_rel_age[sex])

fig4a <- ggplot(fig4a_dat, aes(biological_age, med, colour = component, fill = component)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.18, colour = NA) +
  geom_line(data = filter(fig4a_dat, observed), linewidth = 0.9) +
  geom_line(data = filter(fig4a_dat, relative_age >= max_rel_age[sex]), linewidth = 0.6, linetype = "22") +
  facet_wrap(~ sex_lab, ncol = 2) +
  scale_colour_manual(values = comp_cols, name = NULL) +
  scale_fill_manual(values = comp_cols, name = NULL) +
  scale_x_continuous(breaks = seq(3, 13, 2)) +
  scale_y_continuous(labels = scales::comma) +
  labs(tag = "a", x = "Biological age (years; first detected at age 3)",
       y = "Increase in salmon consumed since age 3\n(kg per bear per year)") +
  theme_ms + theme(legend.position = "top")

fig4b_dat <- decomp_table %>%
  mutate(window_label = factor(window_label, levels = anchor_windows$window_label),
         sex_lab = ifelse(sex == "female", "Female", "Male"))

fig4b <- ggplot(fig4b_dat, aes(window_label, pct_growth.50., colour = sex_lab)) +
  geom_hline(yintercept = 50, colour = "grey88") +
  geom_errorbar(aes(ymin = pct_growth.2.5., ymax = pct_growth.97.5.),
                width = 0.15, linewidth = 0.7, position = position_dodge(0.4)) +
  geom_point(size = 2.5, position = position_dodge(0.4)) +
  scale_colour_manual(values = setNames(sex_cols, c("Female", "Male")), name = NULL) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(tag = "b", x = "Age-at-first-detection scenario",
       y = "Share of 10-year increase\ndue to body growth (%)") +
  theme_ms + theme(legend.position = "top")

fig4 <- fig4a + fig4b + plot_layout(widths = c(2, 1))
print(fig4)
ggsave(file.path(out_dir, "Fig4_decomposition_propagated_2026-11-04.png"), fig4, width = 10, height = 4, dpi = 300)
ggsave(file.path(out_dir, "Fig4_decomposition_propagated_2026-11-04.pdf"), fig4, width = 10, height = 4)

## --- WITHIN- vs BETWEEN-BEAR effect of relative age ----
# splitting relative age into within and between parts 
# For each bear, between = that bear's average relative age (how long it was followed overall).
# within = how far each sample's relative age is from that bear's own average.
df_wb <- df %>%
  group_by(individual) %>%                        
  mutate(n_years     = n_distinct(year),          
         age_between = mean(relative_age),        
         age_within  = relative_age - age_between) %>%  
  ungroup() %>%                                  
  mutate(age_between_z = as.numeric(scale(age_between)),  
         age_within_z  = as.numeric(scale(age_within))) 


#  fit to whole dataset
m_wb_all <- glmmTMB(Salmon.med ~ sex + LogBM_z + age_between_z + age_within_z +
                      (1 | WSHDFID) + (1 | year) + (1 | individual),
                    data = df_wb, family = beta_family(link = "logit"))

#  Sensitivity check: only bears seen in more than one year
m_wb_multi <- update(m_wb_all, data = filter(df_wb, n_years > 1))

m1_chk <- update(m1, data = df_wb)

summary(m_wb_all)
summary(m_wb_multi)
r.squaredGLMM(m_wb_all)

# How many bears and observations carry the within-bear estimate?
n_multi <- df_wb %>% filter(n_years > 1) %>%
  summarise(bears = n_distinct(individual), obs = n())
print(n_multi)

# m1 is the special case where within and between slopes are equal, so the models are nested.
lrt_wb <- anova(m1_chk, m_wb_all)
print(lrt_wb)
p_wb <- lrt_wb$`Pr(>Chisq)`[2]                   # p-value from the second row

## ---slopes per year (model-only) 
sd_w <- sd(df_wb$age_within)
sd_b <- sd(df_wb$age_between)

get_slope <- function(model, term, sd_term) {
  est <- fixef(model)$cond[term]
  se  <- sqrt(vcov(model)$cond[term, term])
  c(est = est / sd_term,
    lwr = (est - 1.96 * se) / sd_term,
    upr = (est + 1.96 * se) / sd_term)
}

wb_table <- rbind(
  within_all      = get_slope(m_wb_all,   "age_within_z",  sd_w),
  between_all     = get_slope(m_wb_all,   "age_between_z", sd_b),
  within_multiyr  = get_slope(m_wb_multi, "age_within_z",  sd(filter(df_wb, n_years > 1)$age_within)),
  between_multiyr = get_slope(m_wb_multi, "age_between_z", sd(filter(df_wb, n_years > 1)$age_between)))
wb_table <- data.frame(term = rownames(wb_table), wb_table, row.names = NULL)
wb_table$LRT_p_slopes_differ <- p_wb              # same p-value repeated on each row for reference
print(wb_table, digits = 3)
write.csv(wb_table, file.path(out_dir, "table_age_within_between_2025-11-04.csv"), row.names = FALSE)

## ---- Refit the within/between model to the MixSIAR iterations (builds B_wb) ---- ##
stopifnot(identical(df_wb$individual, df$individual), identical(df_wb$year, df$year))
fit_glmm_wb <- function(y) {
  d <- df_wb
  d$y <- y                                     # one MixSIAR iteration as the response
  m <- glmmTMB(y ~ sex + LogBM_z + age_between_z + age_within_z +
                 (1 | WSHDFID) + (1 | year) + (1 | individual),
               data = d, family = beta_family(link = "logit"),
               control = glmmTMBControl(parallel = 1))
  list(beta = fixef(m)$cond, vcov = vcov(m)$cond,
       ok = isTRUE(m$sdr$pdHess) && all(is.finite(vcov(m)$cond)))
}

fits_wb_file <- file.path(out_dir, paste0("fits_mixsiar_draws_wb_", n_draws, ".rds"))
if (file.exists(fits_wb_file)) {
  fits_wb <- readRDS(fits_wb_file)
} else {
  fits_wb <- mclapply(draw_ids, function(j) fit_glmm_wb(Y[, j]), mc.cores = n_cores)
  saveRDS(fits_wb, fits_wb_file)
}
ok_wb <- vapply(fits_wb, function(f) is.list(f) && isTRUE(f$ok), logical(1))
cat("Within/between refits converged:", sum(ok_wb), "of", length(ok_wb), "\n")
fits_wb <- fits_wb[ok_wb]
stopifnot(length(fits_wb) == length(fits))     # same draws in both, so B_wb pairs with B_prop

set.seed(20261027)
B_wb <- do.call(rbind, lapply(fits_wb, function(f) mvrnorm(n_sims_per_draw, f$beta, f$vcov)))
colnames(B_wb)                                 # should include age_within_z and age_between_z

# Pooled draws were made above; slope_draws puts the pooled, within and between slopes on the same per-year scale
stopifnot(nrow(B_wb) == nrow(B_prop))      # same draws, so rows pair up
slope_draws <- cbind(pooled  = B_prop[, "relative_age_z"] / age_scale,
                     within  = B_wb[, "age_within_z"]  / sd_w,
                     between = B_wb[, "age_between_z"] / sd_b)
wb_prop <- data.frame(
  term = c("relative age (pooled, per year)", "within bear (per year)", "between bears (per year)"),
  est  = apply(slope_draws, 2, median),
  lwr  = apply(slope_draws, 2, quantile, 0.025),
  upr  = apply(slope_draws, 2, quantile, 0.975),
  row.names = NULL)
print(wb_prop, digits = 3)                 # propagated values reported in the text (Table S4)
write.csv(wb_prop, file.path(out_dir, "table_age_within_between_propagated_2026-11-04.csv"), row.names = FALSE)
mean(slope_draws[, "within"] >= slope_draws[, "between"])   # share of draws with within >= between

## ///////
# Re-create the draws here 
# slope_draws: columns "within" and "between", one row per draw, in years.
slope_draws <- cbind(within  = B_wb[, "age_within_z"]  / sd_w,
                     between = B_wb[, "age_between_z"] / sd_b)
sapply(as.data.frame(slope_draws), quantile, probs = c(0.025, 0.5, 0.975))
mean(slope_draws[, "within"] >= slope_draws[, "between"])   # share of draws where within >= between
##////////

## fold-change if only the within-bear slope is used (Table SX)
start_grid <- data.frame(sex = c("female", "male"), LogBM_z = 0,
                         relative_age_z = (1 - age_center) / age_scale)
L0 <- qlogis(predict_draws(B_prop, start_grid))   # rows = sex, cols = draws; qlogis() = proportion -> logit

# Within-bear slope, one value per draw. Use the same draws for both,
# so each draw's starting diet and slope come from the same model refit.
n_d <- min(ncol(L0), nrow(slope_draws))
b_w <- slope_draws[1:n_d, "within"]

ra <- 1:11                                         # relative ages 1 to 11 (a 10-year window)

# For each sex: diet -> intake -> kg, then fold-change from relative age 1 to 11
fold_within <- sapply(c("female", "male"), function(s) {
  r <- which(c("female", "male") == s)             # row of L0 for this sex
  logit_P <- outer(ra - 1, b_w) +                  # years since first detection x slope, every draw
    matrix(L0[r, 1:n_d], length(ra), n_d, byrow = TRUE)
  P <- plogis(logit_P)                             # logit -> proportion (rows = ages, cols = draws)
  I <- exp((100 * P - intake_param[[s]]["a"]) / intake_param[[s]]["b"])  # kg salmon per kg bear
  M <- mass_at_age(3 + ra - 1, s)                  # body mass at biological ages 3 to 13
  K <- I * M                                       # kg salmon per bear per year
  quantile(K[length(ra), ] / K[1, ], probs)        # age 13 vs age 3: 2.5%, 50%, 97.5%
})
print(round(fold_within, 2))
write.csv(fold_within, file.path(out_dir, "table_fold_change_within_only_2026-11-04.csv"))

## ---- FIGURE S1: number of bear-years at each relative age, by sex ---- ##
sex_cols <- c(female = "#0072B2", male = "#E69F00")
# 1. Count bear-years for every combination of relative age and sex
n_by_age <- df %>%
  count(sex, relative_age) %>%
  complete(sex, relative_age = 1:11, fill = list(n = 0))

n_sex <- table(df$sex)
sex_labs <- c(female = paste0("Female (n = ", n_sex[["female"]], " bear-years)"),
              male   = paste0("Male (n = ",   n_sex[["male"]],   " bear-years)"))

figS1 <- ggplot(n_by_age, aes(x = relative_age, y = n, fill = sex)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(data = filter(n_by_age, n > 0), aes(label = n),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = sex_cols, labels = sex_labs, name = NULL) +
  scale_x_continuous(breaks = 1:11) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Relative age (years since first detection)",
       y = "Number of bear-years") +
  theme_ms + theme(legend.position = c(0.75, 0.85))

print(figS1)
ggsave(file.path(out_dir, "FigS1_bear_years_by_age.png"), figS1, width = 6, height = 3.6, dpi = 300)
ggsave(file.path(out_dir, "FigS1_bear_years_by_age.pdf"), figS1, width = 6, height = 3.6)

## ---- FIGURE S2: within- vs between-individual effect of relative age ---- ##
# 1. Take the propagated slopes and give each row a readable label
slopes <- wb_prop
slopes$component <- factor(c("Pooled\n(relative age, m1)",
                             "Within-\nindividual",
                             "Between-\nindividual"),
                           levels = c("Between-\nindividual",
                                      "Within-\nindividual",
                                      "Pooled\n(relative age, m1)"))

# 2. One colour per component
comp_cols_wb <- c("Pooled\n(relative age, m1)" = "grey30",
                  "Within-\nindividual"        = "#009E73",
                  "Between-\nindividual"       = "#CC79A7")

# 3. Plot the slopes with their 95% intervals
figS2 <- ggplot(slopes, aes(x = est, y = component, colour = component)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lwr, xmax = upr), orientation = "y",
                width = 0.15, linewidth = 0.7) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = comp_cols_wb, guide = "none") +
  labs(x = "Effect of relative age\n(logit scale, per year, 95% CI)", y = NULL) +
  theme_ms

print(figS2)
ggsave(file.path(out_dir, "FigS2_within_between_propagated.png"), figS2, width = 4.5, height = 3.2, dpi = 300)
ggsave(file.path(out_dir, "FigS2_within_between_propagated.pdf"), figS2, width = 4.5, height = 3.2)
