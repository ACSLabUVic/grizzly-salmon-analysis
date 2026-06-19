#################
#
#  Here we estimate the amount of salmon consumed by grizzly bears by age.
#  Use predicted salmon consumption from m1 and sex-specific bear mass estimates
#  
#################
## sampling_age  = years since first detection (1–10, used in GLMM)
## biological_age = assumed true age = relative_age + anchor_offset

## ---- Libraries ----
library(glmmTMB)
library(dplyr)
library(tidyr)
library(ggplot2)


## ---- KINGSLEY MASS FUNCTION & DATA FRAME ---- ##

# 1.  Define von Bertalanffy mass function
# USING GROWTH CURVE PARAMETER ESTIMATES BESED ON FALL VALUES FOR MALE AND FEMALE
# USING VALUES FROM TABLE 1 in KINGSLEY 1983
## must convert mass from tonnes to kg (mtlpy by 10000)

#W_inf = asymptotic proportional mass

vbgm_king83 <- function(age, W_inf, k, x) {
  W_inf * (1 - exp(-k * age - x))^3
}

# 2.  Set model parameters of Spring mass from Kingsley et al. (1983) from table 2
kingsley83_param <- list(
  female   = list(W_inf = 105.0, k = 0.409, x = 0.506),
  male = list(W_inf = 190.6, k = 0.283, x = 0.407)
)

## Checks:
# female at age 20 should approach ~105 kg (near asymptote)
vbgm_king83(20, kingsley83_param$female$W_inf,
            kingsley83_param$female$k,
            kingsley83_param$female$x)
# female at age 3 should be well below asympote (~35-50 kg)
vbgm_king83(3, kingsley83_param$female$W_inf,
            kingsley83_param$female$k,
            kingsley83_param$female$x)
# Male at age 20 should approach ~190 kg
vbgm_king83(20, kingsley83_param$male$W_inf,
            kingsley83_param$male$k,
            kingsley83_param$male$x)

# Compute mass for biological ages 1–20
ages_bio <- 1:20

kingsleymass_df <- data.frame(
  biological_age = ages_bio,
  Female_mass_kg = vbgm_king83(ages_bio,
                               kingsley83_param$female$W_inf,
                               kingsley83_param$female$k,
                               kingsley83_param$female$x),
  Male_mass_kg   = vbgm_king83(ages_bio,
                               kingsley83_param$male$W_inf,
                               kingsley83_param$male$k,
                               kingsley83_param$male$x)
)

print(kingsleymass_df)

# Reshape to long format for joining
kingsley_long <- kingsleymass_df %>%
  pivot_longer(cols = -biological_age,
               names_to  = "sex",
               values_to = "kingsley_mass_kg") %>%
  mutate(sex = ifelse(grepl("Female", sex), "female", "male"))

## ---- LOAD GB DATA & REFIT MODEL ---- ##
df <- read.csv(here("data", "processed", "all_cov_scaled_no_locations.csv"))


m1 <- glmmTMB(Salmon.med ~ sex + LogBM_z + relative_age_z + 
                (1 | WSHDFID ) + (1 | year) + (1 | individual), 
              data = df, family = beta_family(link = "logit"))

# Scaling constants from the ORIGINAL data (must match what model was fitted on)
age_center <- mean(df$relative_age, na.rm = TRUE)
age_scale  <- sd(df$relative_age,   na.rm = TRUE)
## -------------------------------------

## ---- CONVERSION FUNCTIONS ---- ##
# inverted from: %Salmon = 26.73 * ln(X) + 14.61  
convert_female <- function(pct_salmon) exp((pct_salmon - 14.61) / 26.73)

# inverted from: %Salmon = 16.50 * ln(X) + 43.26  
convert_male   <- function(pct_salmon) exp((pct_salmon - 43.26) / 16.50)

## create anchor windows ---
# Three scenarios for unknown age at first detection
# offset = bio_age_start - 1
# so that at relative_age = 1, biological_age = bio_age_start
anchor_windows <- list(
  list(label = "Window 1: Bio age 3-13",  bio_age_start = 3,  offset = 2),
  list(label = "Window 2: Bio age 5-15",  bio_age_start = 5,  offset = 4),
  list(label = "Window 3: Bio age 10-20", bio_age_start = 10, offset = 9)
)


##  BUILD PREDICTION GRIDS FOR EACH ANCHOR WINDOW ##
build_window_grid <- function(window, kingsley_long, age_center, age_scale) {
  
  expand.grid(
    sex          = c("female", "male"),
    relative_age = 1:11,         # observed sampling range
    LogBM_z      = 0             # hold salmon biomass at mean
  ) %>%
    mutate(
      window_label   = window$label,
      bio_age_start  = window$bio_age_start,
      biological_age = relative_age + window$offset,  
      relative_age_z = (relative_age - age_center) / age_scale
    ) %>%
    # Join Kingsley mass by BIOLOGICAL age (not sampling age)
    left_join(kingsley_long, by = c("sex" = "sex", "biological_age" = "biological_age"))
}

# Build grids for all windows and stack into one data frame
all_grids <- bind_rows(lapply(anchor_windows, build_window_grid,
                              kingsley_long = kingsley_long,
                              age_center    = age_center,
                              age_scale     = age_scale))


## --- PREDICT SALMON PROPORTION & CONVERT TO KG ---

all_grids$predict_salm.med <- predict(
  m1,
  newdata = all_grids,
  type    = "response",
  re.form = NA
)

all_grids <- all_grids %>%
  mutate(
    salmon_percent        = predict_salm.med * 100,
    kg_salmon_per_kg_bear = case_when(
      sex == "female" ~ convert_female(salmon_percent),
      sex == "male"   ~ convert_male(salmon_percent)
    ),
    kg_salmon_per_bear = kg_salmon_per_kg_bear * kingsley_mass_kg
  )

## SALMON RATIO FUNCTION (START vs END OF EACH AGE WINDOW)

salmon_ratio_by_window <- function(grid = all_grids) {
  grid %>%
    group_by(window_label, sex) %>%
    summarise(
      bio_age_start   = min(biological_age),
      bio_age_end     = max(biological_age),
      kg_start        = kg_salmon_per_bear[which.min(biological_age)],
      kg_end          = kg_salmon_per_bear[which.max(biological_age)],
      diff_kg         = kg_end - kg_start,
      ratio           = kg_end / kg_start,
      .groups         = "drop"
    )
}

ratio_table <- salmon_ratio_by_window()
print(ratio_table)
# A tibble: 6 × 8
#window_label            sex    bio_age_start bio_age_end kg_start kg_end diff_kg ratio
#<chr>                   <chr>          <dbl>       <dbl>    <dbl>  <dbl>   <dbl> <dbl>
#1 Window 1: Bio age 3-13  female             3          13     372.  1462.   1090.  3.93
#2 Window 1: Bio age 3-13  male               3          13     477.  3209.   2732.  6.73
#3 Window 2: Bio age 5-15  female             5          15     523.  1469.    947.  2.81
#4 Window 2: Bio age 5-15  male               5          15     768.  3280.   2512.  4.27
#5 Window 3: Bio age 10-20 female            10          20     647.  1475.    828.  2.28
#6 Window 3: Bio age 10-20 male              10          20    1156.  3353.   2196.  2.90

## create a summary table of ages, masses, and estimated kg salmon
summary_table <- all_grids %>%
  dplyr::select(
    window_label,
    sex,
    relative_age,
    biological_age,
    kingsley_mass_kg,
    predict_salm.med,
    salmon_percent,
    kg_salmon_per_kg_bear,
    kg_salmon_per_bear
  ) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 2))
  ) %>%
  arrange(window_label, sex, biological_age)


## ---- PLOT ---- ##
all_grids <- all_grids %>%
  mutate(anchor_label = case_when(
    bio_age_start == 3  ~ "Age 3-13",
    bio_age_start == 5  ~ "Age 5-15",
    bio_age_start == 10 ~ "Age 10-20"
  )) %>%
  mutate(anchor_label = factor(anchor_label,
                               levels = c("Age 3-13", "Age 5-15", "Age 10-20")))

ggplot(filter(all_grids, relative_age <= 10),
       aes(x = relative_age, y = kg_salmon_per_bear, colour = anchor_label)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ sex, ncol = 2,
             labeller = labeller(sex = c("female" = "Female (F)",
                                         "male"   = "Male (M)"))) +
  scale_colour_manual(
    values = c("Age 3-13"  = "#66c2a5",
               "Age 5-15"  = "#fc8d62",
               "Age 10-20" = "#7570b3")
  ) +
  scale_x_continuous(
    breaks    = seq(1, 10, by = 1),
    expand    = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x      = "Relative age (years since first detection)",
    y      = "Predicted salmon consumption (kg/bear)",
    colour = "Biological Age Scenario"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text       = element_text(face = "bold", size = 13),
    legend.position  = "bottom"
  )


print(kingsleymass_df)


## --- GETTING VALUES FOR PAPER --- ##
# Predicted proportions at biological age 3 and 13
all_grids %>%
  filter(bio_age_start == 3,
         biological_age %in% c(3, 13)) %>%
  dplyr::select(sex, biological_age, predict_salm.med) %>%
  arrange(sex, biological_age)

#sex biological_age predict_salm.med
#1 female              3        0.6375994
#2 female             13        0.8561782
#3   male              3        0.7494754
#4   male             13        0.9100921


## -----------------------------------------------------------------------------
## June 12th test
library(ggplot2)
library(dplyr)
library(patchwork)

# ---- Define equation curves ---- #
x_seq <- seq(0.01, 17, by = 0.1)  # kg salmon per kg bear

# Female equation curve (Deacy et al. 2018)
female_curve <- data.frame(
  x          = x_seq,
  pct_salmon = 26.73 * log(x_seq) + 14.61,
  sex        = "Female"
)

# Male equation curve (Robbins unpublished)
male_curve <- data.frame(
  x          = x_seq,
  pct_salmon = 16.50 * log(x_seq) + 43.26,
  sex        = "Male"
)

curves_df <- rbind(female_curve, male_curve) %>%
  filter(pct_salmon >= 0, pct_salmon <= 100)

# ---- Your predicted values from agegrid (Window 1: bio age 3-13) ---- #
# These are the predicted % salmon and kg/kg bear from your model
predicted_points <- all_grids %>%
  filter(bio_age_start == 3) %>%
  dplyr::select(sex, biological_age, 
                predict_salm.med,
                kg_salmon_per_kg_bear,
                kg_salmon_per_bear) %>%
  mutate(
    sex        = ifelse(sex == "female", "Female", "Male"),
    pct_salmon = predict_salm.med * 100
  )

# ---- Published male empirical values from literature ---- #
# From your spreadsheet — adult males plotted at approximate kg/kg bear
# Using mean body mass ~180kg for adult males to convert kg/year to kg/kg bear
male_lit <- data.frame(
  study   = c(
    "Van Daele et al. 2013 (SW Kodiak)",
    "Van Daele et al. 2013 (NE Kodiak)",
    "Van Daele et al. 2013 (E Kodiak)",
    "Van Daele et al. 2013 (N Islands)",
    "Van Daele et al. 2013 (NW Kodiak)"
  ),
  mean_kg_year = c(3429, 2182, 2273, 2370, 2834),
  sd_kg_year   = c(2234, 1799, 1390, 1470, 1923),
  body_mass_kg = 180  # approximate adult male mass
) %>%
  mutate(
    kg_per_kg_bear     = mean_kg_year / body_mass_kg,
    kg_per_kg_bear_low = pmax((mean_kg_year - sd_kg_year) / body_mass_kg, 0),
    kg_per_kg_bear_hi  = (mean_kg_year + sd_kg_year) / body_mass_kg,
    # Back-calculate % salmon using male equation
    pct_salmon = 16.50 * log(kg_per_kg_bear) + 43.26
  )

# ================================================================
# PANEL 1 — Female equation with curve + your predicted values
# ================================================================
p_female <- ggplot() +
  
  # Equation curve
  geom_line(
    data = female_curve,
    aes(x = x, y = pct_salmon),
    colour    = "#FF83FA",
    linewidth = 1.2
  ) +
  
  # Your predicted values from model
  geom_point(
    data = predicted_points %>% filter(sex == "Female"),
    aes(x = kg_salmon_per_kg_bear, y = pct_salmon,
        colour = factor(biological_age)),
    size  = 2.5,
    shape = 16
  ) +
  
  # Annotations
  annotate("text", x = 12, y = 15,
           label = "Female\n(Deacy et al. 2018)",
           colour = "#FF83FA", size = 3.5, fontface = "italic") +
  
  scale_colour_viridis_d(name = "Biological age", option = "plasma") +
  scale_x_continuous(limits = c(0, 17),
                     breaks = seq(0, 17, by = 2)) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, by = 20)) +
  labs(
    tag    = "a",
    x      = "Salmon intake (kg salmon kg bear⁻¹)",
    y      = "Assimilated diet (% salmon)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    axis.text       = element_text(size = 11),
    axis.title      = element_text(size = 12),
    plot.tag        = element_text(face = "bold", size = 14)
  )

# ================================================================
# PANEL 2 — Male equation curve + your predicted values
#           + published literature points
# ================================================================
p_male <- ggplot() +
  
  # Equation curve
  geom_line(
    data = male_curve,
    aes(x = x, y = pct_salmon),
    colour    = "#7A67EE",
    linewidth = 1.2
  ) +
  
  # Your predicted values from model
  geom_point(
    data = predicted_points %>% filter(sex == "Male"),
    aes(x = kg_salmon_per_kg_bear, y = pct_salmon,
        colour = factor(biological_age)),
    size  = 2.5,
    shape = 16
  ) +
  
  # Published literature points with error bars
  geom_errorbarh(
    data = male_lit,
    aes(y        = pct_salmon,
        xmin     = kg_per_kg_bear_low,
        xmax     = kg_per_kg_bear_hi),
    colour    = "grey40",
    height    = 1.5,
    linewidth = 0.6
  ) +
  geom_point(
    data  = male_lit,
    aes(x = kg_per_kg_bear, y = pct_salmon,
        shape = study),
    colour = "grey30",
    size   = 3
  ) +
  
  # Annotations
  annotate("text", x = 12, y = 15,
           label = "Male\n(Robbins unpublished)",
           colour = "#7A67EE", size = 3.5, fontface = "italic") +
  annotate("text", x = 12, y = 25,
           label = "● Published estimates\n(Van Daele et al. 2013)",
           colour = "grey30", size = 3, hjust = 0) +
  
  scale_colour_viridis_d(name = "Biological age", option = "plasma") +
  scale_shape_manual(
    values = c(15, 16, 17, 18, 8),
    guide  = "none"  # too many to show clearly
  ) +
  scale_x_continuous(limits = c(0, 35),
                     breaks = seq(0, 35, by = 5)) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, by = 20)) +
  labs(
    tag    = "b",
    x      = "Salmon intake (kg salmon kg bear⁻¹)",
    y      = NULL  # shared with panel a
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",  # shared legend from panel a
    axis.text       = element_text(size = 11),
    axis.title      = element_text(size = 12),
    plot.tag        = element_text(face = "bold", size = 14)
  )

# ================================================================
# COMBINE PANELS
# ================================================================
validation_plot <- p_female + p_male +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

validation_plot

# ---- Save ---- #
ggsave("SI_equation_validation.png",
       plot  = validation_plot,
       width = 12, height = 5,
       dpi   = 300)

### ATTEMPT 2:##################################################################
library(ggplot2)
library(dplyr)
library(patchwork)

# ---- Published empirical values from literature ---- #
# ---- Literature data ---- #
lit_female <- data.frame(
  study   = c("Hilderbrand et al. 1999", "Deacy et al. 2018",
              "Gende & Quinn 2004",
              "Van Daele et al. 2013\n(SW Kodiak)",
              "Van Daele et al. 2013\n(NE Kodiak)",
              "Van Daele et al. 2013\n(E Kodiak)",
              "Van Daele et al. 2013\n(N Islands)",
              "Van Daele et al. 2013\n(NW Kodiak)"),
  mean_kg = c(1800, 1099, 2160, 1939, 980, 1139, 985, 1143),
  sd_kg   = c(NA, 560, NA, 1416, 561, 1212, 462, 1285),
  sex     = "Female",
  bio_age = 10
)

lit_male <- data.frame(
  study   = c("Hilderbrand et al. 1999",
              "Van Daele et al. 2013\n(SW Kodiak)",
              "Van Daele et al. 2013\n(NE Kodiak)",
              "Van Daele et al. 2013\n(E Kodiak)",
              "Van Daele et al. 2013\n(N Islands)",
              "Van Daele et al. 2013\n(NW Kodiak)"),
  mean_kg = c(950, 3429, 2182, 2273, 2370, 2834),
  sd_kg   = c(NA, 2234, 1799, 1390, 1470, 1923),
  sex     = "Male",
  bio_age = 12
)

lit_values <- rbind(lit_female, lit_male) %>%
  mutate(
    ci_low  = pmax(mean_kg - sd_kg, 0, na.rm = TRUE),
    ci_high = mean_kg + sd_kg
  )

# ---- Get predicted consumption from all_grids (Window 1: bio age 3-13) ---- #
pred_data <- all_grids %>%
  filter(bio_age_start == 3) %>%
  dplyr::select(sex, biological_age, kg_salmon_per_bear) %>%
  mutate(sex = ifelse(sex == "female", "Female", "Male"))

# ---- Approximate biological age for literature points ---- #
# Adult females ~ age 10, adult males ~ age 12
# Adjust these if you have better information
lit_values <- lit_values %>%
  mutate(bio_age = ifelse(sex == "Female", 10, 12))

# ---- Colour scale ---- #
sex_colours <- c("Female" = "#FF83FA", "Male" = "#7A67EE")

# ================================================================
# PANEL A — FEMALES
# ================================================================
p_female <- ggplot() +
  
  # Predicted consumption line
  geom_line(
    data = pred_data %>% filter(sex == "Female"),
    aes(x = biological_age, y = kg_salmon_per_bear),
    colour    = "#FF83FA",
    linewidth = 1.2
  ) +
  
  # Shaded region showing published range
  annotate(
    "rect",
    xmin  = 3, xmax = 13,
    ymin  = 250,    # min of published female range
    ymax  = 2224,   # max of published female range
    fill  = "grey80",
    alpha = 0.3
  ) +
  annotate(
    "text", x = 3.2, y = 2300,
    label  = "Published adult female range\n(Van Daele et al. 2013)",
    colour = "grey40", size = 3, hjust = 0
  ) +
  
  # Literature error bars
  geom_errorbar(
    data = lit_values %>%
      filter(sex == "Female", !is.na(sd_kg)),
    aes(x    = bio_age,
        ymin = ci_low,
        ymax = ci_high),
    colour    = "grey40",
    width     = 0.3,
    linewidth = 0.7
  ) +
  
  # Literature points
  geom_point(
    data = lit_values %>% filter(sex == "Female"),
    aes(x = bio_age, y = mean_kg, shape = study),
    colour = "grey30",
    size   = 3
  ) +
  
  scale_shape_manual(
    values = c(15, 16, 17, 18, 8, 11, 12, 13),
    name   = "Literature source"
  ) +
  scale_x_continuous(
    limits = c(3, 13),
    breaks = seq(3, 13, by = 1)
  ) +
  scale_y_continuous(
    limits = c(0, 2500),
    breaks = seq(0, 2500, by = 500),
    labels = scales::comma
  ) +
  labs(
    tag   = "a",
    title = "Female",
    x     = "Biological age (years)",
    y     = "Predicted salmon consumption (kg year⁻¹)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold", size = 13),
    plot.tag        = element_text(face = "bold", size = 14),
    axis.text       = element_text(size = 11),
    axis.title      = element_text(size = 12)
  )

# ================================================================
# PANEL B — MALES
# ================================================================
p_male <- ggplot() +
  
  # Predicted consumption line
  geom_line(
    data = pred_data %>% filter(sex == "Male"),
    aes(x = biological_age, y = kg_salmon_per_bear),
    colour    = "#7A67EE",
    linewidth = 1.2
  ) +
  
  # Shaded region showing published range
  annotate(
    "rect",
    xmin  = 3, xmax = 13,
    ymin  = 854,    # min of published male range
    ymax  = 4381,   # max of published male range
    fill  = "grey80",
    alpha = 0.3
  ) +
  annotate(
    "text", x = 3.2, y = 4500,
    label  = "Published adult male range\n(Van Daele et al. 2013)",
    colour = "grey40", size = 3, hjust = 0
  ) +
  
  # Literature error bars
  geom_errorbar(
    data = lit_values %>%
      filter(sex == "Male", !is.na(sd_kg)),
    aes(x    = bio_age,
        ymin = ci_low,
        ymax = ci_high),
    colour    = "grey40",
    width     = 0.3,
    linewidth = 0.7
  ) +
  
  # Literature points
  geom_point(
    data = lit_values %>% filter(sex == "Male"),
    aes(x = bio_age, y = mean_kg, shape = study),
    colour = "grey30",
    size   = 3
  ) +
  
  scale_shape_manual(
    values = c(15, 16, 17, 18, 8, 11),
    name   = "Literature source"
  ) +
  scale_x_continuous(
    limits = c(3, 13),
    breaks = seq(3, 13, by = 1)
  ) +
  scale_y_continuous(
    limits = c(0, 5000),
    breaks = seq(0, 5000, by = 500),
    labels = scales::comma
  ) +
  labs(
    tag   = "b",
    title = "Male",
    x     = "Biological age (years)",
    y     = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 13),
    plot.tag        = element_text(face = "bold", size = 14),
    axis.text       = element_text(size = 11),
    axis.title      = element_text(size = 12),
    legend.text     = element_text(size = 9),
    legend.title    = element_text(size = 10, face = "bold")
  )

# ================================================================
# COMBINE AND SAVE
# ================================================================
validation_plot <- p_female + p_male +
  plot_layout(widths = c(1, 1.4))  # extra width for legend

validation_plot

ggsave("SI_validation_figure.png",
       plot  = validation_plot,
       width = 12, height = 5,
       dpi   = 300)

# Alternative ways to view pred_data
pred_data %>%
  filter(sex == "Female") %>%
  as.data.frame() %>%
  print()

# Check your prediction data
pred_data %>%
  group_by(sex) %>%
  summarise(
    min_age = min(biological_age),
    max_age = max(biological_age),
    min_kg  = min(kg_salmon_per_bear),
    max_kg  = max(kg_salmon_per_bear),
    n       = n()
  )
