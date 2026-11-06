####################
#
# This script spatially joins grizzly bear hair snag locations to focal watersheds.
# Salmon biomass density, salmon species diversity, and mean EVI
# are joined to each bear location.Filters to individual bears with salmon consumption
# estimates (Salmon.med), selects first detection per individual per year,
# calculates relative age, and scales all predictors for modelling.
#
# Raw bear location data cannot be shared due to data sharing agreements. 
# 
####################

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)

## ----  Watershed IDs for each detection ---- ##
# we want to join covariates to merged watersheds using the focal watershed ID

# Each large (merged) watershed was built around one focal watershed that contains site
# information, and the large watersheds overlap each other. So a bear is assigned to the
# large watershed built AROUND ITS OWN focal watershed: the shapefile's FOCALWSHD is that
# centre, and FOCALWSHD -> WSHDFID is a one-to-one lookup. (We do NOT use any WSHDFID that is
# already stored in the bear file or in focal_clipped.shp: with overlapping watersheds that
# value depends on how the spatial join resolved the overlap.)
shp_df <- st_read(here("data", "processed", "LRG_Watersheds_for_all_sites", "LRG_Watersheds_all_sites_FINAL.shp")) %>%
  st_drop_geometry()

wshd_lookup <- shp_df %>%
  dplyr::select(FOCALWSHD, WSHDFID) %>%
  distinct(FOCALWSHD, .keep_all = TRUE)


gbs_in_wshds <- read.csv(here("data", "processed", "gbs_in_wshds_no_location_2026-11-04.csv"))

# focal 12820 appears in two large watersheds (WSHDFID 129 and 131), so the lookup cannot say
# which one to use. No detection WITH a diet estimate is in 12820 (checked below, after
# first_detection is built), so this only affects detections that are not used in the analysis.
#sum(gbs_in_wshds$FOCALWSHD == 12820, na.rm = TRUE)   # detections in 12820 (not used in analysis)

# focal 12845 has no large watershed built around it; assigned by hand (verified in QGIS)
wshd_lookup <- bind_rows(wshd_lookup, data.frame(FOCALWSHD = 12845, WSHDFID = 63))
stopifnot(!anyDuplicated(wshd_lookup$FOCALWSHD))


## Assign each detection to the large watershed centred on its focal watershed.
siteind_merged <- gbs_in_wshds %>%
  dplyr::select(-any_of("WSHDFID")) %>%       # drop the WSHDFID that came with the file
  left_join(wshd_lookup, by = "FOCALWSHD")
stopifnot(nrow(siteind_merged) == nrow(gbs_in_wshds))   # the join must not add or remove rows

# any individuals still without a merged watershed? (manual fixes in section 3 should cover these)
siteind_merged %>%
  filter(!is.na(individual)) %>%
  filter(is.na(WSHDFID)) %>%
  dplyr::select(individual, year, FOCALWSHD, WSHDFID)
nrow(siteind_merged)

# ----  Read the watershed covariote tables (one row per watershed and year) ---- ##

# --- Salmon biomass --- (biomass of year t-1 goes with hair collected in year t)
salm_merg <- read.csv(here("data", "processed", "merg_salm_biomass_per_year_w_density.csv")) %>%
  dplyr::select(WSHDFID, Salmon_year = year, YearlyBMinMWS, YearlyBMDensityinMWS, LogBMDensityinMWS) %>%
  filter(!is.na(Salmon_year)) %>%
  mutate(year = Salmon_year + 1)

# --- Salmon diversity (Shannon) ---
div_merge <- read.csv(here("data", "processed", "merg_spp_diversity_per_year_with_Hill.csv")) %>%
  dplyr::select(WSHDFID, SalmonDiv_year = year, Salmon_diversity_MWS = Shannon, Hill_diversity_MWS = Hill) %>%
  mutate(year = SalmonDiv_year + 1)

# --- EVI --- (output of 03)
# 03 calls its first column WTRSHD_FID, but it holds terra's 'ID', which is the ROW NUMBER of
# the polygon in the watershed shapefile (1, 2, 3 ...), NOT the watershed ID. So here we
# convert that row number to the real WSHDFID, using the same shapefile that 03 read.
wshd_order <- st_read(here("data", "processed", "LRG_Watersheds_for_all_sites", "LRG_Watersheds_all_sites_FINAL.shp")) %>%
  st_drop_geometry() %>%
  pull(WSHDFID)                      # WSHDFID in shapefile row order (row 1 first)
stopifnot(!anyDuplicated(wshd_order))

evi_raw <- read.csv(here("data", "processed", "EVI", "EVI_by_year_in_merged_WSHDs.csv"))
stopifnot(nrow(evi_raw) == length(wshd_order),                  # one row per watershed polygon
          all(evi_raw$WTRSHD_FID %in% seq_along(wshd_order)))   # row numbers, not watershed IDs

evi <- evi_raw %>%
  mutate(WSHDFID = wshd_order[WTRSHD_FID]) %>%       # row number -> real watershed ID
  dplyr::select(WSHDFID, starts_with("mean_EVI")) %>%
  pivot_longer(starts_with("mean_EVI"), names_to = "EVI_year", values_to = "EVI_Mean_MWS") %>%
  mutate(EVI_year = as.integer(substr(EVI_year, 9, 12)),
         year = EVI_year + 1)

# EVI should be on the 0-1 scale (MODIS stores it x 10,000; if you see values in the
# thousands, multiply EVI_Mean_MWS by 0.0001 here)
stopifnot("EVI is not on the 0-1 scale" = max(evi$EVI_Mean_MWS, na.rm = TRUE) <= 1.5)

# each watershed-year must appear once, otherwise the joins below would duplicate bear rows
stopifnot(!anyDuplicated(salm_merg[, c("WSHDFID", "year")]),
          !anyDuplicated(div_merge[, c("WSHDFID", "year")]),
          !anyDuplicated(evi[, c("WSHDFID", "year")]))

## ----  Function that attaches the covariates using a data frame's CURRENT WSHDFID ---- ##
# We use a function so that it can be run again after the manual watershed corrections
# in section 5 (otherwise corrected bears keep the covariates of their old watershed).

cov_cols <- c("Salmon_year", "YearlyBMinMWS", "YearlyBMDensityinMWS", "LogBMDensityinMWS",
              "SalmonDiv_year", "Salmon_diversity_MWS", "Hill_diversity_MWS", "EVI_year", "EVI_Mean_MWS")

attach_covariates <- function(df) {
  n_before <- nrow(df)
  out <- df %>%
    ungroup() %>%
    dplyr::select(-any_of(cov_cols)) %>%
    left_join(salm_merg, by = c("WSHDFID", "year")) %>%
    left_join(div_merge, by = c("WSHDFID", "year")) %>%
    left_join(evi,       by = c("WSHDFID", "year")) %>%
    # watersheds with no salmon streams have no biomass or diversity rows: assign 0
    mutate(across(c(YearlyBMinMWS, YearlyBMDensityinMWS, LogBMDensityinMWS, Salmon_diversity_MWS),
                  ~ replace_na(.x, 0)),
        Hill_diversity_MWS = replace_na(Hill_diversity_MWS, 1))   # exp(0) = 1, as for Shannon = 0
  stopifnot(nrow(out) == n_before)   # a join must never add or remove bear rows
  out
}

## ----  All detections, with covariates (saved without locations) ---- ##
siteind_EVI_salm <- siteind_merged %>%
  filter(!is.na(revisitid)) %>%
  attach_covariates()

colnames(siteind_EVI_salm)

## NOW HAVE SALMON.MED, SALMON BIOMASS, DIVERSITY, EVI AT THE LARGE WATERSHED SCALE FOR ALL SNAG SITES
write.csv(siteind_EVI_salm %>% dplyr::select(-any_of(c("X", "...1", "X.1"))),
          here("data", "processed", "sitedata_salmon_EVI_no_location_2026-11-04.csv"), row.names = FALSE)

## -----------------------------------------------------------------------------
## site data has all individuals, even those without a Salmon.med; we only need
## individuals that have a proportion of salmon

gbconsumption <- siteind_EVI_salm %>%
  filter(!is.na(Salmon.med)) %>%
  dplyr::select(individual, year, project, Salmon.med, WSHDFID, FOCALWSHD, site_id,
                revisitid, revisit_date, revisit, all_of(cov_cols))

# Remove exact duplicate rows
gb_consump_filtered <- gbconsumption %>%
  distinct() %>%
  filter(!is.na(Salmon.med))

# some individuals have multiple detections in a year; select only the first detection date
# for each individual-year combination.
# revisit_date is text (MM/DD/YYYY), so convert it to a real date before sorting
gb_consump_filtered <- gb_consump_filtered %>%
  mutate(revisit_date_parsed = as.Date(revisit_date, format = "%m/%d/%Y"))
stopifnot(!anyNA(gb_consump_filtered$revisit_date_parsed))

# Sort by date to ensure the first detection comes first
gb_consump_filtered <- gb_consump_filtered %>%
  arrange(individual, revisit_date_parsed)

# Keep only the first detection per individual per year
first_detection <- gb_consump_filtered %>%
  group_by(individual, year) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(-revisit_date_parsed)

# No analysis row may sit in focal 12820 (it belongs to two large watersheds)
stopifnot("a detection with a diet estimate is in focal 12820: choose its large watershed" =
            !(12820 %in% first_detection$FOCALWSHD))

# Check for NA's
first_detection %>%
  filter(is.na(WSHDFID)) %>%
  dplyr::select(individual, year, WSHDFID, FOCALWSHD)

## ----  Manual watershed corrections, THEN re-attach the covariates ---- ##
# Some individual-years were assigned incorrect focal watersheds in the spatial join.
# I manually checked all locations in QGIS to make sure that they are in the
# correct watershed. Here I manually fix the 19 that were incorrect.
manual_fixes <- data.frame(
  individual = c(11243, 14483, 25721, 25777, 27978, 30560, 32633, 56190,
                 139903, 139903, 149554, 197067, 197067, 197067,
                 197098, 197098, 197098, 197159, 197212),
  year       = c(2012, 2014, 2014, 2017, 2015, 2018, 2018, 2018,
                 2015, 2018, 2016, 2011, 2016, 2018,
                 2012, 2017, 2019, 2017, 2016),
  FOCALWSHD_correct = c(12380, 12373, 9625, 5979, 11814, 12358, 12372, 12728,
                        11780, 11816, 12358, 12369, 12369, 12382,
                        11814, 11815, 11814, 12375, 12380),
  WSHDFID_correct   = c(18, 21, 192, 80, 92, 106, 23, 98,
                        71, 118, 106, 22, 22, 16,
                        92, 95, 92, 33, 18)
)

# Apply corrections: fix the focal watershed, then take WSHDFID from the lookup again
first_detection <- first_detection %>%
  left_join(manual_fixes, by = c("individual", "year")) %>%
  mutate(FOCALWSHD = ifelse(!is.na(FOCALWSHD_correct), FOCALWSHD_correct, FOCALWSHD)) %>%
  dplyr::select(-WSHDFID) %>%
  left_join(wshd_lookup, by = "FOCALWSHD")

# Verify fixes applied correctly: the lookup must give the WSHDFID you chose by hand
fixed_rows <- first_detection %>%
  filter(!is.na(FOCALWSHD_correct)) %>%
  dplyr::select(individual, year, FOCALWSHD, FOCALWSHD_correct, WSHDFID, WSHDFID_correct)
fixed_rows
stopifnot("not all 19 manual fixes were found" = nrow(fixed_rows) == nrow(manual_fixes),
          all(fixed_rows$FOCALWSHD == fixed_rows$FOCALWSHD_correct),
          all(fixed_rows$WSHDFID   == fixed_rows$WSHDFID_correct))

first_detection <- first_detection %>%
  dplyr::select(-FOCALWSHD_correct, -WSHDFID_correct)

# Re-attach biomass, diversity and EVI using the CORRECTED WSHDFID
first_detection <- attach_covariates(first_detection)

# Check: everyone in the same watershed and year must have the same covariate values
# (this should print 0 rows)
first_detection %>%
  group_by(WSHDFID, year) %>%
  filter(n_distinct(YearlyBMDensityinMWS) > 1 |
           n_distinct(Salmon_diversity_MWS) > 1 |
           n_distinct(EVI_Mean_MWS) > 1) %>%
  dplyr::select(individual, year, WSHDFID, YearlyBMDensityinMWS, Salmon_diversity_MWS, EVI_Mean_MWS)

# Check: how many rows are missing EVI? (should be 0)
sum(is.na(first_detection$EVI_Mean_MWS))

## NOW WE HAVE DF WITH ONLY GBS WITH CONSUMPTION DATA

## ---- Relative age ---- ##
# relative age = number of years since the individual's first year WITH A DIET ESTIMATE
# (Salmon.med), counting that year as 1. This is the first detection with an isotope
# sample, which can be later than the first detection of the bear in the hair-snag survey.
gbcons_age <- first_detection %>%
  arrange(individual, year) %>%
  group_by(individual) %>%
  mutate(relative_age = year - min(year) + 1) %>%
  ungroup()

# Remove rows where Salmon.med is NA
gbcons_age <- gbcons_age %>%
  filter(!is.na(Salmon.med))

### Lets get an idea of age classes:
age_class_table <- gbcons_age %>%
  count(relative_age, name = "n_observations") %>%
  arrange(relative_age)

print(age_class_table)

## ---- Add sex ---- ##
# Sex is added from a lookup (individual, year, sex) created in the private script
# 04a_individual_locations_PRIVATE.R.
sex_lookup <- read.csv(here("data", "processed", "individual_year_sex_lookup.csv")) %>%
  dplyr::select(individual, year, sex) %>%
  distinct()
stopifnot(!anyDuplicated(sex_lookup[, c("individual", "year")]))

gbcons_age <- gbcons_age %>%
  left_join(sex_lookup, by = c("individual", "year"))
stopifnot(!anyNA(gbcons_age$sex))

# age classes by sex:
age_class_table_sex <- gbcons_age %>%
  count(relative_age, sex, name = "n_observations") %>%
  tidyr::pivot_wider(names_from = sex, values_from = n_observations, values_fill = 0) %>%
  arrange(relative_age)

print(age_class_table_sex)

sample_size_summary <- gbcons_age %>%
  group_by(relative_age) %>%
  summarise(
    n_individuals_total  = n_distinct(individual),
    n_female = sum(sex == "female"),
    n_male   = sum(sex == "male"),
    .groups = "drop"
  )

print(sample_size_summary)


## ---- Scaling of predictors ---- ##
# Shannon diversity is non-linear, so we also calculate Hill's number (q = 1), the
# exponential of Shannon entropy, and use it as the diversity metric in modelling.
scaled_df <- gbcons_age %>%
  rename(hill_diversity = Hill_diversity_MWS)  

scaled_df$EVI_Mean_z      <- as.numeric(scale(scaled_df$EVI_Mean_MWS))
scaled_df$BMDensity_z     <- as.numeric(scale(scaled_df$YearlyBMDensityinMWS))
scaled_df$LogBM_z         <- as.numeric(scale(scaled_df$LogBMDensityinMWS))
scaled_df$diversity_z     <- as.numeric(scale(scaled_df$Salmon_diversity_MWS))
scaled_df$hill_diversity_z <- as.numeric(scale(scaled_df$hill_diversity))
scaled_df$relative_age_z  <- as.numeric(scale(scaled_df$relative_age))

# Verify scaling worked (should be ~0 and ~1)
mean(scaled_df$relative_age_z)
sd(scaled_df$relative_age_z)


## ---- Clean data frame and save ---- ##
all_covariotes <- scaled_df %>%
  dplyr::select(year, individual, sex, Salmon.med, project, WSHDFID, FOCALWSHD, site_id, revisitid, revisit, revisit_date,
                Salmon_year, YearlyBMinMWS, YearlyBMDensityinMWS, LogBMDensityinMWS,
                SalmonDiv_year, Salmon_diversity_MWS, hill_diversity, EVI_year, EVI_Mean_MWS, relative_age,
                EVI_Mean_z, BMDensity_z, LogBM_z, diversity_z, relative_age_z,
                hill_diversity_z) %>%
  rename(
    BM        = YearlyBMinMWS,
    BMDensity = YearlyBMDensityinMWS,
    LogBM     = LogBMDensityinMWS,
    diversity = Salmon_diversity_MWS
  )

nrow(all_covariotes)   # should be 550 individual-years
n_distinct(all_covariotes$individual)   # should be 261

#### ---- save df for modelling ---- ###
write.csv(all_covariotes, here("data", "processed", "covariotes_no_locations_2026-11-04.csv"), row.names = FALSE)

# Save the package versions used
sessionInfo()

## ---------------------------------------------
# recaptured individuals (seen in more than one year)
years_per_bear <- all_covariotes %>% count(individual, name = "n_years")
sum(years_per_bear$n_years > 1)                                   # recaptured individuals
range(years_per_bear$n_years[years_per_bear$n_years > 1] - 1)     # recaptures after the first

# multiple detections within a year
# all detections, one row per detection (siteind_EVI_salm was built in section 4)
det_all <- siteind_EVI_salm %>%
  dplyr::select(individual, year, revisitid, site_id) %>%
  distinct()

# keep only detections from the 550 bear-years used in the models
det_used <- det_all %>%
  semi_join(all_covariotes, by = c("individual", "year"))

by_bearyear <- det_used %>%
  group_by(individual, year) %>%
  summarise(n_detections = n(), n_sites = n_distinct(site_id), .groups = "drop")

nrow(by_bearyear)                  # should be 550
sum(by_bearyear$n_detections > 1)  # bear-years with more than one detection
sum(by_bearyear$n_sites > 1)       # bear-years detected at more than one site

# Save the package versions used
sessionInfo()


