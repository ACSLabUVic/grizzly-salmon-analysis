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

library(readr)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)

# ==============================================================================
## This section has been commented out as it uses data with location info
## see line 114 for reproducable code

## NOT SHARED: contains raw bear locations and consumption data
#siteind1 <- read_csv("data/raw/sensitivity data/sitedata_grizzly_salmonconsumpt_to2019.csv") 
#plot(siteind1$utme, siteind1$utmn)
# ==============================================================================
## STANDARDIZE COORDINATES
## original file has utm's in zone 10 and 9, we just want UTM 9

#siteindZ10 <- subset(siteind1, siteind1$utm_zone==10)
#siteindZ9 <- subset(siteind1, siteind1$utm_zone==9)

#siteindZ9$utme9 <- siteindZ9$utme
#siteindZ9$utmn9 <- siteindZ9$utmn

## make the utm 10 sataset into a sf to get coordinates
#siteindZ10SF <- sf::st_as_sf(siteindZ10, coords = c("utme", "utmn"), crs=3157)
## convert to UTM zone 9
#siteindZ10_9SF <- sf::st_transform(siteindZ10SF, crs=6338)
## get coordinates for UTM9
#coordsZ10_9SF <- data.frame(sf::st_coordinates(siteindZ10_9SF))
#names(coordsZ10_9SF) <- c("utme9", "utmn9")
## bind back into original df
#siteindZ109 <- cbind(siteindZ10, coordsZ10_9SF)

## Combine zones and remove Gitga'at
#shared_cols <- intersect(names(siteindZ9), names(siteindZ109))
#siteind <- rbind(siteindZ9[, shared_cols], siteindZ109[, shared_cols]) %>%
  #filter(project != "gt")

## Make a spatial version of siteind
#siteind_sf <- st_as_sf(siteind, coords = c("utme9", "utmn9"), crs = 6338)


## first we are going to join individuals to the focal watersheds
#focal <- st_read(here("data", "processed", "Focal Watersheds in Study Area", "focal_clipped.shp")) %>%
#  st_transform(crs = 6338)  # match siteind_sf CRS


#st_crs(focal)
#st_crs(siteind_sf)

## visual check
#plot(focal$geometry)
#plot(siteind_sf$geometry, pch = 20, col = "red", add = TRUE)

## spatial join to focal watersheds
#siteind_focal <- st_join(siteind_sf, focal, left = TRUE)


## there are some snag locations that fall just outside of watersheds (e.g., in the water due to gps error)
## to deal with these snag sites, we will assign them to the nearest watershed
#siteind_focal_na <- siteind_focal %>% filter(is.na(WTRSHD_FID))
#siteind_focal_ok <- siteind_focal %>% filter(!is.na(WTRSHD_FID))

#if(nrow(siteind_focal_na) > 0){
  
  ## Find nearest focal watershed polygon for each unmatched point
#  nearest_idx <- st_nearest_feature(siteind_focal_na, focal)
#  siteind_focal_na$WTRSHD_FID <- focal$WTRSHD_FID[nearest_idx]
  
  ## Apply manual corrections from QGIS investigation:
  ## sites 201204 and 201453 were reassigned based on spatial configuration of
  ## nearby points and polygons
#  siteind_focal_na <- siteind_focal_na %>%
#    mutate(WTRSHD_FID = case_when(
#      site_id == "201204" ~ 9624,
#      site_id == "201453" ~ 12426,
#      TRUE ~ WTRSHD_FID
#    ))
  
  ## Rejoin with matched points
#  siteind_focal <- rbind(siteind_focal_ok, siteind_focal_na)
#}

## Drop geometry and extract coordinates
#siteind_focal_joined <- siteind_focal %>%
#  mutate(utme9 = st_coordinates(.)[,1],
#         utmn9 = st_coordinates(.)[,2]) %>%
#  st_drop_geometry()

#siteind_focal_joined <- siteind_focal_joined%>%
#  rename(FOCALWSHD = WTRSHD_FID)

## NOW WE HAVE A DATAFRAME THAT HAS HAIR SNAG LOCATIONS AND THEIR ASSOCIATED FOCAL WATERSHED ID

## Save file without location info for sharing:
#sitedata_nolocation <- siteind_focal_joined %>%
#  select(-...1, -X, -utme, -utmn, -utm_zone, -utme9, -utmn9, -AREA_HA, -F_CODE, -OBJECTID, -AREA_SQM, -FEAT_LEN, -usblAREA))

#write.csv(sitedata_nolocation, here("data", "processed", "gbs_in_wshds_no_location.csv"))
# ==============================================================================

## ---- Now join our covariates:----- ##
# we want to join covariates to merged watersheds using focal watershed ID

wshd <- st_read(here("data", "processed", "LRG_Watersheds_for_all_sites", "LRG_Watersheds_all_sites_FINAL.shp")) %>%
  st_drop_geometry() %>%
  dplyr::select(FOCALWSHD, WSHDFID) %>%
  distinct(FOCALWSHD, .keep_all = TRUE)

gbs_in_wshds <- read.csv(here("data", "processed", "gbs_in_wshds_no_location.csv"))


## Fix WSHDFID for focal watersheds that fall outside merged watershed boundaries
# Verified these in QGIS. here we assigned to nearest correct merged watershed
siteind_merged <- gbs_in_wshds %>%
  mutate(WSHDFID = case_when(
    FOCALWSHD == 12415 ~ 213L,
    FOCALWSHD == 12414 ~ 212L,
    FOCALWSHD == 12426 ~ 66L,
    FOCALWSHD == 11783 ~ 200L,
    FOCALWSHD == 11780 ~ 71L,
    FOCALWSHD == 12845 ~ 63L,
    FOCALWSHD == 9624  ~ 191L,  
    FOCALWSHD == 11833 ~ 44L,
    FOCALWSHD == 11825 ~ 42L,
    FOCALWSHD == 12403 ~ 58L,
    TRUE ~ WSHDFID
  ))

siteind_merged %>%
  filter(!is.na(individual)) %>%
  filter(is.na(WSHDFID)) %>%
  dplyr::select(individual, year, FOCALWSHD, WSHDFID)
nrow(siteind_merged)  


# --- Salmon biomass ---
salm_merg <- read.csv(here("data", "processed", "merg_salm_biomass_per_year_w_density.csv")) %>%
  dplyr::select(-FOCALWSHD)
names(salm_merg)[names(salm_merg) == "year"] <- "Salmon_year"
salm_merg$year <- salm_merg$Salmon_year+1  

siteind_salmon <- merge(siteind_merged, salm_merg, by = c("WSHDFID", "year"), all = TRUE)

siteind_salmon <- siteind_salmon[which(is.na(siteind_salmon$revisitid)== FALSE),]

# Assign 0 for watersheds that have no salmon streams (NA)
siteind_salmon <- siteind_salmon %>%
  mutate(
    YearlyBMinMWS        = ifelse(is.na(YearlyBMinMWS), 0, YearlyBMinMWS),
    YearlyBMDensityinMWS = ifelse(is.na(YearlyBMDensityinMWS), 0, YearlyBMDensityinMWS),
    LogBMDensityinMWS    = ifelse(is.na(LogBMDensityinMWS), 0, LogBMDensityinMWS)
  )

# --- diversity ---
div_merge <- read.csv(here("data", "processed", "merg_spp_diversity_per_year.csv")) %>%
  dplyr::select(WSHDFID, year, Shannon) %>%
  rename(Salmon_diversity_MWS = Shannon) 
names(div_merge)[names(div_merge) == "year"] <- "SalmonDiv_year"
div_merge$year <- div_merge$SalmonDiv_year+1

siteind_salmon_div <- merge(siteind_salmon, div_merge, by = c("WSHDFID", "year"), all = TRUE)

siteind_salmon_div <- siteind_salmon_div[which(is.na(siteind_salmon_div$revisitid)== FALSE),]

#assign 0 for the watersheds that have no salmon streams (NA)
siteind_salmon_div <- siteind_salmon_div %>% 
  mutate(Salmon_diversity_MWS = ifelse(is.na(Salmon_diversity_MWS), 0, Salmon_diversity_MWS))

# --- EVI ---
evi <- read.csv(here("data", "processed", "EVI", "EVI_by_year_in_merged_WSHDs.csv")) %>%
  dplyr::select(WTRSHD_FID, starts_with("mean_EVI")) %>%
  pivot_longer(starts_with("mean_EVI"), names_to = "EVI_year", values_to = "EVI_Mean_MWS") %>%
  mutate(EVI_year = as.integer(substr(EVI_year, 9, 12))) %>%
  rename(WSHDFID = WTRSHD_FID)

names(evi)[names(evi) == "year"] <- "EVI_year"
evi$year <- evi$EVI_year+1

siteind_EVI_salm <- merge(siteind_salmon_div, evi, by = c("WSHDFID", "year"), all = TRUE)

siteind_EVI_salm <- siteind_EVI_salm[which(is.na(siteind_EVI_salm$revisitid)== FALSE),]
colnames(siteind_EVI_salm)

siteind_EVI_salm <- siteind_EVI_salm[,c(1, 2, 4:15, 19:26)]

## NOW HAVE SALMON.MED, SALMON, DIVERSITY, EVI AT THE LARGE WATERSHED SCALE FOR ALL SNAG SITES 

write.csv(siteind_EVI_salm, here("data", "processed", "sitedata_salmon_EVI_no_location.csv"))

## -----------------------------------------------------------------------------

## #site data has all individuals, even those without a salm.med. we only need
# individuals that have a proportion of salmon 

gbconsumption <- siteind_EVI_salm %>%
  filter(!is.na(Salmon.med)) %>%
  dplyr::select(individual, year, project, Salmon.med, WSHDFID, FOCALWSHD, site_id,
                revisitid, revisit_date, revisit, Salmon_year, YearlyBMinMWS, YearlyBMDensityinMWS, LogBMDensityinMWS,
                SalmonDiv_year, Salmon_diversity_MWS, EVI_year, EVI_Mean_MWS)

# Remove exact duplicate rows
gb_consump_filtered <- gbconsumption %>%
  distinct() %>%
  filter(!is.na(Salmon.med))  

# some individuals have multiple detections in a year; select only the first detection date
# for each individual-year combination

# Sort by date to ensure first detection comes first
gb_consump_filtered <- gb_consump_filtered %>%
  arrange(individual, revisit_date)

# Keep only the first detection per individual per year
first_detection <- gb_consump_filtered %>%
  group_by(individual, year) %>%
  dplyr::slice(1)

# Check for NA's 
first_detection %>%
  filter(is.na(WSHDFID)) %>%
  dplyr::select(individual, year, WSHDFID, FOCALWSHD)

# Some individual-years were assigned incorrect focal watersheds in the spatial join
# I manually checked all locations in QGIS to make sure that they are all in the 
# correct watershed. Here I manually fix the 19 that were incorrect
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

# Apply corrections
first_detection <- first_detection %>%
  left_join(manual_fixes, by = c("individual", "year")) %>%
  mutate(
    FOCALWSHD = ifelse(!is.na(FOCALWSHD_correct), FOCALWSHD_correct, FOCALWSHD),
    WSHDFID   = ifelse(!is.na(WSHDFID_correct),   WSHDFID_correct,   WSHDFID)
  ) %>%
  dplyr::select(-FOCALWSHD_correct, -WSHDFID_correct)

# Verify fixes applied correctly
first_detection %>%
  inner_join(manual_fixes, by = c("individual", "year")) %>%
  dplyr::select(individual, year, FOCALWSHD, FOCALWSHD_correct, WSHDFID, WSHDFID_correct)

## NOW WE HAVE DF WITH ONLY GBS WITH CONSUMPTION DATA 

# --- CALCULATE RELATIVE AGE --- 
# relative age = number of years since first detection for each individual
gbcons_age <- first_detection %>%
  arrange(individual, year) %>%
  group_by(individual) %>%
  mutate(relative_age = year - min(year) + 1)

# Remove rows where Salmon.med is NA
gbcons_age <- gbcons_age %>%
  filter(!is.na(Salmon.med))

# ---- scale predictors for modelling ---
scaled_df <- gbcons_age

scaled_df$EVI_Mean_z     <- as.numeric(scale(scaled_df$EVI_Mean_MWS))
scaled_df$BMDensity_z    <- as.numeric(scale(scaled_df$YearlyBMDensityinMWS))
scaled_df$LogBM_z        <- as.numeric(scale(scaled_df$LogBMDensityinMWS))
scaled_df$diversity_z    <- as.numeric(scale(scaled_df$Salmon_diversity_MWS))
scaled_df$relative_age_z <- as.numeric(scale(scaled_df$relative_age))

# Verify scaling worked (should be ~0 and ~1)
mean(scaled_df$relative_age_z)
sd(scaled_df$relative_age_z)

# ==============================================================================
## OUR DATA FRAME DOESNT HAVE A SEX COLUMN YET. ADD THIS HERE!
## this section has been commented out due to sensitive location info
## see output on line 327

## NOT SHARED: individual bear sex data:
#sex <- read.csv("data/raw/sensitivity data/Gizzly_consumption_2009_to_2019.csv")

## Join sex data
scaled_df_with_sex <- scaled_df %>%
  left_join(
    sex %>% dplyr::select(individual, year, sex),
    by = c("individual", "year")
  )

## clean df
all_covariotes <- scaled_df_with_sex %>%
  dplyr::select(year, individual, sex, Salmon.med, project, WSHDFID, FOCALWSHD, site_id, revisitid, revisit, revisit_date,
                Salmon_year, YearlyBMinMWS, YearlyBMDensityinMWS, LogBMDensityinMWS,
                SalmonDiv_year, Salmon_diversity_MWS, EVI_year, EVI_Mean_MWS, relative_age, EVI_Mean_z, BMDensity_z, LogBM_z, diversity_z, relative_age_z) %>%
  rename(
    BM        = YearlyBMinMWS,
    BMDensity = YearlyBMDensityinMWS,
    LogBM     = LogBMDensityinMWS,
    diversity = Salmon_diversity_MWS
  )
# ==============================================================================

# --- SAVE FILE FOR ANALYSIS --- #
#write.csv(all_covariotes, here("data", "processed", "all_cov_scaled_MWS_V3.csv"))

# Save stripped version for publication — removes UTM coordinates

write.csv(all_covariotes, here("data", "processed", "all_cov_scaled_no_locations.csv"))
          