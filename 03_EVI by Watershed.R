####################
#
# This code is for attributing EVI at watersheds
# Extracts mean, median, and maximum EVI from MODIS satellite imagery
# (MOD13Q1, 250m, 16-day composites) for each merged watershed. Mean EVI is used 
# as a covariate in GLMM models.
#
# 
####################

library(here)
library(terra)
library(tidyverse)
library(sf)
library(dplyr)

###################

# Load Raw EVI layers folder 
evi_folder <- here("data", "raw", "RAW_EVI_APRIL22_JULY27")

# List all EVI TIFF files for each year (filtering by year in filename) --> call in specific year
evi_files_2008 <- list.files(evi_folder, 
                             pattern = "MOD13Q1.*doy2008.*\\.tif$", 
                             full.names = TRUE)
evi_files_2009 <- list.files(evi_folder, 
                             pattern = "MOD13Q1.*doy2009.*\\.tif$", 
                             full.names = TRUE)
evi_files_2010 <- list.files(evi_folder, 
                             pattern = "MOD13Q1.*doy2010.*\\.tif$", 
                             full.names = TRUE)
evi_files_2011 <- list.files(evi_folder, 
                             pattern = "MOD13Q1.*doy2011.*\\.tif$", 
                             full.names = TRUE)
evi_files_2012 <- list.files(evi_folder, 
                             pattern = "MOD13Q1.*doy2012.*\\.tif$", 
                             full.names = TRUE)
evi_files_2013 <- list.files(evi_folder, 
                             pattern = "MOD13Q1.*doy2013.*\\.tif$", 
                             full.names = TRUE)
evi_files_2014 <- list.files(evi_folder, 
                             pattern = "MOD13Q1.*doy2014.*\\.tif$", 
                             full.names = TRUE)
evi_files_2015 <- list.files(evi_folder, 
                             pattern = "MOD13Q1.*doy2015.*\\.tif$", 
                             full.names = TRUE)
evi_files_2016 <- list.files(evi_folder, 
                             pattern = "MOD13Q1.*doy2016.*\\.tif$", 
                             full.names = TRUE)
evi_files_2017 <- list.files(evi_folder, 
                             pattern = "MOD13Q1.*doy2017.*\\.tif$", 
                             full.names = TRUE)
evi_files_2018 <- list.files(evi_folder, 
                             pattern = "MOD13Q1.*doy2018.*\\.tif$", 
                             full.names = TRUE)


r2008 = rast(evi_files_2008)  # load rasters into a single multi-layer SpatRaster2008 folder of all tiff files
r2009 = rast(evi_files_2009)
r2010 = rast(evi_files_2010)
r2011 = rast(evi_files_2011)
r2012 = rast(evi_files_2012)
r2013 = rast(evi_files_2013)
r2014 = rast(evi_files_2014)
r2015 = rast(evi_files_2015)
r2016 = rast(evi_files_2016)
r2017 = rast(evi_files_2017)
r2018 = rast(evi_files_2018)

# Load LRGwatershed shapefile 
WSHEDS <- vect(here("data", "processed", "LRG_Watersheds_for_all_sites", "LRG_Watersheds_all_sites_FINAL.shp"))  

# Check CRS - reproject the EVI data to CRS of watersheds:
r2008 <- project(r2008, crs(WSHEDS))
r2009 <- project(r2009, crs(WSHEDS))
r2010 <- project(r2010, crs(WSHEDS))
r2011 <- project(r2011, crs(WSHEDS))
r2012 <- project(r2012, crs(WSHEDS))
r2013 <- project(r2013, crs(WSHEDS))
r2014 <- project(r2014, crs(WSHEDS))
r2015 <- project(r2015, crs(WSHEDS))
r2016 <- project(r2016, crs(WSHEDS))
r2017 <- project(r2017, crs(WSHEDS))
r2018 <- project(r2018, crs(WSHEDS))

#extract the raw values from the rasters. The 'ID' column is a number, here 1 to 12 which is how many polygons there are in the layer 'v' above.
#ext_r20xx is just a big table with all of the extracted pixel values, we can use this to calculate the median and max
ext_r2008 = terra::extract(r2008, WSHEDS, method = 'simple', exact = TRUE, touches = TRUE) %>%
  mutate(FID = WSHEDS$WTRSHD_FID[ID]) %>% 
  select(-fraction) #the exact method gives you the fraction that is covered by or touches a polygon line. You could just select all those where fraction = 1 if you want to be very exact, but it cuts off the edges
ext_r2009 = terra::extract(r2009, WSHEDS, method = 'simple', exact = TRUE, touches = TRUE) %>%
  mutate(FID = WSHEDS$WTRSHD_FID[ID]) %>% 
  select(-fraction) 
ext_r2010 = terra::extract(r2010, WSHEDS, method = 'simple', exact = TRUE, touches = TRUE) %>%
  mutate(FID = WSHEDS$WTRSHD_FID[ID]) %>% 
  select(-fraction) 
ext_r2011 = terra::extract(r2011, WSHEDS, method = 'simple', exact = TRUE, touches = TRUE) %>%
  mutate(FID = WSHEDS$WTRSHD_FID[ID]) %>% 
  select(-fraction) 
ext_r2012 = terra::extract(r2012, WSHEDS, method = 'simple', exact = TRUE, touches = TRUE) %>%
  mutate(FID = WSHEDS$WTRSHD_FID[ID]) %>% 
  select(-fraction) 
ext_r2013 = terra::extract(r2013, WSHEDS, method = 'simple', exact = TRUE, touches = TRUE) %>%
  mutate(FID = WSHEDS$WTRSHD_FID[ID]) %>% 
  select(-fraction) 
ext_r2014 = terra::extract(r2014, WSHEDS, method = 'simple', exact = TRUE, touches = TRUE) %>%
  mutate(FID = WSHEDS$WTRSHD_FID[ID]) %>% 
  select(-fraction) 
ext_r2015 = terra::extract(r2015, WSHEDS, method = 'simple', exact = TRUE, touches = TRUE) %>%
  mutate(FID = WSHEDS$WTRSHD_FID[ID]) %>% 
  select(-fraction) 
ext_r2016 = terra::extract(r2016, WSHEDS, method = 'simple', exact = TRUE, touches = TRUE) %>%
  mutate(FID = WSHEDS$WTRSHD_FID[ID]) %>% 
  select(-fraction) 
ext_r2017 = terra::extract(r2017, WSHEDS, method = 'simple', exact = TRUE, touches = TRUE) %>%
  mutate(FID = WSHEDS$WTRSHD_FID[ID]) %>% 
  select(-fraction) 
ext_r2018 = terra::extract(r2018, WSHEDS, method = 'simple', exact = TRUE, touches = TRUE) %>%
  mutate(FID = WSHEDS$WTRSHD_FID[ID]) %>% 
  select(-fraction) 


#calculate mean, median and max for each year
val_2008 <- ext_r2008 %>% #this is our extracted data frame
  group_by(ID) %>% #
  summarize( #we want to summarize the data
    median_EVI2008 = median(unlist(across(where(is.numeric))), na.rm = TRUE), #we want the median of all the images, so we unlist them and calculate across all values, where(is .numeric()) is just a helper function
    max_EVI2008 = max(unlist(across(where(is.numeric))), na.rm = TRUE), #same for max
    mean_EVI2008 = mean(unlist(across(where(is.numeric))), na.rm = TRUE) #same for mean
  )

val_2009 <- ext_r2009 %>% #this is our extracted data frame
  group_by(ID) %>% 
  summarize(
    median_EVI2009 = median(unlist(across(where(is.numeric))), na.rm = TRUE),
    max_EVI2009 = max(unlist(across(where(is.numeric))), na.rm = TRUE),
    mean_EVI2009 = mean(unlist(across(where(is.numeric))), na.rm = TRUE)
  )

val_2010 <- ext_r2010 %>% #this is our extracted data frame
  group_by(ID) %>% 
  summarize(
    median_EVI2010 = median(unlist(across(where(is.numeric))), na.rm = TRUE),
    max_EVI2010 = max(unlist(across(where(is.numeric))), na.rm = TRUE),
    mean_EVI2010 = mean(unlist(across(where(is.numeric))), na.rm = TRUE))

val_2011 <- ext_r2011 %>% #this is our extracted data frame
  group_by(ID) %>% 
  summarize(
    median_EVI2011 = median(unlist(across(where(is.numeric))), na.rm = TRUE),
    max_EVI2011 = max(unlist(across(where(is.numeric))), na.rm = TRUE),
    mean_EVI2011 = mean(unlist(across(where(is.numeric))), na.rm = TRUE))

val_2012 <- ext_r2012 %>% #this is our extracted data frame
  group_by(ID) %>% 
  summarize(
    median_EVI2012 = median(unlist(across(where(is.numeric))), na.rm = TRUE),
    max_EVI2012 = max(unlist(across(where(is.numeric))), na.rm = TRUE),
    mean_EVI2012 = mean(unlist(across(where(is.numeric))), na.rm = TRUE))

val_2013 <- ext_r2013 %>% #this is our extracted data frame
  group_by(ID) %>% 
  summarize(
    median_EVI2013 = median(unlist(across(where(is.numeric))), na.rm = TRUE),
    max_EVI2013 = max(unlist(across(where(is.numeric))), na.rm = TRUE),
    mean_EVI2013 = mean(unlist(across(where(is.numeric))), na.rm = TRUE))

val_2014 <- ext_r2014 %>% #this is our extracted data frame
  group_by(ID) %>% 
  summarize(
    median_EVI2014 = median(unlist(across(where(is.numeric))), na.rm = TRUE),
    max_EVI2014 = max(unlist(across(where(is.numeric))), na.rm = TRUE),
    mean_EVI2014 = mean(unlist(across(where(is.numeric))), na.rm = TRUE))

val_2015 <- ext_r2015 %>% #this is our extracted data frame
  group_by(ID) %>% 
  summarize(
    median_EVI2015 = median(unlist(across(where(is.numeric))), na.rm = TRUE),
    max_EVI2015 = max(unlist(across(where(is.numeric))), na.rm = TRUE),
    mean_EVI2015 = mean(unlist(across(where(is.numeric))), na.rm = TRUE))

val_2016 <- ext_r2016 %>% #this is our extracted data frame
  group_by(ID) %>% 
  summarize(
    median_EVI2016 = median(unlist(across(where(is.numeric))), na.rm = TRUE),
    max_EVI2016 = max(unlist(across(where(is.numeric))), na.rm = TRUE),
    mean_EVI2016 = mean(unlist(across(where(is.numeric))), na.rm = TRUE))

val_2017 <- ext_r2017 %>% #this is our extracted data frame
  group_by(ID) %>% 
  summarize(
    median_EVI2017 = median(unlist(across(where(is.numeric))), na.rm = TRUE),
    max_EVI2017 = max(unlist(across(where(is.numeric))), na.rm = TRUE),
    mean_EVI2017 = mean(unlist(across(where(is.numeric))), na.rm = TRUE))

val_2018 <- ext_r2018 %>% #this is our extracted data frame
  group_by(ID) %>% 
  summarize(
    median_EVI2018 = median(unlist(across(where(is.numeric))), na.rm = TRUE),
    max_EVI2018 = max(unlist(across(where(is.numeric))), na.rm = TRUE),
    mean_EVI2018 = mean(unlist(across(where(is.numeric))), na.rm = TRUE))


# print the values, you'll have the mean, median and max for each of the polygons/watersheds
print(val_2008)

## Merge all the data frames by ID
yearly_vals <- list(val_2008, val_2009, val_2010, val_2011, val_2012, 
                    val_2013, val_2014, val_2015, val_2016, val_2017, val_2018)

# Use reduce() to iteratively apply full_join on all data frames
EVI_by_year_in_merged_WSHDs <- reduce(yearly_vals, full_join, by = "ID") %>%
  rename(WTRSHD_FID = ID)

##
EVI_mean <- EVI_by_year_in_merged_WSHDs[,c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34)]
EVI_median <- EVI_by_year_in_merged_WSHDs[,c(1, 2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32)]
EVI_max <- EVI_by_year_in_merged_WSHDs[,c(1, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33)]

merged_EVI_mean <-select(EVI_mean, WTRSHD_FID, mean_EVI2008, mean_EVI2009, mean_EVI2010, mean_EVI2011, mean_EVI2012, mean_EVI2013, mean_EVI2014, mean_EVI2015, mean_EVI2016, mean_EVI2017, mean_EVI2018)%>%
  gather(year, mean , -WTRSHD_FID)%>%
  mutate(year=as.numeric(substr(year, 9, 12)))

merged_EVI_median <-select(EVI_median, WTRSHD_FID, median_EVI2008, median_EVI2009, median_EVI2010, median_EVI2011, median_EVI2012, median_EVI2013, median_EVI2014, median_EVI2015, median_EVI2016, median_EVI2017, median_EVI2018)%>%
  gather(year, median , -WTRSHD_FID)%>%
  mutate(year=as.numeric(substr(year, 11, 14)))

merged_EVI_max <-select(EVI_max, WTRSHD_FID, max_EVI2008, max_EVI2009, max_EVI2010, max_EVI2011, max_EVI2012, max_EVI2013, max_EVI2014, max_EVI2015, max_EVI2016, max_EVI2017, max_EVI2018)%>%
  gather(year, max , -WTRSHD_FID)%>%
  mutate(year=as.numeric(substr(year, 8, 11)))


## SAVE FILE 
write.csv(EVI_by_year_in_merged_WSHDs, here("data", "processed", "EVI", "EVI_by_year_in_merged_WSHDs.csv"))
write.csv(merged_EVI_mean, here("data", "processed", "EVI", "Mean_EVI_by_year_in_merged_WSHDs.csv"))
write.csv(merged_EVI_median, here("data", "processed", "EVI", "Median_EVI_by_year_in_merged_WSHDs.csv"))
write.csv(merged_EVI_max, here("data", "processed", "EVI", "Max_EVI_by_year_in_merged_WSHDs.csv"))

## SAVE AS SHP
#merge EVI data with watersheds using "ID" column
WSHEDS <- st_as_sf(WSHEDS)
watersheds_evi <- merge(WSHEDS, EVI_by_year_in_merged_WSHDs, by.x = "WSHDFID", by.y = "WTRSHD_FID")
mean_evi <- merge(WSHEDS, merged_EVI_mean, by.x = "WSHDFID", by.y = "WTRSHD_FID")
median_evi <- merge(WSHEDS, merged_EVI_median,by.x = "WSHDFID", by.y = "WTRSHD_FID")
max_evi <- merge(WSHEDS, merged_EVI_max, by.x = "WSHDFID", by.y = "WTRSHD_FID")

st_write(watersheds_evi, here("data", "processed", "EVI", "evi_shps", "EVI_by_year_in_merged_WSHDs.shp"))
st_write(mean_evi, here("data", "processed", "EVI", "evi_shps", "Mean_EVI_by_year_in_merged_WSHDs.shp"))
st_write(median_evi, here("data", "processed", "EVI", "evi_shps", "Median_EVI_by_year_in_merged_WSHDs.shp"))
st_write(max_evi, here("data", "processed", "EVI", "evi_shps", "Max_EVI_by_year_in_merged_WSHDs.shp"))


