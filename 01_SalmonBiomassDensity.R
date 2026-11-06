#####################
#
# # This script calculates salmon biomass and salmon biomass density per watershed
# and year. It starts from the focal watershed layer (focal_clipped.shp) in
# data/processed/, which was created from hair-snag locations
#
# Data that went into this is the imputed salmon biomass data from Heather Bryan
# Code written by Monica Short
#
# 
#####################

library(here)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)

#####################

#Intersect salmon with shapefiles

#read in the focal watersheds (made in 00_sensitive.R) and the merged (large) watersheds
#the focal layer must have ONE row per focal watershed: stream lengths are summed over it below,
wsh_clip <- st_read(here("data", "processed", "Focal Watersheds in Study Area", "focal_clipped.shp")) %>%
  distinct(WTRSHD_FID, .keep_all = TRUE)
stopifnot(!anyDuplicated(wsh_clip$WTRSHD_FID))

# merged watersheds:
usable_WSH <- st_read(here("data", "processed", "LRG_Watersheds_for_all_sites", "LRG_Watersheds_all_sites_FINAL.shp")) %>%
  st_transform(crs = st_crs(wsh_clip))

#read in salmon data
salm <- st_read(here("data","processed", "combinedStreamsBMTotal_NoZeros2019", "combinedStreamsBMTotal_NoZeros2019.shp"))

#remove columns we don't need
salm <- salm[,c(56:70)]

#ensure the crs is the same as our focal watershed layer
salm <- st_transform(salm, crs = st_crs(wsh_clip))

#intersect salmon streams with the focal watersheds
focal_salm <- st_intersection(salm, wsh_clip)

#calculate the length of the stream that is in each focal watershed
focal_salm$Shape_Leng <- st_length(focal_salm)
focal_salm$Shape_Leng <- as.numeric(focal_salm$Shape_Leng)

focal_salm_df <- as.data.frame(focal_salm)

#intersect salmon streams with the merged watersheds
merg_salm <- st_intersection(salm, usable_WSH)

#calculate the length of each stream that is in each merged watershed
merg_salm$Shape_Leng <- st_length(merg_salm)
merg_salm$Shape_Leng <- as.numeric(merg_salm$Shape_Leng)

merg_salm_df <- as.data.frame(merg_salm)

#######################################################################################
# Now summarize salmon biomass by watershed: 
# Biomass is calculated at the focal and merged watershed scale for completeness but is not used in 
# final analysis
#--------------------------------------------------------------------------------------

## Calculate the total stream length for each stream using the values from the focal watershed layer
StreamLengths <- focal_salm_df %>% 
  group_by(Waterbody) %>% 
  summarise(Stream_Length = sum(Shape_Leng))

StreamLengths <- StreamLengths[which(StreamLengths$Waterbody != ""),]

## also need the stream lengths for unnamed water bodies
unNamed <- focal_salm_df[which(is.na(focal_salm_df$Waterbody) == TRUE),]
unNamedLengths <- unNamed %>% 
  group_by(BM2008) %>% 
  summarise(Stream_Length1 = sum(Shape_Leng))

## add the stream lengths to the dataframes
#focal_salm_df$Stream_Length0 <- StreamLengths$Stream_Length[match(focal_salm_df$Waterbody, StreamLengths$Waterbody)]
#focal_salm_df$Stream_Length1 <- unNamedLengths$Stream_Length1[match(focal_salm_df$BM2008, unNamedLengths$BM2008)]
merg_salm_df$Stream_Length0 <- StreamLengths$Stream_Length[match(merg_salm_df$Waterbody, StreamLengths$Waterbody)]
merg_salm_df$Stream_Length1 <- unNamedLengths$Stream_Length1[match(merg_salm_df$BM2008, unNamedLengths$BM2008)]

#focal_salm_df$Stream_Length <- rowSums(focal_salm_df[,c("Stream_Length0", "Stream_Length1")], na.rm = TRUE)
merg_salm_df$Stream_Length <- rowSums(merg_salm_df[,c("Stream_Length0", "Stream_Length1")], na.rm = TRUE)

## remove columns that we don't need
#focal_salm_df <- focal_salm_df%>%
  #select(-Stream_Length0, -Stream_Length1)
merg_salm_df <- merg_salm_df %>%
  select(-Stream_Length0, -Stream_Length1)

## Rotate the biomass columns and have a biomass value for each year for each watershed. 
#focal_salm_tidy<-select(focal_salm_df, which(names(focal_salm_df)=="BM2008"):which(names(focal_salm_df)=="BM2019"), nSalmPops, Waterbody, WTRSHD_FID, AREA_HA, F_CODE, OBJECTID, AREA_SQM, FEAT_LEN, geometry, Shape_Leng, Stream_Length)%>%
  #gather(year, biomass, -nSalmPops, -Waterbody, -WTRSHD_FID, -AREA_HA, -F_CODE, -OBJECTID, -AREA_SQM, -FEAT_LEN, -geometry, -Shape_Leng, -Stream_Length)%>%
  #mutate(year=as.numeric(substr(year, 3, 6)))

merg_salm_tidy<-select(merg_salm_df, which(names(merg_salm_df)=="BM2008"):which(names(merg_salm_df)=="BM2019"), nSalmPops, Waterbody, WSHDFID, AREAKM2, FOCALWSHD, rociceKM2, usblAREA, geometry, Shape_Leng, Stream_Length)%>%
  gather(year, biomass, -nSalmPops, -Waterbody, -WSHDFID, -AREAKM2, -FOCALWSHD, -rociceKM2, -usblAREA, -geometry, -Shape_Leng, -Stream_Length)%>%
  mutate(year=as.numeric(substr(year, 3, 6)))

## calculate proportion of each stream length that is in each watershed
#focal_salm_tidy$ProportioninFWS <- focal_salm_tidy$Shape_Leng/focal_salm_tidy$Stream_Length
merg_salm_tidy$ProportioninMWS <- merg_salm_tidy$Shape_Leng/merg_salm_tidy$Stream_Length

## multiply the salmon biomass by the proportion of the stream length in that watershed
#focal_salm_tidy$biomassinFWS <- focal_salm_tidy$biomass * focal_salm_tidy$ProportioninFWS
merg_salm_tidy$biomassinMWS <- merg_salm_tidy$biomass * merg_salm_tidy$ProportioninMWS

## Sum biomass across streams within a watershed and year
#bm_focal_salm <- focal_salm_tidy %>%
  #group_by(WTRSHD_FID, year) %>%
  #summarise(YearlyBMinFWS=sum(biomassinFWS))

bm_merg_salm <- merg_salm_tidy %>%
  group_by(WSHDFID, year) %>%
  summarise(YearlyBMinMWS=sum(biomassinMWS))

## check that it worked
#length(unique(bm_focal_salm$year))
length(unique(bm_merg_salm$year))
#[1] 12 years
#length(unique(bm_focal_salm$WTRSHD_FID))
#[1] 270 watersheds
length(unique(bm_merg_salm$WSHDFID))
#[1] 216 large watersheds

## biomass estimate for each focal watershed and year
#write.csv(bm_focal_salm, here("data", "processed", "focal_salm_biomass_per_year.csv"))

## we have one biomass estimate for each larger watershed and year
write.csv(bm_merg_salm, here("data","processed", "merg_salm_biomass_per_year.csv"))

##--------------------------------------------------------------------------------------------------------
## Merge the spatial component back into the biomass data
#bm_focal_spatial <- merge(wsh_clip, bm_focal_salm, by = "WTRSHD_FID")
#bm_focal_spatial$AREA_KM2 <- bm_focal_spatial$AREA_SQM/1000000
#bm_focal_spatial <- bm_focal_spatial[,c(1, 7:12)]

## Now lets calculate salmon biomass density 
#this is salmon biomass per km
#bm_focal_spatial$YearlyBMDensityinFWS <- bm_focal_spatial$YearlyBMinFWS/ bm_focal_spatial$AREA_KM2

## Lets also calculate a log transformed value for salmon biomass density
#bm_focal_spatial$LogBMDensityinFWS <- log(bm_focal_spatial$YearlyBMDensityinFWS +1) # +1 because there are lots of zeros

## now have one biomass density estimate (and log transformed version) for each focal watershed and year
#st_write(bm_focal_spatial, here("data", "processed", "focal_salm_biomass_per_year_w_density.csv"))

#--------------------------------------------------------------------------------------------------------
#Merge the spatial data back into the biomass data
bm_merge_spatial <- merge(usable_WSH, bm_merg_salm, by = "WSHDFID", all = TRUE)

#Now lets calculate salmon biomass density 
#this is salmon biomass per km
bm_merge_spatial$YearlyBMDensityinMWS <- bm_merge_spatial$YearlyBMinMWS/ bm_merge_spatial$usblAREA

#Lets also calculate a log transformed value for salmon biomass density
bm_merge_spatial$LogBMDensityinMWS <- log(bm_merge_spatial$YearlyBMDensityinMWS +1) # +1 because there are lots of zeros


#we also have one biomass density estimate (and log transformed version) for each larger watershed and year
st_write(bm_merge_spatial, here("data", "processed", "merg_salm_biomass_per_year_w_density.csv"), delete_dsn = TRUE)

#------------------------------------------------------------------------------------------------
#Plotting the larger watersheds

#I want to plot the salmon biomass across the watersheds each year
#for ease of plotting, I'm going to split this into years.
#I'm keeping NA watersheds in here so we can plot all the larger watershed outlines
bm_2008MWS <- bm_merge_spatial[which(bm_merge_spatial$year == "2008" | is.na(bm_merge_spatial$year) == TRUE),]
bm_2009MWS <- bm_merge_spatial[which(bm_merge_spatial$year == "2009" | is.na(bm_merge_spatial$year) == TRUE),]
bm_2010MWS <- bm_merge_spatial[which(bm_merge_spatial$year == "2010" | is.na(bm_merge_spatial$year) == TRUE),]
bm_2011MWS <- bm_merge_spatial[which(bm_merge_spatial$year == "2011" | is.na(bm_merge_spatial$year) == TRUE),]
bm_2012MWS <- bm_merge_spatial[which(bm_merge_spatial$year == "2012" | is.na(bm_merge_spatial$year) == TRUE),]
bm_2013MWS <- bm_merge_spatial[which(bm_merge_spatial$year == "2013" | is.na(bm_merge_spatial$year) == TRUE),]
bm_2014MWS <- bm_merge_spatial[which(bm_merge_spatial$year == "2014" | is.na(bm_merge_spatial$year) == TRUE),]
bm_2015MWS <- bm_merge_spatial[which(bm_merge_spatial$year == "2015" | is.na(bm_merge_spatial$year) == TRUE),]
bm_2016MWS <- bm_merge_spatial[which(bm_merge_spatial$year == "2016" | is.na(bm_merge_spatial$year) == TRUE),]
bm_2017MWS <- bm_merge_spatial[which(bm_merge_spatial$year == "2017" | is.na(bm_merge_spatial$year) == TRUE),]
bm_2018MWS <- bm_merge_spatial[which(bm_merge_spatial$year == "2018" | is.na(bm_merge_spatial$year) == TRUE),]
bm_2019MWS <- bm_merge_spatial[which(bm_merge_spatial$year == "2019" | is.na(bm_merge_spatial$year) == TRUE),]

#Plot salmon biomass across LARGER, MERGED watersheds
ggplot()+
  geom_sf(data = bm_2008MWS, aes(fill = YearlyBMinMWS))+
  scale_fill_viridis_c(limits = c(0, 3276000), na.value = "white")+
  labs(x = "", y = "", fill = "Total Salmon Biomass \nin Merged Watershed\n2008") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2009MWS, aes(fill = YearlyBMinMWS))+
  scale_fill_viridis_c(limits = c(0, 3276000), na.value = "white")+
  labs(x = "", y = "", fill = "Total Salmon Biomass \nin Merged Watershed\n2009") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2010MWS, aes(fill = YearlyBMinMWS))+
  scale_fill_viridis_c(limits = c(0, 3276000), na.value = "white")+
  labs(x = "", y = "", fill = "Total Salmon Biomass \nin Merged Watershed\n2010") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2011MWS, aes(fill = YearlyBMinMWS))+
  scale_fill_viridis_c(limits = c(0, 3276000), na.value = "white")+
  labs(x = "", y = "", fill = "Total Salmon Biomass \nin Merged Watershed\n2011") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2012MWS, aes(fill = YearlyBMinMWS))+
  scale_fill_viridis_c(limits = c(0, 3276000), na.value = "white")+
  labs(x = "", y = "", fill = "Total Salmon Biomass \nin Merged Watershed\n2012") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2013MWS, aes(fill = YearlyBMinMWS))+
  scale_fill_viridis_c(limits = c(0, 3276000), na.value = "white")+
  labs(x = "", y = "", fill = "Total Salmon Biomass \nin Merged Watershed\n2013") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2014MWS, aes(fill = YearlyBMinMWS))+
  scale_fill_viridis_c(limits = c(0, 3276000), na.value = "white")+
  labs(x = "", y = "", fill = "Total Salmon Biomass \nin Merged Watershed\n2014") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2015MWS, aes(fill = YearlyBMinMWS))+
  scale_fill_viridis_c(limits = c(0, 3276000), na.value = "white")+
  labs(x = "", y = "", fill = "Total Salmon Biomass \nin Merged Watershed\n2015") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2016MWS, aes(fill = YearlyBMinMWS))+
  scale_fill_viridis_c(limits = c(0, 3276000), na.value = "white")+
  labs(x = "", y = "", fill = "Total Salmon Biomass \nin Merged Watershed\n2016") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2017MWS, aes(fill = YearlyBMinMWS))+
  scale_fill_viridis_c(limits = c(0, 3276000), na.value = "white")+
  labs(x = "", y = "", fill = "Total Salmon Biomass \nin Merged Watershed\n2017") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2018MWS, aes(fill = YearlyBMinMWS))+
  scale_fill_viridis_c(limits = c(0, 3276000), na.value = "white")+
  labs(x = "", y = "", fill = "Total Salmon Biomass \nin Merged Watershed\n2018") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2019MWS, aes(fill = YearlyBMinMWS))+
  scale_fill_viridis_c(limits = c(0, 3276000), na.value = "white")+
  labs(x = "", y = "", fill = "Total Salmon Biomass \nin Merged Watershed\n2019") +
  theme_minimal()

#Plot salmon biomass density across MERGED watersheds
ggplot()+
  geom_sf(data = bm_2008MWS, aes(fill = YearlyBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 10310), na.value = "white")+
  labs(x = "", y = "", fill = "Salmon Biomass \nDensity (kg/km2) \nin Merged Watershed\n2008") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2009MWS, aes(fill = YearlyBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 10310), na.value = "white")+
  labs(x = "", y = "", fill = "Salmon Biomass \nDensity (kg/km2) \nin Merged Watershed\n2009") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2010MWS, aes(fill = YearlyBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 10310), na.value = "white")+
  labs(x = "", y = "", fill = "Salmon Biomass \nDensity (kg/km2) \nin Merged Watershed\n2010") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2011MWS, aes(fill = YearlyBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 10310), na.value = "white")+
  labs(x = "", y = "", fill = "Salmon Biomass \nDensity (kg/km2) \nin Merged Watershed\n2011") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2012MWS, aes(fill = YearlyBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 10310), na.value = "white")+
  labs(x = "", y = "", fill = "Salmon Biomass \nDensity (kg/km2) \nin Merged Watershed\n2012") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2013MWS, aes(fill = YearlyBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 10310), na.value = "white")+
  labs(x = "", y = "", fill = "Salmon Biomass \nDensity (kg/km2) \nin Merged Watershed\n2013") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2014MWS, aes(fill = YearlyBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 10310), na.value = "white")+
  labs(x = "", y = "", fill = "Salmon Biomass \nDensity (kg/km2) \nin Merged Watershed\n2014") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2015MWS, aes(fill = YearlyBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 10310), na.value = "white")+
  labs(x = "", y = "", fill = "Salmon Biomass \nDensity (kg/km2) \nin Merged Watershed\n2015") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2016MWS, aes(fill = YearlyBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 10310), na.value = "white")+
  labs(x = "", y = "", fill = "Salmon Biomass \nDensity (kg/km2) \nin Merged Watershed\n2016") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2017MWS, aes(fill = YearlyBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 10310), na.value = "white")+
  labs(x = "", y = "", fill = "Salmon Biomass \nDensity (kg/km) \nin Merged Watershed\n2017") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2018MWS, aes(fill = YearlyBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 10310), na.value = "white")+
  labs(x = "", y = "", fill = "Salmon Biomass \nDensity (kg/km) \nin Merged Watershed\n2018") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2019MWS, aes(fill = YearlyBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 10310), na.value = "white")+
  labs(x = "", y = "", fill = "Salmon Biomass \nDensity (kg/km) \nin Merged Watershed\n2019") +
  theme_minimal()


#Plot log of salmon biomass density across Merged watersheds
ggplot()+
  geom_sf(data = bm_2008MWS, aes(fill = LogBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 9.25), na.value = "white")+
  labs(x = "", y = "", fill = "Log of Salmon \nBiomass Density  \n(kg/km2) in Merged \nWatershed 2008") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2009MWS, aes(fill = LogBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 9.25), na.value = "white")+
  labs(x = "", y = "", fill = "Log of Salmon \nBiomass Density  \n(kg/km2) in Merged \nWatershed 2009") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2010MWS, aes(fill = LogBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 9.25), na.value = "white")+
  labs(x = "", y = "", fill = "Log of Salmon \nBiomass Density  \n(kg/km2) in Merged \nWatershed 2010") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2011MWS, aes(fill = LogBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 9.25), na.value = "white")+
  labs(x = "", y = "", fill = "Log of Salmon \nBiomass Density  \n(kg/km2) in Merged \nWatershed 2011") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2012MWS, aes(fill = LogBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 9.25), na.value = "white")+
  labs(x = "", y = "", fill = "Log of Salmon \nBiomass Density  \n(kg/km2) in Merged \nWatershed 2012") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2013MWS, aes(fill = LogBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 9.25), na.value = "white")+
  labs(x = "", y = "", fill = "Log of Salmon \nBiomass Density  \n(kg/km2) in Merged \nWatershed 2013") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2014MWS, aes(fill = LogBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 9.25), na.value = "white")+
  labs(x = "", y = "", fill = "Log of Salmon \nBiomass Density  \n(kg/km2) in Merged \nWatershed 2014") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2015MWS, aes(fill = LogBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 9.25), na.value = "white")+
  labs(x = "", y = "", fill = "Log of Salmon \nBiomass Density  \n(kg/km2) in Merged \nWatershed 2015") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2016MWS, aes(fill = LogBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 9.25), na.value = "white")+
  labs(x = "", y = "", fill = "Log of Salmon \nBiomass Density  \n(kg/km2) in Merged \nWatershed 2016") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2017MWS, aes(fill = LogBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 9.25), na.value = "white")+
  labs(x = "", y = "", fill = "Log of Salmon \nBiomass Density  \n(kg/km2) in Merged \nWatershed 2017") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2018MWS, aes(fill = LogBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 9.25), na.value = "white")+
  labs(x = "", y = "", fill = "Log of Salmon \nBiomass Density  \n(kg/km2) in Merged \nWatershed 2018") +
  theme_minimal()

ggplot()+
  geom_sf(data = bm_2019MWS, aes(fill = LogBMDensityinMWS))+
  scale_fill_viridis_c(limits = c(0, 9.25), na.value = "white")+
  labs(x = "", y = "", fill = "Log of Salmon \nBiomass Density  \n(kg/km2) in Merged \nWatershed 2019") +
  theme_minimal()

