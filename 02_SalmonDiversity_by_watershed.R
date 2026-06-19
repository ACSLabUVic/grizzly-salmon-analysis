#####################
#
# calculates Shannon diversity of salmon species for merged watersheds.
# 
#####################

library(readr)
library(here)
library(dplyr)
library(sf)

#read in salmon data that is separated by salmon species. 
salm_spp <- st_read(here("data", "processed",  "combinedStreams2019", "combinedStreams2019.shp"))

#select only the columns we want
salm_spp <- salm_spp[,c(1:3, 128:139, 141)]

# read in shapefiles for watersheds
focal_wsh <- st_read(here("data", "processed", "Focal Watersheds in Study Area", "focal_clipped.shp"))
merged_wsh <- st_read(here("data", "processed", "LRG_Watersheds_for_all_sites", "LRG_Watersheds_all_sites_FINAL.shp"))

#ensure the crs is same for all layers
merged_wsh <- st_transform(merged_wsh, crs = st_crs(focal_wsh))
salm_spp <- st_transform(salm_spp, crs = st_crs(focal_wsh))

## Intersect salmon with shapefiles
#--------------------------------------------------------------------------------------
#intersect salmon streams with the focal watersheds
focal_salm_spp <- st_intersection(salm_spp, focal_wsh)

#calculate the length of the stream that is in each focal watershed
focal_salm_spp$Shape_Leng <- st_length(focal_salm_spp)
focal_salm_spp$Shape_Leng <- as.numeric(focal_salm_spp$Shape_Leng)

focal_salm_spp_df <- as.data.frame(focal_salm_spp)

#intersect salmon streams with the merged watersheds
merg_salm_spp <- st_intersection(salm_spp, merged_wsh)

#calculate the length of each stream that is in each merged watershed
merg_salm_spp$Shape_Leng <- st_length(merg_salm_spp)
merg_salm_spp$Shape_Leng <- as.numeric(merg_salm_spp$Shape_Leng)

merg_salm_spp_df <- as.data.frame(merg_salm_spp)
                                  
#--------------------------------------------------------------------------------------
#the salmon dataset including the species also has a bunch of zeros we want to remove. 
 
focal_spp <- focal_salm_spp_df[which(focal_salm_spp_df$BM2008 != 0 |
                           focal_salm_spp_df$BM2009 != 0 |
                           focal_salm_spp_df$BM2010 != 0 |
                           focal_salm_spp_df$BM2011 != 0 |
                           focal_salm_spp_df$BM2012 != 0 |
                           focal_salm_spp_df$BM2013 != 0 |
                           focal_salm_spp_df$BM2014 != 0 |
                           focal_salm_spp_df$BM2015 != 0 |
                           focal_salm_spp_df$BM2016 != 0 |
                           focal_salm_spp_df$BM2017 != 0 |
                           focal_salm_spp_df$BM2018 != 0 |
                           focal_salm_spp_df$BM2019 != 0),]

merg_spp <- merg_salm_spp_df[which(merg_salm_spp_df$BM2008 != 0 |
                           merg_salm_spp_df$BM2009 != 0 |
                           merg_salm_spp_df$BM2010 != 0 |
                           merg_salm_spp_df$BM2011 != 0 |
                           merg_salm_spp_df$BM2012 != 0 |
                           merg_salm_spp_df$BM2013 != 0 |
                           merg_salm_spp_df$BM2014 != 0 |
                           merg_salm_spp_df$BM2015 != 0 |
                           merg_salm_spp_df$BM2016 != 0 |
                           merg_salm_spp_df$BM2017 != 0 |
                           merg_salm_spp_df$BM2018 != 0 |
                           merg_salm_spp_df$BM2019 != 0),]

# calculate the total stream length in each watershed
StreamLengths <- focal_spp %>% 
  group_by(Waterbody) %>% 
  summarise(Stream_Length = sum(Shape_Leng))

StreamLengths <- StreamLengths[which(is.na(StreamLengths$Waterbody) == FALSE),]

#also need the stream lengths for unnamed water bodies
unNamed <- focal_spp[which(is.na(focal_spp$Waterbody) == TRUE),]
unNamedLengths <- unNamed %>%
  group_by(BM2008) %>%
  summarise(Stream_Length1 = sum(Shape_Leng))

unNamedLengths <- unNamedLengths[which(unNamedLengths$BM2008 != 0),]

focal_spp$Stream_Length0 <- StreamLengths$Stream_Length[match(focal_spp$Waterbody, StreamLengths$Waterbody)]
focal_spp$Stream_Length1 <- unNamedLengths$Stream_Length1[match(focal_spp$BM2008, unNamedLengths$BM2008)]

#there were some unknown waterbodies with no salmon in 2008 so we want to fill in those as well
unNamed2 <- focal_spp[which(is.na(focal_spp$Stream_Length0) == TRUE & is.na(focal_spp$Stream_Length1)== TRUE),]
unNamed2Lengths <- unNamed2 %>%
  group_by(BM2009) %>%
  summarise(Stream_Length2 = sum(Shape_Leng))

unNamed2Lengths <- unNamed2Lengths[which(unNamed2Lengths$BM2009 != 0),]

focal_spp$Stream_Length2 <- unNamed2Lengths$Stream_Length2[match(focal_spp$BM2009, unNamed2Lengths$BM2009)]
merg_spp$Stream_Length0 <- StreamLengths$Stream_Length[match(merg_spp$Waterbody, StreamLengths$Waterbody)]
merg_spp$Stream_Length1 <- unNamedLengths$Stream_Length1[match(merg_spp$BM2008, unNamedLengths$BM2008)]
merg_spp$Stream_Length2 <- unNamed2Lengths$Stream_Length2[match(merg_spp$BM2009, unNamed2Lengths$BM2009)]

focal_spp$Stream_Length <- rowSums(focal_spp[,c("Stream_Length0", "Stream_Length1", "Stream_Length2")], na.rm = TRUE)
merg_spp$Stream_Length <- rowSums(merg_spp[,c("Stream_Length0", "Stream_Length1", "Stream_Length2")], na.rm = TRUE)

focal_spp_lengths <- focal_spp[, c(2:23, 27)]
focal_spp_lengths$AREA_KM2 <- focal_spp_lengths$AREA_SQM/1000000
merg_spp_lengths <- merg_spp[, c(2:18, 20:22, 26)]

#-----------------------------------------------------------------------------------------------
## Calculate diversity in focal watersheds
## Diversity is calculated at the focal and merged watershed scale for completeness but only merged
## is used in final analysis. focal have been commmented out here

## rotate the biomass columns and have a biomass value for each species and year in watershed
#sppFWS_tidy<-select(focal_spp_lengths, which(names(focal_spp_lengths)=="BM2008"):which(names(focal_spp_lengths)=="BM2019"), Waterbody, SpeciesId, WTRSHD_FID, AREA_HA, F_CODE, OBJECTID, AREA_SQM, FEAT_LEN, AREA_KM2, geometry, Shape_Leng, Stream_Length)%>%
  #gather(year, biomass, -Waterbody, -SpeciesId, -WTRSHD_FID, -AREA_HA, -F_CODE, -OBJECTID, -AREA_SQM, -FEAT_LEN, -AREA_KM2, -geometry, -Shape_Leng, -Stream_Length)%>%
  #mutate(year=as.numeric(substr(year, 3, 6)))

sppMWS_tidy<-select(merg_spp_lengths, which(names(merg_spp_lengths)=="BM2008"):which(names(merg_spp_lengths)=="BM2019"), Waterbody, SpeciesId, WSHDFID, AREAKM2, FOCALWSHD, usblAREA, geometry, Shape_Leng, Stream_Length)%>%
  gather(year, biomass, -Waterbody, -SpeciesId, -WSHDFID, -AREAKM2, -FOCALWSHD, -usblAREA, -geometry, -Shape_Leng, -Stream_Length)%>%
  mutate(year=as.numeric(substr(year, 3, 6)))

## calculate proportion of each stream length that is in each watershed
#sppFWS_tidy$ProportioninFWS <- sppFWS_tidy$Shape_Leng/sppFWS_tidy$Stream_Length
sppMWS_tidy$ProportioninMWS <- sppMWS_tidy$Shape_Leng/sppMWS_tidy$Stream_Length

## multiply the salmon biomass by the proportion of the stream length in that watershed
#sppFWS_tidy$biomassinFWS <- sppFWS_tidy$biomass * sppFWS_tidy$ProportioninFWS
sppMWS_tidy$biomassinMWS <- sppMWS_tidy$biomass * sppMWS_tidy$ProportioninMWS

#sppFWS_tidy <- sppFWS_tidy[which(is.na(sppFWS_tidy$biomassinFWS) == FALSE),]
sppMWS_tidy <- sppMWS_tidy[which(is.na(sppMWS_tidy$biomassinMWS) == FALSE),]

## sum biomass across streams within a watershed and year for each species
#bm_sppFWS <- sppFWS_tidy %>%
  #group_by(WTRSHD_FID, year, SpeciesId) %>%
  #summarise(YearlyBMinFWS=sum(biomassinFWS))

bm_sppMWS <- sppMWS_tidy %>%
  group_by(WSHDFID, year, SpeciesId) %>%
  summarise(YearlyBMinMWS=sum(biomassinMWS))

#bm_FWS <- bm_sppFWS %>% 
  #group_by(WTRSHD_FID, year) %>%
  #summarise(totalYrlyBMinFWS = sum(YearlyBMinFWS))

bm_MWS <- bm_sppMWS %>% 
  group_by(WSHDFID, year) %>%
  summarise(totalYrlyBMinMWS = sum(YearlyBMinMWS))

#bm_sppFWS <- merge(bm_sppFWS, bm_FWS, by = c("WTRSHD_FID", "year"))
#bm_sppFWS$proportion <- bm_sppFWS$YearlyBMinFWS / bm_sppFWS$totalYrlyBMinFWS

bm_sppMWS <- merge(bm_sppMWS, bm_MWS, by = c("WSHDFID", "year"))
bm_sppMWS$proportion <- bm_sppMWS$YearlyBMinMWS / bm_sppMWS$totalYrlyBMinMWS

library(vegan)

#bm_diversityFWS <- bm_sppFWS %>%
#  group_by(WTRSHD_FID, year) %>%
#  mutate(Shannon = diversity(x = YearlyBMinFWS, index = 'shannon'))

bm_diversityMWS <- bm_sppMWS %>%
  group_by(WSHDFID, year) %>%
  mutate(Shannon = diversity(x = YearlyBMinMWS, index = 'shannon'))

#bm_diversityFWS <- bm_diversityFWS[,c(1,2,7)]
bm_diversityMWS <- bm_diversityMWS[,c(1,2,7)]

#bm_diversityFWS <- unique(bm_diversityFWS)
bm_diversityMWS <- unique(bm_diversityMWS)

## now have a shannon diversity index for each focal watershed and year
#write.csv(bm_diversityFWS, "C:/Users/monic/Documents/Documents/PhD/focal_spp_diversity_per_year.csv")

#have a shannon diversity index for each larger watershed and year
write.csv(bm_diversityMWS, here("data", "processed", "merg_spp_diversity_per_year.csv"))

          