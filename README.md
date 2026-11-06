Manuscript Title: Increases in salmon biomass consumption by grizzly bears and its ecological implications

**Description of the data and code files**

Run the scripts in numerical order. Scripts 05 to 07 run from the processed data files provided. Hair-snag locations of individual bears are sensitive and are not shared. The scripts that use them (00_sensitive.R and 04a_gb_to_watershed.R) were run privately and are not included; the location-free files they produced are provided below.

01_SalmonBiomassDensity.R: calculates annual salmon biomass density per merged watershed (2008-2019); reproduces methods described in section '2.31 Salmon biomass density' in the manuscript

02_SalmonDiversity_by_watershed.R: calculates annual Shannon diversity (and the Hill number) of salmon species per merged watershed (2008-2019); reproduces methods described in section '2.32 Salmon diversity' in the manuscript

03_EVI_by_Watershed.R: extracts mean spring-summer EVI from MODIS imagery per merged watershed (2008-2018); reproduces methods described in section '2.33 Terrestrial plant productivity (vegetation)' in the manuscript

04_attribute_covariotes.R: attributes all covariates to each bear-year with a one-year lag, calculates relative age and scales the predictors; reproduces methods described in section '2.3 Covariate attribution' in the manuscript

05_glmm_analysis.R: fits and compares beta regression GLMMs, with model diagnostics; reproduces analyses described in section '2.4 Data analysis' and results in section '3' in the manuscript

06_propagations_and_bioenergetics.R: propagates MixSIAR uncertainty through the GLMM and estimates absolute salmon consumption (kg) as bears age using m1 model predictions and Kingsley et al. (1983) growth curves; reproduces analyses described in section '2.41 Estimating increase in salmon biomass consumption with age' in the manuscript

07_SI_full_prop_realizations.R: sensitivity analysis adding uncertainty in age at first detection, the growth curve and the intake equation; reproduces analyses described in [SI section] of the Supplementary Information

gbs_in_wshds_no_location.csv: hair-snag detections with their focal watershed and MixSIAR summaries of the salmon proportion; coordinates removed (see section '2.3 Covariate attribution' in the manuscript)

individual_year_sex_lookup.csv: sex of each bear-year (see section '2.3 Covariate attribution' in the manuscript)

all_cov_scaled_no_locations.csv: final modelling dataset (550 individual-years) with scaled covariates; coordinates removed (see sections '2.3 Covariate attribution' and '2.4 Data analysis' in the manuscript)

gb_all_posteriors_no_location_2026-10-20.csv: posterior draws of the dietary salmon proportion for each bear-year from the MixSIAR model; coordinates removed (see [MixSIAR methods section] in the manuscript)

Deacy_Robbins_intake_data.csv: intake data used for the intake equation in script 07 [Charlie Robbins, unpublished]

Data for scripts 01-03 [salmon escapement, watershed layers, MODIS imagery]: 
NUSEDS: 31.	Fisheries and Oceans Canada. 2019 NuSEDS regional adult salmon escapement database 1950–2010 [data set]. Nanaimo, BC: Pacific Biological Station, Fisheries and Oceans Canada.<img width="468" height="69" alt="image" src="https://github.com/user-attachments/assets/3df20ed5-e9a6-492f-83ce-01302c7e4c35" />

Watershed Layers: 25.	Government of British Columbia. 2011 Freshwater Atlas assessment watersheds [data set]. BC Data Catalogue, GeoBC, Ministry of Water, Land and Resource Stewardship. See https://catalogue.data.gov.bc.ca/dataset/freshwater-atlas-assessment-watersheds (accessed 20 April 2024).<img width="468" height="85" alt="image" src="https://github.com/user-attachments/assets/d3803252-686d-4168-aec6-0fd9b9ce1f72" />

MODIS IMAGERY:	Didan K. 2021 MODIS/Terra vegetation indices 16-day L3 global 250m SIN grid V061 [data set]. NASA Land Processes Distributed Active Archive Center. (doi:10.5067/MODIS/MOD13Q1.061) (accessed 7 February 2025).<img width="468" height="69" alt="image" src="https://github.com/user-attachments/assets/31dbea1b-ecfb-4c24-87fa-a6a8d40b9453" />


## Data folder layout

Zenodo stores files without folders. After downloading, recreate the
structure below inside the project folder so the scripts can find each file.
All paths in the scripts are relative to the project root via `here()`.

```
data/
├── raw/
│   ├── Deacy_Robbins_intake_data.csv
│   ├── FWA_ASSESSMENT_WATERSHEDS_POLY.zip
│   ├── gb_all_posteriors_no_location.csv
│   └── RAW_EVI_APRIL22_JULY27.zip
└── processed/
    ├── all_cov_scaled_no_locations_FINAL_20261104.csv
    ├── combinedStreams2019.zip
    ├── combinedStreamsBMTotal_NoZeros2019.zip
    ├── covariotes_no_loactions_6-11-04.csv
    ├── Focal Watersheds in Study Area.zip
    ├── gbs_in_wshds_no_loaactions_6-11-04.csv
    ├── individual_year_sex_lookup.csv
    ├── LRG_Watersheds_for_all_sites.zip
    ├── merg_salm_biomass_per_year_w_density.csv
    ├── merg_spp_diversity_per_year_with_Hill.csv
    └── sitedata_salmon_EVI_no_location_2026-11-04.csv
```

**raw/** holds source data as received, before processing. 
**processed/** holds files created or cleaned by the analysis scripts 
(or prepared by co-authors), which later scripts read in.


**Code/Software**

Analyses were conducted in R. R Core Team (2026). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/. Version [R version, e.g. 4.5.1]

Packages: here, dplyr, tidyr, ggplot2, patchwork, sf, terra, data.table, glmmTMB, bbmle, MuMIn, DHARMa, performance, ggeffects, MASS. Versions are listed by sessionInfo() at the end of each script.
