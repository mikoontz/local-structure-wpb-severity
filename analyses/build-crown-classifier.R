# This script will build a classifier to predict the categorical response of 
# tree species from the reflectance data within each crown polygon.
# The classifier could then be used to classify all the trees in a given scene

# First step is to hand-classify a bunch of trees within a few scenes. We'll
# then build a model to predict the hand-classified species as a response
# using the reflectance data within that crown as predictors.

# We'll use the merged sites as a starting point, because they also include the
# more spatially resolved imagery from the X3 camera. That might make it easier
# to determine the species from the air.

# RandomForest
# Support Vector Machines

library(sf)
library(tidyverse)
library(purrr)
library(raster)
library(tictoc)
library(here)
library(caret)

source(here::here("data/data_carpentry/extract-reflectance-within-crowns-from-mosaics.R"))

# These sites had X3 and RedEdge photos merged into the same project, so we look in a different place for some of the relevant files. These are also the ones we'll start with for hand classifying
merged_sites <- c("eldo_3k_2",
                  "eldo_3k_3",
                  "eldo_4k_2")

other_sites_to_hand_classify <-
  c("eldo_3k_1", "sequ_4k_1", "sequ_5k_1", "sequ_6k_2", "stan_3k_1", "stan_4k_1", "stan_5k_1", "sier_3k_1", "sier_4k_1", "sier_5k_1", 
    "eldo_4k_1", "eldo_5k_2", "eldo_5k_3", "sequ_5k_2", "sequ_6k_3", "sier_5k_3", "stan_5k_2")

unusable_sites <- c("eldo_4k_3", # too many blocks
                    "stan_4k_3", # too many blocks
                    "stan_5k_3", # too many blocks
                    "sequ_4k_2") # middle section flown on a separate day and the stitch looks terrible

sites_to_hand_classify <-
  c(merged_sites, other_sites_to_hand_classify)


# Copy segmented crowns shapefile to hand-classified directory ------------

sites_to_hand_classify %>%
  walk(.f = function(current_site) {
    
    if(!dir.exists(here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns")))) {
      dir.create(here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns")))
    } 
    
    if(!file.exists(here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns/", current_site, "_hand-classified-crowns.shp")))) {  
      crowns <- sf::st_read(dsn = here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_crowns/", current_site, "_crowns.shp")))
      
      sf::st_write(obj = crowns, dsn = here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns/", current_site, "_hand-classified-crowns.shp")))
    } else (print(paste0("Shapefile for hand classified crowns of ", current_site, " already exists!")))
    
  })


# Pause to go hand classify the crowns ------------------------------------

### Pause!! Go use QGIS to overlay these crowns shapefiles on top of the index mosaic
### and add two extra attributes: live (1/0) and species (pipo/pila/cade/abco). Use the
### known tree locations from the ground plot to pick some obvious examples that fit
### into these categories and manually add the appropriate data to the new attributes
### for a couple hundred trees.

# Extract the reflectance data from within hand-classified crowns ---------
tic()
crowns_with_reflectance <-
  map(.x = sites_to_hand_classify, .f = function(current_site) {
    
    # target_path <- here::here(paste0("data/data_output/classified/model-classified/crowns-with-reflectance/", current_site, "_crowns-with-reflectance/", current_site, "_crowns-with-reflectance.shp"))
    crowns_path <- here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns/", current_site, "_hand-classified-crowns.shp"))
    ttops_path <- here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_ttops/", current_site, "_ttops.shp"))
    crowns <- sf::st_read(crowns_path) %>% filter(!is.na(live) | !is.na(species))
    ttops <- sf::st_read(ttops_path)
    
    index_path <- here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_index.tif"))
    index <- velox::velox(index_path)
    index_copy <- index$copy()
    
    current_crowns <- extract_reflectance_from_crowns(index = index_copy,
                                                      crowns = crowns,
                                                      ttops = ttops)
    
    current_crowns <- 
      current_crowns %>%
      dplyr::select(treeID, height, ch_area, live, species, x, y, b_mean, g_mean, r_mean, re_mean, nir_mean, ndvi_mean, rgi_mean, gbi_mean, ndre_mean) %>% 
      dplyr::mutate(treeID = paste(current_site, treeID, sep = "_"),
                    crs = st_crs(.)$proj4string) %>% 
      st_drop_geometry()
    
    return(current_crowns)
    
  }) %>% 
  do.call("rbind", .)
toc()
  
  glimpse(crowns_with_reflectance)
  write_csv(x = crowns_with_reflectance, path = here::here("data/data_output/classified/hand-classified/hand-classified-crowns.csv"))


# Use the hand classified data to build supervised classifier -------------

crowns_with_reflectance <- readr::read_csv(here::here("data/data_output/classified/hand-classified/hand-classified-crowns.csv"))

live_or_dead_idx <- caret::createDataPartition(crowns_with_reflectance$live, p = 0.8, list = FALSE)

live_or_dead_training <- crowns_with_reflectance[live_or_dead_idx, ]
live_or_dead_testing <- crowns_with_reflectance[-live_or_dead_idx, ]

live_or_dead_fit <- train(
  as.factor(live) ~ b_mean + g_mean + r_mean + re_mean + nir_mean + ndvi_mean + rgi_mean + gbi_mean + ndre_mean,
  data = live_or_dead_training,
  method = "LogitBoost"
)

live_or_dead_fit

live_crowns <- 
  crowns_with_reflectance %>% 
  dplyr::filter(live == 1) %>% 
  dplyr::mutate(functional_group = case_when(species == "pila" ~ "pinus",
                                             species == "pipo" ~ "pinus",
                                             TRUE ~ species))


species_idx <-
  live_crowns %>%
  pull(species) %>% 
  caret::createDataPartition(p = 0.8, list = FALSE)

species_training <- live_crowns[species_idx, ]
species_testing <- live_crowns[-species_idx, ]


# collapsing PIPO and PILA into Pinus -------------------------------------


functional_group_idx <-
  live_crowns %>%
  pull(functional_group) %>% 
  caret::createDataPartition(p = 0.8, list = FALSE)

functional_group_training <- live_crowns[functional_group_idx, ]
functional_group_testing <- live_crowns[-functional_group_idx, ]


# Some different classification approaches --------------------------------


ldaFit <- train(
  species ~ b_mean + g_mean + r_mean + re_mean + nir_mean + ndvi_mean + rgi_mean + gbi_mean + ndre_mean,
  data = species_training,
  method = "lda"
)

ldaFit
ldaProbs <- predict(ldaFit, newdata = species_testing, type = "prob")
head(ldaProbs)
species_testing

# ldaFit_library <- train(
#   species ~ b_mean + g_mean + r_mean + re_mean + nir_mean + ndvi_mean + rgi_mean + gbi_mean + ndre_mean,
#   data = spectral_library %>% filter(source == "queally"),
#   method = "lda"
# )
# 
# 
# lda_library_Probs <- predict(ldaFit_library, newdata = species_testing, type = "prob")
# test <-
#   species_testing %>% 
#   dplyr::mutate(predicted_species = predict(ldaFit_library, newdata = species_testing)) %>% 
#   dplyr::select(species, predicted_species)
# head(test, 20)


rfFit <- train(
  species ~ b_mean + g_mean + r_mean + re_mean + nir_mean + ndvi_mean + rgi_mean + gbi_mean + ndre_mean,
  data = species_training,
  method = "rf"
)

rfFit

rfProbs <- predict(rfFit, newdata = species_testing, type = "prob")
head(rfProbs)
species_testing



cforestFit <- train(
  species ~ b_mean + g_mean + r_mean + re_mean + nir_mean + ndvi_mean + rgi_mean + gbi_mean + ndre_mean,
  data = species_training,
  method = "cforest"
)

cforestFit
cforestProbs <- predict(cforestFit, newdata = species_testing, type = "prob")
head(cforestProbs)
species_testing


live_crowns %>% 
  count(species)

ggplot(plsFit)

plsClasses <- predict(plsFit, newdata = testing)
str(plsClasses)
#>  Factor w/ 2 levels "M","R": 2 1 1 2 1 2 2 2 2 2 ...
plsProbs <- predict(plsFit, newdata = testing, type = "prob")
head(plsProbs)


# RedEdge3 band information
# http://www.leptron.com/manuals/RedEdge_User_Manual.pdf

# Spectral library information --------------------------------------------
# NASA HyspIRI Airborne Campaign Leaf and Canopy Spectra and Leaf Traits
# https://ecosis.org/#result/dd94f09c-1794-44f4-82e9-a7ca707a1ec0
# nasa-hyspiri-airborne-campaign-leaf-and-canopy-spectra-and-leaf-traits


# Susan Meerdink
# fresh-leaf-spectra-to-estimate-leaf-traits-for-california-ecosystems.csv
# https://www.sciencedirect.com/science/article/pii/S0034425716303066?via%3Dihub
# https://ecosis.org/#result/0fadcc45-f79e-4fd3-a6ca-8afaf26ae299

# california_species_data.csv
# Natalie Queally
# https://ecosis.org/#result/c1a5b651-9c46-4e06-a07f-38a2e7b4faf4

micasense_rededge <- 
  data_frame(band = 1:5,
             band_name = c("b", "g", "r", "nir", "re"),
             center_wavelength = c(475, 560, 668, 840, 717),
             bandwidth = c(20, 20, 10, 40, 10)) %>% 
  dplyr::mutate(start_wavelength = center_wavelength - (bandwidth / 2),
                end_wavelength = center_wavelength + (bandwidth / 2))

species_present <- c("abco", "cade", "pipo", "pila", "quke")

serbin <- read_csv(here::here("data/data_raw/nasa-hyspiri-airborne-campaign-leaf-and-canopy-spectra-and-leaf-traits.csv")) %>% 
  dplyr::select(`Latin Genus`, `Latin Species`, Leaf_or_Canopy, Sample_Name, 21:ncol(.)) %>% 
  dplyr::filter(Leaf_or_Canopy == "leaf") %>% 
  dplyr::mutate(genus = substr(`Latin Genus`, start = 1, stop = 2), species = substr(`Latin Species`, start = 1, stop = 2)) %>% 
  dplyr::mutate(species = tolower(paste0(genus, species))) %>% 
  dplyr::select(-`Latin Genus`, -`Latin Species`, -Leaf_or_Canopy) %>% 
  dplyr::rename(replicate = Sample_Name) %>% 
  dplyr::select(species, replicate, everything()) %>% 
  dplyr::mutate(replicate = paste("serbin", species, replicate, sep = "_")) %>% 
  dplyr::filter(species %in% species_present) %>% 
  tidyr::gather(key = "wavelength", value = "reflectance", -(1:2)) %>% 
  dplyr::mutate(wavelength = as.numeric(wavelength)) %>% 
  dplyr::mutate(reflectance = as.numeric(reflectance)) %>% 
  dplyr::mutate(band_name = case_when(  (wavelength >= micasense_rededge$start_wavelength[1] & wavelength <= micasense_rededge$end_wavelength[1]) ~ "b",
                                        (wavelength >= micasense_rededge$start_wavelength[2] & wavelength <= micasense_rededge$end_wavelength[2]) ~ "g",
                                        (wavelength >= micasense_rededge$start_wavelength[3] & wavelength <= micasense_rededge$end_wavelength[3]) ~ "r",
                                        (wavelength >= micasense_rededge$start_wavelength[4] & wavelength <= micasense_rededge$end_wavelength[4]) ~ "nir",
                                        (wavelength >= micasense_rededge$start_wavelength[5] & wavelength <= micasense_rededge$end_wavelength[5]) ~ "re")) %>% 
  dplyr::filter(!is.na(band_name)) %>% 
  dplyr::group_by(species, replicate, band_name) %>% 
  dplyr::summarize(reflectance = mean(reflectance)) %>% 
  dplyr::mutate(source = "serbin")

    
# meerdink <- read_csv(here::here("data/data_raw/fresh-leaf-spectra-to-estimate-leaf-traits-for-california-ecosystems.csv")) %>% 
#   dplyr::mutate(species = tolower(species)) %>% 
#   dplyr::mutate(age = ifelse(str_detect(`sample name`, pattern = "Old"), yes = "old", no = "young")) %>% 
#   dplyr::filter(age == "young") %>% 
#   dplyr::mutate(replicate = paste("meerdink", species, Replicate, sep = "_")) %>% 
#   dplyr::filter(species %in% species_present) %>% 
#   dplyr::select(species, replicate, 35:ncol(.)) %>% 
#   tidyr::gather(key = "wavelength", value = "reflectance", -(1:2)) %>% 
#   dplyr::mutate(wavelength = as.numeric(wavelength)) %>% 
#   dplyr::mutate(band_name = case_when(  (wavelength >= micasense_rededge$start_wavelength[1] & wavelength <= micasense_rededge$end_wavelength[1]) ~ "b",
#                                         (wavelength >= micasense_rededge$start_wavelength[2] & wavelength <= micasense_rededge$end_wavelength[2]) ~ "g",
#                                         (wavelength >= micasense_rededge$start_wavelength[3] & wavelength <= micasense_rededge$end_wavelength[3]) ~ "r",
#                                         (wavelength >= micasense_rededge$start_wavelength[4] & wavelength <= micasense_rededge$end_wavelength[4]) ~ "nir",
#                                         (wavelength >= micasense_rededge$start_wavelength[5] & wavelength <= micasense_rededge$end_wavelength[5]) ~ "re")) %>% 
#   dplyr::filter(!is.na(band_name)) %>% 
#   dplyr::group_by(species, replicate, band_name) %>% 
#   dplyr::summarize(reflectance = mean(reflectance)) %>% 
#   dplyr::mutate(source = "meerdink")
  
queally <- read_csv(here::here("data/data_raw/california_species_data.csv")) %>% 
  dplyr::mutate(species = USDA_symbol) %>% 
  dplyr::filter(species %in% species_present) %>% 
  dplyr::mutate(replicate = paste("queally", species, 1:nrow(.), sep = "_")) %>% 
  dplyr::select(species, replicate, 23:ncol(.)) %>% 
  tidyr::gather(key = "wavelength", value = "reflectance", -(1:2)) %>% 
  dplyr::mutate(wavelength = as.numeric(wavelength)) %>% 
  dplyr::mutate(band_name = case_when(  (wavelength >= micasense_rededge$start_wavelength[1] & wavelength <= micasense_rededge$end_wavelength[1]) ~ "b",
                                        (wavelength >= micasense_rededge$start_wavelength[2] & wavelength <= micasense_rededge$end_wavelength[2]) ~ "g",
                                        (wavelength >= micasense_rededge$start_wavelength[3] & wavelength <= micasense_rededge$end_wavelength[3]) ~ "r",
                                        (wavelength >= micasense_rededge$start_wavelength[4] & wavelength <= micasense_rededge$end_wavelength[4]) ~ "nir",
                                        (wavelength >= micasense_rededge$start_wavelength[5] & wavelength <= micasense_rededge$end_wavelength[5]) ~ "re")) %>% 
  dplyr::filter(!is.na(band_name)) %>% 
  dplyr::group_by(species, replicate, band_name) %>% 
  dplyr::summarize(reflectance = mean(reflectance)) %>% 
  dplyr::mutate(source = "queally")

spectral_library <-
  rbind(serbin %>% filter(species %in% c("pila", "cade")), queally)

spectral_library <-
  spectral_library %>% 
  spread(key = band_name, value = reflectance) %>% 
  ungroup()

colnames(spectral_library)[4:8] <- paste(colnames(spectral_library)[4:8], "mean", sep = "_")

spectral_library <-
  spectral_library %>% 
  dplyr::mutate(ndvi_mean = (nir_mean - r_mean) / (nir_mean + r_mean),
                rgi_mean = r_mean / g_mean,
                gbi_mean = g_mean / b_mean,
                ndre_mean = (re_mean - r_mean) / (re_mean + r_mean))


spectral_library %>% dplyr::select(-replicate) %>% group_by(species, source) %>% 
  summarize_all(funs(mean))



