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
library(here)
library(caret)
library(caTools)
library(klaR)

source(here::here("workflow/21a_extract-reflectance-within-crowns-from-mosaics-function.R"))

sites_to_hand_classify <-
  c("eldo_3k_1", "eldo_3k_2", "eldo_3k_3", "eldo_4k_1", "eldo_4k_2", "eldo_5k_2", "eldo_5k_3",
    "stan_3k_1", "stan_4k_1", "stan_5k_1", "stan_5k_2",
    "sier_3k_1", "sier_4k_1", "sier_5k_1", "sier_5k_3",
    "sequ_4k_1", "sequ_5k_1", "sequ_5k_2", "sequ_6k_2", "sequ_6k_3")

# Copy segmented crowns shapefile to hand-classified directory ------------

if(!dir.exists("data/data_drone/L3b/hand-classified-trees")) {
  dir.create("data/data_drone/L3b/hand-classified-trees", recursive = TRUE)
}

sites_to_hand_classify %>%
  walk(.f = function(current_site) {
    
    if(!file.exists(here::here(paste0("data/data_drone/L3b/hand-classified-trees/", current_site, "_hand-classified-trees.gpkg")))) {  
      
      file.copy(from = here::here(paste0("data/data_drone/L3a/crowns/", current_site, "_crowns.gpkg")),
                to = here::here(paste0("data/data_drone/L3b/hand-classified-trees/", current_site, "_hand-classified-trees.gpkg")))
      
    } else (print(paste0("GPKG for hand classified trees of ", current_site, " already exists!")))
    
  })

# Pause to go hand classify the trees ------------------------------------

### Pause!! Go use QGIS to overlay these crowns shapefiles on top of the index mosaic
### and add two extra attributes: live (1/0) and species (pipo/pila/cade/abco). Use the
### known tree locations from the ground plot to pick some obvious examples that fit
### into these categories and manually add the appropriate data to the new attributes
### for a couple hundred trees.

# Extract the reflectance data from within hand-classified crowns ---------
if(!file.exists(here::here("data/data_drone/L3b/hand-classified-trees.csv"))) {
  
  crowns_with_reflectance <-
    map(.x = sites_to_hand_classify, .f = function(current_site) {
      
      ttops <- sf::st_read(here::here(paste0("data/data_drone/L3a/ttops/", current_site, "_ttops.gpkg")))

      crowns <- 
        sf::st_read(here::here(paste0("data/data_drone/L3b/hand-classified-trees/", current_site, "_hand-classified-crowns.gpkg"))) %>% 
        filter(!is.na(live) | !is.na(species))
      
      index <- velox::velox(here::here(paste0("data/data_drone/L2/index/", current_site, "_index.tif")))
      index_copy <- index$copy()
      
      current_crowns <- extract_reflectance_from_crowns(index = index_copy,
                                                        crowns = crowns,
                                                        ttops = ttops)
      
      current_crowns <- 
        current_crowns %>%
        dplyr::select(treeID, height, ch_area, live, species, x, y, b_mean, g_mean, r_mean, re_mean, nir_mean, ndvi_mean, rgi_mean, cire_mean, cig_mean, ndre_mean) %>% 
        dplyr::mutate(treeID = paste(current_site, treeID, sep = "_"),
                      crs = st_crs(.)$proj4string) %>% 
        st_drop_geometry()
      
      return(current_crowns)
      
    }) %>% 
    do.call("rbind", .)

  write_csv(x = crowns_with_reflectance, path = here::here("data/data_drone/L3b/hand-classified-trees_all.csv"))
}


# read in the hand-classified crowns --------------------------------------

crowns_with_reflectance <- readr::read_csv(here::here("data/data_drone/L3b/hand-classified-trees_all.csv"))

# create the live/dead classifier -----------------------------------------
set.seed(1409)
live_or_dead_idx <- caret::createDataPartition(crowns_with_reflectance$live, p = 0.8, list = FALSE)

live_or_dead_training <- crowns_with_reflectance[live_or_dead_idx, ]
live_or_dead_testing <- crowns_with_reflectance[-live_or_dead_idx, ]

live_or_dead_classifier <- caret::train(
  as.factor(live) ~ b_mean + g_mean + r_mean + re_mean + nir_mean + ndvi_mean + rgi_mean + cire_mean + cig_mean + ndre_mean,
  data = live_or_dead_training,
  method = "LogitBoost"
)

# Classifier is called `live_or_dead_classifier`
live_or_dead_classifier


# subset to only the live crowns ------------------------------------------

live_crowns <- 
  crowns_with_reflectance %>% 
  dplyr::filter(live == 1) %>% 
  dplyr::mutate(functional_group = case_when(species == "pila" ~ "pinus",
                                             species == "pipo" ~ "pinus",
                                             TRUE ~ species))


# create the species classifier -------------------------------------------

species_idx <-
  live_crowns %>%
  pull(species) %>% 
  caret::createDataPartition(p = 0.8, list = FALSE)

species_training <- live_crowns[species_idx, ]
species_testing <- live_crowns[-species_idx, ]

# Best classification method after tests (see docs/interim-reports/)
species_classifier <- caret::train(
  species ~ b_mean + g_mean + r_mean + re_mean + nir_mean + ndvi_mean + rgi_mean + cire_mean + cig_mean + ndre_mean,
  data = species_training,
  method = "rda",
  preProcess = c("center", "scale"),
  tuneGrid = expand.grid(gamma = 0, lambda = c(seq(0.2, 1.0, by = 0.1)))
)

# species classifier is called 'species_classifier
species_classifier


if(!dir.exists("data/data_drone/L3b/classifier-models")) {
  dir.create("data/data_drone/L3b/classifier-models", recursive = TRUE)
}

readr::write_rds(x = live_or_dead_classifier, path = "data/data_drone/L3b/classifier-models/live-or-dead-classifier.rds")
readr::write_rds(x = species_classifier, path = "data/data_drone/L3b/classifier-models/species-classifier.rds")
