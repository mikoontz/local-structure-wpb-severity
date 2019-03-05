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

source(here::here("data/data_carpentry/extract-reflectance-within-crowns-from-mosaics-function.R"))

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
if(!file.exists(here::here("data/data_output/classified/hand-classified/hand-classified-crowns.csv"))) {
  
  tic()
  crowns_with_reflectance <-
    map(.x = sites_to_hand_classify, .f = function(current_site) {
      
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
        dplyr::select(treeID, height, ch_area, live, species, x, y, b_mean, g_mean, r_mean, re_mean, nir_mean, ndvi_mean, rgi_mean, cire_mean, cig_mean, ndre_mean) %>% 
        dplyr::mutate(treeID = paste(current_site, treeID, sep = "_"),
                      crs = st_crs(.)$proj4string) %>% 
        st_drop_geometry()
      
      return(current_crowns)
      
    }) %>% 
    do.call("rbind", .)
  toc()
  
  glimpse(crowns_with_reflectance)
  write_csv(x = crowns_with_reflectance, path = here::here("data/data_output/classified/hand-classified/hand-classified-crowns.csv"))
}


# read in the hand-classified crowns --------------------------------------

crowns_with_reflectance <- readr::read_csv(here::here("data/data_output/classified/hand-classified/hand-classified-crowns.csv"))


# create the live/dead classifier -----------------------------------------
live_or_dead_idx <- caret::createDataPartition(crowns_with_reflectance$live, p = 0.8, list = FALSE)

live_or_dead_training <- crowns_with_reflectance[live_or_dead_idx, ]
live_or_dead_testing <- crowns_with_reflectance[-live_or_dead_idx, ]

live_or_dead_fit <- train(
  as.factor(live) ~ b_mean + g_mean + r_mean + re_mean + nir_mean + ndvi_mean + rgi_mean + cire_mean + cig_mean + ndre_mean,
  data = live_or_dead_training,
  method = "LogitBoost"
)

# Classifier is called `live_or_dead_fit`
live_or_dead_fit


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
rdaFit <- train(
  species ~ b_mean + g_mean + r_mean + re_mean + nir_mean + ndvi_mean + rgi_mean + cire_mean + cig_mean + ndre_mean,
  data = species_training,
  method = "rda",
  preProcess = c("center", "scale"),
  tuneGrid = expand.grid(gamma = 0, lambda = c(seq(0.2, 1.0, by = 0.1)))
)

# species classifier is called rdaFit
rdaFit

