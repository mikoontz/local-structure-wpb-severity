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
library(lidR)
library(viridis)
library(raster)
library(ForestTools)
library(gstat)
library(nngeo)
# devtools::install_github("Jean-Romain/lidRplugins")
library(lidRplugins)
library(furrr)
library(stars)
library(tictoc)
library(here)
library(units)

source(here::here("data/data_carpentry/make-processing-checklist.R"))
source(here::here("data/data_carpentry/extract-reflectance-within-crowns-from-mosaics.R"))

sites_checklist
# These sites had X3 and RedEdge photos merged into the same project, so we look in a different place for some of the relevant
# files.
merged_sites <- c("eldo_3k_2",
                  "eldo_3k_3",
                  "eldo_4k_2")

other_sites_to_hand_classify <-
  c("eldo_3k_1")

sites_to_hand_classify <-
  c(merged_sites, other_sites_to_hand_classify)

tic()
index_mosaics_to_hand_classify <-
  sites_to_hand_classify %>%
  map(.f = function(current_site) {
    crowns <- sf::st_read(dsn = here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_crowns/", current_site, "_crowns.shp")))
    
    if(!dir.exists(here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns")))) {
      dir.create(here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns")))
    } 
    
    if(!file.exists(here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns/", current_site, "_hand-classified-crowns.shp")))) {  
      
      sf::st_write(obj = crowns, dsn = here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns/", current_site, "_hand-classified-crowns.shp")))
    } else (print(paste0("Shapefile for hand classified crowns of ", current_site, " already exists!")))
    
    index_path <- here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_index.tif"))
    
    index <- velox::velox(index_path)
    
    return(index)
    
  })
toc()

### Pause!! Go use QGIS to overlay these crowns shapefiles on top of the index mosaic
### and add two extra attributes: live (1/0) and species (pipo/pila/cade/abco). Use the
### known tree locations from the ground plot to pick some obvious examples that fit
### into these categories and manually add the appropriate data to the new attributes
### for a couple hundred trees.

# pick one to use as an example for testing
# current_site <- "eldo_3k_1"

tic()
crowns_with_reflectance <-
  map2(.x = sites_to_hand_classify, .y = index_mosaics_to_hand_classify, .f = function(current_site, index) {
    
    # target_path <- here::here(paste0("data/data_output/classified/model-classified/crowns-with-reflectance/", current_site, "_crowns-with-reflectance/", current_site, "_crowns-with-reflectance.shp"))
    crowns_path <- here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns/", current_site, "_hand-classified-crowns.shp"))
    ttops_path <- here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_ttops/", current_site, "_ttops.shp"))
    crowns <- sf::st_read(crowns_path) %>% filter(!is.na(live) | !is.na(species))
    ttops <- sf::st_read(ttops_path)
    
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
    
  }) %>% 
  do.call("rbind", .)
toc()

glimpse(crowns_with_reflectance)
