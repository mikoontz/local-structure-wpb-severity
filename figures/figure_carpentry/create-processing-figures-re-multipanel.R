# Purpose: create a multi-panel processing figure that includes:

# Pix4D outputs (DSM; Point Cloud)
# DTM
# CHM
# Tree detection
# Segmentation
# Classification (live/dead; species)
# Two columns, 4 rows I think.

library(tidyverse)
library(raster)
library(viridis)
library(lidR)

source("data/data_carpentry/make-processing-checklist.R")
example_site <- 
  sites_checklist %>% 
  dplyr::filter(rasterized_trees_check) %>% 
  slice(1) %>% 
  pull(site)
"data/data_output/site_data/eldo_3k_1/"
dsm <- raster::raster(paste0("data/data_output/site_data/", example_site, "))


dtm <- raster::raster(paste0("data/data_output/site_data/", site, "/", site, "_2m-dtm.tif"))

site_bounds <- 
  sf::st_read(paste0("data/data_output/site_data/", site, "/", site, "_mission-footprint/", site, "_site-bounds.geoJSON")) %>% 
  sf::st_transform(proj4string(dsm))

dsm <- 
  dsm %>% 
  raster::crop(site_bounds)

dtm <-
  dtm %>% 
  raster::crop(site_bounds) %>% 
  raster::resample(dsm)

chm <- dsm - dtm
chm <-
  raster::focal(chm, w = matrix(1, 3, 3), mean)
chm[chm < 0] <- 0

ttops <- sf::st_read(paste0("data/data_output/site_data/", site, "/", site, "_ttops/", site, "_ttops.shp"))
crowns <- sf::st_read(paste0("data/data_output/site_data/", site, "/", site, "_crowns/", site, "_crowns.shp"))

cc <- sf::st_read(paste0("data/data_output/classified/model-classified/crown-shapefiles/", site, "_classified-crowns/", site, "_classified-crowns.shp"))
