library(sf)
library(raster)
library(tabularaster)
library(tidyverse)
library(viridis)
library(purrr)
library(velox)
library(viridis)

# Source in the classifier code.
# The live/dead classifier is an R object called `live_or_dead_fit`
# The species classifier is an R object called `rdafit`

source(here::here("analyses/build-crown-classifier.R"))
live_or_dead_fit
rdaFit

# crowns_with_reflectance <- readr::read_csv(file = here::here("data/data_output/classified/model-classified/crowns-with-reflectance_all.csv"))
crowns_with_reflectance <- readr::read_csv(file = here::here("data/data_output/classified/model-classified/crowns-with-reflectance_35m-buffer.csv"))

classified_trees  <-
  crowns_with_reflectance %>% 
  dplyr::mutate(live = live_or_dead_fit$levels[predict(live_or_dead_fit, newdata = .)]) %>% 
  dplyr::mutate(live = as.numeric(as.character(live))) %>%
  dplyr::mutate(species = ifelse(live == 1, yes = rdaFit$levels[predict(rdaFit, newdata = .)], no = NA))

classified_trees_3310 <-
  classified_trees %>% 
  split(f = .$crs) %>% 
  lapply(FUN = function(trees) {
    current_crs <- unique(trees$crs)
    
    trees3310 <-
      trees %>% 
      st_as_sf(coords = c("x", "y"), crs = current_crs, remove = TRUE) %>% 
      st_transform(3310) %>% 
      tidyr::separate(col = treeID, into = c("forest", "elev", "rep", "id"), sep = "_", remove = FALSE) %>% 
      dplyr::mutate(site = paste(forest, elev, rep, sep = "_")) %>% 
      dplyr::select(-crs, -id) %>% 
      dplyr::select(treeID, site, forest, elev, rep, live, species, everything())
  }) %>% 
  do.call("rbind", .)


sf::st_write(classified_trees_3310, dsn = here::here("analyses/analyses_output/classified-trees/classified-trees.shp"))
