library(sf)
library(raster)
library(tidyverse)
library(viridis)
library(purrr)
library(viridis)
library(modelr)

# The live/dead classifier is an R object called `live_or_dead_classifier`
# The species classifier is an R object called `species_classifier`
if (!file.exists("data/data_drone/L3b/classifier-models/live-or-dead-classifier.rds") | !file.exists("data/data_drone/L3b/classifier-models/species-classifier.rds")) {
  source(here::here("workflow/22_build-crown-classifier.R"))
}

live_or_dead_classifier <- readr::read_rds("data/data_drone/L3b/classifier-models/live-or-dead-classifier.rds")
species_classifier <- readr::read_rds("data/data_drone/L3b/classifier-models/species-classifier.rds")

# Source in the allometric scaling models based on the ground data
# object is called `allometry_models`
source(here::here("workflow/23_allometric-scaling-models.R"))
allometry_models

crowns_with_reflectance <- readr::read_csv(file = here::here("data/data_drone/L3b/crowns-with-reflectance_35m-buffer.csv"))

classified_trees_nonspatial  <-
  crowns_with_reflectance %>% 
  dplyr::mutate(live = live_or_dead_classifier$levels[predict(live_or_dead_classifier, newdata = .)]) %>% 
  dplyr::mutate(live = as.numeric(as.character(live))) %>%
  dplyr::mutate(species = ifelse(live == 1, yes = species_classifier$levels[predict(species_classifier, newdata = .)], no = NA)) 

allometry_models %>% filter(species == "pipo") %>% pull(model) %>% magrittr::extract2(1) %>% predict(newdata = data.frame(height = 20))

classified_and_allometried_trees <-
  classified_trees_nonspatial %>% 
  dplyr::left_join(allometry_models, by = "species") %>% 
  dplyr::mutate(model = ifelse(live == 0, 
                               yes = allometry_models %>% 
                                 dplyr::filter(species == "pipo") %>% 
                                 pull(model), 
                               no = model)) %>% 
  dplyr::do(modelr::add_predictions(., model = first(.$model), var = "estimated_dbh")) %>% 
  dplyr::select(-model) %>% 
  dplyr::mutate(estimated_ba = (estimated_dbh / 2)^2 * pi / 10000)

classified_trees_4326 <-
  classified_and_allometried_trees %>% 
  split(f = .$crs) %>% 
  lapply(FUN = function(trees) {
    current_crs <- unique(trees$crs)
    
    trees4326 <-
      trees %>% 
      st_as_sf(coords = c("x", "y"), crs = current_crs, remove = FALSE) %>% 
      st_transform(4326) %>% 
      dplyr::rename(local_crs = crs,
                    local_x = x,
                    local_y = y) %>% 
      tidyr::separate(col = treeID, into = c("forest", "elev", "rep", "id"), sep = "_", remove = FALSE) %>% 
      dplyr::mutate(site = paste(forest, elev, rep, sep = "_")) %>% 
      dplyr::select(-id) %>% 
      dplyr::select(treeID, site, forest, elev, rep, live, species, height, estimated_dbh, estimated_ba, everything())
  }) %>% 
  do.call("rbind", .)

sf::st_write(classified_trees_4326, dsn = here::here("data/data_drone/L3b/classified-trees_all.gpkg"), delete_dsn = TRUE)
