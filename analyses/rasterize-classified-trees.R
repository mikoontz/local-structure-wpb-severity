# Purpose: rasterize the classified trees to 20m x 20m cells in order to mimic the footprint of the field plots

library(tidyverse)
library(raster)
library(sf)
library(velox)
library(here)
library(tictoc)

if(file.exists(here::here("analyses/analyses_output/classified-trees.geojson"))) {
  
  classified_trees <- 
    sf::st_read(here::here("analyses/analyses_output/classified-trees.geojson")) %>% 
    sf::st_transform(3310)
  
} else {
  stop("Trees haven't been classified yet! See the analyses/classify-crown-segments.R script.")
}

# Now there is an R object in the environment called "sites_checklist" that has
# infomation about how far along all processing steps are.
source("data/data_carpentry/make-processing-checklist.R")

source(here::here("data/data_carpentry/extract-cwd-from-locations.R"))
# R object is called `cwd_data`


unusable_sites <- c("eldo_4k_3", # too many blocks
                    "stan_4k_3", # too many blocks
                    "stan_5k_3", # too many blocks
                    "sequ_4k_2") # middle section flown on a separate day and the stitch looks terrible

# These sites were processed with their X3 and RedEdge imagery combined so some of their
# output products will be in a slightly different place in the project directory
merged_sites <- c("eldo_3k_2",
                  "eldo_3k_3",
                  "eldo_4k_2")

# This is where I can put in sites that need their processing redone. An empty 
# string means that no already-processed site output will be overwritten
# (but sites that have yet to be processed will still have their processing done)
sites_to_overwrite <- "all"
sites_checklist$overwrite <- ifelse(sites_to_overwrite == "all", yes = TRUE, no = FALSE)

sites_checklist[sites_checklist$site %in% sites_to_overwrite, "overwrite"] <- TRUE

sites_to_process <- 
  sites_checklist %>% 
  dplyr::filter(!(site %in% unusable_sites)) %>%
  dplyr::filter(overwrite | !rasterized_trees_check) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()

results_list <- vector(mode = "list", length = length(sites_to_process))

for(i in seq_along(sites_to_process)) {
  
  current_site <- sites_to_process[i]
  
  current_trees <- 
    classified_trees %>%
    dplyr::filter(site == current_site)
  
  # Build the template raster using the site bounds as an outer border
  # One site has additional restrictions on its flight bounds in order to avoid any chance of visible private property in the processed imagery.
  if(current_site == "stan_3k_2") {
    site_bounds <- sf::st_read(here::here("data/data_output/site_data/stan_3k_2/stan_3k_2_mission-footprint/stan_3k_2_site-bounds_privacy.geoJSON"))
  } else {
    site_bounds <- sf::st_read(here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_mission-footprint/", current_site, "_site-bounds.geoJSON")))
  }
  
  buffered_bounds <-
    site_bounds %>% 
    st_transform(3310) %>% 
    st_buffer(-40) # buffer in a little further than the trees were to reduce edge effects
  
  raster_template <- raster::raster(buffered_bounds, res = 20)
  
  live_and_dead <- raster::rasterize(x = current_trees, 
                                     y = raster_template, 
                                     field = "live", 
                                     background = 0, 
                                     fun = function(x, ...) {
                                       c(length(which(x == 1)), 
                                         length(which(x == 0)))
                                     })
  
  pipo_count <- raster::rasterize(x = current_trees %>% 
                                    dplyr::filter((species == "pipo") | live == 0), 
                                  y = raster_template, 
                                  field = "species", 
                                  background = 0, 
                                  fun = "count")
  
  non_pipo_count <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 1 & species != "pipo")), 
                                      y = raster_template, 
                                      field = "species", 
                                      background = 0, 
                                      fun = "count")
  
  pipo_basal_area <- raster::rasterize(x = current_trees %>% 
                                         dplyr::filter((species == "pipo") | live == 0), 
                                       y = raster_template, 
                                       field = "estimated_ba", 
                                       background = 0, 
                                       fun = sum)
  
  non_pipo_basal_area <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 1 & species != "pipo")), 
                                      y = raster_template, 
                                      field = "estimated_ba", 
                                      background = 0, 
                                      fun = sum)
  
  
  results_raster <- raster::stack(live_and_dead, pipo_count, non_pipo_count, pipo_basal_area, non_pipo_basal_area)
  names(results_raster) <- c("live_count", "dead_count", "pipo_count", "non_pipo_count", "pipo_ba", "non_pipo_ba")
  
  writeRaster(x = results_raster, filename = here::here(paste0("analyses/analyses_output/rasterized-trees/", current_site, "_rasterized-trees.tif")), overwrite = TRUE)
  
  results_df <- 
    results_raster %>% 
    as.data.frame(xy = TRUE) %>% 
    dplyr::mutate(site = current_site) %>% 
    tidyr::separate(col = site, into = c("forest", "elev", "rep"), remove = FALSE) %>% 
    dplyr::mutate(crs = 3310)
  
  results_list[[i]] <- results_df
}

final_results <- do.call("rbind", results_list)

readr::write_csv(final_results, here::here("analyses/analyses_output/data-from-rasterized-classified-trees.csv"))

