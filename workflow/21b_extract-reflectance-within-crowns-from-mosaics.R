library(sf)
library(raster)
library(tabularaster)
library(tidyverse)
library(viridis)
library(purrr)
library(velox)
library(viridis)

source(here::here("workflow/01_make-processing-checklist.R"))
source(here::here("workflow/21a_extract-reflectance-within-crowns-from-mosaics-function.R"))

sites_checklist

if(!dir.exists("data/data_drone/L3b/crowns-with-reflectance/")) {
  dir.create("data/data_drone/L3b/crowns-with-reflectance/", recursive = TRUE)
}

# This is where I can put in sites that need their processing redone. An empty 
# string means that no already-processed site output will be overwritten
# (but sites that have yet to be processed will still have their processing done)
sites_to_overwrite <- "all"
sites_checklist$overwrite <- ifelse(sites_to_overwrite == "all", yes = TRUE, no = FALSE)

sites_checklist[sites_checklist$site %in% sites_to_overwrite, "overwrite"] <- TRUE

sites_to_process <- 
  sites_checklist %>% 
  dplyr::filter(overwrite | !reflectance_extraction_check) %>% 
  dplyr::pull(site)

crowns_with_reflectance <-
  lapply(seq_along(sites_to_process), FUN = function(i) {
    current_site <- sites_to_process[i]
    
    ttops <- sf::st_read(here::here(paste0("data/data_drone/L3a/ttops/", current_site, "_ttops.gpkg")))
    crowns <- sf::st_read(here::here(paste0("data/data_drone/L3a/crowns/", current_site, "_crowns.gpkg")))
    
    index <- velox::velox(here::here(paste0("data/data_drone/L2/index/", current_site, "_index.tif")))
    index_copy <- index$copy()
    
    current_crowns <- extract_reflectance_from_crowns(index = index_copy,
                                                      crowns = crowns,
                                                      ttops = ttops)
    
    current_crowns <- 
      current_crowns %>%
      dplyr::select(treeID, height, ch_area, x, y, b_mean, g_mean, r_mean, re_mean, nir_mean, ndvi_mean, rgi_mean, cire_mean, cig_mean, ndre_mean) %>% 
      dplyr::mutate(treeID = paste(current_site, treeID, sep = "_"),
                    crs = st_crs(.)$proj4string)
    
    sf::st_write(obj = current_crowns, 
                 dsn = here::here(paste0("data/data_drone/L3b/crowns-with-reflectance/", current_site, "_crowns-with-reflectance.gpkg")), delete_dsn = TRUE)
    
    return(st_drop_geometry(current_crowns))
    
  }) %>% 
  do.call("rbind", .)

# Write the full (nonspatial) set of trees to a file
write_csv(crowns_with_reflectance, path = here::here(paste0("data/data_drone/L3b/crowns-with-reflectance_all.csv")))

# Do some spatial subsetting to remove trees too close to the border of the study area
surveyed_area <- sf::st_read("data/data_drone/L0/surveyed-area-3310.gpkg")

# Subset the crowns here
crowns <-
  sites_checklist$site %>% 
  map(.f = function(current_site) {
    
    current_crowns <- sf::st_read(here::here(paste0("data/data_drone/L3b/crowns-with-reflectance/", current_site, "_crowns-with-reflectance.gpkg")))
    
    site_bounds <- dplyr::filter(surveyed_area, site == current_site)

    buffered_bounds <-
      site_bounds %>% 
      st_transform(st_crs(current_crowns)) %>% 
      st_buffer(-35)

    buffered_trees <-
      current_crowns %>%
      st_drop_geometry() %>% 
      st_as_sf(coords = c("x", "y"), crs = st_crs(buffered_bounds)) %>% 
      st_intersection(buffered_bounds) %>% 
      dplyr::mutate(x = st_coordinates(.)[, 1],
                    y = st_coordinates(.)[, 2]) %>% 
      st_drop_geometry()

    return(buffered_trees)
    
  }) %>% 
  do.call("rbind", .)
  
# write the buffered set of crowns to a file
write_csv(crowns, path = here::here(paste0("data/data_drone/L3b/crowns-with-reflectance_35m-buffer.csv")))
