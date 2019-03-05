library(sf)
library(raster)
library(tabularaster)
library(tidyverse)
library(viridis)
library(purrr)
library(velox)
library(viridis)

source(here::here("data/data_carpentry/make-processing-checklist.R"))
source(here::here("data/data_carpentry/extract-reflectance-within-crowns-from-mosaics-function.R"))

sites_checklist

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
  dplyr::filter(overwrite | !reflectance_extraction_check) %>% 
  dplyr::pull(site)

tic()
crowns_with_reflectance <-
  lapply(seq_along(sites_to_process), FUN = function(i) {
    current_site <- sites_to_process[i]
    
    crowns_path <- here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_crowns/", current_site, "_crowns.shp"))
    ttops_path <- here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_ttops/", current_site, "_ttops.shp"))
    crowns <- sf::st_read(crowns_path)
    ttops <- sf::st_read(ttops_path)
    
    index_path <- here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_index.tif"))
    index <- velox::velox(index_path)
    index_copy <- index$copy()
    
    current_crowns <- extract_reflectance_from_crowns(index = index_copy,
                                                      crowns = crowns,
                                                      ttops = ttops)
    
    current_crowns <- 
      current_crowns %>%
      dplyr::select(treeID, height, ch_area, x, y, b_mean, g_mean, r_mean, re_mean, nir_mean, ndvi_mean, rgi_mean, cire_mean, cig_mean, ndre_mean) %>% 
      dplyr::mutate(treeID = paste(current_site, treeID, sep = "_"),
                    crs = st_crs(.)$proj4string)
    
    suppressWarnings(dir.create(here::here(paste0("data/data_output/classified/model-classified/crowns-with-reflectance/", current_site, "_crowns-with-reflectance"))))
    
    sf::st_write(current_crowns, dsn = here::here(paste0("data/data_output/classified/model-classified/crowns-with-reflectance/", current_site, "_crowns-with-reflectance/", current_site, "_crowns-with-reflectance.shp")), delete_dsn = TRUE)
    
    return(st_drop_geometry(current_crowns))
    
  }) %>% 
  do.call("rbind", .)
toc()

# Write the full set of trees to a file
write_csv(crowns_with_reflectance, path = here::here(paste0("data/data_output/classified/model-classified/crowns-with-reflectance_all.csv")))

# Do some spatial subsetting to remove trees too close to the border of the study area

# With a given buffer size, how big will the total area surveyed be?
survey_area <-
  list.files(here::here("data/data_output/classified/model-classified/crowns-with-reflectance")) %>% 
  map_dbl(.f = function(current_crowns_path) {
    
    current_site <- substr(current_crowns_path, start = 1, stop = 9)
    
    # One site has additional restrictions on its flight bounds in order to avoid any chance of visible private property in the processed imagery.
    if(current_site == "stan_3k_2") {
      site_bounds <- sf::st_read(here::here("data/data_output/site_data/stan_3k_2/stan_3k_2_mission-footprint/stan_3k_2_site-bounds_privacy.geoJSON"))
    } else {
      site_bounds <- sf::st_read(here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_mission-footprint/", current_site, "_site-bounds.geoJSON")))
    }
    
    buffered_bounds <-
      site_bounds %>% 
      st_transform(3310) %>% 
      st_buffer(-35)
    
    return(st_area(buffered_bounds))
    
  })

current_crowns_path <-  list.files(here::here("data/data_output/classified/model-classified/crowns-with-reflectance"))[27]

# Subset the crowns here
crowns <-
  list.files(here::here("data/data_output/classified/model-classified/crowns-with-reflectance")) %>% 
  map(.f = function(current_crowns_path) {
    
    current_site <- substr(current_crowns_path, start = 1, stop = 9)
    current_crowns <- sf::st_read(here::here(paste0("data/data_output/classified/model-classified/crowns-with-reflectance/", current_crowns_path, "/", current_crowns_path, ".shp")))
    
    # One site has additional restrictions on its flight bounds in order to avoid any chance of visible private property in the processed imagery.
    if(current_site == "stan_3k_2") {
      site_bounds <- sf::st_read(here::here("data/data_output/site_data/stan_3k_2/stan_3k_2_mission-footprint/stan_3k_2_site-bounds_privacy.geoJSON"))
    } else {
      site_bounds <- sf::st_read(here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_mission-footprint/", current_site, "_site-bounds.geoJSON")))
    }
    
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
write_csv(crowns, path = here::here(paste0("data/data_output/classified/model-classified/crowns-with-reflectance_35m-buffer.csv")))
