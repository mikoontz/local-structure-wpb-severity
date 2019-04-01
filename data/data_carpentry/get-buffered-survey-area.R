# Get the surveyed area for each site so we can convert site-level metrics of density and basal area
# to equivalent units

library(sf)
library(tidyverse)
library(here)

source(here::here("data/data_carpentry/make-processing-checklist.R"))
sites_checklist

sites_to_process <-
  sites_checklist %>% 
  dplyr::filter(crowns_check) %>% 
  dplyr::pull(site)

# With a given buffer size, how big will the total area surveyed be?
survey_area <-
  sites_checklist %>% 
  dplyr::mutate(geometry = map(.x = site, .f = function(current_site) {
    
    # One site has additional restrictions on its flight bounds in order to avoid any chance of visible private property in the processed imagery.
    if(current_site == "stan_3k_2") {
      site_bounds <- sf::st_read(here::here("data/data_output/site_data/stan_3k_2/stan_3k_2_mission-footprint/stan_3k_2_site-bounds_privacy.geoJSON"))
    } else {
      site_bounds <- sf::st_read(here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_mission-footprint/", current_site, "_site-bounds.geoJSON")))
    }
    
    geometry <-
      site_bounds %>% 
      st_transform(3310) %>% 
      st_geometry()

    # buffered_bounds <-
    #   site_bounds %>% 
    #   st_transform(3310) %>% 
    #   st_buffer(-35)
    # 
    return(geometry)
    # return(st_area(buffered_bounds))
    
  })) %>% 
  dplyr::mutate(geometry = st_as_sfc(purrr::flatten(geometry))) %>% 
  dplyr::select(site, geometry) %>% 
  sf::st_as_sf(crs = 3310) %>% 
  dplyr::mutate(survey_area = st_area(geometry)) %>% 
  dplyr::mutate(buffered_survey_area = st_area(st_buffer(geometry, -35)))

units(survey_area$survey_area) <- "ha"
units(survey_area$buffered_survey_area) <- "ha"

survey_area

sf::st_write(obj = survey_area, dsn = "data/data_output/surveyed-area-3310.gpkg")