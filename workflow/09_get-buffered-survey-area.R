# Get the surveyed area for each site so we can convert site-level metrics of density and basal area
# to equivalent units

library(sf)
library(tidyverse)
library(here)

source(here::here("workflow/01_make-processing-checklist.R"))
sites_checklist

# With a given buffer size, how big will the total area surveyed be?
survey_area <-
  sites_checklist %>% 
  dplyr::mutate(geometry = map(.x = site, .f = function(current_site) {
    
    # One site has additional restrictions on its flight bounds in order to avoid any chance of visible private property in the processed imagery.
    if(current_site == "stan_3k_2") {
      site_bounds <- sf::st_read(here::here("data/data_drone/L0/mission-footprint/site-bounds/stan_3k_2_site-bounds_privacy.geoJSON"))
    } else {
      site_bounds <- sf::st_read(here::here(paste0("data/data_drone/L0/mission-footprint/site-bounds/", current_site, "_site-bounds.geoJSON")))
    }
    
    geometry <-
      site_bounds %>% 
      st_transform(3310) %>% 
      st_geometry()

    return(geometry)

  })) %>% 
  dplyr::mutate(geometry = st_as_sfc(purrr::flatten(geometry))) %>% 
  dplyr::select(site, geometry) %>% 
  sf::st_as_sf(crs = 3310) %>% 
  dplyr::mutate(survey_area = st_area(geometry)) %>% 
  dplyr::mutate(buffered_survey_area = st_area(st_buffer(geometry, -35)))

units(survey_area$survey_area) <- "ha"
units(survey_area$buffered_survey_area) <- "ha"

survey_area

sf::st_write(obj = survey_area, dsn = "data/data_drone/L0/surveyed-area-3310.gpkg", delete_dsn = TRUE)
