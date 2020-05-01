# Purpose: make shapefile of all identifiable ground plots from aerial orthomosaic

library(tidyverse)
library(sf)

source(here::here("workflow", "01_make-processing-checklist.R"))
sites_checklist

sites <- sites_checklist$site

plot_locations <-
  sites %>%
  map(.f = function(current_site) {
    
      current_site_plot_locations <-
        sf::st_read(here::here("data", "data_drone", "L1", "plot-locations", paste0(current_site, "_plot-locations.gpkg"))) %>%
        mutate(plot = paste(current_site, id, sep = "_")) %>%
        dplyr::arrange(id) %>%
        sf::st_zm() %>% 
        dplyr::mutate(local_x = st_coordinates(.)[, "X"],
                      local_y = st_coordinates(.)[, "Y"]) %>% 
        dplyr::mutate(local_crs = st_crs(.)$proj4string) %>%
        sf::st_transform(3310)
  }) %>% 
  do.call("rbind", .) %>% 
  dplyr::mutate(site = substr(plot, start = 1, stop = 9)) %>% 
  dplyr::rename(plot_id = id) %>% 
  dplyr::select(site, plot_id, plot, local_x, local_y, local_crs)

plot_locations

sf::st_write(obj = plot_locations, dsn = here::here("data", "data_drone", "L1", "plot-centers-identifiable-from-air_3310.gpkg"), delete_dsn = TRUE)
