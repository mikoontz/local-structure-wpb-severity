# Purpose: make shapefile of all identifiable ground plots from aerial orthomosaic

library(tidyverse)
library(sf)

sites <- list.files("data/data_output/site_data")
plot_locations <-
  sites %>%
  map(.f = function(current_site) {
    
    if (current_site %in% merged_sites) {
      current_site_plot_locations <-
        sf::st_read(paste0("data/data_output/site_data/", current_site, "/", current_site, "_plot-locations/", current_site, "_plot-locations.shp")) %>%
        mutate(plot = paste(current_site, id, sep = "_")) %>%
        dplyr::arrange(id) %>%
        sf::st_zm() %>% 
        dplyr::mutate(local_x = st_coordinates(.)[, "X"],
                      local_y = st_coordinates(.)[, "Y"]) %>% 
        dplyr::mutate(local_crs = st_crs(.)$proj4string) %>%
        sf::st_transform(3310)
    } else {
      current_site_plot_locations <-
        sf::st_read(paste0("data/data_output/site_data/", current_site, "/", current_site, "_plot-locations_re/", current_site, "_plot-locations_re.shp")) %>%
        mutate(plot = paste(current_site, id, sep = "_")) %>%
        dplyr::arrange(id) %>%
        sf::st_zm() %>% 
        dplyr::mutate(local_x = st_coordinates(.)[, "X"],
                      local_y = st_coordinates(.)[, "Y"]) %>% 
        dplyr::mutate(local_crs = st_crs(.)$proj4string) %>%
        sf::st_transform(3310)
    }
  }) %>% 
  do.call("rbind", .) %>% 
  dplyr::mutate(site = substr(plot, start = 1, stop = 9)) %>% 
  dplyr::rename(plot_id = id) %>% 
  dplyr::select(site, plot_id, plot, local_x, local_y, local_crs)

plot_locations

sf::st_write(obj = plot_locations, dsn = "data/data_output/plot-centers-identifiable-from-air_3310.gpkg")
