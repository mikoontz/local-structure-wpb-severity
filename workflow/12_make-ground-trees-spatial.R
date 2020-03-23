library(sf)
library(tidyverse)

# This script will generate a shapefile of all the locations of the trees that
# were measured in the ~110 field plots that are visible from the aerial imagery.
# We can only use this subset of plots because we base the center point location
# on exactly where it can be found in the aerial imagery, and then the tree
# locations are displacements from that center. With no center visible in an 
# aerial image, the only location data we have for plot centers comes from a
# handheld GPS, which could be off by several meters (and thus the displacement
# of all the trees would also be off)

source(here::here("workflow/01_make-processing-checklist.R"))

if(!file.exists("data/data_output/formatted-ground-data.csv")) {
  source(here::here("workflow/11_format-ground-data.R"))
}

d <- readr::read_csv("data/data_output/formatted-ground-data.csv")

# The merged versus the unmerged sites will have different numbers of bands
# Because the merged sites will have the X3 imagery incorporated, there will
# be an extra 3 bands for the ortho and the index outputs (the R, G, and B
# from the X3 camera)
# The R, G, and B bands from the X3 images will always be the final three bands
# if they exist.
# For the RedEdge-derived products, the bands go in the order of wavelength,
# from shortest to longest (B, G, R, RE, NIR)
# There is one Pix4D derived index (NDVI), which will go after the NIR band
# for the index mosaic

sites_to_process <- 
  sites_checklist %>% 
  dplyr::select(site) %>% 
  dplyr::pull() 

ground_trees <- 
  sites_to_process %>% 
  map(.f = function(current_site) {
    
    current_site_plot_locations <- sf::st_read(paste0("data/data_drone/L1/plot-locations/", current_site, "_plot-locations.gpkg"), stringsAsFactors = FALSE)
    
    current_site_ground_trees <-
      map(current_site_plot_locations$plot, .f = function(current_plot) {
        
        current_plot_ground_trees <- 
          d %>% 
          filter(plot == current_plot) %>% 
          dplyr::mutate(delta_x = cospi((1 / 2) - (azm / 180)) * dist) %>% 
          dplyr::mutate(delta_y = sinpi((1 / 2) - (azm / 180)) * dist) %>%
          left_join(current_site_plot_locations, by = "plot") %>%
          st_as_sf() %>%
          dplyr::mutate(x = st_coordinates(.)[, "X"],
                        y = st_coordinates(.)[, "Y"]) %>% 
          dplyr::mutate(new_x = x + delta_x,
                        new_y = y + delta_y) %>% 
          sf::st_drop_geometry() %>% 
          st_as_sf(coords = c("new_x", "new_y")) %>% 
          st_set_crs(st_crs(current_site_plot_locations))
        
        nn <- 
          nngeo::st_nn(x = current_plot_ground_trees, 
                       y = current_plot_ground_trees, 
                       k = min(c(4, nrow(current_plot_ground_trees))), 
                       returnDist = TRUE, sparse = FALSE, progress = FALSE)$dist %>% 
          do.call("rbind", .) %>% 
          as_tibble() %>% 
          dplyr::select(-1)
        
        if (ncol(nn) == 1) {nn <- nn %>% rename(nn1 = V2) %>% mutate(nn2 = NA, nn3 = NA)} else
          if (ncol(nn) == 2) {nn <- nn %>% rename(nn1 = V2, nn2 = V3) %>% mutate(nn3 = NA)} else
            if (ncol(nn) == 3) {nn <- nn %>% rename(nn1 = V2, nn2 = V3, nn3 = V4) }
        
        current_plot_ground_trees <- 
          current_plot_ground_trees %>% 
          bind_cols(nn)  %>% 
          st_transform(3310)
      }) %>% 
      do.call("rbind", .)
    
  }) %>%
  do.call("rbind", .)


sf::st_write(obj = ground_trees, 
             dsn = here::here(paste0("data/data_drone/L1/ground-trees.gpkg")), 
             delete_dsn = TRUE)
