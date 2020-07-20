# Size distribution of our classified trees

library(tidyverse)
library(sf)
library(here)

#### Create circle polygons representing the plot boundaries
# Can only really use the plot centers identfiable from the air to do this
# since the geolocation accuracy of the plots that can't be tied directly
# to the location of the orange cloth "X" in the orthophotos is not good
# enough
plot_radius <- sqrt((66*66) / pi) * (12 * 2.54 / 100) 
plot_locations <- 
  sf::st_read(here::here("data", "data_drone", "L1", "plot-centers-identifiable-from-air_3310.gpkg")) %>% 
  sf::st_transform(3310) %>% 
  sf::st_buffer(plot_radius) %>% 
  dplyr::select(-site, -local_x, -local_y, -local_crs, -plot_id)

#### Get trees identified (segmented and classified) from the air
d <- 
  sf::st_read(here::here("data", "data_drone", "L3b", "model-classified-trees_all.gpkg"), 
              stringsAsFactors = FALSE)  

# subset classified trees to just those within plots that are identifiable
# from the air
air_trees <-
  d %>% 
  sf::st_intersection(plot_locations) %>% 
  dplyr::mutate(species = ifelse(is.na(species), yes = "pipo", no = species),
                plot = as.character(plot)) %>% 
  dplyr::select(forest, elev, rep, site, plot, everything())

# Write the tree locations and attributes to plots
sf::st_write(obj = air_trees, dsn = here::here("data", "data_drone", "L3b", "model-classified-trees-within-ground-plots.gpkg"), append = FALSE)

