# What trends materialize in the drone data if we only use aerially detected and
# classified trees that fall within the plot footprints? Note that this requires
# that we only use the plots that are visible from the orthomosaic (the orange
# X was visible) so that we can get a precise boundary of the ground plots.

library(tidyverse)
library(sf)
library(here)

cwd <-
  read_csv(here::here("data", "data_output", "cwd-data.csv"))

# Get classified trees
# d <- 
#   sf::st_read(here::here("data/data_drone/L3b/model-classified-trees_all_height-corrected.gpkg"), 
#               stringsAsFactors = FALSE)  
# 
# # Join classified trees with CWD data
# dd <-
#   d %>%
#   dplyr::left_join(cwd, by = "site")

# get plot locations
# plot_radius <- sqrt((66*66) / pi) * (12 * 2.54 / 100) 
# plot_locations <- 
#   sf::st_read(here::here("data", "data_drone", "L1", "plot-centers-identifiable-from-air_3310.gpkg")) %>% 
#   sf::st_buffer(plot_radius) %>% 
#   dplyr::select(-site, -local_x, -local_y, -local_crs)

# subset classified trees to just those within plots that are identifiable
# from the air
air_trees <- sf::st_read(here::here("data", "data_drone", "L3b", "model-classified-trees-within-ground-plots_height-corrected.gpkg"), stringsAsFactors = FALSE)

air_trees_by_plot <-
  air_trees %>% 
  dplyr::mutate(live = as.numeric(as.character(live))) %>% 
  st_drop_geometry() %>% 
  group_by(plot) %>% 
  summarize(mean_height = mean(height[species == "pipo"]),
            prop_host = length(which(species == "pipo")) / n(),
            prop_dead = 1 - mean(live[species == "pipo"]),
            overall_density = n(),
            cwd = mean(site_cwd_zscore)) %>% 
  dplyr::mutate(plot = as.character(plot)) %>% 
  dplyr::mutate(air_ground = "drone data over field plots")

# Get random pixels stratified by each site that describe mean host height
# and proportion of host mortality
r_data <- readr::read_csv("data/data_drone/L4/data-from-rasterized-classified-trees.csv")

set.seed(123)
set.seed(1609)
r_data_sample <-
  r_data %>% 
  dplyr::filter(pipo_and_dead_count > 0) %>% 
  group_by(site) %>% 
  sample_n(5) %>%
  dplyr::mutate(prop_dead = dead_count / pipo_and_dead_count) %>% 
  dplyr::ungroup()

air_trees_by_plot_from_raster <-
  r_data_sample %>% 
  dplyr::select(plot = site, mean_height = pipo_and_dead_mean_height, prop_dead) %>% 
  dplyr::mutate(plot = 1:n())

# Get the ground trees from the 110 plots that are identifiable by air
ground_trees <-
  read_csv(here::here("data", "data_output", "formatted-ground-data.csv")) %>%
  filter(is.na(year_fall)) %>%
  left_join(cwd, by = "site") %>% 
  filter(plot %in% unique(plot_locations$plot))

ground_trees_by_plot <-
  ground_trees %>% 
  group_by(plot) %>% 
  summarize(mean_height = mean(height[species == "PIPO"]),
            prop_host = length(which(species == "PIPO")) / n(),
            prop_dead = 1 - mean(live[species == "PIPO"]),
            overall_density = n(),
            cwd = mean(site_cwd_zscore)) %>% 
  dplyr::mutate(plot = as.character(plot)) %>% 
  dplyr::mutate(air_ground = "field plot data")

# Plot air and ground assessment

air_ground_raster <-
  dplyr::select(air_trees_by_plot, plot, mean_height, prop_dead, air_ground) %>% 
  rbind(dplyr::select(ground_trees_by_plot, plot, mean_height, prop_dead, air_ground)) %>% 
  rbind(dplyr::mutate(air_trees_by_plot_from_raster, air_ground = "drone data randomly selected\nat each site"))

ggplot(air_ground_raster, aes(x = mean_height, y = prop_dead, color = air_ground)) +
  geom_point() +
  geom_smooth() +
  scale_color_viridis_d() +
  theme_bw() +
  labs(x = "Mean height of host trees in plot (m)",
       y = "Proportion of host trees dead",
       color = "Data source")

ggsave(filename = "figures/supplemental-figure_height-mortality-relationship-depends-on-data-source.png")
