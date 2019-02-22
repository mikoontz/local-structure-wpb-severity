# Purpose: extract climatic water deficit values for each tree or each site

# Data are from the California Climate Commons

library(raster)
library(sf)
library(tidyverse)
library(here)

cwd <- raster::raster(here::here("data/features/cwd1981_2010_ave_HST_1550861123/cwd1981_2010_ave_HST_1550861123.tif"))

# The .prj file doesn't seem to be reading in properly with the .tif, but we can look at it in a text editor and see that it is EPSG3310

crs(cwd) <- st_crs(3310)$proj4string

plot(cwd)

# There are approximately 4 CWD pixels (at 270m spatial resolution) per
# 40ha site (if square, about 625m on a side).
# First pass, let's just extract the CWD values at the centers of each 
# study area.

# Remember to put them in the same coordinate reference system
# as the CWD raster
site_centers <- st_read(here::here("data/features/plot-centers_ground-gps-measured.kml")) %>% 
  st_transform(3310) %>% 
  tidyr::separate(col = Name, into = c("forest", "elevation_band", "rep", "nickname", "plot_id"), sep = "_") %>% 
  dplyr::select(-Description) %>% 
  dplyr::group_by(forest, elevation_band, rep) %>% 
  dplyr::summarize() %>% 
  sf::st_centroid() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(cwd = raster::extract(cwd, ., method = "bilinear")) %>% 
  dplyr::mutate(elevation_band = as.numeric(as.character(elevation_band)))

plot(cwd)
plot(site_centers$geometry, add = TRUE, pch = 19)

site_centers

ggplot(site_centers, aes(x = elevation_band, y = cwd, col = forest)) +
  geom_smooth(method = "lm") +
  scale_color_viridis_d()

ggplot(site_centers, aes(col = cwd, shape = forest, size = elevation_band)) +
  geom_sf() +
  scale_color_viridis_c()

