library(sf)
library(tidyverse)
library(stringr)
library(purrr)
library(viridis)
library(raster)

flight_logs <- 
  list.files("data/data_output/eldo_3k_1_flight-path/", pattern = "[0-9].csv", full.names = TRUE) %>% 
  purrr::map(read_csv) %>% 
  bind_rows() %>% 
  st_as_sf(coords = c("Longitude", "Latitude")) %>%
  st_set_crs(4326)

photo_points <-
  flight_logs %>% 
  filter(isFlying == 1) %>%
  filter(flyState == 14) %>% 
  filter(Images > lag(Images))

dem <- raster::raster("data/srtm_30m.tif")

takeoff_elev <-
  flight_logs %>% 
  filter(`Altitude (m)` == 0) %>% 
  st_geometry() %>% 
  as("Spatial") %>% 
  raster::extract(dem, ., method = "bilinear") %>% 
  mean()

photo_points <-
  photo_points %>% 
  dplyr::mutate(photo_elev = raster::extract(dem, as(geometry, "Spatial"), method = "bilinear")) %>% 
  dplyr::mutate(agl = `Altitude (m)` + takeoff_elev - photo_elev) %>% 
  filter(agl > 100)

flight_bounds <- st_as_sfc(st_bbox(photo_points))
site_dem <- crop(x = dem, y = as(flight_bounds, "Spatial"), snap = "out")

plot(site_dem)
plot(photo_points$geometry, pch = 19, add = TRUE)
plot(flight_bounds, add = TRUE)