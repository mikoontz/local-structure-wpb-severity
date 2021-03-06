# Purpose: extract climatic water deficit values for each tree or each site

# Data are from the California Climate Commons

library(raster)
library(sf)
library(tidyverse)
library(here)
library(lubridate)

cwd <- raster::raster(here::here("data/data_raw/cwd1981_2010_ave_HST_1550861123/cwd1981_2010_ave_HST_1550861123.tif"))

# The .prj file doesn't seem to be reading in properly with the .tif, but we can look at it in a text editor and see that it is EPSG3310

# Note that using a PROJ4 string probably won't work going forward, given the
# shift to PROJ6
crs(cwd) <- sf::st_crs(3310)$proj4string

# There are approximately 4 CWD pixels (at 270m spatial resolution) per
# 40ha site (if square, about 625m on a side).
# First pass, let's just extract the CWD values at the centers of each 
# study area.

# Remember to put them in the same coordinate reference system
# as the CWD raster
site_centers <- 
  st_read(here::here("data/data_raw/plot-centers_ground-gps-measured.kml")) %>% 
  st_transform(3310) %>% 
  tidyr::separate(col = Name, into = c("forest", "elevation_band", "rep", "nickname", "plot_id"), sep = "_") %>% 
  dplyr::select(-Description) %>% 
  dplyr::group_by(forest, elevation_band, rep) %>% 
  dplyr::summarize() %>% 
  sf::st_centroid() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(site_cwd = raster::extract(cwd, ., method = "bilinear")) %>% 
  dplyr::mutate(elevation_band = as.numeric(as.character(elevation_band))) %>% 
  dplyr::mutate(forest = substr(forest, start = 1, stop = 4)) %>% 
  dplyr::mutate(elevation_band = paste0(substr(elevation_band, start = 1, stop = 1), "k")) %>% 
  dplyr::mutate(site = paste(forest, elevation_band, rep, sep = "_"))

# plot(cwd)
# plot(site_centers$geometry, add = TRUE, pch = 19)

site_centers

# Quick exploratory plots
# ggplot(site_centers, aes(x = elevation_band, y = cwd, col = forest)) +
#   geom_smooth(method = "lm") +
#   scale_color_viridis_d()
# 
# ggplot(site_centers, aes(col = cwd, shape = forest, size = elevation_band)) +
#   geom_sf() +
#   scale_color_viridis_c()

#### CWD data per plot
plot_centers <- 
  st_read(here::here("data/data_raw/plot-centers_ground-gps-measured.kml")) %>% 
  st_transform(3310) %>% 
  tidyr::separate(col = Name, into = c("forest", "elevation_band", "rep", "nickname", "plot_id"), sep = "_") %>% 
  dplyr::select(-Description) %>% 
  dplyr::mutate(plot_cwd = raster::extract(cwd, ., method = "bilinear")) %>% 
  dplyr::mutate(elevation_band = as.numeric(as.character(elevation_band))) %>% 
  dplyr::mutate(forest = substr(forest, start = 1, stop = 4)) %>% 
  dplyr::mutate(elevation_band = paste0(substr(elevation_band, start = 1, stop = 1), "k")) %>% 
  dplyr::mutate(site = paste(forest, elevation_band, rep, sep = "_"))

# This is the raw CWD data. 
cwd_site_data <- 
  site_centers %>% 
  dplyr::select(site, site_cwd) %>% 
  sf::st_drop_geometry()

cwd_plot_data <- 
  plot_centers %>% 
  dplyr::select(site, plot_id, rep, plot_cwd) %>% 
  dplyr::mutate(plot_id = as.numeric(plot_id) - (5 * (as.numeric(rep) - 1))) %>% 
  dplyr::mutate(plot = paste(site, plot_id, sep = "_")) %>% 
  dplyr::select(site, plot, plot_cwd) %>% 
  sf::st_drop_geometry()

# What do these CWD values mean for the overall PIPO distribution
# in the Sierra Nevada?

sn <- sf::st_read(here::here("data/data_output/sierra-nevada-jepson.gpkg"))

sn_pipo <-
  data.table::fread(here::here("data/data_raw/California_Species_clean_All_epsg_3310.csv")) %>% 
  dplyr::as_tibble() %>% 
  dplyr::mutate(current_genus = tolower(current_genus),
                current_species = tolower(current_species)) %>% 
  dplyr::filter(current_genus == "pinus") %>% 
  dplyr::filter(current_species == "ponderosa") %>% 
  sf::st_as_sf(coords = c("x_epsg_3310", "y_epsg_3310"), crs = 3310) %>% 
  dplyr::select(id, early_julian_day, late_julian_day, verbatim_date, elevation) %>% 
  sf::st_intersection(sn) %>% 
  dplyr::mutate(date = parse_date_time(early_julian_day, c("mdy", "ymd", "ymdHM"))) %>% 
  dplyr::mutate(year = year(date)) %>% 
  dplyr::mutate(cwd = raster::extract(cwd, ., method = "bilinear"))

# sn_pipo has 179 rows and thus represents 179 herbaria records
# of ponderosa pine for which the CWD was extracted (using CWD as defined # by the Basin Characterization Model)
sf::st_write(obj = sn_pipo, dsn = "data/data_output/sierra-nevada-pipo-cwd.gpkg")

mean_cwd_sn_pipo <- mean(sn_pipo$cwd, na.rm = TRUE)
sd_cwd_sn_pipo <- sd(sn_pipo$cwd, na.rm = TRUE)

cwd_site_data <-
  cwd_site_data %>% 
  dplyr::mutate(site_cwd_zscore = (site_cwd - mean_cwd_sn_pipo) / sd_cwd_sn_pipo)

readr::write_csv(cwd_site_data, here::here("data/data_output/cwd-data.csv"))

cwd_plot_data <-
  cwd_plot_data %>% 
  dplyr::mutate(plot_cwd_zscore = (plot_cwd - mean_cwd_sn_pipo) / sd_cwd_sn_pipo)

readr::write_csv(cwd_plot_data, here::here("data/data_output/cwd-plot-data.csv"))
