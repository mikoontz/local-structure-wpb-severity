# Purpose is to take the raw flight logs from a DJI drone flown using Map Pilot for iOS
# and create the "mission footprint" represented by 3 geospatial products:
# 1) the 30m elevation raster layer that was used by the Map Pilot flight software to
# adjust flight altitude to stay a constant altitude above ground level. This will also
# be useful for "normalizing" .las point clouds using the lidR package later
# 2) a convex hull geoJSON representing the extent of the area where pictures were taken
# This can be used to subset the .las point cloud
# 3) The photo points themselves

# Load libraries
library(sf)
library(tidyverse)
library(purrr)
library(raster)

source("workflow/01_make-processing-checklist.R")

# overwrite variable if user wants to rewrite all exported geospatial files (dem of each site, bounding box around
# each site, photo points for each site)
overwrite <- FALSE

# character vector of all study sites that don't have the mission footprint figured out yet
sites_to_process <-
  sites_checklist %>% 
  dplyr::filter(!mission_footprint_check) %>% 
  dplyr::select(site) %>% 
  pull()

# the 30m resolution SRTM digital elevation model for the Sierra Nevada region
# No direct source for this exists in R (not raster, elevatr, or FedData), so
# I downloaded the whole thing from Google Earth Engine and subset it using
# this code

dem <- raster::raster("data/data_raw/srtm_30m.tif")


# Create directories if necessary
if (!dir.exists("data/data_drone/L0/mission-footprint/srtm30m")) {
  dir.create("data/data_drone/L0/mission-footprint/srtm30m")
}

if (!dir.exists("data/data_drone/L0/mission-footprint/photo-points")) {
  dir.create("data/data_drone/L0/mission-footprint/photo-points")
}

if (!dir.exists("data/data_drone/L0/mission-footprint/site-bounds")) {
  dir.create("data/data_drone/L0/mission-footprint/site-bounds")
}

# Iterate through all the available sites
# For loop is much more inuitive to use here (IMO)
for (i in seq_along(sites_to_process)) {
  # get the character string representing the ith site
  current_site <- sites_to_process[i]
  
  # read all the individual raw flight logs from the specified directory
  # turn the csv files into shape files by assigning the appropriate columns
  # to be the coordinates
  
  flight_logs_list <- 
    list.files(paste0("data/data_drone/L0/flight-logs/", current_site), pattern = "[0-9].csv", full.names = TRUE) %>% 
    purrr::map(read_csv) %>% 
    purrr::map(.f = function(x) {
      x %>% 
        sf::st_as_sf(coords = c("Longitude", "Latitude")) %>%
        sf::st_set_crs(4326)
      })

  # assume the first row of each log file is the location of the takeoff point
  # All altitude calculations are relative to this point, so getting the
  # elevation of this point from the DEM tells us the offset
  takeoff_points <-
    purrr::map(flight_logs_list, .f = function(x) {
      x %>% 
        dplyr::slice(1)
    })

  # filter the spatial flight log to just the rows where images incremented
  # and the drone is in a "flying" condition
  photo_points_list <-
    purrr::map(flight_logs_list, .f = function(x) {
      x %>% 
        dplyr::filter(isFlying == 1) %>% 
        dplyr::filter(flyState == 14) %>% 
        dplyr::filter(Images > lag(Images))
    })
  
  # Iterate over the list elements representing the photo points for each flight
  # and the take off point for each flight to calculate the elevation for each
  # photo point (on the ground) and the agl (above ground level) measure for
  # each photo point (by taking the ground elevation, the takeoff elevation, and
  # the relative altitude offset from the takeoff elevation into account)
  # Note: the agl should be fairly consistent across the whole mission.
  photo_points <- 
    purrr::map2(.x = photo_points_list, .y = takeoff_points, .f = function(x, y, ...) {
      x %>% 
        dplyr::mutate(photo_elev = raster::extract(dem, as(geometry, "Spatial"), method = "bilinear")) %>% 
        dplyr::mutate(agl = `Altitude (m)` + raster::extract(dem, as(st_geometry(y), "Spatial")) - photo_elev) %>% 
        dplyr::filter(agl > 105)
    }) %>% 
    do.call(rbind, .)
  
  # create a vector object representing the convex hull of the photopoints
  site_bounds <- 
    photo_points %>% 
    st_union() %>% 
    st_convex_hull()
  
  site_dem <- raster::crop(x = dem, y = as(site_bounds, "Spatial"), snap = "out")
  
  raster::writeRaster(x = site_dem, 
                      filename = paste0("data/data_drone/L0/mission-footprint/srtm30m/", current_site, "_srtm30m.tif"), 
                      overwrite = overwrite)

  sf::st_write(obj = photo_points, 
               dsn = paste0("data/data_drone/L0/mission-footprint/photo-points/", current_site, "_photo-points.geoJSON"), 
               delete_dsn = overwrite)

  sf::st_write(obj = site_bounds, 
               dsn = paste0("data/data_drone/L0/mission-footprint/site-bounds/", current_site, "_site-bounds.geoJSON"), 
               delete_dsn = overwrite)
  
}
