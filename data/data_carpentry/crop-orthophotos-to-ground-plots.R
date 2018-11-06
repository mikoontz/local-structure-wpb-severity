# This script will create small, cropped rasters of the Pix4D orthophoto outputs that overlap the plots (where the orange X's over the plot centers were visible)
# By saving the much smaller rasters as separte files, they will be much faster to work with and test different crown segmentation methods on the in situ data.
# Before chopping up the original orthophoto, each raster (of the 40ha area; ~625m x 625m) is ~600MB, but we only use between 1 and 5, ~25m x 25m sections of the orthophoto that overlap the ground plots. Performing multiple tree top detection and crown segmentation algorithms on the 144 plots with known coordinates in the orthophoto (because their orange X's were visible in the orthophoto) will be much faster using the smaller rasters.

library(sf)
library(tidyverse)
library(purrr)
library(raster)
library(lidR)
# library(viridis)
# library(ForestTools)
# library(gstat)

# Now there is an R object in the environment called "sites_checklist" that has
# infomation about how far along all processing steps are.

source("data/data_carpentry/make_processing-checklist.R")
source("data/data_carpentry/format_ground-data.R")

sites_to_process <- 
  sites_checklist %>% 
  dplyr::filter(!plot_remote_data_check) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()

for (i in seq_along(sites_to_process)) {
   # get the character string representing the ith site
  current_site <- sites_to_process[i]
  
  current_dir <- paste0("data/data_output/site_data/", current_site, "/")
  
  # Create a new directory to put the plot-specific remotely sensed data
  dir.create(paste0(current_dir, current_site, "_plot-remote-data"))
  
  # The Digital Terrain Model (dtm) is the 2m resolution "ground" underneath
  # the current site. Created using CloudCompare and the Cloth Simulator Filter
  dtm <- raster::raster(paste0(current_dir, current_site, "_2m-dtm.tif"))
  
  # The Digital Surface Model (dsm) is the ~5cm resolution raster representing
  # the surface (ground + objects on top) that the drone flew over
  dsm <- raster::raster(paste0(current_dir, "3_dsm_ortho/1_dsm/", current_site, "_x3_dsm.tif"))
  
  # Get the orthophoto also
  ortho <- raster::brick(paste0("data/data_output/site_data/", current_site, "/3_dsm_ortho/2_mosaic/", current_site, "_x3_transparent_mosaic_group1.tif"))
  
  # Make sure the Coordinate Reference Systems (crs) are the same
  raster::crs(dtm) <- raster::crs(ortho)
  
  # The flight bounds is the polygon surrounding the images that the drone took. Determined
  # in the "convert_flight-path-to-site-boundary.R script using the elevation model used
  # during flight planning, the takeoff altitude, the ground altitude underneath each photo,
  # and the altitude difference between each photo and the takeoff point. Essentially, we 
  # crop all the products (e.g., the dsm, the dtm, the lidar point cloud) to just be 
  # area within the flight path, rather than some of the spillover information beyond the
  # area that the drone directly flew over.
  # Could consider buffering this further (inward-- so using a negative buffer value) to 
  # reduce edge effects
  site_bounds <- 
    sf::st_read(paste0(current_dir, current_site, "_mission-footprint/", current_site, "_site-bounds.geoJSON")) %>% 
    sf::st_transform(sp::proj4string(dtm))
  
  # Cropping the dtm and dsm to the site bounds
  site_dtm <- raster::crop(dtm, site_bounds)
  site_dsm <- raster::crop(dsm, site_bounds)
  site_ortho <- raster::crop(ortho, site_bounds)
  
  # Using bilinear interpolation to downsample the 2m resolution DTM to have the
  # same resolution as the dsm (~5cm, but slightly different for each site)
  dtm_resamp <- raster::resample(site_dtm, site_dsm, method = "bilinear")
  
  # The Canopy Height Model (chm) is the dsm (vegetation + ground) minus the dtm (ground)
  # to give just the height of the vegetation.
  chm <- site_dsm - dtm_resamp
  
  # Smooth out the chm and set any negative values to 0 (meaning "ground") following
  # advice from Zagalikis, Cameron, and Miller (2004) and references therein
  # More recently, a 3x3 pixel smoothing filter was specifically suggested as ideal
  # for sUAS derived chm by Mohan et al. (2017)
  # w <- focalWeight(x = chm, d = 0.5, type = "circle")
  # w
  chm_smooth <- raster::focal(chm, w = matrix(1, nrow = 3, ncol = 3))
  chm_smooth[raster::getValues(chm_smooth) < 0] <- 0
  
  plot_radius <- sqrt((66*66) / pi) * 12 *2.54 / 100
  
  # Get plot locations as determined by visually inspecting orthophotos and finding orange X's laid
  # across the plot centers
  current_site_plot_locations <- 
    sf::st_read(paste0(current_dir, current_site, "_plot-locations/", current_site, "_plot-locations.shp")) %>% 
    mutate(plot = paste(current_site, id, sep = "_")) %>%
    dplyr::arrange(id) %>% 
    st_zm()
  
  plots_visible_from_ortho <- unique(current_site_plot_locations$id)
  
  
  for (j in seq_along(plots_visible_from_ortho)) {
    current_plot_id <- plots_visible_from_ortho[j]  
    
    current_plot_location <-
      current_site_plot_locations %>% 
      filter(plot == paste(current_site, current_plot_id, sep = "_"))
    
    # We want enough canopy beyond the ground plots that we can avoid edge effects in our crown segmentation efforts during validation. That is, we want the crown segmentation that we perform on the cropped rasters just over each plot to essentially be the same crown segmentations that occurs when we work the whole site's CHM
    # So we buffer around the plot boundary by another 5 meters (from 11.35m radius to 16.35m radius)
    current_plot_boundary <-
      current_plot_location %>% 
      st_buffer(plot_radius + 5)
    
    # current_plot_dtm <- raster::crop(dtm_resamp, current_plot_boundary)
    # current_plot_dsm <- raster::crop(dsm, current_plot_boundary)
    current_plot_chm <- raster::crop(chm_smooth, current_plot_boundary)
    current_plot_ortho <- raster::crop(site_ortho, current_plot_boundary)
    
    current_site_las_catalog <- lidR::catalog(paste0(current_dir, current_site, "_vegetation-from-csf.las"))
    current_plot_las <- lidR::lasclip(current_site_las_catalog, as(current_plot_boundary, "Spatial")@polygons[[1]]@Polygons[[1]])
    lidR::lasnormalize(las = current_plot_las, dtm = site_dtm)
    
    # Write the cropped geospatial data to files for rapid recall and use in segmentation algorithm validation against ground data
    # chm
    writeRaster(x = current_plot_chm, filename = paste0(current_dir, current_site, "_plot-remote-data/", current_site, "_", current_plot_id, "_chm.tif"))
    # las
    writeLAS(.las = current_plot_las, file = paste0(current_dir, current_site, "_plot-remote-data/", current_site, "_", current_plot_id, "_normalized-point-cloud.las"))
    # ortho
    writeRaster(x = current_plot_ortho, filename = paste0(current_dir, current_site, "_plot-remote-data/", current_site, "_", current_plot_id, "_ortho.tif"))

  } # End loop working through the different plots within a site
} # End loop working through the different sites
