# This script will create small, cropped rasters of the Pix4D orthophoto outputs that overlap the plots (where the orange X's over the plot centers were visible)
# By saving the much smaller rasters as separte files, they will be much faster to work with and test different crown segmentation methods on the in situ data.
# Before chopping up the original orthophoto, each raster (of the 40ha area; ~625m x 625m) is ~600MB, but we only use between 1 and 5, ~25m x 25m sections of the orthophoto that overlap the ground plots. Performing multiple tree top detection and crown segmentation algorithms on the 144 plots with known coordinates in the orthophoto (because their orange X's were visible in the orthophoto) will be much faster using the smaller rasters.

library(sf)
library(tidyverse)
library(purrr)
library(raster)
library(lidR)

# Now there is an R object in the environment called "sites_checklist" that has
# infomation about how far along all processing steps are.

source("workflow/01_make-processing-checklist.R")

# Create directories if needed
if(!dir.exists("data/data_drone/L1/dsm/cropped-to-plot/")) {
  dir.create("data/data_drone/L1/dsm/cropped-to-plot/")
}

if(!dir.exists("data/data_drone/L1/ortho/cropped-to-plot/")) {
  dir.create("data/data_drone/L1/ortho/cropped-to-plot/")
}

if(!dir.exists("data/data_drone/L2/chm/cropped-to-plot/")) {
  dir.create("data/data_drone/L2/chm/cropped-to-plot/")
}

if(!dir.exists("data/data_drone/L2/classified-point-cloud/cropped-to-plot/")) {
  dir.create("data/data_drone/L2/classified-point-cloud/cropped-to-plot/")
}

if(!dir.exists("data/data_drone/L2/dtm/cropped-to-plot/")) {
  dir.create("data/data_drone/L2/dtm/cropped-to-plot/")
}

if(!dir.exists("data/data_drone/L2/index/cropped-to-plot/")) {
  dir.create("data/data_drone/L2/index/cropped-to-plot/")
}

# This is where I can put in sites that need their processing redone. An empty 
# string means that no already-processed site output will be overwritten
# (but sites that have yet to be processed will still have their processing done)
sites_to_overwrite <- ""
sites_checklist$overwrite <- ifelse(sites_to_overwrite == "all", yes = TRUE, no = FALSE)

sites_checklist[sites_checklist$site %in% sites_to_overwrite, "overwrite"] <- TRUE

sites_to_process <- 
  sites_checklist %>% 
  dplyr::filter(overwrite | (!L1_plot_remote_data_check & chm_check)) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()

(start <- Sys.time())

for (i in seq_along(sites_to_process)) {
  # get the character string representing the ith site
  current_site <- sites_to_process[i]
  
  # Import all of the products that we want to crop to the inidivdual ground plots
  # DTM, DSM, point cloud, index tifs, ortho tifs
  
  # The Digital Surface Model (dsm) is the ~5cm resolution raster representing
  # the surface (ground + objects on top) that the drone flew over
  dsm <- raster::raster(x = here::here(paste0("data/data_drone/L1/dsm/", current_site, "_dsm.tif")))
  
  # Orthomosaic (could be 5-band or 8-band depending on whether it includes the 
  # X3-derived R, G, and B bands)
  ortho <- raster::brick(here::here(paste0("data/data_drone/L1/ortho/", current_site, "_ortho.tif")))
  
  # The Digital Terrain Model (dtm) represents the "ground" underneath
  # the current site. Created using the cloth simulation filter implemented in 
  # the lidR package. See the "create-canopy-height-models.R" script.
  dtm <- raster::raster(here::here(paste0("data/data_drone/L2/dtm/", current_site, "_dtm.tif")))
  
  # The Canopy Height Model (chm) represents the height of the canopy above
  # the ground level. It is derived by subtracting the digital terrain
  # model from the digital surface model (chm = dsm - dtm). See the
  # "create-canopy-height-models.R" script
  chm <- raster::raster(here::here(paste0("data/data_drone/L2/chm/", current_site, "_chm.tif")))
  
  # Index mosaic (could be 6-band or 9-band depending on whether it includes the
  # X3-derived R, G, and B bands)
  index <- raster::brick(here::here(paste0("data/data_drone/L2/index/", current_site, "_index.tif")))
  
  # The densified point cloud can be read in and cropped using the LAScatalog functionality from {lidR}
  current_site_las_catalog <- lidR::readLAScatalog(here::here(paste0("data/data_drone/L2/classified-point-cloud/", current_site, "_classified-point-cloud.las")))
  
  plot_radius <- sqrt((66*66) / pi) * 12 * 2.54 / 100
  
  # Get plot locations as determined by visually inspecting orthophotos and finding orange X's laid
  # across the plot centers
  current_site_plot_locations <- sf::st_read(here::here(paste0("data/data_drone/L1/plot-locations/", current_site, "_plot-locations.gpkg")), stringsAsFactors = FALSE)
  
  for (j in seq_along(current_site_plot_locations$plot)) {
    current_plot <- current_site_plot_locations[j, ]
    
    # We want enough canopy beyond the ground plots that we can avoid edge effects 
    # in our crown segmentation efforts during validation. That is, we want the 
    # crown segmentation that we perform on the cropped rasters just over each plot 
    # to essentially be the same crown segmentations that occurs when we work the 
    # whole site's CHM. So we buffer around the plot boundary by another 5 meters 
    # (from 11.35m radius to 16.35m radius)
    current_plot_boundary <-
      current_plot %>% 
      st_buffer(plot_radius + 5)
    
    # Until patch is integrated into {lidR}, we have to rename the geometry column to be "geometry"
    current_plot_boundary <- st_sf(sf::st_drop_geometry(current_plot_boundary), geometry = sf::st_geometry(current_plot_boundary))
    
    current_plot_dtm <- raster::crop(dtm, current_plot_boundary)
    current_plot_dsm <- raster::crop(dsm, current_plot_boundary)
    current_plot_chm <- raster::crop(chm, current_plot_boundary)
    current_plot_ortho <- raster::crop(ortho, current_plot_boundary)
    current_plot_index <- raster::crop(index, current_plot_boundary)
    current_plot_las <- lidR::clip_roi(current_site_las_catalog, current_plot_boundary)
    
    # Write the cropped geospatial data to files for rapid recall and use in segmentation algorithm validation against ground data
    # Level 1
    # dsm
    writeRaster(x = current_plot_dsm, filename = here::here(paste0("data/data_drone/L1/dsm/cropped-to-plot/", current_plot$plot, "_dsm.tif")), overwrite = TRUE)
    # ortho
    writeRaster(x = current_plot_ortho, filename = here::here(paste0("data/data_drone/L1/ortho/cropped-to-plot/", current_plot$plot, "_ortho.tif")), overwrite = TRUE)

    # Level 2
    # dtm
    writeRaster(x = current_plot_dtm, filename = here::here(paste0("data/data_drone/L2/dtm/cropped-to-plot/", current_plot$plot, "_dtm.tif")), overwrite = TRUE)
    # chm
    writeRaster(x = current_plot_chm, filename = here::here(paste0("data/data_drone/L2/chm/cropped-to-plot/", current_plot$plot, "_chm.tif")), overwrite = TRUE)
    # las
    writeLAS(las = current_plot_las, file = here::here(paste0("data/data_drone/L2/classified-point-cloud/cropped-to-plot/", current_plot$plot, "_classified-point-cloud.las")))
    # index
    writeRaster(x = current_plot_index, filename = here::here(paste0("data/data_drone/L2/index/cropped-to-plot/", current_plot$plot, "_index.tif")), overwrite = TRUE)
  } # End loop working through the different plots within a site
  
  print(Sys.time() - start)
  print(paste0("...", current_site, " complete..."))
} # End loop working through the different sites

(end <- Sys.time())
(end - start)