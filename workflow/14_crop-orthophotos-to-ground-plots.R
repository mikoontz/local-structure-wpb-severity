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

source("data/data_carpentry/make-processing-checklist.R")

unusable_sites <- c("eldo_4k_3", # too many blocks
                    "stan_4k_3", # too many blocks
                    "stan_5k_3", # too many blocks
                    "sequ_4k_2") # middle section flown on a separate day and the stitch looks terrible

G# These sites were processed with their X3 and RedEdge imagery combined so some of their
# output products will be in a slightly different place in the project directory
merged_sites <- c("eldo_3k_2",
                  "eldo_3k_3",
                  "eldo_4k_2")

# This is where I can put in sites that need their processing redone. An empty 
# string means that no already-processed site output will be overwritten
# (but sites that have yet to be processed will still have their processing done)
sites_to_overwrite <- ""
sites_checklist$overwrite <- ifelse(sites_to_overwrite == "all", yes = TRUE, no = FALSE)

sites_checklist[sites_checklist$site %in% sites_to_overwrite, "overwrite"] <- TRUE


sites_to_process <- 
  sites_checklist %>% 
  dplyr::filter(!(site %in% unusable_sites)) %>%
  dplyr::filter(overwrite | (!plot_remote_data_check & chm_check)) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()

(start <- Sys.time())

for (i in seq_along(sites_to_process)) {
  # get the character string representing the ith site
  current_site <- sites_to_process[i]
  
  current_dir <- paste0("data/data_output/site_data/", current_site, "/")
  
  if (sites_checklist$overwrite[sites_checklist$site == current_site]) {
    if (dir.exists(paste0(current_dir, current_site, "_plot-remote-data"))) {
      unlink(paste0(current_dir, current_site, "_plot-remote-data"), recursive = TRUE)
      message(paste("...Erasing", paste0(current_dir, current_site, "_plot-remote-data", "...")))
    }
  }
  
  # Create a new directory to put the plot-specific remotely sensed data
  dir.create(paste0(current_dir, current_site, "_plot-remote-data"))
  
  # Import all of the products that we want to crop to the inidivdual ground plots
  # DTM, DSM, point cloud, index tifs, ortho tifs
  
  # The Digital Terrain Model (dtm) represents the "ground" underneath
  # the current site. Created using the cloth simulation filter implemented in 
  # the lidR package. See the "create-canopy-height-models.R" script.
  dtm <- raster::raster(paste0(current_dir, current_site, "_dtm.tif"))
  
  # The Digital Surface Model (dsm) is the ~5cm resolution raster representing
  # the surface (ground + objects on top) that the drone flew over
  if (current_site %in% merged_sites) {
    dsm <- raster::raster(x = here::here(paste0("data/data_output/site_data/", current_site, "/", "3_dsm_ortho/1_dsm/", current_site, "_dsm.tif")))
  } else {
    dsm <- raster::raster(x = here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_re/3_dsm_ortho/1_dsm/", current_site, "_re_dsm.tif")))
  }
  
  # The Canopy Height Model (chm) represents the height of the canopy above
  # the ground level. It is derived by subtracting the digital terrain
  # model from the digital surface model (chm = dsm - dtm). See the
  # "create-canopy-height-models.R" script
  chm <- raster::raster(paste0(current_dir, current_site, "_chm.tif"))
  
  # Orthomosaic (could be 5-band or 8-band depending on whether it includes the 
  # X3-derived R, G, and B bands)
  ortho <- raster::brick(here::here(paste0(current_dir, current_site, "_ortho.tif")))
  
  # Index mosaic (could be 6-band or 9-band depending on whether it includes the
  # X3-derived R, G, and B bands)
  index <- raster::brick(here::here(paste0(current_dir, current_site, "_index.tif")))
  
  # # Make sure the Coordinate Reference Systems (crs) are the same
  # raster::crs(dtm) <- raster::crs(ortho)
  # 
  
  plot_radius <- sqrt((66*66) / pi) * 12 * 2.54 / 100
  
  # Get plot locations as determined by visually inspecting orthophotos and finding orange X's laid
  # across the plot centers
  if (current_site %in% merged_sites) {
    current_site_plot_locations <- 
      sf::st_read(paste0(current_dir, current_site, "_plot-locations/", current_site, "_plot-locations.shp")) %>% 
      mutate(plot = paste(current_site, id, sep = "_")) %>%
      dplyr::arrange(id) %>% 
      st_zm()
  } else {
    current_site_plot_locations <- 
      sf::st_read(paste0(current_dir, current_site, "_plot-locations_re/", current_site, "_plot-locations_re.shp")) %>% 
      mutate(plot = paste(current_site, id, sep = "_")) %>%
      dplyr::arrange(id) %>% 
      st_zm()
  }
  
  
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
    
    current_plot_dtm <- raster::crop(dtm, current_plot_boundary)
    current_plot_dsm <- raster::crop(dsm, current_plot_boundary)
    current_plot_chm <- raster::crop(chm, current_plot_boundary)
    current_plot_ortho <- raster::crop(ortho, current_plot_boundary)
    current_plot_index <- raster::crop(index, current_plot_boundary)
    
    current_site_las_catalog <- lidR::catalog(paste0(current_dir, current_site, "_classified_point_cloud.las"))
    current_plot_las <- lidR::lasclip(current_site_las_catalog, as(current_plot_boundary, "Spatial")@polygons[[1]]@Polygons[[1]])
    
    # Write the cropped geospatial data to files for rapid recall and use in segmentation algorithm validation against ground data
    # dsm
    writeRaster(x = current_plot_dsm, filename = paste0(current_dir, current_site, "_plot-remote-data/", pull(current_plot, plot), "_dsm.tif"), overwrite = TRUE)
    # dtm
    writeRaster(x = current_plot_dtm, filename = paste0(current_dir, current_site, "_plot-remote-data/", pull(current_plot, plot), "_dtm.tif"), overwrite = TRUE)
    # chm
    writeRaster(x = current_plot_chm, filename = paste0(current_dir, current_site, "_plot-remote-data/", pull(current_plot, plot), "_chm.tif"), overwrite = TRUE)
    # las
    writeLAS(las = current_plot_las, file = paste0(current_dir, current_site, "_plot-remote-data/", pull(current_plot, plot), "_classified-point-cloud.las"))
    # ortho
    writeRaster(x = current_plot_ortho, filename = paste0(current_dir, current_site, "_plot-remote-data/", pull(current_plot, plot), "_ortho.tif"), overwrite = TRUE)
    # index
    writeRaster(x = current_plot_index, filename = paste0(current_dir, current_site, "_plot-remote-data/", pull(current_plot, plot), "_index.tif"), overwrite = TRUE)
  } # End loop working through the different plots within a site
  
  print(Sys.time() - start)
  print(paste0("...", current_site, " complete..."))
} # End loop working through the different sites

(end <- Sys.time())
(end - start)