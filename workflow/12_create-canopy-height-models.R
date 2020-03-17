# Load packages
library(RCSF)
library(lidR)
library(sf)
library(tidyverse)
library(purrr)
library(raster)
library(rgl)

source("data/data_carpentry/make-processing-checklist.R")

sites_checklist

unusable_sites <- c("eldo_4k_3", # too many blocks
                    "stan_4k_3", # too many blocks
                    "stan_5k_3", # too many blocks
                    "sequ_4k_2") # middle section flown on a separate day and the stitch looks terrible

# These sites were processed with their X3 and RedEdge imagery combined so some of their
# output products will be in a slightly different place in the project directory
merged_sites <- c("eldo_3k_2",
                  "eldo_3k_3",
                  "eldo_4k_2")

# Set the Cloth Simulation Filter processing parameters for different sites
# These parameters are the defaults
csf_parameters <- data_frame(site = sites_checklist$site,
                             sloop_smooth = TRUE,
                             class_threshold = 0.25,
                             cloth_resolution = 1.0,
                             rigidness = 1,
                             iterations = 500,
                             time_step = 0.65)

# Sites that need a slightly finer resolution cloth (found through trial and error)
cloth_res_0.9_sites <- c("eldo_5k_2", "sequ_5k_1", "sequ_5k_2")
csf_parameters[csf_parameters$site %in% cloth_res_0.9_sites, "cloth_resolution"] <- 0.9

cloth_res_0.75_sites <- c("eldo_3k_1", "stan_3k_1", "sequ_6k_1")
csf_parameters[csf_parameters$site %in% cloth_res_0.75_sites, "cloth_resolution"] <- 0.75

cloth_res_0.6_sites <- c("stan_3k_2", "sier_3k_2", "sier_3k_1", "sier_4k_2", "sier_5k_1", "sier_5k_3")
csf_parameters[csf_parameters$site %in% cloth_res_0.6_sites, "cloth_resolution"] <- 0.6

cloth_res_0.5_sites <- c("eldo_4k_1", "stan_3k_3", "sier_3k_3", "sequ_4k_1", "sequ_6k_2")
csf_parameters[csf_parameters$site %in% cloth_res_0.5_sites, "cloth_resolution"] <- 0.5

cloth_res_0.4_sites <- c("eldo_5k_1", "stan_4k_1", "sequ_4k_3", "sequ_6k_3")
csf_parameters[csf_parameters$site %in% cloth_res_0.4_sites, "cloth_resolution"] <- 0.4

sites_to_overwrite <- ""
sites_checklist$overwrite <- FALSE

sites_checklist[sites_checklist$site %in% sites_to_overwrite, "overwrite"] <- TRUE

sites_to_process <- 
  sites_checklist %>% 
  dplyr::filter(!(site %in% unusable_sites)) %>%
  dplyr::filter(overwrite | !classified_point_cloud_check | !dtm_check | !chm_check) %>% 
  dplyr::pull(site)


(start <- Sys.time())

for (i in seq_along(sites_to_process)) {
  current_site <- sites_to_process[i]
  
  if (current_site %in% merged_sites) {
    current_point_cloud <- lidR::readLAS(files = here::here(paste0("data/data_output/site_data/", current_site, "/", "2_densification/point_cloud/", current_site, "_densified_point_cloud.las")))
  } else {
    current_point_cloud <- lidR::readLAS(files = here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_re/2_densification/point_cloud/", current_site, "_re_Green_densified_point_cloud.las")))
  }
  
  
  # Wipe clean all of Pix4D classification
  current_point_cloud@data[, Classification := 0L]

  # Use the cloth simulation filter by Zhang et al. (2016) [http://www.mdpi.com/2072-4292/8/6/501/htm]
  # implemented in the lidR package to classify points in the point cloud as ground or non-ground
  
  # cloth resolution set to be ~5 times the average spacing between points
  # when point cloud density is ~30 pts per m^2
  current_idx <- csf_parameters$site == current_site
  
  current_point_cloud <- lidR::lasground(las = current_point_cloud, 
                                         algorithm = csf(sloop_smooth = csf_parameters$sloop_smooth[current_idx], 
                                                         class_threshold = csf_parameters$class_threshold[current_idx], 
                                                         cloth_resolution = csf_parameters$cloth_resolution[current_idx],  
                                                         rigidness = csf_parameters$rigidness[current_idx], 
                                                         iterations = csf_parameters$iterations[current_idx], 
                                                         time_step = csf_parameters$time_step[current_idx]))
  
  
  # Plot the classification of the point cloud for inspection
  plot(current_point_cloud, color = "Classification")
  legend3d("topright", legend = current_site, bty = "n")
  
  
  # Export the classified point cloud to disk so we can use the vegetation points for tree segmentation if
  # we want
  lidR::writeLAS(las = current_point_cloud, file = here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_classified-point-cloud.las")))
  
  
  # Create a 1m resolution digital terrain model using the classified ground points
  # and interpolation between those ground points using the tin() function from lidR which 
  # implements a Delaunay triangulation method
  dtm <- lidR::grid_terrain(las = current_point_cloud,
                            res = 1.0,
                            algorithm = tin())
  
  
  # Write the dtm file to disk so we can use it later.
  # Note that all of these outputs generated using R get written to the same place, regardless of whether the
  # output is derived from merged X3+RedEdge imagery versus just being derived from RedEdge imagery
  raster::writeRaster(x = dtm, filename = here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_dtm.tif")), overwrite = TRUE)
  
  # Import the digital surface model Pix4D output in order to generate a canopy height model
  if (current_site %in% merged_sites) {
    dsm <- raster::raster(x = here::here(paste0("data/data_output/site_data/", current_site, "/", "3_dsm_ortho/1_dsm/", current_site, "_dsm.tif")))
  } else {
    dsm <- raster::raster(x = here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_re/3_dsm_ortho/1_dsm/", current_site, "_re_dsm.tif")))
  }

  
  # Using bilinear interpolation to downsample the 1m resolution DTM to have the
  # same resolution as the dsm (~5cm, but slightly different for each site)
  dtm_resamp <- raster::resample(dtm, dsm, method = "bilinear")
  
  # The Canopy Height Model (chm) is the dsm (vegetation + ground) minus the dtm (ground)
  # to give just the height of the vegetation.
  chm <- dsm - dtm_resamp
  
  # Need to decide whether to do CHM smoothing now or just before the segmentation process.
  # # Smooth out the chm and set any negative values to 0 (meaning "ground") following
  # # advice from Zagalikis, Cameron, and Miller (2004) and references therein
  # # More recently, a 3x3 pixel smoothing filter was specifically suggested as ideal
  # # for sUAS derived chm by Mohan et al. (2017)
  # chm_smooth <- raster::focal(chm, w = matrix(1, 3, 3), mean)
  # chm_smooth[raster::getValues(chm_smooth) < 0] <- 0
  
  # Write the chm file to disk so we can use it later
  # Note that all of these outputs generated using R get written to the same place, regardless of whether the
  # output is derived from merged X3+RedEdge imagery versus just being derived from RedEdge imagery
  raster::writeRaster(x = chm, filename = here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_chm.tif")), overwrite = TRUE)

  print(paste0("...", current_site, " complete..."))
}  

(end <- Sys.time())
(end - start)

# # Testing the snag detection algorithms. Might be another way to find dead trees.
# bbpr_thresholds <- matrix(c(0.80, 0.80, 0.70,
#                             0.85, 0.85, 0.60,
#                             0.80, 0.80, 0.60,
#                             0.90, 0.90, 0.55),
#                           nrow =3, ncol = 4)
# 
# # Run snag classification and assign classes to each point
# snags <- lassnags(las = current_point_cloud, algorithm = wing2015(neigh_radii = c(1.5, 1, 2), BBPRthrsh_mat = bbpr_thresholds))
