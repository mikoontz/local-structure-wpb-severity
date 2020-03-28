# Purpose is to stack the orthomosaics and the index mosaics for both the
# merged and the RedEdge-only data outputs from Pix4D. This will keep
# all the data conveniently stacked and available for cross-band calculations
# It also serves the purpose to re-home these data products from the
# default (and convoluted) Pix4D file structure to one that we can work with
# a bit more easily.

# Recall that the processing happened in a separate folder on a solid-state hard drive with the
# following file structure (see also: workflow/07_integrate-x3-and-re-photos-for-one-project.txt)

# H:/ (a solid-state local drive)
# |---- my_project/
#       |---- readme.md
#       |---- my_project.Rproj
#       |---- data/
#             |---- data_drone/
#                   |---- L0/
#                   |---- L1/
#                   |---- L2/
#                   |---- L3a/
#                   |---- L3b/
#                   |---- L4/
#             |---- data_working/
#                   |---- "current_site".p4d
#                   |---- "current_site"_re_photos/
#                         |---- calibration/
#                   |---- "current_site"_x3_photos/
#                   |____ "current_site"/
#                         |---- "current_site"_photos/
#                         |---- "current_site"_re/
#                               |---- 1_initial/
#                               |---- 2_densification/
#                               |---- 3_dsm_ortho/
#                               |____ 4_index/
#                         |---- "current_site"_x3/
#                               |---- 1_initial/
#                               |---- 2_densification/
#                               |____ 3_dsm_ortho/
#                         |---- "current_site"_re.p4d  
#                         |---- "current_site"_x3.p4d
#                         |---- 1_initial/
#                         |---- 2_densification/
#                         |---- 3_dsm_ortho/
#                         |____ 4_index/  

library(sf)
library(tidyverse)
library(purrr)
library(raster)
library(lidR)

# Now there is an R object in the environment called "sites_checklist" that has
# infomation about how far along all processing steps are.
source("workflow/01_make-processing-checklist.R")

# These sites were processed with their X3 and RedEdge imagery combined so some of their
# output products will be in a slightly different place in the project directory
merged_sites <- c("eldo_3k_2",
                  "eldo_3k_3",
                  "eldo_4k_2")

# The merged versus the unmerged sites will have different numbers of bands
# Because the merged sites will have the X3 imagery incorporated, there will
# be an extra 3 bands for the ortho and the index outputs (the R, G, and B
# from the X3 camera)
# The R, G, and B bands from the X3 images will always be the final three bands
# if they exist.
# For the RedEdge-derived products, the bands go in the order of wavelength,
# from shortest to longest (B, G, R, RE, NIR)
# There is one Pix4D derived index (NDVI), which will go after the NIR band
# for the index mosaic

# This is where I can put in sites that need their processing redone. An empty 
# string means that no already-processed site output will be overwritten
# (but sites that have yet to be processed will still have their processing done)
sites_to_overwrite <- ""
sites_checklist$overwrite <- FALSE

sites_checklist[sites_checklist$site %in% sites_to_overwrite, "overwrite"] <- TRUE

if (sites_to_overwrite == "all") {
  sites_checklist$overwrite <- TRUE
}

sites_to_process <- 
  sites_checklist %>% 
  dplyr::filter(overwrite | ((!stacked_ortho_check | !stacked_index_check))) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()

if(!dir.exists("data/data_drone/L1/ortho/")) {
  dir.create("data/data_drone/L1/ortho/")
}

if(!dir.exists("data/data_drone/L2/index/")) {
  dir.create("data/data_drone/L2/index/")
}

(start <- Sys.time())

for (i in seq_along(sites_to_process)) {
  # get the character string representing the ith site
  current_site <- sites_to_process[i]
  
  # The orthophoto mosaics for each site include some efforts by Pix4D to "retouch"
  # so their colors match to each other better.
  if (current_site %in% merged_sites) {
    ortho_blue <- raster::raster(here::here(paste0("data/data_working/", current_site, "/3_dsm_ortho/2_mosaic/", current_site, "_transparent_mosaic_blue.tif")))
    ortho_green <- raster::raster(here::here(paste0("data/data_working/", current_site, "/3_dsm_ortho/2_mosaic/", current_site, "_transparent_mosaic_green.tif")))
    ortho_red <- raster::raster(here::here(paste0("data/data_working/", current_site, "/3_dsm_ortho/2_mosaic/", current_site, "_transparent_mosaic_red.tif")))
    ortho_re <- raster::raster(here::here(paste0("data/data_working/", current_site, "/3_dsm_ortho/2_mosaic/", current_site, "_transparent_mosaic_red edge.tif")))
    ortho_nir <- raster::raster(here::here(paste0("data/data_working/", current_site, "/3_dsm_ortho/2_mosaic/", current_site, "_transparent_mosaic_nir.tif")))
    ortho_rgb <- raster::brick(here::here(paste0("data/data_working/", current_site, "/3_dsm_ortho/2_mosaic/", current_site, "_transparent_mosaic_group1.tif")))
    ortho <- raster::brick(ortho_blue, ortho_green, ortho_red, ortho_re, ortho_nir, ortho_rgb)  
  } else {
    ortho_blue <- raster::raster(x = here::here(paste0("data/data_working/", current_site, "/", current_site, "_re/3_dsm_ortho/2_mosaic/", current_site, "_re_transparent_mosaic_blue.tif")))
    ortho_green <- raster::raster(x = here::here(paste0("data/data_working/", current_site, "/", current_site, "_re/3_dsm_ortho/2_mosaic/", current_site, "_re_transparent_mosaic_green.tif")))
    ortho_red <- raster::raster(x = here::here(paste0("data/data_working/", current_site, "/", current_site, "_re/3_dsm_ortho/2_mosaic/", current_site, "_re_transparent_mosaic_red.tif")))
    ortho_re <- raster::raster(x = here::here(paste0("data/data_working/", current_site, "/", current_site, "_re/3_dsm_ortho/2_mosaic/", current_site, "_re_transparent_mosaic_red edge.tif")))
    ortho_nir <- raster::raster(x = here::here(paste0("data/data_working/", current_site, "/", current_site, "_re/3_dsm_ortho/2_mosaic/", current_site, "_re_transparent_mosaic_nir.tif")))
    ortho <- raster::brick(ortho_blue, ortho_green, ortho_red, ortho_re, ortho_nir)  
  }
  
  # The index mosaics for each site include no color correction by Pix4D and
  # (for the RedEdge imagery) also include corrections to the calibration panels
  if (current_site %in% merged_sites) {
    index_blue <- raster::raster(here::here(paste0("data/data_working/", current_site, "/4_index/reflectance/", current_site, "_transparent_reflectance_blue.tif")))
    index_green <- raster::raster(here::here(paste0("data/data_working/", current_site, "/4_index/reflectance/", current_site, "_transparent_reflectance_green.tif")))
    index_red <- raster::raster(here::here(paste0("data/data_working/", current_site, "/4_index/reflectance/", current_site, "_transparent_reflectance_red.tif")))
    index_re <- raster::raster(here::here(paste0("data/data_working/", current_site, "/4_index/reflectance/", current_site, "_transparent_reflectance_red edge.tif")))
    index_nir <- raster::raster(here::here(paste0("data/data_working/", current_site, "/4_index/reflectance/", current_site, "_transparent_reflectance_nir.tif")))
    index_ndvi <- raster::raster(here::here(paste0("data/data_working/", current_site, "/4_index/indices/ndvi/", current_site, "_index_ndvi.tif")))
    index_rgb <- raster::brick(here::here(paste0("data/data_working/", current_site, "/4_index/reflectance/", current_site, "_transparent_reflectance_group1.tif")))
    index <- raster::brick(index_blue, index_green, index_red, index_re, index_nir, index_ndvi, index_rgb)
  } else {
    index_blue <- raster::raster(x = here::here(paste0("data/data_working/", current_site, "/", current_site, "_re/4_index/reflectance/", current_site, "_re_transparent_reflectance_blue.tif")))
    index_green <- raster::raster(x = here::here(paste0("data/data_working/", current_site, "/", current_site, "_re/4_index/reflectance/", current_site, "_re_transparent_reflectance_green.tif")))
    index_red <- raster::raster(x = here::here(paste0("data/data_working/", current_site, "/", current_site, "_re/4_index/reflectance/", current_site, "_re_transparent_reflectance_red.tif")))
    index_re <- raster::raster(x = here::here(paste0("data/data_working/", current_site, "/", current_site, "_re/4_index/reflectance/", current_site, "_re_transparent_reflectance_red edge.tif")))
    index_nir <- raster::raster(x = here::here(paste0("data/data_working/", current_site, "/", current_site, "_re/4_index/reflectance/", current_site, "_re_transparent_reflectance_nir.tif")))
    index_ndvi <- raster::raster(x = here::here(paste0("data/data_working/", current_site, "/", current_site, "_re/4_index/indices/ndvi/", current_site, "_re_index_ndvi.tif")))
    index <- raster::brick(index_blue, index_green, index_red, index_re, index_nir, index_ndvi)
  }
  
  raster::writeRaster(x = ortho, 
                      filename = here::here(paste0("data/data_drone/L1/ortho/", current_site, "_ortho.tif")), 
                      overwrite = sites_checklist[sites_checklist$site == current_site, "overwrite"])
  raster::writeRaster(x = index, 
                      filename = here::here(paste0("data/data_drone/L2/index/", current_site, "_index.tif")), 
                      overwrite = sites_checklist[sites_checklist$site == current_site, "overwrite"])
  
  print(paste0("...", current_site, " complete..."))
}  

(end <- Sys.time())
(end - start)
