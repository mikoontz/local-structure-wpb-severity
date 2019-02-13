# Load libraries
library(RCSF)
library(lidR)
library(sf)
library(tidyverse)
library(purrr)
library(raster)

source("data/data_carpentry/make_processing-checklist.R")

sites_checklist

unusable_sites <- c("eldo_4k_3", # too many blocks
                    "stan_4k_3", # too many blocks
                    "stan_5k_3", # too many blocks
                    "sequ_4k_2") # middle section flown on a separate day and the stitch looks terrible

usable_sites <- 
  sites_checklist %>% 
  dplyr::filter(!(site %in% unusable_sites)) %>% 
  dplyr::pull(site)

i = 1
for (i in seq_along(usable_sites)) {
  current_site <- usable_sites[i]
  # Wipe clean all of Pix4D classification
  
  # Get the file system right. Is it time to move all of the data_workign stuff to site_data??
  
  las <- readLAS(files = here::here(paste0("data/data_working/", current_site, "/", current_site, stan_3k_2_re/2_densification/point_cloud/stan_3k_2_re_Green_densified_point_cloud.las")
  las@data[, Classification := 0L]
  
}  

(start <- Sys.time())

# cloth resolution set to be ~5 times the average spacing between points
# when point cloud density is ~30 pts per m^2

las <- lasground(las, csf(sloop_smooth = TRUE, 
                          class_threshold = 0.25, 
                          cloth_resolution = 1.0,  
                          rigidness = 1, 
                          iterations = 500, 
                          time_step = 0.65))

plot(las, color = "Classification")

(end <- Sys.time())
print(end - start)

# DTM's partially using Pix4D point classification
dtm1a = grid_terrain(las = las,
                     res = 1.0,
                     algorithm = knnidw(k = 100, p = 2))

(start <- Sys.time())

raster::plot(dtm1a)

dtm1b = grid_terrain(las = las,
                     res = 1.0,
                     algorithm = tin())

dtm <- raster::raster("data/data_output/site_data/stan_3k_2/stan_3k_2_re_1m-dtm.tif")
dtm[dtm == -999] <- NA

raster::plot(dtm)
raster::plot(dtm1b)
?tin
