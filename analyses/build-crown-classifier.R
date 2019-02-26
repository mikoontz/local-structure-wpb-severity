# This script will build a classifier to predict the categorical response of 
# tree species from the reflectance data within each crown polygon.
# The classifier could then be used to classify all the trees in a given scene

# First step is to hand-classify a bunch of trees within a few scenes. We'll
# then build a model to predict the hand-classified species as a response
# using the reflectance data within that crown as predictors.

# We'll use the merged sites as a starting point, because they also include the
# more spatially resolved imagery from the X3 camera. That might make it easier
# to determine the species from the air.

# RandomForest
# Support Vector Machines


library(sf)
library(tidyverse)
library(purrr)
library(lidR)
library(viridis)
library(raster)
library(ForestTools)
library(gstat)
library(nngeo)
# devtools::install_github("Jean-Romain/lidRplugins")
library(lidRplugins)
library(furrr)
library(stars)
library(tictoc)
library(here)
library(units)

source(here::here("data/data_carpentry/make-processing-checklist.R"))

sites_checklist
# These sites had X3 and RedEdge photos merged into the same project, so we look in a different place for some of the relevant
# files.
merged_sites <- c("eldo_3k_2",
                  "eldo_3k_3",
                  "eldo_4k_2")

merged_sites %>%
  walk(.f = function(current_site) {
    crowns <- sf::st_read(dsn = here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_crowns/", current_site, "_crowns.shp")))
    
    if(!dir.exists(here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns")))) {
      dir.create(here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns")))
    }
    
    sf::st_write(obj = crowns, dsn = here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns/", current_site, "_hand-classified-crowns.shp")))
  })

other_sites_to_hand_classify <-
  c("eldo_3k_1")

other_sites_to_hand_classify %>%
  walk(.f = function(current_site) {
    crowns <- sf::st_read(dsn = here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_crowns/", current_site, "_crowns.shp")))
    
    if(!dir.exists(here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns")))) {
      dir.create(here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns")))
    }
    
    sf::st_write(obj = crowns, dsn = here::here(paste0("data/data_output/classified/hand-classified/", current_site, "_hand-classified-crowns/", current_site, "_hand-classified-crowns.shp")))
  })

# 
# 
# current_site <- "eldo_3k_1"
# sites_to_process <- current_site
# 
# ground_trees <- 
#   sites_to_process %>% 
#   map(.f = function(current_site) {
#     
#     if (current_site %in% merged_sites) {
#       current_site_plot_locations <- 
#         sf::st_read(paste0("data/data_output/site_data/", current_site, "/", current_site, "_plot-locations/", current_site, "_plot-locations.shp")) %>% 
#         mutate(plot = paste(current_site, id, sep = "_")) %>%
#         dplyr::arrange(id) %>% 
#         st_zm()
#     } else {
#       current_site_plot_locations <- 
#         sf::st_read(paste0("data/data_output/site_data/", current_site, "/", current_site, "_plot-locations_re/", current_site, "_plot-locations_re.shp")) %>% 
#         mutate(plot = paste(current_site, id, sep = "_")) %>%
#         dplyr::arrange(id) %>% 
#         st_zm()
#     }
#     
#     current_plot_ground_trees <- 
#       d %>% 
#       filter(site == current_site) %>% 
#       left_join(current_site_plot_locations, by = "plot") %>% 
#       dplyr::mutate(delta_x = cospi((1 / 2) - (azm / 180)) * dist) %>% 
#       dplyr::mutate(delta_y = sinpi((1 / 2) - (azm / 180)) * dist) %>%
#       st_as_sf() %>%
#       dplyr::mutate(x = st_coordinates(.)[, "X"],
#                     y = st_coordinates(.)[, "Y"]) %>% 
#       as.data.frame() %>% 
#       dplyr::mutate(new_x = x + delta_x,
#                     new_y = y + delta_y) %>% 
#       filter(!is.na(x)) %>% 
#       st_as_sf(coords = c("new_x", "new_y")) %>% 
#       st_set_crs(st_crs(current_site_plot_locations))
#     
#     nn <- 
#       st_nn(current_plot_ground_trees, current_plot_ground_trees, k = min(c(4, nrow(current_plot_ground_trees))), returnDist = TRUE, sparse = FALSE, progress = FALSE)$dist %>% 
#       as.data.frame()
#     
#     nn <- nn %>% dplyr::select(-1)
#     
#     if (ncol(nn) == 1) {nn <- nn %>% rename(nn1 = V2) %>% mutate(nn2 = NA, nn3 = NA)} else
#       if (ncol(nn) == 2) {nn <- nn %>% rename(nn1 = V2, nn2 = V3) %>% mutate(nn3 = NA)} else
#         if (ncol(nn) == 3) {nn <- nn %>% rename(nn1 = V2, nn2 = V3, nn3 = V4) }
#     
#     current_plot_ground_trees <- 
#       current_plot_ground_trees %>% 
#       bind_cols(nn)
#   }) %>%
#   do.call("rbind", .)
# 
# current_plot <- "eldo_3k_1_4"
# current_index <- raster::brick(here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_plot-remote-data/", current_plot, "_index.tif")))
# names(current_index) <- c("b", "g", "r", "re", "nir", "ndvi")
# 
# 
# 
# current_plot_trees <- 
#   ground_trees %>% 
#   dplyr::filter(is.na(year_fall)) %>% 
#   dplyr::filter(plot == current_plot)
# 
# # plot(current_index[["ndvi"]], col = viridis(100))
# # plotRGB(current_index, r = "r", g = "g", b = "b", scale = 0.22)
# plotRGB(current_index, r = "r", g = "g", b = "b", scale = 0.2)
# plot(current_plot_trees[, "species"], add = TRUE, pch = 19, pal = rainbow(5))
# 
# ggplot(current_plot_trees, aes(color = species)) +
#   geom_sf()
# 
# 
# 











unusable_sites <- c("eldo_4k_3", # too many blocks
                    "stan_4k_3", # too many blocks
                    "stan_5k_3", # too many blocks
                    "sequ_4k_2") # middle section flown on a separate day and the stitch looks terrible

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

sites_to_process <- 
  sites_checklist %>% 
  dplyr::filter(!(site %in% unusable_sites)) %>%
  dplyr::filter(overwrite | !ttops_check | !crowns_check) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()
