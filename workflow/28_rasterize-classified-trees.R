# Purpose: rasterize the classified trees to 20m x 20m cells in order to mimic the footprint of the field plots

library(tidyverse)
library(raster)
library(sf)
library(velox)
library(here)
library(tictoc)
library(lubridate)
library(nngeo)

# Now there is an R object in the environment called "sites_checklist" that has
# infomation about how far along all processing steps are.
source("workflow/01_make-processing-checklist.R")

# Converts a per-cell value (20 x 20m) to a per acre value
perCell_to_perAc <- function(r) {
  r / (prod(res(r)) / 4046.856)
} # returns trees per acre when r represents a raster with counts of trees per cell and raster has units of meters

perCell_to_perHa <- function(r) {
  r / (prod(res(r)) / 10000)
} # returns trees per hectare when r represents a raster with counts of trees per cell


if(!file.exists(here::here("data/data_drone/L3b/model-classified-trees_all.gpkg"))) {
  
  source("workflow/24_classify-all-trees.R")
  
} 

classified_trees <- 
  sf::st_read(here::here("data/data_drone/L3b/model-classified-trees_all.gpkg"))

# This is where I can put in sites that need their processing redone. An empty 
# string means that no already-processed site output will be overwritten
# (but sites that have yet to be processed will still have their processing done)
sites_to_overwrite <- "all"
sites_checklist$overwrite <- ifelse(sites_to_overwrite == "all", yes = TRUE, no = FALSE)

sites_checklist[sites_checklist$site %in% sites_to_overwrite, "overwrite"] <- TRUE

sites_to_process <- 
  sites_checklist %>% 
  dplyr::filter(overwrite | !rasterized_trees_check) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()

surveyed_area <- 
  sf::st_read("data/data_drone/L0/surveyed-area-3310.gpkg")

if(!dir.exists(here::here(paste0("data/data_drone/L4/rasterized-trees/")))) {
  dir.create(here::here(paste0("data/data_drone/L4/rasterized-trees/")), recursive = TRUE)
}

results_list <- vector(mode = "list", length = length(sites_to_process))

for(i in seq_along(sites_to_process)) {
  
  current_site <- sites_to_process[i]

  buffered_bounds <-
    surveyed_area %>% 
    dplyr::filter(site == current_site) %>% 
    st_buffer(-35) # buffer in a little further than the trees were to reduce edge effects (-35 worked)
  
  raster_template <- raster::raster(buffered_bounds, res = 20)
  
  # Get the classified trees for this particular site
  # Calculate the voronoi polygons around each tree (and the voronoi polygon area)
  current_trees <- 
    classified_trees %>%
    dplyr::filter(site == current_site)
  
  current_trees <-
    current_trees %>% 
    st_geometry() %>% 
    st_union() %>% 
    st_voronoi() %>% 
    st_cast() %>% 
    st_intersection(st_geometry(buffered_bounds)) %>% 
    sf::st_sf(geometry = .) %>% 
    dplyr::mutate(voronoi_area = st_area(.)) %>% 
    dplyr::mutate(voronoi_poly = geometry) %>% 
    st_join(current_trees, .) %>% 
    dplyr::mutate(nn = nngeo::st_nn(., ., k = 4, returnDist = TRUE, sparse = FALSE, progress = FALSE)$dist) %>% 
    dplyr::mutate(nn1 = map_dbl(nn, magrittr::extract, 2),
                  nn2 = map_dbl(nn, magrittr::extract, 3),
                  nn3 = map_dbl(nn, magrittr::extract, 4)) %>% 
    dplyr::select(-nn)
  
  current_cwd <- raster::resample(cwd, raster_template, method = "bilinear")

  # See the workflow/03_extract-cwd-from-locations.R script to get this vector file
  sn_pipo <- sf::st_read("data/data_output/sierra-nevada-pipo-cwd.gpkg")
  
  mean_cwd_sn_pipo <- mean(sn_pipo$cwd, na.rm = TRUE)
  sd_cwd_sn_pipo <- sd(sn_pipo$cwd, na.rm = TRUE)
  
  current_cwd_zscore <- (current_cwd - mean_cwd_sn_pipo) / sd_cwd_sn_pipo
  
  # # count of trees per cell -------------------------------------------------
  
  live_count <- raster::rasterize(x = current_trees, 
                                  y = raster_template, 
                                  field = "live", 
                                  background = 0,
                                  fun = function(x, ...) {
                                    length(which(na.omit(x) == 1))
                                  })
  
  dead_count <- raster::rasterize(x = current_trees, 
                                  y = raster_template, 
                                  field = "live", 
                                  background = 0, 
                                  fun = function(x, ...) {
                                    length(which(na.omit(x) != 1))
                                  })
  
  
  pipo_count <- raster::rasterize(x = current_trees, 
                                  y = raster_template, 
                                  field = "species",
                                  background = 0, 
                                  fun = function(x, ...) {
                                    length(which(na.omit(x) == "pipo"))
                                  })
  
  non_pipo_count <- raster::rasterize(x = current_trees, 
                                      y = raster_template, 
                                      field = "species",
                                      background = 0, 
                                      fun = function(x, ...) {
                                        length(which(na.omit(x) != "pipo"))
                                      })
  
  pipo_and_dead_count <- dead_count + pipo_count
  
  total_count <- raster::rasterize(x = current_trees, 
                                   y = raster_template, 
                                   field = "live", 
                                   background = 0, 
                                   fun = function(x, ...) {
                                     length(na.omit(x))
                                   })
  
  # r1 <- live_count + dead_count
  # compareRaster(r1, total_count, values = TRUE)
  # 
  # r2 <- pipo_count + non_pipo_count
  # compareRaster(r2, live_count, values = TRUE)
  # 
  
  # convert counts to trees per hectare -------------------------------------
  
  live_tpha <- perCell_to_perHa(live_count)  
  dead_tpha <- perCell_to_perHa(dead_count)
  pipo_tpha <- perCell_to_perHa(pipo_count)
  non_pipo_tpha <- perCell_to_perHa(non_pipo_count)
  pipo_and_dead_tpha <- perCell_to_perHa(pipo_and_dead_count)
  overall_tpha <- perCell_to_perHa(total_count)
  

  # mean heights of trees per cell ------------------------------------------

  live_mean_height <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 1)), 
                                                           y = raster_template, 
                                                           field = "height", 
                                                           background = NA, 
                                                           fun = mean)
  
  dead_mean_height <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 0)), 
                                       y = raster_template, 
                                       field = "height", 
                                       background = NA, 
                                       fun = mean)
  
  pipo_mean_height <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo")), 
                                       y = raster_template, 
                                       field = "height", 
                                       background = NA, 
                                       fun = mean)
  
  non_pipo_mean_height <- raster::rasterize(x = current_trees %>% dplyr::filter((species != "pipo")), 
                                           y = raster_template, 
                                           field = "height", 
                                           background = NA, 
                                           fun = mean)
  
  pipo_and_dead_mean_height <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 0) | (species == "pipo")), 
                                                 y = raster_template, 
                                                 field = "height", 
                                                 background = NA, 
                                                 fun = mean)
  
  overall_mean_height <- raster::rasterize(x = current_trees, 
                                        y = raster_template, 
                                        field = "height", 
                                        background = NA, 
                                        fun = mean)
  
  # total basal area per cell -----------------------------------------------
  
  live_basal_area <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 1)), 
                                       y = raster_template, 
                                       field = "estimated_ba", 
                                       background = NA, 
                                       fun = sum)
  
  dead_basal_area <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 0)), 
                                       y = raster_template, 
                                       field = "estimated_ba", 
                                       background = NA, 
                                       fun = sum)
  
  pipo_basal_area <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo")), 
                                       y = raster_template, 
                                       field = "estimated_ba", 
                                       background = NA, 
                                       fun = sum)
  
  non_pipo_basal_area <- raster::rasterize(x = current_trees %>% dplyr::filter((species != "pipo")), 
                                           y = raster_template, 
                                           field = "estimated_ba", 
                                           background = NA, 
                                           fun = sum)
  
  pipo_and_dead_basal_area <- dead_basal_area + pipo_basal_area
  
  total_basal_area <- raster::rasterize(x = current_trees, 
                                        y = raster_template, 
                                        field = "estimated_ba", 
                                        background = NA, 
                                        fun = sum)
  
  # r3 <- live_basal_area + dead_basal_area
  # compareRaster(r3, total_basal_area, values = TRUE)
  # 
  # r4 <- pipo_basal_area + non_pipo_basal_area
  # compareRaster(r4, live_basal_area, values = TRUE)
  
  # convert to basal area per hectare ---------------------------------------
  
  # We can actually use the same equations as for trees per hectare
  
  live_bapha <- perCell_to_perHa(live_basal_area)  
  dead_bapha <- perCell_to_perHa(dead_basal_area)
  pipo_bapha <- perCell_to_perHa(pipo_basal_area)
  non_pipo_bapha <- perCell_to_perHa(non_pipo_basal_area)
  pipo_and_dead_bapha <- perCell_to_perHa(pipo_and_dead_basal_area)
  overall_bapha <- perCell_to_perHa(total_basal_area)
  
  # mean basal area per tree ------------------------------------------------
  
  
  live_mean_ba <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 1)), 
                                    y = raster_template, 
                                    field = "estimated_ba", 
                                    background = NA, 
                                    fun = mean)
  
  dead_mean_ba <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 0)), 
                                    y = raster_template, 
                                    field = "estimated_ba", 
                                    background = NA, 
                                    fun = mean)
  
  pipo_mean_ba <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo")), 
                                    y = raster_template, 
                                    field = "estimated_ba", 
                                    background = NA, 
                                    fun = mean)
  
  non_pipo_mean_ba <- raster::rasterize(x = current_trees %>% dplyr::filter((species != "pipo")), 
                                        y = raster_template, 
                                        field = "estimated_ba", 
                                        background = NA, 
                                        fun = mean)
  
  pipo_and_dead_mean_ba <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo") | (live == 0)), 
                                             y = raster_template, 
                                             field = "estimated_ba", 
                                             background = NA, 
                                             fun = mean)
  
  overall_mean_ba <- raster::rasterize(x = current_trees, 
                                       y = raster_template, 
                                       field = "estimated_ba", 
                                       background = NA, 
                                       fun = mean)
  

# mean diameter -----------------------------------------------------------

  live_mean_dbh <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 1)), 
                                    y = raster_template, 
                                    field = "estimated_dbh", 
                                    background = NA, 
                                    fun = mean)
  
  dead_mean_dbh <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 0)), 
                                    y = raster_template, 
                                    field = "estimated_dbh", 
                                    background = NA, 
                                    fun = mean)
  
  pipo_mean_dbh <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo")), 
                                    y = raster_template, 
                                    field = "estimated_dbh", 
                                    background = NA, 
                                    fun = mean)
  
  non_pipo_mean_dbh <- raster::rasterize(x = current_trees %>% dplyr::filter((species != "pipo")), 
                                        y = raster_template, 
                                        field = "estimated_dbh", 
                                        background = NA, 
                                        fun = mean)
  
  pipo_and_dead_mean_dbh <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo") | (live == 0)), 
                                             y = raster_template, 
                                             field = "estimated_dbh", 
                                             background = NA, 
                                             fun = mean)
  
  overall_mean_dbh <- raster::rasterize(x = current_trees, 
                                       y = raster_template, 
                                       field = "estimated_dbh", 
                                       background = NA, 
                                       fun = mean)  
  
  # quadratic mean diameter -------------------------------------------------
  
  live_qmd <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 1)), 
                                y = raster_template, 
                                field = "estimated_dbh", 
                                background = NA, 
                                fun = function(x, ...) {sqrt(sum(x^2) / length(x))})
  
  dead_qmd <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 0)), 
                                y = raster_template, 
                                field = "estimated_dbh", 
                                background = NA, 
                                fun = function(x, ...) {sqrt(sum(x^2) / length(x))})
  
  pipo_qmd <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo")), 
                                y = raster_template, 
                                field = "estimated_dbh", 
                                background = NA, 
                                fun = function(x, ...) {sqrt(sum(x^2) / length(x))})
  
  non_pipo_qmd <- raster::rasterize(x = current_trees %>% dplyr::filter((species != "pipo")), 
                                    y = raster_template, 
                                    field = "estimated_dbh", 
                                    background = NA, 
                                    fun = function(x, ...) {sqrt(sum(x^2) / length(x))})
  
  pipo_and_dead_qmd <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo") | (live == 0)), 
                                         y = raster_template, 
                                         field = "estimated_dbh", 
                                         background = NA, 
                                         fun = function(x, ...) {sqrt(sum(x^2) / length(x))})
  
  overall_qmd <- raster::rasterize(x = current_trees, 
                                   y = raster_template, 
                                   field = "estimated_dbh", 
                                   background = NA, 
                                   fun = function(x, ...) {sqrt(sum(x^2) / length(x))})  
  
  
  # stand density index -----------------------------------------------------
  # A transformation to the equivalent number of trees with a 10 inch diameter
  # Original by Reineke 1933 (if qmd is in centimters)
  # sdi = tpa * ( qmd / 25.4)^1.605
  
  # Proposed constant for ponderosa in California by Oliver 1995
  # (if qmd is in centimeters)
  # sdi = tpa * ( qmd / 25.4)^1.77
  
  live_sdi_ha <- perCell_to_perHa(live_count) * (live_qmd / 25.4)^1.77
  live_sdi_ac <- perCell_to_perAc(live_count) * (live_qmd / 25.4)^1.77
  
  dead_sdi_ha <- perCell_to_perHa(dead_count) * (dead_qmd / 25.4)^1.77
  dead_sdi_ac <- perCell_to_perAc(dead_count) * (dead_qmd / 25.4)^1.77
  
  pipo_sdi_ha <- perCell_to_perHa(pipo_count) * (pipo_qmd / 25.4)^1.77
  pipo_sdi_ac <- perCell_to_perAc(pipo_count) * (pipo_qmd / 25.4)^1.77
  
  non_pipo_sdi_ha <- perCell_to_perHa(non_pipo_count) * (non_pipo_qmd / 25.4)^1.77
  non_pipo_sdi_ac <- perCell_to_perAc(non_pipo_count) * (non_pipo_qmd / 25.4)^1.77
  
  pipo_and_dead_sdi_ha <- perCell_to_perHa(pipo_and_dead_count) * (pipo_and_dead_qmd / 25.4)^1.77
  pipo_and_dead_sdi_ac <- perCell_to_perAc(pipo_and_dead_count) * (pipo_and_dead_qmd / 25.4)^1.77
  
  overall_sdi_ha <- perCell_to_perHa(total_count) * (overall_qmd / 25.4)^1.77
  overall_sdi_ac <- perCell_to_perAc(total_count) * (overall_qmd / 25.4)^1.77
  

# mean voronoi areas -----------------------------------------------------------

  live_mean_voronoi <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 1)), 
                                     y = raster_template, 
                                     field = "voronoi_area", 
                                     background = NA, 
                                     fun = mean)
  
  dead_mean_voronoi <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 0)), 
                                     y = raster_template, 
                                     field = "voronoi_area", 
                                     background = NA, 
                                     fun = mean)
  
  pipo_mean_voronoi <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo")), 
                                     y = raster_template, 
                                     field = "voronoi_area", 
                                     background = NA, 
                                     fun = mean)
  
  non_pipo_mean_voronoi <- raster::rasterize(x = current_trees %>% dplyr::filter((species != "pipo")), 
                                         y = raster_template, 
                                         field = "voronoi_area", 
                                         background = NA, 
                                         fun = mean)
  
  pipo_and_dead_mean_voronoi <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo") | (live == 0)), 
                                              y = raster_template, 
                                              field = "voronoi_area", 
                                              background = NA, 
                                              fun = mean)
  
  overall_mean_voronoi <- raster::rasterize(x = current_trees, 
                                        y = raster_template, 
                                        field = "voronoi_area", 
                                        background = NA, 
                                        fun = mean)  
  

# mean distance to first nearest neighbor --------------------------------------------------

  live_mean_nn1 <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 1)), 
                                         y = raster_template, 
                                         field = "nn1", 
                                         background = NA, 
                                         fun = mean)
  
  dead_mean_nn1 <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 0)), 
                                         y = raster_template, 
                                         field = "nn1", 
                                         background = NA, 
                                         fun = mean)
  
  pipo_mean_nn1 <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo")), 
                                         y = raster_template, 
                                         field = "nn1", 
                                         background = NA, 
                                         fun = mean)
  
  non_pipo_mean_nn1 <- raster::rasterize(x = current_trees %>% dplyr::filter((species != "pipo")), 
                                             y = raster_template, 
                                             field = "nn1", 
                                             background = NA, 
                                             fun = mean)
  
  pipo_and_dead_mean_nn1 <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo") | (live == 0)), 
                                                  y = raster_template, 
                                                  field = "nn1", 
                                                  background = NA, 
                                                  fun = mean)
  
  overall_mean_nn1 <- raster::rasterize(x = current_trees, 
                                            y = raster_template, 
                                            field = "nn1", 
                                            background = NA, 
                                            fun = mean)    
  
  # mean distance to second nearest neighbor --------------------------------------------------
  
  live_mean_nn2 <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 1)), 
                                     y = raster_template, 
                                     field = "nn2", 
                                     background = NA, 
                                     fun = mean)
  
  dead_mean_nn2 <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 0)), 
                                     y = raster_template, 
                                     field = "nn2", 
                                     background = NA, 
                                     fun = mean)
  
  pipo_mean_nn2 <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo")), 
                                     y = raster_template, 
                                     field = "nn2", 
                                     background = NA, 
                                     fun = mean)
  
  non_pipo_mean_nn2 <- raster::rasterize(x = current_trees %>% dplyr::filter((species != "pipo")), 
                                         y = raster_template, 
                                         field = "nn2", 
                                         background = NA, 
                                         fun = mean)
  
  pipo_and_dead_mean_nn2 <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo") | (live == 0)), 
                                              y = raster_template, 
                                              field = "nn2", 
                                              background = NA, 
                                              fun = mean)
  
  overall_mean_nn2 <- raster::rasterize(x = current_trees, 
                                        y = raster_template, 
                                        field = "nn2", 
                                        background = NA, 
                                        fun = mean)     
  
  # mean distance to third nearest neighbor --------------------------------------------------
  
  live_mean_nn3 <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 1)), 
                                     y = raster_template, 
                                     field = "nn3", 
                                     background = NA, 
                                     fun = mean)
  
  dead_mean_nn3 <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 0)), 
                                     y = raster_template, 
                                     field = "nn3", 
                                     background = NA, 
                                     fun = mean)
  
  pipo_mean_nn3 <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo")), 
                                     y = raster_template, 
                                     field = "nn3", 
                                     background = NA, 
                                     fun = mean)
  
  non_pipo_mean_nn3 <- raster::rasterize(x = current_trees %>% dplyr::filter((species != "pipo")), 
                                         y = raster_template, 
                                         field = "nn3", 
                                         background = NA, 
                                         fun = mean)
  
  pipo_and_dead_mean_nn3 <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo") | (live == 0)), 
                                              y = raster_template, 
                                              field = "nn3", 
                                              background = NA, 
                                              fun = mean)
  
  overall_mean_nn3 <- raster::rasterize(x = current_trees, 
                                        y = raster_template, 
                                        field = "nn3", 
                                        background = NA, 
                                        fun = mean)    
  
  
  # standard deviation distance to first nearest neighbor --------------------------------------------------
  
  live_sd_nn1 <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 1)), 
                                     y = raster_template, 
                                     field = "nn1", 
                                     background = NA, 
                                     fun = sd)
  
  dead_sd_nn1 <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 0)), 
                                     y = raster_template, 
                                     field = "nn1", 
                                     background = NA, 
                                     fun = sd)
  
  pipo_sd_nn1 <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo")), 
                                     y = raster_template, 
                                     field = "nn1", 
                                     background = NA, 
                                     fun = sd)
  
  non_pipo_sd_nn1 <- raster::rasterize(x = current_trees %>% dplyr::filter((species != "pipo")), 
                                         y = raster_template, 
                                         field = "nn1", 
                                         background = NA, 
                                         fun = sd)
  
  pipo_and_dead_sd_nn1 <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo") | (live == 0)), 
                                              y = raster_template, 
                                              field = "nn1", 
                                              background = NA, 
                                              fun = sd)
  
  overall_sd_nn1 <- raster::rasterize(x = current_trees, 
                                        y = raster_template, 
                                        field = "nn1", 
                                        background = NA, 
                                        fun = sd)    
  
 # coefficient of variation distance to first nearest neighbor --------------------------------------------------
  
  live_cov_nn1 <- live_sd_nn1 / live_mean_nn1
    
  dead_cov_nn1 <- dead_sd_nn1 / dead_mean_nn1
  
  pipo_cov_nn1 <- pipo_sd_nn1 / pipo_mean_nn1
  
  non_pipo_cov_nn1 <- non_pipo_sd_nn1 / non_pipo_mean_nn1
  
  pipo_and_dead_cov_nn1 <- pipo_and_dead_sd_nn1 / pipo_and_dead_mean_nn1
  
  overall_cov_nn1 <- overall_sd_nn1 / overall_mean_nn1  
  
  

# standard deviation of voronoi polygon areas -----------------------------

  live_sd_voronoi <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 1)), 
                                         y = raster_template, 
                                         field = "voronoi_area", 
                                         background = NA, 
                                         fun = sd)
  
  dead_sd_voronoi <- raster::rasterize(x = current_trees %>% dplyr::filter((live == 0)), 
                                         y = raster_template, 
                                         field = "voronoi_area", 
                                         background = NA, 
                                         fun = sd)
  
  pipo_sd_voronoi <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo")), 
                                         y = raster_template, 
                                         field = "voronoi_area", 
                                         background = NA, 
                                         fun = sd)
  
  non_pipo_sd_voronoi <- raster::rasterize(x = current_trees %>% dplyr::filter((species != "pipo")), 
                                             y = raster_template, 
                                             field = "voronoi_area", 
                                             background = NA, 
                                             fun = sd)
  
  pipo_and_dead_sd_voronoi <- raster::rasterize(x = current_trees %>% dplyr::filter((species == "pipo") | (live == 0)), 
                                                  y = raster_template, 
                                                  field = "voronoi_area", 
                                                  background = NA, 
                                                  fun = sd)
  
  overall_sd_voronoi <- raster::rasterize(x = current_trees, 
                                            y = raster_template, 
                                            field = "voronoi_area", 
                                            background = NA, 
                                            fun = sd)    
  

# coefficient of variation voronoi area -----------------------------------
  
  live_cov_voronoi <- live_sd_voronoi / live_mean_voronoi
  
  dead_cov_voronoi <- dead_sd_voronoi / dead_mean_voronoi
  
  pipo_cov_voronoi <- pipo_sd_voronoi / pipo_mean_voronoi
  
  non_pipo_cov_voronoi <- non_pipo_sd_voronoi / non_pipo_mean_voronoi
  
  pipo_and_dead_cov_voronoi <- pipo_and_dead_sd_voronoi / pipo_and_dead_mean_voronoi
  
  overall_cov_voronoi <- overall_sd_voronoi / overall_mean_voronoi 
  
# stack results together --------------------------------------------------

  
  results_raster <- raster::stack(live_count, dead_count, pipo_count, non_pipo_count, pipo_and_dead_count, total_count, 
                                  live_tpha, dead_tpha, pipo_tpha, non_pipo_tpha, pipo_and_dead_tpha, overall_tpha, 
                                  live_mean_height, dead_mean_height, pipo_mean_height, non_pipo_mean_height, pipo_and_dead_mean_height, overall_mean_height,
                                  live_basal_area, dead_basal_area, pipo_basal_area, non_pipo_basal_area, pipo_and_dead_basal_area, total_basal_area,
                                  live_bapha, dead_bapha, pipo_bapha, non_pipo_bapha, pipo_and_dead_bapha, overall_bapha,
                                  live_mean_ba, dead_mean_ba, pipo_mean_ba, non_pipo_mean_ba, pipo_and_dead_mean_ba, overall_mean_ba,
                                  live_mean_dbh, dead_mean_dbh, pipo_mean_dbh, non_pipo_mean_dbh, pipo_and_dead_mean_dbh, overall_mean_dbh,
                                  live_qmd, dead_qmd, pipo_qmd, non_pipo_qmd, pipo_and_dead_qmd, overall_qmd,
                                  live_sdi_ha, dead_sdi_ha, pipo_sdi_ha, non_pipo_sdi_ha, pipo_and_dead_sdi_ha, overall_sdi_ha,
                                  live_sdi_ac, dead_sdi_ac, pipo_sdi_ac, non_pipo_sdi_ac, pipo_and_dead_sdi_ac, overall_sdi_ac,
                                  live_mean_voronoi, dead_mean_voronoi, pipo_mean_voronoi, non_pipo_mean_voronoi, pipo_and_dead_mean_voronoi, overall_mean_voronoi,
                                  live_mean_nn1, dead_mean_nn1, pipo_mean_nn1, non_pipo_mean_nn1, pipo_and_dead_mean_nn1, overall_mean_nn1,
                                  live_mean_nn2, dead_mean_nn2, pipo_mean_nn2, non_pipo_mean_nn2, pipo_and_dead_mean_nn2, overall_mean_nn2,
                                  live_mean_nn3, dead_mean_nn3, pipo_mean_nn3, non_pipo_mean_nn3, pipo_and_dead_mean_nn3, overall_mean_nn3,
                                  live_sd_nn1, dead_sd_nn1, pipo_sd_nn1, non_pipo_sd_nn1, pipo_and_dead_sd_nn1, overall_sd_nn1,
                                  live_cov_nn1, dead_cov_nn1, pipo_cov_nn1, non_pipo_cov_nn1, pipo_and_dead_cov_nn1, overall_cov_nn1,
                                  live_sd_voronoi, dead_sd_voronoi, pipo_sd_voronoi, non_pipo_sd_voronoi, pipo_and_dead_sd_voronoi, overall_sd_voronoi,
                                  live_cov_voronoi, dead_cov_voronoi, pipo_cov_voronoi, non_pipo_cov_voronoi, pipo_and_dead_cov_voronoi, overall_cov_voronoi,
                                  current_cwd, current_cwd_zscore)
  
  names(results_raster) <- c("live_count", "dead_count", "pipo_count", "non_pipo_count", "pipo_and_dead_count", "total_count", 
                             "live_tpha", "dead_tpha", "pipo_tpha", "non_pipo_tpha", "pipo_and_dead_tpha", "overall_tpha",
                             "live_mean_height", "dead_mean_height", "pipo_mean_height", "non_pipo_mean_height", "pipo_and_dead_mean_height", "overall_mean_height",
                             "live_ba", "dead_ba", "pipo_ba", "non_pipo_ba", "pipo_and_dead_ba", "total_ba",
                             "live_bapha", "dead_bapha", "pipo_bapha", "non_pipo_bapha", "pipo_and_dead_bapha", "overall_bapha",
                             "live_mean_ba", "dead_mean_ba", "pipo_mean_ba", "non_pipo_mean_ba", "pipo_and_dead_mean_ba", "overall_mean_ba",
                             "live_mean_dbh", "dead_mean_dbh", "pipo_mean_dbh", "non_pipo_mean_dbh", "pipo_and_dead_mean_dbh", "overall_mean_dbh",
                             "live_qmd", "dead_qmd", "pipo_qmd", "non_pipo_qmd", "pipo_and_dead_qmd", "overall_qmd",
                             "live_sdi_ha", "dead_sdi_ha", "pipo_sdi_ha", "non_pipo_sdi_ha", "pipo_and_dead_sdi_ha", "overall_sdi_ha",
                             "live_sdi_ac", "dead_sdi_ac", "pipo_sdi_ac", "non_pipo_sdi_ac", "pipo_and_dead_sdi_ac", "overall_sdi_ac",
                             "live_mean_voronoi", "dead_mean_voronoi", "pipo_mean_voronoi", "non_pipo_mean_voronoi", "pipo_and_dead_mean_voronoi", "overall_mean_voronoi",
                             "live_mean_nn1", "dead_mean_nn1", "pipo_mean_nn1", "non_pipo_mean_nn1", "pipo_and_dead_mean_nn1", "overall_mean_nn1",
                             "live_mean_nn2", "dead_mean_nn2", "pipo_mean_nn2", "non_pipo_mean_nn2", "pipo_and_dead_mean_nn2", "overall_mean_nn2",
                             "live_mean_nn3", "dead_mean_nn3", "pipo_mean_nn3", "non_pipo_mean_nn3", "pipo_and_dead_mean_nn3", "overall_mean_nn3",
                             "live_sd_nn1", "dead_sd_nn1", "pipo_sd_nn1", "non_pipo_sd_nn1", "pipo_and_dead_sd_nn1", "overall_sd_nn1",
                             "live_cov_nn1", "dead_cov_nn1", "pipo_cov_nn1", "non_pipo_cov_nn1", "pipo_and_dead_cov_nn1", "overall_cov_nn1",
                             "live_sd_voronoi", "dead_sd_voronoi", "pipo_sd_voronoi", "non_pipo_sd_voronoi", "pipo_and_dead_sd_voronoi", "overall_sd_voronoi",
                             "live_cov_voronoi", "dead_cov_voronoi", "pipo_cov_voronoi", "non_pipo_cov_voronoi", "pipo_and_dead_cov_voronoi", "overall_cov_voronoi",
                             "local_cwd", "local_cwd_zscore")
  
  writeRaster(x = results_raster, filename = here::here(paste0("data/data_drone/L4/rasterized-trees/", current_site, "_rasterized-trees.tif")), overwrite = TRUE)
  
  results_df <- 
    results_raster %>% 
    as.data.frame(xy = TRUE) %>% 
    dplyr::mutate(site = current_site) %>% 
    tidyr::separate(col = site, into = c("forest", "elev", "rep"), remove = FALSE) %>% 
    dplyr::mutate(crs = 3310)
  
  results_list[[i]] <- results_df
}

final_results <- do.call("rbind", results_list)

write.csv(final_results, file = here::here("data/data_drone/L4/data-from-rasterized-classified-trees.csv"), row.names = FALSE)
