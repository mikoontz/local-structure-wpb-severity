# Get plot centers as discovered from finding orange X's on the orthophotos.

library(sf)
library(tidyverse)
library(purrr)
library(lidR)
library(viridis)
library(raster)
library(ForestTools)
library(gstat)

# Pseudocode
# Calculate the "validation metrics" for the ground data at the plots that are visible in the drone-derived orthophotos
#
# (Explanation: we define "validation metrics" to be things that we care about in the forest
# By picking 2+ "validation metrics" that are somewhat non-correlated and then matching them to the same metrics calculated
# using the treetop detection/crown segmentation algorithm, we can compare the quality of different algorithms.
# The greater the combined correlation of the "validation metrics" between algorithm-derived metrics and ground data-derived
# metrics, the better the algorithm.)
#
# Candidate "validation metrics" include:
# tree density within a plot
# tree size distribution within a plot
# proportion of live trees
# Voronoi polygon areas might not work due to edge effects, but perhaps distance to nearest tree?

# Establish a method for detecting tree tops and segmenting crowns
# Work through all of the remotely-visible plots at each site (144 in total) applying the ttops and crown segmentation approach
# 

# First get the formatted ground data; the object is called `d`
source("data/data_carpentry/format_ground-data.R")

ground_trees <- 
  list.files("data/data_output/site_data") %>% 
  map(.f = function(current_site) {
    current_plot_location <-
      st_read(paste0("data/data_output/site_data/", current_site, "/", current_site, "_plot-locations/", current_site, "_plot-locations.shp")) %>% 
      st_zm() %>% 
      dplyr::arrange(id) %>% 
      dplyr::mutate(plot = paste(current_site, id, sep = "_"))
    
    current_plot_ground_trees <- 
      d %>% 
      filter(site == current_site) %>% 
      left_join(current_plot_location, by = "plot") %>% 
      dplyr::mutate(delta_x = cospi((1 / 2) - (azm / 180)) * dist) %>% 
      dplyr::mutate(delta_y = sinpi((1 / 2) - (azm / 180)) * dist) %>%
      st_as_sf() %>%
      dplyr::mutate(x = st_coordinates(.)[, "X"],
                    y = st_coordinates(.)[, "Y"]) %>% 
      as.data.frame() %>% 
      dplyr::mutate(new_x = x + delta_x,
                    new_y = y + delta_y) %>% 
      filter(!is.na(x)) %>% 
      st_as_sf(coords = c("new_x", "new_y")) %>% 
      st_set_crs(st_crs(current_plot_location)) %>% 
      st_transform(4326)
    
    current_plot_ground_trees    
  }) %>%
  do.call("rbind", .)

# Use the plot radius to determine the area of the plot and thus the tree density
plot_radius <- sqrt((66*66) / pi) * 12 *2.54 / 100

# Get the validation metrics for the ground plot trees (that hadn't fallen by the time of the flights)
ground_tree_summary <-
  ground_trees %>% 
  dplyr::filter(is.na(year_fall)) %>% 
  dplyr::group_by(plot) %>% 
  dplyr::summarize(total_tree_count = n(),
                   live_tree_count = sum(live),
                   dead_tree_count = n() - sum(live),
                   total_density_tph = (n() * 10000) / (pi * plot_radius ^ 2),
                   live_density_tph = (sum(live) * 10000) / (pi * plot_radius ^ 2),
                   dead_density_tph = ((n() - sum(live)) * 10000) / (pi * plot_radius ^ 2),
                   live_proportion = sum(live) / n(),
                   height_lwr_25 = quantile(height, prob = 0.25),
                   height_mean = mean(height),
                   height_upr_25 = quantile(height, prob = 0.75))


# Now there is an R object in the environment called "sites_checklist" that has
# infomation about how far along all processing steps are.
source("data/data_carpentry/make_processing-checklist.R")

all_sites <- list.files("data/data_output/site_data/")
all_validation_plots <- ground_tree_summary$plot

for (i in seq_along(all_sites)) {
  current_site <- all_sites[i]
  current_dir <- paste0("data/data_output/site_data/", current_site, "/")
  
  current_plot_location <-
    st_read(paste0("data/data_output/site_data/", current_site, "/", current_site, "_plot-locations/", current_site, "_plot-locations.shp")) %>% 
    st_zm() %>% 
    dplyr::arrange(id) %>% 
    dplyr::mutate(plot = paste(current_site, id, sep = "_"))
  
  # The plots_at_current_site object represents the plots that are visible on the orthophoto
  plots_at_current_site <- 
    ground_tree_summary %>% 
    dplyr::mutate(site = substr(all_validation_plots, start = 1, stop = 9)) %>% 
    dplyr::filter(site == current_site) %>% 
    pull(plot)
  
  for (j in seq_along(plots_at_current_site)) {
    current_plot <- plots_at_current_site[j]
    
    current_ortho <- raster::brick(paste0(current_dir, current_site, "_plot-remote-data/", current_plot, "_ortho.tif"))
    current_dsm <- raster::raster(paste0(current_dir, current_site, "_plot-remote-data/", current_plot, "_dsm.tif"))
    current_dtm <- raster::raster(paste0(current_dir, current_site, "_plot-remote-data/", current_plot, "_dtm.tif"))
    current_las <- lidR::readLAS(paste0(current_dir, current_site, "_plot-remote-data/", current_plot, "_normalized-point-cloud.las"))

    current_plot_boundary <-
      current_plot_location %>% 
      dplyr::filter(plot == current_plot) %>% 
      st_buffer(plot_radius)

    # First make the Canopy Height Model using the DSM and the DTM
    # Using bilinear interpolation to downsample the 2m resolution DTM to have the
    # same resolution as the dsm (~5cm, but slightly different for each site)
    current_dtm_resamp <- raster::resample(current_dtm, current_dsm, method = "bilinear")
    
    # The Canopy Height Model (chm) is the dsm (vegetation + ground) minus the dtm (ground)
    # to give just the height of the vegetation.
    current_chm_rough <- current_dsm - current_dtm_resamp
    
    # Smooth out the chm and set any negative values to 0 (meaning "ground") following
    # advice from Zagalikis, Cameron, and Miller (2004) and references therein
    # More recently, a 3x3 pixel smoothing filter was specifically suggested as ideal
    # for sUAS derived chm by Mohan et al. (2017)
    # w <- focalWeight(x = chm, d = 0.5, type = "circle")
    # w
    current_chm <- raster::focal(current_chm_rough, w = matrix(1/49, nrow = 7, ncol = 7))
    current_chm[raster::getValues(current_chm) < 0] <- 0
    
    
    # First step for many segmentation algorithms is to detect tree tops
    # which can be done in a number of ways
    
    # We only measured trees that were greater than 2.5 inches DBH and focused on ponderosa pine
    # which translated to a height of approximately 6 meters (Wonn and O'Hara, 2001) so we ignored
    # trees less than this height
    # Empirically, over 99% of trees were greater than 3 meters in height... (1 - ecdf(ground_trees$height)(3) = 0.99001)
    # ... and over 80% of trees measured in the ground plots were greater than 6 meters in height (1 - ecdf(ground_trees$height)(6) = 0.8085248)
    # ... and over 90% of trees measured in the ground plots were greater than 4.5 meters in height (1 - ecdf(ground_trees$height)(4.5) = 0.9020979)
    # Uses the "variable window filter" algorithm by Popescu and Wynne (2004)
    
    min_height <- 3
    
    # Methods 1 through 3 use a "variable window filter" to search for local maxima within window sizes
    # that are functions of how high a particular pixel is
    
    # Method 1
    # Using the window function defined in the help file
    dynamicWindow_default <- function(x){x * 0.06 + 0.5}
    
    # Method 2
    # Equation and coefficients taken from Popescu and Wynne (2004)
    # Divide the Popescu and Wynne (2004) equations by 2 to convert to *radius* of search window
    # which is what the vwf() function requires
    
    # Equation and coefficients taken from Popescu and Wynne (2004)'s "Pines model"
    dynamicWindow_pines <- function(x) {
      window_radius <- 0.5 * (3.75105 - 0.17919*x + 0.01241 * x^2)
      return(window_radius)
    }
    
    # Method 3
    # Equation and coefficients taken from Popescu and Wynne (2004)
    # Divide the Popescu and Wynne (2004) equations by 2 to convert to *radius* of search window
    # which is what the vwf() function requires
    
    # Equation and coefficients taken from Popescu and Wynne (2004)'s conifer + deciduous model: "Combined model"
    dynamicWindow_combined <- function(x) {
      window_radius <- 0.5 * (2.51503 + 0.00901 * x^2)
      return(window_radius)
    }
    
    ttops_vwf_default <- 
      ForestTools::vwf(CHM = current_chm, winFun = dynamicWindow_default, minHeight = min_height, maxWinDiameter = NULL) %>% 
      st_as_sf() %>% 
      st_intersection(y = current_plot_boundary)
    
    ttops_vwf_pines <- 
      ForestTools::vwf(CHM = current_chm, winFun = dynamicWindow_pines, minHeight = min_height, maxWinDiameter = NULL) %>% 
      st_as_sf() %>% 
      st_intersection(y = current_plot_boundary)
    
    ttops_vwf_combined <- 
      ForestTools::vwf(CHM = current_chm, winFun = dynamicWindow_combined, minHeight = min_height, maxWinDiameter = NULL) %>% 
      st_as_sf() %>% 
      st_intersection(y = current_plot_boundary)
    

    # Method 4
    # Using a fixed square window for detecting local maxima
    # 1.5 meters per side in pixels
    window_size <- 1.5
    ws_in_pixels <- round(window_size * (1 / res(current_chm)[1]))
    ws_in_pixels <- ifelse(ws_in_pixels %% 2 == 0, yes = ws_in_pixels + 1, no = ws_in_pixels)

    # On a raster-like object, the window size is in pixels
    # The returned value when a raster-like object is used is also a raster-like object, so we transform it into a
    # sf POINT object
    ttops_localMaxima_chm_1.5 <- 
      current_chm %>% 
      lidR::tree_detection(ws = ws_in_pixels) %>% 
      xyFromCell(which(.[] > 0)) %>% 
      as.data.frame() %>% 
      st_as_sf(coords = c("x", "y"),
               crs = proj4string(current_chm)) %>% 
      mutate(id = 1:nrow(.)) %>% 
      mutate(height = extract(current_chm, as(., "Spatial"))) %>% 
      st_intersection(y = current_plot_boundary)
    
    # Method 5
    # Using a fixed window for detecting local maxima
    # On a raw point cloud object, the window size is in meters
    # 1.5 las units per side = 1.5 meters
    ttops_localMaxima_las_1.5 <- 
      current_las %>% 
      lidR::tree_detection(ws = window_size) %>% 
      as.data.frame() %>% 
      st_as_sf(coords = c("X", "Y"),
               crs = proj4string(current_chm)) %>% 
      mutate(treeID = 1:nrow(.)) %>% 
      rename(height = Z) %>% 
      st_intersection(y = current_plot_boundary)
    
    # Method 6
    # Using a fixed window for detecting local maxima
    # 2 meters on a side in pixels
    window_size <- 2
    ws_in_pixels <- round(window_size * (1 / res(current_chm)[1]))
    ws_in_pixels <- ifelse(ws_in_pixels %% 2 == 0, yes = ws_in_pixels + 1, no = ws_in_pixels)
    
    # On a raster-like object, the window size is in pixels
    # The returned value when a raster-like object is used is also a raster-like object, so we transform it into a
    # sf POINT object
    ttops_localMaxima_chm_2 <- 
      current_chm %>% 
      lidR::tree_detection(ws = ws_in_pixels) %>% 
      xyFromCell(which(.[] > 0)) %>% 
      as.data.frame() %>% 
      st_as_sf(coords = c("x", "y"),
               crs = proj4string(current_chm)) %>% 
      mutate(id = 1:nrow(.)) %>% 
      mutate(height = extract(current_chm, as(., "Spatial"))) %>% 
      st_intersection(y = current_plot_boundary)
    
    # Method 7
    # Using a fixed window for detecting local maxima
    # On a raw point cloud object, the window size is in meters
    # 2 las units on a side = 2 meters
    ttops_localMaxima_las_2 <- 
      current_las %>% 
      lidR::tree_detection(ws = window_size) %>% 
      as.data.frame() %>% 
      st_as_sf(coords = c("X", "Y"),
               crs = proj4string(current_chm)) %>% 
      mutate(treeID = 1:nrow(.)) %>% 
      rename(height = Z) %>% 
      st_intersection(y = current_plot_boundary)
    
    # Method 8
    # Using a fixed window for detecting local maxima
    # 2.5 meters on a side in pixels
    window_size <- 2.5
    ws_in_pixels <- round(window_size * (1 / res(current_chm)[1]))
    ws_in_pixels <- ifelse(ws_in_pixels %% 2 == 0, yes = ws_in_pixels + 1, no = ws_in_pixels)
    
    # On a raster-like object, the window size is in pixels
    # The returned value when a raster-like object is used is also a raster-like object, so we transform it into a
    # sf POINT object
    ttops_localMaxima_chm_2.5 <- 
      current_chm %>% 
      lidR::tree_detection(ws = ws_in_pixels) %>% 
      xyFromCell(which(.[] > 0)) %>% 
      as.data.frame() %>% 
      st_as_sf(coords = c("x", "y"),
               crs = proj4string(current_chm)) %>% 
      mutate(id = 1:nrow(.)) %>% 
      mutate(height = extract(current_chm, as(., "Spatial"))) %>% 
      st_intersection(y = current_plot_boundary)
    
    # Method 9
    # Using a fixed window for detecting local maxima
    # On a raw point cloud object, the window size is in meters
    # 2.5 las units = 2.5 meters
    ttops_localMaxima_las_2.5 <- 
      current_las %>% 
      lidR::tree_detection(ws = window_size) %>% 
      as.data.frame() %>% 
      st_as_sf(coords = c("X", "Y"),
               crs = proj4string(current_chm)) %>% 
      mutate(treeID = 1:nrow(.)) %>% 
      rename(height = Z) %>% 
      st_intersection(y = current_plot_boundary)
    
    # Method 10
    # Using the approach developed by Li et al. (2012)
    # Using values from Shin et al. (2018) Remote Sensing.
    # "The segmentation method relies on four user-defined parameters: minimum height of a tree, maximum crown radius, and two numeric distances, which define horizontal distance (in meters) thresholds between all points above 15 m in height, and below 15 m in height. These thresholds values are hereafter referred to as distance thresholds (DT). The DT value is site specific [42]. A low DT value generally results in over-segmentation with many additional trees identified in the point cloud, whereas a high DT value causes under-segmentation, where many tree canopies are merged together into single large canopies. In this study, 16 different iterations of varying DT values were tested. We used a minimum tree height of 2 m, and a maximum canopy diameter of 7 m, based on the ranges observed in our field data."
    # "The parameters used in a given tree segmentation algorithm must be “tuned” to match the specific site and user’s needs. We used a point-based algorithm [42] to segment individual trees from the point cloud. The main parameter that affected the segmentation was the DT parameter—a distance threshold between points that determined whether a point was or was not part of a particular tree. Within the Li et al. [42] segmentation algorithm, there are two different DT values, both of which were set as equal in our study. In future studies, these values can be set differently to potentially achieve better segmentation results."
    # "The optimized iteration contains two DT values: 1.4 m for areas with more than 50% canopy cover, and 1.7 m for areas of 50% or less canopy cover."
    lastrees_li2(current_las, dt1 = 1.4, dt2 = 1.4, R = 2, Zu = 15, hmin = min_height, speed_up = 10)

    crowns_li2 <- 
      current_las %>% 
      as.spatial() %>% 
      st_as_sf(coords = c("X", "Y"),
               crs = proj4string(current_chm)) %>% 
      mutate(treeID = 1:nrow(.)) %>% 
      rename(height = Z) %>% 
      st_intersection(y = current_plot_boundary)
    
    ttops_li2 <-
      current_las@data %>% 
      dplyr::group_by(treeID) %>% 
      dplyr::filter(Z == max(Z)) %>% 
      ungroup() %>% 
      st_as_sf(coords = c("X", "Y"),
               crs = proj4string(current_chm)) %>% 
      rename(height = Z) %>% 
      st_intersection(y = current_plot_boundary)
    
    
    # Method 11
    # tolerance: The minimum height of the object in the units of image intensity between its highest point (seed) and the point where it contacts another object (checked for every contact pixel). If the height is smaller than the tolerance, the object will be combined with one of its neighbors, which is the highest. Tolerance should be chosen according to the range of x. Default value is 1, which is a reasonable value if x comes from distmap.
    # ext: Radius of the neighborhood in pixels for the detection of neighboring objects. Higher value smoothes out small objects.

    crowns_watershed <- 
      lidR::lastrees_watershed(chm = current_chm, th_tree = 2, tol = 3, ext = round((1 / res(current_chm)[1])), extra = TRUE) %>% 
      APfun::APpolygonize() %>% 
      st_as_sf() %>% 
      rename(treeID = DN) %>% 
      group_by(treeID) %>%
      summarize() %>% 
      st_set_crs(proj4string(current_chm))
    
    ttops_watershed <- 
      current_chm %>% 
      raster::extract(y = crowns_watershed, cellnumbers = TRUE, df = TRUE) %>% 
      rename(treeID = ID, height = layer) %>% 
      group_by(treeID) %>% 
      filter(height == max(height)) %>% 
      ungroup() %>% 
      bind_cols(xyFromCell(current_chm, .$cell) %>% 
                  as.data.frame()) %>% 
      st_as_sf(coords = c("x", "y"),
               crs = proj4string(current_chm)) %>% 
      st_intersection(y = current_plot_boundary)

    
    # All the different ttops
    
    ttops_vwf_default
    ttops_vwf_pines
    ttops_vwf_combined
    ttops_localMaxima_chm_1.5
    ttops_localMaxima_las_1.5
    ttops_localMaxima_chm_2
    ttops_localMaxima_las_2
    ttops_localMaxima_chm_2.5
    ttops_localMaxima_las_2.5
    ttops_li2
    ttops_watershed
    
    # Crown Segmentation as a 2nd step
    # Method 1
    # EBImage watershed
    # crowns_watershed
    
    # Method 2
    # Li algorithm
    # crowns_li2
    
    # # Methods 3 through 5 rely on tree top detection ahead of time
    # # Method 3
    # # Dalponte
    # lastrees_dalponte(current_las, current_chm, ttops, th_tree = 2, th_seed = 0.45,
    #                   th_cr = 0.55, max_cr = 10, extra = TRUE)
    # 
    # # Method 4
    # # Silva
    # lastrees_silva(current_las, current_chm, ttops, max_cr_factor = 0.6, exclusion = 0.3,
    #                extra = TRUE)
    # 
    # # Method 5
    # # Marker controlled watershed segmentation
    # crowns <- 
    #   ForestTools::mcws(treetops = as(ttops, "Spatial"), 
    #                     CHM = current_chm, 
    #                     minHeight = 2, 
    #                     format = "polygons", 
    #                     OSGeoPath = "C:\\OSGeo4W64") %>% 
    #   st_as_sf() %>% 
    #   mutate(treeID = 1:nrow(.))
    }
}




# lastrees_li2(current_plot_las, dt1 = 1, dt2 = 1.5, R = 1, Zu = 15, hmin = 5, speed_up = 25)
# 
# r <- 
#   lastrees_dalponte(las = current_plot_las, 
#                     chm = current_plot_chm, 
#                     treetops = as.data.frame(ttops3)[, -5], 
#                     th_tree = 4, 
#                     th_seed = .9, 
#                     th_cr = 0.95, 
#                     max_cr = 1 / res(current_plot_chm)[1] * 10,
#                     extra = TRUE)


