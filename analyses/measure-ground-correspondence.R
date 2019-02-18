library(foreach)
library(doParallel)
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

# begin main program ------------------------------------------------------


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
# live tree density
# dead tree density
# tree size distribution within a plot (25th and 75th percentile height; mean height)
# proportion of live trees
# Voronoi polygon areas might not work due to edge effects, but let's try mean distance to 3 nearest trees

# Establish a method for detecting tree tops and segmenting crowns
# Work through all of the remotely-visible plots at each site (144 in total) applying the ttops and crown segmentation approach

# These sites were processed with their X3 and RedEdge imagery combined so some of their
# output products will be in a slightly different place in the project directory
merged_sites <- c("eldo_3k_2",
                  "eldo_3k_3",
                  "eldo_4k_2")

unusable_sites <- c("eldo_4k_3", # too many blocks
                    "stan_4k_3", # too many blocks
                    "stan_5k_3", # too many blocks
                    "sequ_4k_2") # middle section flown on a separate day and the stitch looks terrible

# First get the formatted ground data; the object is called `d`
source("data/data_carpentry/format_ground-data.R")
source("data/data_carpentry/make_processing-checklist.R")

# This is where I can put in sites that need their processing redone. An empty 
# string means that no already-processed site output will be overwritten
# (but sites that have yet to be processed will still have their processing done)

sites_to_overwrite <- ""
sites_checklist$overwrite <- FALSE

sites_checklist[sites_checklist$site %in% sites_to_overwrite, "overwrite"] <- TRUE

sites_to_process <- 
  sites_checklist %>% 
  dplyr::filter(!(site %in% unusable_sites)) %>%
  dplyr::filter(overwrite | plot_remote_data_check) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()

ground_trees <- 
  sites_to_process %>% 
  map(.f = function(current_site) {
    
    if (current_site %in% merged_sites) {
      current_site_plot_locations <- 
        sf::st_read(paste0("data/data_output/site_data/", current_site, "/", current_site, "_plot-locations/", current_site, "_plot-locations.shp")) %>% 
        mutate(plot = paste(current_site, id, sep = "_")) %>%
        dplyr::arrange(id) %>% 
        st_zm()
    } else {
      current_site_plot_locations <- 
        sf::st_read(paste0("data/data_output/site_data/", current_site, "/", current_site, "_plot-locations_re/", current_site, "_plot-locations_re.shp")) %>% 
        mutate(plot = paste(current_site, id, sep = "_")) %>%
        dplyr::arrange(id) %>% 
        st_zm()
    }
    
    current_plot_ground_trees <- 
      d %>% 
      filter(site == current_site) %>% 
      left_join(current_site_plot_locations, by = "plot") %>% 
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
      st_set_crs(st_crs(current_site_plot_locations))
    
    nn <- 
      st_nn(current_plot_ground_trees, current_plot_ground_trees, k = min(c(4, nrow(current_plot_ground_trees))), returnDist = TRUE, sparse = FALSE, progress = FALSE)$dist %>% 
      as.data.frame()
    
    nn <- nn %>% dplyr::select(-1)
    
    if (ncol(nn) == 1) {nn <- nn %>% rename(nn1 = V2) %>% mutate(nn2 = NA, nn3 = NA)} else
      if (ncol(nn) == 2) {nn <- nn %>% rename(nn1 = V2, nn2 = V3) %>% mutate(nn3 = NA)} else
        if (ncol(nn) == 3) {nn <- nn %>% rename(nn1 = V2, nn2 = V3, nn3 = V4) }
    
    current_plot_ground_trees <- 
      current_plot_ground_trees %>% 
      bind_cols(nn)  %>% 
      st_transform(4326)
  }) %>%
  do.call("rbind", .)

# Use the plot radius to determine the area of the plot and thus the tree density
plot_radius <- sqrt((66*66) / pi) * 12 *2.54 / 100

# Get the validation metrics for the ground plot trees (that hadn't fallen by the time of the flights)
ground_tree_summary <-
  ground_trees %>% 
  dplyr::filter(is.na(year_fall)) %>% 
  dplyr::group_by(plot) %>% 
  dplyr::summarize(ttop_method = "ground",
                   total_tree_count = n(),
                   live_tree_count = sum(live),
                   dead_tree_count = n() - sum(live),
                   total_density_tph = (n() * 10000) / (pi * plot_radius ^ 2),
                   live_density_tph = (sum(live) * 10000) / (pi * plot_radius ^ 2),
                   dead_density_tph = ((n() - sum(live)) * 10000) / (pi * plot_radius ^ 2),
                   live_proportion = sum(live) / n(),
                   height_lwr_25 = quantile(height, prob = 0.25),
                   height_mean = mean(height),
                   height_upr_25 = quantile(height, prob = 0.75),
                   nn_1_mean = mean(nn1),
                   nn_2_mean = mean(nn2),
                   nn_3_mean = mean(nn3),
                   tree_count_above_15m = sum(height >= 15),
                   tree_count_below_15m = sum(height < 15)) %>% 
  as.data.frame() %>% 
  dplyr::select(-geometry)

ground_trees %>% 
  group_by(live) %>% 
  summarize(min = min(height),
            `10` = quantile(height, probs = 0.10),
            `20` = quantile(height, probs = 0.20),
            `30` = quantile(height, probs = 0.30),
            `40` = quantile(height, probs = 0.40),
            mean = mean(height),
            `50` = quantile(height, probs = 0.50),
            `60` = quantile(height, probs = 0.60),
            `70` = quantile(height, probs = 0.70),
            `80` = quantile(height, probs = 0.80),
            `90` = quantile(height, probs = 0.90),
            max = max(height),
            ecdf_05 = ecdf(height)(5),
            ecdf_10 = ecdf(height)(10),
            ecdf_15 = ecdf(height)(15),
            ecdf_20 = ecdf(height)(20),
            ecdf_10_dbh = ecdf(dbh)(10),
            ecdf_15_dbh = ecdf(dbh)(15),
            ecdf_20_dbh = ecdf(dbh)(20),
            ecdf_25_dbh = ecdf(dbh)(25))

ggplot(ground_trees, aes(x = dbh, y = height, col = as.factor(live))) +
  geom_point() +
  geom_smooth(method = "lm")

# Now we have a summary of what all of the ground plots looks like with respsect
# to a few of these plot-level characteristics.
# Now we need to employ a bunch of segmentation algorithms to find out which
# one matches best to the ground data.

all_validation_plots <- ground_tree_summary$plot

start <- Sys.time()

# I revisit this StackOverflow question and answer, like, everytime I want
# to do anything in parallel to remind myself how
# https://stackoverflow.com/questions/38318139/run-a-for-loop-in-parallel-in-r

cores <- detectCores() 
cl <- makeCluster(min(c(cores - 2), 8)) #not to overload your computer
registerDoParallel(cl)

ttops_summary <- foreach(i = seq_along(all_validation_plots), .combine = rbind) %do% {
  
  library(sf)
  library(tidyverse)
  library(purrr)
  library(lidR)
  library(viridis)
  library(raster)
  library(ForestTools)
  library(gstat)
  library(nngeo)
  library(lidRplugins)
  
  source("analyses/vwf2.R")
  
  
  # treetop detection algorithm helper functions ----------------------------
  
  
  st_vwf <- function(CHM, plot_boundary, winFun, minHeight, maxWinDiameter) {
    start <- Sys.time()
    
    ttops_sf <-
      CHM %>% 
      vwf2(winFun = winFun, minHeight = minHeight, maxWinDiameter = maxWinDiameter) %>%
      st_as_sf() %>%
      st_set_agr("constant") %>%
      st_intersection(y = plot_boundary)
    
    end <- Sys.time()
    
    ttops_time <- end - start
    
    return(list(ttops_sf = ttops_sf, ttops_time = ttops_time))
  }
  
  
  st_lmf <- function(obj, plot_boundary, ws) {
    start <- Sys.time()
    ttops_sf <-
      obj %>%
      lidR::tree_detection(algorithm = lmf(ws = ws)) %>%
      st_as_sf() %>% 
      rename(height = Z) %>%
      st_set_agr("constant") %>%
      st_intersection(y = plot_boundary)
    
    end <- Sys.time()
    
    ttops_time <- end - start
    
    return(list(ttops_sf = ttops_sf, ttops_time = ttops_time))
  }
  
  
  st_li2012 <- function(las, plot_boundary, dt1, dt2, R, Zu, hmin, speed_up) {
    start <- Sys.time()
    ttops_las <- 
      las %>% 
      lidR::lastrees(algorithm = li2012(dt1 = dt1, dt2 = dt2, R = R, Zu = Zu, hmin = min_height, speed_up = speed_up))
    
    ttops_sf <-
      ttops_las %>% 
      slot("data") %>% 
      dplyr::group_by(treeID) %>%
      dplyr::filter(Z == max(Z)) %>%
      dplyr::ungroup() %>%
      st_as_sf(coords = c("X", "Y"),
               crs = proj4string(las)) %>%
      rename(height = Z) %>%
      st_set_agr("constant") %>%
      st_intersection(y = plot_boundary)
    
    end <- Sys.time()
    
    ttops_time <- end - start
    
    return(list(ttops_las = ttops_las, ttops_sf = ttops_sf, ttops_time = ttops_time))
  }
  
  st_watershed <- function(las, plot_boundary, chm, th_tree, tol, ext) {
    start <- Sys.time()
    
    ttops_las <-
      las %>% 
      lidR::lastrees(algorithm = watershed(chm = current_chm, th_tree = th_tree, tol = tol, ext = ext))
    
    ttops_sf <-  
      ttops_las %>% 
      slot("data") %>% 
      dplyr::group_by(treeID) %>%
      dplyr::filter(Z == max(Z)) %>%
      dplyr::ungroup() %>%
      st_as_sf(coords = c("X", "Y"),
               crs = proj4string(las)) %>%
      rename(height = Z) %>%
      st_set_agr("constant") %>%
      st_intersection(y = plot_boundary)
    
    end <- Sys.time()
    
    ttops_time <- end - start
    
    return(list(ttops_las = ttops_las, ttops_sf = ttops_sf, ttops_time = ttops_time))
  }
  
  st_ptree <- function(las, plot_boundary, k, hmin, nmax) {
    start <- Sys.time()
    
    ttops_las <-
      las %>%
      lidR::lastrees(algorithm = ptrees(k = k, hmin = hmin, nmax = nmax))
    
    ttops_sf <-
      ttops_las %>%
      slot("data") %>% 
      dplyr::group_by(treeID) %>%
      dplyr::filter(Z == max(Z)) %>%
      dplyr::ungroup() %>%
      st_as_sf(coords = c("X", "Y"),
               crs = proj4string(las)) %>%
      rename(height = Z) %>%
      st_set_agr("constant") %>%
      st_intersection(y = plot_boundary)
    
    end <- Sys.time()
    ttops_time <- end - start
    
    return(list(ttops_las = ttops_las, ttops_sf = ttops_sf, ttops_time = ttops_time))
    
  }
  
  st_multichm <- function(las, plot_boundary, res, layer_thickness, dist_2d, dist_3d, use_max, ws) {
    
    start <- Sys.time()
    
    ttops_sf <-
      las %>% 
      lidR::tree_detection(algorithm = multichm(res = res, layer_thickness = layer_thickness, dist_2d = dist_2d, dist_3d = dist_3d, use_max = use_max, ws = ws)) %>% 
      st_as_sf() %>% 
      rename(height = Z) %>%
      st_set_agr("constant") %>%
      st_intersection(y = plot_boundary)
    
    end <- Sys.time()
    ttops_time <- end - start
    
    return(list(ttops_sf = ttops_sf, ttops_time = ttops_time))
  }
  
  
  st_lmfx <- function(las, plot_boundary, hmin, dist_2d, ws) {
    
    start <- Sys.time()
    
    ttops_sf <- 
      las %>% 
      lidR::tree_detection(algorithm = lmfx(hmin = hmin, dist_2d = dist_2d, ws = ws)) %>% 
      st_as_sf() %>% 
      rename(height = Z) %>%
      st_set_agr("constant") %>%
      st_intersection(y = plot_boundary)
    
    end <- Sys.time()
    
    ttops_time <- end - start
    
    return(list(ttops_sf = ttops_sf, ttops_time = ttops_time))
    
  }
  
  # get data for particular plot --------------------------------------------
  
  
  current_plot <- all_validation_plots[i]
  current_site <- substr(current_plot, start = 1, stop = 9)
  current_dir <- paste0("data/data_output/site_data/", current_site, "/")
  
  if (current_site %in% merged_sites) {
    current_site_plot_locations <- 
      sf::st_read(paste0("data/data_output/site_data/", current_site, "/", current_site, "_plot-locations/", current_site, "_plot-locations.shp"), quiet = TRUE) %>% 
      mutate(plot = paste(current_site, id, sep = "_")) %>%
      dplyr::arrange(id) %>% 
      st_zm()
  } else {
    current_site_plot_locations <- 
      sf::st_read(paste0("data/data_output/site_data/", current_site, "/", current_site, "_plot-locations_re/", current_site, "_plot-locations_re.shp"), quiet = TRUE) %>% 
      mutate(plot = paste(current_site, id, sep = "_")) %>%
      dplyr::arrange(id) %>% 
      st_zm()
  }
  
  # The Canopy Height Model (chm) is the dsm (vegetation + ground) minus the dtm (ground)
  # to give just the height of the vegetation.
  current_chm_rough <- raster::raster(paste0(current_dir, current_site, "_plot-remote-data/", current_plot, "_chm.tif"))
  
  # The point cloud is directly used for some segmentation algorithms, so we import that too
  current_las <- suppressWarnings(suppressMessages(lidR::readLAS(paste0(current_dir, current_site, "_plot-remote-data/", current_plot, "_classified-point-cloud.las"))))
  current_las_normalized <- lidR::lasnormalize(las = current_las, algorithm = tin())
  
  plot_radius <- sqrt((66*66) / pi) * 12 *2.54 / 100
  
  current_plot_boundary <-
    current_site_plot_locations %>% 
    dplyr::filter(plot == current_plot) %>% 
    st_buffer(plot_radius) %>% 
    st_set_agr("constant")
  
  
  # Smooth out the chm and set any negative values to 0 (meaning "ground") following
  # advice from Zagalikis, Cameron, and Miller (2004) and references therein
  # More recently, a 3x3 pixel smoothing filter was specifically suggested as ideal
  # for sUAS derived chm by Mohan et al. (2017)
  # w <- focalWeight(x = chm, d = 0.5, type = "circle")
  current_chm <- raster::focal(current_chm_rough, w = matrix(1/9, nrow = 3, ncol = 3))
  # current_chm <- current_chm_rough
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
  
  
  # vwf treetop detection ---------------------------------------------------
  
  
  # These algorithms use a "variable window filter" to search for local maxima within window sizes
  # that are functions of how high a particular pixel is
  
  # Using the window function defined in the help file
  dynamicWindow_default <- function(x){x * 0.06 + 0.5}
  
  ttops_vwf_default <-
    current_chm %>% 
    st_vwf(plot_boundary = current_plot_boundary, winFun = dynamicWindow_default, minHeight = min_height, maxWinDiameter = NULL)
  
  # Equation and coefficients taken from Popescu and Wynne (2004)
  # Divide the Popescu and Wynne (2004) equations by 2 to convert to *radius* of search window
  # which is what the vwf() function requires
  
  # Equation and coefficients taken from Popescu and Wynne (2004)'s "Pines model"
  dynamicWindow_pines <- function(x) {
    window_radius <- 0.5 * (3.75105 - 0.17919*x + 0.01241 * x^2)
    return(window_radius)
  }
  
  ttops_vwf_pines <-
    current_chm %>% 
    st_vwf(plot_boundary = current_plot_boundary, winFun = dynamicWindow_pines, minHeight = min_height, maxWinDiameter = NULL)
  
  # Equation and coefficients taken from Popescu and Wynne (2004)
  # Divide the Popescu and Wynne (2004) equations by 2 to convert to *radius* of search window
  # which is what the vwf() function requires
  
  # Equation and coefficients taken from Popescu and Wynne (2004)'s conifer + deciduous model: "Combined model"
  dynamicWindow_combined <- function(x) {
    window_radius <- 0.5 * (2.51503 + 0.00901 * x^2)
    return(window_radius)
  }
  
  ttops_vwf_combined <-
    current_chm %>% 
    st_vwf(plot_boundary = current_plot_boundary, winFun = dynamicWindow_combined, minHeight = min_height, maxWinDiameter = NULL)
  
  
  
  # lmf treetop detection ---------------------------------------------------
  
  
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
    st_lmf(plot_boundary = current_plot_boundary, ws = ws_in_pixels)
  
  # Method 5
  # Using a fixed window for detecting local maxima
  # On a raw point cloud object, the window size is in meters
  # 1.5 las units per side = 1.5 meters
  ttops_localMaxima_las_1.5 <-
    current_las_normalized %>% 
    st_lmf(plot_boundary = current_plot_boundary, ws = window_size)
  
  
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
    st_lmf(plot_boundary = current_plot_boundary, ws = ws_in_pixels)
  
  # Using a fixed window for detecting local maxima
  # On a raw point cloud object, the window size is in meters
  # 2 las units on a side = 2 meters
  ttops_localMaxima_las_2 <-
    current_las_normalized %>%
    st_lmf(plot_boundary = current_plot_boundary, ws = window_size)
  
  
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
    st_lmf(plot_boundary = current_plot_boundary, ws = ws_in_pixels)
  
  # Using a fixed window for detecting local maxima
  # On a raw point cloud object, the window size is in meters
  # 2.5 las units = 2.5 meters
  
  ttops_localMaxima_las_2.5 <-
    current_las_normalized %>%
    st_lmf(plot_boundary = current_plot_boundary, ws = window_size)
  
  
  # li2012 treetop detection ------------------------------------------------
  
  
  # Using the approach developed by Li et al. (2012)
  # Using values from Shin et al. (2018) Remote Sensing.
  # "The segmentation method relies on four user-defined parameters: minimum height of a tree, 
  # maximum crown radius, and two numeric distances, which define horizontal distance (in meters) 
  # thresholds between all points above 15 m in height, and below 15 m in height. These 
  # thresholds values are hereafter referred to as distance thresholds (DT). The DT value is site 
  # specific [42]. A low DT value generally results in over-segmentation with many additional 
  # trees identified in the point cloud, whereas a high DT value causes under-segmentation, where 
  # many tree canopies are merged together into single large canopies. In this study, 16 different 
  # iterations of varying DT values were tested. We used a minimum tree height of 2 m, and a 
  # maximum canopy diameter of 7 m, based on the ranges observed in our field data."
  
  # "The parameters used in a given tree segmentation algorithm must be “tuned” to match the 
  # specific site and user’s needs. We used a point-based algorithm [42] to segment individual 
  # trees from the point cloud. The main parameter that affected the segmentation was the DT 
  # parameter—a distance threshold between points that determined whether a point was or was not 
  # part of a particular tree. Within the Li et al. [42] segmentation algorithm, there are two 
  # different DT values, both of which were set as equal in our study. In future studies, these 
  # values can be set differently to potentially achieve better segmentation results."
  # "The optimized iteration contains two DT values: 1.4 m for areas with more than 50% canopy 
  # cover, and 1.7 m for areas of 50% or less canopy cover."
  
  # Could systematically iterate through many parameter combinations...
  # dt1_vals <- c(1.0, 1.4)
  # dt2_vals <- c(1.0, 1.4)
  # R_vals <- 0
  # Zu_vals <- c(15, 25)
  # min_height <- 3
  # speed_up <- 10
  # 
  # li2012_params <- as_data_frame(expand.grid(dt1 = dt1_vals, dt2 = dt2_vals, R = R_vals, Zu = Zu_vals, hmin = min_height, speed_up = speed_up))
  # 
  # li2012_summary <-
  #   li2012_params %>% 
  #   pmap(.f = st_li2012, las = current_las_normalized, plot_boundary = current_plot_boundary)
  
  # Default parameters
  ttops_li2012_dt1_1.4_dt2_1.4_zu_15_R_0_speedUp_10 <- 
    current_las_normalized %>% 
    st_li2012(plot_boundary = current_plot_boundary, dt1 = 1.4, dt2 = 1.4, R = 0, Zu = 15, hmin = min_height, speed_up = 10)
  
  # Using approach developed by Li et al. (2012)
  # Using values from Jakubowski et al. (2013), a study done in mixed-conifer forests in the Sierra Nevada (Tahoe National Forest, mostly)
  ttops_li2012_dt1_1.5_dt2_2.0_zu_15_R_0_speedUp_10 <- 
    current_las_normalized %>% 
    st_li2012(plot_boundary = current_plot_boundary, dt1 = 1.5, dt2 = 2, R = 0, Zu = 15, hmin = min_height, speed_up = 10)
  
  ttops_li2012_dt1_1.0_dt2_2.0_zu_15_R_0_speedUp_10 <-
    current_las_normalized %>% 
    st_li2012(plot_boundary = current_plot_boundary, dt1 = 1.0, dt2 = 2, R = 0, Zu = 15, hmin = min_height, speed_up = 10)
  
  ttops_li2012_dt1_1.0_dt2_2.0_zu_20_R_0_speedUp_10 <-
    current_las_normalized %>% 
    st_li2012(plot_boundary = current_plot_boundary, dt1 = 1.0, dt2 = 2, R = 0, Zu = 20, hmin = min_height, speed_up = 10)
  
  ttops_li2012_dt1_1.0_dt2_1.5_zu_25_R_0_speedUp_10 <-
    current_las_normalized %>% 
    st_li2012(plot_boundary = current_plot_boundary, dt1 = 1.0, dt2 = 1.5, R = 0, Zu = 25, hmin = min_height, speed_up = 10)
  
  ttops_li2012_dt1_1.0_dt2_1.5_zu_20_R_0_speedUp_10 <-
    current_las_normalized %>% 
    st_li2012(plot_boundary = current_plot_boundary, dt1 = 1.0, dt2 = 1.5, R = 0, Zu = 20, hmin = min_height, speed_up = 10) 
  
  ttops_li2012_dt1_1.0_dt2_1.4_zu_20_R_0_speedUp_10 <-
    current_las_normalized %>% 
    st_li2012(plot_boundary = current_plot_boundary, dt1 = 1.0, dt2 = 1.4, R = 0, Zu = 20, hmin = min_height, speed_up = 10)
  
  ttops_li2012_dt1_1.0_dt2_1.4_zu_25_R_0_speedUp_10 <-
    current_las_normalized %>% 
    st_li2012(plot_boundary = current_plot_boundary, dt1 = 1.0, dt2 = 1.4, R = 0, Zu = 25, hmin = min_height, speed_up = 10)
  
  ttops_li2012_dt1_1.0_dt2_1.3_zu_20_R_0_speedUp_10 <-
    current_las_normalized %>% 
    st_li2012(plot_boundary = current_plot_boundary, dt1 = 1.0, dt2 = 1.3, R = 0, Zu = 20, hmin = min_height, speed_up = 10)
  
  ttops_li2012_dt1_1.0_dt2_1.0_zu_20_R_0_speedUp_10 <-
    current_las_normalized %>% 
    st_li2012(plot_boundary = current_plot_boundary, dt1 = 1.0, dt2 = 1.0, R = 0, Zu = 20, hmin = min_height, speed_up = 10)
  
  ttops_li2012_dt1_1.0_dt2_1.3_zu_20_R_0_speedUp_20 <-
    current_las_normalized %>% 
    st_li2012(plot_boundary = current_plot_boundary, dt1 = 1.0, dt2 = 1.3, R = 0, Zu = 20, hmin = min_height, speed_up = 20)
  
  
  
  # watershed treetop detection ---------------------------------------------
  
  
  # th_tree: Threshold below which a pixel cannot be a tree
  # tolerance: The minimum height of the object in the units of image intensity between 
  # its highest point (seed) and the point where it contacts another object (checked for 
  # every contact pixel). If the height is smaller than the tolerance, the object will 
  # be combined with one of its neighbors, which is the highest. Tolerance should be 
  # chosen according to the range of x. Default value is 1, which is a reasonable value 
  # if x comes from distmap.
  # ext: Radius of the neighborhood in pixels for the detection of neighboring objects. 
  # Higher value smoothes out small objects.
  
  ttops_watershed_tol_3_ext_0.5 <-
    current_las_normalized %>% 
    st_watershed(plot_boundary = current_plot_boundary, chm = current_chm, th_tree = 3, tol = 3, ext = 0.5 * round((1 / res(current_chm)[1])))
  
  
  ttops_watershed_tol_1_ext_0.5 <-
    current_las_normalized %>% 
    st_watershed(plot_boundary = current_plot_boundary, chm = current_chm, th_tree = 3, tol = 1, ext = 0.5 * round((1 / res(current_chm)[1])))
  
  
  ttops_watershed_tol_1_ext_1.0 <-
    current_las_normalized %>% 
    st_watershed(plot_boundary = current_plot_boundary, chm = current_chm, th_tree = 3, tol = 1, ext = 1 * round((1 / res(current_chm)[1])))
  
  
  
  ## lidRplugins experimental tree detection
  
  # ptree treetop detection -------------------------------------------------
  
  # ptrees
  # importantly, we have to use the non-normalized point cloud for this algorithm
  
  ttops_ptrees_hmin_2_k_30_15_nmax_7 <-
    current_las %>%
    st_ptree(plot_boundary = current_plot_boundary, k = c(30, 15), hmin = min_height, nmax = 7)
  
  
  # multichm treetop detection ----------------------------------------------
  # Eysn, L., Hollaus, M., Lindberg, E., Berger, F., Monnet, J. M., Dalponte, M., … Pfeifer, N. (2015). A benchmark of lidar-based single tree detection methods using heterogeneous forest data from the Alpine Space. Forests, 6(5), 1721–1747. https://doi.org/10.3390/f6051721
  
  ttops_multichm_res_1.0_thickness_0.5_dist_2d_3.0_dist_3d_5_usemax_FALSE_ws_5 <-
    current_las_normalized %>% 
    st_multichm(plot_boundary = current_plot_boundary, res = 1.0, layer_thickness = 0.5, dist_2d = 3.0, dist_3d = 5.0, use_max = FALSE, ws = 5)
  
  
  # lmfx treetop detection --------------------------------------------------
  # lmfx
  

  ttops_lmfx_dist_2d_3.0_ws_3.0 <- 
    current_las_normalized %>% 
    st_lmfx(plot_boundary = current_plot_boundary, hmin = min_height, dist_2d = 3.0, ws = 3.0)
  
  # combine treetop detection algorithm outputs -----------------------------
  
  
  # All the different ttops
  
  ttops <-
    list(
      vwf_default = ttops_vwf_default,
      vwf_pines = ttops_vwf_pines,
      vwf_combined = ttops_vwf_combined,
      localMaxima_chm_1.5 = ttops_localMaxima_chm_1.5,
      localMaxima_las_1.5 = ttops_localMaxima_las_1.5,
      localMaxima_chm_2 = ttops_localMaxima_chm_2,
      localMaxima_las_2 = ttops_localMaxima_las_2,
      localMaxima_chm_2.5 = ttops_localMaxima_chm_2.5,
      localMaxima_las_2.5 = ttops_localMaxima_las_2.5,
      li2012_dt1_1.4_dt2_1.4_zu_15_R_0_speedUp_10 = ttops_li2012_dt1_1.4_dt2_1.4_zu_15_R_0_speedUp_10,
      li2012_dt1_1.5_dt2_2.0_zu_15_R_0_speedUp_10 = ttops_li2012_dt1_1.5_dt2_2.0_zu_15_R_0_speedUp_10,
      li2012_dt1_1.0_dt2_2.0_zu_15_R_0_speedUp_10 = ttops_li2012_dt1_1.0_dt2_2.0_zu_15_R_0_speedUp_10,
      li2012_dt1_1.0_dt2_2.0_zu_20_R_0_speedUp_10 = ttops_li2012_dt1_1.0_dt2_2.0_zu_20_R_0_speedUp_10,
      li2012_dt1_1.0_dt2_1.5_zu_25_R_0_speedUp_10 = ttops_li2012_dt1_1.0_dt2_1.5_zu_25_R_0_speedUp_10,
      li2012_dt1_1.0_dt2_1.5_zu_20_R_0_speedUp_10 = ttops_li2012_dt1_1.0_dt2_1.5_zu_20_R_0_speedUp_10,
      li2012_dt1_1.0_dt2_1.4_zu_20_R_0_speedUp_10 = ttops_li2012_dt1_1.0_dt2_1.4_zu_20_R_0_speedUp_10,
      li2012_dt1_1.0_dt2_1.4_zu_25_R_0_speedUp_10 = ttops_li2012_dt1_1.0_dt2_1.4_zu_25_R_0_speedUp_10,
      li2012_dt1_1.0_dt2_1.3_zu_20_R_0_speedUp_10 = ttops_li2012_dt1_1.0_dt2_1.3_zu_20_R_0_speedUp_10,
      li2012_dt1_1.0_dt2_1.0_zu_20_R_0_speedUp_10 = ttops_li2012_dt1_1.0_dt2_1.0_zu_20_R_0_speedUp_10,
      li2012_dt1_1.0_dt2_1.3_zu_20_R_0_speedUp_20 = ttops_li2012_dt1_1.0_dt2_1.3_zu_20_R_0_speedUp_20,
      watershed_tol_3_ext_0.5 = ttops_watershed_tol_3_ext_0.5,
      watershed_tol_1_ext_0.5 = ttops_watershed_tol_1_ext_0.5,
      watershed_tol_1_ext_1.0 = ttops_watershed_tol_1_ext_1.0,
      ptrees_hmin_2_k_30_15_nmax_7 = ttops_ptrees_hmin_2_k_30_15_nmax_7,
      multichm_res_1.0_thickness_0.5_dist_2d_3.0_dist_3d_5_usemax_FALSE_ws_5 = ttops_multichm_res_1.0_thickness_0.5_dist_2d_3.0_dist_3d_5_usemax_FALSE_ws_5,
      lmfx_dist_2d_3.0_ws_3.0 = ttops_lmfx_dist_2d_3.0_ws_3.0)
  
  
  # Combine all of the sf objects from the different treetop detection algorithms
  # into a single sf object
  # Also add a column that says which algorithm detected the tree
  # Further, calculate the nearest neighbor distance between the
  # three nearest neighbors for each tree
  ttops_sf <-
    lapply(seq_along(ttops), FUN = function(j) {

      # Each list element in the ttops list is itself a list of two or three pieces:
      # the sf object with the location and height of each tree segmented
      # the elapsed time that the algorithm took to complete
      # potentially the returned las point cloud if the tree detection algorithm
      # returns that as well
      
      # which of those 2-3 list element in the current ttop list is the sf object?
      ttops_sf_idx <- which(names(ttops[[j]]) == "ttops_sf")

      # just select the height and geometry columns of the sf object
      ttops[[j]][[ttops_sf_idx]] <-
        ttops[[j]][[ttops_sf_idx]] %>%
        dplyr::select(height, geometry)
      
      # If no trees were segmented, create an "empty" sf object so that the rest of
      # the objects will bind together nicely since they'll have the same column names, etc.
      if (nrow(ttops[[j]][[ttops_sf_idx]]) == 0) {
        ttops[[j]][[ttops_sf_idx]] <- 
          st_sf(height = NA, geometry = st_sfc(st_multipoint()), crs = st_crs(ttops[[j]][[ttops_sf_idx]]))}
      
      # Add the treetop detection algorithm name to the data frame row
      ttops[[j]][[ttops_sf_idx]] <-
        ttops[[j]][[ttops_sf_idx]] %>% 
        dplyr::mutate(ttop_method = names(ttops[j]))
      
      # If there are any trees in the sf object, get the nearest neighbor distance between each tree
      # and its three nearest neighbors
      if(nrow(ttops[[j]][[ttops_sf_idx]]) > 0) {
        nn <- st_nn(ttops[[j]][[ttops_sf_idx]], ttops[[j]][[ttops_sf_idx]], k = min(c(4, nrow(ttops[[j]][[ttops_sf_idx]]))), returnDist = TRUE, sparse = FALSE, progress = FALSE)$dist %>% 
          as.data.frame() } else {nn <- data.frame()}
      
      # remove nuisance column
      nn <- nn %>% dplyr::select(-1)
      
      # build up a nearest neighbors dataframe, but ensure that the data frame will always have 
      # three columns (even if there are not enough trees segmented to result in 3
      # nearest neighbor distances for each tree). If there aren't enough trees, fill the
      # first, first and second, or first, second and third columns of this data frame with NA
      if (ncol(nn) == 0) {nn <- nn %>% mutate(nn1 = NA, nn2 = NA, nn3 = NA)} else
        if (ncol(nn) == 1) {nn <- nn %>% rename(nn1 = V2) %>% mutate(nn2 = NA, nn3 = NA)} else
          if (ncol(nn) == 2) {nn <- nn %>% rename(nn1 = V2, nn2 = V3) %>% mutate(nn3 = NA)} else
            if (ncol(nn) == 3) {nn <- nn %>% rename(nn1 = V2, nn2 = V3, nn3 = V4) }
      
      # bind the nearest neighbor dataframe we just created with the sf object
      ttops[[j]][[ttops_sf_idx]] <-
        ttops[[j]][[ttops_sf_idx]] %>%
        bind_cols(nn)
      
      # transform the crs to epsg4326 to make it compatible across all segmented sf objects
      # if the sf object is empty, set its crs to 4326
      if (nrow(ttops[[j]][[ttops_sf_idx]]) > 0) {
        ttops[[j]][[ttops_sf_idx]] <- 
          ttops[[j]][[ttops_sf_idx]] %>% 
          st_transform(4326)} else 
            {ttops[[j]][[ttops_sf_idx]] <- suppressWarnings(ttops[[j]][[ttops_sf_idx]] %>% st_set_crs(4326))}
      
    }) %>%
    do.call("rbind", .) # bind the rows for all the sf objects resulting from different segmentation algorithms together

  # Get the time elapsed for each method-- this might help break any close calls
  # with respect to which algorithm best captured the data collected on the
  # ground
  ttops_time <-
    lapply(seq_along(ttops), FUN = function(k) {
      ttops_time_idx <- which(names(ttops[[k]]) == "ttops_time")
      ttops_time <- ttops[[k]][[ttops_time_idx]]
      
      data_frame(ttops_method = names(ttops[k]), elapsed_time = ttops_time)
    }) %>% do.call("rbind", .)
  
  
  current_plot_ttops_summary <-
    ttops %>%
    dplyr::group_by(ttop_method) %>%
    dplyr::summarize(plot = current_plot,
                     total_tree_count = length(which(!is.na(height))),
                     live_tree_count = NA,
                     dead_tree_count = NA,
                     total_density_tph = (total_tree_count * 10000) / (pi * plot_radius ^ 2),
                     live_density_tph = NA,
                     dead_density_tph = NA,
                     live_proportion = NA,
                     height_lwr_25 = quantile(height, prob = 0.25, na.rm = TRUE),
                     height_mean = mean(height),
                     height_upr_25 = quantile(height, prob = 0.75, na.rm = TRUE),
                     nn_1_mean = mean(nn1),
                     nn_2_mean = mean(nn2),
                     nn_3_mean = mean(nn3),
                     tree_count_above_15m = sum(height >= 15),
                     tree_count_below_15m = sum(height < 15)) %>%
    as.data.frame() %>%
    dplyr::select(-geometry) %>%
    dplyr::select(plot, everything())
  
  
  rbind(ground_tree_summary %>% filter(plot == current_plot), current_plot_ttops_summary)
  
  # cat(paste0("...", all_validation_plots[i], " (", i, " of ", length(all_validation_plots), ") complete.\n"))
  # 
  # # Crown Segmentation as a 2nd step
  # # Method 1
  # EBImage watershed
  # crowns_watershed
  # 
  # # Method 2
  # # Li algorithm
  # # Currently not good final crown objects
  # crowns_li2_dt1_1.4_dt2_1.4
  # crowns_li2_dt1_1.5_dt2_2
  # 
  # # Methods 3 through 5 rely on tree top detection ahead of time
  # # Method 3
  # # Dalponte
  # lastrees_dalponte(current_las_normalized, current_chm, ttops, th_tree = 2, th_seed = 0.45,
  #                   th_cr = 0.55, max_cr = 10, extra = TRUE)
  # 
  # # Method 4
  # # Silva
  # lastrees_silva(current_las_normalized, current_chm, ttops, max_cr_factor = 0.6, exclusion = 0.3,
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
stopCluster(cl)
end <- Sys.time()
end - start

glimpse(ttops_summary)

if (!file.exists("data/data_output/ttops-summary.csv")) {
  write_csv(ttops_summary, path = "data/data_output/ttops-summary.csv")
}

ttops_summary_from_file <- read_csv("data/data_output/ttops-summary.csv")

ttops_summary_to_save <-
  ttops_summary %>%
  filter(!(ttop_method %in% ttops_summary_from_file$ttop_method)) %>% 
  rbind(ttops_summary_from_file)

# Overwrite the originally saved data with the data.frame representing the old data plus the new data
# write_csv(ttops_summary_to_save, path = "data/data_output/ttops-summary.csv")

ttops_summary <- ttops_summary_from_file
# ttops_summary <- ttops_summary_to_save
ground <- ttops_summary %>% filter(ttop_method == "ground") %>% dplyr::select(-live_tree_count, -dead_tree_count, -total_density_tph, -live_density_tph, -dead_density_tph, -live_proportion)
air <- ttops_summary %>% filter(ttop_method != "ground") %>% dplyr::select(-live_tree_count, -dead_tree_count, -total_density_tph, -live_density_tph, -dead_density_tph, -live_proportion)

ground <- gather(ground, key = forest_metric, value = ground_value, -plot, -ttop_method) %>% dplyr::select(-ttop_method)
air <- gather(air, key = forest_metric, value = air_value, -plot, -ttop_method)

air_ground <- left_join(air, ground, by = c("plot", "forest_metric"))
air_ground$diff <- air_ground$air_value - air_ground$ground_value
air_ground

air_ground %>%
  filter(forest_metric == "total_tree_count") %>% 
  filter(diff > 10 | diff < -10) %>% 
  group_by(plot, "pos_or_neg" = ifelse(diff > 0, yes = "positive", no = "negative")) %>% 
  summarize(n = n()) %>% 
  as.data.frame()

ggplot(air_ground, aes(x = ground_value, y = diff)) +
  geom_point() +
  geom_smooth(aes(col = ttop_method)) +
  facet_wrap(~forest_metric, scales = "free") +
  scale_color_viridis_d()

air_ground_summary <-
  air_ground %>%
  group_by(ttop_method, forest_metric) %>% 
  summarize(ground_correlation = cor(air_value, ground_value, use = "complete.obs"),
            rmse = sqrt(mean((air_value - ground_value)^2)),
            me = mean(air_value - ground_value),
            med_error = median(air_value - ground_value),
            min = min(air_value - ground_value, na.rm = TRUE),
            max = max(air_value - ground_value, na.rm = TRUE)) %>% 
  as.data.frame()

unique(air_ground_summary$forest_metric)
air_ground_summary %>% filter(forest_metric == "total_tree_count") %>% arrange(desc(ground_correlation))
air_ground_summary %>% filter(forest_metric == "total_tree_count") %>% arrange(rmse)

air_ground_summary %>% filter(forest_metric == "height_mean")
air_ground_summary %>% filter(forest_metric == "nn_1_mean") %>% arrange(desc(ground_correlation))
air_ground_summary %>% filter(forest_metric == "tree_count_above_15m") %>% arrange(rmse)
air_ground_summary %>% filter(forest_metric == "tree_count_below_15m") %>% arrange(rmse)
