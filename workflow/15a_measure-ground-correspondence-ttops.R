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
library(furrr)

source(here::here("data/data_carpentry/segmentation-helper-functions.R"))

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

# These sites had X3 and RedEdge photos merged into the same project, so we look in a different place for some of the relevant
# files.
merged_sites <- c("eldo_3k_2",
                  "eldo_3k_3",
                  "eldo_4k_2")

# This is where I can put in sites that need their processing redone. An empty 
# string means that no already-processed site output will be overwritten
# (but sites that have yet to be processed will still have their processing done)
if(file.exists(here::here("analyses/analyses_output/ground-tree-summary.csv"))) {
  ground_tree_summary <- readr::read_csv(here::here(paste0("analyses/analyses_output/ground-tree-summary.csv")))
} else {
  source(here::here(paste0("analyses/summarize-ground-data.R")))
}

all_validation_plots <- ground_tree_summary$plot

# Start the timer
start <- Sys.time()
# Set up the parallelization
num_cores_to_use <- availableCores() - 2
plan(multiprocess, workers = num_cores_to_use)

suppressWarnings(suppressMessages( # Suppress all the outputs from lidR functions. Too much red text!
ttops_summary <- 
  all_validation_plots %>% 
  furrr::future_map(.f = function(current_plot) {

  # get data for particular plot --------------------------------------------
  print(current_plot)
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
  
  # The dtm is the terrain model for a particular site. We need it to normalize the objects returned from the
  # ptrees algorithm, which acts on a non-normalized point cloud.
  current_dtm <- raster::raster(paste0(current_dir, current_site, "_plot-remote-data/", current_plot, "_dtm.tif"))
  
  # The Canopy Height Model (chm) is the dsm (vegetation + ground) minus the dtm (ground)
  # to give just the height of the vegetation.
  current_chm_rough <- raster::raster(paste0(current_dir, current_site, "_plot-remote-data/", current_plot, "_chm.tif"))
  
  # The point cloud is directly used for some segmentation algorithms, so we import that too
  current_las <- lidR::readLAS(paste0(current_dir, current_site, "_plot-remote-data/", current_plot, "_classified-point-cloud.las"))
  # current_las_normalized <- lidR::lasnormalize(las = current_las, algorithm = tin())
  current_las_normalized <- lidR::lasnormalize(las = current_las, algorithm = current_dtm, na.rm = TRUE)
  
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
  dynamicWindow_default <- function(x){
    window_radius <- x * 0.06 + 0.5
    window_radius
  }

  ttops_vwf_default <-
    current_chm %>%
    st_vwf(plot_boundary = current_plot_boundary, winFun = dynamicWindow_default, minHeight = min_height, maxWinDiameter = NULL)

  # Equation and coefficients taken from Popescu and Wynne (2004)
  # Divide the Popescu and Wynne (2004) equations by 2 to convert to *radius* of search window
  # which is what the vwf() function requires

  # Equation and coefficients taken from Popescu and Wynne (2004)'s "Pines model"
  dynamicWindow_pines <- function(x) {
    window_radius <- 0.5 * (3.75105 - 0.17919*x + 0.01241 * x^2)
    window_radius
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
    window_radius
  }

  ttops_vwf_combined <-
    current_chm %>%
    st_vwf(plot_boundary = current_plot_boundary, winFun = dynamicWindow_combined, minHeight = min_height, maxWinDiameter = NULL)

  vwf_list <-
  list(#vwf_default = ttops_vwf_default,
       vwf_pines = ttops_vwf_pines,
       vwf_combined = ttops_vwf_combined)
  
  # lmf on a chm treetop detection ---------------------------------------------------
  ws_vals <- c(1.5, 2.0, 2.5)
  ws_in_pixels_vals <- round(ws_vals * (1 / res(current_chm)[1]))
  ws_in_pixels_vals <- ifelse(ws_in_pixels_vals %% 2 == 0, yes = ws_in_pixels_vals + 1, no = ws_in_pixels_vals)
  ws_names <- paste("ws", ws_vals, sep = "_")

  # First row is default values for li2012 algorithm
  # second row come from Jakubowski et al. (2013) [mixed conifer forest near Tahoe-- pretty comparable to our study]
  lmf_params <- data_frame(ws = ws_vals, ws_in_pixels_vals, ws_names)

  lmf_chm_list <-
    lmf_params %>%
    pmap(.f = function(ws_in_pixels_vals, ...) {

      current_chm %>%
        st_lmf(plot_boundary = current_plot_boundary, ws = ws_in_pixels_vals)

    })

  lmf_chm_names <-
    lmf_params %>%
    tidyr::unite(col = "ttops_method", ends_with("names")) %>%
    pull(ttops_method) %>%
    paste("localMaxima_chm", ., sep = "_")

  names(lmf_chm_list) <- lmf_chm_names



# lmf on a las treetop detection ------------------------------------------


  lmf_las_list <-
    lmf_params %>%
    pmap(.f = function(ws_vals, ...) {

      current_las_normalized %>%
        st_lmf(plot_boundary = current_plot_boundary, ws = ws_vals)

    })

  lmf_las_names <-
    lmf_params %>%
    tidyr::unite(col = "ttops_method", ends_with("names")) %>%
    pull(ttops_method) %>%
    paste("localMaxima_las", ., sep = "_")

  names(lmf_las_list) <- lmf_las_names


  
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
  
  # First round used a few deliberate combinations of parameters (one from the literature)
  dt1_vals <- c(1.4, 1.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
  dt2_vals <- c(1.4, 2.0, 2.0, 2.0, 1.5, 1.5, 1.4, 1.4, 1.3, 1.5, 1.3)
  R_vals <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  Zu_vals <- c(15, 15, 15, 20, 25, 20, 20, 25, 20, 25, 20)
  speed_up <- c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 20)

  dt1_names <- paste("dt1", dt1_vals, sep = "_")
  dt2_names <- paste("dt2", dt2_vals, sep = "_")
  R_names <- paste("R", R_vals, sep = "_")
  Zu_names <- paste("zu", Zu_vals, sep = "_")
  speed_up_names <- paste("speedUp", speed_up, sep = "_")

  # First row is default values for li2012 algorithm
  # second row come from Jakubowski et al. (2013) [mixed conifer forest near Tahoe-- pretty comparable to our study]
  li2012_params <- data_frame(dt1 = dt1_vals, dt2 = dt2_vals, R = R_vals, Zu = Zu_vals, hmin = min_height, speed_up = speed_up,
                              dt1_names, dt2_names, Zu_names, R_names, speed_up_names)


  # Second round was systematically exploring the parameter space
  dt1_vals <- c(1, 1.5, 2, 2.5)
  dt2_vals <- c(dt1_vals + 0.5)
  R_vals <- c(0, 2)
  Zu_vals <- c(15, 20, 25)
  speed_up <- c(10, 15)

  li2012_params <-
    as_data_frame(expand.grid(dt1_vals, dt2_vals, R_vals, Zu_vals, speed_up)) %>%
    setNames(c("dt1", "dt2", "R", "Zu", "speed_up")) %>%
    filter(dt2 > dt1) %>%
    dplyr::mutate(dt1_names = paste("dt1", dt1, sep = "_")) %>%
    dplyr::mutate(dt2_names = paste("dt2", dt2, sep = "_")) %>%
    dplyr::mutate(R_names = paste("R", R, sep = "_")) %>%
    dplyr::mutate(Zu_names = paste("zu", Zu, sep = "_")) %>%
    dplyr::mutate(speed_up_names = paste("speedUp", speed_up, sep = "_"))

    li2012_list <-
    li2012_params %>%
    pmap(.f = function(dt1, dt2, R, Zu, speed_up, ...) {

      current_las_normalized %>%
      st_li2012(plot_boundary = current_plot_boundary, dt1 = dt1, dt2 = dt2, R = R, Zu = Zu, hmin = min_height, speed_up = speed_up)

    })

  li2012_names <-
    li2012_params %>%
    tidyr::unite(col = "ttops_method", ends_with("names")) %>%
    pull(ttops_method) %>%
    paste0("li2012_", .)

  names(li2012_list) <- li2012_names

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
  th_tree_vals <- min_height
  tol_vals <- c(3, 1, 1)
  tol_names <- paste("tol", tol_vals, sep = "_")
  ext_vals <- c(0.5, 0.5, 1.0)
  ext_names <- paste("ext", ext_vals, sep = "_")

  watershed_params <- data_frame(th_tree = th_tree_vals, tol = tol_vals, ext = ext_vals, tol_names, ext_names)
  # First row is default values for li2012 algorithm
  # second row come from Jakubowski et al. (2013)

  watershed_list <-
    watershed_params %>%
    pmap(.f = function(th_tree, tol, ext, ...) {
       current_las_normalized %>%
        st_watershed(plot_boundary = current_plot_boundary, chm = current_chm, th_tree = th_tree, tol = tol, ext = ext * round((1 / res(current_chm)[1])))
    })

  watershed_names <-
    watershed_params %>%
    tidyr::unite(col = "ttops_method", ends_with("names")) %>%
    dplyr::pull(ttops_method) %>%
    paste0("watershed_", .)

  names(watershed_list) <- paste0("watershed_", watershed_names)

  # lidRplugins experimental tree detection
  
  # ptree treetop detection -------------------------------------------------
  
  # ptrees
  # importantly, we have to use the non-normalized point cloud for this algorithm

  ttops_ptrees_k_30_15_nmax_7 <-
    current_las %>%
    st_ptree(plot_boundary = current_plot_boundary, dtm = current_dtm, k = c(30, 15), algorithm_hmin = -Inf, post_hmin = min_height, nmax = 7)

  ttops_ptrees_k_vega2014_nmax_7 <-
    current_las %>%
    st_ptree(plot_boundary = current_plot_boundary, dtm = current_dtm, k = c(100, 80, 60, 40, 30, 25, 20, 15, 12, 10, 8, 7, 6, 5), algorithm_hmin = -Inf, post_hmin = min_height, nmax = 7)


  ttops_ptrees_k_vega2014_trunc_nmax_7 <-
    current_las %>%
    st_ptree(plot_boundary = current_plot_boundary, dtm = current_dtm, k = c(12, 10, 8, 7, 6, 5), algorithm_hmin = -Inf, post_hmin = min_height, nmax = 7)


  ptrees_list <-
  list(ptrees_k_30_15_nmax_7 = ttops_ptrees_k_30_15_nmax_7,
       ptrees_k_vega2014_nmax_7 = ttops_ptrees_k_vega2014_nmax_7,
       ptrees_k_vega2014_trunc_nmax_7 = ttops_ptrees_k_vega2014_trunc_nmax_7)

  # multichm treetop detection ----------------------------------------------
  # Eysn, L., Hollaus, M., Lindberg, E., Berger, F., Monnet, J. M., Dalponte, M., … Pfeifer, N. (2015). A benchmark of lidar-based single tree detection methods using heterogeneous forest data from the Alpine Space. Forests, 6(5), 1721–1747. https://doi.org/10.3390/f6051721
  
  ttops_multichm_res_1.0_thickness_0.5_dist2d_3_dist3d_5_usemax_FALSE_ws_5 <-
    current_las_normalized %>%
    st_multichm(plot_boundary = current_plot_boundary, res = 1.0, layer_thickness = 0.5, dist_2d = 3.0, dist_3d = 5.0, use_max = FALSE, ws = 5)

  multichm_list <-
    list(multichm_res_1.0_thickness_0.5_dist2d_3_dist3d_5_usemax_FALSE_ws_5 = ttops_multichm_res_1.0_thickness_0.5_dist2d_3_dist3d_5_usemax_FALSE_ws_5)

  # lmfx treetop detection --------------------------------------------------
  
  dist_2d_vals <- 3.0
  ws_vals <- 3.0

  dist_2d_names <- paste("dist2d", dist_2d_vals, sep = "_")
  ws_names <- paste("ws", ws_vals, sep = "_")

  lmfx_params <- dplyr::data_frame(dist2d = dist_2d_vals, ws = ws_vals, dist_2d_names, ws_names)

  # Recall the dynamic window from when we used the vwf() function on its own
  # dynamicWindow_pines <- function(x) {
  #   window_radius <- 0.5 * (3.75105 - 0.17919*x + 0.01241 * x^2)
  #   return(window_radius)
  # }

  dist_2d_vals <- c(1.0, 1.5, 2.0, 2.5, 3)
  ws_vals <- list("ws_1.5" = 1.5, "ws_2" = 2, "ws_2.5" = 2.5, "ws_3" = 3, ws_dynamicWindow_pines = dynamicWindow_pines, ws_dynamicWindow_combined = dynamicWindow_combined)

  lmfx_params <-
    dplyr::as_data_frame(expand.grid(dist_2d_vals, ws_vals)) %>%
    setNames(c("dist2d", "ws")) %>%
    dplyr::mutate(dist2d_names = paste("dist2d", dist2d, sep = "_")) %>%
    dplyr::mutate(ws_names = names(ws))

   lmfx_list <-
    lmfx_params %>%
    purrr::pmap(.f = function(dist2d, ws, ...) {
      current_las_normalized %>%
        st_lmfx(plot_boundary = current_plot_boundary, hmin = min_height, dist_2d = dist2d, ws = ws)
    })

  lmfx_names <-
    lmfx_params %>%
    tidyr::unite(col = "ttops_method", ends_with("names")) %>%
    dplyr::pull(ttops_method) %>%
    paste0("lmfx_", .)

  names(lmfx_list) <- lmfx_names
  
  # combine treetop detection algorithm outputs -----------------------------
  
  
  # All the different ttops
  
  ttops <-
    list(
      vwf_list,
      lmf_chm_list,
      lmf_las_list,
      li2012_list,
      watershed_list,
      ptrees_list,
      multichm_list,
      lmfx_list
      ) %>%
    purrr::flatten()
  
  
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
        dplyr::mutate(ttops_method = names(ttops[j]))
      
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
      if (ncol(nn) == 0) {nn <- nn %>% dplyr::mutate(nn1 = NA, nn2 = NA, nn3 = NA)} else
        if (ncol(nn) == 1) {nn <- nn %>% dplyr::rename(nn1 = V2) %>% dplyr::mutate(nn2 = NA, nn3 = NA)} else
          if (ncol(nn) == 2) {nn <- nn %>% dplyr::rename(nn1 = V2, nn2 = V3) %>% dplyr::mutate(nn3 = NA)} else
            if (ncol(nn) == 3) {nn <- nn %>% dplyr::rename(nn1 = V2, nn2 = V3, nn3 = V4) }
      
      # bind the nearest neighbor dataframe we just created with the sf object
      ttops[[j]][[ttops_sf_idx]] <-
        ttops[[j]][[ttops_sf_idx]] %>%
        dplyr::bind_cols(nn)
      
      # transform the crs to epsg4326 to make it compatible across all segmented sf objects
      # if the sf object is empty, set its crs to 4326
      if (nrow(ttops[[j]][[ttops_sf_idx]]) > 0) {
        ttops[[j]][[ttops_sf_idx]] <- 
          ttops[[j]][[ttops_sf_idx]] %>% 
          sf::st_transform(4326)} else 
            {ttops[[j]][[ttops_sf_idx]] <- suppressWarnings(ttops[[j]][[ttops_sf_idx]] %>% sf::st_set_crs(4326))}
      
    }) %>%
    do.call("rbind", .) # bind the rows for all the sf objects resulting from different segmentation algorithms together

  # Get the time elapsed for each method-- this might help break any close calls
  # with respect to which algorithm best captured the data collected on the
  # ground
  ttops_time <-
    lapply(seq_along(ttops), FUN = function(k) {
      ttops_time_idx <- which(names(ttops[[k]]) == "ttops_time")
      ttops_time <- ttops[[k]][[ttops_time_idx]]
      
      dplyr::data_frame(ttops_method = names(ttops[k]), elapsed_time = ttops_time)
    }) %>% do.call("rbind", .)
  
  
  current_plot_ttops_summary <-
    ttops_sf %>%
    dplyr::group_by(ttops_method) %>%
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
    dplyr::left_join(ttops_time, by = "ttops_method") %>% 
    dplyr::select(plot, ttops_method, elapsed_time, everything())
})
))
 
ttops_summary <- 
  do.call("rbind", ttops_summary) %>% 
  rbind(ground_tree_summary)

end <- Sys.time()
end - start

glimpse(ttops_summary)

if (!file.exists("analyses/analyses_output/ttops-summary.csv")) {
  write_csv(ttops_summary, path = "analyses/analyses_output/ttops-summary.csv")
}

ttops_summary_from_file <- read_csv("analyses/analyses_output/ttops-summary.csv")

ttops_summary_to_save <-
  ttops_summary %>%
  filter(!(ttops_method %in% ttops_summary_from_file$ttops_method)) %>%
  rbind(ttops_summary_from_file) %>% 
  dplyr::mutate(elapsed_time = as.numeric(elapsed_time))

# Overwrite the originally saved data with the data.frame representing the old data plus the new data
write_csv(ttops_summary_to_save, path = "analyses/analyses_output/ttops-summary.csv")
