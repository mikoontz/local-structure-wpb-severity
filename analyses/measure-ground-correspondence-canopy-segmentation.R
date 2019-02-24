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
library(stars)

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
current_plot = all_validation_plots[1]

# Start the timer
start <- Sys.time()
# Set up the parallelization
num_cores_to_use <- availableCores() - 2
plan(multiprocess, workers = num_cores_to_use)

# suppressWarnings(suppressMessages( # Suppress all the outputs from lidR functions. Too much red text!
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
      
      
      # Equation and coefficients taken from Popescu and Wynne (2004)'s "Pines model"
      dynamicWindow_pines <- function(x) {
        window_radius <- 0.5 * (3.75105 - 0.17919*x + 0.01241 * x^2)
        window_radius
      }
      
      # Li2012 algorithm for treetop detection ----------------------------------
      
      
      dt1_vals <- c(1, 1, 1, 1, 1, 1, 1, 1.5)
      dt2_vals <- c(2, 1.5, 2, 1.3, 1.3, 1.5, 1.4, 2)
      R_vals <- 0
      Zu_vals <- c(20, 20, 15, 20, 20, 25, 20, 25)
      speed_up <- c(10, 10, 10, 10, 20, 10, 10, 10)
      
      dt1_names <- paste("dt1", dt1_vals, sep = "_")
      dt2_names <- paste("dt2", dt2_vals, sep = "_")
      R_names <- paste("R", R_vals, sep = "_")
      Zu_names <- paste("zu", Zu_vals, sep = "_")
      speed_up_names <- paste("speedUp", speed_up, sep = "_")
      
      # First row is default values for li2012 algorithm
      # second row come from Jakubowski et al. (2013) [mixed conifer forest near Tahoe-- pretty comparable to our study]
      li2012_params <- data_frame(dt1 = dt1_vals, dt2 = dt2_vals, R = R_vals, Zu = Zu_vals, hmin = min_height, speed_up = speed_up,
                                  dt1_names, dt2_names, Zu_names, R_names, speed_up_names)
      
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
      
      # lmfx treetop detection --------------------------------------------------
      
      dist_2d_vals <- c(1.5, 1, 2.5, 1, 2.5, 3, 2.5, 2.5, 2, 2, 3, 3)
      ws_vals <- list("ws_2.5" = 2.5, "ws_2.5" = 2.5, "ws_2.5" = 2.5, ws_dynamicWindow_pines = dynamicWindow_pines, "ws_3" = 3, "ws_2.5" = 2.5, "ws_2" = 2, ws_dynamicWindow_pines = dynamicWindow_pines, "ws_1.5" = 1.5, "ws_3" = 3, "ws_1.5" = 1.5, "ws_2" = 2)
      
      dist2d_names <- paste("dist2d", dist_2d_vals, sep = "_")
      
      lmfx_params <- dplyr::data_frame(dist2d = dist_2d_vals, ws = ws_vals, dist2d_names, ws_names = names(ws_vals))
      
      # Recall the dynamic window from when we used the vwf() function on its own
      # dynamicWindow_pines <- function(x) {
      #   window_radius <- 0.5 * (3.75105 - 0.17919*x + 0.01241 * x^2)
      #   return(window_radius)
      # }
      
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
      
      ttops <-
        list(
          li2012_list,
          lmfx_list
        ) %>%
        purrr::flatten()
      
      ttops_sf <-
      lapply(seq_along(ttops), FUN = function(j) {
        # which of those 2-3 list element in the current ttop list is the sf object?
        ttops_sf_idx <- which(names(ttops[[j]]) == "ttops_sf")
        
        # just select the height and geometry columns of the sf object
        ttops[[j]][[ttops_sf_idx]] <-
          ttops[[j]][[ttops_sf_idx]] %>%
          dplyr::select(treeID, height, geometry) %>% 
          dplyr::mutate(treeID = 1:nrow(.))
      })
      
      ttops_las <- ttops[which(str_detect(names(ttops), pattern = "li2012"))] 
      
      ttops_las <-
        lapply(seq_along(ttops_las), FUN = function(j) {
          # which of those 2-3 list element in the current ttop list is the sf object?
          ttops_las_idx <- which(names(ttops_las[[j]]) == "ttops_las")
          
          ttops_las[[j]][[ttops_las_idx]] <-
            ttops_las[[j]][[ttops_las_idx]]
        })
      
      tic()
      mcws_segmentation1 <-
        ttops_sf %>% 
        purrr::map(.f = function(ttops) {
          ttops %>%
            as("Spatial") %>% 
            ForestTools::mcws(CHM = current_chm, minHeight = 1.5, format = "raster") %>% 
            st_as_stars() %>% 
            st_as_sf(merge = TRUE)
        })
      toc()
      
      tic()
      mcws_segmentation2 <-
        ttops_sf %>% 
        purrr::map(.f = function(ttops) {
          ttops %>%
            as("Spatial") %>% 
            ForestTools::mcws(CHM = current_chm, minHeight = 1.5, format = "polygons", OSGeoPath = "C:/OSGeo4W64") %>% 
            st_as_sf()
        })
      toc()
      
      mcws_segmentation3 <- 
        ForestTools::mcws(treetops = ttops_sf[[1]], CHM = current_chm, minHeight = 1.5, format = "raster") %>% 
        setNames(nm = "treeID") %>% 
        st_as_stars() %>% 
        st_as_sf(as_points = TRUE)
      
      mcwatershed_segmentation <- 
        lapply(ttops_sf, FUN = function(ttops) {
          lidR::mcwatershed(chm = current_chm, treetops = as(ttops, "Spatial"), th_tree = 1.5)() %>% 
            APfun::APpolygonize() %>% 
            st_as_sf()
        })
      
      # tolerance: The minimum height of the object in the units of image intensity between its highest point (seed) and the point where it contacts another object (checked for every contact pixel). If the height is smaller than the tolerance, the object will be combined with one of its neighbors, which is the highest. Tolerance should be chosen according to the range of x. Default value is 1, which is a reasonable value if x comes from distmap.
      # ext: Radius of the neighborhood in pixels for the detection of neighboring objects. Higher value smoothes out small objects.
      watershed_segmentation <- 
        lidR::watershed(chm = current_chm, th_tree = 1.5, tol = 2, ext = ceiling(2 * (1 / res(current_chm))[1]))() %>% 
        APfun::APpolygonize() %>% 
        st_as_sf() %>% 
        st_set_crs(proj4string(current_chm))
      
      watershed_segmentation <- 
        lidR::watershed(chm = current_chm, th_tree = 1.5, tol = 2, ext = ceiling(2 * (1 / res(current_chm))[1]))() 
      
      crs(watershed_segmentation) <- crs(current_chm) 
      
      watershed_segmentation <-
        watershed_segmentation %>% 
        st_as_stars() %>% 
        st_as_sf(merge = TRUE)

      plot(watershed_segmentation[, "layer"])
      
      
      
      
      
      dalponte2016_segmentation <-
        lapply(ttops_sf, FUN = function(ttops) {
          lidR::dalponte2016(chm = current_chm, treetops = as(ttops, "Spatial"), th_tree = 1.5, th_seed = 0.45, th_cr = 0.55, max_cr = ceiling((10 * 1 / res(current_chm))[1]))() %>% 
            APfun::APpolygonize() %>% 
            st_as_sf() %>% 
            st_set_crs(proj4string(current_chm))
        })
      
      
      silva2016_segmentation <-
        lapply(ttops_sf, FUN = function(ttops) {
          lidR::silva2016(chm = current_chm, treetops = as(ttops, "Spatial"), max_cr_factor = 0.6, exclusion = 0.3)() %>% 
            APfun::APpolygonize() %>% 
            st_as_sf() %>% 
            st_set_crs(proj4string(current_chm))
        })
      
      li2012_segmentation <-
        lapply(ttops_las, FUN = function(las) {
          crowns <- lidR::tree_hulls(las, type = "concave", concavity = 2) %>% st_as_sf()
        })

      
      plot(current_chm)
      plot(current_plot_boundary$geometry, add = TRUE, lwd = 2)
      plot(ttops_sf[[1]]$geometry, add = TRUE, pch = 19)
      plot(li2012_segmentation[[1]]$geometry, add = TRUE, col = "red")
      plot(mcws_segmentation[[1]]$geometry, add = TRUE)
      plot(mcwatershed_segmentation[[1]]$geometry, add = TRUE)
      
      current_ortho <- raster::brick(paste0(current_dir, current_site, "_plot-remote-data/", current_plot, "_ortho.tif"))
      
      source("data/data_carpentry/format_ground-data.R")
      current_plot_ground_trees <- 
        d %>% 
        filter(plot == current_plot) %>% 
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
      
      plotRGB(current_ortho, r = 3, g = 2, b = 1)
      plot(current_plot_ground_trees$geometry[is.na(current_plot_ground_trees$year_fall)], add = TRUE, col = "red", pch = 19)
      plot(current_plot_ground_trees$geometry[!is.na(current_plot_ground_trees$year_fall)], add = TRUE, col = "blue", pch = 19)
      plot(current_plot_boundary$geometry, add = TRUE)
      plot(canopy_segmentation[[1]]$geometry, add = TRUE)
      plot(ttops_sf[[2]]$geometry, add = TRUE, col = "orange", pch = 19)
      
      plot(current_chm)
      plot(current_plot_boundary$geometry, add = TRUE)
      plot(canopy_segmentation[[1]]$geometry, add = TRUE)
      plot(ttops_sf[[1]]$geometry, add = TRUE, col = "orange", pch = 19)
      
            }) # end future_map # end suppressWarnings # end suppressMessages
  

  # First get the tree tops using the 20 best treetop detection algorithms
  # Then try each segmentation algorithm:
  
  # lidR::dalponte2016
  # lidR::watershed # same as EBImage::watershed
  # lidR::mcwatershed # same as imager::watershed with tree tops as a priority map
  # lidR::li2012
  # lidR::silva2016
  # ForestTools::mcws
  
  
  # cat(paste0("...", all_validation_plots[i], " (", i, " of ", length(all_validation_plots), ") complete.\n"))
  # 
  # # Crown Segmentation as a 2nd step
  # # Method 1
  # EBImage watershed
  # crowns_watershed
  # 
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
  
  
  
      
      
      
      
    #   # Combine all of the sf objects from the different treetop detection algorithms
    #   # into a single sf object
    #   # Also add a column that says which algorithm detected the tree
    #   # Further, calculate the nearest neighbor distance between the
    #   # three nearest neighbors for each tree
    #   ttops_sf <-
    #     lapply(seq_along(ttops), FUN = function(j) {
    #       
    #       # Each list element in the ttops list is itself a list of two or three pieces:
    #       # the sf object with the location and height of each tree segmented
    #       # the elapsed time that the algorithm took to complete
    #       # potentially the returned las point cloud if the tree detection algorithm
    #       # returns that as well
    #       
    #       # which of those 2-3 list element in the current ttop list is the sf object?
    #       ttops_sf_idx <- which(names(ttops[[j]]) == "ttops_sf")
    #       
    #       # just select the height and geometry columns of the sf object
    #       ttops[[j]][[ttops_sf_idx]] <-
    #         ttops[[j]][[ttops_sf_idx]] %>%
    #         dplyr::select(height, geometry)
    #       
    #       # If no trees were segmented, create an "empty" sf object so that the rest of
    #       # the objects will bind together nicely since they'll have the same column names, etc.
    #       if (nrow(ttops[[j]][[ttops_sf_idx]]) == 0) {
    #         ttops[[j]][[ttops_sf_idx]] <- 
    #           st_sf(height = NA, geometry = st_sfc(st_multipoint()), crs = st_crs(ttops[[j]][[ttops_sf_idx]]))}
    #       
    #       # Add the treetop detection algorithm name to the data frame row
    #       ttops[[j]][[ttops_sf_idx]] <-
    #         ttops[[j]][[ttops_sf_idx]] %>% 
    #         dplyr::mutate(ttops_method = names(ttops[j]))
    #       
    #       # If there are any trees in the sf object, get the nearest neighbor distance between each tree
    #       # and its three nearest neighbors
    #       if(nrow(ttops[[j]][[ttops_sf_idx]]) > 0) {
    #         nn <- st_nn(ttops[[j]][[ttops_sf_idx]], ttops[[j]][[ttops_sf_idx]], k = min(c(4, nrow(ttops[[j]][[ttops_sf_idx]]))), returnDist = TRUE, sparse = FALSE, progress = FALSE)$dist %>% 
    #           as.data.frame() } else {nn <- data.frame()}
    #       
    #       # remove nuisance column
    #       nn <- nn %>% dplyr::select(-1)
    #       
    #       # build up a nearest neighbors dataframe, but ensure that the data frame will always have 
    #       # three columns (even if there are not enough trees segmented to result in 3
    #       # nearest neighbor distances for each tree). If there aren't enough trees, fill the
    #       # first, first and second, or first, second and third columns of this data frame with NA
    #       if (ncol(nn) == 0) {nn <- nn %>% dplyr::mutate(nn1 = NA, nn2 = NA, nn3 = NA)} else
    #         if (ncol(nn) == 1) {nn <- nn %>% dplyr::rename(nn1 = V2) %>% dplyr::mutate(nn2 = NA, nn3 = NA)} else
    #           if (ncol(nn) == 2) {nn <- nn %>% dplyr::rename(nn1 = V2, nn2 = V3) %>% dplyr::mutate(nn3 = NA)} else
    #             if (ncol(nn) == 3) {nn <- nn %>% dplyr::rename(nn1 = V2, nn2 = V3, nn3 = V4) }
    #       
    #       # bind the nearest neighbor dataframe we just created with the sf object
    #       ttops[[j]][[ttops_sf_idx]] <-
    #         ttops[[j]][[ttops_sf_idx]] %>%
    #         dplyr::bind_cols(nn)
    #       
    #       # transform the crs to epsg4326 to make it compatible across all segmented sf objects
    #       # if the sf object is empty, set its crs to 4326
    #       if (nrow(ttops[[j]][[ttops_sf_idx]]) > 0) {
    #         ttops[[j]][[ttops_sf_idx]] <- 
    #           ttops[[j]][[ttops_sf_idx]] %>% 
    #           sf::st_transform(4326)} else 
    #           {ttops[[j]][[ttops_sf_idx]] <- suppressWarnings(ttops[[j]][[ttops_sf_idx]] %>% sf::st_set_crs(4326))}
    #       
    #     })) %>%
    # do.call("rbind", .) # bind the rows for all the sf objects resulting from different segmentation algorithms together