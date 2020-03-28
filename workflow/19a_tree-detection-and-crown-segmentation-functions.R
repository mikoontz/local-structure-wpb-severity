# This set of functions helps me to produce comparable outputs from various
# tree detection algorithms so that I can run them all in a different script
# and know exactly what the returned values will be.
# I also include some extra functionality like including the runtime for each algorithm

library(sf)
library(tidyverse)
library(raster)
library(lidR)

# We want to use the ForestTools version that has our patch applied in order to avoid
# even dimensions of the focal matrix
# It is not yet on CRAN (and no telling when it will be), so we install from GitHub
# See PR here: https://github.com/andrew-plowright/ForestTools/pull/6
remotes::install_github("andrew-plowright/ForestTools@01a789ee37b10430eff4fb1df841d5298aa91f25")

library(ForestTools)

if (class(try(lidRplugins::lmfx)) != "function" | try(packageVersion("lidRplugins")) != "0.2.0") {
  remotes::install_github("mikoontz/lidRplugins@master")
}

library(lidRplugins)

# This is the "variable window filter" function from the {ForestTools} package (which is an awesome package), 
# but with this patch applied to prevent even-sided windows in the focal operation: 
# https://github.com/AndyPL22/ForestTools/pull/6
# As of 2020-03-24, this patch has been merged into the AndyPL22/ForestTools repo, but hasn't been 
# submitted to CRAN so I chose to recreate the relevant function here rather than using remotes::install_github()
# on the unsubmitted version.
# To be clear: most of this function was written by Andrew Plowright (thank you!) except for the small
# patch. For more clarity about which part I wrote, see the link above.

# treetop detection algorithm helper functions ----------------------------

st_vwf <- function(CHM, plot_boundary = NULL, winFun, minHeight, maxWinDiameter) {
  start <- Sys.time()
  
  ttops_sf <-
    CHM %>% 
    ForestTools::vwf(winFun = winFun, minHeight = minHeight, maxWinDiameter = maxWinDiameter) %>%
    st_as_sf() %>%
    st_set_agr("constant") 
  
  if(!is.null(plot_boundary)) {
    ttops_sf <- ttops_sf %>%
      st_intersection(y = plot_boundary)
  }
  
  
  end <- Sys.time()
  
  ttops_time <- as.numeric(base::difftime(time1 = end, time2 = start, units = "secs"))
  
  return(list(ttops_sf = ttops_sf, ttops_time = ttops_time))
}


st_lmf <- function(obj, plot_boundary = NULL, ws) {
  start <- Sys.time()
  ttops_sf <-
    obj %>%
    lidR::find_trees(algorithm = lmf(ws = ws)) %>%
    st_as_sf() %>% 
    rename(height = Z) %>%
    st_set_agr("constant")
  
  if(!is.null(plot_boundary)) {
    ttops_sf <- ttops_sf %>%
      st_intersection(y = plot_boundary)
  }
  
  end <- Sys.time()
  
  ttops_time <- as.numeric(base::difftime(time1 = end, time2 = start, units = "secs"))
  
  return(list(ttops_sf = ttops_sf, ttops_time = ttops_time))
}


st_li2012 <- function(las, plot_boundary = NULL, dt1, dt2, R, Zu, hmin, speed_up) {
  start <- Sys.time()
  ttops_las <- 
    las %>% 
    lidR::segment_trees(algorithm = li2012(dt1 = dt1, dt2 = dt2, R = R, Zu = Zu, hmin = hmin, speed_up = speed_up))
  
  ttops_sf <-
    ttops_las %>% 
    slot("data") %>% 
    dplyr::group_by(treeID) %>%
    dplyr::filter(Z == max(Z)) %>%
    dplyr::ungroup() %>%
    st_as_sf(coords = c("X", "Y"),
             crs = proj4string(las)) %>%
    rename(height = Z) %>%
    st_set_agr("constant") 
  
  if(!is.null(plot_boundary)) {
    ttops_sf <- ttops_sf %>%
      st_intersection(y = plot_boundary)
  }
  
  end <- Sys.time()
  
  ttops_time <- as.numeric(base::difftime(time1 = end, time2 = start, units = "secs"))
  
  return(list(ttops_las = ttops_las, ttops_sf = ttops_sf, ttops_time = ttops_time))
}

st_watershed <- function(las, plot_boundary = NULL, chm, th_tree, tol, ext) {
  start <- Sys.time()
  
  ttops_las <-
    las %>% 
    lidR::segment_trees(algorithm = watershed(chm = current_chm, th_tree = th_tree, tol = tol, ext = ext))
  
  ttops_sf <-  
    ttops_las %>% 
    slot("data") %>% 
    dplyr::group_by(treeID) %>%
    dplyr::filter(Z == max(Z)) %>%
    dplyr::ungroup() %>%
    st_as_sf(coords = c("X", "Y"),
             crs = proj4string(las)) %>%
    rename(height = Z) %>%
    st_set_agr("constant")
  
  if(!is.null(plot_boundary)) {
    ttops_sf <- ttops_sf %>%
      st_intersection(y = plot_boundary)
  }
  
  
  end <- Sys.time()
  
  ttops_time <- as.numeric(base::difftime(time1 = end, time2 = start, units = "secs"))
  
  return(list(ttops_las = ttops_las, ttops_sf = ttops_sf, ttops_time = ttops_time))
}

st_ptree <- function(las, plot_boundary = NULL, dtm, k, algorithm_hmin, post_hmin = NULL, nmax = 7, normalize_las = TRUE) {
  
  start <- Sys.time()
  
  ttops_las <-
    las %>%
    lidR::segment_trees(algorithm = ptrees(k = k, hmin = algorithm_hmin, nmax = nmax))
  
  ttops_sf <-
    ttops_las %>%
    slot("data") %>% 
    dplyr::group_by(treeID) %>%
    dplyr::filter(Z == max(Z)) %>%
    dplyr::ungroup() %>%
    st_as_sf(coords = c("X", "Y"),
             crs = proj4string(las)) %>%
    rename(height = Z) %>%
    st_set_agr("constant")
  
  if(!is.null(plot_boundary)) {
    ttops_sf <- ttops_sf %>%
      st_intersection(y = plot_boundary)
  }
  
  ttops_sf <-
    ttops_sf %>%
    dplyr::mutate(elev = ifelse(nrow(.) > 0, yes = raster::extract(dtm, .), no = NA_real_)) %>% 
    dplyr::mutate(height = ifelse(nrow(.) > 0, yes = height - elev, no = NA_real_))
  
  if(!is.null(post_hmin)) {
    ttops_sf <-
      ttops_sf %>% 
      dplyr::filter(height > post_hmin)
  }
  
  if(normalize_las) {
    ttops_las <- 
      ttops_las %>% 
      lidR::lasnormalize(algorithm = tin())
  }
  
  end <- Sys.time()
  ttops_time <- as.numeric(base::difftime(time1 = end, time2 = start, units = "secs"))
  
  return(list(ttops_las = ttops_las, ttops_sf = ttops_sf, ttops_time = ttops_time))
  
}

st_multichm <- function(las, plot_boundary = NULL, res, layer_thickness, dist_2d, dist_3d, use_max, ws) {
  
  start <- Sys.time()
  
  ttops_sf <-
    las %>% 
    lidR::find_trees(algorithm = multichm(res = res, layer_thickness = layer_thickness, dist_2d = dist_2d, dist_3d = dist_3d, use_max = use_max, ws = ws)) %>% 
    st_as_sf() %>% 
    rename(height = Z) %>%
    st_set_agr("constant")
  
  if(!is.null(plot_boundary)) {
    ttops_sf <- ttops_sf %>%
      st_intersection(y = plot_boundary)
  }
  
  
  end <- Sys.time()
  ttops_time <- as.numeric(base::difftime(time1 = end, time2 = start, units = "secs"))
  
  return(list(ttops_sf = ttops_sf, ttops_time = ttops_time))
}


st_lmfx <- function(las, plot_boundary = NULL, hmin, dist_2d, ws) {
  
  start <- Sys.time()
  
  ttops_sf <- 
    las %>% 
    lidR::find_trees(algorithm = lmfx(hmin = hmin, dist_2d = dist_2d, ws = ws)) %>% 
    st_as_sf() %>% 
    rename(height = Z) %>%
    st_set_agr("constant")
  
  if(!is.null(plot_boundary)) {
    ttops_sf <- ttops_sf %>%
      st_intersection(y = plot_boundary)
  }
  
  
  end <- Sys.time()
  
  ttops_time <- as.numeric(base::difftime(time1 = end, time2 = start, units = "secs"))
  
  return(list(ttops_sf = ttops_sf, ttops_time = ttops_time))
  
}