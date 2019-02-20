library(sf)
library(tidyverse)
library(lidR)
library(raster)
library(ForestTools)
# devtools::install_github("Jean-Romain/lidRplugins")
library(lidRplugins)

source("analyses/vwf2.R")

# treetop detection algorithm helper functions ----------------------------

st_vwf <- function(CHM, plot_boundary = NULL, winFun, minHeight, maxWinDiameter) {
  start <- Sys.time()
  
  ttops_sf <-
    CHM %>% 
    vwf2(winFun = winFun, minHeight = minHeight, maxWinDiameter = maxWinDiameter) %>%
    st_as_sf() %>%
    st_set_agr("constant") 
  
  if(!is.null(plot_boundary)) {
    ttops_sf <- ttops_sf %>%
      st_intersection(y = plot_boundary)
  }
  
  
  end <- Sys.time()
  
  ttops_time <- end - start
  
  return(list(ttops_sf = ttops_sf, ttops_time = ttops_time))
}


st_lmf <- function(obj, plot_boundary = NULL, ws) {
  start <- Sys.time()
  ttops_sf <-
    obj %>%
    lidR::tree_detection(algorithm = lmf(ws = ws)) %>%
    st_as_sf() %>% 
    rename(height = Z) %>%
    st_set_agr("constant")
  
  if(!is.null(plot_boundary)) {
    ttops_sf <- ttops_sf %>%
      st_intersection(y = plot_boundary)
  }
  
  end <- Sys.time()
  
  ttops_time <- end - start
  
  return(list(ttops_sf = ttops_sf, ttops_time = ttops_time))
}


st_li2012 <- function(las, plot_boundary = NULL, dt1, dt2, R, Zu, hmin, speed_up) {
  start <- Sys.time()
  ttops_las <- 
    las %>% 
    lidR::lastrees(algorithm = li2012(dt1 = dt1, dt2 = dt2, R = R, Zu = Zu, hmin = hmin, speed_up = speed_up))
  
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
  
  ttops_time <- end - start
  
  return(list(ttops_las = ttops_las, ttops_sf = ttops_sf, ttops_time = ttops_time))
}

st_watershed <- function(las, plot_boundary = NULL, chm, th_tree, tol, ext) {
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
    st_set_agr("constant")
  
  if(!is.null(plot_boundary)) {
    ttops_sf <- ttops_sf %>%
      st_intersection(y = plot_boundary)
  }
  
  
  end <- Sys.time()
  
  ttops_time <- end - start
  
  return(list(ttops_las = ttops_las, ttops_sf = ttops_sf, ttops_time = ttops_time))
}

st_ptree <- function(las, plot_boundary = NULL, dtm, k, algorithm_hmin, post_hmin = NULL, nmax = 7, normalize_las = TRUE) {
  
  start <- Sys.time()
  
  ttops_las <-
    las %>%
    lidR::lastrees(algorithm = ptrees(k = k, hmin = algorithm_hmin, nmax = nmax))
  
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
    dplyr::mutate(elev = raster::extract(dtm, .)) %>% 
    dplyr::mutate(height = height - elev)
  
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
  ttops_time <- end - start
  
  return(list(ttops_las = ttops_las, ttops_sf = ttops_sf, ttops_time = ttops_time))
  
}

st_multichm <- function(las, plot_boundary = NULL, res, layer_thickness, dist_2d, dist_3d, use_max, ws) {
  
  start <- Sys.time()
  
  ttops_sf <-
    las %>% 
    lidR::tree_detection(algorithm = multichm(res = res, layer_thickness = layer_thickness, dist_2d = dist_2d, dist_3d = dist_3d, use_max = use_max, ws = ws)) %>% 
    st_as_sf() %>% 
    rename(height = Z) %>%
    st_set_agr("constant")
  
  if(!is.null(plot_boundary)) {
    ttops_sf <- ttops_sf %>%
      st_intersection(y = plot_boundary)
  }
  
  
  end <- Sys.time()
  ttops_time <- end - start
  
  return(list(ttops_sf = ttops_sf, ttops_time = ttops_time))
}


st_lmfx <- function(las, plot_boundary = NULL, hmin, dist_2d, ws) {
  
  start <- Sys.time()
  
  ttops_sf <- 
    las %>% 
    lidR::tree_detection(algorithm = lmfx(hmin = hmin, dist_2d = dist_2d, ws = ws)) %>% 
    st_as_sf() %>% 
    rename(height = Z) %>%
    st_set_agr("constant")
  
  if(!is.null(plot_boundary)) {
    ttops_sf <- ttops_sf %>%
      st_intersection(y = plot_boundary)
  }
  
  
  end <- Sys.time()
  
  ttops_time <- end - start
  
  return(list(ttops_sf = ttops_sf, ttops_time = ttops_time))
  
}