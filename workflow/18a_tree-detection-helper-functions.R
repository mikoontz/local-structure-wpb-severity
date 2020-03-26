# This set of functions helps me to produce comparable outputs from various
# tree detection algorithms so that I can run them all in a different script
# and know exactly what the returned values will be.
# I also include some extra functionality like including the runtime for each algorithm

library(sf)
library(tidyverse)
library(raster)
library(ForestTools)
library(lidR)

if (try(packageVersion("lidRplugins")) != "0.2.0") {
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

vwf2 <- function(CHM, winFun, minHeight = NULL, maxWinDiameter = 99, minWinNeib = "queen", verbose = FALSE){
  
  ### CHECK INPUTS
  
  if(verbose) cat("Checking inputs", "\n")
  
  # Check for valid inputs for 'minWinNeib'
  if(!minWinNeib %in% c("queen", "rook")) stop("Invalid input for 'minWinNeib'. Set to 'queen' or 'rook'")
  
  # Check for unprojected rasters
  CHM.crs <- as.character(raster::crs(CHM))
  CHM.prj <- regmatches(CHM.crs, regexpr("(?<=proj=).*?(?=\\s)", CHM.crs, perl = TRUE))
  if(length(CHM.prj) > 0 && CHM.prj %in% c(c("latlong", "latlon", "longlat", "lonlat"))){
    warning("'CHM' map units are in degrees. Ensure that 'winFun' outputs values in this unit.")
  }
  
  # Round out CHM resolution to fifth decimal and check that CHM has square cells.
  # Rounding is necessary since a lack of precision in CHM cell size call cause the
  # 'focalWeight' function to misbehave
  roundRes <- round(raster::res(CHM), 5)
  if(roundRes[1] != roundRes[2]) stop("Input 'CHM' does not have square cells")
  if(roundRes[1] == 0)           stop("The map units of the 'CHM' are too small")
  
  # Ensure that 'minHeight' argument is given a positive value
  if(!is.null(minHeight) && minHeight <= 0) stop("Minimum canopy height must be set to a positive value.")
  
  # Get range of CHM values
  CHM.rng <- suppressWarnings(raster::cellStats(CHM, range))
  names(CHM.rng) <- c("min", "max")
  
  # Check if CHM has usable values
  if(is.infinite(CHM.rng["max"]) | is.infinite(CHM.rng["min"])){stop("Input 'CHM' does not contain any usable values. Check input data.")}
  
  ### APPLY MINIMUM CANOPY HEIGHT ----
  
  if(!is.null(minHeight)){
    
    if(minHeight >= CHM.rng["max"]) stop("'minHeight' is set to a value higher than the highest cell value in 'CHM'")
    
    # Mask sections of CHM that are lower than 'minHeight'
    if(minHeight > CHM.rng["min"]){
      
      CHM[CHM < minHeight] <- NA
      CHM.rng["min"] <- minHeight
      
    }
  }
  
  ### CREATE WINDOWS ----
  
  if(verbose) cat("Creating windows", "\n")
  
  # Generate a list of radii
  seqFloor   <- floor(  winFun(CHM.rng["min"]))
  if(is.infinite(seqFloor)) seqFloor <- 0 # Watch out for parabola!
  seqCeiling <- ceiling(winFun(CHM.rng["max"]))
  radii <- seq(seqFloor, seqCeiling, by = roundRes[1])
  
  # Remove radii that are smaller than the CHM's resolution
  radii <- radii[radii >= roundRes[1]]
  if(length(radii) == 0){
    warning("The maximum window radius computed with 'winFun' is smaller than the CHM's resolution",
            "\nA 3x3 cell search window will be uniformly applied",
            "\nUse a higher resolution 'CHM' or adjust 'winFun' to produce wider dynamic windows")
    radii <- roundRes[1]
  }
  
  # Calculate the dimensions of the largest matrix to be created from the generated list of radii
  maxDimension <- ceiling((max(radii) / roundRes[1]) * 2 + 1)
  if (maxDimension %% 2 == 0) {maxDimension <- maxDimension + 1}
  
  # Check if input formula will yield a window size bigger than the maximum set by 'maxWinDiameter'
  if(!is.null(maxWinDiameter) && maxDimension > maxWinDiameter){
    
    stop("Input function for 'winFun' yields a window diameter of ",  maxDimension, " cells, which is wider than the maximum allowable value in \'maxWinDiameter\'.",
         "\nChange the 'winFun' function or set 'maxWinDiameter' to a higher value (or to NULL).")
  }
  
  # Convert radii into windows
  windows <- lapply(radii, function(radius){
    
    # Based on the unit size of the input CHM and a given radius, this function will create a matrix whose non-zero
    # values will form the shape of a circular window
    win.mat <- raster::focalWeight(raster::raster(resolution = roundRes), radius, type = "circle")
    
    # Apply Queen's neighborhood if circle is 3x3
    if(nrow(win.mat) == 3 && minWinNeib == "queen") win.mat[] <- 1
    
    # Pad the window to the size of the biggest matrix created from the list of radii
    win.pad <- raster::extend(raster::raster(win.mat), (maxDimension - ncol(win.mat)) /2, value = 0)
    
    # The matrix values are then transformed into a vector
    win.vec <- as.vector(win.pad != 0)
    
    return(win.vec)
    
  })
  names(windows) <- radii
  
  ### CREATE VWF FUNCTION ----
  
  .vwMax <- function(x, ...){
    
    # Locate central value in the moving window.
    centralValue <- x[length(x) / 2 + 0.5]
    
    # If central value is NA, then return NA.
    if(is.na(centralValue)){
      
      return(NA)
      
    }else{
      
      # Calculate the expected crown radius.
      radius <- winFun(centralValue)
      
      # Retrieve windows size closest to radius
      window <- windows[[which.min(abs(as.numeric(names(windows)) -  radius))]]
      
      # If the central value is the highest value within the variably-sized window (i.e.: local maxima), return 1. If not, return 0.
      return(if(max(x[window], na.rm = TRUE) == centralValue) 1 else 0)
      
    }
  }
  
  
  ### APPLY VWF FUNCTION ----
  
  # Apply local maxima-finding function to raster
  lm.pts <- raster::rasterToPoints(
    raster::focal(CHM,
                  w = matrix(1, maxDimension, maxDimension),
                  fun = .vwMax,
                  pad = TRUE, padValue = NA),
    fun = function(x) x == 1,
    spatial = TRUE)
  
  lm.pts[["height"]]    <- raster::extract(CHM, lm.pts)
  lm.pts[["winRadius"]] <- winFun(lm.pts[["height"]])
  lm.pts[["treeID"]]    <- 1:length(lm.pts)
  
  lm.pts[["layer"]] <- NULL
  
  ### RETURN OUTPUT ----
  
  return(lm.pts)
  
}


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