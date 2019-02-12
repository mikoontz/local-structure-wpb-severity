# This is the "variable window filter" function from ForestTools, but with this patch applied to prevent even-sided
# windows in the focal operation: https://github.com/AndyPL22/ForestTools/pull/6

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
