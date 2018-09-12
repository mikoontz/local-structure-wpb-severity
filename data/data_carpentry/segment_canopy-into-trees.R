# Load parallelization packages
library(foreach)
library(doParallel)
library(tidyverse)

# character vector of all study sites
all_sites <- list.files("data/data_output")

sites_checklist <- 
  expand.grid(
    c("eldo", "stan", "sier", "sequ"), 
    c("3k", "4k", "5k", "6k"), 
    1:3) %>% 
  as.data.frame() %>% 
  setNames(c("forest", "elev", "rep")) %>% 
  dplyr::filter(!(forest == "sequ" & elev == "3k")) %>%
  dplyr::filter(!(forest != "sequ" & elev == "6k")) %>% 
  dplyr::arrange(forest, elev, rep) %>% 
  dplyr::mutate(site = paste(forest, elev, rep, sep = "_")) %>% 
  dplyr::select(-forest, -elev, -rep) %>% 
  dplyr::mutate(ttops_check = file.exists(paste0("data/data_output/", site, "/", site, "_ttops/", site, "_ttops.shp"))) %>% 
  dplyr::mutate(crowns_check = file.exists(paste0("data/data_output/", site, "/", site, "_crowns/", site, "_crowns.shp")))

sites_to_process <-
  sites_checklist %>% 
  dplyr::filter(!ttops_check & site %in% all_sites) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()

# I revisit this StackOverflow question and answer, like, everytime I want
# to do anything in parallel to remind myself how
# https://stackoverflow.com/questions/38318139/run-a-for-loop-in-parallel-in-r

cores <- detectCores() 
cl <- makeCluster(min(c(cores - 2), 4)) #not to overload your computer
registerDoParallel(cl)

print(Sys.time())
foreach (i = seq_along(sites_to_process)) %dopar% {

  # Load all packages within the foreach() construct for parallelization to work
  library(raster)
  library(sf)
  library(tidyverse)
  library(ForestTools)
  
  # get the character string representing the ith site
  current_site <- sites_to_process[i]
  
  # convenience string to set file paths for input/output
  current_dir <- paste0("data/data_output/", current_site, "/")
  
  # The Digital Terrain Model (dtm) is the 2m resolution "ground" underneath
  # the current site. Created using CloudCompare and the Cloth Simulator Filter
  dtm <- raster::raster(paste0(current_dir, current_site, "_2m-dtm.tif"))
  
  # The Digital Surface Model (dsm) is the ~5cm resolution raster representing
  # the surface (ground + objects on top) that the drone flew over
  dsm <- raster::raster(paste0(current_dir, "3_dsm_ortho/1_dsm/", current_site, "_x3_dsm.tif"))
  
  # Make sure the Coordinate Reference Systems (crs) are the same
  raster::crs(dtm) <- raster::crs(dsm)

  # The flight bounds is the polygon surrounding the images that the drone took. Determined
  # in the "convert_flight-path-to-site-boundary.R script using the elevation model used
  # during flight planning, the takeoff altitude, the ground altitude underneath each photo,
  # and the altitude difference between each photo and the takeoff point. Essentially, we 
  # crop all the products (e.g., the dsm, the dtm, the lidar point cloud) to just be 
  # area within the flight path, rather than some of the spillover information beyond the
  # area that the drone directly flew over.
  # Could consider buffering this further (inward-- so using a negative buffer value) to 
  # reduce edge effects
  site_bounds <- 
    sf::st_read(paste0(current_dir, current_site, "_mission-footprint/", current_site, "_site-bounds.geoJSON")) %>% 
    sf::st_transform(sp::proj4string(dtm))
  
  # Cropping the dtm and dsm to the fligth bounds
  site_dtm <- raster::crop(dtm, site_bounds)
  site_dsm <- raster::crop(dsm, site_bounds)
  
  # Using bilinear interpolation to downsample the 2m resolution DTM to have the
  # same resolution as the dsm (~5cm, but slightly different for each site)
  dtm_resamp <- raster::resample(site_dtm, site_dsm, method = "bilinear")
  
  # The Canopy Height Model (chm) is the dsm (vegetation + ground) minus the dtm (ground)
  # to give just the height of the vegetation.
  chm <- site_dsm - dtm_resamp
  
  # Smooth out the chm and set any negative values to 0 (meaning "ground") following
  # advice from Zagalikis, Cameron, and Miller (2004) and references therein
  # More recently, a 3x3 pixel smoothing filter was specifically suggested as ideal
  # for sUAS derived chm by Mohan et al. (2017)
  chm_smooth <- raster::focal(chm, w = matrix(1, 3, 3), mean)
  chm_smooth[raster::getValues(chm_smooth) < 0] <- 0

  # Using ForestTools package written by Andrew Plowright and the chm
  # (https://cran.r-project.org/web/packages/ForestTools/vignettes/treetopAnalysis.html)
  # define a function that sets the variable window size for finding maxima
  # This linear relationship is the one used in the vignette, which seems pretty appropriate for our
  # forests too. Our maximum tree sizes are about 60 meters, which would mean the window size to search
  # for a maximum height is 3.6 meters x 3.6 meters
  # This is a point where the workflow could be dialed in to improve accuracy given the on-the-ground
  # data we have.
  
  dynamicWindow <- function(x) {
    meters_per_side <- x * 0.05 + 0.6
    return(meters_per_side)
  }
  
  # We only measured trees that were greater than 2.5 inches DBH and focused on ponderosa pine
  # which translated to a height of approximately 6 meters (Wonn and O'Hara, 2001) so we ignored
  # trees less than this height
  # Uses the "variable window filter" algorithm by Popescu and Wynne (2004)
  ttops <- try(ForestTools::vwf(CHM = chm_smooth, winFun = dynamicWindow, minHeight = 6, maxWinDiameter = NULL))
  
  if (class(ttops) != "try-error") {
    
    # Convert ttops sp object to an sf object for easier writing
    sf_ttops <- 
      ttops %>% 
      st_as_sf()
    
    # output the tree tops to a file; create a new directory to puth the .shp
    # related files in it. This will preserve the CRS that the points came from
    # originally
    
    if (!dir.exists(paste0(current_dir, current_site, "_ttops"))) {
      dir.create(paste0(current_dir, current_site, "_ttops"))
    }
    
    st_write(obj = sf_ttops, dsn = paste0(current_dir, current_site, "_ttops/", current_site, "_ttops.shp"))
    
  }
  
  # Do the segmentation constrained by the "known" tree top locations
  # Called "Marker-controlled watershed segmentation" following Beucher & Meyer, 1993
  
  # Note this function masks all raster points lower than the "minHeight" variable, so it
  # is a different kind of parameter than the minHeight parameter from the vwf() function
  # Function returns a SpatialPolygonsDataFrame with tree attributes included
  
  # Uses the python version (in OSGeo4w64) of polygonize because the R version is bad.
  # Andrew Plowright is brilliant. I spent days trying to get the gdal_polygonize.py script
  # to work or trying to get the gdal.polygonize() method to work directly in python and
  # could not make it happen for the life of me (tried command line, OSGeo4w64 shell, reticulate
  # in R... nothing; finally got it to work in QGIS, but that would have been an unfortunate 
  # pointing and clicking step...). So thank you, Andrew Plowright, because your implementation
  # using the APpolygonize() function in the APfun package you wrote just *works*
  
  if (class(ttops) != "try-error") {
    
    crowns <- 
      try({
        ForestTools::mcws(treetops = ttops, 
                          CHM = chm_smooth, 
                          minHeight = 3.5, 
                          format = "polygons", 
                          OSGeoPath = "C:\\OSGeo4W64") %>% 
          st_as_sf()
      })
    
    if (all(class(crowns) != "try-error")) {
      
      # Write the crowns polygons to a file
      if (!dir.exists(paste0(current_dir, current_site, "_crowns"))) {
        dir.create(paste0(current_dir, current_site, "_crowns"))
      }
      st_write(obj = crowns, dsn = paste0(current_dir, current_site, "_crowns/", current_site, "_crowns.shp"))
    
  }
  paste0("...", current_site, " completed.")
}
}
print(Sys.time())
stopCluster(cl)
