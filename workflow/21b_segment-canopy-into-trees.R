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
library(tictoc)
library(here)
library(units)

source(here::here("data/data_carpentry/segmentation-helper-functions.R"))
source(here::here("data/data_carpentry/make-processing-checklist.R"))

sites_checklist
# These sites had X3 and RedEdge photos merged into the same project, so we look in a different place for some of the relevant
# files.
merged_sites <- c("eldo_3k_2",
                  "eldo_3k_3",
                  "eldo_4k_2")


unusable_sites <- c("eldo_4k_3", # too many blocks
                    "stan_4k_3", # too many blocks
                    "stan_5k_3", # too many blocks
                    "sequ_4k_2") # middle section flown on a separate day and the stitch looks terrible

# The merged versus the unmerged sites will have different numbers of bands
# Because the merged sites will have the X3 imagery incorporated, there will
# be an extra 3 bands for the ortho and the index outputs (the R, G, and B
# from the X3 camera)
# The R, G, and B bands from the X3 images will always be the final three bands
# if they exist.
# For the RedEdge-derived products, the bands go in the order of wavelength,
# from shortest to longest (B, G, R, RE, NIR)
# There is one Pix4D derived index (NDVI), which will go after the NIR band
# for the index mosaic

# This is where I can put in sites that need their processing redone. An empty 
# string means that no already-processed site output will be overwritten
# (but sites that have yet to be processed will still have their processing done)
sites_to_overwrite <- "all"
sites_checklist$overwrite <- ifelse(sites_to_overwrite == "all", yes = TRUE, no = FALSE)

sites_checklist[sites_checklist$site %in% sites_to_overwrite, "overwrite"] <- TRUE

sites_to_process <- 
  sites_checklist %>% 
  dplyr::filter(!(site %in% unusable_sites)) %>%
  dplyr::filter(overwrite | !ttops_check | !crowns_check) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()

# Get an example for testing
# current_site <- sites_to_process[1]

sites_to_process %>% 
  walk(.f = function(current_site) {
    current_dir <- paste0("data/data_output/site_data/", current_site, "/")
    
    if (dir.exists(here::here(paste0(current_dir, current_site, "_ttops")))) {
      unlink(here::here(paste0(current_dir, current_site, "_ttops")), recursive = TRUE)
    }
    
    if (dir.exists(here::here(paste0(current_dir, current_site, "_crowns")))) {
      unlink(here::here(paste0(current_dir, current_site, "_crowns")), recursive = TRUE)
    }
  })

tic()
# Set up the parallelization
num_cores_to_use <- availableCores() - 6 # Took 1996.05 sec on Alienware using 6 workers + 64GB RAM
plan(multiprocess, workers = num_cores_to_use)

crowns <-
  sites_to_process %>% 
  furrr::future_map(.f = function(current_site) {
 
    current_dir <- paste0("data/data_output/site_data/", current_site, "/")
 
    # The dtm is the terrain model for a particular site. We need it to normalize the objects returned from the
    # ptrees algorithm, which acts on a non-normalized point cloud.
    current_dtm <- raster::raster(here::here(paste0(current_dir, current_site, "_dtm.tif")))
    
    # The Canopy Height Model (chm) is the dsm (vegetation + ground) minus the dtm (ground)
    # to give just the height of the vegetation.
    current_chm_rough <- raster::raster(here::here(paste0(current_dir, current_site, "_chm.tif")))
    
    # The point cloud is directly used for some segmentation algorithms, so we import that too
    current_las <- lidR::readLAS(paste0(current_dir, current_site, "_classified-point-cloud.las"))
    current_las_normalized <- lidR::lasnormalize(las = current_las, algorithm = current_dtm, na.rm = TRUE)
    
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
    
    dist2d <- 1.5
    ws <- 2.5
    
    max_crown_area <- pi * 10^2 %>% set_units(m^2)
    
    ttops <-
      current_las_normalized %>%
      st_lmfx(hmin = min_height, dist_2d = dist2d, ws = ws)
    
    ttops_sf <-
      ttops$ttops_sf 

    # Here, we may want to subset the tree tops to buffer inside from the edge of
    # the outermost treetops that were detected. We would do the same thing for
    # the crown segments.
    # survey_area <- ttops_sf %>% st_combine() %>% st_convex_hull() %>% st_buffer(-25)
    # ttops_sf <- ttops_sf %>% st_intersection(survey_area)
    # 
    max_crown_area <- pi * 10^2 %>% set_units(m^2)
    
    # Write the ttops point geometries to a file
     dir.create(here::here(paste0(current_dir, current_site, "_ttops")))
    
    st_write(obj = ttops_sf, dsn = here::here(paste0(current_dir, current_site, "_ttops/", current_site, "_ttops.shp")), delete_dsn = TRUE)
    

    non_spatial_ttops <-
      ttops_sf %>%
      dplyr::mutate(x = st_coordinates(.)[, 1],
                    y = st_coordinates(.)[, 2]) %>%
      sf::st_drop_geometry()
  
    crowns <-
      ttops_sf %>%
      as("Spatial") %>%
      ForestTools::mcws(CHM = current_chm, minHeight = 1.5, format = "raster") %>%
      setNames(nm = "treeID") %>%
      st_as_stars() %>%
      st_as_sf(merge = TRUE) %>%
      st_join(ttops_sf) %>% 
      dplyr::filter(!is.na(treeID.y)) %>% 
      dplyr::rename(treeID = treeID.x) %>% 
      dplyr::select(-treeID.y) %>% 
      dplyr::mutate(ch_area = st_area(st_convex_hull(.))) %>% 
      dplyr::filter(ch_area < max_crown_area)
    
    points_to_crowns <-
      non_spatial_ttops %>%
      dplyr::anti_join(crowns, by = "treeID") %>%
      st_as_sf(coords = c("x", "y"), crs = st_crs(ttops_sf)) %>%
      st_buffer(dist = 0.5) %>% 
      dplyr::mutate(ch_area = st_area(st_convex_hull(.)))

    crowns <-
      crowns %>%
      rbind(points_to_crowns)

    # Write the crowns polygons to a file
    dir.create(here::here(paste0(current_dir, current_site, "_crowns")))
    
    st_write(obj = crowns, dsn = here::here(paste0(current_dir, current_site, "_crowns/", current_site, "_crowns.shp")), delete_dsn = TRUE)
    
    return(list(ttops = ttops_sf, crowns = crowns))
    
  }) # end unnamed, mapped function; end map
toc()