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

# Establish a method for detecting tree tops and segmenting crowns
# Work through all of the remotely-visible plots at each site (144 in total) applying the ttops and crown segmentation approach
# 













# Now there is an R object in the environment called "sites_checklist" that has
# infomation about how far along all processing steps are.






source("data/data_carpentry/make_processing-checklist.R")
source("data/data_carpentry/format_ground-data.R")

plot_locations_files <- list.files("data/data_output/site_data/", pattern = "plot-locations.shp", recursive = TRUE, full.names = TRUE)
# character vector of all study sites

all_sites <- list.files("data/data_output/site_data")

sites_to_process <-
  sites_checklist %>% 
  dplyr::filter(!crowns_check) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()

i = 1
# get the character string representing the ith site
current_site <- all_sites[i]

current_dir <- paste0("data/data_output/site_data/", current_site, "/")

# The Digital Terrain Model (dtm) is the 2m resolution "ground" underneath
# the current site. Created using CloudCompare and the Cloth Simulator Filter
dtm <- raster::raster(paste0(current_dir, current_site, "_2m-dtm.tif"))

# The Digital Surface Model (dsm) is the ~5cm resolution raster representing
# the surface (ground + objects on top) that the drone flew over
dsm <- raster::raster(paste0(current_dir, "3_dsm_ortho/1_dsm/", current_site, "_x3_dsm.tif"))

# Make sure the Coordinate Reference Systems (crs) are the same
raster::crs(dtm) <- raster::crs(dsm)

# Get the orthophoto also
ortho <- raster::brick(paste0("data/data_output/site_data/", current_site, "/3_dsm_ortho/2_mosaic/", current_site, "_x3_transparent_mosaic_group1.tif"))

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

# Cropping the dtm and dsm to the site bounds
site_dtm <- raster::crop(dtm, site_bounds)
site_dsm <- raster::crop(dsm, site_bounds)
site_ortho <- raster::crop(ortho, site_bounds)

# Using bilinear interpolation to downsample the 2m resolution DTM to have the
# same resolution as the dsm (~5cm, but slightly different for each site)
dtm_resamp <- raster::resample(site_dtm, site_dsm, method = "bilinear")

# The Canopy Height Model (chm) is the dsm (vegetation + ground) minus the dtm (ground)
# to give just the height of the vegetation.
chm <- site_dsm - dtm_resamp
chm2 <- chm
chm2[raster::getValues(chm2) < 0] <- 0

# Smooth out the chm and set any negative values to 0 (meaning "ground") following
# advice from Zagalikis, Cameron, and Miller (2004) and references therein
# More recently, a 3x3 pixel smoothing filter was specifically suggested as ideal
# for sUAS derived chm by Mohan et al. (2017)
w <- focalWeight(x = chm, d = 0.5, type = "circle")
w

chm_smooth <- raster::focal(chm, w = w)
chm_smooth[raster::getValues(chm_smooth) < 0] <- 0

plot_radius <- sqrt((66*66) / pi) * 12 *2.54 / 100

# Get plot locations as determined by visually inspecting orthophotos and finding orange X's laid
# across the plot centers
current_site_plot_locations <- 
  sf::st_read(paste0(current_dir, current_site, "_plot-locations/", current_site, "_plot-locations.shp")) %>% 
  mutate(plot = paste(current_site, id, sep = "_")) %>% 
  st_zm()

current_plot_id <- 2 

current_plot_location <-
  current_site_plot_locations %>% 
  filter(plot == paste(current_site, current_plot_id, sep = "_"))

current_plot_boundary <-
  current_plot_location %>% 
  st_buffer(plot_radius)

current_plot_dtm <- raster::crop(dtm_resamp, current_plot_boundary)
current_plot_dsm <- raster::crop(dsm, current_plot_boundary)
current_plot_chm <- raster::crop(chm_smooth, current_plot_boundary)
current_plot_chm2 <- raster::crop(chm2, current_plot_boundary)

current_plot_ortho <- raster::crop(site_ortho, current_plot_boundary)

las_catalog <- lidR::catalog(paste0(current_dir, current_site, "_vegetation-from-csf.las"))
current_plot_las <- lidR::lasclip(las_catalog, as(current_plot_boundary, "Spatial")@polygons[[1]]@Polygons[[1]])
lidR::lasnormalize(las = current_plot_las, dtm = site_dtm)

# current_plot_chm3 <- lidR::grid_canopy(current_plot_las, subcircle = 0.05, na.fill = "knnidw", k = 10, p = 2, res = 0.05) %>% 
#   as.raster()
current_plot_chm3 <- lidR::grid_canopy(current_plot_las, subcircle = 0.05, na.fill = "kriging", k = 10, res = 0.05, model = vgm(10, "Exp", 300, 4.5)) %>% 
  as.raster()

par(mfrow = c(1, 2))
plot(current_plot_chm)
plot(current_plot_boundary$geometry, add = TRUE, lwd = 3)
plot(current_plot_chm3)
dev.off()
current_plot_ground_trees <- 
  d %>% 
  filter(plot == paste(current_site, current_plot_id, sep = "_")) %>% 
  left_join(current_plot_location, by = "plot") %>% 
  mutate(delta_x = cospi((1 / 2) - (azm / 180)) * dist) %>% 
  mutate(delta_y = sinpi((1 / 2) - (azm / 180)) * dist) %>%
  st_as_sf() %>%
  mutate(x = st_coordinates(.)[, "X"],
         y = st_coordinates(.)[, "Y"]) %>% 
  as.data.frame() %>% 
  mutate(new_x = x + delta_x,
         new_y = y + delta_y) %>% 
  filter(!is.na(x)) %>% 
  st_as_sf(coords = c("new_x", "new_y")) %>% 
  st_set_crs(st_crs(current_plot_location))

dynamicWindow <- function(x) {
  # meters_per_side <- x * 0.05 + 0.5
  # Divide the Popescu and Wynne (2004) equations by 2 to convert to *radius* of search window
  # which is what the vwf() function requires
  # Equation and coefficients taken from Popescu and Wynne (2004)'s "Pines model"
  # meters_per_side <- 0.5 * (3.75105 + 0.17919*x + 0.01241 * x^2)
  meters_per_side <- x * 0.025 + 1
  # Equation and coefficients taken from Popescu and Wynne (2004)'s conifer + deciduous model: "Combined model"
  # meters_per_side <- 0.5 * (2.51503 + 0.00901 * x^2)
  
  return(meters_per_side)
}

data.frame(height = 0:60, radius = dynamicWindow(0:60))
# We only measured trees that were greater than 2.5 inches DBH and focused on ponderosa pine
  # which translated to a height of approximately 6 meters (Wonn and O'Hara, 2001) so we ignored
  # trees less than this height
  # Uses the "variable window filter" algorithm by Popescu and Wynne (2004)
  # ttops <- try(ForestTools::vwf(CHM = chm_cropped, winFun = dynamicWindow, minHeight = 4, maxWinDiameter = NULL))
  
  ttops <- 
    current_plot_chm %>% 
    ForestTools::vwf(winFun = dynamicWindow, minHeight = 4, maxWinDiameter = NULL) %>% 
    st_as_sf() %>% 
    st_intersection(current_plot_boundary)
  
  lastrees_li2(current_plot_las, dt1 = 1, dt2 = 1.5, R = 1, Zu = 15, hmin = 5, speed_up = 25)
  
  ttops2 <-
    current_plot_las %>% 
    as.spatial() %>% 
    st_as_sf() %>% 
    group_by(treeID) %>% 
    filter(Z == max(Z)) %>% 
    st_intersection(current_plot_boundary)
  
  ttops3 <-
    tree_detection(current_plot_las, ws = 2, hmin = 4) %>% 
    st_as_sf(coords = c("X", "Y"), remove = FALSE) %>% 
    st_set_crs(st_crs(current_plot_chm)) %>% 
    mutate(treeID = 1:nrow(.))
  
  r <- 
    lastrees_dalponte(las = current_plot_las, 
                    chm = current_plot_chm, 
                    treetops = as.data.frame(ttops3)[, -5], 
                    th_tree = 4, 
                    th_seed = .9, 
                    th_cr = 0.95, 
                    max_cr = 1 / res(current_plot_chm)[1] * 10,
                    extra = TRUE)
  
  crowns <- 
      ForestTools::mcws(treetops = as(ttops, "Spatial"), 
                        CHM = current_plot_chm, 
                        minHeight = 2, 
                        format = "polygons", 
                        OSGeoPath = "C:\\OSGeo4W64") %>% 
        st_as_sf() %>% 
    mutate(id = 1:nrow(.))
  
  crowns2 <- 
    ForestTools::mcws(treetops = as(ttops2, "Spatial"), 
                      CHM = current_plot_chm, 
                      minHeight = 3.5, 
                      format = "polygons", 
                      OSGeoPath = "C:\\OSGeo4W64") %>% 
    st_as_sf() %>% 
    mutate(id = 1:nrow(.))
  
  crowns3 <- 
    ForestTools::mcws(treetops = as(ttops3, "Spatial"), 
                      CHM = current_plot_chm, 
                      minHeight = 3.5, 
                      format = "polygons", 
                      OSGeoPath = "C:\\OSGeo4W64") %>% 
    st_as_sf() %>% 
    mutate(treeID = 1:nrow(.))
  
  crowns4 <- 
    ForestTools::mcws(treetops = as(current_plot_ground_trees, "Spatial"), 
                      CHM = current_plot_chm, 
                      minHeight = 4, 
                      format = "polygons", 
                      OSGeoPath = "C:\\OSGeo4W64") %>% 
    st_as_sf() %>% 
    mutate(id = 1:nrow(.))
  
  s <- 
    current_plot_las %>% 
    as.spatial() %>% 
    st_as_sf()
  
  length(unique(s$treeID))
  plot(s$geometry, col = viridis(length(unique(s$treeID)))[s$treeID])
  # r <- current_plot_chm
  # r <-
  #   raster(nrow = 25, ncol = 25) %>%
  #   setExtent(current_plot_chm)
  # crs(r) <- crs(current_plot_chm)
  # 
  # ttops_r <-
  #   current_plot_las %>%
  #   rasterize(r, field = "treeID")
  # 
  # ttops_poly <-
  #   ttops_r %>% 
  #   APfun::APpolygonize() %>% 
  #   st_as_sf()
  # 
  # raster::plot(ttops_r)
  # plot(ttops_poly$geometry)
  # # raster::plot(current_plot_chm)
  
  fw <- 
    focalWeight(current_plot_chm, d = 1.25, type = "circle")
  fw[fw > 0] <- 1
  
  local_max <-
    current_plot_chm %>% 
    focal(w = fw, fun = max, pad = TRUE, na.rm = TRUE)
  
  local_max[] <- ifelse(local_max[] == current_plot_chm[] & local_max[] != 0, yes = TRUE, no = FALSE)
  plot(local_max)
  
  local_max_xy <- 
    xyFromCell(local_max, which(local_max[])) %>% 
    as.data.frame() %>% 
    st_as_sf(coords = c("x", "y")) %>% 
    st_set_crs(st_crs(current_plot_chm)) %>% 
    mutate(id = 1:nrow(.))

  raster::plotRGB(current_plot_ortho)
  # plot(local_max_xy$geometry, col = "red", pch = 19, cex = 2, add = TRUE)
  # raster::plot(current_plot_chm)
  plot(current_plot_boundary$geometry, add = TRUE, lwd = 3)
  plot(current_plot_ground_trees$geometry, pch = 19, add = TRUE)
  plot(crowns$geometry, col = viridis(nrow(crowns))[crowns$id], add = TRUE)
  
  plot(crowns2$geometry, col = viridis(nrow(crowns2))[crowns2$id], add = TRUE)
  plot(crowns3$geometry, col = viridis(nrow(crowns3))[crowns3$treeID], add = TRUE)
  plot(crowns4$geometry, col = rainbow(nrow(crowns4))[crowns4$id], add = TRUE)
  
  
  plot(ttops$geometry, pch = 19, col = "red", add = TRUE)
  plot(ttops2$geometry, pch = 19, cex = 2, col = "blue", add = TRUE)
  plot(ttops3$geometry, pch = 19, cex = 2, col = "purple", add = TRUE)
  
  plot(current_plot_location$geometry, add = TRUE, col = "red", pch = 19)
  
  text(st_coordinates(current_plot_ground_trees)[,"X"] + 1, st_coordinates(current_plot_ground_trees)[, "Y"], labels = round(current_plot_ground_trees$height, 2))
  text(st_coordinates(ttops)[,"X"] + 1, st_coordinates(ttops)[, "Y"], labels = round(ttops$height, 2), col = "red")

  