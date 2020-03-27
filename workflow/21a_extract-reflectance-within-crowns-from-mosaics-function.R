# Extract reflectance values from within the crown polygons

library(sf)
library(raster)
library(tidyverse)
library(viridis)
library(purrr)
library(viridis)

# Until velox goes back on CRAN, use the GitHub version
if(!require(velox)) {
  remotes::install_github('hunzikp/velox')
  library(velox)
}

# If the crowns vector shapefile represents a result of hand classifying live/dead and
# species, set the hand_classified= argument to TRUE and the crowns will be subset
# to just ones with some hand classification. This should greatly speed up the process
# and eliminate all the crowns that are just going to be excluded for building a 
# model anyhow.

extract_reflectance_from_crowns <- function(index, crowns, ttops) {

  non_spatial_ttops <-
    ttops %>% 
    dplyr::mutate(x = st_coordinates(.)[,1],
                  y = st_coordinates(.)[,2]) %>% 
    st_drop_geometry()
  
  
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
  # give the velox object bands for its different raster layers
  if(length(index$rasterbands) == 6) {
    names(index$rasterbands) <- c("b", "g", "r", "re", "nir", "ndvi")
  } else if(length(index$rasterbands) == 9) {
    names(index$rasterbands) <- c("b", "g", "r", "re", "nir", "ndvi", "x3_r", "x3_g", "x3_b")
  }
  
  # create new bands in the velox object using raster algebra on a pixel-by-pixel basis 
  # to get our vegetation indices that might be helpful (RGI and GBI)
  index$rasterbands$rgi <- index$rasterbands$r / index$rasterbands$g
  index$rasterbands$cire <- (index$rasterbands$nir / index$rasterbands$re) - 1
  index$rasterbands$cig <- (index$rasterbands$nir / index$rasterbands$g) - 1
  index$rasterbands$ndvi <- (index$rasterbands$nir - index$rasterbands$r) / (index$rasterbands$nir + index$rasterbands$r)
  index$rasterbands$ndre <- (index$rasterbands$re - index$rasterbands$r) / (index$rasterbands$re + index$rasterbands$r)
  
  # When using the velox objects and adding new bands, you have to update the $nbands list elemnent
  # manually in order for other velox methods to work properly on it
  index$nbands <- length(names(index$rasterbands))
  
  # # read in the corresponding crown segment polygon data for the current site and add some extra information
  # # to it
  crowns <-
    crowns %>% 
    dplyr::left_join(non_spatial_ttops) %>% 
    dplyr::mutate(object_ = 1:nrow(.))
  
  # I have since learned that the {exactextractr} package benchmarks even better than the {velox}
  # package for polygon-based cell value extraction from rasters.
  # Learned from here: https://github.com/hunzikp/velox/issues/43#issuecomment-604663649
  # It might be worth refactoring this whole function to use exactextractr instead of velox.
  # https://isciences.gitlab.io/exactextractr/
  
  # I learned from here: https://gis.stackexchange.com/questions/286409/fastest-way-to-extract-a-raster-in-r-improve-the-time-of-my-reproducible-code
  # that the velox$extract() method works fastest if you crop the raster to the polygons first.
  index$crop(crowns)
  
  # this is the extract operation; I summarize all pixel values in each crown segment using mean()
  extracted_vals <- index$extract(crowns, fun = function(x) mean(x, na.rm = TRUE), df = TRUE)
  
  # add column names to the resulting object and add a treeID column
  colnames(extracted_vals) <- c("object_", paste0(names(index$rasterbands), "_mean"))
  
  crowns <-
    crowns %>%
    dplyr::left_join(extracted_vals, by = "object_")
  
  return(crowns)
}

# # Example use:
# index_path <- "data/data_drone/L2/index/eldo_3k_1_index.tif"
# crowns_path <- "data/data_drone/L3b/hand-classified-crowns/eldo_3k_1_hand-classified-crowns.gpkg"
# ttops_path <- "data/data_drone/L3a/ttops/eldo_3k_1_ttops.gpkg"
# 
# index <- velox::velox(here::here(index_path))
# crowns <- sf::st_read(here::here(crowns_path)) %>% dplyr::filter(!is.na(live) | !is.na(species))
# ttops <- sf::st_read(here::here(ttops_path))
# 
# crowns <- extract_reflectance_from_crowns(index = index,
#                                           crowns = crowns,
#                                           ttops = ttops)

