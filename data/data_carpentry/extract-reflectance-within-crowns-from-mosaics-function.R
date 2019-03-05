# Extract reflectance values from within the crown polygons

library(sf)
library(raster)
library(tidyverse)
library(viridis)
library(purrr)
library(velox)
library(viridis)

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
# index_path <- "data/data_output/site_data/eldo_3k_1/eldo_3k_1_index.tif"
# crowns_path <- "data/data_output/classified/hand-classified/eldo_3k_1_hand-classified-crowns/eldo_3k_1_hand-classified-crowns.shp"
# ttops_path <- "data/data_output/site_data/eldo_3k_1/eldo_3k_1_ttops/eldo_3k_1_ttops.shp"
# 
# index <- velox::velox(here::here(index_path))
# crowns <- sf::st_read(here::here(crowns_path)) %>% dplyr::filter(!is.na(live) | !is.na(species))
# ttops <- sf::st_read(here::here(ttops_path))
# 
# crowns <- extract_reflectance_from_crowns(index = index,
#                                           crowns = crowns,
#                                           ttops = ttops)

