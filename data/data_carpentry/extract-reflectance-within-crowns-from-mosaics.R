# Extract reflectance values from within the crown polygons

library(sf)
library(raster)
library(tabularaster)
library(tidyverse)
library(viridis)
library(purrr)
library(velox)
library(viridis)

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
  dplyr::filter(overwrite | !classified_check) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()

# Create a list for the classified crowns to live in.
classified_crowns <- vector(mode = "list", length = length(sites_to_process))
i = 1

# # loop through all the sites_to_process
# for (i in seq_along(sites_to_process)) {
#   
#   # `current_site` is the current site in the loop
#   current_site <- sites_to_process[i]
#   
#   # Use the velox package to read in the current site's index mosaic raster
#   # Index mosaic (could be 6-band or 9-band depending on whether it includes the
#   # X3-derived R, G, and B bands)
#   index <- velox::velox(here::here(paste0(current_dir, current_site, "_index.tif")))
#   
index_path <- "data/data_output/site_data/eldo_3k_1/eldo_3k_1_index.tif"
crowns_path <- "data/data_output/classified/hand-classified/eldo_3k_1_hand-classified-crowns/eldo_3k_1_hand-classified-crowns.shp"
ttops_path <- "data/data_output/site_data/eldo_3k_1/eldo_3k_1_ttops/eldo_3k_1_ttops.shp"

extract_reflectance_from_crowns <- function(index_path, crowns_path, ttops_path) {
  index <- velox::velox(here::here(index_path))
  crowns <- sf::st_read(here::here(crowns_path))
  ttops <- sf::st_read(here::here(ttops_path))
  non_spatial_ttops <-
    ttops %>% 
    dplyr::mutate(x = st_coordinates(.)[,1],
                  y = st_coordinates(.)[,2]) %>% 
    st_drop_geometry()
  
  # give the velox object bands for its different raster layers
  if(length(index$rasterbands) == 6) {
    names(index$rasterbands) <- c("b", "g", "r", "re", "nir", "ndvi")
  } else if(length(index$rasterbands) == 6) {
    names(index$rasterbands) <- c("b", "g", "r", "re", "nir", "ndvi", "x3_r", "x3_g", "x3_b")
  }
  
  # create new bands in the velox object using raster algebra on a pixel-by-pixel basis 
  # to get our vegetation indices that might be helpful (RGI and GBI)
  index$rasterbands$rgi <- index$rasterbands$r / index$rasterbands$g
  index$rasterbands$gbi <- index$rasterbands$g / index$rasterbands$b
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
  
  # I copy the velox object, because velox methods alter the original object even without using
  # the assignment operator
  index_copy <- index$copy()
  
  # I learned from here: https://gis.stackexchange.com/questions/286409/fastest-way-to-extract-a-raster-in-r-improve-the-time-of-my-reproducible-code
  # that the velox$extract() method works fastest if you crop the raster to the polygons first.
  index_copy$crop(crowns)
  
  # this is the extract operation; I summarize all pixel values in each crown segment using mean()
  extracted_vals <- index_copy$extract(crowns, fun = function(x) mean(x, na.rm = TRUE), df = TRUE)
  
  # add column names to the resulting object and add a treeID column
  colnames(extracted_vals) <- c("object_", paste0(names(index_copy$rasterbands), "_mean"))
  
  crowns <-
    crowns %>%
    dplyr::left_join(extracted_vals, by = "object_")
  
  return(crowns)
}
