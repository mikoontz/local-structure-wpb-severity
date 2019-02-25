library(sf)
library(raster)
library(tabularaster)
library(tidyverse)
library(viridis)
library(purrr)
library(velox)
library(viridis)


# character vector of all study sites
all_sites <- list.files("data/data_output/site_data")

overwrite <- FALSE

### This first block of code reads in all the hand classified crown shapefiles, extracts all pixel values from
### the orthomosaics which the crown polygons overlap, and summarizes those values into a single value representing
### each tree's reflectance signature (red band, green band, blue band, red / green, green / blue)

### We use the raster and tabularaster package for this operation. It turns out this approach is pretty slow, and
### so we use a different approach for the second block of code to extract reflectance data for *all* tree crowns
### at *all* sites. We could refactor this code to make it match, if we'd like.

if (!file.exists("data/data_output/classified/hand-classified/hand-classified-and-extracted.rda") | overwrite) {
  
  # the list of site names that we spent some time hand classifying
  hand_classified_site_names <- substr(x = list.files("data/data_output/classified/hand-classified/crown-shapefiles"), start = 1, stop = 9)
  
  # ttop data for those sites that we hand classified
  hand_classified_ttops <- 
    hand_classified_site_names %>% 
    map(.f = function(x) {
      sf::st_read(paste0("data/data_output/site_data/", x, "/", x, "_ttops/", x, "_ttops.shp"))
    })
  
  # read the crown segment polygon data for the hand classified sites, add some information from the ttop data
  hand_classified_crowns <-
    hand_classified_site_names %>% 
    map2(.y = hand_classified_ttops, .f = function(x, y) {
      crowns <- 
        sf::st_read(paste0("data/data_output/classified/hand-classified/crown-shapefiles/", x, "_classified-crowns/", x, "_classified-crowns.shp")) %>% 
        dplyr::mutate(site = x) %>% 
        dplyr::mutate(crs = st_crs(.)$proj4string) %>%
        dplyr::mutate(treeID = y$treeID) %>% 
        dplyr::mutate(ttop_x = st_coordinates(y)[, "X"]) %>% 
        dplyr::mutate(ttop_y = st_coordinates(y)[, "Y"]) %>% 
        dplyr::filter(!is.na(live)) %>% 
        dplyr::mutate(object_ = 1:nrow(.))
    })

  # read the orthomosaic data corresponding to the sites that we hand classified
  hand_classified_orthos <-
    hand_classified_site_names %>% 
    map(.f = function(x) {
      raster::brick(paste0("data/data_output/site_data/", x, "/3_dsm_ortho/2_mosaic/", x, "_x3_transparent_mosaic_group1.tif"))
    })
  
  # extract pixel values for each crown segment that was classified, do some algebra on them to get additional reflectance
  # values that might be helpful (RGI and GBI), get the mean for each spectral band across all pixels in each crown segment
  # We also got the 75th percentile data, but didn't end up using it.
  
  hand_classified_crowns <- 
    hand_classified_crowns %>% 
    map2(.y = hand_classified_orthos, .f = function(crowns, ortho) {
      live_or_dead <- 
        cellnumbers(ortho[[1]], crowns) %>% 
        dplyr::mutate(r = raster::extract(ortho[[1]], cell_)) %>% 
        dplyr::mutate(g = raster::extract(ortho[[2]], cell_)) %>% 
        dplyr::mutate(b = raster::extract(ortho[[3]], cell_)) %>% 
        dplyr::mutate(rgi = r/g) %>% 
        dplyr::mutate(rgvi = (r - g) / (r + g)) %>% 
        dplyr::mutate(gbi = g/b) %>% 
        dplyr::mutate(gbvi = (g - b) / (g + b)) %>% 
        dplyr::group_by(object_) %>% 
        dplyr::summarize(r_mean = mean(r, na.rm = TRUE),
                         g_mean = mean(g, na.rm = TRUE),
                         b_mean = mean(b, na.rm = TRUE),
                         rgi_mean = mean(rgi, na.rm = TRUE),
                         rgvi_mean = mean(rgvi, na.rm = TRUE),
                         gbi_mean = mean(gbi, na.rm = TRUE),
                         gbvi_mean = mean(gbvi, na.rm = TRUE),
                         r_hi = quantile(r, probs = 0.75, na.rm = TRUE),
                         g_hi = quantile(g, probs = 0.75, na.rm = TRUE),
                         b_hi = quantile(b, probs = 0.75, na.rm = TRUE),
                         rgi_hi = quantile(rgi, probs = 0.75, na.rm = TRUE),
                         rgvi_hi = quantile(rgvi, probs = 0.75, na.rm = TRUE),
                         gbi_hi = quantile(gbi, probs = 0.75, na.rm = TRUE),
                         gbvi_hi = quantile(gbvi, probs = 0.75, na.rm = TRUE))
      
      # add the spectral data to the original crown segment polygon object
      crowns <- 
        crowns %>% 
        left_join(live_or_dead, by = "object_")
    })
  
  # save the hand classified + spectral data object for later use.
  save(hand_classified_crowns, file = "data/data_output/classified/hand-classified/hand-classified-and-extracted.rda")
}

# Object name is "hand_classified_crowns"
load("data/data_output/classified/hand-classified/hand-classified-and-extracted.rda")
hand_classified_crowns

# pull out just the data (ignoring all spatial components) in order to model the probability of a tree being alive
# given its spectral data
hcc_data <- lapply(hand_classified_crowns, FUN = as.data.frame) %>% 
  do.call(rbind, .)

# fit a logistic regression using the different reflectance values available for each tree
m1 <- glm(live ~ r_mean + g_mean + b_mean + rgi_mean + gbi_mean, family = "binomial", data = hcc_data)
summary(m1)

### This second block of code extracts pixel values from the corresponding orthomosaic for all of the crown segments
### at all of the study sites.

### first create a checklist to see what data_output products have been completed for each site.
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
  dplyr::mutate(ttops_check = file.exists(paste0("data/data_output/site_data/", site, "/", site, "_ttops/", site, "_ttops.shp"))) %>% 
  dplyr::mutate(crowns_check = file.exists(paste0("data/data_output/site_data/", site, "/", site, "_crowns/", site, "_crowns.shp"))) %>% 
  dplyr::mutate(classified_check = file.exists(paste0("data/data_output/classified/model-classified/crown-shapefiles/", site, "_classified-crowns/", site, "_classified-crowns.shp")))

# the sites to process is the intersection between the sites that haven't been model classified and the
# sites for which we have Pix4D data to work with
sites_to_process <-
  sites_checklist %>% 
  dplyr::filter(!classified_check & site %in% all_sites) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()

overwrite <- FALSE

  if (overwrite) {
    
    # First, create shapefile directories for each site in the "sites_to_process" vector
    sites_to_process %>% 
      map(.f = function(site) {
        if (dir.exists(paste0("data/data_output/classified/model-classified/crown-shapefiles/", site, "_classified-crowns"))) {
          unlink(paste0("data/data_output/classified/model-classified/crown-shapefiles/", site, "_classified-crowns"), recursive = TRUE)
        }
        dir.create(paste0("data/data_output/classified/model-classified/crown-shapefiles/", site, "_classified-crowns"))
      })

    # Create a list for the classified crowns to live in.
    classified_crowns <- vector(mode = "list", length = length(sites_to_process))
    
    # loop through all the sites_to_process
    for (i in seq_along(sites_to_process)) {
      
      # `site` is the current site in the loop
      site <- sites_to_process[i]
      
      # Use the velox package to read in the current site's orthomosaic raster
      ortho <- velox::velox(paste0("data/data_output/site_data/", site, "/3_dsm_ortho/2_mosaic/", site, "_x3_transparent_mosaic_group1.tif"))
      # give the velox object bands for its different raster layers
      names(ortho$rasterbands) <- c("r", "g", "b", "x")
      
      # create new bands in the velox object using raster algebra on a pixel-by-pixel basis 
      # to get our vegetation indices that might be helpful (RGI and GBI)
      ortho$rasterbands$rgi <- ortho$rasterbands$r / ortho$rasterbands$g
      ortho$rasterbands$gbi <- ortho$rasterbands$g / ortho$rasterbands$b
      
      # When using the velox objects and adding new bands, you have to update the $nbands list elemnent
      # manually in order for other velox methods to work properly on it
      ortho$nbands <- length(names(ortho$rasterbands))
      
      # read in the corresponding ttop data for the current site
      ttops <- sf::st_read(paste0("data/data_output/site_data/", site, "/", site, "_ttops/", site, "_ttops.shp"))
      
      # read in the corresponding crown segment polygon data for the current site and add some extra information
      # to it
      crowns <- 
        sf::st_read(paste0("data/data_output/site_data/", site, "/", site, "_crowns/", site, "_crowns.shp")) %>% 
        dplyr::mutate(site = site) %>% 
        dplyr::mutate(crs = st_crs(.)$proj4string) %>% 
        dplyr::mutate(treeID = ttops$treeID) %>% 
        dplyr::mutate(ttop_x = st_coordinates(ttops)[, "X"]) %>% 
        dplyr::mutate(ttop_y = st_coordinates(ttops)[, "Y"]) %>% 
        dplyr::mutate(object_ = 1:nrow(.))
      
      # I copy the velox object, because velox methods alter the original object even without using
      # the assignment operator
      ortho_copy <- ortho$copy()
      
      # I learned from here: https://gis.stackexchange.com/questions/286409/fastest-way-to-extract-a-raster-in-r-improve-the-time-of-my-reproducible-code
      # that the velox$extract() method works fastest if you crop the raster to the polygons first.
      ortho_copy$crop(crowns)
      
      # this is the extract operation; I summarize all pixel values in each crown segment using mean()
      extracted_vals <- ortho_copy$extract(crowns, fun = function(x) mean(x, na.rm = TRUE))
      
      # add column names to the resulting object and add a treeID column
      colnames(extracted_vals) <- paste0(names(ortho_copy$rasterbands), "_mean")
      extracted_vals <- as.data.frame(extracted_vals)
      extracted_vals$treeID <- 1:nrow(extracted_vals)
      
      # join the result with the original crowns shapefile, and add a column representing the model prediction of the 
      # probability that each tree is alive given its spectral data
      # store the object in a list
      classified_crowns[[i]] <-
        crowns %>% 
        left_join(extracted_vals, by = "treeID") %>% 
        dplyr::mutate(live_prob = plogis(predict(m1, newdata = ., type = "link")))
      
      # write the spatial object, with the spectral data and its classification result, to a file.
      sf::st_write(classified_crowns[[i]], paste0("data/data_output/classified/model-classified/crown-shapefiles/", site, "_classified-crowns/", site, "_classified-crowns.shp"))
    }
  }

# A simple plot
ggplot(classified_crowns[[i]], aes(x = ttop_x, y = ttop_y, col = live_prob)) +
  geom_point() +
  scale_color_viridis_c(name = "Pr(live)", option = "E") +
  coord_equal() +
  theme_bw() +
  xlab(label = "longitude") +
  ylab(label = "latitude") +
  ggtitle(names(classified_crowns[i])) +
  theme(plot.title = element_text(margin = margin(b = 10)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

