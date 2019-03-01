library(sf)
library(tidyverse)

# This script will generate a shapefile of all the locations of the trees that
# were measured in the ~110 field plots that are visible from the aerial imagery.
# We can only use this subset of plots because we base the center point location
# on exactly where it can be found in the aerial imagery, and then the tree
# locations are displacements from that center. With no center visible in an 
# aerial image, the only location data we have for plot centers comes from a
# handheld GPS, which could be off by several meters (and thus the displacement
# of all the trees would also be off)

source(here::here("data/data_carpentry/make-processing-checklist.R"))
source(here::here("data/data_carpentry/format-ground-data.R"))

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
sites_to_overwrite <- ""
sites_checklist$overwrite <- ifelse(sites_to_overwrite == "all", yes = TRUE, no = FALSE)

sites_checklist[sites_checklist$site %in% sites_to_overwrite, "overwrite"] <- TRUE

sites_to_process <- 
  sites_checklist %>% 
  dplyr::filter(!(site %in% unusable_sites)) %>%
  dplyr::filter(overwrite | !ground_trees_check) %>% 
  dplyr::select(site) %>% 
  dplyr::pull() 

current_site = sites_to_process[1]

ground_trees <-
  sites_to_process %>%
  map(.f = function(current_site) {

    if (current_site %in% merged_sites) {
      current_site_plot_locations <-
        sf::st_read(paste0("data/data_output/site_data/", current_site, "/", current_site, "_plot-locations/", current_site, "_plot-locations.shp")) %>%
        mutate(plot = paste(current_site, id, sep = "_")) %>%
        dplyr::arrange(id) %>%
        st_zm()
    } else {
      current_site_plot_locations <-
        sf::st_read(paste0("data/data_output/site_data/", current_site, "/", current_site, "_plot-locations_re/", current_site, "_plot-locations_re.shp")) %>%
        mutate(plot = paste(current_site, id, sep = "_")) %>%
        dplyr::arrange(id) %>%
        st_zm()
    }
    
    current_plot_ground_trees <-
      d %>%
      left_join(current_site_plot_locations, by = "plot") %>%
      dplyr::rename(forst_code = forest_code) %>% 
      dplyr::mutate(delta_x = cospi((1 / 2) - (azm / 180)) * dist) %>%
      dplyr::mutate(delta_y = sinpi((1 / 2) - (azm / 180)) * dist) %>%
      st_as_sf() %>%
      dplyr::mutate(x = st_coordinates(.)[, "X"],
                    y = st_coordinates(.)[, "Y"]) %>%
      as.data.frame() %>%
      dplyr::mutate(new_x = x + delta_x,
                    new_y = y + delta_y) %>%
      filter(!is.na(x)) %>%
      st_as_sf(coords = c("new_x", "new_y")) %>%
      st_set_crs(st_crs(current_site_plot_locations))
    
    nn <-
      st_nn(current_plot_ground_trees, 
            current_plot_ground_trees, 
            k = min(c(4, nrow(current_plot_ground_trees))), 
            returnDist = TRUE, 
            sparse = FALSE, 
            progress = FALSE)$dist %>%
      as.data.frame()
    
    nn <- nn %>% dplyr::select(-1)
    
    if (ncol(nn) == 1) {
      nn <- nn %>% rename(nn1 = V2) %>% mutate(nn2 = NA, nn3 = NA)
    } else if (ncol(nn) == 2) {
      nn <- nn %>% rename(nn1 = V2, nn2 = V3) %>% mutate(nn3 = NA)
    } else if (ncol(nn) == 3) {
      nn <- nn %>% rename(nn1 = V2, nn2 = V3, nn3 = V4) 
    }
    
    current_plot_ground_trees <-
      current_plot_ground_trees %>%
      bind_cols(nn)
    
    if(!dir.exists(here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_ground-trees/")))) {
      dir.create(here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_ground-trees/")))
    }
    
    sf::st_write(current_plot_ground_trees, dsn = here::here(paste0("data/data_output/site_data/", current_site, "/", current_site, "_ground-trees/", current_site, "_ground-trees.shp")), delete_dsn = TRUE)
    
    return(current_plot_ground_trees_3310 = st_transform(current_plot_ground_trees, 3310))
    
}) %>% 
  do.call("rbind", .)
6
