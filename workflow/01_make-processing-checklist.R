library(tidyverse)
library(purrr)

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
  dplyr::filter(!(site %in% c("eldo_4k_3", "stan_4k_3", "stan_5k_3", "sequ_4k_2"))) %>% # These sites did not process well
  dplyr::select(-forest, -elev, -rep) %>% 
  dplyr::mutate(plot_locations_check = file.exists(paste0("data/data_drone/L1/plot-locations/", site, "_plot-locations.gpkg"))) %>% 
  dplyr::mutate(mission_footprint_check = 
                  file.exists(paste0("data/data_drone/L0/mission-footprint/srtm30m/", site, "_srtm30m.tif")) &
                  file.exists(paste0("data/data_drone/L0/mission-footprint/photo-points/", site, "_photo-points.geoJSON")) &
                  file.exists(paste0("data/data_drone/L0/mission-footprint/site-bounds/", site, "_site-bounds.geoJSON"))) %>% 
  dplyr::mutate(stacked_ortho_check = file.exists(paste0("data/data_drone/L1/ortho/", site, "_ortho.tif"))) %>% 
  dplyr::mutate(stacked_index_check = file.exists(paste0("data/data_drone/L2/index/", site, "_index.tif"))) %>% 
  dplyr::mutate(classified_point_cloud_check = file.exists(paste0("data/data_drone/L2/classified-point-cloud/", site, "_classified-point-cloud.las"))) %>% 
  dplyr::mutate(dtm_check = file.exists(paste0("data/data_drone/L2/dtm/", site, "_dtm.tif"))) %>% 
  dplyr::mutate(chm_check = file.exists(paste0("data/data_drone/L2/chm/", site, "_chm.tif"))) %>% 
  dplyr::mutate(L1_plot_remote_data_check = map_lgl(site, .f = function(current_site) {
    x <- paste0("data/data_drone/L1/dsm/cropped-to-plot/")
    length(list.files(x)) > 0
    })) %>% 
  dplyr::mutate(ground_trees_check = file.exists("data/data_drone/L1/ground-trees.gpkg")) %>% 
  dplyr::mutate(ttops_check = file.exists(paste0("data/data_drone/L3a/ttops/", site, "_ttops.gpkg"))) %>% 
  dplyr::mutate(crowns_check = file.exists(paste0("data/data_drone/L3a/crowns/", site, "_crowns.gpkg"))) %>% 
  dplyr::mutate(reflectance_extraction_check = file.exists(paste0("data/data_drone/L3b/crowns-with-reflectance/", site, "_crowns-with-reflectance.gpkg"))) %>% 
  dplyr::mutate(rasterized_trees_check = file.exists(paste0("data/data_drone/L4/rasterized-trees/", site, "_rasterized-trees.tif")))

sites_checklist

