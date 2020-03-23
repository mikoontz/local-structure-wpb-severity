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
  dplyr::mutate(plot_locations_check = ifelse(site %in%  c("eldo_3k_2", "eldo_3k_3", "eldo_4k_2"),
                                              yes = file.exists(paste0("data/data_output/site_data/", site, "/", site, "_plot-locations/", site, "_plot-locations.shp")),
                                              no = file.exists(paste0("data/data_output/site_data/", site, "/", site, "_plot-locations_re/", site, "_plot-locations_re.shp")))) %>% 
  dplyr::mutate(mission_footprint_check = 
                  file.exists(paste0("data/data_drone/L0/mission-footprint/srtm30m/", site, "_srtm30m.tif")) &
                  file.exists(paste0("data/data_drone/L0/mission-footprint/photo-points/", site, "_photo-points.geoJSON")) &
                  file.exists(paste0("data/data_drone/L0/mission-footprint/site-bounds/", site, "_site-bounds.geoJSON"))) %>% 
  dplyr::mutate(classified_point_cloud_check = file.exists(paste0("data/data_output/site_data/", site, "/", site, "_classified_point_cloud.las"))) %>% 
  dplyr::mutate(dtm_check = file.exists(paste0("data/data_output/site_data/", site, "/", site, "_dtm.tif"))) %>% 
  dplyr::mutate(chm_check = file.exists(paste0("data/data_output/site_data/", site, "/", site, "_chm.tif"))) %>% 
  dplyr::mutate(stacked_ortho_check = file.exists(paste0("data/data_output/site_data/", site, "/", site, "_ortho.tif"))) %>% 
  dplyr::mutate(stacked_index_check = file.exists(paste0("data/data_output/site_data/", site, "/", site, "_index.tif"))) %>% 
  dplyr::mutate(plot_remote_data_check = map_lgl(paste0("data/data_output/site_data/", site, "/", site, "_plot-remote-data"), .f = function(x) length(list.files(x)) > 0)) %>% 
  dplyr::mutate(ground_trees_check = file.exists(paste0("data/data_output/site_data/", site, "/", site, "_ground-trees/", site, "_ground-trees.shp"))) %>% 
  dplyr::mutate(ttops_check = file.exists(paste0("data/data_output/site_data/", site, "/", site, "_ttops/", site, "_ttops.shp"))) %>% 
  dplyr::mutate(crowns_check = file.exists(paste0("data/data_output/site_data/", site, "/", site, "_crowns/", site, "_crowns.shp"))) %>% 
  dplyr::mutate(reflectance_extraction_check = file.exists(paste0("data/data_output/classified/model-classified/crowns-with-reflectance/", site, "_crowns-with-reflectance/", site, "_crowns-with-reflectance.shp"))) %>% 
  dplyr::mutate(rasterized_trees_check = file.exists(paste0("analyses/analyses_output/rasterized-trees/", site, "_rasterized-trees.tif")))

sites_checklist

