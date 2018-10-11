# General and basic summary statistics for project

library(sf)
library(tidyverse)
library(purrr)

all_sites <- substr(list.files("data/data_output/classified/model-classified/crown-shapefiles/"), start = 1, stop = 9)

# all_site_bounds <- 
#   all_sites %>% 
#   map(.f = function(site) {
#     site_bounds <-
#       sf::st_read(paste0("data/data_output/site_data/", site, "/", site, "_mission-footprint/", site, "_site-bounds.geoJSON"))
#   }) %>% 
#   setNames(all_sites)
# 
# all_site_areas <-
#   all_site_bounds %>% 
#   map2_df(.y = all_sites, .f = function(x, y) {
#     site_area <-
#       data.frame(site = y, area = st_area(x), stringsAsFactors = FALSE)
#   })

classified_forest <- 
  all_sites %>% 
  map(.f = function(site) {
    crowns <-
      sf::st_read(paste0("data/data_output/classified/model-classified/crown-shapefiles/", site, "_classified-crowns/", site, "_classified-crowns.shp"))
      
    forest <-
      crowns %>% 
      mutate(forest = substr(site, start = 1, stop = 4)) %>% 
      mutate(elev = as.numeric(paste0(substr(site, start = 6, stop = 6), "000"))) %>% 
      mutate(rep = as.numeric(substr(site, start = 9, stop = 9))) %>% 
      mutate(site = as.character(site)) %>%
      mutate(crown_poly = geometry) %>% 
      as.data.frame() %>% 
      dplyr::select(-geometry) %>% 
      st_as_sf(coords = c("ttop_x", "ttop_y"), remove = FALSE, sf_column_name = "geometry") %>% 
      st_set_crs(st_crs(crowns))
    
    site_bounds <-
      sf::st_read(paste0("data/data_output/site_data/", site, "/", site, "_mission-footprint/", site, "_site-bounds.geoJSON")) %>% 
      sf::st_transform(st_crs(crowns)) %>% 
      sf::st_geometry()
    
    forest <-
      forest %>% 
      st_geometry() %>% 
      st_union() %>% 
      st_voronoi() %>% 
      st_cast() %>% 
      st_intersection(st_geometry(site_bounds)) %>%  
      data.frame(geometry = .) %>% 
      st_sf() %>% 
      dplyr::mutate(voronoi_area = st_area(.)) %>% 
      dplyr::mutate(voronoi_poly = geometry) %>% 
      st_join(forest, .) %>% 
      st_set_geometry(value = "geometry")
    
    if (forest$forest[1] == "sequ") {
      forest <-
        forest %>% 
        mutate(elev_relative = ifelse(elev == 4000, yes = "lo", no = ifelse(elev == 5000, yes = "mi", no = "hi")))
    } else {
      forest <-
        forest %>% 
        mutate(elev_relative = ifelse(elev == 3000, yes = "lo", no = ifelse(elev == 4000, yes = "mi", no = "hi")))
    }
    
    return(forest)
  
    }) %>% 
  setNames(all_sites)

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
  dplyr::mutate(classified_check = file.exists(paste0("data/data_output/classified/model-classified/crown-shapefiles/", site, "_classified-crowns/", site, "_classified-crowns.shp"))) %>% 
  dplyr::mutate(augmented_classified_check = file.exists(paste0("data/data_output/classified/model-classified/augmented-crowns/", site, "_augmented-crowns.rds")))

# the sites to process is the intersection between the sites that haven't been model classified and the
# sites for which we have Pix4D data to work with
sites_to_process <-
  sites_checklist %>% 
  dplyr::filter(!augmented_classified_check & site %in% all_sites) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()

overwrite <- FALSE

if (overwrite) {
  if (dir.exists("data/data_output/classified/model-classified/augmented-crowns")) {
    unlink("data/data_output/classified/model-classified/augmented-crowns", recursive = TRUE)
  }
  dir.create("data/data_output/classified/model-classified/augmented-crowns/") 
}

for (i in seq_along(sites_to_process)) {
  idx <- which(names(classified_forest) == sites_to_process[i])
  assign(x = sites_to_process[i], value = classified_forest[[idx]])
  saveRDS(object = get(sites_to_process[i]), file = paste0("data/data_output/classified/model-classified/augmented-crowns/", sites_to_process[i], "_augmented-crowns.rds"))
}

