# General and basic summary statistics for project

library(sf)
library(tidyverse)

cwd <- read_csv("data/data_output/cwd-data.csv")

surveyed_area <- 
  sf::st_read("data/data_output/surveyed-area-3310.gpkg", stringsAsFactors = FALSE)

ground_trees <- read_csv("data/data_output/ground-data-for-modeling-summarized-by-site.csv")

trees <- 
  sf::st_read("analyses/analyses_output/classified-trees.geojson", stringsAsFactors = FALSE) %>% 
  st_transform(3310)

trees_split <-
  trees %>% 
  split(f = trees$site)

trees_buffered <-
  trees_split %>% 
  lapply(FUN = function(x) {
    
    site <- unique(x$site)
    site_bounds <- surveyed_area[surveyed_area$site == site, ]
    site_bounds_buffer <- sf::st_buffer(site_bounds, -35)
    
    return(x[site_bounds_buffer, ])
  })

trees_3310 <-
  do.call("rbind", trees_buffered)

summary_table_by_site_ground <-
  ground_trees %>% 
  dplyr::mutate(total_ground_trees = live_count + dead_count,
                ground_mortality = dead_count / total_ground_trees) %>% 
  dplyr::select(site, ground_mortality, total_ground_trees)
  
summary_table_by_site <-
  trees_3310 %>%
  st_drop_geometry() %>%
  group_by(site) %>%
  summarize(prop_mortality = 1 - mean(live),
            n = n())

# summary_table_by_forest <-
#   trees_3310 %>% 
#   st_drop_geometry() %>% 
#   group_by(forest) %>% 
#   summarize(prop_mortality = 1 - mean(live))
# 
# summary_table_by_elev <-
#   trees_3310 %>% 
#   st_drop_geometry() %>% 
#   group_by(elev) %>% 
#   summarize(prop_mortality = 1 - mean(live))
# 
# summary_table_by_forest_elev <-
#   trees_3310 %>% 
#   st_drop_geometry() %>% 
#   group_by(forest, elev) %>% 
#   summarize(prop_mortality = 1 - mean(live))
# 
# summary_table_by_forest
# summary_table_by_elev
# summary_table_by_forest_elev

summary_print_table <-
  summary_table_by_site %>% 
  dplyr::left_join(summary_table_by_site_ground, by = "site") %>% 
  dplyr::left_join(surveyed_area %>% st_drop_geometry(), by = "site") %>% 
  dplyr::left_join(cwd, by = "site") %>% 
  dplyr::mutate(ground_tpha = total_ground_trees / 0.202343,
                air_tpha = n / buffered_survey_area) %>% 
  dplyr::select(site, prop_mortality, ground_mortality, air_tpha, ground_tpha, survey_area, buffered_survey_area, site_cwd, site_cwd_zscore)

write_csv(summary_print_table, path = "analyses/analyses_output/summarized-non-spatial-site-data.csv")
