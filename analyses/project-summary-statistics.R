# General and basic summary statistics for project

library(sf)
library(tidyverse)
library(purrr)
library(lme4)
library(effects)

all_sites <- list.files("data/data_output/classified/model-classified/augmented-crowns", full.names = TRUE)
all_site_names <- substr(list.files("data/data_output/classified/model-classified/augmented-crowns"), start = 1, stop = 9)

augmented_crowns <- 
  lapply(all_sites, FUN = readRDS) %>% 
  setNames(all_site_names)

all_site_bounds <-
  all_site_names %>% 
  map2(.y = augmented_crowns, .f = function(site, forest) {
    site_bounds <-
      sf::st_read(paste0("data/data_output/site_data/", site, "/", site, "_mission-footprint/", site, "_site-bounds.geoJSON")) %>% 
      sf::st_transform(st_crs(forest)) %>% 
      sf::st_geometry()
  }) %>% 
  setNames(all_site_names)

all_site_areas <-
  all_site_bounds %>% 
  map2_df(.y = all_site_names, .f = function(site_bounds, site_name) {
    data.frame(site = site_name, site_area = as.numeric(st_area(site_bounds)), stringsAsFactors = FALSE)
  })

# pull out just the data (ignoring all spatial components) in order to model the probability of a tree being alive
# given its spectral data
hcc_data <- lapply(augmented_crowns, FUN = as.data.frame) %>%
  do.call(rbind, .)

head(hcc_data)

hcc_data <- 
  hcc_data %>% 
  mutate(live = ifelse(live_prob > 0.5, yes = 1, no = 0))

summarized_hcc_data <-
  hcc_data %>% 
  group_by(forest, elev, rep, site, live) %>% 
  summarize(count = n()) %>% 
  spread(key = live, value = count) %>% 
  rename(dead = `0`, live = `1`) %>% 
  mutate(total_trees = dead + live) %>% 
  mutate(survivorship = live / total_trees) %>% 
  mutate(mortality = dead / total_trees) %>% 
  left_join(all_site_areas, by = "site") %>% 
  mutate(density = total_trees / (site_area / 10000))

ggplot(summarized_hcc_data, aes(x = density, y = mortality)) + 
  geom_point() +
  geom_smooth(method = "lm")

summarized_hcc_data
sum(summarized_hcc_data$total_trees)