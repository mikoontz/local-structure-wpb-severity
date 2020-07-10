# Size distribution of our classified trees

library(tidyverse)
library(sf)
library(modelr)
library(here)

cwd <-
  read_csv(here::here("data", "data_output", "cwd-data.csv"))

#### Read in the trees that fall within the visible plots
air_trees <- 
  sf::st_read(here::here("data", "data_drone", "L3b", "model-classified-trees-within-ground-plots.gpkg"), 
              stringsAsFactors = FALSE)

#### Get ground trees
# Subset ground trees to the 110 plots that are identifiable by air
ground_trees <- 
  sf::st_read(here::here("data", "data_drone", "L1", "ground-trees.gpkg"),
              stringsAsFactors = FALSE) %>% 
  dplyr::filter(is.na(year_fall)) %>%
  dplyr::left_join(cwd, by = "site") %>%
  dplyr::filter(plot %in% unique(air_trees$plot))


#### Drop geometries from air and ground trees
atrees <-
  air_trees %>% 
  sf::st_drop_geometry() %>% 
  dplyr::select(site, plot, live, species, height) %>% 
  dplyr::mutate(ag = "a")

gtrees <-
  ground_trees %>% 
  sf::st_drop_geometry() %>% 
  dplyr::select(site, plot, live, species, height) %>% 
  dplyr::mutate(species = tolower(species)) %>% 
  dplyr::mutate(ag = "g")

trees <- 
  rbind(atrees, gtrees)

mean_heights <-
  trees %>% 
  # dplyr::filter(species == "pipo") %>% 
  dplyr::group_by(plot, live, ag) %>% 
  dplyr::summarize(mean_height = mean(height),
                   n = n())

mean_heights_wide <- 
  mean_heights %>% 
  tidyr::pivot_wider(names_from = "ag", values_from = c("mean_height", "n")) %>% 
  dplyr::mutate(a_minus_g = mean_height_a - mean_height_g,
                a_div_g = mean_height_a / mean_height_g) %>% 
  dplyr::mutate(site = substr(plot, start = 1, stop = 9))%>% 
  dplyr::left_join(cwd, by = "site") %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(live = ifelse(live == 1, yes = "live", no = "dead")) %>% 
  dplyr::mutate(live = as.factor(live))

#### Plot to show any interacting effect of CWD on the live vs. dead height estimates

ggplot(mean_heights_wide[complete.cases(mean_heights_wide), ], aes(x = site_cwd_zscore, y = a_div_g, color = live)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"))

#### Model the effect of the site CWD and live/dead status on the ratio of drone-derived
height_corrections <-
  mean_heights_wide %>% 
  dplyr::group_by(live) %>% 
  dplyr::summarize(a_div_g = mean(a_div_g, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(a_div_g_reciprocal = 1 / a_div_g) %>% 
  dplyr::mutate(live = ifelse(live == "live", yes = 1, no = 0))

height_corrections

readr::write_csv(x = preds, path = here::here("analyses", "analyses_output", "tree-height-correction-factor-by-live-dead.csv"))

#### Correct the tree heights
#### Get trees identified (segmented and classified) from the air
trees <- 
  sf::st_read(here::here("data", "data_drone", "L3b", "model-classified-trees_all.gpkg"), 
              stringsAsFactors = FALSE)  

# Source in the allometric scaling models based on the ground data
# object is called `allometry_models`
source(here::here("workflow/14_allometric-scaling-models.R"))
allometry_models

# Use the corrected height to estimate basal area and dbh with allometry equations
trees <- 
  trees %>% 
  dplyr::left_join(height_corrections, by = "live") %>% 
  dplyr::mutate(height_raw = height,
                height = height_raw * a_div_g_reciprocal) %>% 
  dplyr::left_join(allometry_models, by = "species") %>% 
  dplyr::mutate(model = ifelse(live == 0, 
                               yes = allometry_models %>% 
                                 dplyr::filter(species == "pipo") %>% 
                                 pull(model), 
                               no = model)) %>% 
  dplyr::do(modelr::add_predictions(., model = first(.$model), var = "estimated_dbh_corrected")) %>% 
  dplyr::select(-model) %>% 
  dplyr::mutate(estimated_ba_corrected = (estimated_dbh_corrected / 2)^2 * pi / 10000,
                height_corrected = height,
                height = height_raw)

trees

sf::st_write(obj = trees, dsn = here::here("data", "data_drone", "L3b", "model-classified-trees_all_height-corrected.gpkg"))  
