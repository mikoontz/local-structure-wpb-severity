# Size distribution of our classified trees

library(tidyverse)
library(sf)
library(mgcv)
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
  dplyr::select(treeID, plot, live, height) %>% 
  dplyr::mutate(ag = "a")

gtrees <-
  ground_trees %>% 
  sf::st_drop_geometry() %>% 
  dplyr::select(tree, plot, live, height) %>% 
  dplyr::rename(treeID = tree) %>% 
  dplyr::mutate(ag = "g")

trees <- 
  rbind(atrees, gtrees)

mean_heights <-
  trees %>% 
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
  dplyr::mutate(live = as.factor(live)) %>% 
  dplyr::mutate(height_raw = mean_height_a,
                height = mean_height_g) %>% 
  dplyr::mutate(height_raw_s = scale(height_raw),
                height_s = scale(height))

#### Plot to show effect of live/dead on tree height estimates

ggplot(mean_heights_wide, aes(x = mean_height_a, y = mean_height_g, color = live)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ 1 - s(x, bs = "cs")) +
  geom_abline(slope = 1, intercept = 0)

ggplot(mean_heights_wide, aes(x = mean_height_a, y = mean_height_g, color = live)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_abline(slope = 1, intercept = 0)

ggplot(mean_heights_wide, aes(x = mean_height_a, y = mean_height_g, color = live)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0)

fm1 <- mgcv::gam(formula = height ~ live + s(height_raw, by = live), data = mean_heights_wide)
summary(fm1)
gratia::draw(fm1)

fm2 <- stats::lm(formula = height ~ live * height_raw + 0, data = mean_heights_wide)
summary(fm2)

fm3 <- stats::lm(formula = height ~ height_raw + 0, data = mean_heights_wide[mean_heights_wide$live == 0, ])
summary(fm3)

fm4 <- stats::lm(formula = height ~ height_raw + 0, data = mean_heights_wide[mean_heights_wide$live == 1, ])
summary(fm4)

preds0 <- 
  data.frame(height_raw = 0:50, live = 0) %>%  
  modelr::add_predictions(model = fm3, var = "height")

preds1 <- 
  data.frame(height_raw = 0:50, live = 1) %>%  
  modelr::add_predictions(model = fm4, var = "height")

preds <- rbind(preds0, preds1)

plot(x = mean_heights_wide$mean_height_a, y = mean_heights_wide$mean_height_g, pch = 19, col = c("brown", "darkgreen")[as.factor(mean_heights_wide$live)])
lines(x = preds0$height_raw, y = preds0$height, col = "brown")
lines(x = preds1$height_raw, y = preds1$height, col = "darkgreen")
abline(a = 0, b = 1, lty = 2, lwd = 3)

#### Correct the individual tree heights using the plot level averages;

#### Get trees identified (segmented and classified) from the air
all_trees <- 
  sf::st_read(here::here("data", "data_drone", "L3b", "model-classified-trees_all.gpkg"), 
              stringsAsFactors = FALSE)  

# Source in the allometric scaling models based on the ground data
# object is called `allometry_models`
source(here::here("workflow/14_allometric-scaling-models.R"))
allometry_models

# Use the corrected height to estimate basal area and dbh with allometry equations
all_trees <- 
  all_trees %>% 
  dplyr::mutate(height_raw = height,
                estimated_dbh_raw = estimated_dbh) %>% 
  modelr::add_predictions(model = fm1, var = "height") %>% 
  dplyr::mutate(height = ifelse(live == 1, yes = height_raw, no = height)) %>% 
  dplyr::left_join(allometry_models, by = "species") %>% 
  dplyr::mutate(model = ifelse(live == 0, 
                               yes = allometry_models %>% 
                                 dplyr::filter(species == "pipo") %>% 
                                 pull(model), 
                               no = model)) %>% 
  dplyr::do(modelr::add_predictions(., model = first(.$model), var = "estimated_dbh")) %>% 
  dplyr::select(-model) %>% 
  dplyr::mutate(estimated_ba = (estimated_dbh / 2)^2 * pi / 10000)

sf::st_write(obj = all_trees, dsn = here::here("data", "data_drone", "L3b", "model-classified-trees_all_height-corrected.gpkg"), append = FALSE)  
