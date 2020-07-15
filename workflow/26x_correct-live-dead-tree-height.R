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

### Using the individual-tree level corrections
# corrections <- 
#   readr::read_csv("data/data_drone/L3b/height-corrections-individual-tree-level.csv") %>% 
#   dplyr::mutate(atreeID = paste0(site, "_", air_tree_id),
#                 gtreeID = paste0(plot, "_", ground_tree_id)) %>% 
#   dplyr::select(atreeID, gtreeID)
# 
# correction_table <-
#   atrees %>% 
#   dplyr::inner_join(corrections, by = c(treeID = "atreeID"), live) %>% 
#   dplyr::left_join(dplyr::select(gtrees, -live), by = c(gtreeID = "treeID"), live)
# 
# ggplot(correction_table, aes(x = height_g, y = height_a, color = as.factor(live))) +
#   geom_point() +
#   geom_smooth() +
#   geom_abline(slope = 1, intercept = 0)
# 
# correction_table %>% 
#   mutate(ratio = height_a / height_g) %>% 
#   group_by(live) %>% 
#   summarize(mean_ratio = mean(ratio))

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
  dplyr::mutate(live = ifelse(live == 1, yes = "live", no = "dead")) %>% 
  dplyr::mutate(live = as.factor(live))

#### Plot to show any interacting effect of CWD on the live vs. dead height estimates

ggplot(mean_heights_wide, aes(x = mean_height_a, y = mean_height_g, color = live)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_abline(slope = 1, intercept = 0)

ggplot(mean_heights_wide, aes(x = mean_height_a, y = mean_height_g, color = live)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0)

ggplot(mean_heights_wide, aes(x = mean_height_a, y = mean_height_g, color = live)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x)) +
  geom_abline(slope = 1, intercept = 0)

fm1 <- mgcv::gam(formula = mean_height_g ~ live + s(mean_height_a, by = live), data = mean_heights_wide)
summary(fm1)

gratia::draw(fm1)


readr::write_csv(x = preds, path = here::here("analyses", "analyses_output", "tree-height-correction-factor-by-live-dead.csv"))

#### Correct the tree heights
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
