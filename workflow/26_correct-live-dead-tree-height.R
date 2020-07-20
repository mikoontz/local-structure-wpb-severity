# Size distribution of our classified trees

library(tidyverse)
library(sf)
library(modelr)
library(here)

# Read in the data frame of site/cwd/cwd z-score information
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
  dplyr::select(treeID, live, height, ch_area) %>% 
  dplyr::rename(height_raw = height)

gtrees <-
  ground_trees %>% 
  sf::st_drop_geometry() %>% 
  dplyr::select(tree, live, height) %>% 
  dplyr::rename(treeID = tree)

### Read in the ~450 trees on the ground positively identified as being also
### identified from the drone in order to match the drone-derived height
### with the ground-derived height
corrections <- 
  readr::read_csv("data/data_drone/L3b/height-corrections-individual-tree-level.csv") %>% 
  dplyr::mutate(atreeID = paste0(site, "_", air_tree_id),
                gtreeID = paste0(plot, "_", ground_tree_id)) %>% 
  dplyr::select(atreeID, gtreeID)

# Join the drone- and ground-derived tree information using the `corrections`
# dataframe as the in-between connection between the treeID of a ground-derived
# tree and the treeID of a drone-derived tree
correction_table <-
  corrections %>% 
  dplyr::left_join(atrees, by = c("atreeID" = "treeID")) %>% 
  dplyr::left_join(dplyr::select(gtrees, -live), by = c("gtreeID" = "treeID"), live) %>% 
  dplyr::mutate(live = as.factor(ifelse(live == 0, yes = "dead", no = "live")),
                site = substr(atreeID, start = 1, stop = 9)) %>% 
  dplyr::left_join(cwd, by = "site") %>% 
  dplyr::as_tibble()

## Plot of all the data shows 4 key things:
## 1) We measure live tree height pretty darn well. (Fairly remarkable 
##    considering the point cloud has to be normalized by an interpolated
##    raster layer of what is considered the ground which is itself derived
##    from the point cloud and a cloth simulator filter algorithm)
## 2) Dead trees measured as pretty short from the drone may in fact be quite 
##    tall. This is almost certainly due to the snag itself being much less
##    voluminous and thus harder for the photogrammetry process to pinpoint
##    where the top really is
## 3) Dead trees measured as pretty tall from the drone tend to be measured
##    reasonably well (one more meter of drone-measured height seems to equal
##    one more meter of ground-measured height)
## 4) Even when the slope of ground- vs. drone-measured dead tree heights is 1 
##    (above drone-measured heights of ~18), the drone still *underestimates* 
##    dead tree height relative to ground-measured heights by a constant amount 
##    (~5 meters) as a function of the drone-measured heights. That is, we
##    should probably add ~5 meters to the heights of the drone-measured dead 
##    trees
## 5) The drone seems to *overestimate* the height of the live trees relative
##    to their ground measurements, but only slightly. We should probably 
##    subtract a little bit from the drone-measured heights of the live trees.

ggplot(correction_table, aes(x = height_raw, y = height, color = live)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Height measured from drone (m)",
       y = "Height measured from ground (m)",
       title = "Comparison of drone-derived and field-derived tree heights",
       color = "Live or dead") +
  theme_bw()

## As a *conservative estimate* of the calibrations needed to make the drone-
## measured heights comparable to the ground-measured heights, we can use the
## relationship between the drone-measured and ground-measured heights of the
## tallest of the live and dead trees. Both slopes (for both live and dead 
## trees) is very close to 1, suggesting 1 more meter of drone-measured height
## is associated with 1 more meter of ground-measured height. 

## Here, we subset to the trees greater than 18 meters in height
conservative_height_calibration <-
  correction_table %>% 
  dplyr::filter(height_raw > 15)

## Check the plot and the GAM fits look good!
ggplot(conservative_height_calibration, aes(x = height_raw, y = height, color = live)) +
  geom_point() +
  geom_smooth(method = "gam") +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Height measured from drone (m)",
       y = "Height measured from ground (m)",
       title = "Comparison of drone-derived and field-derived tree heights",
       color = "Live or dead") +
  theme_bw()

## The linear fits also look good! We will use these for simplicity when 
## calibrating the drone-measured tree heights
ggplot(conservative_height_calibration, aes(x = height_raw, y = height, color = live)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Height measured from drone (m)",
       y = "Height measured from ground (m)",
       title = "Comparison of drone-derived and field-derived tree heights",
       color = "Live or dead") +
  theme_bw()


## The linear slopes are mostly parallel and track the 1:1 line. However, there
## is a persistent bias in the drone-measured heights such that dead tree 
## heights are underestimated and live tree heights are overestimated relative
## to the ground-measured trees.
## How should we adjust the drone-measured tree heights? Let's build a model
## to investigate how the live/dead status and the drone-measured tree height
## predict the *difference* between the ground- and drone-derived tree heights
## Then we'll be able to use this model, the live/dead status of each drone-
## detected tree, and the drone-measured height of each tree to give us a 
## "correction" amount-- the value that we should add to each drone-detected
## tree.
## Again, it is important to recall that this is a *conservative* calibration
## More than likely, the drone-measured tree heights should be augmented further
## though that is outside the scope of this particular study. And would be an
## excellent avenue for new research! (i.e., how can we better detect and 
## measure dead trees using Structure from Motion techniques?)
fm1 <- lm(formula = (height - height_raw) ~ height_raw * live, 
           data = correction_table)
summary(fm1)

fm2 <- lm(formula = (height - height_raw) ~ height_raw * live, 
          data = conservative_height_calibration)
summary(fm2)

ggplot(conservative_height_calibration, aes(x = height_raw, y = height - height_raw, color = live)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Height measured from drone (m)",
       y = "Height measured from ground (m)",
       title = "Comparison of drone-derived and field-derived tree heights",
       color = "Live or dead") +
  theme_bw()

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
  dplyr::mutate(live = ifelse(live == 1, yes = "live", no = "dead"),
                height_raw = height,
                estimated_dbh_raw = estimated_dbh) %>% 
  modelr::add_predictions(model = fm2, var = "height_correction") %>% 
  dplyr::mutate(height = height_raw + height_correction) %>% 
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
