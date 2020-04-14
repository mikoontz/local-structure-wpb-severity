library(sf)
library(tidyverse)

# Pseudocode
# Calculate the "validation metrics" for the ground data at the plots that are visible in the drone-derived orthophotos
#
# (Explanation: we define "validation metrics" to be things that we care about in the forest
# By picking 2+ "validation metrics" that are somewhat non-correlated and then matching them to the same metrics calculated
# using the treetop detection/crown segmentation algorithm, we can compare the quality of different algorithms.
# The greater the combined correlation of the "validation metrics" between algorithm-derived metrics and ground data-derived
# metrics, the better the algorithm.)
#
# Candidate "validation metrics" include:
# tree density within a plot
# live tree density
# dead tree density
# tree size distribution within a plot (25th and 75th percentile height; mean height)
# proportion of live trees
# Voronoi polygon areas might not work due to edge effects, but let's try mean distance to 3 nearest trees

# Establish a method for detecting tree tops and segmenting crowns
# Work through all of the remotely-visible plots at each site (144 in total) applying the ttops and crown segmentation approach

# First get the formatted ground data; the object is called `d`
if(!file.exists("data/data_drone/L1/ground-trees.gpkg")) {
  source("workflow/12_make-ground-trees-spatial.R")
}

ground_trees <- sf::st_read("data/data_drone/L1/ground-trees.gpkg")

# Use the plot radius to determine the area of the plot and thus the tree density
plot_radius <- sqrt((66*66) / pi) * 12 *2.54 / 100

# Get the validation metrics for the ground plot trees (that hadn't fallen by the time of the flights)
ground_tree_summary <-
  ground_trees %>% 
  dplyr::filter(is.na(year_fall)) %>% 
  dplyr::group_by(plot) %>% 
  dplyr::summarize(ttops_method = "ground",
                   total_tree_count = n(),
                   live_tree_count = sum(live),
                   dead_tree_count = n() - sum(live),
                   total_density_tph = (n() * 10000) / (pi * plot_radius ^ 2),
                   live_density_tph = (sum(live) * 10000) / (pi * plot_radius ^ 2),
                   dead_density_tph = ((n() - sum(live)) * 10000) / (pi * plot_radius ^ 2),
                   live_proportion = sum(live) / n(),
                   height_lwr_25 = quantile(height, prob = 0.25),
                   height_mean = mean(height),
                   height_upr_25 = quantile(height, prob = 0.75),
                   nn_1_mean = mean(nn1),
                   nn_2_mean = mean(nn2),
                   nn_3_mean = mean(nn3),
                   tree_count_above_15m = sum(height >= 15),
                   tree_count_below_15m = sum(height < 15)) %>% 
  sf::st_drop_geometry() %>% 
  dplyr::mutate(elapsed_time = NA) %>% 
  dplyr::select(plot, ttops_method, elapsed_time, everything())

ground_trees %>% 
  group_by(live) %>% 
  summarize(min = min(height),
            `10` = quantile(height, probs = 0.10),
            `20` = quantile(height, probs = 0.20),
            `30` = quantile(height, probs = 0.30),
            `40` = quantile(height, probs = 0.40),
            mean = mean(height),
            `50` = quantile(height, probs = 0.50),
            `60` = quantile(height, probs = 0.60),
            `70` = quantile(height, probs = 0.70),
            `80` = quantile(height, probs = 0.80),
            `90` = quantile(height, probs = 0.90),
            max = max(height),
            ecdf_05 = ecdf(height)(5),
            ecdf_10 = ecdf(height)(10),
            ecdf_15 = ecdf(height)(15),
            ecdf_20 = ecdf(height)(20),
            ecdf_10_dbh = ecdf(dbh)(10),
            ecdf_15_dbh = ecdf(dbh)(15),
            ecdf_20_dbh = ecdf(dbh)(20),
            ecdf_25_dbh = ecdf(dbh)(25))

ggplot(ground_trees, aes(x = dbh, y = height, color = species, lty = as.factor(live))) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(color = "Species",
       lty = "Live (1) or dead (0)?")

ggplot(ground_trees %>% filter(species == "PIPO"), aes(x = dbh, y = height, color = as.factor(live))) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(color = "Live (1) or dead (0)?")

# Now we have a summary of what all of the ground plots looks like with respsect
# to a few of these plot-level characteristics.
# Now we need to employ a bunch of segmentation algorithms to find out which
# one matches best to the ground data.

readr::write_csv(ground_tree_summary, here::here("analyses/analyses_output/ground-tree-summary.csv"))
