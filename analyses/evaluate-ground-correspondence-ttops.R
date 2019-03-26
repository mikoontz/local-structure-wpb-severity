# Assess which tree top detection algorithms to further test with segmenting the 
# canopies themselves

library(sf)
library(tidyverse)
library(purrr)
library(viridis)
library(raster)

ttops_summary <- read.csv(here::here("analyses/analyses_output/ttops-summary.csv"),
                          stringsAsFactors = FALSE)

ground <- 
  ttops_summary %>% 
  filter(ttops_method == "ground") %>% 
  dplyr::select(-live_tree_count, -dead_tree_count, -total_density_tph, -live_density_tph, -dead_density_tph, -live_proportion) %>% 
  gather(key = forest_metric, value = ground_value, -plot, -ttops_method) %>% 
  dplyr::select(-ttops_method)

air <- 
  ttops_summary %>% 
  filter(ttops_method != "ground") %>% 
  dplyr::select(-live_tree_count, -dead_tree_count, -total_density_tph, -live_density_tph, -dead_density_tph, -live_proportion) %>% 
  gather(key = forest_metric, value = air_value, -plot, -ttops_method)

air_ground <- 
  left_join(air, ground, by = c("plot", "forest_metric")) %>% 
  dplyr::mutate(diff = air_value - ground_value)

elapsed_time <-
  air_ground %>% 
  dplyr::filter(forest_metric == "elapsed_time") %>% 
  group_by(ttops_method) %>% 
  summarize(mean_time = mean(air_value)) %>% 
  dplyr::rename(mean_algorithm_time_s = mean_time)

air_ground_summary <-
  air_ground %>%
  dplyr::filter(forest_metric != "elapsed_time") %>% 
  group_by(ttops_method, forest_metric) %>% 
  summarize(ground_correlation = cor(air_value, ground_value, use = "na.or.complete"),
            complete_cases = length(which(!is.na(air_value))),
            rmse = sqrt(mean((air_value - ground_value)^2, na.rm = TRUE)),
            me = mean(air_value - ground_value, na.rm = TRUE),
            med_error = median(air_value - ground_value, na.rm = TRUE),
            min = min(air_value - ground_value, na.rm = TRUE),
            max = max(air_value - ground_value, na.rm = TRUE)) %>% 
  dplyr::left_join(elapsed_time)

# Some example comparisons of correlations between air and ground data as well as RMSE between
# air and ground data.
air_ground_summary %>% filter(forest_metric == "total_tree_count") %>% arrange(desc(ground_correlation))
air_ground_summary %>% filter(forest_metric == "total_tree_count") %>% arrange(rmse)

air_ground_summary %>% filter(forest_metric == "height_mean") %>% arrange(desc(ground_correlation))
air_ground_summary %>% filter(forest_metric == "height_mean") %>% arrange(rmse)

air_ground_summary %>% filter(forest_metric == "nn_1_mean") %>% arrange(desc(ground_correlation))
air_ground_summary %>% filter(forest_metric == "nn_1_mean") %>% arrange(rmse)


air_ground_summary %>% filter(forest_metric == "tree_count_above_15m") %>% arrange(rmse)
air_ground_summary %>% filter(forest_metric == "tree_count_below_15m") %>% arrange(rmse)


# The percentage threshold for assessing whether the ground correlation or RMSE of a 
# particular algorithm is within the top `threshold` percentile of all models tested 
# and whether the ground correlationa or RMSE values of a particular algorithm
# are within `threshold` percent of the best performing algorithm
threshold <- 5

# Generate the comparisons across algorithms using ground correlation and RMSE as the key
# metrics. See which algorithms fall into the Top `threshold` percentile or whether
# their response metrics are within `threshold` percent of the best performing algorithm
comps <-
  air_ground_summary %>% 
  dplyr::filter(complete_cases > 100) %>% 
  dplyr::group_by(forest_metric) %>% 
  dplyr::mutate(ground_cor_top_threshold_pct = ifelse(ground_correlation > quantile(ground_correlation, probs = (100 - threshold) / 100, na.rm = TRUE), yes = 1, no = 0)) %>% 
  dplyr::mutate(ground_cor_within_threshold_pct = ifelse(ground_correlation > (((100 - threshold) / 100) * max(ground_correlation, na.rm = TRUE)), yes = 1, no = 0)) %>% 
  dplyr::mutate(rmse_top_threshold_pct = ifelse(rmse < quantile(rmse, probs = threshold / 100, na.rm = TRUE), yes = 1, no = 0)) %>% 
  dplyr::mutate(rmse_within_threshold_pct = ifelse(rmse < ((100 + threshold) / 100) * min(rmse, na.rm = TRUE), yes = 1, no = 0))

# To summarize the models, count the number of forest metrics in each algorithm implementation
# that are rated highly. Can the tree top detection method perform well for several
# different forest metrics that we care about?
best_ttops_detection <-
  comps %>% 
  dplyr::filter(forest_metric %in% c("total_tree_count", "tree_count_above_15m", "height_mean", "height_upr_25", "height_lwr_25", "nn_1_mean", "nn_2_mean")) %>% 
  dplyr::group_by(ttops_method) %>% 
  summarize(ground_cor_top_threshold_pct = sum(ground_cor_top_threshold_pct),
            ground_cor_within_threshold_pct = sum(ground_cor_within_threshold_pct),
            rmse_top_threshold_pct = sum(rmse_top_threshold_pct),
            rmse_within_threshold_pct = sum(rmse_within_threshold_pct)) %>% 
  dplyr::mutate(both_within_threshold_pct = ground_cor_within_threshold_pct + rmse_within_threshold_pct,
                both_top_threshold_pct = ground_cor_top_threshold_pct + rmse_top_threshold_pct) %>% 
  dplyr::arrange(desc(both_within_threshold_pct)) %>% 
  dplyr::left_join(elapsed_time)

single_best_algorithm <- best_ttops_detection %>% slice(1) %>% pull(ttops_method)

best_algorithm_deets <-
  comps %>% 
  dplyr::filter(ttops_method == single_best_algorithm)

algorithm_details <-
  best_ttops_detection %>% 
  separate(col = ttops_method, into = c("algorithm", "var1", "var1_val", "var2", "var2_val", "zu", "zu_val", "R", "R_val", "speedUp", "speedUp_val"), sep = "_") %>% 
  dplyr::select(1:11)

algorithm_summary <-
  best_ttops_detection %>% 
  separate(col = ttops_method, into = c("algorithm", "var1", "var1_val", "var2", "var2_val", "zu", "zu_val", "R", "R_val", "speedUp", "speedUp_val"), sep = "_") %>% 
  dplyr::select(1:11) %>% 
  dplyr::group_by(algorithm) %>% 
  tally()

