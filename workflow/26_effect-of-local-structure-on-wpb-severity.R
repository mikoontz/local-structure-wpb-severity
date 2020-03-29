# Purpose: build the primary paper analysis

library(tidyverse)
library(raster)
library(sf)
library(here)
library(mgcv)
library(brms)
library(future)

if(file.exists(here::here("data/data_drone/L4/data-from-rasterized-classified-trees.csv"))) {
  
  data_from_rasterized_trees <- 
    readr::read_csv(here::here("data/data_drone/L4/data-from-rasterized-classified-trees.csv"))
} else {
  stop("You need to extract the data from the rasterized version of the classified trees! See the workflow/25_rasterize-classified-trees.R script.")
}

glimpse(data_from_rasterized_trees)

cwd_data <- readr::read_csv("data/data_output/cwd-data.csv")

center_param <- TRUE
scale_param <- TRUE

# Center and scale predictors
analysis_df <-
  data_from_rasterized_trees %>% 
  dplyr::left_join(cwd_data, by = "site") %>% 
  as_tibble() %>% 
  dplyr::mutate_at(.vars = vars(-x, -y, -local_cwd, -local_cwd_zscore, -site, -forest, -elev, -rep, -crs, -site_cwd, -site_cwd_zscore), .funs = list(s = function(x) {scale(x, center = center_param, scale = scale_param)[, 1]}))
                     
# Implement a zero-inflated binomial with an exact Gaussian process to
# estimate the spatial autocorrelation

# exact GP 200 samples ----------------------------------------------------
# 9.8 hours
# Use height as a predictor instead of DBH, which requires incorporating
# more error into the estimate (because we translate height to DBH using
# data-derived allometric relationships that have plenty of variability both
# in using the species-specific equations, in using the live-only equation
# for PIPO, and then in the allometric relationship itself)

set.seed(1356) # Only meaningful for when randomly subsetting data
adf <-
  analysis_df %>% 
  dplyr::mutate(prop_host = pipo_and_dead_count / total_count,
                prop_host_s = as.numeric(scale(prop_host))) %>% 
  dplyr::filter(pipo_and_dead_count != 0) %>% 
  dplyr::mutate(site = as.factor(site)) %>% 
  dplyr::group_by(site) %>%
  dplyr::sample_n(200)

(start <- Sys.time())
fm1_brms <- brm(dead_count | trials(pipo_and_dead_count) ~ 
                   site_cwd_zscore*prop_host_s*pipo_and_dead_mean_height_s +
                   prop_host_s*overall_tpha_s +
                   site_cwd_zscore*overall_tpha_s +
                   gp(x, y, by = site, scale = FALSE),
                 data = adf,
                 family = zero_inflated_binomial(),
                 iter = 4000,
                 chains = 4,
                 cores = 4,
                 control = list(adapt_delta = 0.90))
summary(fm1_brms)
(elapsed <- Sys.time() - start)
pp_check(fm1_brms, nsamples = 50)

saveRDS(fm1_brms, here::here('analyses/analyses_output/fitted-model_zibinomial_site-cwdZscore_prop-host_pipo-height_overall-tpha_exact-gp-per-site_200-samples.rds'))

