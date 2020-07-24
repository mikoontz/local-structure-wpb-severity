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
  dplyr::mutate(prop_host_count = pipo_and_dead_count / total_count,
                prop_host_count_s = as.numeric(scale(prop_host_count)),
                prop_host_ba = pipo_and_dead_ba / total_ba,
                prop_host_ba_s = as.numeric(scale(prop_host_ba))) %>% 
  dplyr::filter(pipo_and_dead_count != 0) %>% 
  dplyr::mutate(site = as.factor(site)) %>% 
  dplyr::group_by(site) %>%
  dplyr::sample_n(200)

(start <- Sys.time())
fm1_brms <- brm(dead_count | trials(pipo_and_dead_count) ~ 
                  site_cwd_zscore +
                  prop_host_count_s +
                  pipo_and_dead_mean_height_s +
                  overall_tpha_s +
                  overall_bapha_s +
                  site_cwd_zscore:prop_host_count_s +
                  site_cwd_zscore:pipo_and_dead_mean_height_s +
                  site_cwd_zscore:overall_tpha_s +
                  site_cwd_zscore:overall_bapha_s +
                  prop_host_count_s:overall_tpha_s +
                  prop_host_count_s:pipo_and_dead_mean_height_s +
                  site_cwd_zscore:prop_host_count_s:pipo_and_dead_mean_height_s +
                  gp(x, y, by = site, scale = FALSE),
                data = adf,
                family = zero_inflated_binomial(),
                iter = 4000,
                chains = 4,
                cores = 4,
                control = list(adapt_delta = 0.80))
summary(fm1_brms)
(elapsed <- difftime(Sys.time(), start, units = "hours"))
pp_check(fm1_brms, nsamples = 50)

# The final model to use, I think.
readr::write_rds(x = fm1_brms, path = here::here('analyses/analyses_output/fitted-model_zibinomial_site-cwdZscore_prop-host_pipo-height_overall-tpha_overall-bapha_height-corrected_exact-gp-per-site_200-samples.rds'))


###
# (start <- Sys.time())
# fm1b_brms <- brm(dead_count | trials(pipo_and_dead_count) ~ 
#                   site_cwd_zscore +
#                   prop_host_count_s +
#                   prop_host_ba_s +
#                   total_count_s +
#                   total_ba_s +
#                   site_cwd_zscore:prop_host_count_s +
#                   site_cwd_zscore:prop_host_ba_s +
#                   site_cwd_zscore:total_count_s +
#                   site_cwd_zscore:total_ba_s +
#                   prop_host_count_s:total_count_s +
#                   prop_host_ba_s:total_ba_s +
#                   prop_host_count_s:prop_host_ba_s +
#                   site_cwd_zscore:prop_host_ba_s:total_ba_s +
#                   site_cwd_zscore:prop_count_ba_s:total_count_s +
#                   gp(x, y, by = site, scale = FALSE),
#                 data = adf,
#                 family = zero_inflated_binomial(),
#                 iter = 5000,
#                 chains = 4,
#                 cores = 4,
#                 control = list(adapt_delta = 0.80))
# summary(fm1b_brms)
# (elapsed <- difftime(Sys.time(), start, units = "hours"))
# pp_check(fm1b_brms, nsamples = 50)

(start <- Sys.time())
fm1b_brms <- brm(dead_count | trials(pipo_and_dead_count) ~ 
                  site_cwd_zscore +
                  prop_host_count_s +
                  pipo_and_dead_ba_s +
                  overall_tpha_s +
                  non_pipo_ba_s +
                  site_cwd_zscore:prop_host_count_s +
                  site_cwd_zscore:pipo_and_dead_ba_s +
                  site_cwd_zscore:overall_tpha_s +
                  site_cwd_zscore:non_pipo_ba_s +
                  pipo_and_dead_ba_s:non_pipo_ba_s +
                  prop_host_count_s:overall_tpha_s +
                  prop_host_count_s:pipo_and_dead_ba_s +
                  site_cwd_zscore:prop_host_count_s:pipo_and_dead_ba_s +
                  gp(x, y, by = site, scale = FALSE),
                data = adf,
                family = zero_inflated_binomial(),
                iter = 5000,
                chains = 4,
                cores = 4,
                control = list(adapt_delta = 0.80))
summary(fm1b_brms)
(elapsed <- difftime(Sys.time(), start, units = "hours"))
pp_check(fm1b_brms, nsamples = 50)






(start <- Sys.time())
fm2_brms <- brm(dead_ba / total_ba ~ 
                  site_cwd_zscore +
                  prop_host_count_s +
                  pipo_and_dead_mean_height_s +
                  overall_tpha_s +
                  overall_bapha_s +
                  site_cwd_zscore:prop_host_count_s +
                  site_cwd_zscore:pipo_and_dead_mean_height_s +
                  site_cwd_zscore:overall_tpha_s +
                  site_cwd_zscore:overall_bapha_s +
                  prop_host_count_s:overall_tpha_s +
                  pipo_and_dead_mean_height_s:prop_host_count_s +
                  pipo_and_dead_mean_height_s:overall_bapha_s +
                  site_cwd_zscore:prop_host_count_s:pipo_and_dead_mean_height_s +
                  gp(x, y, by = site, scale = TRUE),
                data = adf,
                family = zero_one_inflated_beta(),
                iter = 4000,
                chains = 4,
                cores = 4,
                control = list(adapt_delta = 0.80))
summary(fm2_brms)
(elapsed <- difftime(Sys.time(), start, units = "hours"))
pp_check(fm2_brms, nsamples = 50)



