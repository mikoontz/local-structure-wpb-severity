# Purpose: build the primary paper analysis

library(tidyverse)
library(raster)
library(sf)
library(here)
library(tictoc)
library(mgcv)
library(brms)
library(lme4)
library(effects)
library(future)


if(file.exists(here::here("analyses/analyses_output/data-from-rasterized-classified-trees.csv"))) {
  
  data_from_rasterized_trees <- 
    readr::read_csv(here::here("analyses/analyses_output/data-from-rasterized-classified-trees.csv"))
} else {
  stop("You need to extract the data from the rasterized version of the classified trees! See the analyses/rasterize-classified-trees.R script.")
}

glimpse(data_from_rasterized_trees)

# Get the CWD data
cwd_data <- readr::read_csv(here::here("data/data_output/cwd-data.csv"))

center_param <- TRUE
scale_param <- TRUE

analysis_df <-
  data_from_rasterized_trees %>% 
  dplyr::left_join(cwd_data, by = "site") %>% 
  as_tibble() %>% 
  dplyr::mutate(pipo_and_dead_tpha_s = scale(pipo_and_dead_tpha, center = center_param, scale = scale_param),
                overall_tpha_s = scale(overall_tpha, center = center_param, scale = scale_param),
                pipo_and_dead_bapha_s = scale(pipo_and_dead_bapha, center = center_param, scale = scale_param),
                overall_bapha_s = scale(overall_bapha, center = center_param, scale = scale_param),
                pipo_and_dead_qmd_s = scale(pipo_and_dead_qmd, center = center_param, scale = scale_param),
                overall_qmd_s = scale(overall_qmd, center = center_param, scale = scale_param),
                live_sdi_ac_s = scale(live_sdi_ac, center = center_param, scale = scale_param),
                pipo_and_dead_sdi_ac_s = scale(pipo_and_dead_sdi_ac, center = center_param, scale = scale_param),
                overall_sdi_ac_s = scale(overall_sdi_ac, center = center_param, scale = scale_param))

set.seed(0314) # Only meaningful for when randomly subsetting data
adf <-
  analysis_df %>% 
  dplyr::filter(pipo_and_dead_count != 0) %>% 
  dplyr::mutate(site = as.factor(site))

# dplyr::group_by(site) %>% 
# dplyr::sample_n(50) # ...OR...
# dplyr::sample_frac(0.5)

# Exact Gaussian Process
# n = 10; total = 320; iter = 2000: 217 seconds (0.06 hours)
# n = 50; total = 1600; iter = 2000: 1454/1570 seconds (0.4 hours)
# n = 100; total = 3200; iter = 2000: 1.6 hours
# n = 200; total = 6400; iter = 2500 = 8.9 hours

# Approximate Gaussian Process
# n = 10; total = 320; iter = 2000; k = 3; adapt_delta = 0.95; 212 seconds (0.059 hours) *7 divergent transitions
# n = 10; total = 320; iter = 2000; k = 10; adapt_delta = 0.95; 342 seconds (0.095 hours) *393 divergent transitions
# n = 50; total = 1600; iter = 2000; k = 3; adapt_delta = 0.95; 1084 seconds (0.3 hours) *39 divergent transitions
# n = 50; total = 1600; iter = 2000; k = 10; adapt_delta = 0.95; 989 seconds (0.275 hours) *48 divergent transitions
# n = 100; total = 3200; iter = 2000; k = 3; adapt_delta = 0.95; 1945 seconds (0.54 hours) *132 divergent transitions
# n = 100; total = 3200; iter = 2000; k = 3; adapt_delta = 0.99; 3816 seconds (1.06 hours) *40 divergent transitions **1598 transitions exceeding max_treedepth
# n = 100; total = 3200; iter = 2000; k = 5; adapt_delta = 0.95; 2807.55 seconds (0.78 hours) *28 divergent transitions
# n = 100; total = 3200; iter = 2000; k = 10; adapt_delta = 0.95; 2125 seconds (0.59 hours)
# frac = 0.5; total = 11145; iter = 2000; k = 10; adapt_delta = 0.95; (4.07 hours)
# all; total = 22289; iter = 2000; k = 10; adapt_delta = 0.95; (15.08 hours)

# (start <- Sys.time())
# fm1_brms <- brm(dead_count | trials(pipo_and_dead_count) ~ 
#                   site_cwd_zscore*pipo_and_dead_tpha_s*pipo_and_dead_qmd_s +
#                   site_cwd_zscore*overall_tpha_s*overall_qmd_s +
#                   gp(x, y, by = site, k = 10, c = 5/4, scale = FALSE),
#                 data = adf,
#                 family = binomial(link = "logit"),
#                 chains = 4,
#                 cores = 4,
#                 control = list(adapt_delta = 0.95))
# summary(fm1_brms)
# (elapsed <- Sys.time() - start)
# 
# prior_summary(fm1_brms)

# saveRDS(fm1_brms, file = here::here("analyses/analyses_output/fitted-model_binomial_site-cwdZscore_pipo-tpha-qmd_overall-tpha-qmd_approx-gp-per-site.rds"))

# Posterior predictive checks
# Indicate a pretty bad model fit. Way underpredicting zeros

# fm1 <- readRDS(here::here("analyses/analyses_output/fitted-model_site-cwdZscore_pipo-tpha-qmd_overall-tpha-qmd_approx-gp-per-site.rds"))
# pp_check(fm1, nsamples = 50)

# Implement a zero-inflated binomial instead
# exact GP 200 samples ----------------------------------------------------
# 7.7 hours
set.seed(0314) # Only meaningful for when randomly subsetting data
adf <-
  analysis_df %>% 
  dplyr::filter(pipo_and_dead_count != 0) %>% 
  dplyr::mutate(site = as.factor(site)) %>% 
  dplyr::group_by(site) %>%
  dplyr::sample_n(200)

(start <- Sys.time())
fm16_brms <- brm(dead_count | trials(pipo_and_dead_count) ~ 
                   site_cwd_zscore*pipo_and_dead_tpha_s*pipo_and_dead_qmd_s +
                   site_cwd_zscore*overall_tpha_s*overall_qmd_s +
                   gp(x, y, by = site, scale = FALSE),
                 data = adf,
                 family = zero_inflated_binomial(),
                 chains = 4,
                 cores = 4,
                 control = list(adapt_delta = 0.95))
summary(fm16_brms)
(elapsed <- Sys.time() - start)
pp_check(fm16_brms, nsamples = 50)

# saveRDS(fm16_brms, here::here('analyses/analyses_output/fitted-model_zibinomial_site-cwdZscore_pipo-tpha-qmd_overall-tpha-qmd_exact-gp-per-site_200-samples.rds'))
fm16_brms <- readRDS(here::here('analyses/analyses_output/fitted-model_zibinomial_site-cwdZscore_pipo-tpha-qmd_overall-tpha-qmd_exact-gp-per-site_200-samples.rds'))

for (i in seq_along(unique(adf$site))) {
  current_site <- unique(adf$site)[i]
  current_site_idx <- which(adf$site == current_site)
  
  resids <- resid(fm16_brms)[current_site_idx, "Estimate"]
  coords <- fm16_brms$data[current_site_idx, c("x", "y")]
  
  spatial_autocor <- Variogram(object = resids, dist(coords))
  # plot(spatial_autocor)
  
}
