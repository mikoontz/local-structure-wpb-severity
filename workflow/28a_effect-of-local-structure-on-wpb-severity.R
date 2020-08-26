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
                prop_host_ba_s = as.numeric(scale(prop_host_ba)),
                prop_dead_ba = dead_ba / total_ba,
                prop_dead_ba_transformed = (prop_dead_ba * (nrow(.) - 1) + 0.5) / nrow(.),
                prop_dead_pipo_count = dead_count / pipo_and_dead_count,
                prop_dead_pipo_ba = dead_ba / pipo_and_dead_ba,
                prop_dead_pipo_ba_transformed = (prop_dead_pipo_ba * (nrow(.) - 1) + 0.5) / nrow(.),
                prop_dead_count = dead_count / pipo_and_dead_count) %>% 
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
                  pipo_and_dead_mean_height_s:prop_host_count_s +
                  pipo_and_dead_mean_height_s:overall_bapha_s +
                  site_cwd_zscore:prop_host_count_s:pipo_and_dead_mean_height_s +
                  gp(x, y, by = site, scale = FALSE),
                data = adf,
                family = zero_inflated_binomial(),
                iter = 5000,
                warmup = 2000,
                chains = 4,
                cores = 4,
                control = list(adapt_delta = 0.80))
summary(fm1_brms)
(elapsed <- difftime(Sys.time(), start, units = "hours"))
pp_check(fm1_brms, nsamples = 50)

# The final model to use, I think.
# readr::write_rds(x = fm1_brms, path = here::here('analyses/analyses_output/fitted-model_zibinomial_site-cwdZscore_prop-host_pipo-height_overall-tpha_overall-bapha_height-corrected_exact-gp-per-site_200-samples.rds'))

### Basal area consequences as a function of proportion of dead trees
adf <-
  adf %>% 
  dplyr::ungroup() %>% 
  # mutate(dead_count_binned = ggplot2::cut_width(dead_count, width = 5, ordered_results = TRUE, boundary = 0)) %>% 
  dplyr::mutate(dead_count_binned = ggplot2::cut_number(dead_count, n = 4, ordered_results = TRUE, boundary = 0)) %>% 
  dplyr::mutate(pipo_count_binned = ggplot2::cut_number(pipo_and_dead_count, n = 6, ordered_results = TRUE, boundary = 0)) %>% 
  dplyr::mutate(prop_host_count_binned = ggplot2::cut_number(prop_host_count, n = 4, ordered_results = TRUE, boundary = 0)) %>% 
  dplyr::mutate(prop_host_ba_binned = ggplot2::cut_number(prop_host_ba, n = 6, ordered_results = TRUE, boundary = 0))

stacked <- 
  adf %>% 
  dplyr::select(prop_dead_pipo_count, dead_ba, pipo_and_dead_ba, total_ba) %>% 
  tidyr::pivot_longer(cols = c(dead_ba, pipo_and_dead_ba, total_ba), names_to = "type", values_to = "basal_area") %>% 
  dplyr::mutate(type = case_when(type == "dead_ba" ~ "dead tree basal area",
                                 type == "pipo_and_dead_ba" ~ "host tree basal area",
                                 TRUE ~ "total basal area"))

ggplot(stacked, aes(prop_dead_pipo_count, basal_area, color = type)) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(formula = y ~ s(x, bs = "cs", k = 4)) +
  scale_color_viridis_d() +
  labs(x = "Proportion of dead hosts\n(#dead / #hosts)",
       y = "Basal area\n(m^2)",
       color = "Basal area type") +
  theme_bw()

### Zoomed in version so smoothed lines are more legible
# Huge thanks to this Twitter thread for pointing out behavior of 
# ylim() and scale_y_continuous(limits = c()) [i.e., dropping data
# prior to performing geom_smooth() fits]
ba_v_prop_dead_hosts_gg <-
  ggplot(stacked, aes(prop_dead_pipo_count, basal_area, color = type)) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(formula = y ~ s(x, bs = "cs", k = 4)) +
  scale_color_viridis_d() +
  labs(x = "Proportion of dead hosts\n(#dead / #hosts)",
       y = "Basal area\n(m^2)",
       color = "Basal area type") +
  theme_bw() +
  coord_cartesian(ylim = c(0, 2.5))

ba_v_prop_dead_hosts_gg

ggsave(filename = "figures/basal-area-consequences.png", plot = ba_v_prop_dead_hosts_gg)


ggplot(adf, aes(x = prop_dead_pipo_count, y = prop_dead_ba, color = pipo_count_binned)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 7)) +
  scale_color_viridis_d() +
  labs(x = "Proportion of dead hosts\n(#dead / #hosts)",
       y = "Proportion of dead basal area\n(dead basal area / total basal area)",
       color = "Number\nof host trees") +
  theme_bw()

