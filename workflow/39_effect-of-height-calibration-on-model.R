# Show effect of height calibration on model output

library(tidyverse)
library(raster)
library(sf)
library(here)
library(mgcv)
library(brms)
library(future)
library(tidybayes)


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


set.seed(1356) # Only meaningful for when randomly subsetting data
adf <-
  analysis_df %>% 
  dplyr::mutate(prop_host = pipo_and_dead_count / total_count,
                prop_host_s = as.numeric(scale(prop_host)),
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

# Set the raw height to be the mean height in order to keep the model parameter
# names the same
adf$pipo_and_dead_mean_height_s <- adf$pipo_and_dead_mean_height_raw_s

(start <- Sys.time())
fm1_brms <- brm(dead_count | trials(pipo_and_dead_count) ~ 
                  site_cwd_zscore +
                  prop_host_s +
                  pipo_and_dead_mean_height_s +
                  overall_tpha_s +
                  overall_bapha_s +
                  site_cwd_zscore:prop_host_s +
                  site_cwd_zscore:pipo_and_dead_mean_height_s +
                  site_cwd_zscore:overall_tpha_s +
                  site_cwd_zscore:overall_bapha_s +
                  prop_host_s:overall_tpha_s +
                  prop_host_s:pipo_and_dead_mean_height_s +
                  site_cwd_zscore:prop_host_s:pipo_and_dead_mean_height_s +
                  gp(x, y, by = site, scale = FALSE),
                data = adf,
                family = zero_inflated_binomial(),
                iter = 4000,
                chains = 4,
                cores = 4,
                control = list(adapt_delta = 0.90))
summary(fm1_brms)
(elapsed <- difftime(Sys.time(), start, units = "hours"))
pp_check(fm1_brms, nsamples = 50)

readr::write_rds(x = fm1_brms, path = here::here('analyses/analyses_output/fitted-model_zibinomial_site-cwdZscore_prop-host_pipo-height_overall-tpha_overall-bapha_exact-gp-per-site_200-samples_uncalibrated.rds'))

### Build a halfeye plot
fm1_brms <- readr::read_rds(path = here::here('analyses/analyses_output/fitted-model_zibinomial_site-cwdZscore_prop-host_pipo-height_overall-tpha_overall-bapha_exact-gp-per-site_200-samples_uncalibrated.rds'))

covariates <- rownames(fixef(fm1_brms))

samps <- 
  posterior_samples(fm1_brms) %>% 
  dplyr::select(tidyselect::all_of(paste0("b_", covariates))) %>% 
  dplyr::sample_n(size = 4000)

original_names <- colnames(samps)
replacement_names <- c("Intercept", 
                       "Site CWD", 
                       "Proportion host (ponderosa)", 
                       "Ponderosa mean height", 
                       "Overall density", 
                       "Overall basal area",
                       "Site CWD : Proportion host (ponderosa)",
                       "Site CWD : Ponderosa mean height",
                       "Site CWD : Overall density",
                       "Site CWD : Overall basal area",
                       "Proportion host (ponderosa) : Overall density",
                       "Proportion host (ponderosa) : Ponderosa mean height",
                       "Site CWD : Proportion host (ponderosa) : Ponderosa mean size")

replacement_name_order <- 1:13

colnames(samps) <- replacement_names

long_samps <-
  samps %>% 
  tidyr::pivot_longer(cols = 1:ncol(.), names_to = "variable", values_to = "samps") %>% 
  dplyr::mutate(variable = factor(variable,
                                  levels = replacement_names[rev(replacement_name_order)]))

effect_sizes_halfeye <-
  ggplot(long_samps, aes(x = samps, y = variable)) +
  geom_halfeyeh() +
  theme_bw(base_size = 9) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Effect size\nLog odds change in Pr(Ponderosa mortality)\nfor a 1 standard deviation increase in covariate", y = NULL)

effect_sizes_halfeye

ggsave(plot = effect_sizes_halfeye, filename = "figures/effect-sizes-halfeye.png", width = 6, height = 3.5, units = "in")



replacement_names <- c("Intercept", 
                       "Site CWD", 
                       "Proportion host (ponderosa)", 
                       "Ponderosa mean height", 
                       "Overall density", 
                       "Overall basal area",
                       "Site CWD : Proportion host (ponderosa)",
                       "Site CWD : Ponderosa mean height",
                       "Site CWD : Overall density",
                       "Site CWD : Overall basal area",
                       "Proportion host (ponderosa) : Overall density",
                       "Proportion host (ponderosa) : Ponderosa mean height",
                       "Ponderosa mean height : Overall basal area",
                       "Site CWD : Proportion host (ponderosa) : Ponderosa mean size")

replacement_name_order <- 1:14

colnames(samps) <- replacement_names

long_samps <-
  samps %>% 
  tidyr::gather(key = "variable", value = "samps") %>% 
  dplyr::mutate(variable = factor(variable,
                                  levels = replacement_names[rev(replacement_name_order)]))

effect_sizes_halfeye <-
  ggplot(long_samps, aes(x = samps, y = variable)) +
  geom_halfeyeh() +
  theme_bw(base_size = 9) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Effect size\nLog odds change in Pr(Ponderosa mortality)\nfor a 1 standard deviation increase in covariate", y = NULL)

effect_sizes_halfeye

ggsave(plot = effect_sizes_halfeye, filename = "figures/effect-sizes-halfeye_uncalibrated.png", width = 6, height = 3.5, units = "in")

