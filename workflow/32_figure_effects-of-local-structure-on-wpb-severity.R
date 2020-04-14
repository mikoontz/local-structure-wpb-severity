# Purpose: visualize the results of the primary analysis
# Triptych of low, average, and high CWD in different facets
# pondo qmd vs. pondo density with color of point representing probability of
# pondo mortality
# all tree qmd vs. all tree density with color of point representing Pr(pipo mortality)

library(tidyverse)
library(brms)
library(tidybayes)

fm1 <- readr::read_rds(path =  here::here('analyses', 'analyses_output', 'fitted-model_zibinomial_site-cwdZscore_prop-host_pipo-height_overall-tpha_overall-bapha_exact-gp-per-site_200-samples.rds'))

fitted_compact <- read_csv("analyses/analyses_output/model-predictions.csv")

step_size <- 0.1
tree_size_small <- -0.7
tree_size_large <- 0.7

prop_host_height_cwd_interaction_df <-
  fitted_compact %>% 
  dplyr::filter(overall_tpha_s == 0 & overall_bapha_s == 0) %>% 
  dplyr::mutate(pipo_and_dead_mean_height_s = round(pipo_and_dead_mean_height_s, 1)) %>% 
  dplyr::filter(pipo_and_dead_mean_height_s %in% c(tree_size_small, tree_size_large)) %>% 
  dplyr::mutate(pipo_and_dead_mean_height = ifelse(pipo_and_dead_mean_height_s == tree_size_small, yes = "Smaller trees", no = "Larger trees")) %>% 
  dplyr::mutate(cwd = case_when(site_cwd_zscore == min(site_cwd_zscore) ~ "cool/wet site",
                                site_cwd_zscore == max(site_cwd_zscore) ~ "hot/dry site",
                                TRUE ~ "average site")) %>% 
  dplyr::mutate(cwd = factor(cwd, levels = c("cool/wet site", "average site", "hot/dry site")))

prop_host_height_cwd_interaction_gg <-
  ggplot(prop_host_height_cwd_interaction_df, 
         aes(x = prop_host_s, 
             y = est_mn, 
             color = pipo_and_dead_mean_height)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, group = pipo_and_dead_mean_height, fill = pipo_and_dead_mean_height), alpha = 0.5, color = NA) +
  geom_line() +
  facet_grid(~ cwd) +
  theme_bw() +
  labs(x = "Proportion host trees (Ponderosa pine); scaled",
       y = "Pr(ponderosa mortality)",
       color = "Mean height",
       fill = "Mean height") +
  # scale_y_continuous(limits = c(0, 1)) +
  scale_fill_viridis_d(option = "E") +
  scale_color_viridis_d(option = "E")

prop_host_height_cwd_interaction_gg

ggsave(filename = "figures/prop-host_pipo-height_cwd_interaction.png", width = 6, height = 3.5, units = "in", plot = prop_host_height_cwd_interaction_gg)

# halfeye plots of model coefficients -------------------------------------
samps <- readr::read_csv(file = here::here("analyses", "analyses_output", "final-model-posterior-samples.csv"))
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

ggsave(plot = effect_sizes_halfeye, filename = "figures/effect-sizes-halfeye.png", width = 6, height = 3.5, units = "in")
