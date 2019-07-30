# Purpose: visualize the results of the primary analysis
# Triptych of low, average, and high CWD in different facets
# pondo qmd vs. pondo density with color of point representing probability of
# pondo mortality
# all tree qmd vs. all tree density with color of point representing Pr(pipo mortality)

library(tidyverse)
library(brms)
library(tidybayes)

fm1  <- readRDS("analyses/analyses_output/fitted-model_zibinomial_site-cwdZscore_pipo-tpha-qmd_overall-tpha-qmd_exact-gp-per-site_200-samples.rds")
samps <- posterior_samples(fm1)[, 1:14]

newdata <- read_csv("analyses/analyses_output/newdata-for-model-predictions.csv")
fitted_compact <- read_csv("analyses/analyses_output/model-predictions.csv")
newdata <- read.csv("analyses/analyses_output/newdata-for-model-predictions_2.csv", stringsAsFactors = FALSE)
fitted_compact <- read.csv("analyses/analyses_output/model-predictions_2.csv", stringsAsFactors = FALSE)

step_size <- 0.1

fitted_just_estimates <- newdata
fitted_just_estimates$prob_pipo_mortality <- plogis(as.matrix(newdata) %*% fixef(fm1)[, "Estimate"])

fitted_subset <-
  fitted_just_estimates %>% 
  dplyr::filter(overall_tpha_s == 0,
                overall_qmd_s == 0) %>% 
  dplyr::mutate(cwd = case_when(site_cwd_zscore == -1 ~ "cool/wet site",
                                site_cwd_zscore == 0 ~ "average site",
                                site_cwd_zscore == 1 ~ "hot/dry site")) %>% 
  dplyr::mutate(cwd = factor(cwd, levels = c("cool/wet site", "average site", "hot/dry site")))

pipo_tpha_qmd_cwd_interaction_gg_raster <-
  ggplot(fitted_subset, aes(x = pipo_and_dead_tpha_s, y = pipo_and_dead_qmd_s, fill = prob_pipo_mortality)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(~ cwd) +
  theme_bw() +
  labs(x = "Ponderosa pine density\n(trees per hectare; scaled)",
       y = "Ponderosa pine quadratic mean diameter\n(cm; scaled)",
       fill = "Probability of mortality\nof ponderosa pine")

pipo_tpha_qmd_cwd_interaction_df <-
  fitted_compact %>% 
  dplyr::filter(site_cwd_zscore %in% c(-1, 0, 1)) %>% 
  dplyr::filter(overall_tpha_s == 0) %>% 
  dplyr::filter(overall_qmd_s == 0) %>% 
  dplyr::filter(abs(abs(pipo_and_dead_qmd_s) - 0.7) < 0.05) %>% 
  dplyr::mutate(pipo_and_dead_qmd = ifelse(pipo_and_dead_qmd_s == -0.7, yes = "Smaller trees", no = "Larger trees")) %>% 
  dplyr::mutate(cwd = case_when(site_cwd_zscore == -1 ~ "cool/wet site",
                                site_cwd_zscore == 0 ~ "average site",
                                site_cwd_zscore == 1 ~ "hot/dry site")) %>% 
  dplyr::mutate(cwd = factor(cwd, levels = c("cool/wet site", "average site", "hot/dry site")))

pipo_tpha_qmd_cwd_interaction_gg <-
  ggplot(pipo_tpha_qmd_cwd_interaction_df, 
         aes(x = pipo_and_dead_tpha_s, 
             y = est_mn, 
             color = pipo_and_dead_qmd)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, group = as.factor(pipo_and_dead_qmd_s)), fill = "lightgray", color = NA) +
  geom_line() +
  facet_grid(~ cwd) +
  theme_bw() +
  labs(x = "Ponderosa pine density\n(trees per hectare; scaled)",
       y = "Pr(ponderosa mortality)",
       color = "Quadratic\nmean diameter")

ggsave(plot = pipo_tpha_qmd_cwd_interaction_gg, filename = "figures/pipo_tpha_qmd_cwd_interaction.png", width = 6, height = 3.5, units = "in")

ggsave(plot = pipo_tpha_qmd_cwd_interaction_gg_raster, filename = "figures/pipo_tpha_qmd_cwd_interaction_raster.png")


# halfeye plots of model coefficients -------------------------------------
original_names <- colnames(samps)
replacement_names <- c("Intercept", 
                       "Site CWD", 
                       "Ponderosa density", 
                       "Ponderosa mean size", 
                       "Overall density", 
                       "Overall mean size", 
                       "Site CWD : Ponderosa density",
                       "Site CWD : Ponderosa mean size",
                       "Ponderosa density : Ponderosa mean size",
                       "Site CWD : Overall density",
                       "Site CWD : Overall mean size",
                       "Overall density : Overall mean size",
                       "Site CWD : Ponderosa density : Ponderosa mean size",
                       "Site CWD : Overall density : Overall mean size")

replacement_name_order <- c(1, 2, 3, 5, 4, 6, 9, 12, 7, 10, 8, 11, 13, 14)
colnames(samps) <- replacement_names
long_samps <-
  samps %>% 
  tidyr::gather(key = "variable", value = "samps") %>% 
  dplyr::mutate(variable = factor(variable,
                                  levels = replacement_names[rev(replacement_name_order)]))

effect_sizes_halfeye <-
  ggplot(long_samps, aes(x = samps, y = variable)) +
  geom_halfeyeh() +
  theme_bw(base_size = 10) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Effect size\nLog odds change in Pr(Ponderosa mortality)\nfor a 1 standard deviation increase in covariate", y = NULL)

ggsave(plot = effect_sizes_halfeye, filename = "figures/effect-sizes-halfeye.png", width = 6, height = 3.5, units = "in")



# newdata2 ----------------------------------------------------------------

newdata <- read.csv("analyses/analyses_output/newdata-for-model-predictions_2.csv", stringsAsFactors = FALSE)
fitted_compact <- read.csv("analyses/analyses_output/model-predictions_2.csv", stringsAsFactors = FALSE)

step_size <- 0.1

fitted_just_estimates <- newdata
fitted_just_estimates$prob_pipo_mortality <- plogis(as.matrix(newdata) %*% fixef(fm1)[, "Estimate"])

fitted_subset <-
  fitted_just_estimates %>% 
  dplyr::filter(non_pipo_tpha_s == 0,
                overall_bapha_s == 0) %>% 
  dplyr::mutate(cwd = case_when(site_cwd_zscore == -1 ~ "cool/wet site",
                                site_cwd_zscore == 0 ~ "average site",
                                site_cwd_zscore == 1 ~ "hot/dry site")) %>% 
  dplyr::mutate(cwd = factor(cwd, levels = c("cool/wet site", "average site", "hot/dry site")))

pipo_tpha_qmd_cwd_interaction_gg_raster <-
  ggplot(fitted_subset, aes(x = pipo_and_dead_tpha_s, y = pipo_and_dead_qmd_s, fill = prob_pipo_mortality)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(~ cwd) +
  theme_bw() +
  labs(x = "Ponderosa pine density\n(trees per hectare; scaled)",
       y = "Ponderosa pine quadratic mean diameter\n(cm; scaled)",
       fill = "Probability of mortality\nof ponderosa pine")

pipo_tpha_qmd_cwd_interaction_df <-
  fitted_compact %>% 
  dplyr::filter(site_cwd_zscore %in% c(-1, 0, 1)) %>% 
  dplyr::filter(non_pipo_tpha_s == 0) %>% 
  dplyr::filter(overall_bapha_s == 0) %>% 
  dplyr::filter(abs(abs(pipo_and_dead_qmd_s) - 0.7) < 0.05) %>% 
  dplyr::mutate(pipo_and_dead_qmd = ifelse(pipo_and_dead_qmd_s == -0.7, yes = "Smaller trees", no = "Larger trees")) %>% 
  dplyr::mutate(cwd = case_when(site_cwd_zscore == -1 ~ "cool/wet site",
                                site_cwd_zscore == 0 ~ "average site",
                                site_cwd_zscore == 1 ~ "hot/dry site")) %>% 
  dplyr::mutate(cwd = factor(cwd, levels = c("cool/wet site", "average site", "hot/dry site")))

pipo_tpha_qmd_cwd_interaction_gg <-
  ggplot(pipo_tpha_qmd_cwd_interaction_df, 
         aes(x = pipo_and_dead_tpha_s, 
             y = est_mn, 
             color = pipo_and_dead_qmd)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, group = as.factor(pipo_and_dead_qmd_s)), fill = "lightgray", color = NA) +
  geom_line() +
  facet_grid(~ cwd) +
  theme_bw() +
  labs(x = "Ponderosa pine density\n(trees per hectare; scaled)",
       y = "Pr(ponderosa mortality)",
       color = "Quadratic\nmean diameter")
