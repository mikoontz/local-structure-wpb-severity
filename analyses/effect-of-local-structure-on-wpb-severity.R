# Purpose: build the primary paper analysis

library(tidyverse)
library(raster)
library(sf)
library(here)
library(tictoc)
library(brms)
library(lme4)
library(effects)
library(future)


if(file.exists(here::here("analyses/analyses_output/data-from-rasterized-classified-trees.csv"))) {
  
  final_results <- 
    readr::read_csv(here::here("analyses/analyses_output/data-from-rasterized-classified-trees.csv"))
} else {
  stop("You need to extract the data from the rasterized version of the classified trees! See the analyses/rasterize-classified-trees.R script.")
}

glimpse(final_results)

# Get the CWD data
source(here::here("data/data_carpentry/extract-cwd-from-locations.R"))
# R object is called `cwd_data`
cwd_data

analysis_df <-
  final_results %>% 
  dplyr::left_join(cwd_data, by = "site") %>% 
  as_tibble() %>% 
  dplyr::mutate(live_and_dead_pipo_count = pipo_count + dead_count,
                live_and_dead_pipo_ba = pipo_ba + dead_ba,
                live_and_dead_pipo_mean_ba = pipo_mean_ba + dead_mean_ba) %>% 
  dplyr::mutate(live_and_dead_pipo_count_s = scale(live_and_dead_pipo_count, center = TRUE, scale = FALSE),
                non_pipo_count_s = scale(non_pipo_count, center = TRUE, scale = FALSE),
                total_count_s = scale(total_count, center = TRUE, scale = FALSE),
                live_and_dead_pipo_ba_s = scale(live_and_dead_pipo_ba, center = TRUE, scale = FALSE),
                non_pipo_ba_s = scale(non_pipo_ba, center = TRUE, scale = FALSE),
                total_ba_s = scale(total_ba, center = TRUE, scale = FALSE),
                live_and_dead_pipo_mean_ba_s = scale(live_and_dead_pipo_mean_ba, center = TRUE, scale = FALSE)) %>% 
  dplyr::mutate(unique_cellID = 1:nrow(.))
  
  
  
  
  
  # summarized_df <-
  #   analysis_df %>% 
  #   group_by(site) %>% 
  #   summarize(live_count = sum(live_count),
  #             dead_count = sum(dead_count),
  #             pipo_count = sum(pipo_count),
#             non_pipo_count = sum(non_pipo_count),
#             total_count = sum(total_count),
#             pipo_ba = sum(pipo_ba),
#             non_pipo_ba = sum(non_pipo_ba),
#             total_ba = sum(total_ba),
#             cwd_zscore = mean(cwd_zscore))
# 
# fm0a <- glm(cbind(dead_count, live_count) ~ cwd_zscore*total_ba, data = analysis_df, family = "binomial")
# summary(fm0a)
# 
# fm0b <- glm(cbind(dead_count, pipo_count) ~ cwd_zscore*total_ba, data = analysis_df, family = "binomial")
# summary(fm0b)
# 
# fm0c <- glm(cbind(dead_count, pipo_count) ~ cwd_zscore*live_and_dead_pipo_ba, data = analysis_df, family = "binomial")
# summary(fm0c)
# 
# fm1 <- glm(cbind(dead_count, pipo_count) ~ cwd_zscore*pipo_count*pipo_ba + cwd_zscore*total_count*total_ba, data = analysis_df, family = "binomial")
# summary(fm1)
# analysis_df$pipo_count - analysis_df$dead_count
# analysis_df$total_count - (analysis_df$dead_count)
# 
# data_frame(beta = names(coef(fm1)), estimate = coef(fm1))
# 
# e <- Effect(c("cwd_zscore", "pipo_count", "pipo_ba", "total_count", "total_ba"), fm1, xlevels = list(cwd_zscore = c(-1, 0, 1)))
# e_gg <- data.frame(e)
# 
# ggplot(e_gg %>% filter(total_count == 20 & total_ba == 9.8), aes(x = pipo_count, y = fit, color = as.factor(pipo_ba))) +
#   geom_line() +
#   facet_wrap(~ cwd_zscore) +
#   scale_color_viridis_d()
# 
# fm2 <- glm(cbind(dead_count, pipo_count - dead_count) ~ cwd_zscore*pipo_ba*total_ba, data = analysis_df, family = "binomial")
# 
# fm3 <- glmer(cbind(dead_count, pipo_count - dead_count) ~ cwd_zscore*pipo_ba*non_pipo_ba + (1 | site), data = analysis_df, family = "binomial")
# summary(fm3)
# 
# future::plan(strategy = multiprocess)
# fm2_brms <- brm(dead_count | trials(pipo_count + dead_count) ~ cwd_zscore*pipo_ba*total_ba + (1 | site), 
#                 data = analysis_df, 
#                 family = binomial(link = "logit"),
#                 prior=c(set_prior("normal (0, 8)")),
#                 chains = 3,
#                 cores = 3)
# 
# summary(fm2_brms)
# 
# future::plan(strategy = multiprocess)
# fm3_brms <- brm(dead_count | trials(pipo_count + dead_count) ~ cwd_zscore*pipo_ba*total_ba + (1 | unique_cellID), 
#                 data = analysis_df, 
#                 family = binomial(link = "logit"),
#                 prior=c(set_prior("normal (0, 8)")),
#                 chains = 4,
#                 future = TRUE)
# 
# summary(fm3_brms)
# saveRDS(fm3_brms, here::here("analyses/analyses_output/fitted-model_cwd-zscore_pipo-ba_total-ba_uniqueCellID.rds"))
# 
# 
# future::plan(strategy = multiprocess)
# fm4_brms <- brm(dead_count | trials(pipo_count + dead_count) ~ cwd_zscore*pipo_count*pipo_ba + cwd_zscore*non_pipo_count*non_pipo_ba + (1 | unique_cellID), 
#                 data = analysis_df, 
#                 family = binomial(link = "logit"),
#                 prior=c(set_prior("normal (0, 8)")),
#                 chains = 4,
#                 future = TRUE)
# summary(fm4_brms)


tic()
future::plan(strategy = multiprocess)
fm5_brms <- brm(dead_count | trials(live_and_dead_pipo_count) ~ cwd_zscore*live_and_dead_pipo_count_s*live_and_dead_pipo_ba_s + cwd_zscore*non_pipo_count_s*non_pipo_ba_s + (1 | unique_cellID), 
                data = analysis_df, 
                family = binomial(link = "logit"),
                prior=c(set_prior("normal (0, 8)")),
                chains = 4,
                future = TRUE)
summary(fm5_brms)
toc()

tic()
future::plan(strategy = multiprocess)
fm6_brms <- brm(dead_count | trials(live_and_dead_pipo_count) ~ cwd_zscore*live_and_dead_pipo_count_s*live_and_dead_pipo_ba_s + cwd_zscore*total_count_s*total_ba_s + (1 | unique_cellID), 
                data = analysis_df, 
                family = binomial(link = "logit"),
                prior=c(set_prior("normal (0, 8)")),
                chains = 4,
                future = TRUE)
summary(fm6_brms)
toc()

saveRDS(fm6_brms, here::here("analyses/analyses_output/fitted-model_cwdZscore_pipo-count-ba_total-count-ba_uniqueCellID.rds"))
pp_check(fm6_brms)
pp_check(fm6_brms, type = "error_hist", nsamples = 11)
pp_check(fm6_brms, type = "scatter_avg", nsamples = 100)
pp_check(fm6_brms, type = "stat_2d")
pp_check(fm6_brms, type = "rootogram")
pp_check(fm6_brms, type = "loo_pit")


library(GGally)
ggpairs(analysis_df %>% dplyr::select(cwd_zscore, pipo_count, pipo_ba))
coef(fm4_brms)
ggpairs(analysis_df %>% dplyr::select(cwd_zscore, non_pipo_count, non_pipo_ba))
