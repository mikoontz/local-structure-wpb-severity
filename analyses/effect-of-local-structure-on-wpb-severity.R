# Purpose: build the primary paper analysis

library(tidyverse)
library(raster)
library(sf)
library(here)
library(tictoc)
library(brms)
library(lme4)
library(effects)


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
  dplyr::mutate(total_count = pipo_count + non_pipo_count) %>% 
  dplyr::mutate(total_ba = pipo_ba + non_pipo_ba) %>% 
  dplyr::mutate(unique_cellID = 1:nrow(.))

summarized_df <-
  analysis_df %>% 
  group_by(site) %>% 
  summarize(live_count = sum(live_count),
            dead_count = sum(dead_count),
            pipo_count = sum(pipo_count),
            non_pipo_count = sum(non_pipo_count),
            total_count = sum(total_count),
            pipo_ba = sum(pipo_ba),
            non_pipo_ba = sum(non_pipo_ba),
            total_ba = sum(total_ba),
            cwd_zscore = mean(cwd_zscore))

fm0 <- glm(cbind(dead_count, pipo_count) ~ cwd_zscore*pipo_count*pipo_ba + cwd_zscore*total_count*total_ba, data = summarized_df, family = "binomial")
summary(fm0)

fm1 <- glm(cbind(dead_count, pipo_count) ~ cwd_zscore*pipo_count*pipo_ba + cwd_zscore*total_count*total_ba, data = analysis_df, family = "binomial")
summary(fm1)

data_frame(beta = names(coef(fm1)), estimate = coef(fm1))

e <- Effect(c("cwd_zscore", "pipo_count", "pipo_ba", "total_count", "total_ba"), fm1, xlevels = list(cwd_zscore = c(-1, 0, 1)))
e_gg <- data.frame(e)

ggplot(e_gg %>% filter(total_count == 20 & total_ba == 9.8), aes(x = pipo_count, y = fit, color = as.factor(pipo_ba))) +
  geom_line() +
  facet_wrap(~ cwd_zscore) +
  scale_color_viridis_d()

fm2 <- glm(cbind(dead_count, pipo_count) ~ cwd_zscore*pipo_ba*total_ba, data = analysis_df, family = "binomial")

summary(fm2)

fm2_brms <- brm(dead_count | trials(dead_count + pipo_count) ~ cwd_zscore*pipo_ba*total_ba, data = analysis_df, family = binomial(link = "logit"))
