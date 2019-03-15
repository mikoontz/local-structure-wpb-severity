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

set.seed(0314)
adf <-
  analysis_df %>% 
  dplyr::filter(pipo_and_dead_count != 0) %>% 
  dplyr::mutate(site = as.factor(site)) %>% 
  dplyr::group_by(site) %>% 
  dplyr::sample_n(200)


# n = 10: 217 seconds (0.06 hours)
# n = 50: 1454 seconds (0.4 hours)
# n = 100: 1.6 hours
# n = 200 = 8.9 hours

# tic()
# fm1_mgcv <- bam(cbind(dead_count, pipo_count) ~ 
#                   site_cwd_zscore*pipo_and_dead_tpha_s*pipo_and_dead_qmd_s +
#                   site_cwd_zscore*overall_tpha_s*overall_qmd_s +
#                   s(x, y, by = site, bs = "gp", k = 5), 
#                 data = adf, 
#                 family = binomial(link = "logit"))
# toc()
# summary(fm1_mgcv)
# plot(fm1_mgcv)

# tic()
# fm2_mgcv <- bam(cbind(dead_count, pipo_count) ~ 
#                   site_cwd_zscore*pipo_and_dead_tpha_s*pipo_and_dead_qmd_s +
#                   site_cwd_zscore*overall_tpha_s*overall_qmd_s +
#                   s(x, y, by = site, bs = "gp", k = 10, m = c(3, 5)), 
#                 data = adf, 
#                 family = binomial(link = "logit"))
# toc()
# summary(fm2_mgcv)


# 
# tic()
# fm3_mgcv <- bam(cbind(dead_count, pipo_count) ~ 
#                   local_cwd_zscore*pipo_and_dead_tpha_s*pipo_and_dead_qmd_s +
#                   local_cwd_zscore*overall_tpha_s*overall_qmd_s +
#                   s(x, y, by = site, k = 10), 
#                 data = adf, 
#                 family = binomial(link = "logit"))
# toc()
# summary(fm3_mgcv)
# 
# 
# tic()
# fm4_mgcv <- bam(cbind(dead_count, pipo_count) ~ 
#                   site_cwd_zscore*pipo_and_dead_tpha_s*pipo_and_dead_qmd_s +
#                   site_cwd_zscore*overall_tpha_s*overall_qmd_s +
#                   s(x, y, by = site, k = 10), 
#                 data = adf, 
#                 family = binomial(link = "logit"))
# toc()
# summary(fm4_mgcv)

# library(parallel)  
# nc <- 6   ## cluster size, set for example portability
# cl <- makeCluster(nc) 
# tic()
# fm5_mgcv <- bam(cbind(dead_count, pipo_count) ~ 
#                   site_cwd_zscore*pipo_and_dead_tpha_s*pipo_and_dead_qmd_s +
#                   site_cwd_zscore*overall_tpha_s*overall_qmd_s +
#                   te(x, y, by = site), 
#                 data = adf, 
#                 family = binomial(link = "logit"))
# # ,
# #                 cluster = cl)
# toc()
# summary(fm5_mgcv)
# stopCluster(cl)

(start <- Sys.time())
fm1_brms <- brm(dead_count | trials(pipo_and_dead_count) ~ 
                  site_cwd_zscore*pipo_and_dead_tpha_s*pipo_and_dead_qmd_s +
                  site_cwd_zscore*overall_tpha_s*overall_qmd_s +
                  gp(x, y, by = site, k = 3, scale = FALSE),
                data = adf,
                family = binomial(link = "logit"),
                chains = 4,
                cores = 4,
                iter = 2500,
                control = list(adapt_delta = 0.95))
summary(fm1_brms)
(elapsed <- Sys.time() - start)

prior_summary(fm1_brms)

saveRDS(fm1_brms, file = here::here("analyses/analyses_output/fitted-model_site-cwdZscore_pipo-tpha-qmd_overall-tpha-qmd_gp-per-site_200-samples.rds"))

# if(!file.exists(here::here("analyses/analyses_output/fitted-model_cwdZscore_pipo-count-ba_total-count-ba_uniqueCellID.rds"))) {
#   # tic()
#   future::plan(strategy = multiprocess)
#   fm1_brms <- brm(dead_count | trials(live_and_dead_pipo_count) ~ 
#                     cwd_zscore*live_and_dead_pipo_count_s*live_and_dead_pipo_ba_s +
#                     cwd_zscore*total_count_s*total_ba_s + 
#                     (1 | unique_cellID), 
#                   data = analysis_df, 
#                   family = binomial(link = "logit"),
#                   prior=c(set_prior("normal (0, 8)")),
#                   chains = 4,
#                   future = TRUE)
#   summary(fm1_brms)
#   # toc() # 7200 seconds
#   
#   saveRDS(fm1_brms, here::here("analyses/analyses_output/fitted-model_cwdZscore_pipo-count-ba_total-count-ba_uniqueCellID.rds"))
# } else {
#   fm1_brms <- readRDS(here::here("analyses/analyses_output/fitted-model_cwdZscore_pipo-count-ba_total-count-ba_uniqueCellID.rds"))
# }
# 

