# Purpose: Use the fitted model to get model predictions for a range of newdata
# Will be useful for plotting, and "zeroing out" the spatial Gaussian Process 
# in order to focus on the fixed effects of the forest structure/environment
# covariates

library(tidyverse)
library(brms)

cwd_data <- read_csv("data/data_output/cwd-data.csv")

fm1 <- readRDS("analyses/analyses_output/fitted-model_zibinomial_site-cwdZscore_pipo-tpha-qmd_overall-tpha-qmd_exact-gp-per-site_200-samples.rds")

step_size <- 0.1

covariates <- rownames(fixef(fm1))
interaction_idx <- grep(pattern = ":", x = covariates)
main_effect_covariates <- covariates[-interaction_idx]

interaction_covariates <- 
  tibble(interaction_names = covariates[interaction_idx],
         interaction_maths = str_replace_all(interaction_names, pattern = ":", replacement = "*"))

lapply(main_effect_covariates, FUN = function(x) assign(x, seq(-1, 1, by = step_size), pos = 1))

Intercept <- 1
site_cwd_zscore <- c(min(cwd_data$site_cwd_zscore), 
                     mean(cwd_data$site_cwd_zscore), 
                     max(cwd_data$site_cwd_zscore))

newdata  <-
  expand.grid(lapply(main_effect_covariates, FUN = get)) %>% 
  as_tibble() %>% 
  setNames(nm = main_effect_covariates)

for (i in seq_along(1:nrow(interaction_covariates))) {
 
  interaction_name <- (interaction_covariates[i, "interaction_names"] %>% pull())
  interaction_math <- (interaction_covariates[i, "interaction_maths"] %>% pull())
  
  newdata <-
    newdata %>% 
    dplyr::mutate(!!interaction_name := !!str2lang(interaction_math))
   
}

# newdata <- 
#   expand.grid(Intercept = 1,
#               site_cwd_zscore = -1:1,
#               pipo_and_dead_tpha_s = seq(-1, 1, by = step_size),
#               pipo_and_dead_qmd_s = seq(-1, 1, by = step_size),
#               overall_tpha_s = seq(-1, 1, by = step_size),
#               overall_qmd_s = seq(-1, 1, by = step_size)) %>% 
#   as_tibble() %>% 
#   dplyr::mutate(`site_cwd_zscore:pipo_and_dead_tpha_s` = site_cwd_zscore * pipo_and_dead_tpha_s,
#                 `site_cwd_zscore:pipo_and_dead_qmd_s` = site_cwd_zscore * pipo_and_dead_qmd_s,
#                 `pipo_and_dead_tpha_s:pipo_and_dead_qmd_s` = pipo_and_dead_tpha_s * pipo_and_dead_qmd_s,
#                 `site_cwd_zscore:overall_tpha_s` = site_cwd_zscore * overall_tpha_s,
#                 `site_cwd_zscore:overall_qmd_s` = site_cwd_zscore * overall_qmd_s,
#                 `overall_tpha_s:overall_qmd_s` = overall_tpha_s * overall_qmd_s,
#                 `site_cwd_zscore:pipo_and_dead_tpha_s:pipo_and_dead_qmd_s` = site_cwd_zscore * pipo_and_dead_tpha_s * pipo_and_dead_qmd_s,
#                 `site_cwd_zscore:overall_tpha_s:overall_qmd_s` = site_cwd_zscore * overall_tpha_s * overall_qmd_s) %>% 
#   dplyr::select(Intercept,
#                 site_cwd_zscore,
#                 pipo_and_dead_tpha_s,
#                 pipo_and_dead_qmd_s,
#                 overall_tpha_s,
#                 overall_qmd_s,
#                 `site_cwd_zscore:pipo_and_dead_tpha_s`,
#                 `site_cwd_zscore:pipo_and_dead_qmd_s`,
#                 `pipo_and_dead_tpha_s:pipo_and_dead_qmd_s`,
#                 `site_cwd_zscore:overall_tpha_s`,
#                 `site_cwd_zscore:overall_qmd_s`,
#                 `overall_tpha_s:overall_qmd_s`,
#                 `site_cwd_zscore:pipo_and_dead_tpha_s:pipo_and_dead_qmd_s`,
#                 `site_cwd_zscore:overall_tpha_s:overall_qmd_s`)
# 
# write_csv(newdata, path = "analyses/analyses_output/newdata-for-model-predictions.csv")

samps <- 
  posterior_samples(fm1) %>% 
  dplyr::select(one_of(paste0("b_", covariates)))

# samps <- posterior_samples(fm1)[, 1:14]

(start <- Sys.time())
lwr <- purrr::map_dbl(1:nrow(newdata), .f = function(row) {quantile(plogis(as.matrix(samps) %*% t(newdata[row, ])), prob = 0.025)})
(Sys.time() - start)

# (start <- Sys.time())
# est_med <- purrr::map_dbl(1:nrow(newdata), .f = function(row) {median(plogis(as.matrix(samps) %*% t(newdata[row, ])))})
# (Sys.time() - start)

(start <- Sys.time())
est_mn <- purrr::map_dbl(1:nrow(newdata), .f = function(row) {mean(plogis(as.matrix(samps) %*% t(newdata[row, ])))})
(Sys.time() - start)

# (start <- Sys.time())
# est_sd <- purrr::map_dbl(1:nrow(newdata), .f = function(row) {sd(plogis(as.matrix(samps) %*% t(newdata[row, ])))})
# (Sys.time() - start)

(start <- Sys.time())
upr <- purrr::map_dbl(1:nrow(newdata), .f = function(row) {quantile(plogis(as.matrix(samps) %*% t(newdata[row, ])), prob = 0.975)})
(Sys.time() - start)

fitted <- newdata
fitted$lwr <- lwr
# fitted$est_med <- est_med
fitted$est_mn <- est_mn
fitted$upr <- upr
# fitted$est_sd <- est_sd

object.size(fitted)
head(fitted)

fitted_compact <-
  fitted %>% 
  dplyr::select(one_of(main_effect_covariates), 
                lwr, est_mn, upr)
# fitted_compact <- 
#   fitted %>% 
#   dplyr::select(site_cwd_zscore,
#                 pipo_and_dead_tpha_s,
#                 pipo_and_dead_qmd_s,
#                 overall_tpha_s,
#                 overall_qmd_s,
#                 lwr, est_med, est_mn, upr, est_sd)

object.size(fitted_compact)
write_csv(fitted_compact, path = "analyses/analyses_output/model-predictions.csv")



newdata2 <- 
  expand.grid(Intercept = 1,
              site_cwd_zscore = -1:1,
              pipo_and_dead_tpha_s = seq(-1, 1, by = step_size),
              pipo_and_dead_qmd_s = seq(-1, 1, by = step_size),
              overall_bapha_s = seq(-1, 1, by = step_size),
              non_pipo_tpha_s = seq(-1, 1, by = step_size)) %>% 
  as_tibble() %>% 
  dplyr::mutate(`site_cwd_zscore:pipo_and_dead_tpha_s` = site_cwd_zscore * pipo_and_dead_tpha_s,
                `site_cwd_zscore:pipo_and_dead_qmd_s` = site_cwd_zscore * pipo_and_dead_qmd_s,
                `pipo_and_dead_tpha_s:pipo_and_dead_qmd_s` = pipo_and_dead_tpha_s * pipo_and_dead_qmd_s,
                `site_cwd_zscore:non_pipo_tpha_s` = site_cwd_zscore * non_pipo_tpha_s,
                `site_cwd_zscore:overall_bapha_s` = site_cwd_zscore * overall_bapha_s,
                `non_pipo_tpha_s:overall_bapha_s` = non_pipo_tpha_s * overall_bapha_s,
                `site_cwd_zscore:pipo_and_dead_tpha_s:pipo_and_dead_qmd_s` = site_cwd_zscore * pipo_and_dead_tpha_s * pipo_and_dead_qmd_s,
                `site_cwd_zscore:non_pipo_tpha_s:overall_bapha_s` = site_cwd_zscore * non_pipo_tpha_s * overall_bapha_s) %>% 
  dplyr::select(Intercept,
                site_cwd_zscore,
                pipo_and_dead_tpha_s,
                pipo_and_dead_qmd_s,
                non_pipo_tpha_s,
                overall_bapha_s,
                `site_cwd_zscore:pipo_and_dead_tpha_s`,
                `site_cwd_zscore:pipo_and_dead_qmd_s`,
                `pipo_and_dead_tpha_s:pipo_and_dead_qmd_s`,
                `site_cwd_zscore:non_pipo_tpha_s`,
                `site_cwd_zscore:overall_bapha_s`,
                `non_pipo_tpha_s:overall_bapha_s`,
                `site_cwd_zscore:pipo_and_dead_tpha_s:pipo_and_dead_qmd_s`,
                `site_cwd_zscore:non_pipo_tpha_s:overall_bapha_s`)

write_csv(newdata2, path = "analyses/analyses_output/newdata-for-model-predictions_2.csv")

samps <- posterior_samples(fm1)[, 1:14]

(start <- Sys.time())
lwr <- purrr::map_dbl(1:nrow(newdata2), .f = function(row) {quantile(plogis(as.matrix(samps) %*% t(newdata2[row, ])), prob = 0.025)})
(Sys.time() - start)

(start <- Sys.time())
est_med <- purrr::map_dbl(1:nrow(newdata2), .f = function(row) {median(plogis(as.matrix(samps) %*% t(newdata2[row, ])))})
(Sys.time() - start)

(start <- Sys.time())
est_mn <- purrr::map_dbl(1:nrow(newdata2), .f = function(row) {mean(plogis(as.matrix(samps) %*% t(newdata2[row, ])))})
(Sys.time() - start)

(start <- Sys.time())
est_sd <- purrr::map_dbl(1:nrow(newdata2), .f = function(row) {sd(plogis(as.matrix(samps) %*% t(newdata2[row, ])))})
(Sys.time() - start)

(start <- Sys.time())
upr <- purrr::map_dbl(1:nrow(newdata2), .f = function(row) {quantile(plogis(as.matrix(samps) %*% t(newdata2[row, ])), prob = 0.975)})
(Sys.time() - start)

fitted <- newdata2
fitted$lwr <- lwr
fitted$est_med <- est_med
fitted$est_mn <- est_mn
fitted$upr <- upr
fitted$est_sd <- est_sd

object.size(fitted)
head(fitted)

fitted_compact <- 
  fitted %>% 
  dplyr::select(site_cwd_zscore,
                pipo_and_dead_tpha_s,
                pipo_and_dead_qmd_s,
                overall_bapha_s,
                non_pipo_tpha_s,
                lwr, est_med, est_mn, upr, est_sd)

object.size(fitted_compact)
write_csv(fitted_compact, path = "analyses/analyses_output/model-predictions_2.csv")
