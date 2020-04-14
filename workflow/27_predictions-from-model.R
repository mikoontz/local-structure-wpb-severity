# Purpose: Use the fitted model to get model predictions for a range of newdata
# Will be useful for plotting, and "zeroing out" the spatial Gaussian Process 
# in order to focus on the fixed effects of the forest structure/environment
# covariates

library(tidyverse)
library(brms)

cwd_data <- read_csv("data/data_output/cwd-data.csv")

fm1 <- readr::read_rds(path =  here::here('analyses', 'analyses_output', 'fitted-model_zibinomial_site-cwdZscore_prop-host_pipo-height_overall-tpha_overall-bapha_exact-gp-per-site_200-samples.rds'))

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

samps <- 
  posterior_samples(fm1) %>% 
  dplyr::select(tidyselect::all_of(paste0("b_", covariates))) %>% 
  dplyr::sample_n(size = 4000)

(start <- Sys.time())
lwr <- purrr::map_dbl(1:nrow(newdata), .f = function(row) {quantile(plogis(as.matrix(samps) %*% t(newdata[row, ])), prob = 0.025)})
(Sys.time() - start)

(start <- Sys.time())
est_mn <- purrr::map_dbl(1:nrow(newdata), .f = function(row) {mean(plogis(as.matrix(samps) %*% t(newdata[row, ])))})
(Sys.time() - start)

(start <- Sys.time())
upr <- purrr::map_dbl(1:nrow(newdata), .f = function(row) {quantile(plogis(as.matrix(samps) %*% t(newdata[row, ])), prob = 0.975)})
(Sys.time() - start)

fitted <- newdata
fitted$lwr <- lwr
fitted$est_mn <- est_mn
fitted$upr <- upr

fitted_compact <-
  fitted %>% 
  dplyr::select(tidyselect::all_of(main_effect_covariates), 
                lwr, est_mn, upr)

write_csv(fitted_compact, path = "analyses/analyses_output/model-predictions.csv")
