library(tidyverse)
library(brms)

fm1 <- readr::read_rds(path =  here::here('analyses', 'analyses_output', 'fitted-model_zibinomial_site-cwdZscore_prop-host_pipo-height_overall-tpha_overall-bapha_exact-gp-per-site_200-samples.rds'))
pp_check_gg <- 
  brms::pp_check(fm1, nsamples = 50)

pp_check_gg <-
  pp_check_gg +
  labs(x = "Actual number of dead trees in a cell",
       y = "Model prediction of number of dead trees in a cell") +
  theme_bw(base_size = 9)

ggplot2::ggsave(filename = here::here("figures", "posterior-prediction-model-checking.png"), width = 6, height = 3.5, units = "in")
