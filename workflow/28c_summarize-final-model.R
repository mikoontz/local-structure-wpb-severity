# Purpose: summarize the effect sizes of the final model

library(brms)
library(tidyverse)

fm1 <- readr::read_rds(path =  here::here('analyses', 'analyses_output', 'fitted-model_zibinomial_site-cwdZscore_prop-host_pipo-height_overall-tpha_overall-bapha_exact-gp-per-site_200-samples.rds'))

fm1
summary(fm1)

out <- 
  fm1 %>% 
  brms::posterior_summary() %>% 
  as_tibble() %>% 
  mutate(beta = rownames(posterior_summary(fm1))) %>% 
  rename(lwr0_025 = Q2.5, 
         upr0_975 = Q97.5) %>% 
  dplyr::mutate(beta = ifelse(stringi::stri_startswith_fixed(str = beta, pattern = "b_"), yes = stringi::stri_sub(beta, from = 3, to = -1), no = beta) )

write_csv(x = out, path = here::here("analyses", "analyses_output", "final-model-summary.csv"))