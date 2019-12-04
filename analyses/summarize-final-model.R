# Purpose: summarize the effect sizes of the final model

library(brms)
library(tidyverse)

fm1 <- readRDS("analyses/analyses_output/fitted-model_zibinomial_site-cwdZscore_prop-host_pipo-height_overall-tpha_exact-gp-per-site_200-samples.rds")

fm1
summary(fm1)

out <- 
  fm1 %>% 
  fixef() %>% 
  as_tibble() %>% 
  mutate(beta = rownames(fixef(fm1))) %>% 
  rename(lwr0_025 = Q2.5, 
         upr0_975 = Q97.5)

write_csv(x = out, path = "analyses/analyses_output/final-model-summary.csv")
