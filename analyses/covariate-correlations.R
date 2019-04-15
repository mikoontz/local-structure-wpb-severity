# Purpose: show correlations between covariates

library(tidyverse)
library(GGally)

data <- read_csv("analyses/analyses_output/data-from-rasterized-classified-trees.csv")

small_data <-
  data %>% 
  dplyr::select(pipo_and_dead_qmd, overall_qmd, pipo_and_dead_tpha, overall_tpha)

ggcorr(small_data, label = TRUE)

ggplot(small_data, aes(pipo_and_dead_qmd, overall_qmd)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red")
