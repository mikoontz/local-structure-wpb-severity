# Purpose: get classification details including number of hand classified trees and accuracy/kappa of the live/dead as well as species classifications

library(tidyverse)
library(caret)

crowns_with_reflectance <- readr::read_csv(here::here("data/data_output/classified/hand-classified/hand-classified-crowns.csv"))

# create the live/dead classifier -----------------------------------------
live_or_dead_idx <- caret::createDataPartition(crowns_with_reflectance$live, p = 0.8, list = FALSE)

live_or_dead_training <- crowns_with_reflectance[live_or_dead_idx, ]
live_or_dead_testing <- crowns_with_reflectance[-live_or_dead_idx, ]

live_or_dead_fit <- train(
  as.factor(live) ~ b_mean + g_mean + r_mean + re_mean + nir_mean + ndvi_mean + rgi_mean + cire_mean + cig_mean + ndre_mean,
  data = live_or_dead_training,
  method = "LogitBoost"
)

# Classifier is called `live_or_dead_fit`
live_or_dead_fit


# subset to only the live crowns ------------------------------------------

live_crowns <- 
  crowns_with_reflectance %>% 
  dplyr::filter(live == 1) %>% 
  dplyr::mutate(functional_group = case_when(species == "pila" ~ "pinus",
                                             species == "pipo" ~ "pinus",
                                             TRUE ~ species))


# create the species classifier -------------------------------------------

species_idx <-
  live_crowns %>%
  pull(species) %>% 
  caret::createDataPartition(p = 0.8, list = FALSE)

species_training <- live_crowns[species_idx, ]
species_testing <- live_crowns[-species_idx, ]

# Best classification method after tests (see docs/interim-reports/)
rdaFit <- train(
  species ~ b_mean + g_mean + r_mean + re_mean + nir_mean + ndvi_mean + rgi_mean + cire_mean + cig_mean + ndre_mean,
  data = species_training,
  method = "rda",
  preProcess = c("center", "scale"),
  tuneGrid = expand.grid(gamma = 0, lambda = c(seq(0.2, 1.0, by = 0.1)))
)

# species classifier is called rdaFit
rdaFit