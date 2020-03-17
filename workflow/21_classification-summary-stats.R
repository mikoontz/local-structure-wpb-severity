# Purpose: get classification details including number of hand classified trees and accuracy/kappa of the live/dead as well as species classifications

library(tidyverse)
library(caret)

crowns_with_reflectance <- readr::read_csv(here::here("data/data_output/classified/hand-classified/hand-classified-crowns.csv"))
crowns_with_reflectance <-
  crowns_with_reflectance %>% 
  dplyr::mutate(live = ifelse(live == 1, yes = "live", no = "dead"))

# create the live/dead classifier -----------------------------------------
live_or_dead_idx <- caret::createDataPartition(crowns_with_reflectance$live, p = 0.8, list = FALSE)

live_or_dead_training <- crowns_with_reflectance[live_or_dead_idx, ]
live_or_dead_testing <- crowns_with_reflectance[-live_or_dead_idx, ]

live_or_dead_fit <- train(
  as.factor(live) ~ b_mean + g_mean + r_mean + re_mean + nir_mean + ndvi_mean + rgi_mean + cire_mean + cig_mean + ndre_mean,
  data = live_or_dead_training,
  method = "LogitBoost",
  trControl = trainControl(classProbs = TRUE, 
                            summaryFunction = twoClassSummary),
  metric = "ROC"
)

# Classifier is called `live_or_dead_fit`
live_or_dead_fit

live_or_dead_predict <- predict(live_or_dead_fit, newdata = live_or_dead_testing)

live_or_dead_pred_accuracy <- 
  data.frame(data = live_or_dead_testing$live, pred = live_or_dead_predict) %>% 
  dplyr::mutate(match = data == pred)

mean(live_or_dead_pred_accuracy$match)

# subset to only the live crowns ------------------------------------------

live_crowns <- 
  crowns_with_reflectance %>% 
  dplyr::filter(live == "live") %>% 
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

species_predict <- predict(rdaFit, newdata = species_testing)

species_pred_accuracy <- 
  data.frame(data = species_testing$species, pred = species_predict) %>% 
  dplyr::mutate(match = data == pred) %>% 
  dplyr::mutate(data_host = ifelse(data == "pipo", yes = "host", no = "non-host"),
                pred_host = ifelse(pred == "pipo", yes = "host", no = "non-host")) %>% 
  dplyr::mutate(match_host = data_host == pred_host)

mean(species_pred_accuracy$match)
mean(species_pred_accuracy$match_host)

classification_summary_stats <-
  tibble(type = c("live/dead", "species", "host/non-host"),
         accuracy = c(mean(live_or_dead_pred_accuracy$match), mean(species_pred_accuracy$match), mean(species_pred_accuracy$match_host)))

write_csv(classification_summary_stats, path = "analyses/analyses_output/classification-summary-stats.csv")
