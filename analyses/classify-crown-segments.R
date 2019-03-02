library(sf)
library(raster)
library(tabularaster)
library(tidyverse)
library(viridis)
library(purrr)
library(velox)
library(viridis)

# Source in the classifier code.
# The live/dead classifier is an R object called `live_or_dead_fit`
# The species classifier is an R object called `rdafit`

source(here::here("analyses/build-crown-classifier.R"))
live_or_dead_fit
rdaFit

# crowns_with_reflectance <- readr::read_csv(file = here::here("data/data_output/classified/model-classified/crowns-with-reflectance_all.csv"))
crowns_with_reflectance <- readr::read_csv(file = here::here("data/data_output/classified/model-classified/crowns-with-reflectance_35m-buffer.csv"))

classified_trees  <-
  crowns_with_reflectance %>% 
  dplyr::mutate(live = live_or_dead_fit$levels[predict(live_or_dead_fit, newdata = .)]) %>% 
  dplyr::mutate(live = as.numeric(as.character(live))) %>%
  dplyr::mutate(species = ifelse(live == 1, yes = rdaFit$levels[predict(rdaFit, newdata = .)], no = NA))

readr::write_csv(classified_trees, path = here::here("analyses/analyses_output/classified-trees.csv"))
