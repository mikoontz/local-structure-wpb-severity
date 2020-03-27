# Purpose: derive allometric equations to convert height to basal area for the study region

library(tidyverse)
library(mgcv)
library(here)

if(!file.exists("data/data_output/formatted-ground-data.csv")) {
 source(here::here("workflow/11_format-ground-data.R"))
}

d <- readr::read_csv("data/data_output/formatted-ground-data.csv")

key_species <- c("PIPO", "PILA", "CADE", "QUKE", "ABCO")

dd <- 
  d %>% 
  dplyr::filter(species %in% key_species) %>% 
  dplyr::filter(live == 1) %>% 
  dplyr::mutate(species = tolower(species),
                ba = pi * (dbh / 2)^2,
                height_squared = height^2) 
  
allometry_models <-
  dd %>% 
  dplyr::group_by(species) %>% 
  dplyr::do(model = lm(dbh ~ height, data = .))

