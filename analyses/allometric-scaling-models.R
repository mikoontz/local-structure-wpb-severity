# Purpose: derive allometric equations to convert height to basal area for the study region

library(tidyverse)
library(here)

source(here::here("data/data_carpentry/format-ground-data.R"))

glimpse(d) 

key_species <- c("PIPO", "PILA", "CADE", "QUKE", "ABCO")

dd <- 
  d %>% 
  dplyr::filter(species %in% key_species) %>% 
  dplyr::filter(live == 1)

# Seems like a linear model predicting dbh from height is a reasonable fit
ggplot(dd, aes(x = height, y = dbh)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ species)

allometry_models <-
  dd %>% 
  dplyr::mutate(species = tolower(species)) %>% 
  dplyr::group_by(species) %>% 
  do(model = lm(dbh ~ height, data = .))

