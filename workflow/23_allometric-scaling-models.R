# Purpose: derive allometric equations to convert height to basal area for the study region

library(tidyverse)
library(here)

source(here::here("workflow/11_format-ground-data.R"))

glimpse(d) 

key_species <- c("PIPO", "PILA", "CADE", "QUKE", "ABCO")

dd <- 
  d %>% 
  dplyr::filter(species %in% key_species) %>% 
  dplyr::filter(live == 1) %>% 
  dplyr::mutate(species = tolower(species),
                ba = pi * (dbh / 2)^2,
                height_squared = height^2) 
  


fm_ba <- gam(ba ~ height*species, data = dd)
fm_dbh <- gam(dbh ~ height*species, data = dd)

allometry_models <-
  dd %>% 
  dplyr::mutate(predicted_ba = predict(fm_ba),
                predicted_dbh = predict(fm_dbh),
                predicted_ba_from_dbh = pi * (predicted_dbh / 2)^2)

ggplot(allometry_models, aes(x = ba, y = predicted_ba, color = species)) +
  geom_point() +
  geom_smooth() +
  geom_abline(slope = 1, intercept = 0)

ggplot(allometry_models, aes(x = ba, y = predicted_ba_from_dbh, color = species)) +
  geom_point() +
  geom_smooth() +
  geom_abline(slope = 1, intercept = 0)

ggplot(allometry_models, aes(x = dbh, y = predicted_dbh, color = species)) +
  geom_point() +
  geom_smooth() +
  geom_abline(slope = 1, intercept = 0)


plot(allometry_models$height, allometry_models$dbh, pch = 19, col = factor(allometry_models$species))
points(allometry_models$height, allometry_models$predicted_dbh, col = factor(allometry_models$species))

grouped_dd <-
  dd %>% 
  dplyr::group_by(species, plot) %>% 
  dplyr::summarize(total_height = sum(height),
                   total_height_squared = sum(height_squared),
                   total_ba = sum(ba),
                   mean_ba = mean(ba),
                   mean_height = mean(height))

# Seems like a linear model predicting dbh from height is a reasonable fit
ggplot(dd, aes(x = height, y = dbh)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ species)

ggplot(dd, aes(x = height, y = dbh)) +
  geom_point() +
  geom_smooth(method = "auto") +
  facet_wrap(~ species)

ggplot(dd, aes(x = height, y = ba)) +
  geom_point() +
  geom_smooth(method = "auto") +
  facet_wrap(~ species)

ggplot(dd, aes(x = height_squared, y = ba)) +
  geom_point() +
  geom_smooth(method = "auto") +
  facet_wrap(~ species)

ggplot(grouped_dd, aes(x = total_height, y = total_ba)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ species)

ggplot(grouped_dd, aes(x = total_height_squared, y = total_ba)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ species)

ggplot(grouped_dd, aes(x = mean_height, y = total_ba)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ species)



