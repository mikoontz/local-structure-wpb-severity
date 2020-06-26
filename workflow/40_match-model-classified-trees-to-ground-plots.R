# Size distribution of our classified trees

library(tidyverse)
library(sf)
library(here)

cwd <-
  read_csv(here::here("data", "data_output", "cwd-data.csv"))

cwd_plot <-
  read_csv(here::here("data", "data_output", "cwd-plot-data.csv"))

#### Create circle polygons representing the plot boundaries
# Can only really use the plot centers identfiable from the air to do this
# since the geolocation accuracy of the plots that can't be tied directly
# to the location of the orange cloth "X" in the orthophotos is not good
# enough
plot_radius <- sqrt((66*66) / pi) * (12 * 2.54 / 100) 
plot_locations <- 
  sf::st_read(here::here("data", "data_drone", "L1", "plot-centers-identifiable-from-air_3310.gpkg")) %>% 
  sf::st_buffer(plot_radius) %>% 
  dplyr::select(-site, -local_x, -local_y, -local_crs)

#### Get ground trees
# Subset ground trees to the 110 plots that are identifiable by air
ground_trees <- sf::st_read(here::here(paste0("data/data_drone/L1/ground-trees.gpkg")))

ground_trees <-
  read_csv(here::here("data", "data_output", "formatted-ground-data.csv")) %>%
  filter(is.na(year_fall)) %>%
  left_join(cwd, by = "site") %>% 
  filter(plot %in% unique(plot_locations$plot))

#### Get trees identified (segmented and classified) from the air
# Also join with CWD data
d <- 
  sf::st_read(here::here("data", "data_drone", "L3b", "model-classified-trees_all.gpkg"), 
              stringsAsFactors = FALSE)  

# Join classified trees with CWD data
dd <-
  d %>%
  dplyr::left_join(cwd, by = "site")

# Compare to Stovall et al., 2019 size categories
# dd <-
#   dd %>% 
#   dplyr::mutate(height_cat = case_when(height < 5 ~ "very small",
#                                        height >= 5 & height < 15 ~ "small",
#                                        height >= 15 & height < 30 ~ "medium",
#                                        height >= 30 ~ "large"))
# 
# 
# dd %>% 
#   sf::st_drop_geometry() %>% 
#   dplyr::group_by(height_cat) %>% 
#   summarize(n = n()) %>% 
#   dplyr::mutate(pct = n / sum(n))
# 
# ggplot(dd, aes(x = height)) + 
#   geom_histogram(bins = 200)


# subset classified trees to just those within plots that are identifiable
# from the air
air_trees <-
  dd %>% 
  sf::st_intersection(plot_locations) %>% 
  dplyr::mutate(species = ifelse(is.na(species), yes = "pipo", no = species),
                plot = as.character(plot))

atrees <-
  air_trees %>% 
  sf::st_drop_geometry() %>% 
  dplyr::select(site, plot, live, species, height) %>% 
  dplyr::mutate(ag = "a")

gtrees <-
  ground_trees %>% 
  sf::st_drop_geometry() %>% 
  dplyr::select(site, plot, live, species, height) %>% 
  dplyr::mutate(species = tolower(species)) %>% 
  dplyr::mutate(ag = "g")

trees <- 
  rbind(atrees, gtrees) %>% 
  dplyr::left_join(cwd_plot, by = "plot")

pipos <-
  trees %>% 
  dplyr::filter(species == "pipo") %>% 
  dplyr::group_by(plot, plot_cwd_zscore, live, ag) %>% 
  dplyr::summarize(mean_height = mean(height),
                   n = n())

pipos_wide <- 
  pipos %>% 
  tidyr::pivot_wider(names_from = "ag", values_from = c("mean_height", "n")) %>% 
  dplyr::mutate(hgt_diff = mean_height_a - mean_height_g,
                hgt_ratio = mean_height_a / mean_height_g) %>% 
  dplyr::mutate(site = substr(plot, start = 1, stop = 9)) %>% 
  dplyr::left_join(cwd, by = "site")

ggplot(pipos_wide[complete.cases(pipos_wide), ], aes(x = site_cwd_zscore, y = hgt_ratio, size = n_g, color = as.factor(live))) +
  geom_point() +
  geom_smooth(method = "lm")

summary(lm(hgt_ratio ~ site_cwd_zscore * live, data = pipos_wide))
ggplot(pipos_wide[complete.cases(pipos_wide), ], aes(x = plot_cwd_zscore, y = hgt_ratio, size = n_g, color = as.factor(live))) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(pipos_wide[complete.cases(pipos_wide), ], aes(x = site_cwd_zscore, y = hgt_diff, size = n_g, color = as.factor(live))) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "CWD z-score",
       y = "Height difference (drone minus ground)")

ggplot(pipos_wide[complete.cases(pipos_wide), ], aes(x = plot_cwd_zscore, y = hgt_ratio, size = n_g, color = as.factor(live))) +
  geom_point() +
  geom_smooth(method = "lm")


ggplot(pipos, aes(x = plot_cwd_zscore, y = mean_height, color = ag)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(vars(live))

ggplot(pipos_wide[complete.cases(pipos_wide), ], aes(x = plot_cwd_zscore, y = hgt_diff, size = n_g, color = as.factor(live))) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "CWD z-score",
       y = "Height difference (drone minus ground)")
  
ggplot(pipos_wide[complete.cases(pipos_wide), ], aes(x = plot_cwd_zscore, y = hgt_ratio, size = n_g, color = as.factor(live))) +
  geom_point() +
  geom_smooth(method = "lm")
  

amort <- 
  atrees %>% 
  group_by(plot) %>% 
  summarize(pct_mortality = 1 - mean(live)) %>% 
  dplyr::mutate(ag = "a")

gmort <- 
  gtrees %>% 
  group_by(plot) %>% 
  summarize(pct_mortality = 1 - mean(live)) %>% 
  dplyr::mutate(ag = "g")

aheight <- 
  atrees %>% 
  dplyr::filter(species == "pipo") %>% 
  group_by(plot, live) %>% 
  summarize(mn_hgt = mean(height)) %>% 
  dplyr::mutate(ag = "a")

ggplot(aheight, aes(x = mn_hgt, color = as.factor(live))) +
  geom_histogram()

gheight <- 
  gtrees %>% 
  dplyr::filter(species == "pipo") %>% 
  group_by(plot, live) %>% 
  summarize(mn_hgt = mean(height)) %>% 
  dplyr::mutate(ag = "g")

ggplot(gheight, aes(x = mn_hgt, color = as.factor(live))) +
  geom_histogram()


mort <- 
  rbind(amort, gmort) %>% 
  pivot_wider(names_from = ag, values_from = pct_mortality) %>% 
  dplyr::mutate(mort_diff = a - g)

ggplot(mort, aes(mort_diff)) + 
  geom_histogram()

range(mort$mort_diff)

