# Make a map of the study area

library(tidyverse)
library(sf)
library(tmap)
library(ggrepel)
library(spData)
library(grid)
library(tmaptools)

sn <- 
  st_read("data/data_output/sierra-nevada-jepson/sierra-nevada-jepson.shp") %>% 
  st_transform(3310)
usfs <- 
  st_read("data/features/S_USA.AdministrativeForest/S_USA.AdministrativeForest.shp") %>% 
  st_transform(3310)
veg_plots <- 
  st_read("data/features/plot-centers_ground-gps-measured.kml") %>% 
  st_transform(3310) %>% 
  group_by(forest, elevation_band, site) %>% 
  summarize() %>% 
  st_centroid() %>% 
  dplyr::mutate(elev = case_when(elevation_band == 3000 ~ "1219-1524",
                                 elevation_band == 4000 ~ "1219-1524",
                                 elevation_band == 5000 ~ "1524-1828",
                                 elevation_band == 6000 ~ "1828-2133"))

western_states <- 
  us_states %>% 
  dplyr::filter(NAME %in% c("California", "Nevada", "Oregon", "Idaho")) %>% 
  st_transform(3310)

relevant_forests <-
  usfs %>% 
  dplyr::filter(FORESTNAME %in% c("Eldorado National Forest", 
                                  "Stanislaus National Forest",
                                  "Sierra National Forest",
                                  "Sequoia National Forest")) %>% 
  dplyr::mutate(FORESTNAME = case_when(FORESTNAME == "Eldorado National Forest" ~ "Eldorado NF",
                                       FORESTNAME == "Stanislaus National Forest" ~ "Stanislaus NF",
                                       FORESTNAME == "Sierra National Forest" ~ "Sierra NF",
                                       FORESTNAME == "Sequoia National Forest" ~ "Sequoia NF")) %>% 
  dplyr::mutate(x = st_coordinates(st_centroid(st_geometry(.)))[, 1],
                y = st_coordinates(st_centroid(st_geometry(.)))[, 2])

study_extent_bbox <- st_bbox(sn) + c(0, 0, 90000, 0)

study_extent <-
  tm_shape(sn, bbox = study_extent_bbox) +
  tm_borders() +
  tm_shape(relevant_forests) +
  tm_fill(col = "FORESTNAME", title = "Forest name", palette = get_brewer_pal("Accent", n = 4)) +
  tm_shape(veg_plots) +
  tm_dots(shape = "elev", title.shape = "Elevation band (m)", size = 0.15) +
  tm_layout(legend.position = c("left", "bottom"), legend.width = 1, outer.margins = c(0.02, 0.02, 0, 0)) +
  tm_ylab(text = "Latitude (km)", space = 2) +
  tm_xlab(text = "Longitude (km)", space = 1.5) +
  tm_grid(alpha = 0.2, labels.inside.frame = FALSE, labels.format = list(fun = function(x) as.character(x / 1000))) +
  tm_compass(type = "arrow", position = c(0.85, 0.05))

broader_extent_bbox <- st_bbox(western_states) + c(0, 0, -450000, -750000)

broader_extent <-
  tm_shape(western_states, projection = 3310, bbox = broader_extent_bbox) + 
  tm_borders() +
  tm_shape(sn) +
  tm_fill() +
  tm_shape(relevant_forests) +
  tm_borders()

  # tm_compass(type = "arrow", position = c("left", "bottom"))
  # tm_scale_bar(breaks = c(0, 100, 200), size = 0.5, position = c("left", "bottom"))

broader_extent

# panel_plot <- tmap_arrange(broader_extent, study_extent)
# panel_plot
# tmap_save(tm = panel_plot, filename = "figures/study-geographic-extent.png")

vp <- viewport(x = grid::unit(1.3, "npc"), y = grid::unit(1, "npc"), height = 0.4, just = c(1, 1))

# study_extent
# print(broader_extent, vp = vp)

tmap_save(tm = study_extent, filename = "figures/study-geographic-extent-inset.png", insets_tm = broader_extent, insets_vp = vp, height = 7, units = "in", dpi = 600)
