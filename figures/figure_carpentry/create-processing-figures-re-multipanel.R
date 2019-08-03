# Purpose: create a multi-panel processing figure that includes:

# Pix4D outputs (DSM; Point Cloud)
# DTM
# CHM
# Tree detection
# Segmentation
# Classification (live/dead; species)
# rasterized?

# 3x3

# DSM; Point Cloud; DTM
# CHM; ttops; crowns
# live/dead; species; rasterized

# alternatively, for rasterized
# rasterized live; rasterized dead
# rasterized host; rasterized non-host
# rasterized host qmd; rasterized non-host qmd

library(tidyverse)
library(raster)
library(viridis)
library(lidR)
library(sf)
library(png)
library(magick)
library(cowplot)
library(here)
library(tmap)
library(tictoc)

source("data/data_carpentry/make-processing-checklist.R")

# These sites were processed with their X3 and RedEdge imagery combined so some of their
# output products will be in a slightly different place in the project directory
merged_sites <- c("eldo_3k_2",
                  "eldo_3k_3",
                  "eldo_4k_2")

example_site <- 
  sites_checklist %>% 
  dplyr::filter(rasterized_trees_check) %>% 
  slice(3) %>% 
  pull(site)

# The Digital Surface Model (dsm) is the ~5cm resolution raster representing
# the surface (ground + objects on top) that the drone flew over
if (example_site %in% merged_sites) {
  dsm <- raster::raster(x = here::here(paste0("data/data_output/site_data/", example_site, "/", "3_dsm_ortho/1_dsm/", example_site, "_dsm.tif")))
} else {
  dsm <- raster::raster(x = here::here(paste0("data/data_output/site_data/", example_site, "/", example_site, "_re/3_dsm_ortho/1_dsm/", example_site, "_re_dsm.tif")))
}

# if (example_site %in% merged_sites) {
#   point_cloud <- lidR::readLAS(files = here::here(paste0("data/data_output/site_data/", example_site, "/", "2_densification/point_cloud/", example_site, "_densified_point_cloud.las")))
# } else {
#   point_cloud <- lidR::readLAS(files = here::here(paste0("data/data_output/site_data/", example_site, "/", example_site, "_re/2_densification/point_cloud/", example_site, "_re_Green_densified_point_cloud.las")))
# }
# plot(point_cloud, colorPalette = viridis(100))

dtm <- raster::raster(paste0("data/data_output/site_data/", example_site, "/", example_site, "_dtm.tif"))
current_chm_rough <- raster::raster(paste0("data/data_output/site_data/", example_site, "/", example_site, "_chm.tif"))
chm <- raster::focal(current_chm_rough, w = matrix(1/9, nrow = 3, ncol = 3))
# current_chm <- current_chm_rough
chm[raster::getValues(chm) < 0] <- 0


dsm_gg <-
  dsm %>% 
  as.data.frame(xy = TRUE) %>% 
  setNames(c("x", "y", "z")) %>% 
  ggplot(mapping = aes(x = x, y = y, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c(na.value = "white") +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude", fill = "Meters\nabove\nsea\nlevel") +
  theme_bw()

dtm_gg <-
  dtm %>% 
  as.data.frame(xy = TRUE) %>% 
  setNames(c("x", "y", "z")) %>% 
  ggplot(mapping = aes(x = x, y = y, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c(na.value = "white") +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude", fill = "Meters\nabove\nsea\nlevel") +
  theme_bw()

chm_gg <-
  chm %>% 
  as.data.frame(xy = TRUE) %>% 
  setNames(c("x", "y", "z")) %>% 
  ggplot(mapping = aes(x = x, y = y, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c(na.value = "white") +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude", fill = "Meters\nabove\nground\nlevel") +
  theme_bw()


ttops <- 
  sf::st_read(paste0("data/data_output/site_data/", example_site, "/", example_site, "_ttops/", example_site, "_ttops.shp")) %>% 
  dplyr::mutate(x = st_coordinates(.)[,1],
                y = st_coordinates(.)[,2])

crowns <- sf::st_read(paste0("data/data_output/site_data/", example_site, "/", example_site, "_crowns/", example_site, "_crowns.shp"))

crop_poly <- 
  matrix(c(715200, 4270500,
           715400, 4270500,
           715400, 4270700,
           715200, 4270700,
           715200, 4270500),
         ncol = 2, byrow = TRUE) %>% 
  list() %>% 
  st_polygon()

plot(crop_poly, col = "red", add = TRUE)

ttops_gg <-
  ggplot(ttops, aes(x = x, y = y)) +
  geom_point(cex = 0.2) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

ttops_crop_gg <-
  ggplot(ttops[crop_poly, ], aes(x = x, y = y)) +
  geom_point(cex = 0.2) +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

crowns_gg <-
  ggplot(crowns) +
  geom_sf() +
  coord_sf() +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

crowns_crop_gg <-
  ggplot(crowns[crop_poly, ]) +
  geom_sf() +
  coord_sf() +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

png("figures/eldo_3k_3_ttops.png", width = 6, height = 6, units = "in", res = 2400)
plot(st_geometry(ttops), axes = TRUE, xlab = "Longitude (m)", ylab = "Latitude (m)", pch = 19, cex = 0.2)
dev.off()

png("figures/eldo_3k_3_ttops_with_crop_poly.png", width = 6, height = 6, units = "in", res = 2400)
plot(st_geometry(ttops), axes = TRUE, xlab = "Longitude (m)", ylab = "Latitude (m)", pch = 19, cex = 0.2)
plot(crop_poly, col = adjustcolor("red", alpha.f = 0.2), border = "red", lwd = 3, add = TRUE)
dev.off()

png("figures/eldo_3k_3_ttops_cropped.png", width = 6, height = 6, units = "in", res = 2400)
plot(st_geometry(ttops[crop_poly, ]), axes = TRUE, xlab = "Longitude (m)", ylab = "Latitude (m)", pch = 19, cex = 0.2)
dev.off()

png("figures/eldo_3k_3_crowns.png", width = 6, height = 6, units = "in", res = 2400)
plot(st_geometry(crowns), axes = TRUE, xlab = "Longitude (m)", ylab = "Latitude (m)")
dev.off()

png("figures/eldo_3k_3_crowns_with_crop_poly.png", width = 6, height = 6, units = "in", res = 2400)
plot(st_geometry(crowns), axes = TRUE, xlab = "Longitude (m)", ylab = "Latitude (m)")
plot(crop_poly, col = adjustcolor("red", alpha.f = 0.2), border = "red", lwd = 3, add = TRUE)
dev.off()

png("figures/eldo_3k_3_crowns_cropped.png", width = 6, height = 6, units = "in", res = 2400)
plot(st_geometry(crowns[crop_poly, ]), axes = TRUE, xlab = "Longitude (m)", ylab = "Latitude (m)")
dev.off()

index <- raster::brick(paste0("data/data_output/site_data/", example_site, "/", example_site, "_index.tif"))
index_rgb <- r[[7:9]]
crop_rgb <- raster::crop(index_rgb, y = as(crop_poly, "Spatial"))

png("figures/eldo_3k_3_ortho_rgb.png", width = 6, height = 5.5, units = "in", res = 2400)
plotRGB(index_rgb)
dev.off()

png("figures/eldo_3k_3_ortho_rgb_with_crop_poly.png", width = 6, height = 5.5, units = "in", res = 2400)
plotRGB(index_rgb)
plot(crop_poly, col = adjustcolor("black", alpha.f = 0.2), lwd = 3, add = TRUE)
dev.off()

png("figures/eldo_3k_3_ortho_rgb_cropped.png", width = 6, height = 5.5, units = "in", res = 2400)
plotRGB(crop_rgb)
dev.off()

if (!file.exists(here::here(paste0("data/data_output/site_data/", example_site, "/", example_site, "_classified-crowns.gpkg")))) {
  
  cc <- sf::st_read(here::here("analyses/analyses_output/classified-trees.geojson"))
  
  site_cc <- 
    cc %>% 
    dplyr::filter(site == example_site) %>% 
    sf::st_transform(3310) %>% 
    dplyr::mutate(live = ifelse(live == 1, yes = "live", no = "dead")) %>% 
    dplyr::mutate(host = ifelse((live == "live" & species == "pipo") | live == "dead",
                                yes = "host",
                                no = "non-host"))
  
  sf::st_write(obj = site_cc, dsn = here::here(paste0("data/data_output/site_data/", example_site, "/", example_site, "_classified-crowns.gpkg")))
}

site_cc <- sf::st_read(dsn = here::here(paste0("data/data_output/site_data/", example_site, "/", example_site, "_classified-crowns.gpkg")))

live_dead_gg <-
  ggplot(site_cc, aes(x = local_x, y = local_y, color = factor(live, levels = c("live", "dead")))) +
  geom_point(cex = 0.55) +
  labs(x = "Longitude (m)", y = "Latitude (m)", color = "Status") +
  theme_bw() +
  coord_equal() +
  scale_x_continuous(limits = c(714800, 715400)) +
  scale_color_viridis_d()

host_gg <-
  ggplot(site_cc, aes(x = local_x, y = local_y, color = factor(host, levels = c("non-host", "host")))) +
  geom_point(cex = 0.55) +
  labs(x = "Longitude (m)", y = "Latitude (m)", color = "Species") +
  theme_bw() +
  coord_equal() +
  scale_x_continuous(limits = c(714800, 715400)) +
  scale_color_viridis_d()


ggsave(plot = dsm_gg, filename = "figures/eldo_3k_3_dsm.png", width = 6, units = "in")
ggsave(plot = dtm_gg, filename = "figures/eldo_3k_3_dtm.png", width = 6, units = "in")
ggsave(plot = chm_gg, filename = "figures/eldo_3k_3_chm.png", width = 6, units = "in")
ggsave(plot = ttops_gg, filename = "figures/eldo_3k_3_ttops.png", width = 6, units = "in")
ggsave(plot = crowns_gg, filename = "figures/eldo_3k_3_crowns.png", width = 6, units = "in")
ggsave(plot = live_dead_gg, filename = "figures/eldo_3k_3_live_dead.png", width = 6, height = 5.5, units = "in")
ggsave(plot = host_gg, filename = "figures/eldo_3k_3_host_nonhost.png", width = 6, height = 5.5, units = "in")

raw_photo <- ggdraw() + draw_image("figures/eldo_3k_3_raw_x3_100MEDIA_DJI_0705.JPG")
# point_cloud_png <- ggdraw() + draw_image("figures/eldo_3k_3_point_cloud.png")
# dsm_png <- ggdraw() + draw_image("figures/eldo_3k_3_dsm.png")
# dtm_png <- ggdraw() + draw_image("figures/eldo_3k_3_dtm.png")
# chm_png <- ggdraw() + draw_image("figures/eldo_3k_3_chm.png")
# ttops_png <- ggdraw() + draw_image("figures/eldo_3k_3_ttops.png")
# crowns_png <- ggdraw() + draw_image("figures/eldo_3k_3_crowns.png")
# live_dead_png <- ggdraw() + draw_image("figures/eldo_3k_3_live_dead.png")
# host_nonhost_png <- ggdraw() + draw_image("figures/eldo_3k_3_host_nonhost.png")
# 
# panel_plot <- plot_grid(dsm_png, point_cloud_png,
#                         dtm_png, chm_png,
#                         ttops_png, crowns_png,
#                         live_dead_gg, host_nonhost_png,
#                         nrow = 4, ncol = 2)
# 
# panel_plot
