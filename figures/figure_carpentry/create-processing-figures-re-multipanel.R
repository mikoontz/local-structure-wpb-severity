# Purpose: create a multi-panel processing figure that includes:

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

# figure specs in inches
# height <- 4
# width <- 6

# These sites were processed with their X3 and RedEdge imagery combined so some of their
# output products will be in a slightly different place in the project directory
merged_sites <- c("eldo_3k_2",
                  "eldo_3k_3",
                  "eldo_4k_2")

example_site <- "eldo_3k_3"

####
# Level 1: basic photogrammetric outputs
####

# digital surface model, uncorrected orthomosaic, point cloud (as a screen shot; not implemented here)

if (example_site %in% merged_sites) {
  dsm <- raster::raster(x = here::here(paste0("data/data_output/site_data/", example_site, "/", "3_dsm_ortho/1_dsm/", example_site, "_dsm.tif")))
} else {
  dsm <- raster::raster(x = here::here(paste0("data/data_output/site_data/", example_site, "/", example_site, "_re/3_dsm_ortho/1_dsm/", example_site, "_re_dsm.tif")))
}

dsm_gg <-
  dsm %>% 
  as.data.frame(xy = TRUE) %>% 
  setNames(c("x", "y", "z")) %>%
  ggplot(mapping = aes(x = x, y = y, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c(na.value = "white") +
  coord_equal() +
  labs(x = "Longitude (m)", y = "Latitude (m)", fill = "Meters\nabove\nsea\nlevel") +
  theme_bw()

ggsave(plot = dsm_gg, filename = "figures/L1_eldo_3k_3_dsm.png", width = 6, height = 4.5, units = "in")

index <- raster::brick(paste0("data/data_output/site_data/", example_site, "/", example_site, "_index.tif"))
index_rgb <- index[[7:9]]

png("figures/L1_eldo_3k_3_ortho_rgb.png", width = 6, height = 6, units = "in", res = 2400)
plotRGB(index_rgb)
dev.off()

####
# Level 2: radiometric and/or geometric corrected photogrammetric outputs
####

# radiometric corrected: index (red and near infrared)
# geometric corrected: digital terrain model, canopy height model

dtm <- raster::raster(paste0("data/data_output/site_data/", example_site, "/", example_site, "_dtm.tif"))
current_chm_rough <- raster::raster(paste0("data/data_output/site_data/", example_site, "/", example_site, "_chm.tif"))
chm <- raster::focal(current_chm_rough, w = matrix(1/9, nrow = 3, ncol = 3))
chm[raster::getValues(chm) < 0] <- 0

dtm_gg <-
  dtm %>% 
  as.data.frame(xy = TRUE) %>% 
  setNames(c("x", "y", "z")) %>% 
  ggplot(mapping = aes(x = x, y = y, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c(na.value = "white") +
  coord_equal() +
  labs(x = "Longitude (m)", y = "Latitude (m)", fill = "Meters\nabove\nsea\nlevel") +
  theme_bw()

chm_gg <-
  chm %>% 
  as.data.frame(xy = TRUE) %>% 
  setNames(c("x", "y", "z")) %>% 
  ggplot(mapping = aes(x = x, y = y, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c(na.value = "white") +
  coord_equal() +
  labs(x = "Longitude (m)", y = "Latitude (m)", fill = "Meters\nabove\nground\nlevel") +
  theme_bw()

ndvi <- index[[c(3, 5)]]

r_gg <-
  ndvi[[1]] %>% 
  as.data.frame(xy = TRUE) %>% 
  setNames(c("x", "y", "z")) %>% 
  ggplot(mapping = aes(x = x, y = y, fill = z)) +
  geom_raster() +
  scale_fill_gradient(low = "#4D4D4D", high = "#E6E6E6", na.value = "white") +
  coord_equal() +
  labs(x = "Longitude (m)", y = "Latitude (m)", fill = "Surface\nreflectance\n(red)") +
  theme_bw()

nir_gg <-
  ndvi[[2]] %>% 
  as.data.frame(xy = TRUE) %>% 
  setNames(c("x", "y", "z")) %>% 
  ggplot(mapping = aes(x = x, y = y, fill = z)) +
  geom_raster() +
  scale_fill_gradient(low = "#4D4D4D", high = "#E6E6E6", na.value = "white") +
  coord_equal() +
  labs(x = "Longitude (m)", y = "Latitude (m)", fill = "Surface\nreflectance\n(near\ninfrared)") +
  theme_bw()


ggsave(plot = dtm_gg, filename = "figures/L2_eldo_3k_3_dtm.png", width = 6, height = 4.5, units = "in")
ggsave(plot = chm_gg, filename = "figures/L2_eldo_3k_3_chm.png", width = 6, height = 4.5, units = "in")
ggsave(filename = "figures/L2_eldo_3k_3_ortho_r.png", plot = r_gg, width = 6, height = 4.5, units = "in")
ggsave(filename = "figures/L2_eldo_3k_3_ortho_nir.png", plot = nir_gg, width = 6, height = 4.5, units = "in")

####
# Level 3: domain specific information extraction
####

# Polygon for zooming into an area in the eldo_3k_1 site to show better detail on the tree tops and crowns
crop_poly <- 
  matrix(c(715200, 4270500,
           715400, 4270500,
           715400, 4270700,
           715200, 4270700,
           715200, 4270500),
         ncol = 2, byrow = TRUE) %>% 
  list() %>% 
  st_polygon() %>% 
  st_sfc(crs = 32610) %>% 
  st_transform(3310)

# Level 3a: spectral OR geometric information
# tree tops, crown segments

ndvi_r <- overlay(ndvi[[1]], 
                  ndvi[[2]],
                  fun = function(red, nir) { return( (nir - red) / (nir + red) ) })

ndvi_r_projected <- projectRaster(from = ndvi_r, crs = st_crs(3310)$proj4string)

ndvi_gg <-
  ndvi_r_projected %>% 
  as.data.frame(xy = TRUE) %>% 
  setNames(c("x", "y", "z")) %>% 
  ggplot(mapping = aes(x = x, y = y, fill = z)) +
  geom_raster() +
  scale_fill_gradient(low = "#E6E6E6", high = "#228B22", na.value = "white") +
  # scale_fill_viridis_c(option = "E", na.value = "white") +
  coord_equal() +
  labs(x = "Longitude (m)", y = "Latitude (m)", fill = "NDVI") +
  theme_bw()

ggsave(filename = "figures/L3a_eldo_3k_3_ortho_ndvi_projected.png", plot = ndvi_gg, width = 6, height = 4.5, units = "in")

ttops <- 
  sf::st_read(paste0("data/data_output/site_data/", example_site, "/", example_site, "_ttops/", example_site, "_ttops.shp")) %>% 
  sf::st_transform(3310) %>% 
  dplyr::mutate(x = st_coordinates(.)[,1],
                y = st_coordinates(.)[,2])

crowns <- sf::st_read(paste0("data/data_output/site_data/", example_site, "/", example_site, "_crowns/", example_site, "_crowns.shp")) %>% 
  sf::st_transform(3310)

png("figures/L3a_eldo_3k_3_ttops_with_crop_poly.png", width = 6, height = 6, units = "in", res = 2400)
par(mar = c(5.1, 4.1, 1.6, 1.6))
plot(st_geometry(ttops), axes = TRUE, xlab = "Longitude (m)", ylab = "Latitude (m)", pch = 19, cex = 0.2)
plot(crop_poly, col = adjustcolor("red", alpha.f = 0.2), border = "red", lwd = 3, add = TRUE)
dev.off()

png("figures/L3a_eldo_3k_3_ttops_cropped.png", width = 6, height = 6, units = "in", res = 2400)
par(mar = c(5.1, 4.1, 1.6, 1.6))
plot(st_geometry(ttops[crop_poly, ]), axes = TRUE, xlab = "Longitude (m)", ylab = "Latitude (m)", pch = 19, cex = 0.2)
dev.off()

png("figures/L3a_eldo_3k_3_crowns_cropped.png", width = 6, height = 6, units = "in", res = 2400)
par(mar = c(5.1, 4.1, 1.6, 1.6))
plot(st_geometry(crowns[crop_poly, ]), axes = TRUE, xlab = "Longitude (m)", ylab = "Latitude (m)")
dev.off()

# Level 3b: spectral AND geometric information
# live/dead classified tree tops, host/non-host classified tree tops

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

site_cc <- 
  sf::st_read(dsn = here::here(paste0("data/data_output/site_data/", example_site, "/", example_site, "_classified-crowns.gpkg"))) %>% 
  dplyr::mutate(x = st_coordinates(.)[,1],
                y = st_coordinates(.)[,2])


png(filename = "figures/L3b_eldo_3k_3_live_dead.png", width = 6, height = 6, units = "in", res = 2400)
par(mar = c(5.1, 4.1, 2.1, 2.1))
plot(x = site_cc$x, y = site_cc$y, 
     col = viridis(2)[factor(site_cc$live, levels = c("live", "dead"))], 
     pch = 19, cex = 0.35, asp = 1, bty = "L",
     xaxt = "n", yaxt = "n",
     xlab = "Longitude (m)", ylab = "Latitude (m)")
legend(x = -46275, y = 60450, xpd = NA, legend = c("live", "dead"), xjust = 0.5, pch = 19, col = viridis(2), horiz = TRUE, bty = "n")
axis(side = 1, at = seq(-46500, -46050, length.out = 3))
axis(side = 2, at = seq(59850, 60350, length.out = 3))
dev.off()


png(filename = "figures/L3b_eldo_3k_3_host_nonhost.png", width = 6, height = 6, units = "in", res = 2400)
par(mar = c(5.1, 4.1, 2.1, 2.1))
plot(x = site_cc$x, y = site_cc$y, 
     col = viridis(2)[factor(site_cc$host, levels = c("non-host", "host"))], 
     pch = 19, cex = 0.35, asp = 1, bty = "L",
     xaxt = "n", yaxt = "n",
     xlab = "Longitude (m)", ylab = "Latitude (m)")
legend(x = -46275, y = 60450, xpd = NA, legend = c("non-host", "host"), xjust = 0.5, pch = 19, col = viridis(2), horiz = TRUE, bty = "n")
axis(side = 1, at = seq(-46500, -46050, length.out = 3))
axis(side = 2, at = seq(59850, 60350, length.out = 3))
dev.off()

####
# Level 4: aggregated to regular grids
####

# rasterized version of live/dead classification at 20m x 20m grid cells
# rasterized version of host/non-host classification at 20m x 20m grid cells

r_eldo_3k_3 <- raster::brick(here::here("analyses/analyses_output/rasterized-trees/eldo_3k_3_rasterized-trees.tif"))
names(r_eldo_3k_3) <- c("live_count", "dead_count", "pipo_count", "non_pipo_count", "pipo_and_dead_count", "total_count", 
                        "live_tpha", "dead_tpha", "pipo_tpha", "non_pipo_tpha", "pipo_and_dead_tpha", "overall_tpha",
                        "live_mean_height", "dead_mean_height", "pipo_mean_height", "non_pipo_mean_height", "pipo_and_dead_mean_height", "overall_mean_height",
                        "live_ba", "dead_ba", "pipo_ba", "non_pipo_ba", "pipo_and_dead_ba", "total_ba",
                        "live_bapha", "dead_bapha", "pipo_bapha", "non_pipo_bapha", "pipo_and_dead_bapha", "overall_bapha",
                        "live_mean_ba", "dead_mean_ba", "pipo_mean_ba", "non_pipo_mean_ba", "pipo_and_dead_mean_ba", "overall_mean_ba",
                        "live_mean_dbh", "dead_mean_dbh", "pipo_mean_dbh", "non_pipo_mean_dbh", "pipo_and_dead_mean_dbh", "overall_mean_dbh",
                        "live_qmd", "dead_qmd", "pipo_qmd", "non_pipo_qmd", "pipo_and_dead_qmd", "overall_qmd",
                        "live_sdi_ha", "dead_sdi_ha", "pipo_sdi_ha", "non_pipo_sdi_ha", "pipo_and_dead_sdi_ha", "overall_sdi_ha",
                        "live_sdi_ac", "dead_sdi_ac", "pipo_sdi_ac", "non_pipo_sdi_ac", "pipo_and_dead_sdi_ac", "overall_sdi_ac",
                        "live_mean_voronoi", "dead_mean_voronoi", "pipo_mean_voronoi", "non_pipo_mean_voronoi", "pipo_and_dead_mean_voronoi", "overall_mean_voronoi",
                        "live_mean_nn1", "dead_mean_nn1", "pipo_mean_nn1", "non_pipo_mean_nn1", "pipo_and_dead_mean_nn1", "overall_mean_nn1",
                        "live_mean_nn2", "dead_mean_nn2", "pipo_mean_nn2", "non_pipo_mean_nn2", "pipo_and_dead_mean_nn2", "overall_mean_nn2",
                        "live_mean_nn3", "dead_mean_nn3", "pipo_mean_nn3", "non_pipo_mean_nn3", "pipo_and_dead_mean_nn3", "overall_mean_nn3",
                        "live_sd_nn1", "dead_sd_nn1", "pipo_sd_nn1", "non_pipo_sd_nn1", "pipo_and_dead_sd_nn1", "overall_sd_nn1",
                        "live_cov_nn1", "dead_cov_nn1", "pipo_cov_nn1", "non_pipo_cov_nn1", "pipo_and_dead_cov_nn1", "overall_cov_nn1",
                        "live_sd_voronoi", "dead_sd_voronoi", "pipo_sd_voronoi", "non_pipo_sd_voronoi", "pipo_and_dead_sd_voronoi", "overall_sd_voronoi",
                        "live_cov_voronoi", "dead_cov_voronoi", "pipo_cov_voronoi", "non_pipo_cov_voronoi", "pipo_and_dead_cov_voronoi", "overall_cov_voronoi",
                        "local_cwd", "local_cwd_zscore")

live <- r_eldo_3k_3[["live_count"]]
dead <- r_eldo_3k_3[["dead_count"]]
pipo_height <- r_eldo_3k_3[["pipo_and_dead_mean_height"]]
total_count <- r_eldo_3k_3[["total_count"]]

prop_dead <- dead / (live + dead) * 100

prop_dead_df <-
  prop_dead %>% 
  as.data.frame(xy = TRUE)

prop_dead_gg <-
  ggplot(prop_dead_df, aes(x, y, fill = layer)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c() +
  labs(x = "Longitude (m)",
       y = "Latitude (m)",
       fill = "Fraction\ndead (%)") +
  theme_classic() +
  theme(text = element_text(size = 15), 
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))

nonhost <- r_eldo_3k_3[["non_pipo_count"]]
host <- r_eldo_3k_3[["pipo_and_dead_count"]]
prop_host <- host / (nonhost + host) * 100

prop_host_df <-
  prop_host %>% 
  as.data.frame(xy = TRUE)

prop_host_gg <-
  ggplot(prop_host_df, aes(x, y, fill = layer)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c() +
  labs(x = "Longitude (m)",
       y = "Latitude (m)",
       fill = "Fraction\nhost (%)") +
  theme_classic() +
  theme(text = element_text(size = 15), 
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))

pipo_height_df <-
  pipo_height %>% 
  as.data.frame(xy = TRUE)

pipo_height_gg <-
  ggplot(pipo_height_df, aes(x, y, fill = pipo_and_dead_mean_height)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c() +
  labs(x = "Longitude (m)",
       y = "Latitude (m)",
       fill = "Mean\nhost\nheight\n(m)") +
  theme_classic() +
  theme(text = element_text(size = 15), 
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))

total_count_df <-
  total_count %>% 
  as.data.frame(xy = TRUE)

total_count_gg <-
  ggplot(total_count_df, aes(x, y, fill = total_count)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c() +
  labs(x = "Longitude (m)",
       y = "Latitude (m)",
       fill = "Tree\ncount") +
  theme_classic() +
  theme(text = element_text(size = 15), 
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))

ggsave(plot = prop_dead_gg, filename = "figures/L4_eldo_3k_3_prop-dead-rasterized.png", width = 6, height = 5, units = "in")

ggsave(plot = prop_host_gg, filename = "figures/L4_eldo_3k_3_prop-host-rasterized.png", width = 6, height = 5, units = "in")

ggsave(plot = pipo_height_gg, filename = "figures/L4_eldo_3k_3_pipo-height-rasterized.png", width = 6, height = 5, units = "in")

ggsave(plot = total_count_gg, filename = "figures/L4_eldo_3k_3_total-count-rasterized.png", width = 6, height = 5, units = "in")
