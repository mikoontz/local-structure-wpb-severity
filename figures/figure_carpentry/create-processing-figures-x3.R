# Script to create a Digital Surface Model, Digital Terrain Model, and Canopy Height Model
# with a smaller file size (i.e., coarser resolution)

library(raster)
library(sf)
library(tidyverse)

all_sites <- list.files("data/data_output/site_data")
site <- all_sites[1]

dsm <- raster::raster(paste0("data/data_output/site_data/", site, "/3_dsm_ortho/1_dsm/", site, "_x3_dsm.tif"))
dtm <- raster::raster(paste0("data/data_output/site_data/", site, "/", site, "_2m-dtm.tif"))

site_bounds <- 
  sf::st_read(paste0("data/data_output/site_data/", site, "/", site, "_mission-footprint/", site, "_site-bounds.geoJSON")) %>% 
  sf::st_transform(proj4string(dsm))
  
dsm <- 
  dsm %>% 
  raster::crop(site_bounds)

dtm <-
  dtm %>% 
  raster::crop(site_bounds) %>% 
  raster::resample(dsm)

chm <- dsm - dtm
chm <-
  raster::focal(chm, w = matrix(1, 3, 3), mean)
chm[chm < 0] <- 0

ttops <- sf::st_read(paste0("data/data_output/site_data/", site, "/", site, "_ttops/", site, "_ttops.shp"))
crowns <- sf::st_read(paste0("data/data_output/site_data/", site, "/", site, "_crowns/", site, "_crowns.shp"))

cc <- sf::st_read(paste0("data/data_output/classified/model-classified/crown-shapefiles/", site, "_classified-crowns/", site, "_classified-crowns.shp"))

orig <- par()

### Plot the Digital Surface Model (DSM)
png(paste0("figures/", site, "_dsm.png"), width = 12, height = 10, units = "in", res = 600)
par(mar = c(0.1, 9.1, 0.1, 9.1))
plot(dsm, col = viridis(100),
     xlab = NA,
     ylab = NA,
     legend = FALSE,
     cex.axis = 1.5,
     cex.lab = 1.5,
     axes = FALSE,
     box = FALSE)
plot(dsm, col = viridis(100),
     legend.only = TRUE, 
     add = TRUE,
     legend.args = list(text='surface\nheight\n(m)', side = 3, line = 1.5, cex = 2.5),
     axis.args = list(cex.axis = 1.5),
     smallplot = c(0.9, 0.925, 0.3, 0.7))
mtext(side = 2, text = "latitude", cex = 2.5, line = 1.5, las = 1)
mtext(side = 1, text = "longitude", cex = 2.5, line = -2.5, las = 1)

dev.off()

### Plot the Digital Terrain Model (DTM)
png(paste0("figures/", site, "_dtm.png"), width = 12, height = 10, units = "in", res = 600)
par(mar = c(0.1, 9.1, 0.1, 9.1))
plot(dtm, col = viridis(100),
     xlab = NA,
     ylab = NA,
     legend = FALSE,
     cex.axis = 1.5,
     cex.lab = 1.5,
     axes = FALSE,
     box = FALSE)
plot(dtm, col = viridis(100),
     legend.only = TRUE, 
     add = TRUE,
     legend.args = list(text="ground\nelevation\n(m)", side = 3, line = 1.5, cex = 2.5),
     axis.args = list(cex.axis = 1.5),
     smallplot = c(0.9, 0.925, 0.3, 0.7))
mtext(side = 2, text = "latitude", cex = 2.5, line = 1.5, las = 1)
mtext(side = 1, text = "longitude", cex = 2.5, line = -2.5, las = 1)

dev.off()

### Plot the Canopy Height Model (CHM)
png(paste0("figures/", site, "_chm.png"), width = 12, height = 10, units = "in", res = 600)
par(mar = c(0.1, 9.1, 0.1, 9.1))
plot(chm, col = viridis(100),
     xlab = NA,
     ylab = NA,
     legend = FALSE,
     cex.axis = 1.5,
     cex.lab = 1.5,
     axes = FALSE,
     box = FALSE)
plot(chm, col = viridis(100),
     legend.only = TRUE, 
     add = TRUE,
     legend.args = list(text="canopy\nheight\n(m)", side = 3, line = 1.5, cex = 2.5),
     axis.args = list(cex.axis = 1.5),
     smallplot = c(0.9, 0.925, 0.3, 0.7))
mtext(side = 2, text = "latitude", cex = 2.5, line = 1.5, las = 1)
mtext(side = 1, text = "longitude", cex = 2.5, line = -2.5, las = 1)

dev.off()

### Plot the canopy height model with the tree tops identified
png(paste0("figures/", site, "_chm-with-ttops.png"), width = 12, height = 10, units = "in", res = 600)
par(mar = c(0.1, 9.1, 0.1, 9.1))
plot(chm, col = viridis(100),
     xlab = NA,
     ylab = NA,
     legend = FALSE,
     cex.axis = 1.5,
     cex.lab = 1.5,
     axes = FALSE,
     box = FALSE)
plot(chm, col = viridis(100),
     legend.only = TRUE, 
     add = TRUE,
     legend.args = list(text="canopy\nheight\n(m)", side = 3, line = 1.5, cex = 2.5),
     axis.args = list(cex.axis = 1.5),
     smallplot = c(0.9, 0.925, 0.3, 0.7))
mtext(side = 2, text = "latitude", cex = 2.5, line = 1.5, las = 1)
mtext(side = 1, text = "longitude", cex = 2.5, line = -2.5, las = 1)
plot(ttops$geometry, add = TRUE, pch = 19, cex = 0.5, col = "red")

dev.off()

### Plot just the tree tops
png(paste0("figures/", site, "_ttops.png"), width = 12, height = 10, units = "in", res = 600)
par(mar = c(0.1, 9.1, 0.1, 9.1))

plot(ttops$geometry, pch = 19, cex = 0.5,
     xlab = NA,
     ylab = NA,
     cex.axis = 1.5,
     cex.lab = 1.5,
     axes = FALSE)
mtext(side = 2, text = "latitude", cex = 2.5, line = 1.5, las = 1)
mtext(side = 1, text = "longitude", cex = 2.5, line = -2.5, las = 1)

dev.off()

### Plot the tree tops with the crowns
png(paste0("figures/", site, "_ttops-with-crowns.png"), width = 12, height = 10, units = "in", res = 600)
par(mar = c(0.1, 9.1, 0.1, 9.1))

plot(ttops$geometry, pch = ".",
     xlab = NA,
     ylab = NA,
     cex.axis = 1.5,
     cex.lab = 1.5,
     axes = FALSE)
mtext(side = 2, text = "latitude", cex = 2.5, line = 1.5, las = 1)
mtext(side = 1, text = "longitude", cex = 2.5, line = -2.5, las = 1)
plot(crowns$geometry, add = TRUE)

dev.off()

### Plot the classified tree tops
cc_ggplot <-
  ggplot(cc, aes(x = ttop_x, y = ttop_y, col = live_prob)) +
  geom_point(size = 0.9) +
  coord_equal() +
  xlab(label = "longitude") +
  ylab(label = "latitude") +
  scale_color_viridis_c(name = "Pr (live) ", option = "E") +
  theme_minimal() +
  theme(axis.text = element_blank(), 
        axis.title = element_text(size = 32), 
        legend.title = element_text(size = 32), 
        legend.text = element_text(size = 32), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = margin(t = 0, r = 5.5, b = 0, l = 5.5, "points"),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        legend.key.height = unit(0.075, "npc"))

ggsave(plot = cc_ggplot, filename = paste0("figures/", site, "_classified-ttops.png"), height = 10, width = 12, units = "in", dpi = 600)

# Example stem plots with orthos and mortality numbers
# eldo_4k_1 ~10.8% mortality (low)
# sequ_5k_2 ~32.4% mortality (about average)
# sier_4k_1 ~58.2% mortality (high)

# read summarized data

sd <- read_csv("data/data_output/summarized-non-spatial-site-data.csv")
