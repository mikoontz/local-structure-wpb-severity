# Purpose: create a figure showing the cropped orthomosaic

library(tidyverse)
library(raster)
library(tmap)
library(RStoolbox)

ortho <- raster::brick("figures/eldo_3k_3_remote-plot-data/eldo_3k_3_2_ortho.tif")

names(ortho) <- c("re_b", "re_g", "re_r", "re_re", "re_nir", "x3_r", "x3_g", "x3_b", "x3_trans")

ortho

ortho_rgb <- ortho[[c("x3_r", "x3_g", "x3_b")]]

ortho_rgb_gg <-
  ggRGB(ortho_rgb, r = 1, g = 2, b = 3) +
  labs(x = "Longitude",
       y = "Latitude")

ggsave(plot = ortho_rgb_gg, filename = "figures/eldo_3k_3_2_ortho-rgb.png", width = 6, units = "in")
# tm_shape(ortho_rgb) +
  # tm_rgb()

