# Purpose: use an example site to show the rasterized product

library(tidyverse)
library(raster)
library(cowplot)
library(sf)
library(here)

ttops_eldo_3k_3 <- st_read("data/data_output/site_data/eldo_3k_3/eldo_3k_3_ttops/eldo_3k_3_ttops.shp")

r_eldo_3k_3 <- raster::brick(here::here("analyses/analyses_output/rasterized-trees/eldo_3k_3_rasterized-trees.tif"))
names(r_eldo_3k_3) <- c("live_count", "dead_count", "pipo_count", "non_pipo_count", "pipo_and_dead_count", "total_count", 
                        "live_tpha", "dead_tpha", "pipo_tpha", "non_pipo_tpha", "pipo_and_dead_tpha", "overall_tpha",
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
prop_dead <- dead / (live + dead) * 100

live_df <- 
  live %>% 
  as.data.frame(xy = TRUE) %>% 
  dplyr::rename(count = live_count) %>% 
  dplyr::mutate(live = "live")

dead_df <-
  dead %>%
  as.data.frame(xy = TRUE) %>% 
  dplyr::rename(count = dead_count) %>% 
  dplyr::mutate(live = "dead")

live_and_dead_df <-
  rbind(live_df, dead_df)

prop_dead_df <-
  prop_dead %>% 
  as.data.frame(xy = TRUE)

live_gg <-
  ggplot(live_df, aes(x, y, fill = count)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c()

dead_gg <-
  ggplot(dead_df, aes(x, y, fill = count)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c()

live_and_dead_gg <-
  ggplot(live_and_dead_df, aes(x, y, fill = count)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c() +
  facet_grid(~ live) +
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Tree count")

prop_dead_gg <-
  ggplot(prop_dead_df, aes(x, y, fill = layer)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c() +
  labs(x = "Longitude (m)",
       y = "Latitude (m)",
       fill = "Fraction\ndead (%)")

ggsave(plot = live_and_dead_gg, filename = "figures/live-and-dead-count-rasterized.png")
ggsave(plot = prop_dead_gg, filename = "figures/eldo_3k_3_prop-dead-rasterized.png", width = 6, height = 4.5, units = "in")



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
       fill = "Fraction\nhost (%)")

ggsave(plot = prop_host_gg, filename = "figures/eldo_3k_3_prop-host-rasterized.png", width = 6, height = 4.5, units = "in")
