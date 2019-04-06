# Purpose: use an example site to show the rasterized product

library(tidyverse)
library(raster)
library(cowplot)
library(here)

r_eldo_3k_3 <- raster::brick(here::here("analyses/analyses_output/rasterized-trees/eldo_3k_3_rasterized-trees.tif"))
names(r_eldo_3k_3) <- c("live_count", "dead_count", "pipo_count", "non_pipo_count", "pipo_and_dead_count", "total_count", 
                        "live_tpha", "dead_tpha", "pipo_tpha", "non_pipo_tpha", "pipo_and_dead_tpha", "overall_tpha",
                        "live_ba", "dead_ba", "pipo_ba", "non_pipo_ba", "pipo_and_dead_ba", "total_ba",
                        "live_bapha", "dead_bapha", "pipo_bapha", "non_pipo_bapha", "pipo_and_dead_bapha", "overall_bapha",
                        "live_mean_ba", "dead_mean_ba", "pipo_mean_ba", "non_pipo_mean_ba", "pipo_and_dead_mean_ba", "overall_mean_ba",
                        "live_qmd", "dead_qmd", "pipo_qmd", "non_pipo_qmd", "pipo_and_dead_qmd", "overall_qmd",
                        "live_sdi_ha", "dead_sdi_ha", "pipo_sdi_ha", "non_pipo_sdi_ha", "pipo_and_dead_sdi_ha", "overall_sdi_ha",
                        "live_sdi_ac", "dead_sdi_ac", "pipo_sdi_ac", "non_pipo_sdi_ac", "pipo_and_dead_sdi_ac", "overall_sdi_ac",
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
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Fraction dead (%)")

ggsave(plot = live_and_dead_gg, filename = "figures/live-and-dead-count-rasterized.png")
ggsave(plot = prop_dead_gg, filename = "figures/proportion-dead-rasterized.png")
