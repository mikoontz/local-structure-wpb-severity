library(sf)
library(raster)
library(tabularaster)
library(tidyverse)
library(viridis)


# character vector of all study sites
all_sites <- list.files("data/data_output")

sites_checklist <- 
  expand.grid(
    c("eldo", "stan", "sier", "sequ"), 
    c("3k", "4k", "5k", "6k"), 
    1:3) %>% 
  as.data.frame() %>% 
  setNames(c("forest", "elev", "rep")) %>% 
  dplyr::filter(!(forest == "sequ" & elev == "3k")) %>%
  dplyr::filter(!(forest != "sequ" & elev == "6k")) %>% 
  dplyr::arrange(forest, elev, rep) %>% 
  dplyr::mutate(site = paste(forest, elev, rep, sep = "_")) %>% 
  dplyr::select(-forest, -elev, -rep) %>% 
  dplyr::mutate(ttops_check = file.exists(paste0("data/data_output/", site, "/", site, "_ttops/", site, "_ttops.shp"))) %>% 
  dplyr::mutate(crowns_check = file.exists(paste0("data/data_output/", site, "/", site, "_crowns/", site, "_crowns.shp"))) %>% 
  dplyr::mutate(classified_check = file.exists(paste0("data/data_output/", site, "/", site, "_classified/", site, "_classified.shp")))

sites_to_process <-
  sites_checklist %>% 
  dplyr::filter(!classified_check & site %in% all_sites) %>% 
  dplyr::select(site) %>% 
  dplyr::pull()

hand_classified_sites <-
  list.files("data/data_output")
  sf::st_read("data/data_output/eldo_3k_1/eldo_3k_1_live-or-dead/eldo_3k_1_live-or-dead.shp")

i = 1
for (i in seq_along(sites_to_process)) {
  # get the character string representing the ith site
  current_site <- sites_to_process[i]
  
# convenience string to set file paths for input/output
current_dir <- paste0("data/data_output/", current_site, "/")

ttops <- sf::st_read(paste0(current_dir, "/", current_site, "_ttops/", current_site, "_ttops.shp"))
crowns <- 
  sf::st_read(paste0(current_dir, "/", current_site, "_crowns/", current_site, "_crowns.shp")) %>% 
  mutate(treeID = ttops$treeID)

site_bounds <- 
  sf::st_read(paste0(current_dir, "/", current_site, "_mission-footprint/", current_site, "_site-bounds.geoJSON")) %>% 
  st_transform(st_crs(ttops))

ortho <- 
  raster::brick(paste0(current_dir, 
                       "/3_dsm_ortho/2_mosaic/", 
                       current_site, 
                       "_x3_transparent_mosaic_group1.tif")) %>% 
  raster::crop(site_bounds)

system.time({
  live_or_dead <- 
    cellnumbers(ortho[[1]], crowns) %>% 
    dplyr::mutate(r = raster::extract(ortho[[1]], cell_)) %>% 
    dplyr::mutate(g = raster::extract(ortho[[2]], cell_)) %>% 
    dplyr::mutate(b = raster::extract(ortho[[3]], cell_)) %>% 
    dplyr::mutate(rgi = r/g) %>% 
    dplyr::mutate(rgvi = (r - g) / (r + g)) %>% 
    dplyr::group_by(object_) %>% 
    dplyr::summarize(r_mean = mean(r, na.rm = TRUE),
                     g_mean = mean(g, na.rm = TRUE),
                     b_mean = mean(b, na.rm = TRUE),
                     rgi_mean = mean(rgi, na.rm = TRUE),
                     rgvi_mean = mean(rgvi, na.rm = TRUE),
                     r_hi = quantile(r, probs = 0.75, na.rm = TRUE),
                     g_hi = quantile(g, probs = 0.75, na.rm = TRUE),
                     b_hi = quantile(b, probs = 0.75, na.rm = TRUE),
                     rgi_hi = quantile(rgi, probs = 0.75, na.rm = TRUE),
                     rgvi_hi = quantile(rgvi, probs = 0.75, na.rm = TRUE)) %>% 
    dplyr::rename(treeID = object_)
})

}

crowns_class <- 
  sf::st_read(paste0(current_dir, "/", current_site, "_live-or-dead/", current_site, "_live-or-dead.shp")) %>% 
  mutate(treeID = ttops$treeID)

ttops_extra <- 
  ttops %>% 
  left_join(dplyr::select(as.data.frame(crowns_class), -geometry), by = c("treeID", "layer", "height", "winRadius")) %>% 
  left_join(live_or_dead, by = "treeID")

pca_mean <- prcomp(live_or_dead[, 2:5], scale. = TRUE)
ggplot(as.data.frame(pca_mean$x), aes(x = PC1, y = PC2)) + 
  geom_point()

pca_mean_not_scaled <- prcomp(live_or_dead[, 2:4])
ggplot(as.data.frame(pca_mean_not_scaled$x), aes(x = PC1, y = PC2)) + 
  geom_point()

pca_hi <- prcomp(live_or_dead[, 6:9], scale. = TRUE)
ggplot(as.data.frame(pca_hi$x), aes(x = PC1, y = PC2)) + 
  geom_point()

ggplot(live_or_dead, aes(x = r, y = g)) +
  geom_point()

ggplot(filter(test, !is.na(live)), aes(x = PC1, y = PC2, col = as.factor(live))) +
  geom_point()

pca_full <- prcomp(live_or_dead[, 2:7], scale = TRUE)

ggplot(as.data.frame(pca_full$x), aes(x = PC1, y = PC2)) + 
  geom_point()


ggplot(as.data.frame(pca_full$x), aes(x = PC1, y = PC2, col = k$cluster)) + 
  geom_point()

classified <- 
  ttops_extra %>%
  filter(!is.na(live))

m1 <- glm(live ~ rgi_mean, data = classified, family = "binomial")
m2 <- glm(live ~ rgi_hi + rgi_mean, data = classified, family = "binomial")
m3 <- glm(live ~ rgi_hi + r_mean + g_mean + b_mean, data = classified, family = "binomial")
m4 <- glm(live ~ r_mean + g_mean + b_mean, data = classified, family = "binomial")
m5 <- glm(live ~ r_mean + g_mean + b_mean + rgi_mean, data = classified, family = "binomial")
AIC(m1, m2, m3, m4, m5)

ttops_extra <-
  ttops_extra %>% 
  mutate(predicted_live = plogis(predict(m4, newdata = ., type = "link")))

ggplot(ttops_extra, aes(col = predicted_live)) +
  geom_sf() +
  scale_color_viridis_c()

hist(ttops_extra$predicted_live)

