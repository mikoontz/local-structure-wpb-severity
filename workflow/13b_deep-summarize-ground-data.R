# Purpose: summarized the ground plot data

library(sf)
library(tidyverse)
library(lubridate)
library(raster)
library(here)

if(!file.exists("data/data_output/formatted-ground-data.csv")) {
  source(here::here("workflow/11_format-ground-data.R"))
}
d <- readr::read_csv("data/data_output/formatted-ground-data.csv")

cwd <- raster::raster(here::here("data/features/cwd1981_2010_ave_HST_1550861123/cwd1981_2010_ave_HST_1550861123.tif"))
# The .prj file doesn't seem to be reading in properly with the .tif, but we can look at it in a text editor and see that it is EPSG3310
crs(cwd) <- st_crs(3310)$proj4string

if(!file.exists(here::here("data", "data_output", "sierra-nevada-pipo-cwd.gpkg"))) {
  source(here::here("workflow", "03_extract-cwd-from-locations.R"))
}

sn_pipo <- sf::st_read(here::here("data", "data_output", "sierra-nevada-pipo-cwd.gpkg"))

mean_cwd_sn_pipo <- mean(sn_pipo$cwd, na.rm = TRUE)
sd_cwd_sn_pipo <- sd(sn_pipo$cwd, na.rm = TRUE)

# cwd for individual plots ------------------------------------------------

plot_centers <- 
  st_read(here::here("data", "data_raw", "plot-centers_ground-gps-measured.kml")) %>% 
  st_transform(3310) %>% 
  st_zm() %>% 
  tidyr::separate(col = Name, into = c("forest", "elevation_band", "rep", "nickname", "plot_id"), sep = "_") %>% 
  dplyr::select(-Description) %>% 
  dplyr::mutate(elevation_band = as.numeric(as.character(elevation_band))) %>% 
  dplyr::mutate(forest = substr(forest, start = 1, stop = 4)) %>% 
  dplyr::mutate(elevation_band = paste0(substr(elevation_band, start = 1, stop = 1), "k")) %>%
  dplyr::mutate(plot_id = as.numeric(plot_id),
                rep = as.numeric(rep)) %>% 
  dplyr::mutate(plot_id = plot_id - (5 * (rep - 1))) %>% 
  dplyr::mutate(site = paste(forest, elevation_band, rep, sep = "_"),
                plot = paste(site, plot_id, sep = "_")) %>% 
  dplyr::mutate(plot_cwd = raster::extract(cwd, ., method = "bilinear")) %>% 
  dplyr::mutate(plot_cwd_zscore = (plot_cwd - mean_cwd_sn_pipo) / sd_cwd_sn_pipo) %>% 
  dplyr::select(site, plot, plot_cwd, plot_cwd_zscore)

plot_centers


# summarize data for each plot --------------------------------------------
center_param <- TRUE
scale_param <- TRUE

dd_plot <-
  d %>% 
  dplyr::mutate(ba = (dbh / 2)^2 * pi) %>% 
  group_by(forest, elev, rep, site, plot) %>% 
  summarize(live_count = length(which(live == 1)),
            dead_count = length(which(live == 0)),
            pipo_count = length(which(species == "PIPO" & live == 1)),
            abco_count = length(which(species == "ABCO" & live == 1)),
            cade_count = length(which(species == "CADE" & live == 1)),
            pipo_count_dead = length(which(species == "PIPO" & live == 0)),
            abco_count_dead = length(which(species == "ABCO" & live == 0)),
            cade_count_dead = length(which(species == "CADE" & live == 0)),
            nonhost_count_live = length(which(species != "PIPO" & live == 1)),
            nonhost_count_dead = length(which(species != "PIPO" & live == 0)),
            pipo_and_dead_count = pipo_count + dead_count,
            total_count = n(),
            n_fallen = length(which(!is.na(year_fall))),
            pct_fallen = n_fallen / dead_count,
            qmd_fallen = sqrt(sum(dbh[which(!is.na(year_fall))]^2)),
            n_fallen_pipo = length(which(!is.na(year_fall) & species == "PIPO")),
            pct_fallen_pipo = n_fallen_pipo / pipo_count_dead,
            qmd_fallen_pipo = sqrt(sum(dbh[which(!is.na(year_fall) & species == "PIPO")]^2)),
            live_tpha = live_count / 0.0404686, # there are 0.0404686 hectares in 0.1 acres (the total surveyed area per plot)
            dead_tpha = dead_count / 0.0404686,
            pipo_tpha = pipo_count / 0.0404686,
            abco_tpha = abco_count / 0.0404686,
            cade_tpha = cade_count / 0.0404686,
            pipo_dead_tpha = pipo_count_dead / 0.0404686,
            abco_dead_tpha = abco_count_dead / 0.0404686,
            cade_dead_tpha = cade_count_dead / 0.0404686,
            nonhost_live_tpha = nonhost_count_live / 0.0404686,
            nonhost_dead_tpha = nonhost_count_dead / 0.0404686,
            pipo_and_dead_tpha = pipo_and_dead_count / 0.0404686,
            overall_tpha = total_count / 0.0404686,
            live_ba = sum(ba[live == 1]),
            dead_ba = sum(ba[live == 0]),
            pipo_ba = sum(ba[species == "PIPO" & live == 1]),
            pipo_and_dead_ba = pipo_ba + dead_ba,
            total_ba = sum(ba),
            live_bapha = live_ba / 0.0404686,
            dead_bapha = dead_ba / 0.0404686,
            pipo_bapha = pipo_ba / 0.0404686,
            pipo_and_dead_bapha = pipo_and_dead_ba / 0.0404686,
            overall_bapha = total_ba / 0.0404686,
            live_qmd = sqrt(sum(dbh[live == 1]^2) / live_count),
            dead_qmd = sqrt(sum(dbh[live == 0]^2) / dead_count),
            pipo_qmd = sqrt(sum(dbh[species == "PIPO" & live == 1]^2) / pipo_count),
            pipo_live_and_dead_qmd = sqrt(sum(dbh[species == "PIPO"]^2) / (pipo_count + pipo_count_dead)),
            pipo_dead_qmd = sqrt(sum(dbh[live == 0 & species == "PIPO"]^2) / pipo_count_dead),
            pipo_and_dead_qmd = sqrt(sum(dbh[live == 0 | (species == "PIPO" & live == 1)]^2) / pipo_and_dead_count),
            overall_qmd = sqrt(sum(dbh^2) / total_count),
            live_sdi_ac = live_count * 10 * (live_qmd / 25.4)^1.77, # multiply counts by 10 to get trees per acre
            pipo_and_dead_sdi_ac = pipo_and_dead_count * 10 * (pipo_and_dead_qmd / 25.4)^1.77,
            overall_sdi_ac = total_count * 10 * (overall_qmd / 25.4)^1.77) %>% 
  dplyr::left_join(plot_centers) %>% 
  dplyr::mutate(pipo_and_dead_tpha_s = scale(pipo_and_dead_tpha, center = center_param, scale = scale_param),
                overall_tpha_s = scale(overall_tpha, center = center_param, scale = scale_param),
                pipo_and_dead_bapha_s = scale(pipo_and_dead_bapha, center = center_param, scale = scale_param),
                overall_bapha_s = scale(overall_bapha, center = center_param, scale = scale_param),
                pipo_and_dead_qmd_s = scale(pipo_and_dead_qmd, center = center_param, scale = scale_param),
                overall_qmd_s = scale(overall_qmd, center = center_param, scale = scale_param),
                live_sdi_ac_s = scale(live_sdi_ac, center = center_param, scale = scale_param),
                pipo_and_dead_sdi_ac_s = scale(pipo_and_dead_sdi_ac, center = center_param, scale = scale_param),
                overall_sdi_ac_s = scale(overall_sdi_ac, center = center_param, scale = scale_param)) %>% 
  ungroup() %>% 
  st_as_sf() %>% 
  dplyr::mutate(x_3310 = st_coordinates(.)[, 1],
                y_3310 = st_coordinates(.)[, 2]) %>% 
  st_drop_geometry()

# aggregating up to the site level ----------------------------------------

site_centers <- 
  plot_centers %>% 
  group_by(site) %>% 
  summarize() %>% 
  st_centroid() %>% 
  dplyr::mutate(site_cwd = raster::extract(cwd, ., method = "bilinear")) %>% 
  dplyr::mutate(site_cwd_zscore = (site_cwd - mean_cwd_sn_pipo) / sd_cwd_sn_pipo)

dd_site <-
  d %>% 
  dplyr::mutate(ba = (dbh / 2)^2 * pi) %>% 
  group_by(forest, elev, rep, site) %>% 
  summarize(live_count = length(which(live == 1)),
            dead_count = length(which(live == 0)),
            pipo_count = length(which(species == "PIPO" & live == 1)),
            abco_count = length(which(species == "ABCO")),
            cade_count = length(which(species == "CADE")),
            pipo_count_dead = length(which(species == "PIPO" & live == 0)),
            abco_count_dead = length(which(species == "ABCO" & live == 0)),
            cade_count_dead = length(which(species == "CADE" & live == 0)),
            nonhost_count_live = length(which(species != "PIPO" & live == 1)),
            nonhost_count_dead = length(which(species != "PIPO" & live == 0)),
            pipo_and_dead_count = pipo_count + dead_count,
            total_count = n(),
            n_fallen = length(which(!is.na(year_fall))),
            pct_fallen = n_fallen / dead_count,
            qmd_fallen = sqrt(sum(dbh[which(!is.na(year_fall))]^2)),
            n_fallen_pipo = length(which(!is.na(year_fall) & species == "PIPO")),
            pct_fallen_pipo = n_fallen_pipo / pipo_count_dead,
            qmd_fallen_pipo = sqrt(sum(dbh[which(!is.na(year_fall) & species == "PIPO")]^2)),
            live_tpha = live_count / 0.202343, # there are 0.202343 hectares in 0.5 acres (the total surveyed area per site)
            dead_tpha = dead_count / 0.202343,
            pipo_tpha = pipo_count / 0.202343,
            abco_tpha = abco_count / 0.202343,
            cade_tpha = cade_count / 0.202343,
            pipo_dead_tpha = pipo_count_dead / 0.0404686,
            abco_dead_tpha = abco_count_dead / 0.0404686,
            cade_dead_tpha = cade_count_dead / 0.0404686,
            nonhost_live_tpha = nonhost_count_live / 0.0404686,
            nonhost_dead_tpha = nonhost_count_dead / 0.0404686,
            pipo_and_dead_tpha = pipo_and_dead_count / 0.202343,
            overall_tpha = total_count / 0.202343,
            live_ba = sum(ba[live == 1]),
            dead_ba = sum(ba[live == 0]),
            pipo_ba = sum(ba[species == "PIPO" & live == 1]),
            pipo_and_dead_ba = pipo_ba + dead_ba,
            total_ba = sum(ba),
            live_bapha = live_ba / 0.202343,
            dead_bapha = dead_ba / 0.202343,
            pipo_bapha = pipo_ba / 0.202343,
            pipo_and_dead_bapha = pipo_and_dead_ba / 0.202343,
            overall_bapha = total_ba / 0.202343,
            live_qmd = sqrt(sum(dbh[live == 1]^2) / live_count),
            dead_qmd = sqrt(sum(dbh[live == 0]^2) / dead_count),
            pipo_qmd = sqrt(sum(dbh[species == "PIPO" & live == 1]^2) / pipo_count),
            pipo_live_and_dead_qmd = sqrt(sum(dbh[species == "PIPO"]^2) / (pipo_count + pipo_count_dead)),
            pipo_dead_qmd = sqrt(sum(dbh[live == 0 & species == "PIPO"]^2) / pipo_count_dead),
            pipo_and_dead_qmd = sqrt(sum(dbh[live == 0 | (species == "PIPO" & live == 1)]^2) / pipo_and_dead_count),
            overall_qmd = sqrt(sum(dbh^2) / total_count),
            live_sdi_ac = live_count * 2 * (live_qmd / 25.4)^1.77, # Multiply live count by 2 because that is equivalent to trees per acre (because 0.5 acres surveyed at each site)
            pipo_and_dead_sdi_ac = pipo_and_dead_count * 2 * (pipo_and_dead_qmd / 25.4)^1.77,
            overall_sdi_ac = total_count * 2 * (overall_qmd / 25.4)^1.77) %>% 
  dplyr::left_join(site_centers) %>% 
  ungroup() %>% 
  dplyr::mutate(pipo_and_dead_tpha_s = scale(pipo_and_dead_tpha, center = center_param, scale = scale_param),
                overall_tpha_s = scale(overall_tpha, center = center_param, scale = scale_param),
                pipo_and_dead_bapha_s = scale(pipo_and_dead_bapha, center = center_param, scale = scale_param),
                overall_bapha_s = scale(overall_bapha, center = center_param, scale = scale_param),
                pipo_and_dead_qmd_s = scale(pipo_and_dead_qmd, center = center_param, scale = scale_param),
                overall_qmd_s = scale(overall_qmd, center = center_param, scale = scale_param),
                live_sdi_ac_s = scale(live_sdi_ac, center = center_param, scale = scale_param),
                pipo_and_dead_sdi_ac_s = scale(pipo_and_dead_sdi_ac, center = center_param, scale = scale_param),
                overall_sdi_ac_s = scale(overall_sdi_ac, center = center_param, scale = scale_param)) %>% 
  st_as_sf() %>% 
  dplyr::mutate(x_3310 = st_coordinates(.)[, 1],
                y_3310 = st_coordinates(.)[, 2]) %>% 
  st_drop_geometry()

  
readr::write_csv(dd_site, here::here("analyses", "analyses_output", "deep-ground-tree-summary-by-site.csv"))
readr::write_csv(dd_plot, here::here("analyses", "analyses_output", "deep-ground-tree-summary-by-plot.csv"))


