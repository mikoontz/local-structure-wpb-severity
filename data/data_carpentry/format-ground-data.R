# This script reworks the raw ground data sent to me by Leif to make it
# easier to manipulate in R

library(tidyverse)
library(sf)
library(here)

d <- 
  readr::read_csv(here::here("data/data_raw/ground-data.csv")) %>% 
  dplyr::rename(leif_name = `Plot #`,
         forest = Forest,
         elev = `Elev. Band`,
         rep = Block,
         plot_id = Plot,
         tree_id = `Tree tag #`,
         species = Species,
         dist_ft = `Distance (ft)`,
         azm = Azm,
         dbh_in = `DBH (in)`,
         height_ft = `HT (ft)`,
         live = `Live/dead`,
         year_dead = `Year death (time of attack or time of injury/other cause)`,
         year_fall = `Year Fall (initial)`,
         wpb = `BB/Cause`) %>% 
  dplyr::select(tree_id,
         forest,
         elev,
         rep,
         plot_id,
         species,
         dist_ft,
         azm,
         dbh_in,
         height_ft,
         live,
         year_dead,
         year_fall,
         wpb) %>% 
  dplyr::filter(!is.na(tree_id))

d <- 
  d %>% 
  dplyr::mutate(wpb = tolower(wpb)) %>% 
  dplyr::mutate(wpb = ifelse(wpb == "wbp", yes = "wpb", no = wpb)) %>% 
  dplyr::mutate(wpb = ifelse(str_detect(wpb, pattern = "wpb"), yes = 1, no = 0)) %>% 
  dplyr::mutate(dbh = dbh_in * 2.54) %>% 
  dplyr::mutate(height = height_ft * 12 * 2.54 / 100) %>% 
  dplyr::mutate(dist = dist_ft * 12 * 2.54 / 100) %>% 
  dplyr::mutate(live = ifelse(live == "L", yes = 1, no = 0)) %>% 
  dplyr::mutate(forest = tolower(forest)) %>% 
  dplyr::mutate(forest_code = case_when(forest == "eldorado" ~ "eldo",
                                  forest == "sequoia" ~ "sequ",
                                  forest == "sierra" ~ "sier",
                                  forest == "stanislaus" ~ "stan")) %>% 
  dplyr::mutate(rep = ifelse(forest == "eldorado" & elev == 4000 & rep %in% c(4, 5), yes = 3, no = rep)) %>% 
  dplyr::select(-height_ft, -dist_ft) %>% 
  dplyr::mutate(elev_code = paste0(elev / 1000, "k")) %>% 
  dplyr::mutate(plot_id = plot_id - ((rep - 1) * 5)) %>% 
  tidyr::unite(col = site, forest_code, elev_code, rep, remove = FALSE) %>% 
  tidyr::unite(col = plot, site, plot_id, remove = FALSE) %>% 
  tidyr::unite(col = tree, plot, tree_id, remove = FALSE)
  