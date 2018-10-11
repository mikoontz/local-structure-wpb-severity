# This script reworks the raw ground data sent to me by Leif to make it
# easier to manipulate in R

library(tidyverse)
library(sf)
library(lme4)

summarized_hcc_data <- read_csv("data/data_output/summarized-non-spatial-site-data.csv")

d <- 
  read_csv("data/ground-data.csv") %>% 
  rename(leif_name = `Plot #`,
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
  select(tree_id,
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
  filter(!is.na(tree_id))

d <- 
  d %>% 
  mutate(wpb = tolower(wpb)) %>% 
  mutate(wpb = ifelse(wpb == "wbp", yes = "wpb", no = wpb)) %>% 
  mutate(wpb = ifelse(str_detect(wpb, pattern = "wpb"), yes = 1, no = 0)) %>% 
  mutate(dbh = dbh_in * 2.54) %>% 
  mutate(height = height_ft * 12 * 2.54 / 100) %>% 
  mutate(dist = dist_ft * 12 * 2.54 / 100) %>% 
  mutate(live = ifelse(live == "L", yes = 1, no = 0)) %>% 
  mutate(forest = tolower(forest)) %>% 
  mutate(forest_code = case_when(forest == "eldorado" ~ "eldo",
                                  forest == "sequoia" ~ "sequ",
                                  forest == "sierra" ~ "sier",
                                  forest == "stanislaus" ~ "stan")) %>% 
  select(-height_ft, -dist_ft) %>% 
  mutate(elev_code = paste0(elev / 1000, "k")) %>% 
  unite(col = site, forest_code, elev_code, rep, remove = FALSE) %>% 
  unite(col = plot, site, plot_id, remove = FALSE) %>% 
  unite(col = tree, plot, tree_id, remove = FALSE)
  

glimpse(d)

