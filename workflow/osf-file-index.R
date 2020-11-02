# get files uploaded already to OSF for local-structure-wpb-severity project
# site, re, x3
# 1, eldo_3k_1, 10465, 1957, good, good
# 2, eldo_3k_2, 10650, 1958, good, good
# 3, eldo_3k_3, 11315, 2008, good, good
# 4, eldo_4k_1, 10495, 1871, good, good
# 5, eldo_4k_2, 10505, 1877, good, good
# NA, eldo_4k_3, NA, NA, NA, NA
# 6, eldo_5k_1, 10785, 1922, good, good
# 7, eldo_5k_2, 10915, 1996, good, good
# 8, eldo_5k_3, 10450, 1902, good, good
# 9, stan_3k_1, 12120, 2064, good, good
# 10, stan_3k_2, 7190, 1266, good, good
# 11, stan_3k_3, 10565, 1945, good, good
# 12, stan_4k_1, 11045, 1909, good, bad (3)
# 13, stan_4k_2, 11500, 1990, good, bad (2)
# NA, stan_4k_3, NA, NA, NA, NA
# 14, stan_5k_1, 10665, 1964, good, bad (2)
# 15, stan_5k_2, 11505, 1968, good, bad (1)
# NA, stan_5k_3, NA, NA, NA, NA
# 16, sier_3k_1, 11761, 2006, bad (9), bad (2)
# 17, sier_3k_2, 10630, 1901, good, bad (4)
# 18, sier_3k_3, 10549, 1915, bad (1), bad (3)
# 19, sier_4k_1, 11765, 1568, good, bad (1)
# 20, sier_4k_2, 10436, 1877, bad (4), bad (7)
# 21, sier_4k_3, 9423, 1688, bad (2), bad (3)
# 22, sier_5k_1, 6762, 1095, bad (3), bad (2)
# 23, sier_5k_2, 9776, 1756, bad (3), bad (7)
# 24, sier_5k_3, 10227, 1867, bad (18), bad (4)
# 25, sequ_4k_1, 10499, 1291, bad (1), bad (2)
# NA, sequ_4k_2, NA, NA, NA, NA
# 26, sequ_4k_3, 10785, 1969, good, bad (1)
# 27, sequ_5k_1, 9403, 1742, bad (2), bad (2)
# 28, sequ_5k_2, 10708, 1485, bad (2), bad (4)
# 29, sequ_5k_3, 11518, 1294, bad (7), bad (7)
# 30, sequ_6k_1, 9675, 1796, good, bad (4)
# 31, sequ_6k_2, 5615, 953, good, bad (2)
# 32, sequ_6k_3, 9470, 1726, good, bad (3)

library(osfr)
library(dplyr)
library(glue)

osfr::osf_auth()

local_structure_wpb_severity_project <- osfr::osf_retrieve_node("3cwf9")

# sites_uploaded <- osfr::osf_ls_files(x = local_structure_wpb_severity_project, path = "data/data_drone/L0/photos", type = "folder", n_max = Inf)

sites <-
  readr::read_csv("index, site, re_n, x3_n, re_check, re_missing, x3_check, x3_missing
                 1, eldo_3k_1, 10465, 1957, good, , good, 
                 2, eldo_3k_2, 10650, 1958, good, , good, 
                 3, eldo_3k_3, 11315, 2008, good, , good,
                 4, eldo_4k_1, 10495, 1871, good, , good,
                 5, eldo_4k_2, 10505, 1877, good, , good,
                 NA, eldo_4k_3, NA, NA, NA, , NA,
                 6, eldo_5k_1, 10785, 1922, good, , good,
                 7, eldo_5k_2, 10915, 1996, good, , good,
                 8, eldo_5k_3, 10450, 1902, good, , good,
                 9, stan_3k_1, 12120, 2064, good, , good,
                 10, stan_3k_2, 7190, 1266, good, , good,
                 11, stan_3k_3, 10565, 1945, good, , good,
                 12, stan_4k_1, 11045, 1909, good, , bad, 3
                 13, stan_4k_2, 11500, 1990, good, , bad, 2
                 NA, stan_4k_3, NA, NA, NA, , NA,
                 14, stan_5k_1, 10665, 1964, good, , bad, 2
                 15, stan_5k_2, 11505, 1968, good, , bad, 1
                 NA, stan_5k_3, NA, NA, NA, , NA, 
                 16, sier_3k_1, 11761, 2006, bad, 9, bad, 2
                 17, sier_3k_2, 10630, 1901, good, , bad, 4
                 18, sier_3k_3, 10549, 1915, bad, 1, bad, 3
                 19, sier_4k_1, 11765, 1568, good, , bad, 1
                 20, sier_4k_2, 10436, 1877, bad, 4, bad, 7
                 21, sier_4k_3, 9423, 1688, bad, 2, bad, 3
                 22, sier_5k_1, 6762, 1095, bad, 3, bad, 2
                 23, sier_5k_2, 9776, 1756, bad, 3, bad, 7
                 24, sier_5k_3, 10227, 1867, bad, 18, bad, 4
                 25, sequ_4k_1, 10499, 1291, bad, 1, bad, 2
                 NA, sequ_4k_2, NA, NA, NA, , NA, 
                 26, sequ_4k_3, 10785, 1969, good, , bad, 1
                 27, sequ_5k_1, 9403, 1742, bad, 2, bad, 2
                 28, sequ_5k_2, 10708, 1485, bad, 2, bad, 4
                 29, sequ_5k_3, 11518, 1294, bad, 7, bad, 7
                 30, sequ_6k_1, 9675, 1796, good, , bad, 4
                 31, sequ_6k_2, 5615, 953, good, , bad, 2
                 32, sequ_6k_3, 9470, 1726, good, , bad, 3")

all_photos <- 
  list.files("data/data_drone/L0/photos-metadata/", full.names = TRUE) %>% 
  purrr::map_dfr(.f = function(x) {
    metadata <- readr::read_csv(x)
  })

incomplete_re <-
  sites %>% 
  dplyr::filter(re_check == "bad")

this_site <- incomplete_re$site[1]
purrr::map(incomplete_re$site, .f = function(this_site) {
  (starttime <- Sys.time())
  osf_re <- osfr::osf_ls_files(x = local_structure_wpb_severity_project, path = glue::glue("data/data_drone/L0/photos/{this_site}"), n_max = Inf, pattern = "re_", verbose = TRUE)
  
  readr::write_rds(osf_re, file = glue::glue("data/data_output/osf-files/{this_site}_osf-files.rds"))
  
  osf_re_text <-
    osf_re %>% 
    dplyr::mutate(site = this_site) %>% 
    dplyr::select(site, name, id)
  
  readr::write_csv(osf_re_text, file = glue::glue("data/data_output/osf-files/{this_site}_osf-files.csv"))
  (difftime(Sys.time(), starttime, units = "mins"))
})