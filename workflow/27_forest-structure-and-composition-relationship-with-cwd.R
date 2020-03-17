# Purpose: assess any forest structure/composition differences along the CWD gradient in both ground and aerial data that might offer an alternative explanation for the strong CWD/host size interaction we found in our analysis

library(tidyverse)
library(ggrepel)

site_data <- read_csv("analyses/analyses_output/summarized-non-spatial-site-data.csv")

site_data <-
  site_data %>% 
  dplyr::mutate(air_n_dead = air_n_total * air_prop_mortality,
                ground_n_dead = ground_n_total * ground_prop_mortality,
                air_mortality_rate_tphapyr = air_tpha_total * air_prop_mortality,
                ground_mortality_rate_tphapyr = ground_tpha_total * ground_prop_mortality) %>% 
  dplyr::mutate(forest = factor(substr(site, start = 1, stop = 4), levels = c("eldo", "stan", "sier", "sequ")),
                elev = factor(substr(site, start = 6, stop = 7), levels = c("3k", "4k", "5k", "6k")))


# the overall pattern that we analyzed; CWD predicting pipo mortal; QMD predicting mortality --------

pipo_mortality_cwd <- 
  site_data %>% 
  dplyr::select(site_cwd_zscore, ground_prop_mortality_pipo, air_prop_mortality) %>% 
  dplyr::rename(ground = starts_with("ground"), air = starts_with("air")) %>% 
  pivot_longer(cols = c("ground", "air"), names_to = "ground_or_air", values_to = "prop_pipo_mortality")
  
ggplot(pipo_mortality_cwd, aes(x = site_cwd_zscore, y = prop_pipo_mortality, color = ground_or_air)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  scale_color_viridis_d() +
  theme_bw()

ggplot(site_data, aes(x = site_cwd_zscore, y = ground_prop_mortality_pipo / air_prop_mortality )) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  scale_color_viridis_d() +
  theme_bw()

# Where do the sites fall out on the CWD gradient?

ggplot(site_data, aes(x = forest, y = site_cwd_zscore, color = elev)) + 
  geom_point() +
  scale_color_viridis_d() +
  geom_label_repel(aes(label = site)) +
  labs(x = "Forest",
       y = "Site CWD z-score") +
  ggtitle("Generally higher CWD for lower latitude and lower elevations")

ggplot(site_data, aes(x = elev, y = site_cwd_zscore, color = forest)) + 
  geom_point() +
  scale_color_viridis_d() +
  geom_label_repel(aes(label = site))


# overestimation of mortality in ground plots -----------------------------

ggplot(site_data, aes(x = site_cwd_zscore, y = air_qmd_dead / ground_qmd_dead_pipo)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  # geom_smooth() +
  scale_color_viridis_d() +
  theme_bw() +
  geom_hline(yintercept = 1, color = "red") +
  labs(x = "Site CWD z-score",
       y = "Ratio of aerially-detected mortality to ground-detected mortality") +
  ggtitle("Almost always (31 out of 32 sites) less mortality in surrounding forest compared to in plots\nGreater overestimation of mortality in cool, wet sites")


# relationship with dead PIPO qmd -----------------------------------------

pipo_mortality_qmd <- 
  site_data %>% 
  dplyr::select(air_qmd_pipo_and_dead, ground_prop_mortality_pipo, air_prop_mortality) %>% 
  dplyr::rename(qmd = air_qmd_pipo_and_dead) %>% 
  dplyr::rename(ground = ground_prop_mortality_pipo, air = air_prop_mortality) %>% 
  pivot_longer(cols = c("ground", "air"), names_to = "ground_or_air", values_to = "prop_pipo_mortality")

ggplot(pipo_mortality_qmd, aes(x = qmd, y = prop_pipo_mortality, color = ground_or_air)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  scale_color_viridis_d() +
  theme_bw()

# Did nonhost mortality vary across the CWD gradient? ------------------------

# Site CWD predicting the overall proportion of trees that are dead non-pipos
# slightly increasing trend across the CWD gradient, implying that more non-hosts are dying in hot/dry sites
ggplot(site_data, aes(x = site_cwd_zscore, y = ground_prop_mortality_nonhost)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Site CWD z-score",
       y = "Proportion of all trees that are dead nonhosts")

# Site CWD predicting the proportion of *dead* trees that are non-pipos
# We can think of this as an estimate of our classification error (since we classified all dead trees as PIPO and
# thus this plot would show 0 across the whole CWD gradient)
# Probably a better comparison to the aerial data
# No real trend here, suggesting that the fraction of dead trees that are non-hosts does not change across the
# CWD gradient. This suggests that a CWD by host/nonhost interaction (masked by our classification errors for dead 
# trees) is probably not the cause of the CWD by host size interaction effect on probability of mortality
ggplot(site_data, aes(x = site_cwd_zscore, y = ground_prop_mortality_nonhost / ground_prop_mortality)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Site CWD z-score",
       y = "Proportion of *dead* trees that are nonhosts from ground data") +
  ggtitle("No trend or perhaps *greater* species misclassification in hot, dry sites not cool, wet sites")



# host tree size across the CWD gradient ---------------------------------------

pipo_and_dead_qmd <- 
  site_data %>% 
  dplyr::select(site_cwd_zscore, ground_qmd_pipo_and_dead, air_qmd_pipo_and_dead) %>% 
  dplyr::rename(ground = starts_with("ground"), air = starts_with("air")) %>% 
  pivot_longer(cols = c("ground", "air"), names_to = "ground_or_air", values_to = "qmd_pipo_and_dead")

ggplot(pipo_and_dead_qmd, aes(x = site_cwd_zscore, y = qmd_pipo_and_dead, color = ground_or_air)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  scale_color_viridis_d() +
  theme_bw() +
  labs(x = "Site CWD z-score",
       y = "QMD at each site using only living PIPO and *all* dead trees") +
  ggtitle("Similar trend between air & ground surveys QMD of living PIPO and all dead trees across CWD gradient")

pipo_qmd <- 
  site_data %>% 
  dplyr::select(site_cwd_zscore, ground_qmd_pipo, air_qmd_pipo) %>% 
  dplyr::rename(ground = starts_with("ground"), air = starts_with("air")) %>% 
  pivot_longer(cols = c("ground", "air"), names_to = "ground_or_air", values_to = "qmd_pipo")

ggplot(pipo_qmd, aes(x = site_cwd_zscore, y = qmd_pipo, color = ground_or_air)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  scale_color_viridis_d() +
  theme_bw()

# The QMD of dead trees doesn't change across the gradient;
# This also suggests that the host QMD by CWD interaction on Pr(host mortality) is not
# driven by our classification error.

dead_qmd <- 
  site_data %>% 
  dplyr::select(site_cwd_zscore, ground_qmd_dead, air_qmd_dead) %>% 
  dplyr::rename(ground = starts_with("ground"), air = starts_with("air")) %>% 
  pivot_longer(cols = c("ground", "air"), names_to = "ground_or_air", values_to = "qmd_dead")

ggplot(dead_qmd, aes(x = site_cwd_zscore, y = qmd_dead, color = ground_or_air)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  scale_color_viridis_d() +
  theme_bw() +
  labs(x = "Site CWD z-score",
       y = "QMD of dead trees per site") +
  ggtitle("No difference in QMD vs. CWD trend between aerial and ground surveys")

# How did we classify forest composition across CWD gradient? ----------------

prop_nonhost_live <-
  site_data %>% 
  dplyr::select(site_cwd_zscore, air_prop_nonhost_live, ground_prop_nonhost_live) %>% 
  dplyr::rename(ground = starts_with("ground"), air = starts_with("air")) %>% 
  pivot_longer(cols = c("ground", "air"), names_to = "ground_or_air", values_to = "prop_nonhost_live")

ggplot(prop_nonhost_live, aes(x = site_cwd_zscore, y = prop_nonhost_live, color = ground_or_air)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_viridis_d() + 
  theme_bw() +
  labs(x = "Site CWD (z score)",
       y = "Proportion of live trees that are nonhosts")

prop_nonhost_all <-
  site_data %>% 
  dplyr::select(site_cwd_zscore, air_prop_nonhost_all, ground_prop_nonhost_all) %>% 
  dplyr::rename(ground = starts_with("ground"), air = starts_with("air")) %>% 
  pivot_longer(cols = c("ground", "air"), names_to = "ground_or_air", values_to = "prop_nonhost_all")

ggplot(prop_nonhost_all, aes(x = site_cwd_zscore, y = prop_nonhost_all, color = ground_or_air)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_viridis_d() + 
  theme_bw() +
  labs(x = "Site CWD (z score)",
       y = "Proportion of all trees that are nonhosts")

ggplot(site_data, aes(x = site_cwd_zscore, y = air_prop_nonhost_all)) +
  geom_point() +
  geom_smooth(method = "lm")


# size relationship with ground pipo qmd and mortality? -------------------

ggplot(site_data, aes(x = ground_qmd_pipo_live_and_dead, y = ground_prop_mortality_pipo)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Ponderosa (live and dead) QMD from ground data",
       y = "Proportion of ponderosa mortality from ground data") +
  ggtitle("Strong, possibly negative relationship between ponderosa QMD and proportion of ponderosa mortality")



ground_trees <- read_csv("data/data_output/ground-data-for-modeling-summarized-by-plot.csv")

summary_table_by_plot_ground <-
  ground_trees %>% 
  dplyr::mutate(ground_n_total = live_count + dead_count,
                ground_tpha_total = ground_n_total / 0.0404686,
                ground_tpha_abco = abco_count / 0.0404686,
                ground_tpha_cade = cade_count / 0.0404686,
                ground_tpha_pipo_live = pipo_count / 0.0404686,
                ground_tpha_pipo_and_dead = pipo_and_dead_count / 0.0404686,
                ground_prop_mortality = dead_count / ground_n_total,
                ground_prop_mortality_pipo = pipo_count_dead / ground_n_total,
                ground_prop_mortality_nonhost = nonhost_count_dead / ground_n_total,
                ground_prop_mortality_abco = abco_count_dead / ground_n_total,
                ground_prop_pipo_live = pipo_count / ground_n_total,
                ground_prop_pipo_and_dead = pipo_and_dead_count / ground_n_total,
                ground_prop_abco = abco_count / ground_n_total,
                ground_prop_cade = cade_count / ground_n_total,
                ground_prop_nonhost_live = 1 - ground_prop_pipo_live,
                ground_prop_nonhost_all = 1 - ground_prop_pipo_and_dead) %>% 
  dplyr::rename(ground_qmd_pipo = pipo_qmd,
                ground_qmd_pipo_live_and_dead = pipo_live_and_dead_qmd,
                ground_qmd_pipo_and_dead = pipo_and_dead_qmd,
                ground_qmd_overall = overall_qmd,
                ground_qmd_dead = dead_qmd,
                ground_qmd_dead_pipo = pipo_dead_qmd,
                ground_qmd_live = live_qmd) 

ggplot(summary_table_by_plot_ground, aes(x = ground_qmd_pipo_live_and_dead, y = ground_prop_mortality_pipo)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Ponderosa (live and dead) QMD from ground data",
       y = "Proportion of ponderosa mortality from ground data") +
  ggtitle("No, or perhaps negative relationship between ponderosa QMD and proportion of ponderosa mortality\nwhen aggregated to plot")


  
# proportion of live PIPO for ground vs. air measurements -----------------

ggplot(site_data, aes(x = site_cwd_zscore, y = ground_prop_pipo_live)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(site_data, aes(x = site_cwd_zscore, y = air_prop_pipo_live)) +
  geom_point() +
  geom_smooth(method = "lm")

# ground measurements: live pipo + dead pipo ------------------------------

ggplot(site_data, aes(x = site_cwd_zscore, y = ground_prop_pipo_live + ground_prop_mortality_pipo)) +
  geom_point() +
  geom_smooth(method = "lm")


ggplot(site_data, aes(x = site_cwd_zscore, y = air_prop_abco)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(site_data, aes(x = site_cwd_zscore, y = air_prop_cade)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(site_data, aes(x = site_cwd_zscore, y = air_prop_nonhost_live)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(site_data, aes(x = site_cwd_zscore, y = air_prop_nonhost_all)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(site_data, aes(x = site_cwd_zscore, y = ground_prop_nonhost_live)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(site_data, aes(x = site_cwd_zscore, y = ground_qmd_pipo_and_dead)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(site_data, aes(x = site_cwd_zscore, y = ground_qmd_pipo_and_dead)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(site_data, aes(x = site_cwd_zscore, y = ground_qmd_overall)) +
  geom_point() +
  geom_smooth(method = "lm")
