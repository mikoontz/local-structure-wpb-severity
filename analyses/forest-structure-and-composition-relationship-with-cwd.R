# Purpose: assess any forest structure/composition differences along the CWD gradient in both ground and aerial data that might offer an alternative explanation for the strong CWD/host size interaction we found in our analysis

library(tidyverse)

site_data <- read_csv("analyses/analyses_output/summarized-non-spatial-site-data.csv")

site_data <-
  site_data %>% 
  dplyr::mutate(air_n_dead = air_n_total * air_prop_mortality,
                ground_n_dead = ground_n_total * ground_prop_mortality,
                air_mortality_rate_tphapyr = air_tpha_total * air_prop_mortality,
                ground_mortality_rate_tphapyr = ground_tpha_total * ground_prop_mortality)


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


# relationship with dead PIPO qmd -----------------------------------------

ggplot(site_data, aes(x = site_cwd_zscore, y = air_qmd_dead / ground_qmd_dead_pipo)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  scale_color_viridis_d() +
  theme_bw()

pipo_mortality_qmd <- 
  site_data %>% 
  dplyr::select(ground_qmd_pipo_dead, ground_prop_mortality_pipo, air_prop_mortality) %>% 
  dplyr::rename(qmd = ground_qmd_dead) %>% 
  dplyr::rename(ground = starts_with("ground"), air = starts_with("air")) %>% 
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
  labs(x = "Site CWD (z score)",
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
  labs(x = "Site CWD (z score)",
       y = "Proportion of *dead* trees that are nonhosts from ground data")


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
  theme_bw()

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
  theme_bw()

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
