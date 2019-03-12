# We will use a non-spatial model to determine the stricter priors to use for the spatial model

library(tidyverse)
library(here)
library(brms)
library(future)

# Get the summarized ground data
dd_plot <- readr::read_csv(here::here("data/data_output/ground-data-for-modeling-summarized-by-plot.csv"))
dd_site <- readr::read_csv(here::here("data/data_output/ground-data-for-modeling-summarized-by-site.csv"))

dd_plot <-
  dd_plot %>% 
  dplyr::rename(local_cwd_zscore = plot_cwd_zscore) # This will make the prior modeling easier


# basal area model/plot level----------------------------------------------

# Model using total basal area in each plot
fm_ba_plot_nonspatial <- glm(cbind(dead_count, pipo_count) ~ 
                                local_cwd_zscore*pipo_and_dead_tpha_s*pipo_and_dead_bapha_s +
                                local_cwd_zscore*overall_tpha_s*overall_bapha_s,
                              data = dd_plot,
                              family = "binomial")

priors_df_for_ba_plot_model <-
  bind_cols(
    tibble(coef = names(coef(fm_ba_plot_nonspatial))), 
    as_tibble(summary(fm_ba_plot_nonspatial)$coefficients)) %>% 
  dplyr::rename(estimate = Estimate,
                stderr = `Std. Error`) %>% 
  dplyr::select(coef, estimate, stderr) %>% 
  dplyr::mutate(prior_string = paste0("normal (", estimate, ", ", stderr, ")"))

priors_ba_plot_model <- brms::set_prior(priors_df_for_ba_plot_model$prior_string, coef = priors_df_for_ba_plot_model$coef)
priors_ba_plot_model[priors_ba_plot_model$coef == "(Intercept)", "class"] <- "Intercept"
priors_ba_plot_model[priors_ba_plot_model$coef == "(Intercept)", "coef"] <- ""


# basal area/site level ---------------------------------------------------

# Model using total basal area in each plot
fm_ba_site_nonspatial <- glm(cbind(dead_count, pipo_count) ~ 
                                site_cwd_zscore*pipo_and_dead_tpha_s*pipo_and_dead_bapha_s +
                                site_cwd_zscore*overall_tpha_s*overall_bapha_s,
                              data = dd_site,
                              family = "binomial")

priors_df_for_ba_site_model <-
  bind_cols(
    tibble(coef = names(coef(fm_ba_site_nonspatial))), 
    as_tibble(summary(fm_ba_site_nonspatial)$coefficients)) %>% 
  dplyr::rename(estimate = Estimate,
                stderr = `Std. Error`) %>% 
  dplyr::select(coef, estimate, stderr) %>% 
  dplyr::mutate(prior_string = paste0("normal (", estimate, ", ", stderr, ")"))

priors_ba_site_model <- brms::set_prior(priors_df_for_ba_site_model$prior_string, coef = priors_df_for_ba_site_model$coef)
priors_ba_site_model[priors_ba_site_model$coef == "(Intercept)", "class"] <- "Intercept"
priors_ba_site_model[priors_ba_site_model$coef == "(Intercept)", "coef"] <- ""


# quadratic mean diameter/plot level --------------------------------------

fm_qmd_plot_nonspatial <- glm(cbind(dead_count, pipo_count) ~ 
                        local_cwd_zscore*pipo_and_dead_tpha_s*pipo_and_dead_qmd_s +
                        local_cwd_zscore*overall_tpha_s*overall_qmd_s,
                      data = dd_plot,
                      family = "binomial")

priors_df_for_qmd_plot_model <-
  bind_cols(
    tibble(coef = names(coef(fm_qmd_plot_nonspatial))), 
    as_tibble(summary(fm_qmd_plot_nonspatial)$coefficients)) %>% 
  dplyr::rename(estimate = Estimate,
                stderr = `Std. Error`) %>% 
  dplyr::select(coef, estimate, stderr) %>% 
  dplyr::mutate(prior_string = paste0("normal (", estimate, ", ", stderr, ")"))

priors_qmd_plot_model <- brms::set_prior(priors_df_for_qmd_plot_model$prior_string, class = "b", coef = priors_df_for_qmd_plot_model$coef)
priors_qmd_plot_model[priors_qmd_plot_model$coef == "(Intercept)", "class"] <- "Intercept"
priors_qmd_plot_model[priors_qmd_plot_model$coef == "(Intercept)", "coef"] <- ""



# quadratic mean diameter/site level --------------------------------------

fm_qmd_site_nonspatial <- glm(cbind(dead_count, pipo_count) ~ 
                                site_cwd_zscore*pipo_and_dead_tpha_s*pipo_and_dead_qmd_s +
                                site_cwd_zscore*overall_tpha_s*overall_qmd_s,
                              data = dd_site,
                              family = "binomial")

priors_df_for_qmd_site_model <-
  bind_cols(
    tibble(coef = names(coef(fm_qmd_site_nonspatial))), 
    as_tibble(summary(fm_qmd_site_nonspatial)$coefficients)) %>% 
  dplyr::rename(estimate = Estimate,
                stderr = `Std. Error`) %>% 
  dplyr::select(coef, estimate, stderr) %>% 
  dplyr::mutate(prior_string = paste0("normal (", estimate, ", ", stderr, ")"))

priors_qmd_site_model <- brms::set_prior(priors_df_for_qmd_site_model$prior_string, class = "b", coef = priors_df_for_qmd_site_model$coef)
priors_qmd_site_model[priors_qmd_site_model$coef == "(Intercept)", "class"] <- "Intercept"
priors_qmd_site_model[priors_qmd_site_model$coef == "(Intercept)", "coef"] <- ""

priors_ba_plot_model
priors_ba_site_model
priors_qmd_plot_model
priors_qmd_site_model

