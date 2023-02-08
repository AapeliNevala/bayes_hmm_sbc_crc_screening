library(cmdstanr)
library(tidyverse)
library(data.table)
library(bayesplot)
library(posterior)
library(tidybayes)
library(ggthemes)
library(Epi)
source("R/model_conf_utils.R")
source("R/visualisation_utils.R")
source("R/cutoffs_dev.R")


# cutoffs <- c(55, 65)
# 
# attender_data <- add_cutoffs(attender_data, cutoffs)
# refuser_data <- add_cutoffs(refuser_data, cutoffs)
# control_data <- add_cutoffs(control_data, cutoffs)


stan.data_men_const <- list(
  mu1_prior = c(1,85)*5,
  mu2_prior = c(1,30)*10,
  mu3_prior = c(1,30)*10,
  mu4_prior = c(1,4)*10,
  mu_ref_prior_1 = c(.1, .1),
  mu_ref_prior_2 = c(.1, .1),
  mu_ref_prior_3 = c(.1, .1),
  entry = 30,
  entry_int = 30,
  exit = 80,
  S_c_prior = 0.45,
  S_aa_prior = 0.13,
  S_a_prior = 0.05,
  screen_counts = c(43765, 52890, 49386, 38215, 25627)
)


stan.data_women_const <- list(
  mu1_prior = c(1,120)*5,
  mu2_prior = c(1,30)*10,
  mu3_prior = c(1,30)*10,
  mu4_prior = c(1,4)*10,
  mu_ref_prior_1 = c(.1, .1),
  mu_ref_prior_2 = c(.1, .1),
  mu_ref_prior_3 = c(.1, .1),
  entry = 30,
  entry_int = 30,
  exit = 80,
  S_c_prior = 0.35,
  S_aa_prior = 0.05,
  S_a_prior = 0.025
)

prior_pred_model_const <- cmdstan_model(
  "stan_code/generation_models/prior_pred_const.stan",
  compile = TRUE,
  cpp_options = list(stan_threads = TRUE),
  include_paths = c("./stan_code"),
  dir = file.path(here::here(), "model_files"),
  force_recompile = TRUE
)


pp_m_const <- prior_pred_model_const$sample(
  data = stan.data_men_const,
  threads_per_chain = 1,
  fixed_param = TRUE
)

pp_w_const <- prior_pred_model_const$sample(
  data = stan.data_women_const,
  threads_per_chain = 1,
  fixed_param = TRUE
)



posterior_state_probs_from_draws(pp_m_const, time_up = 76, time_lo = 55) -> prior_pred_m_const
posterior_state_probs_from_draws(pp_w_const, time_up = 76, time_lo = 55) -> prior_pred_w_const

ggpubr::ggarrange(prior_pred_m_const, prior_pred_w_const) -> prior_pred_comb



ggsave("final_plots/prior_pred_m_const.png", prior_pred_m_const, device = "png", dpi = 600, scale = 1.5)
ggsave("final_plots/prior_pred_w_const.png", prior_pred_w_const, device = "png", dpi = 600, scale = 1.5)
ggsave("final_plots/prior_pred_comb.png", prior_pred_comb, device = "png", dpi = 600, scale = 1.5)

