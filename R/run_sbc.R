source("R/5sm_data_generating_functions.R")
library(parallel)
library(cmdstanr)



# Fix priors --------------------------------------------------------------

priors_const_m <- list(
  mu1_prior = c(1,85)*5,
  mu2_prior = c(1,40)*10,
  mu3_prior = c(1,40)*10,
  mu4_prior = c(1,4)*10,
  mu_ref_prior_1 = c(0.1, .1),
  mu_ref_prior_2 = c(0.1, .1),
  mu_ref_prior_3 = c(0.1, .1),
  entry = 30,
  entry_int = 30,
  exit = 80
)

priors_const_w <- list(
  mu1_prior = c(1,120)*5,
  mu2_prior = c(1,40)*10,
  mu3_prior = c(1,40)*10,
  mu4_prior = c(1,4)*10,
  mu_ref_prior_1 = c(0.1, .1),
  mu_ref_prior_2 = c(0.1, .1),
  mu_ref_prior_3 = c(0.1, .1),
  entry = 30,
  entry_int = 30,
  exit = 80
)




# Add other configurations ------------------------------------------------



sbc_init_data_const_m <- list(
  N = 5000,
  S_c_prior = 0.45,
  S_aa_prior = 0.13,
  S_a_prior = 0.05,
  death_rates_and_cutoffs = death_rates_and_cutoffs,
  "priors" = priors_const_m,
  cutoffs = NULL,
  hide_perc = 0.3,
  parametrization = "constant"
)


sbc_init_data_const_w <- list(
  N = 5000,
  S_c_prior = 0.35,
  S_aa_prior = 0.05,
  S_a_prior = 0.025,
  death_rates_and_cutoffs = death_rates_and_cutoffs,
  "priors" = priors_const_w,
  cutoffs = NULL,
  hide_perc = 0.3,
  parametrization = "constant"
)

# Build models ------------------------------------------------------------



five_sm_const_rate <- cmdstan_model(
  "stan_code/5sm/5_sm_constant_rates.stan",
  compile = TRUE,
  cpp_options = list(stan_threads = TRUE),
  include_paths = c("./stan_code/5sm"),
  dir = file.path(here::here(), "model_files"),
  force_recompile = F
)


cl_1 <- makeCluster(detectCores())
clusterExport(cl_1, varlist = ls())
clusterEvalQ(cl_1, expr = (source("R/5sm_data_generating_functions.R")))


run_one_sim(50,  sbc_init_data_tv_m, five_sm_tv_rate)

# Women, constant rates model ---------------------------------------------

sbc_results_const_advi_w <- parLapply(cl_1, 1:500, function(i) {
  try(run_one_sim_advi(50,  sbc_init_data_const_w, five_sm_const_rate), TRUE)
})

saveRDS(sbc_results_const_advi_w, "data/sbc_results_constant_w.RDS")


sbc_results_const_hmc_w <- parLapply(cl_1, 1:500, function(i) {
  try(run_one_sim_hmc(50,  sbc_init_data_const_w, five_sm_const_rate), TRUE)
})

saveRDS(sbc_results_const_advi_w, "data/sbc_results_constant_w_hmc.RDS")


# Men, constant rates model -----------------------------------------------

sbc_results_const_advi_m <- parLapply(cl_1, 1:500, function(i) {
  try(run_one_sim_advi(50,  sbc_init_data_const_m, five_sm_const_rate), TRUE)
})

saveRDS(sbc_results_const_advi_m, "data/sbc_results_constant_m.RDS")






