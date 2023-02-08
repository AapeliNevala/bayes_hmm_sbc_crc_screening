library(tidyverse)
library(dplyr)
library(tidyr)
library(data.table)
library(cmdstanr)
library(latex2exp)
source("R/model_conf_utils.R")
setwd(here::here())


# Load and select proper population mortality data --------------------
popmort <- read.csv2("data/popmorts.csv")
popmort <- tibble(popmort)
popmort %>% 
  filter(per > 2000, sex == 0, age > 20) %>% 
  group_by(age) %>% 
  mutate(rate = mean(haz)) %>% 
  distinct(age, rate) ->
  death_rates_and_cutoffs






# Function that generates a single transition time i -> j with hazards
# l1, ..., ln.
generate_single_from_piecewise_exp <- function(hazards, cutoffs, start) {
  
  stopifnot(!is.null(hazards))
  
  ddd <- data.frame(
    "A" = c(cutoffs, Inf),
    "B" = c(start, cutoffs)
  )
  
  ddd <- data.frame(
    "diff" = ddd$A - ddd$B,
    "rexp" = rexp(length(hazards), rate = hazards),
    "start" = ddd$B
  )
  ddd$res <- ddd$rexp < ddd$diff
  `[`(ddd$start + ddd$rexp, which(ddd$res == TRUE)[1])
}

# Helper functions that subsets only relevant cutoffs given
# current time (eg. if i->j happened at t, then hazards l(s), s<t are ignored)
modify_cutoffs <- function(start, cutoffs) {
  cutoffs[cutoffs > start]
}

subset_rates <- function(start, cutoffs, rates) {
  c(rates[c(which(cutoffs >= start), length(rates))])
}

# Do the same thing many times
rpexp <- function(n, hazards, cutoffs, start) {
  replicate(n, generate_single_from_piecewise_exp(hazards, cutoffs, start))
}

# Extract screening observations from latent data
create_event_array <- function(
    data_list
) {
  screening_events <- data_list$events %>% select(-id)
  screening_times <- data_list$times %>% select(-id)
  
  data_list
  
  
  screening_events <- as.matrix(screening_events)
  screening_times <- as.matrix(screening_times)
  screening_events[is.na(screening_events)] <- -1
  screening_events[screening_events == 2] <- 1
  screening_times[is.na(screening_times)] <- -1
  
  data_array <- array(NA, dim = c(nrow(screening_events), ncol(screening_events), 2))
  data_array[,,2] <- screening_events
  data_array[,,1] <- screening_times
  
  data_array
}

## Function to generate observations y | state, sensitivities -- these are
## the observed data in the screening process.
random_from_matrix <- function(state, a_sens = S_a, aa_sens = S_aa, c_sens = S_c) {
  confusion_list <- list(
    "Healthy" = "Healthy",
    "Adenoma" = sample(c("Healthy", "Adenoma"), 1, prob = c(1-a_sens, a_sens)),
    "AA" = sample(c("Healthy", "AA"), 1, prob = c(1-aa_sens, aa_sens)),
    "Pre-clin" = sample(c("Healthy", "Pre-clin"), 1, prob = c(1-c_sens, c_sens)),
    "Out of follow-up" = "Out of follow-up"
  )
  
  confusion_list[[state]]
}

## Function for generating data using smaller functions from piecewise constant/constant
## rates MSM with absorbing state observation exactly observed.
generate_msm <- function(N, mu1, mu2, mu3, mu4, cutoffs, entry_time, exit_time, death_rates_and_cutoffs) {
  stopifnot(
    length(mu4) == 1,
    !is.null(mu1),
    !is.null(mu4),
    !is.null(mu2),
    !is.null(mu3)
  )
  tibble(entry = entry_time, first = rpexp(N, mu1, cutoffs, entry_time)) %>% 
    rowwise() %>% 
    mutate(
      second = rpexp(1, subset_rates(first, cutoffs, mu2), modify_cutoffs(first, cutoffs), first),
      third = rpexp(1, subset_rates(second, cutoffs, mu3), modify_cutoffs(second, cutoffs), second),
      fourth = third + rexp(1, mu4),
      death = min(rpexp(1, death_rates_and_cutoffs$rate, death_rates_and_cutoffs$age[-length(death_rates_and_cutoffs$age)], entry_time), exit_time)
    ) %>% 
    mutate(
      dtime = death,
      death = case_when(
        is.na(fourth) ~ death,
        death < fourth ~ death,
        TRUE ~ 0
      ),
      clinical = case_when(
        fourth < death ~ fourth,
        death == 0 ~ fourth,
        TRUE ~ 0)
    ) 
}

## Prior generation function for piecewise-constant rates model
generate_ridge_params <- function(
    mu1_prior, mu2_prior, mu3_prior, mu4_prior,
    mu_ref_prior_1, mu_ref_prior_2, mu_ref_prior_3,
    mu_raw_prior_1, mu_raw_prior_2, mu_raw_prior_3,
    entry, exit, cutoff_places, entry_int
) {
  
  params <- list()
  
  params$mu1_mean <- mu1_prior[1]/mu1_prior[2]#rgamma(1, shape = mu1_prior[1], rate = mu1_prior[2])
  params$mu2_mean <- mu2_prior[1]/mu2_prior[2]#rgamma(1, shape = mu2_prior[1], rate = mu2_prior[2])
  params$mu3_mean <- mu3_prior[1]/mu3_prior[2]#rgamma(1, shape = mu3_prior[1], rate = mu3_prior[2])
  params$mu4 <- rgamma(1, shape = mu4_prior[1], rate = mu4_prior[2])
  
  params$mu_raw_1 <- rnorm(length(cutoff_places)+1, mu_raw_prior_1[1], mu_raw_prior_1[2])
  params$mu_raw_2 <- rnorm(length(cutoff_places)+1, mu_raw_prior_2[1], mu_raw_prior_2[2])
  params$mu_raw_3 <- rnorm(length(cutoff_places)+1, mu_raw_prior_3[1], mu_raw_prior_3[2])
  
  params$mu1 = exp(params$mu_raw_1) * params$mu1_mean
  params$mu2 = exp(params$mu_raw_2) * params$mu2_mean
  params$mu3 = exp(params$mu_raw_3) * params$mu3_mean
  
  params$mu_raw_ref_1 <- rnorm(1, mu_ref_prior_1[1], mu_ref_prior_1[2])
  params$mu_raw_ref_2 <- rnorm(1, mu_ref_prior_2[1], mu_ref_prior_2[2])
  params$mu_raw_ref_3 <- rnorm(1, mu_ref_prior_3[1], mu_ref_prior_3[2])
  
  params$mu1_ref = exp(params$mu_raw_ref_1) * params$mu1
  params$mu2_ref = exp(params$mu_raw_ref_2) * params$mu2
  params$mu3_ref = exp(params$mu_raw_ref_3) * params$mu3
  
  params
  
} 

## Prior generating function for constant rates model
generate_const_params <- function(
    mu1_prior, mu2_prior, mu3_prior, mu4_prior, mu_ref_prior_1,mu_ref_prior_2, mu_ref_prior_3,
    entry, exit, cutoff_places, entry_int
) {
  params <- list()
  
  params$mu1 <- rgamma(1, shape = mu1_prior[1], rate = mu1_prior[2])
  params$mu2 <- rgamma(1, shape = mu2_prior[1], rate = mu2_prior[2])
  params$mu3 <- rgamma(1, shape = mu3_prior[1], rate = mu3_prior[2])
  params$mu4 <- rgamma(1, shape = mu4_prior[1], rate = mu4_prior[2])
  
  
  params$mu_raw_ref_1 <- rnorm(1, mu_ref_prior_1[1], mu_ref_prior_1[2])
  params$mu_raw_ref_2 <- rnorm(1, mu_ref_prior_2[1], mu_ref_prior_2[2])
  params$mu_raw_ref_3 <- rnorm(1, mu_ref_prior_3[1], mu_ref_prior_3[2])
  
  
  params$mu1_ref <- exp(params$mu_raw_ref_1) * params$mu1
  params$mu2_ref <- exp(params$mu_raw_ref_2) * params$mu2
  params$mu3_ref <- exp(params$mu_raw_ref_3) * params$mu3
  
  params
  
}

## Event history data generation for those who have
## also screening observations
times_to_events_attenders <- function(dt_att) {
  dt_att %>% 
    mutate(
      status_at_60 = case_when(
        !is.na(third) & third < 60 ~ "Pre-clin",
        !is.na(second) & second < 60 ~ "AA",
        !is.na(first) & first < 60 ~ "Adenoma",
        (death > 0 & death < 60) | (!is.na(fourth) & fourth < 60) ~ "Out of follow-up",
        TRUE ~ "Healthy"
      ),
      status_at_62 = case_when(
        !is.na(third) & third < 62 ~ "Pre-clin",
        !is.na(second) & second < 62 ~ "AA",
        !is.na(first) & first < 62 ~ "Adenoma",
        (death > 0 & death < 62) | (!is.na(fourth) & fourth < 62) ~ "Out of follow-up",
        TRUE ~ "Healthy"
      ),
      status_at_64 = case_when(
        !is.na(third) & third < 64 ~ "Pre-clin",
        !is.na(second) & second < 64 ~ "AA",
        !is.na(first) & first < 64 ~ "Adenoma",
        (death > 0 & death < 64) | (!is.na(fourth) & fourth < 64) ~ "Out of follow-up",
        TRUE ~ "Healthy"
      ),
      status_at_66 = case_when(
        !is.na(third) & third < 66 ~ "Pre-clin",
        !is.na(second) & second < 66 ~ "AA",
        !is.na(first) & first < 66 ~ "Adenoma",
        (death > 0 & death < 66 )| (!is.na(fourth) & fourth < 66) ~ "Out of follow-up",
        TRUE ~ "Healthy"
      ),
      status_at_68 = case_when(
        !is.na(third) & third < 68 ~ "Pre-clin",
        !is.na(second) & second < 68 ~ "AA",
        !is.na(first) & first < 68 ~ "Adenoma",
        (death > 0 & death < 68) | (!is.na(fourth) & fourth < 68) ~ "Out of follow-up",
        TRUE ~ "Healthy"
      ) 
    )
}


## Event history data generation for 
## population without screening observations
create_event_time_data_ref_pop <- function(dt_pop, id_start) {
  
  
  dt_pop %>% 
    ungroup() %>% 
    mutate(id = row_number() + id_start) %>% 
    select(id, entry, death, clinical) %>% 
    pivot_longer(cols = entry:clinical, names_to = "event", values_to = "time") %>% 
    filter(time != 0) ->
    exact_events
  
  exact_events %>% 
    group_by(id) %>%
    arrange(time) %>%
    filter(time >= 30) %>%
    mutate(rank = row_number()) %>%
    mutate(time = round(time, digits = 1),
           time_prev = lag(time, 1L, order_by = time)) %>%
    group_by(id, time) %>%
    mutate(id_count = n()) %>%
    mutate(
      time_count = n(),
      has_clinical = "clinical" %in% event
    ) %>%
    ungroup(time) %>%
    filter((time_count == 1 | event != 0) & (has_clinical == FALSE | event == "clinical")) %>%
    filter((id_count == 1 | event != "death")) %>% 
    mutate(time = round(time, digits = 1),
           time_prev = lag(time, 1L, order_by = time)) %>% 
    group_by(id, time) %>% 
    mutate(id_count = n()) %>% 
    mutate(
      time_count = n(),
      has_clinical = "clinical" %in% event
    ) %>% 
    ungroup(time) %>% 
    filter((time_count == 1 | event != 0) & (has_clinical == FALSE | event == "clinical")) %>% 
    filter((id_count == 1 | event != "death"))  ->
    events_and_times_ref_pop
  
  
  events_and_times_ref_pop %>% 
    select(id, time, rank) %>%
    ungroup() %>% 
    pivot_wider(id_cols = id, names_from = rank, values_from = time) ->
    times_ref_pop
  
  events_and_times_ref_pop %>% select(id, event, rank) %>% 
    ungroup() %>% 
    mutate(event = case_when(
      event %in% c("entry", "death", "0") ~ 0,
      event == "Healthy" ~ 1,
      event == "Adenoma" ~ 3,
      event == "AA" ~ 4,
      event == "Pre-clin" ~ 5,
      event == "clinical" ~ 6
    )) %>%
    pivot_wider(id_cols = id, names_from = rank, values_from = event) ->
    events_ref_pop
  
  list("times" = times_ref_pop,
       "events" = events_ref_pop)
}

## Small function that defines some adenomas to be 
## unclassifiable between advanced and non-advanced
hide_adenomas <- function(
    dd_set, 
    perc
) {
  dd_set[,,2][dd_set[,,2] %in% c(3,4) & rbinom(length(dd_set[,,2] %in% c(3,4)), 1, perc)] <- 7L
  dd_set
}

## Combiner function that samples priors
## and generates data for all the subgroups
generate_data_for_sbc <- function(
    N,
    priors,
    parametrization,
    S_a_prior, S_aa_prior, S_c_prior,
    cutoffs,
    death_rates_and_cutoffs, entry_time = 30, exit_time = 80,
    entry_int = entry_time,
    hide_perc = 0.3
) {
  
  N_att <- round((1-0.33)*N, 0)
  N_ref <- N-N_att
  priors$cutoff_places <- cutoffs
  
  
  if(parametrization == "piecewise_ridge") {
    params <- do.call(generate_ridge_params, priors)
  }
  
  if(parametrization == "constant") {
    params <- do.call(generate_const_params, priors)
  }
  
  
  params$S_c <- ifelse(length(S_c_prior) < 2, S_c_prior, rbeta(1, S_c_prior[1], S_c_prior[2]))
  params$S_a <- ifelse(length(S_a_prior) < 2, S_a_prior, rbeta(1, S_a_prior[1], S_a_prior[2]))
  params$S_aa <- ifelse(length(S_aa_prior) < 2, S_aa_prior, rbeta(1, S_aa_prior[1], S_aa_prior[2]))
  
  dt_att <- generate_msm(N_att, params$mu1, params$mu2, params$mu3, params$mu4, cutoffs, entry_time, exit_time, death_rates_and_cutoffs)
  dt_ref <- generate_msm(N_ref, params$mu1_ref, params$mu2_ref, params$mu3_ref, params$mu4, cutoffs, entry_time, exit_time, death_rates_and_cutoffs)
  
  
  dt_att <- dt_att %>% ungroup() %>% mutate(id = row_number())
  dt_att <- times_to_events_attenders(dt_att)
  
  
  dt_att %>%
    rowwise() %>% 
    mutate(
      screen_60 = random_from_matrix(status_at_60, params$S_a, params$S_aa, params$S_c),
      screen_62 = random_from_matrix(status_at_62, params$S_a, params$S_aa, params$S_c),
      screen_64 = random_from_matrix(status_at_64, params$S_a, params$S_aa, params$S_c),
      screen_66 = random_from_matrix(status_at_66, params$S_a, params$S_aa, params$S_c),
      screen_68 = random_from_matrix(status_at_68, params$S_a, params$S_aa, params$S_c)
    ) ->
    dt_att
  
  dt_att %>%
    pivot_longer(cols = screen_60:screen_68, names_to = "tp", values_to = "event") %>% 
    filter(!event %in% c("Healthy", "Healthy_C")) %>% 
    mutate(
      time = case_when(
        tp == "screen_60" ~ 60,
        tp == "screen_62" ~ 62,
        tp == "screen_64" ~ 64,
        tp == "screen_66" ~ 66,
        tp == "screen_68" ~ 68
      )
    ) %>% 
    group_by(id) %>% 
    mutate(min_obs = min(time)) %>% 
    select(min_obs, id) %>% 
    distinct() ->
    min_obs
  
  
  
  dt_att %>% 
    select(id, entry, death, clinical) %>% 
    pivot_longer(cols = entry:clinical, names_to = "event", values_to = "time") %>% 
    filter(time != 0) ->
    exact_events
  
  
  dt_att %>% 
    select(id, entry, screen_60:screen_68, clinical, death) %>% 
    pivot_longer(cols = screen_60:screen_68, names_to = "tp", values_to = "event") %>% 
    mutate(
      time = case_when(
        tp == "screen_60" ~ 60,
        tp == "screen_62" ~ 62,
        tp == "screen_64" ~ 64,
        tp == "screen_66" ~ 66,
        tp == "screen_68" ~ 68
      )
    ) %>% 
    filter(event != "Out of follow-up") %>% 
    filter(time < death | death == 0, time < clinical | clinical == 0) %>% 
    select(id, event, time) %>% 
    bind_rows(exact_events) %>% 
    group_by(id) %>% 
    arrange(time) %>% 
    filter(time >= 30) %>% 
    left_join(min_obs, by = "id") %>% 
    filter(time <= min_obs | is.na(min_obs)) %>% 
    mutate(rank = row_number()) ->
    events_and_times
  
  
  
  
  
  events_and_times %>% 
    mutate(time = round(time, digits = 1),
           time_prev = lag(time, 1L, order_by = time)) %>% 
    group_by(id, time) %>% 
    mutate(id_count = n()) %>% 
    mutate(
      time_count = n(),
      has_clinical = "clinical" %in% event
    ) %>% 
    filter((time_count == 1 | event != 0), (has_clinical == FALSE | event == "clinical")) %>% 
    filter((id_count == 1 | event != "death")) ->
    events_and_times
  
  
  
  events_and_times %>% 
    select(id, time, rank) %>%
    ungroup() %>% 
    pivot_wider(id_cols = id, names_from = rank, values_from = time) ->
    times
  
  
  
  events_and_times %>% 
    select(id, time, event, rank) %>% 
    ungroup() %>% 
    mutate(event = case_when(
      event %in% c("entry", "death", "0") ~ 0,
      event == "Healthy" ~ 1,
      event == "Adenoma" ~ 3,
      event == "AA" ~ 4,
      event == "Pre-clin" ~ 5,
      event == "clinical" ~ 6
    )) %>%
    pivot_wider(id_cols = id, names_from = rank, values_from = event) ->
    events
  
  
  ref <- create_event_time_data_ref_pop(dt_pop = dt_ref, id_start = N_att)
  
  events <- replace_na(events, list(NA, NA, rep(-1, ncol(events)-1)))
  events_ref <- replace_na(ref$events, list(NA, NA, rep(-1, ncol(ref$events)-1)))
  
  times <- replace_na(times, list(NA, NA, rep(-1, ncol(times)-1)))
  times_ref <- replace_na(ref$times, list(NA, NA, rep(-1, ncol(ref$times)-1)))
  
  if(!is.null(cutoffs)) {
    observed_data_screening <- add_cutoffs(list(events = events, times = times), cutoffs, id_cols = "id")
    observed_data_refusers <- add_cutoffs(list(events = events_ref, times = times_ref), cutoffs, id_cols = "id")
  } else {
    observed_data_screening <- list(events = events, times = times)
    observed_data_refusers <- list(events = events_ref, times = times_ref)
  }
  
  
  observed_data_screening <- create_event_array(observed_data_screening)
  observed_data_refusers <- create_event_array(observed_data_refusers)
  
  
  ref_proportions <- ref_props_from_data(observed_data_screening, observed_data_refusers, cutoffs, entry_time)$ref_prop
  
  params$mu1_total_pop = params$mu1_ref * ref_proportions + params$mu1 * (1-ref_proportions)
  params$mu2_total_pop = params$mu2_ref * ref_proportions + params$mu2 * (1-ref_proportions)
  params$mu3_total_pop = params$mu3_ref * ref_proportions + params$mu3 * (1-ref_proportions)
  
  dt_pop <- generate_msm(
    N, 
    params$mu1_total_pop, params$mu2_total_pop, params$mu3_total_pop, params$mu4, 
    cutoffs, entry_time, exit_time, death_rates_and_cutoffs
  )
  
  pop <- create_event_time_data_ref_pop(dt_pop = dt_pop, id_start = N_att + N_ref)
  
  events_pop <- replace_na(pop$events, list(NA, NA, rep(-1, ncol(pop$events)-1)))
  times_pop <- replace_na(pop$times, list(NA, NA, rep(-1, ncol(pop$times)-1)))
  
  if(!is.null(cutoffs)) {
    observed_data_control <- add_cutoffs(list(events = events_pop, times = times_pop), cutoffs, id_cols = "id")
  } else {
    observed_data_control <- list(events = events_pop, times = times_pop)
  }
  observed_data_control <- create_event_array(observed_data_control)
  
  if(hide_perc > 0) {
    observed_data_screening <- hide_adenomas(observed_data_screening, hide_perc)
  }
  
  ret_list <- c(list("params" = params), 
                list(
                  "S_a_prior" = S_a_prior,
                  "S_aa_prior" = S_aa_prior,
                  "S_c_prior" = S_c_prior
                ),
                priors,
                list(
                  N_screen = nrow(observed_data_screening[,,2]), 
                  K_screen = ncol(observed_data_screening[,,2]), 
                  observed_data_screening = observed_data_screening,
                  K_refusers = ncol(observed_data_refusers[,,2]),
                  N_refusers = nrow(observed_data_refusers[,,2]),
                  observed_data_refusers = observed_data_refusers,
                  K_control = ncol(observed_data_control[,,2]),
                  N_control = nrow(observed_data_control[,,2]),
                  observed_data_control = observed_data_control,
                  ref_proportions = ref_proportions,
                  screen_counts = c(43765, 52890, 49386, 38215, 25627)
                )
  )
  
  if(!is.null(cutoffs)) {
    ret_list$cutoff_places <- cutoffs
    ret_list$N_cutoffs <- length(cutoffs)
  }
  
  ret_list
  
}

run_one_sim_hmc <- function(
    L,
    priors,
    bayesian_model,
    thin = 5,
    warmup = 2000,
    threads = 1
) {
  
  dataset <- do.call(generate_data_for_sbc, priors)
  
  
  mm <- bayesian_model$sample(
    data = dataset[-1], 
    output_dir = file.path(here::here(), "model_files"), 
    chains = 1,
    threads = threads,
    save_warmup = FALSE,
    iter_warmup = warmup,
    adapt_delta = 0.95,
    iter_sampling = L,
    thin = thin,
  )
  
  list(
    dataset,
    mm
  )
}


## Function that generates data and runs ADVI once.
run_one_sim_advi <- function(
    L,
    sbc_init_data,
    bayesian_model,
    algorithm
) {
  
  dataset <- do.call(generate_data_for_sbc, sbc_init_data)
  
  mm <- bayesian_model$variational(
    data = dataset[-1], 
    output_dir = file.path(here::here(), "model_files"), 
    algorithm = algorithm,
    threads = 1,
    output_samples = L,
    tol_rel_obj = 0.001
  )
  
  list(
    dataset,
    mm
  )
}

extract_results_tv_rates <- function(x, n) {
  x[[2]]$draws() %>% posterior::as_draws_matrix() %>% as_tibble() %>% select(-c("lp__")) %>% mutate_all(as.numeric) -> x2
  colnames(x2) <- str_remove(colnames(x2), "\\[") %>% str_remove(., "\\]")
  x2 %>% pivot_longer(cols = names(.), names_to = "par", values_to = "post_val")  -> x2
  x1 <- x[[1]]$params %>% unlist() %>% bind_rows() %>% 
    pivot_longer(cols = everything(), names_to = c("par"), values_to = "prior_vals") %>% 
    sample_n(size = n)
  left_join(x1, x2, by = "par")
}

extract_results_const_rates <- function(x, n) {
  x[[2]]$draws() %>% posterior::as_draws_matrix() %>% as_tibble() %>% select(-c("lp__")) %>% mutate_all(as.numeric) -> x2
  colnames(x2) <- str_remove(colnames(x2), "\\[\\d\\]")
  x2 %>% pivot_longer(cols = names(.), names_to = "par", values_to = "post_val") %>% 
    group_by(par) %>% 
    sample_n(size = n) %>% ungroup() -> 
    x2
  x1 <- x[[1]]$params %>% bind_rows() %>% 
    pivot_longer(cols = everything(), names_to = c("par"), values_to = "prior_vals")
  left_join(x1, x2, by = "par")
}

plot_sbc_bars <- function(results_object, rates = "const", remove_pars = "(S_c)|(S_a)|(S_aa)|(raw)", 
                          n = 49, breaks = seq(0, 50, by = 5), latex_pars = read.csv2("data/pars_to_latex.csv")) {
  if(rates == "const") {
    sbc_posts <- lapply(results_object, function(x) try(extract_results_const_rates(x, n)))
    ncols <- 7
  } else {
    sbc_posts <- lapply(results_object, function(x)  try(extract_results_tv_rates(x, n)))
    ncols <- 7
  } 
  
  sbc_posts <- lapply(sbc_posts, function(x) {
    if(is_tibble(x)) return(x)
    if(!is_tibble(x)) return(tibble(par = NA, prior_vals = NA, post_val = NA))
  })
  
  sbc_posts <- rbindlist(sbc_posts)
  
  sbc_posts %>% 
    left_join(latex_pars, by = "par") %>% 
    filter(!is.na(post_val), !str_detect(par, "S|total|raw")) %>% 
    group_by(par, latex_par, prior_vals) %>% summarise(rank_stat = sum(post_val > prior_vals)) %>% 
    ggplot() + geom_histogram(aes(y = rank_stat), breaks = seq(0, 50, by = 5)) + 
    geom_vline(aes(xintercept = 50)) + 
    geom_vline(aes(xintercept = qbinom(0.025, 500, 1/10)), linetype = "dashed") + 
    geom_vline(aes(xintercept = qbinom(0.975, 500, 1/10)), linetype = "dashed") + 
    facet_wrap(~par, ncol = ncols, labeller = jotain <- function(x) {
      latex_pars <- read.csv2("data/pars_to_latex.csv", allowEscapes = TRUE)
      latex_pars %>% right_join(x, by = "par") %>% 
        rowwise() %>% 
        mutate(latex_par = list(TeX(latex_par_2))) %>% 
        select(latex_par)
    }) + coord_flip() + theme_bw() + theme(strip.text.x = element_text(size = 15)) + xlab("Count") + ylab("Rank statistic") -> pp
  
  sbc_posts %>% 
    filter(!is.na(post_val), !str_detect(par, remove_pars)) %>%
    group_by(par) %>% 
    mutate(post_mean = mean(post_val),
           prior_mean = mean(prior_vals)) %>% 
    ggplot() + 
    geom_density(aes(x = post_val), linetype = "dashed") +
    geom_density(aes(x = prior_vals), linetype = "solid") +
    geom_vline(aes(xintercept = prior_mean), col = "red") +
    geom_vline(aes(xintercept = post_mean), col = "green") + 
    facet_wrap(~ par, labeller = jotain <- function(x) {
      latex_pars <- read.csv2("data/pars_to_latex.csv", allowEscapes = TRUE)
      latex_pars %>% right_join(x, by = "par") %>% 
        rowwise() %>% 
        mutate(latex_par = list(TeX(latex_par_2))) %>% 
        select(latex_par)
    }, scales = "free") + ylab("f(x)") + xlab("x") ->
    pp2
  
  list(
    "data" = sbc_posts,
    "plot1" = pp,
    "plot2" = pp2
  )
}







