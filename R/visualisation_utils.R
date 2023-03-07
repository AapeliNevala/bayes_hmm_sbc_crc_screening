## Functions used for creating figures in "Bayesian hidden Markov model for natural history of colorectal cancer: 
## handling misclassified observations, varying observation schemes and unobserved data"

posterior_state_probs_from_draws <- function(model_fit, time_up, time_lo, KM = NULL) {
  model_fit$draws() %>% posterior::summarise_draws() %>% 
    filter(stringr::str_detect(variable, "state_probs")) %>% 
    mutate(process = str_split_fixed(variable, "(\\[)|(\\,)", 3)[,1],
           state = str_split_fixed(variable, "(\\[)|(\\,)", 3)[,2],
           time = as.numeric(str_replace(str_split_fixed(variable, "(\\[)|(\\,)", 3)[,3], "\\]", "")),
           group = str_split_fixed(variable, "(\\_)|(\\[)", 4)[,3]) %>% 
    mutate(state_name = case_when(
      state == 1 ~ "(1) Healthy",
      state == 2 ~ "(2) Non-advanced adenoma",
      state == 3 ~ "(3) Advanced adenoma",
      state == 4 ~ "(4) Pre-clinical CRC",
      state == 5 ~ "(5) Clinical CRC"
    ))  %>% 
    mutate(
      mean_exp = mean,
      q5_exp = q5,
      q95_exp = q95
    ) %>% 
    mutate(time = time + 30) ->
    five_state_model_probs_me
  
  
  
  five_state_model_probs_me %>% 
    filter(time < time_up, time >= time_lo, state %in% c(3,4)) %>% 
    mutate(line = case_when(
      group != "attenders" & state != 5 ~ 0.05,
      state != 5 & time < 60 ~ pmax(exp(-(time-60)^2/5), 0.05),
      state != 5 & time > 68 ~ pmax(exp(-(time-68)^2/5), 0.05),
      state == 5 ~ 1
    )) %>% 
    ggplot(aes(x = time)) +  
    geom_line(aes(y = q5_exp, linetype = state_name), size = 0.5) +
    geom_line(aes(y = q95_exp, linetype = state_name), size = 0.5) +
    geom_line(aes(y = mean_exp, linetype = state_name), size = 1.2) +
    geom_vline(aes(xintercept = c(60))) +
    geom_vline(aes(xintercept = c(68))) +
    geom_label(data = five_state_model_probs_me %>% filter(time %in% c(60, 68, 74), state %in% c(3,4)), 
               aes(x = time, y = mean_exp, label = round(mean_exp, 3)*100), 
               size = 3, position = position_dodge(5), 
               color = "black") + 
    facet_wrap(~ group, scales = "fixed", labeller = label_parsed) + xlab("")  + ylab("") + # scale_color_manual(values = c("#56B4E9", "#F0E442")) + 
    labs(linetype = "State") + guides(alpha = "none") +
    theme_minimal() + theme(legend.position="top") + xlim(c(time_lo, time_up+1)) ->
    adenoma_plot
  
  five_state_model_probs_me %>% 
    filter(time < time_up, time >= time_lo, state %in% c(5)) %>% 
    mutate(line = case_when(
      group != "attenders" & state != 5 ~ 0.05,
      state != 5 & time < 60 ~ pmax(exp(-(time-60)^2/5), 0.05),
      state != 5 & time > 68 ~ pmax(exp(-(time-68)^2/5), 0.05),
      state == 5 ~ 1
    )) -> c_plot_data
  
  if(is.null(KM)) {
    c_plot_data %>%   
      ggplot(aes(x = time)) + 
      geom_line(aes(y = q5_exp, linetype = state_name), size = 0.5) +
      geom_line(aes(y = q95_exp, linetype = state_name), size = 0.5) +
      geom_line(aes(y = mean_exp, linetype = state_name), size = 1.2) +
      geom_vline(aes(xintercept = c(60))) +
      geom_vline(aes(xintercept = c(68))) + 
      geom_label(data = five_state_model_probs_me %>% filter(time %in% c(60, 68, 74), state %in% c(5)), 
                 aes(x = time, y = mean_exp, label = round(mean_exp, 3)*100), 
                 size = 3, position = position_dodge(5), 
                 color = "black") +
      facet_wrap(~ group, scales = "fixed", labeller = label_parsed) + xlab("Age") + # scale_color_manual(values = c("#D55E00", "#000000")) + 
      labs(linetype = "State") + guides(alpha = "none") +
      theme_minimal() + theme(legend.position="top") + xlim(c(time_lo, time_up+1)) + ylab("") ->
      cancer_plot
  } else {
    c_plot_data %>%   
      ggplot(aes(x = time)) + 
      geom_line(aes(y = q5_exp, linetype = state_name), size = 0.5) +
      geom_line(aes(y = q95_exp, linetype = state_name), size = 0.5) +
      geom_line(aes(y = mean_exp, linetype = state_name), size = 1.2) +
      geom_vline(aes(xintercept = c(60))) +
      geom_vline(aes(xintercept = c(68))) + 
      geom_step(data = KM %>% filter(time >= time_lo, time <= time_up), aes(x = time, y = 1-surv), color = "Red") +
      geom_label(data = five_state_model_probs_me %>% filter(time %in% c(60, 68, 74), state %in% c(5)), 
                 aes(x = time, y = mean_exp, label = round(mean_exp, 3)*100), 
                 size = 3, position = position_dodge(5), 
                 color = "black") +
      facet_wrap(~ group, scales = "fixed", labeller = label_parsed) + xlab("Age") + # scale_color_manual(values = c("#D55E00", "#000000")) + 
      labs(linetype = "State") + guides(alpha = "none") +
      theme_minimal() + theme(legend.position="top") + xlim(c(time_lo, time_up+1)) + ylab("") ->
      cancer_plot
  }
  
  
  
  
  five_state_model_probs_me %>% 
    filter(time < time_up, time >= time_lo, state %in% c(1,2)) %>% 
    mutate(line = case_when(
      group != "attenders" & state != 5 ~ 0.05,
      state != 5 & time < 60 ~ pmax(exp(-(time-60)^2/5), 0.05),
      state != 5 & time > 68 ~ pmax(exp(-(time-68)^2/5), 0.05),
      state == 5 ~ 1
    )) %>% 
    ggplot(aes(x = time)) + 
    geom_line(aes(y = q5_exp, linetype = state_name), size = 0.5) +
    geom_line(aes(y = q95_exp, linetype = state_name), size = 0.5) + 
    geom_line(aes(y = mean_exp, linetype = state_name), size = 1.2) +
    geom_vline(aes(xintercept = c(60))) +
    geom_vline(aes(xintercept = c(68))) +
    geom_label(data = five_state_model_probs_me %>% filter(time %in% c(60, 68, 74), state %in% c(1,2)), 
               aes(x = time, y = mean_exp, label = round(mean_exp, 3)*100), 
               size = 3, position = position_dodge(5), 
               color = "black") + 
    facet_wrap(~ group, scales = "fixed", labeller = label_parsed) + # scale_color_manual(values = "#009E73") + 
    labs(linetype = "State") + guides(alpha = "none") + xlab("") +
    theme_minimal() + theme(legend.position="top") + xlim(c(time_lo, time_up+1)) + ylab("") ->
    healthy_plot
  
  ggpubr::ggarrange(healthy_plot, adenoma_plot, cancer_plot, ncol = 1, common.legend = FALSE) -> returned_plot
  ggpubr::annotate_figure(returned_plot, left = "State occupancy probability")
}

table_results <- function(bayesian_model, incl_params, excl_params, tex_pars = read.csv2("data/latex_pars.csv")) {
  bayesian_model$draws() %>% posterior::summarise_draws() %>% 
    filter(str_detect(variable, incl_params), !str_detect(variable, excl_params)) %>% 
    select(variable, mean, median, sd, q5, q95) %>% 
    mutate(variable = str_remove_all(variable, "(\\[)|(\\])")) %>% 
    left_join(tex_pars, by = c("variable" = "par")) %>% 
    select(-variable) %>% 
    xtable::xtable(digits = 4, row.names = FALSE)
}



