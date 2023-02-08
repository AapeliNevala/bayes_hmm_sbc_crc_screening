library(cmdstanr)
library(tidyverse)
library(data.table)

add_cutoffs <- function(
  times_events_datalist,
  cutoffs,
  id_cols = NULL
) {
  
  dd1 <- times_events_datalist$events
  dd2 <- times_events_datalist$times
  
  if(is.null(id_cols)) id_cols <- names(dd1)[1:3]
  
  dd2 %>% 
    rowwise() %>% 
    mutate(
      max_t = pmax(!!!rlang::syms(as.character(1:(ncol(dd2)-length(id_cols)))), na.rm = TRUE)
    ) ->
  dd2
  
  dd2 %>% relocate(max_t, .after = id_cols[length(id_cols)]) -> dd2
  
  
  dd1 %>% cbind(matrix(rep(0, nrow(times_events_datalist$events) * length(cutoffs)), byrow = TRUE, nrow = nrow(times_events_datalist$events))) ->
    dd1
  
  dd2 %>% cbind(matrix(rep(cutoffs, nrow(times_events_datalist$events)), byrow = TRUE, nrow = nrow(times_events_datalist$events))) ->
    dd2
  
  names(dd1) <- c(id_cols, 1:(ncol(dd1)-length(id_cols)))
  names(dd2) <- c(names(dd1)[1:length(id_cols)], "max_t", names(dd1)[(length(id_cols)+1):length(names(dd1))])
  
  dd1 %>% pivot_longer(cols = as.character(1:(ncol(dd1)-length(id_cols))), names_to = "obs", values_to = "val") %>% 
    filter(!is.na(val)) ->
    dd1
  
  dd2 %>% pivot_longer(cols = as.character(1:(ncol(dd2)-length(id_cols)-1)), names_to = "obs", values_to = "val") %>% 
    filter(!is.na(val)) ->
    dd2
  
  dd3 <- left_join(dd1, dd2, by = c(id_cols, "obs"))
  
  setDT(dd3)
  dd3[, obs := as.numeric(obs)]
  dd3 <- unique(dd3, by = names(dd3))
  dd3[ , obs := frank(val.y, ties.method = "average"), keyby = id_cols]
  dd3 <- dd3[val.y <= max_t]
  dd3 <- dd3[obs %% 1 == 0 | val.x != 0]
  dd3[ , obs := frank(val.y, ties.method = "average"), keyby = id_cols]
  dd3 <- dd3[order(obs)]
  
  list(
    times = dd3 %>% select(-val.x) %>% pivot_wider(id_cols = id_cols, names_from = obs, values_from = val.y),
    events = dd3 %>% select(-val.y) %>% pivot_wider(id_cols = id_cols, names_from = obs, values_from = val.x)
  )
}


