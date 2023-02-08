ref_props_from_data <- function(array_a, array_r, cutoffs, entry) {
  pyrs_r <- extract_exit_times(array_r) %>% extract_pyrs(.,cutoffs, entry)
  pyrs_a <- extract_exit_times(array_a) %>% extract_pyrs(.,cutoffs, entry)
  
  left_join(pyrs_r, pyrs_a, by = "age", suffix = c(".ref", ".att")) %>% 
    mutate(ref_prop = pyrs_acc.ref/(pyrs_acc.ref + pyrs_acc.att))
  
}


extract_exit_times <- function(array_gen) {
  array_gen[,,1] %>% 
    as_tibble() %>% 
    mutate(exit_time = do.call(pmax, c(select(., 1:ncol(.)), na.rm = TRUE))) %>% 
    select(exit_time)
  
}

extract_pyrs <- function(exit_times, cutoffs, entry) {
  exit_times %>% mutate(entry_time = entry) %>% 
    Epi::Lexis(
      data = .,
      entry = list("age" = entry_time),
      exit = list("age" = exit_time),
      entry.status = 0,
      exit.status = 0
    ) -> ref_prop_lex
  
  if(!is.null(cutoffs)) {
    Epi::splitLexis(
      ref_prop_lex,
      breaks = cutoffs
    ) %>% 
      group_by(age) %>% 
      summarise(pyrs_acc = sum(lex.dur)) %>% 
      select(age, pyrs_acc) -> ret_lex
  } else {
    
    ret_lex <- ref_prop_lex %>% 
      summarise(pyrs_acc = sum(lex.dur)) %>% 
      mutate(age = entry) %>% 
      select(age, pyrs_acc)
  }
  
  ret_lex
}


