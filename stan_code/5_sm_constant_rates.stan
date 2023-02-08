functions {
    
  #include stan_utils/stan_utility_helpers.stan
  #include 5_state_model_matrix.stan

  matrix generated_state_probs(vector mu1, vector mu2, vector mu3, vector mu4, int entry_int) {
    matrix[5,99-entry_int] state_probs = rep_matrix(negative_infinity(), 5, 99-entry_int);
    matrix[4, 4] P;
    matrix[5,99-entry_int] exp_state_probs;
    vector[4] acc;
    real m1 = mu1[1];
    real m2 = mu2[1];
    real m3 = mu3[1];
    real m4 = mu4[1];

    state_probs[1,1] = 0;
    
    for(t in 2:(99-entry_int)) {

      
      P = transition_matrix(m1, m2, m3, m4, 1.0);
      
      for(i in 1:4) {
        acc = rep_vector(negative_infinity(), 4);
        for(s in 1:i) {
          acc[s] = P[s, i] + state_probs[s, t-1];
        }
        
        
        state_probs[i, t] = log_sum_exp(acc);
      }
      
      state_probs[5, t] = log1m_exp(log_sum_exp(state_probs[1:4, t]));
      
    }
    
    exp_state_probs = exp(state_probs);
    return(exp_state_probs);
  }
  
  
  int[] possible_states(
    int obs
  ) {
    if(obs == 3) return {2};
    if(obs == 4) return {3};
    if(obs == 5) return {4};
    if(obs == 6) return {4};
    if(obs == 2) return {1};
    if(obs == 0) return {1,2,3,4};
    if(obs == 7) return {2,3};
    if(obs == 8) return {1,2};
    return {1,2,3,4};
  }

  vector forward_pass(
    int obs,
    matrix P,
    matrix phii,
    vector gamma_prev
  ) {
    
    vector[4] gamma = rep_vector(negative_infinity(), 4);
    vector[4] acc = rep_vector(negative_infinity(), 4);
    
    for(i in possible_states(obs)) {
      for(j in 1:i) {
        acc[j] = gamma_prev[j] + P[j,i];
      }
      gamma[i] = log_sum_exp(acc);
      if(obs != 0 && obs != 7) {
        gamma[i] += phii[i, obs];
      }
    
    }
    
    return gamma;
  }
  


  real hmm_loop(
    real[] times,
    real[] obs_real,
    vector mu1,
    vector mu2,
    vector mu3,
    vector mu4, 
    matrix phii,
    real miss_mix
  ) {
    int dimension = num_elements(obs_real);
    int obs[dims(obs_real)[1]];
    int j = 1;
    real m1 = mu1[1];
    real m2 = mu2[1];
    real m3 = mu3[1];
    real m4 = mu4[1];

    
    matrix[4,4] P = rep_matrix(negative_infinity(), 4, 4);
    vector[4] gamma = rep_vector(negative_infinity(), 4);
    gamma[1] = 0.0;
    real retval = log_sum_exp(gamma);

    
    
    for(i in 1:dimension) {
      obs[i] = real_2_int(obs_real[i]);
      if(obs[i] < 0 || times[i] < 0) {
        dimension = i-1;
      }
    }
  
    for(m in 2:dimension) {
      
      P = transition_matrix(m1, m2, m3, m4, times[m] - times[m-1]);
      
      gamma = forward_pass(obs[m], P, phii, gamma);

      
      if(obs[m] == 6) {
        gamma += log(m4);
      }
      
      retval = log_sum_exp(gamma);

      if(obs[m] == 7) {
        retval += log_mix(exp(gamma[2] - log_sum_exp(gamma[2:3])), phii[2,7], phii[3,7]);
      }



    }
    
    if(retval == negative_infinity()) {
      reject("WAAAA", times[2]);
    }

    
    
    return retval;
  }
  

  
  real partial_sum(
    real[,,] observed_data, int start, int end,
    vector mu1, vector mu2, vector mu3, vector mu4, matrix phii, real miss_mix
  ) {
    real sum_var = 0;
    
    for(i in 1:dims(observed_data)[1]) {
      sum_var += hmm_loop(
        observed_data[i,,1], 
        observed_data[i,,2],
        mu1, mu2, mu3, mu4, 
        phii,
        miss_mix
      );
    }
    
    return sum_var;
  }
  
  
}

data {
  int<lower = 1> N_screen;
  int<lower = 1> N_refusers;         // num subjects
  int<lower = 1> N_control;   
  
  int<lower = 1> K_screen;
  int<lower = 1> K_refusers;         // num subjects
  int<lower = 1> K_control;
  
  real observed_data_screening[N_screen, K_screen, 2];
  real observed_data_control[N_control, K_control, 2];
  real observed_data_refusers[N_refusers, K_refusers, 2];

  real<lower = 0> mu1_prior[2];
  real<lower = 0> mu2_prior[2];
  real<lower = 0> mu3_prior[2];
  real<lower = 0> mu4_prior[2];
  
  real<lower = 0> mu_ref_prior_1[2];
  real<lower = 0> mu_ref_prior_2[2];
  real<lower = 0> mu_ref_prior_3[2];

  
  real<lower = 0, upper = 1> S_c_prior;
  real<lower = 0, upper = 1> S_a_prior;
  real<lower = 0, upper = 1> S_aa_prior;
  
  int entry_int;
  // int N_cutoffs;
  real ref_proportions;
}


parameters {
  real<lower = 0, upper = 1> miss_mix;
  vector<lower = 0>[1] mu1;                // Transition rate from healthy to progressive pre-invasive stage
  vector<lower = 0>[1] mu2;                // Transition rate for progressive screen-detectable to clinical stage
  vector<lower = 0>[1] mu3;
  vector<lower = 0>[1] mu4;


  real mu_raw_ref_1;
  real mu_raw_ref_2;
  real mu_raw_ref_3;


}

transformed parameters {
    
  real<lower = 0, upper = 1> S_c = S_c_prior;
  real<lower = 0, upper = S_c> S_aa = S_aa_prior;
  real<lower = 0, upper = S_aa> S_a = S_a_prior;

  
  vector[1] mu1_ref = exp(mu_raw_ref_1) * mu1;
  vector[1] mu2_ref = exp(mu_raw_ref_2) * mu2;
  vector[1] mu3_ref = exp(mu_raw_ref_3) * mu3;
  // vector<lower = 0>[N_cutoffs+1] mu4_ref = exp(mu4_raw_ref) * mu4;
  // real<lower = 0, upper = 1> Spec;


  vector[1] mu1_total_pop = mu1_ref * ref_proportions + mu1 * (1-ref_proportions);
  vector[1] mu2_total_pop = mu2_ref * ref_proportions + mu2 * (1-ref_proportions);
  vector[1] mu3_total_pop = mu3_ref * ref_proportions + mu3 * (1-ref_proportions);
  // vector<lower = 0>[N_cutoffs+1] mu4_total_pop = mu4_ref .* ref_proportions + mu4 .* (1-ref_proportions);

}


model {
  
  int grainsize = 1;
  
  matrix[5,7] phii = rep_matrix(negative_infinity(), 5, 7);
  
  mu1 ~ gamma(mu1_prior[1], mu1_prior[2]);
  mu2 ~ gamma(mu2_prior[1], mu2_prior[2]);
  mu3 ~ gamma(mu3_prior[1], mu3_prior[2]);
  mu4 ~ gamma(mu4_prior[1], mu4_prior[2]);

  mu_raw_ref_1 ~ normal(mu_ref_prior_1[1], mu_ref_prior_1[2]);
  mu_raw_ref_2 ~ normal(mu_ref_prior_2[1], mu_ref_prior_2[2]);
  mu_raw_ref_3 ~ normal(mu_ref_prior_3[1], mu_ref_prior_3[2]);


  phii[1,1] = 0.0;
  phii[2,3] = log(S_a);
  phii[2,1] = log1m(S_a);
  phii[3,4] = log(S_aa);
  phii[3,1] = log1m(S_aa); 
  phii[4,5] = log(S_c);
  phii[4,1] = log1m(S_c);
  phii[4,6] = 0.0;
  phii[2,7] = log(S_a);
  phii[3,7] = log(S_aa);

  target += reduce_sum(
    partial_sum, 
    observed_data_screening, 
    grainsize,
    mu1, mu2, mu3, mu4, 
    phii,
    miss_mix
  );
  
  target += reduce_sum(
    partial_sum,
    observed_data_refusers,
    grainsize,
    mu1_ref, mu2_ref, mu3_ref, mu4,
    phii,
    miss_mix
  );

  target += reduce_sum(
    partial_sum,
    observed_data_control,
    grainsize,
    mu1_total_pop,
    mu2_total_pop,
    mu3_total_pop,
    mu4,
    phii,
    miss_mix
  );

}

generated quantities {
  matrix[5,99-entry_int] state_probs_general = generated_state_probs(mu1_total_pop, mu2_total_pop, mu3_total_pop, mu4, entry_int);
  matrix[5,99-entry_int] state_probs_refusers = generated_state_probs(mu1_ref, mu2_ref, mu3_ref, mu4, entry_int);
  matrix[5,99-entry_int] state_probs_attenders = generated_state_probs(mu1, mu2, mu3, mu4, entry_int);

  // matrix[K_screen, 4] state_occ_person[N_screen];
  // matrix[N_screen, K_screen] obs_person;
  // 
  // 
  // for(i in 1:N_screen) {
  //   state_occ_person[i] = prediction_loop(
  //     observed_data[i,,1],
  //     observed_data[i,,2],
  //     mu1, mu2, mu3, mu4,
  //     phii
  //   )
  // }
  // 
  // for(i in 1:N_screen) {
  //   for(k in 1:K_screen) {
  //     if(observed_data[i,,2] == 1 || observed_data[i,,2] == 3 || observed_data[i,,2] == 4 || observed_data[i,,2] == 5 || observed_data[i,,2] == 7) {
  //       obs_person[i,k] = multinomial_rng(exp(state_occ_person[i][k,1:4])/sum(exp(state_occ_person[i][k,1:4])));
  //     } else {
  //       obs_person[i,k] = 0;
  //     }
  //   }
  // }

}
