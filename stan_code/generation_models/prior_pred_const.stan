functions {
    
    
  #include stan_utils/stan_utility_helpers.stan
  #include 5_state_model_matrix.stan

  matrix generated_state_probs(real mu1, real mu2, real mu3, real mu4, int entry_int) {
    matrix[5,99-entry_int] state_probs = rep_matrix(negative_infinity(), 5, 99-entry_int);
    matrix[4, 4] P;
    matrix[5,99-entry_int] exp_state_probs;
    vector[4] acc;

    state_probs[1,1] = 0;
    
    for(t in 2:(99-entry_int)) {

      
      P = transition_matrix(mu1, mu2, mu3, mu4, 1.0);
      
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


}

data {

  
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

}

generated quantities {
  


  real<lower = 0, upper = 1> S_c = 0.4;
  real<lower = 0, upper = S_c> S_aa = 0.1;
  real<lower = 0, upper = S_aa> S_a = 0.05;
  
  real<lower = 0> mu1 = gamma_rng(mu1_prior[1], mu1_prior[2]);
  real<lower = 0> mu2 = gamma_rng(mu2_prior[1], mu2_prior[2]);
  real<lower = 0> mu3 = gamma_rng(mu3_prior[1], mu3_prior[2]);
  real<lower = 0> mu4 = gamma_rng(mu4_prior[1], mu4_prior[2]);
  
  real mu_raw_ref_1 = normal_rng(mu_ref_prior_1[1], mu_ref_prior_1[2]);
  real mu_raw_ref_2 = normal_rng(mu_ref_prior_2[1], mu_ref_prior_2[2]);
  real mu_raw_ref_3 = normal_rng(mu_ref_prior_3[1], mu_ref_prior_3[2]);
  
  real mu1_ref = exp(mu_raw_ref_1) * mu1;
  real mu2_ref = exp(mu_raw_ref_2) * mu2;
  real mu3_ref = exp(mu_raw_ref_3) * mu3;


  matrix[5,99-entry_int] state_probs_refusers = generated_state_probs(mu1_ref, mu2_ref, mu3_ref, mu4, entry_int);
  matrix[5,99-entry_int] state_probs_attenders = generated_state_probs(mu1, mu2, mu3, mu4, entry_int);

}
