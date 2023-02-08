matrix transition_matrix(
  real mu1, real mu2,
  real mu3, real mu4,
  real ts
  ) {
    matrix[4,4] P = rep_matrix(negative_infinity(), 4, 4);
    
    P[1,1] = -mu1 * ts;
    P[1,2] = log(mu1*(exp(-mu2*ts)- exp(-mu1 * ts)) / (mu1-mu2));
    P[1,3] = log(mu1*mu2/((mu2-mu1) * (mu3 - mu2) * (mu1 - mu3)) * 
    ((mu2-mu3) * exp(-mu1 * ts) + (mu3-mu1) * exp(-mu2 * ts) + (mu1-mu2) * exp(-mu3 * ts)));
    P[1,4] = log(mu1*mu2*mu3/((mu2-mu1) * (mu2-mu3) * (mu3-mu1)) * (
      (mu3-mu2)/(mu4-mu1) * (exp(-mu4 * ts)-exp(- mu1 * ts)) + 
      (mu1-mu3)/(mu4-mu2) * (exp(-mu4 * ts)-exp(- mu2 * ts)) +
      (mu2-mu1)/(mu4-mu3) * (exp(-mu4 * ts)-exp(- mu3 * ts))));
      
      P[2,2] = -mu2 * ts;
      P[2,3] = log(mu2 / (mu2-mu3) * (exp(-mu3*ts)- exp(-mu2 * ts)));
      P[2,4] = log(mu2*mu3/((mu3-mu2) * (mu4 - mu3) * (mu2 - mu4)) * 
      ((mu3-mu4) * exp(-mu2 * ts) + (mu4-mu2) * exp(-mu3 * ts) + (mu2-mu3) * exp(-mu4 * ts)));
      
      P[3,3] = -mu3 * ts;
      P[3,4] = log(mu3/(mu3-mu4)*(exp(-mu4*ts)- exp(-mu3 * ts)));
      
      P[4,4] = -mu4 * ts;
      
      return P;
}