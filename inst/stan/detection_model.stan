// 
// This file defines a general glm-like model for species-level detection.
// 

data {
  
  int N;                    // number of plots/sites/quadrats
  int S_obs;                // number of observed species
  int K;                    // number of columns in design matrix
  array[S_obs, N] int W;    // raw incidence matrix
  matrix[K, N] X;           // covariate/design matrix
  
  // --- set design for iNEXT plot comparisons --- //
  int N_new;              // number of design points
  matrix[K, N_new] X_new; // design matrix
  
}

transformed data{
  
  int t_obs = N - K;
  
}

parameters {
  
  matrix[S_obs, K] beta; // regression coefs
  
}

model {
  
  // --- priors --- //
  for(s in 1:S_obs){
    beta[s, ] ~ normal(0, 1);
  }
  
  for(s in 1:S_obs){
    for(i in 1:N){
      target += bernoulli_logit_lpmf(W[s, i] | beta[s, ] * X[, i]);
    }
  }
  
}

generated quantities {
  
  // --- generate the detection probs for propagating uncertainties --- //
  matrix[S_obs, N_new] P_new;
  matrix[t_obs, N_new] Qt_new;
  
  for(s in 1:S_obs){
    for(i in 1:N_new){
      P_new[s, i] = inv_logit(beta[s, ] * X_new[, i]);
    }
  }
  
  // --- construct expected incidence frequency counts --- //
  for(i in 1:N_new){
    for(k in 1:t_obs){
      Qt_new[k, i] = sum(exp(binomial_lpmf(k | t_obs, P_new[, i])));
    }
  }
  
}



