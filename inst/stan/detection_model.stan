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
  int t_obs_max;          // upper bound on effective sample size (used to clamp t_obs_int)

}

// transformed data{
//   
//   int t_obs = N - K;
//   
// }

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

  for(s in 1:S_obs){
    for(i in 1:N_new){
      P_new[s, i] = inv_logit(beta[s, ] * X_new[, i]);
    }
  }

  // --- get information matrices --- //
  vector[N_new] t_obs = rep_vector(0.0, N_new);

  for(s in 1:S_obs){
    vector[N] p_s = inv_logit(beta[s, ] * X)';
    vector[N] weights_s = p_s .* (1 - p_s);
    matrix[K, K] I_beta_s = X * diag_matrix(weights_s) * X';
    matrix[K, K] I_beta_s_inv = inverse(I_beta_s);
    vector[N_new] t_obs_s;
    for(i in 1:N_new){
      // per-design-point gradient: dP/d_beta = p*(1-p)*x_i, shape [K]
      vector[K] grad_s_i = X_new[, i] * (P_new[s, i] * (1 - P_new[s, i]));
      t_obs_s[i] =   P_new[s, i] * (1 - P_new[s, i]) / (quad_form(I_beta_s_inv, grad_s_i));
    }
    t_obs += t_obs_s;
  }

  t_obs = t_obs / S_obs;

  // integer version needed for loop bounds; clamped to avoid unbounded W_new_i allocation
  array[N_new] int t_obs_int;
  for(i in 1:N_new) t_obs_int[i] = min(to_int(t_obs[i]), t_obs_max);

  // --- generate per-species detection counts (iNEXT incidence-freq format) ---
  // incfreq_new[1, i]     = t_obs_int[i]  (effective sample size)
  // incfreq_new[s + 1, i] = detections for species s out of t_obs_int[i] pseudo-sites
  array[S_obs + 1, N_new] int incfreq_new = rep_array(0, S_obs + 1, N_new);

  for(i in 1:N_new){
    incfreq_new[1, i] = t_obs_int[i];
    array[S_obs, t_obs_int[i]] int W_new_i = rep_array(0, S_obs, t_obs_int[i]);
    for(s in 1:S_obs){
      W_new_i[s, ] = bernoulli_rng(rep_array(P_new[s, i], t_obs_int[i]));
      incfreq_new[s + 1, i] = sum(W_new_i[s, ]);
    }
  }

}



