functions {
  /*
  * Alternative to neg_binomial_2_log_rng() that 
  * avoids potential numerical problems during warmup
  */
  int neg_binomial_2_log_safe_rng(real eta, real phi) {
    real gamma_rate = gamma_rng(phi, phi / exp(eta));
    if (gamma_rate >= exp(20.79))
      return -9;
      
    return poisson_rng(gamma_rate);
  }
}
data {
  int<lower=1> N;                     
  int<lower=0> complaints[N];              
  vector<lower=0>[N] traps;                
  
  // 'exposure'
  vector[N] log_sq_foot;  
  
  // building-level data
  int<lower=1> K;
  int<lower=1> J;
  int<lower=1, upper=J> building_idx[N];
  matrix[J,K] building_data;
}
parameters {
  real<lower=0> inv_phi;       // 1/phi (easier to think about prior for 1/phi instead of phi)
  real beta;                   // coefficient on traps
  
  vector[J] alpha;            // buildings-specific intercepts
  real<lower=0> sigma_alpha;  // sd of building-specific intercepts
  real mu;                    // intercept of model for alpha
  vector[K] zeta;             // coefficients on building-level predictors in model for mu 
}
transformed parameters {
  real phi = inv(inv_phi);
}
model {
  alpha ~ normal(mu + building_data * zeta, sigma_alpha);
  sigma_alpha ~ normal(0, 1);
  mu ~ normal(log(4), 1);
  zeta ~ normal(0, 1);  // could also use informative priors on the different elements
  beta ~ normal(-0.25, 1);
  inv_phi ~ normal(0, 1);
  
  complaints ~ neg_binomial_2_log(alpha[building_idx] + beta * traps + log_sq_foot, phi);
} 
generated quantities {
  int y_rep[N];
  for (n in 1:N) {
    real eta_n = alpha[building_idx[n]] + beta * traps[n] + log_sq_foot[n];
    y_rep[n] = neg_binomial_2_log_safe_rng(eta_n, phi);
  }
}
