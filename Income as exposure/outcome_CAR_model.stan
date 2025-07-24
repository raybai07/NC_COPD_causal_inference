// Authors: Evan Funderburg, Emily Sitnik, Mark Ritchie, and Ray Bai
// suitably modified from https://mc-stan.org/learn-stan/case-studies/mbjoseph-CARStan.html

data {
  int<lower=1> n;
  int<lower=1> p;
  // response vector y, confounders X, treatment variable A, generalized propensity score vector gps
  int<lower=0> y[n];
  matrix[n, p] X;
  matrix<lower=0, upper=1>[n, n] W;
  vector[n] A;
  vector<lower=0>[n] gps;
  real a_new;
  vector[n] log_offset;
  
}
transformed data {
  vector[n] zeros;
  matrix<lower=0>[n, n] D;
  {
    vector[n] W_rowsums;
    for (i in 1 : n) {
      W_rowsums[i] = sum(W[i,  : ]);
    }
    D = diag_matrix(W_rowsums);
  }
  zeros = rep_vector(0, n);
  
  // Standardization to improve MCMC convergence
  real mean_A = mean(A);
  real sd_A = sd(A);
  vector[n] A_std = (A - mean_A) / sd_A;
  
  real mean_gps = mean(gps);
  real sd_gps = sd(gps);
  vector[n] gps_std = (gps - mean_gps) / sd_gps;
  
  real anew_std = (a_new - mean_A) / sd_A;

  // confounders X already standardized, don't need to also standardize beta
}  

parameters {
  real beta0_std;
  real betaA_std;
  vector[p] beta;      // confounders X already standardized
  real betar_std;
  vector[n] phi;
  real<lower=0> sigma_phi;
  real<lower=0, upper=1> rho_phi; // will automtically get standard uniform prior
}

model {
  // Samples from multivariate normal with zero mean vector and precision (inverse 
  // covariance) matrix tau*(D-alpha*W)
  phi ~ multi_normal_prec(zeros, (D - rho_phi * W)/square(sigma_phi)); 
  
  // All other priors
  beta0_std ~ normal(0, 5);
  betaA_std ~ normal(0, 5);
  beta ~ normal(0,5);
  betar_std ~ normal(0, 5);
  sigma_phi ~ cauchy(0, 2.5);
  
  // Likelihood
   y ~ poisson_log(beta0_std + betaA_std * A_std + X * beta + betar_std * gps_std + phi + log_offset);
}
generated quantities {

   // Unit-level potential outcomes y_a and
   // average causal effect mu_a for a_new
   int<lower=0> y_a[n];
   real mu_a;
   
   for(i in 1:n) {
   
    // Generate unit-level potential outcome Y(a_new)
    y_a[i] = poisson_rng(exp(beta0_std + betaA_std*anew_std + X[i]*beta + betar_std*gps_std[i] + phi[i] + log_offset[i]));
  
  }
  // Compute average causal effect mu(a) = E[Y(a)]
  mu_a = mean(y_a);
  
  // Get regression coefficient estimates on original scale
  real beta_A = betaA_std / sd_A;
  real beta_r = betar_std / sd_gps;
  real beta0 = beta0_std - (betaA_std * mean_A / sd_A + betar_std * mean_gps / sd_gps);
}