// Authors: Evan Funderburg, Emily Sitnik, Mark Ritchie, and Ray Bai
// Suitably modified from https://mc-stan.org/learn-stan/case-studies/mbjoseph-CARStan.html

data {
  int<lower=1> n;
  int<lower=1> p;
  matrix[n, p] X;
  vector[n] A;
  matrix<lower=0, upper=1>[n, n] W;
}
transformed data {

  matrix<lower=0>[n, n] D;
  {
    vector[n] W_rowsums;
    for (i in 1:n) {
      W_rowsums[i] = sum(W[i, ]);
    }
    D = diag_matrix(W_rowsums);
  }
}
parameters {
  real gamma0;
  vector[p] gamma;
  real<lower=0> sigma_e;
  real<lower=0, upper=1> rho_e;
}
transformed parameters {
  vector[n] mu;
  vector[n] diag_var;
  vector[n] log_pdf;
  vector<lower=0>[n] gps;
  
  diag_var = square(sigma_e) * diagonal(inverse_spd(D - rho_e * W));
  
  // Store the generalized propensity scores
  for (i in 1:n) {
    mu[i] = gamma0 + X[i] * gamma;
    log_pdf[i] = normal_lpdf(A[i] | mu[i], sqrt(diag_var[i]));
    gps[i] = exp(log_pdf[i]);
  }
}
model {
  // Priors 
  gamma0 ~ normal(0, 5);
  gamma ~ normal(0, 5);
  sigma_e ~ cauchy(0, 2.5);

  // Likelihood
  target += multi_normal_prec_lpdf(A | mu, (D - rho_e * W) / square(sigma_e));
}
