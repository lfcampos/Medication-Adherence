data {
  int<lower=0> N;                        // number of training patients
  int<lower=0> N_test;                   // number of test patients
  int<lower=0> P;                        // number of covariate columns
  matrix[N,P] x;                         // covariates
  int total[N];                          // total number of days per patient
  int y[N];                              // number of successes per patient
}

// add intercept column
transformed data {
  vector[N] cons;
  matrix[N, P+1] X;

  cons = rep_vector(1, N);
  X = append_col(cons, x);
}

parameters {
  vector[N] delta;                    // random effects intercepts
  real<lower=0> sigma_delta;          // random effects sd
  vector[P+1] beta;                   // fixed effects coefficients
}

model {
  delta ~ normal(0, sigma_delta);

  for(i in 1:N){
   y[i] ~ binomial(total[i], inv_logit(delta[i] + X[i,] * beta));
  }
}

generated quantities {
  vector[N_test] delta_new;

  for(i in 1:N_test){
    delta_new[i] = normal_rng(0, sigma_delta);
  }

}
