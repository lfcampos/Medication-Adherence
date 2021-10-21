data {
  //data sizes
  int<lower=0> N;         // number of subjects
  int<lower=0> J;         // number of covariates observed
  int<lower=0> T[N];      // number of days followed per subject
  int<lower=0> maxT;      // maximum number of days observed
  int<lower=0> ni[N];     // number of observed outcomes per subject
  int<lower=0> n;         // total number of observed outcomes sum(ni)

  //covariates
  matrix[N, J] X;           // covariate matrix
  matrix[N, maxT] c;        // adherence matrix (ragged), 999 is missing value

  //outcomes
  vector[2*n] y;             // outcome vector
  int<lower=0> id[2*n];      // subject index for outcome vector
  int<lower=0> t[2*n];       // time index for outcome vector
  int<lower=0> yid[2*n];     // time index for outcome vector

}

parameters {
  vector<lower=0,upper=1>[2] rho;          // autoregressive parameter
  vector[2] phi;                           // adherence effect
  vector<lower=0,upper=30>[2] sigma;       // measurement error sd: [sd1, sd2]
  real<lower=-1,upper=1> cor;              // measurement error correlation
  vector<lower=0,upper=10>[2] sigma_nu;    // innovation variance scalings
  vector<lower=0,upper=30>[2] sigma_0;     // innitial state variance scalings
  matrix[J,2] beta;                        // covariate effect matrix
}

model {
  int pos;
  vector[2*n] alpha;
  vector[2*n] mu;
  matrix[2*n,2*n] Sigma;

  // // // // // // // // // // // //
  //            Priors             //
  // // // // // // // // // // // //
  for (i in 1:2){
    rho[i] ~ beta(3,1);            // auto-correlation
    phi[i] ~ normal(0,5);          // daily adherence effect

    sigma[i] ~ cauchy(0,5);      // measurement error standard deviation
    sigma_nu[i] ~ cauchy(0,2.5); // innovation standard deviation
    sigma_0[i] ~ cauchy(0,5);    // initial random effect standard deviation
  }
  cor ~ uniform(-1, 1);            // measurement error correlation

  for (j in 1:J){
    if(j == 1){
      beta[j,1] ~ normal(120, 20); // intercept for SBP
      beta[j,2] ~ normal(80, 20);  // intercept for DBP
    }else{
      beta[j,1] ~ normal(0, 20);   // coefficients
      beta[j,2] ~ normal(0, 20);
    }
  }

  // // // // // // // // // // // //
  // Calculate covariance matrix   //
  // // // // // // // // // // // //

  Sigma = rep_matrix(0, 2*n, 2*n);
  // diagonal
  for (i in 1:(2*n)){
    if(t[i] == 1){
      Sigma[i,i] = sigma_0[yid[i]]^2;
    }else{
      vector[t[i]] s;
      s[1] = rho[yid[i]]^(2*(t[i]-1)) * sigma_0[yid[i]]^2;
      for(j in 2:t[i]){
        s[j] = rho[yid[i]]^(2*(t[i]-j)) * sigma_nu[yid[i]]^2;
      }
      Sigma[i,i] = sum(s);
    }
  }
  // Calculate rest of covariance matrix
  for(i in 1:(2*n)){
    pos = 1;
    for (j in i:(i + ni[id[i]])){
      if(j <= 2*n){
        // covariance within observation vector
        if ((j > i) && (id[i] == id[j]) && (yid[i] == yid[j])){
          Sigma[i,j] = rho[yid[i]]^(pos) * Sigma[i,i];
          Sigma[j,i] = rho[yid[i]]^(pos) * Sigma[i,i];
          pos = pos + 1;
        }
        // Covariance between measurement errors
        if (j == (i + ni[id[i]])){
         Sigma[i,j] = cor*sigma[1]*sigma[2];
         Sigma[j,i] = cor*sigma[1]*sigma[2];
        }
      }
    }
    // Measurement error
    Sigma[i,i] = Sigma[i,i] + (sigma[yid[i]])^2;
  }

  // // // // // // // // // // // //
  //     Calculate mean vector     //
  // // // // // // // // // // // //

  for (i in 1:(2*n)){
    alpha[i] = 0;
    if(t[i]>1){
      vector[t[i]] s;
      s[1] = 0;
      for(j in 2:t[i]){
        s[j] = phi[yid[i]]*c[id[i],j]*rho[yid[i]]^(t[i]-j);
      }
      alpha[i] = sum(s);
    }
  }
  for (i in 1:(2*n)){
    mu[i] = X[id[i],]*beta[,yid[i]] + alpha[i];
  }

  // // // // // // // // // // // //
  //      Model Specification      //
  // // // // // // // // // // // //
  pos = 1;
  for (i in 1:N) {
    segment(y, pos, 2*ni[i]) ~ multi_normal(segment(mu, pos, 2*ni[i]), block(Sigma, pos, pos, 2*ni[i], 2*ni[i]));
    pos = pos + 2*ni[i];
  }
}

