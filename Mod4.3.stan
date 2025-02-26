data{
  int<lower=1> N_z; // Total number of observations of phenotype (Z)
  int<lower=1> N_I; // Total number of unique individuals (I)
  int<lower=1> N_w; // Total number of fitness observations
  int<lower=1> ID[N_z];  //  Individual ID repeated obs
  int<lower=1> ID_w[N_w];  //  Individual ID repeated obs
  vector[N_z] z; // Phenotypic observations
  vector[N_z] x; // Environmental covariate
  vector[N_w] w; // Fitness 
  }
  parameters {
  // Fixed effects
  real B_0; // Population intercept
  real B_1; // covariate slope
  real Bw_0; // Fitness intercept
  vector[5] Bw; // Selection gradient slopes
  // Random effects
  matrix[N_I,2] zI; // intercepts & slopes for each individual
  vector<lower=0>[2] sigma_I; // sd individual intercepts & slopes
  cholesky_factor_corr[2] L; // Factor to estimate correlation int-slopes, is size n-trait  real<lower=0> sigma_e;// Standard deviation of model likelihood
  real<lower=0> sigma_e;// Standard deviation of phenotypic model likelihood
  real<lower=0> sigma_ew;// Standard deviation of fitness model likelihood
  }
  transformed parameters {
  vector[N_z] e_z; // Predicted values for phenotype
  vector[N_w] e_w; // Predicted values for fitness
  matrix[N_I,2] I = zI * diag_pre_multiply(sigma_I, L)'; // get the unscaled value  // model equation 1 phenotypic model
  e_z = B_0 + (B_1 + I[ID,2]) .* x + I[ID, 1];
  // model equation 2 selection gradient
  e_w = Bw_0 + Bw[1] * I[ID_w, 1] +
               Bw[2] * I[ID_w, 2] +
               Bw[2]*(I[ID_w, 1].*I[ID_w, 1]) +
               Bw[2]*(I[ID_w, 2].*I[ID_w, 2]) +
               Bw[2]*(I[ID_w, 1].*I[ID_w, 2]);
  }
  model {
  // Likelihood functions
  z ~ normal(e_z, sigma_e);
  w ~ normal(e_w, sigma_ew);
  // Priors
  B_0 ~ normal(0,1);
  B_1 ~ normal(0,1);
  Bw_0 ~ normal(0,1);
  Bw ~ normal(0,1);
  to_vector(zI) ~ normal(0,1);
  sigma_I ~ exponential(3);
  L ~ lkj_corr_cholesky(3);  sigma_e ~ exponential(3);
  sigma_ew ~ exponential(3);
  }
  generated quantities{
  // Variance
  real<lower=0> var_ID_int = sigma_I[1]^2;
  real<lower=0> var_ID_slopes = sigma_I[2]^2;
  real<lower=0> var_res = sigma_e^2;
  real<lower=0> var_P = var_ID_int + var_ID_slopes + var_res;
  real<lower=0> var_res_w = sigma_ew^2;
  // Correlations & Covariances
  matrix[2,2] Omega_I = L * L'; // Correlation matrix
  matrix[2,2] D_I = diag_matrix(sigma_I); // Diagonal SD matrix
  matrix[2,2] S_I = D_I*Omega_I*D_I; // Covariance matrix
  }
