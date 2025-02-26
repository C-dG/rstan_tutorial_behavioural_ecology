data{
  int<lower=1> N_z; // Total number of observations of phenotype (Z)
  int<lower=1> N_I; // Total number of unique individuals (I)
  int<lower=1> N_pop; // Total number of unique populations (pop)
  int<lower=1> ID[N_z];  //  Individual ID repeated obs
  int<lower=1> pop[N_z];  //  population repeated obs
  vector[N_z] z; // Phenotypic observations
  vector[N_z] x; // Environmental covariate
  }
  parameters {
  // Fixed effects
  real B_0; // Population intercept
  real B_1; // Population slope for covariate x
  // Random effects
  matrix[N_I,2] zI; // intercepts & slopes for each individual
  vector<lower=0>[2] sigma_I; // sd individual intercepts & slopes
  cholesky_factor_corr[2] L; // Factor to estimate correlation int-slopes, is size n-trait
  matrix[N_pop,2] zpop; // intercepts & slopes for each individual
  vector<lower=0>[2] sigma_pop; // sd individual intercepts & slopes
  cholesky_factor_corr[2] L_pop; // Factor to estimate correlation int-slopes, is size n-trait
  real<lower=0> sigma_e;// Standard deviation of model likelihood
  }
  transformed parameters {
  vector[N_z] e_z; // Predicted values for phenotype
  matrix[N_I,2] I  = zI * diag_pre_multiply(sigma_I, L)'; // get the unscaled value
  matrix[N_pop,2] P  = zpop * diag_pre_multiply(sigma_pop, L_pop)'; // get the unscaled value
  // model equation
  e_z = B_0 + (B_1 + I[ID,2] + P[pop,2]) .* x + I[ID, 1] + P[pop,1];
  }
  model {
  // Likelihood function
  z ~ normal(e_z, sigma_e);
  // Priors
  B_0 ~ normal(0,1);
  B_1 ~ normal(0,1);
  to_vector(zI) ~ normal(0,1);
  sigma_I ~ exponential(3);
  L ~ lkj_corr_cholesky(3);
  to_vector(zpop) ~ normal(0,1);
  sigma_pop ~ exponential(3);
  L_pop ~ lkj_corr_cholesky(3);
  sigma_e ~ exponential(3);
  }
  generated quantities{
  // Variance
  real<lower=0> var_ID_int = sigma_I[1]^2;
  real<lower=0> var_ID_slopes = sigma_I[2]^2;
  real<lower=0> var_pop_int = sigma_pop[1]^2;
  real<lower=0> var_pop_slopes = sigma_pop[2]^2;
  real<lower=0> var_res = sigma_e^2;
  real<lower=0> var_P = var_ID_int + var_ID_slopes + var_pop_int + var_pop_slopes + var_res;
  // Correlations & Covariances
  matrix[2,2] Omega_I = L * L'; // Correlation matrix
  matrix[2,2] D_I = diag_matrix(sigma_I); // Diagonal SD matrix
  matrix[2,2] S_I = D_I*Omega_I*D_I; // Covariance matrix
  // For population level effects
  matrix[2,2] Omega_pop = L_pop * L_pop'; // Correlation matrix
  matrix[2,2] D_pop = diag_matrix(sigma_pop); // Diagonal SD matrix
  matrix[2,2] S_pop = D_pop*Omega_pop*D_pop; // Covariance matrix
  }
