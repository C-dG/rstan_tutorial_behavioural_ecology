data{
  int<lower=1> N_z; // Total number of observations of phenotype (Z)
  int<lower=1> N_I; // Total number of unique individuals (I)
  int<lower=1> ID[N_z];  //  Individual ID repeated obs
  vector[N_z] z; // Phenotypic observations
  vector[N_z] x; // Environmental covariate
  }
  parameters {
  // Fixed effects
  real B_0; // Population intercept
  real B_1; // Population slope for covariate x
  // Random effects
  vector[N_I] zI; // (intercepts for each individual)
  real<lower=0> sigma_I; // sd individual intercepts
  real<lower=0> sigma_e;// Standard deviation of model likelihood
  }
  transformed parameters {
  vector[N_z] e_z; // Predicted values for phenotype
  vector[N_I] I  =  zI * sigma_I; // Get the unscaled values for random effect
  // model equation
  e_z = B_0 + B_1 * x + I[ID];
  }
  model {
  // Likelihood function
  z ~ normal(e_z, sigma_e);
  // Priors
  B_0 ~ normal(0,1);
  B_1 ~ normal(0,1);
  to_vector(zI) ~ normal(0,1);
  sigma_I ~ exponential(3);
  sigma_e ~ exponential(3);
  }
  generated quantities{
  // Variance
  real<lower=0> var_ID = sigma_I^2;
  real<lower=0> var_res = sigma_e^2;
  real<lower=0> var_P = var_ID + var_res;
  }
