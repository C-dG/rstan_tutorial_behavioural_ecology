data{
  int<lower=1> N_z; // Total number of observations of phenotype (Z)
  vector[N_z] z; // Phenotypic observations
  vector[N_z] x; // Environmental covariate
  }
  parameters {
  // Fixed effects
  real B_0; // Population intercept
  real B_1; // Population slope for covariate x
  // Random effects
  real<lower=0> sigma_e;// Standard deviation of model likelihood
  }
  transformed parameters {
  }
  model {
    vector[N_z] e_z; // Predicted values for phenotype
  // model equation
  e_z = B_0 + B_1 * x;

  // Likelihood function
  z ~ normal(e_z, sigma_e);
  // Priors
  B_0 ~ normal(0,1);
  B_1 ~ normal(0,1);
  sigma_e ~ exponential(3);
  }
