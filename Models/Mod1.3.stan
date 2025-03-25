data{
  int<lower=1> N_z; // Total number of observations of phenotype (Z)
  vector[N_z] z; // Phenotypic observations
  vector[N_z] x; // Environmental covariate
  vector[N_z] sex; // Sex of every individual
  }
  parameters {
  // Fixed effects
  real B_0; // Population intercept
  vector[2] B; // Vector of betas for fixed effects
  vector[1] B_int; // Vector of n interactions
  // Random effects
  real<lower=0> sigma_e;// Standard deviation of model likelihood
  }
  transformed parameters {
  vector[N_z] e_z; // Predicted values for phenotype
  // model equation
  e_z = B_0 + B[1]* x + B[2] * sex + B_int[1].*x.*sex;
  }
  model {
  // Likelihood function
  z ~ normal(e_z, sigma_e);
  // Priors
  B_0 ~ normal(0,1);
  B ~ normal(0,1);
  B_int ~ normal(0,1);
  sigma_e ~ exponential(3);
  }
