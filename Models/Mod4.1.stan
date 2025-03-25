data{
  int<lower=1> N_z; // Total number of observations of phenotype (Z)
  int<lower=1> N_I; // Total number of unique individuals (I)
  int<lower=1> N_w; // Total number of fitness observations
  int<lower=1> ID[N_z];  //  Individual ID repeated obs
  int<lower=1> ID_w[N_w];  //  Individual ID repeated obs
  vector[N_z] z; // Phenotypic observations
  vector[N_w] w; // Fitness 
  }
  parameters {
  // Fixed effects
  real B_0; // Population intercept
  real Bw_0; // Fitness intercept
  vector[2] Bw; // Selection gradient slopes
  // Random effects
  vector[N_I] zI; // (intercepts for each individual)
  real<lower=0> sigma_I; // sd individual intercepts
  real<lower=0> sigma_e;// Standard deviation of model likelihood
  real<lower=0> sigma_ew;// Standard deviation of fitness model likelihood
  }
  transformed parameters {
  vector[N_z] e_z; // Predicted values for phenotype
  vector[N_w] e_w; // Predicted values for fitness
  vector[N_I] I  =  zI * sigma_I; // Get the unscaled values for random effect
  // model equation 1 phenotypic model
  e_z = B_0 + I[ID];
  // model equation 2 selection gradient
  e_w = Bw_0 + Bw[1] * I[ID_w] + Bw[2]*(I[ID_w].*I[ID_w]);
  }
  model {
  // Likelihood functions
  z ~ normal(e_z, sigma_e);
  w ~ normal(e_w, sigma_ew);
  // Priors
  B_0 ~ normal(0,1);
  Bw_0 ~ normal(0,1);
  Bw ~ normal(0,1);
  to_vector(zI) ~ normal(0,1);
  sigma_I ~ exponential(3);
  sigma_e ~ exponential(3);
  sigma_ew ~ exponential(3);
  }
  generated quantities{
  // Variance
  real<lower=0> var_ID = sigma_I^2;
  real<lower=0> var_res = sigma_e^2;
  real<lower=0> var_P = var_ID + var_res;
  real<lower=0> var_res_w = sigma_ew^2;
  }
