#===============================================================================
# Model 2.2
# Add a random slope 
#===============================================================================

# Prepare the data for the stan model
stan_data_mod2.2 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         x = as.vector(scale(df$Density)),
                         N_I = length(unique(df$Individual)),
                         ID = as.integer(df$Individual))

# Write the stan model
write(
  temp <- "data{
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
  matrix[N_I,2] zI; // intercepts & slopes for each individual
  vector<lower=0>[2] sigma_I; // sd individual intercepts & slopes
  cholesky_factor_corr[2] L; // Factor to estimate correlation int-slopes, is size n-trait
  real<lower=0> sigma_e;// Standard deviation of model likelihood
  }
  transformed parameters {
  vector[N_z] e_z; // Predicted values for phenotype
  matrix[N_I,2] I = zI * diag_pre_multiply(sigma_I, L)'; // get the unscaled value
  // model equation
  e_z = B_0 + (B_1 + I[ID,2]) .* x + I[ID, 1];
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
  sigma_e ~ exponential(3);
  }
  generated quantities{
  // Variance
  real<lower=0> var_ID_int = sigma_I[1]^2;
  real<lower=0> var_ID_slopes = sigma_I[2]^2;
  real<lower=0> var_res = sigma_e^2;
  real<lower=0> var_P = var_ID_int + var_ID_slopes + var_res;
  // Correlations & Covariances
  matrix[2,2] Omega_I = L * L'; // Correlation matrix
  matrix[2,2] D_I = diag_matrix(sigma_I); // Diagonal SD matrix
  matrix[2,2] S_I = D_I*Omega_I*D_I; // Covariance matrix
  }"
  , file = "Mod2.2.stan") 

# Fit the model
fit_mod2.2 <- stan("Mod2.2.stan", data = stan_data_mod2.2, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F)

pars_mod2.2 <- c("B_0", "B_1", 
                 "var_ID_int", "var_ID_slopes", "var_res", "var_P",
                 "Omega_I", "S_I")

# Look at the summary estimates
round(rstan::summary(fit_mod2.2, pars = pars_mod2.2)$summary[,c(1,4,6,8,9,10)],3)