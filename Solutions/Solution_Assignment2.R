#===============================================================================
# Model 2.3
# Assignment 2: Add a random intercept & slope for population
#===============================================================================

# Prepare the data for the stan model
stan_data_mod2.3 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         x = as.vector(scale(df$Density)),
                         N_I = length(unique(df$Individual)),
                         ID = as.integer(as.factor(df$Individual)),
                         N_pop = length(unique(df$Population)),
                         pop = as.integer(as.factor(df$Population)))

# Write the stan model
write(
  temp <- "data{
  int<lower=1> N_z;         // Total number of observations of phenotype (Z)
  int<lower=1> N_I;         // Total number of unique individuals (I)
  int<lower=1> N_pop;       // Total number of unique populations (pop)
  int<lower=1> ID[N_z];     //  Individual ID repeated obs
  int<lower=1> pop[N_z];    //  population repeated obs
  vector[N_z] z;            // Phenotypic observations
  vector[N_z] x;            // Environmental covariate
  }
  parameters {
  // Fixed effects
  real B_0; 
  real B_1; 
  
  // Random effects
  matrix[N_I,2] zI; 
  vector<lower=0>[2] sigma_I; 
  cholesky_factor_corr[2] L; 
  matrix[N_pop,2] zpop; // Population-level intercepts & slopes 
  vector<lower=0>[2] sigma_pop; // Population-level sd's intercepts & slopes
  cholesky_factor_corr[2] L_pop; // Cholesky factor for population-level correlations
  real<lower=0> sigma_e;
  }
  
  transformed parameters {
  vector[N_z] e_z; 
  
  matrix[N_I,2] I  = zI * diag_pre_multiply(sigma_I, L)'; 
  matrix[N_pop,2] P  = zpop * diag_pre_multiply(sigma_pop, L_pop)'; // get the unscaled values for population-level effects
  
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
  to_vector(zpop) ~ normal(0,1);  // Prior for standardised population-level effects
  sigma_pop ~ exponential(3);     // Priors sd population-level effects
  L_pop ~ lkj_corr_cholesky(3);   // Correlation prior population-level effects
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
  matrix[2,2] rho_I = L * L'; 
  matrix[2,2] cov_I = diag_matrix(sigma_I)*rho_I*diag_matrix(sigma_I); 
  
  // For population level effects
  matrix[2,2] rho_pop = L_pop * L_pop'; // Correlation matrix
  matrix[2,2] cov_pop = diag_matrix(sigma_pop)*rho_pop*diag_matrix(sigma_pop); // Covariance matrix
  }"
  , file = "Models/Mod2.3.stan") 

# Fit the model
fit_mod2.3 <- stan("Models/Mod2.3.stan", data = stan_data_mod2.3, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F,
                   seed = 26032025)

pars_mod2.2 <- c("B_0", "B_1", 
                 "var_ID_int", "var_ID_slopes",
                 "var_pop_int", "var_pop_slopes",
                 "var_res", "var_P",
                 "Omega_I", "S_I",
                 "Omega_pop", "S_pop")

# Look at the summary estimates
round(rstan::summary(fit_mod2.3, pars = pars_mod2.2)$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
