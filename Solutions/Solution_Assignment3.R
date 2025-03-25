#===============================================================================
# Model 3.2
# Assignment 3: Add random slopes for individuals and correlate them across traits
#===============================================================================

# Prepare the data for the stan model
stan_data_mod3.2 <- list(N_z = nrow(df),
                         zs = as.matrix(cbind(as.vector(scale(df$Exploration)),
                                              as.vector(scale(df$Aggression)))),
                         x = as.vector(scale(df$Density)),
                         N_I = length(unique(df$Individual)),
                         ID = as.integer(as.factor(df$Individual)))

# Write the stan model
write(
  temp <- "data{
  int<lower=1> N_z;           // Total number of observations of phenotype (Z)
  int<lower=1> N_I;           // Total number of unique individuals (I)
  int<lower=1> ID[N_z];       // Individual ID repeated obs
  array[N_z] vector[2] zs;    // Phenotypic observations for phenotype 1 & 2
  vector[N_z] x;              // Environmental covariate
  }
  
  parameters {
  // Fixed effects
  vector[2] B_0; 
  vector[2] B_1; 
  
  // Random effects
  matrix[N_I,4] zI;             // Note the number of dimension = 4 
  vector<lower=0>[4] sigma_I;   // Note the number of dimension = 4 
  cholesky_factor_corr[4] L;    // Note the number of dimension = 4 
  cholesky_factor_corr[2] LR; 
  vector<lower=0>[2] sd_R; 
  }
 
   transformed parameters {
  vector[N_z] e_z1; 
  vector[N_z] e_z2; 
  
  matrix[N_I,4] I  = zI * diag_pre_multiply(sigma_I, L)'; // Note the number of dimension = 4 
  
  // model equations
  e_z1 = B_0[1] + (B_1[1] + I[ID, 1]) .* x + I[ID, 2];
  e_z2 = B_0[2] + (B_1[2] + I[ID, 3]) .* x + I[ID, 4];
  }
 
  model {
  // Transform expected values to array
  array[N_z] vector[2] mus;
  for (o in 1:N_z) 
  mus[o] = [e_z1[o], e_z2[o]]';
  matrix[2,2] L_sigma = diag_pre_multiply(sd_R, LR);
  
  // Cholesky multinormal likelihood function to estimate residual correlations
  zs ~ multi_normal_cholesky(mus, L_sigma);  
  
  // Priors
  B_0 ~ normal(0,1);
  B_1 ~ normal(0,1);
  to_vector(zI) ~ normal(0,1);
  sigma_I ~ exponential(3);
  L ~ lkj_corr_cholesky(3);
  LR ~ lkj_corr_cholesky(3); //Prior residual correlation
  }
  
  generated quantities{
  // Variance
  real<lower=0> var_ID_int1 = sigma_I[1]^2;
  real<lower=0> var_ID_slopes1 = sigma_I[2]^2;
  real<lower=0> var_res1 = sd_R[1]^2;
  real<lower=0> var_P1 = var_ID_int1 + var_ID_slopes1 + var_res1;
  real<lower=0> var_ID_int2 = sigma_I[3]^2;
  real<lower=0> var_ID_slopes2 = sigma_I[4]^2;
  real<lower=0> var_res2 = sd_R[2]^2;
  real<lower=0> var_P2 = var_ID_int2 + var_ID_slopes2 + var_res2;
  
  // Correlations & Covariances
  matrix[4,4] Omega_I = L * L'; 
  matrix[4,4] D_I = diag_matrix(sigma_I); 
  matrix[4,4] S_I = D_I*Omega_I*D_I; 
  }"
  , file = "Models/Mod3.2.stan") 

# Fit the model
fit_mod3.2 <- stan("Models/Mod3.2.stan", data = stan_data_mod3.2, 
                   chains = 4, iter = 4000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F,
                   seed = 26032025)

pars_mod3.2 <- c("B_0", "B_1", 
                 "var_ID_int1", "var_ID_slopes1", "var_res1", "var_P1",
                 "var_ID_int2", "var_ID_slopes2", "var_res2", "var_P2",
                 "Omega_I", "S_I")

# Look at the summary estimates
round(rstan::summary(fit_mod3.2, pars = pars_mod3.2)$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
