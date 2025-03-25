#===============================================================================
# Model 4.3
# Assignment 4: Add individual random slopes to model 4.1 (not 4.2!)
# Setup up direct, quadratic and correlation selection gradients for individual intercepts and slopes 
#===============================================================================#===============================================================================

# Prepare the data for the stan model
stan_data_mod4.3 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         x = as.vector(scale(df$Density)),
                         N_I = length(unique(df$Individual)),
                         ID = as.integer(as.factor(df$Individual)),
                         w = as.vector(scale(unique(df$LRS))),
                         N_w = length(unique(df$LRS)),
                         ID_w = unique(as.integer(as.factor(df$Individual))))

# Write the stan model
write(
  temp <- "data{
  int<lower=1> N_z;         // Total number of observations of phenotype (Z)
  int<lower=1> N_I;         // Total number of unique individuals (I)
  int<lower=1> N_w;         // Total number of fitness observations
  int<lower=1> ID[N_z];     // Individual ID repeated obs
  int<lower=1> ID_w[N_w];   // Individual ID repeated obs
  vector[N_z] z;            // Phenotypic observations
  vector[N_z] x;            // Environmental covariate
  vector[N_w] w;            // Fitness 
  }
  
  parameters {
  // Fixed effects
  real B_0; 
  real B_1; 
  real Bw_0; 
  vector[5] Bw;     // Note the number of dimensions = 5
  
  // Random effects
  matrix[N_I,2] zI; 
  vector<lower=0>[2] sigma_I; 
  cholesky_factor_corr[2] L; 
  real<lower=0> sigma_e;
  real<lower=0> sigma_ew;
  }
  
  transformed parameters {
  vector[N_z] e_z; // Predicted values for phenotype
  vector[N_w] e_w; // Predicted values for fitness
  
  matrix[N_I,2] I = zI * diag_pre_multiply(sigma_I, L)'; // get the unscaled value  // model equation 1 phenotypic model
  
  e_z = B_0 + (B_1 + I[ID,2]) .* x + I[ID, 1];
  // model equation 2 selection gradient
  e_w = Bw_0 + Bw[1] * I[ID_w, 1] +
               Bw[2] * I[ID_w, 2] +
               Bw[3]*(I[ID_w, 1].*I[ID_w, 1]) +
               Bw[4]*(I[ID_w, 2].*I[ID_w, 2]) +
               Bw[5]*(I[ID_w, 1].*I[ID_w, 2]);
  // Note the estimation of correlational selection Bw[5]*(I[ID_w, 1].*I[ID_w, 2])
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
  L ~ lkj_corr_cholesky(3); 
  sigma_e ~ exponential(3);
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
  matrix[2,2] Omega_I = L * L'; /
  matrix[2,2] D_I = diag_matrix(sigma_I);
  matrix[2,2] S_I = D_I*Omega_I*D_I; 
  }"
  , file = "Models/Mod4.3.stan") 

# Fit the model
fit_mod4.3 <- stan("Models/Mod4.3.stan", data = stan_data_mod4.3, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F,
                   seed = 26032025)

pars_mod4.3 <- c("B_0", "Bw_0", "Bw",
                 "var_ID_int", "var_ID_slopes", "var_res", "var_P", "var_res_w",
                 "Omega_I", "S_I")

# Look at the summary estimates
round(rstan::summary(fit_mod4.3, pars = pars_mod4.3)$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================#===============================================================================
