# Rstan tutorial models

# Load the necessary packages
library(rstan);library(shinystan);library(dplyr)

# Set working directory (sets it where the script is located)
script.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script.dir)

# load the datafile
df <- read.csv("fake_SPI_birds.csv")

# Take a quick look at the structure of the data frame
str(df)
head(df)

# stan does not work with character classes, therefore transform all characters to integers 
# Use "//" to annotate text within rstan files insetad of "#"
# Denote the end of a line with ";"

# Change characters to indices
df <- df %>% mutate(Individual = Individual %>% as.factor() %>% as.numeric())
df <- df %>% mutate(Partner = Partner %>% as.factor() %>% as.numeric())
df <- df %>% mutate(Population = Population %>% as.factor() %>% as.numeric())
df$Sex <- ifelse(df$Sex == "F", 0, 1) 

head(df)

# Stan models have the following chunks
# However, not all are necessary for each model
write(
  temp <- "data{} 
  transformed data{}
  parameters {}
  transformed parameters {}
  model {}
  generated quantities {}"
  , file = "Mod.stan") 
# .stan is necessary to save it as a stan file!

# Model 1 Linear regression 
#===============================================================================
# Model 1.1: Simple regression model
#===============================================================================

# Prepare the data for the stan model
stan_data_mod1.1 <- list(N_z = nrow(df),
                        z = as.vector(scale(df$Exploration)),
                       x = as.vector(scale(df$Density)))

# Take a look at the imported data structure
View(stan_data_mod1.1)

# Write the stan model
write(
  temp <- "data{
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
  vector[N_z] e_z; // Predicted values for phenotype
  // model equation
  e_z = B_0 + B_1 * x;
  }
  model {
  // Likelihood function
  z ~ normal(e_z, sigma_e);
  // Priors
  B_0 ~ normal(0,1);
  B_1 ~ normal(0,1);
  sigma_e ~ exponential(3);
  }"
  , file = "Mod1.1.stan") 

# Fit the model
fit_mod1.1 <- stan("Mod1.1.stan", data = stan_data_mod1.1, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F)

# Look at the summary estimates
round(rstan::summary(fit_mod1.1)$summary[,c(1,4,6,8,9,10)],3)

# We're only interested in the intercept, slope and sd estimation for now
round(rstan::summary(fit_mod1.1, pars = c("B_0", "B_1", "sigma_e"))$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
# Model 1.2
# Assignment 1: Add the fixed effect Sex
#===============================================================================

# Prepare the data for the stan model
stan_data_mod1.2 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         x = as.vector(scale(df$Density)),
                         sex = df$Sex)

# Write the stan model
write(
  temp <- "data{
  int<lower=1> N_z; // Total number of observations of phenotype (Z)
  vector[N_z] z; // Phenotypic observations
  vector[N_z] x; // Environmental covariate
  vector[N_z] sex; // Sex of every individual
  }
  parameters {
  // Fixed effects
  real B_0; // Population intercept
  real B_1; // Population slope for covariate x
  real B_sex; // Population effect of sex
  // Random effects
  real<lower=0> sigma_e;// Standard deviation of model likelihood
  }
  transformed parameters {
  vector[N_z] e_z; // Predicted values for phenotype
  // model equation
  e_z = B_0 + B_1 * x + B_sex * sex;
  }
  model {
  // Likelihood function
  z ~ normal(e_z, sigma_e);
  // Priors
  B_0 ~ normal(0,1);
  B_1 ~ normal(0,1);
  B_sex ~ normal(0,1);
  sigma_e ~ exponential(3);
  }"
  , file = "Mod1.2.stan") 

# Fit the model
fit_mod1.2 <- stan("Mod1.2.stan", data = stan_data_mod1.2, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F)

# Look at the summary estimates
round(rstan::summary(fit_mod1.2, pars = c("B_0", "B_1", "B_sex", "sigma_e"))$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
# Model 1.3 - Add an interaction effect
#===============================================================================

# Prepare the data for the stan model
stan_data_mod1.3 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         x = as.vector(scale(df$Density)),
                         sex = df$Sex)

# used different notation to show alternative option for multiple betas

# Write the stan model
write(
  temp <- "data{
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
  }"
  , file = "Mod1.3.stan") 

# Fit the model
fit_mod1.3 <- stan("Mod1.3.stan", data = stan_data_mod1.3, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F)

# Look at the summary estimates
round(rstan::summary(fit_mod1.3, pars = c("B_0", "B", "B_int", "sigma_e"))$summary[,c(1,4,6,8,9,10)],3)

# Model 2 Generalised linear regression 
#===============================================================================
# Model 2.1: 
#===============================================================================

# Prepare the data for the stan model
stan_data_mod2.1 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         x = as.vector(scale(df$Density)),
                         N_I = length(unique(df$Individual)),
                         ID = df$Individual)

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
  }"
  , file = "Mod2.1.stan") 

# Fit the model
fit_mod2.1 <- stan("Mod2.1.stan", data = stan_data_mod2.1, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F)

pars_mod2.1 <- c("B_0", "B_1", "var_ID", "var_res", "var_P")

# Look at the summary estimates
round(rstan::summary(fit_mod2.1, pars = pars_mod2.1)$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
# Model 2.2
#===============================================================================

# Prepare the data for the stan model
stan_data_mod2.2 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         x = as.vector(scale(df$Density)),
                         N_I = length(unique(df$Individual)),
                         ID = df$Individual)

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

#===============================================================================
# Model 2.3
#===============================================================================

# Prepare the data for the stan model
stan_data_mod2.3 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         x = as.vector(scale(df$Density)),
                         N_I = length(unique(df$Individual)),
                         ID = df$Individual,
                         N_pop = length(unique(df$Population)),
                         pop = df$Population)

# Write the stan model
write(
  temp <- "data{
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
  }"
  , file = "Mod2.3.stan") 

# Fit the model
fit_mod2.3 <- stan("Mod2.3.stan", data = stan_data_mod2.3, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F)

pars_mod2.2 <- c("B_0", "B_1", 
                 "var_ID_int", "var_ID_slopes",
                 "var_pop_int", "var_pop_slopes",
                 "var_res", "var_P",
                 "Omega_I", "S_I",
                 "Omega_pop", "S_pop")

# Look at the summary estimates
round(rstan::summary(fit_mod2.3, pars = pars_mod2.2)$summary[,c(1,4,6,8,9,10)],3)

# Model 3 Mulivariate linear mixed model 
#===============================================================================
# Model 3.1: 
#===============================================================================

# Prepare the data for the stan model
stan_data_mod3.1 <- list(N_z = nrow(df),
                         zs = as.matrix(cbind(as.vector(scale(df$Exploration)),
                                              as.vector(scale(df$Aggression)))),
                         x = as.vector(scale(df$Density)),
                         N_I = length(unique(df$Individual)),
                         ID = df$Individual)

# Write the stan model
write(
  temp <- "data{
  int<lower=1> N_z; // Total number of observations of phenotype (Z)
  int<lower=1> N_I; // Total number of unique individuals (I)
  int<lower=1> ID[N_z]; //  Individual ID repeated obs
  array[N_z] vector[2] zs; // phenotypic observations for phenotype 1 & 2
  vector[N_z] x; // Environmental covariate
  }
  parameters {
  // Fixed effects
  vector[2] B_0; // Population intercept
  vector[2] B_1; // Population slope for covariate x
  // Random effects
  matrix[N_I,2] zI; // intercepts & slopes for each individual
  vector<lower=0>[2] sigma_I; // sd individual intercepts & slopes
  cholesky_factor_corr[2] L; // Factor to estimate correlations
  cholesky_factor_corr[2] LR; // Cholesky corr matrix for residuals
  vector<lower=0>[2] sd_R; // temporary indvidual standard deviation (i.e. residual variance)
  }
  transformed parameters {
  vector[N_z] e_z1; // Predicted values for phenotype
  vector[N_z] e_z2; // Predicted values for phenotype
  matrix[N_I,2] I  = zI * diag_pre_multiply(sigma_I, L)'; // get the unscaled value
  // model equations
  e_z1 = B_0[1] + B_1[1] * x + I[ID, 1];
  e_z2 = B_0[2] + B_1[2] * x + I[ID, 2];
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
  LR ~ lkj_corr_cholesky(3); // Prior residual correlation
  }
  generated quantities{
  // Variance
  real<lower=0> var_ID_int1 = sigma_I[1]^2;
  real<lower=0> var_res1 = sd_R[1]^2;
  real<lower=0> var_P1 = var_ID_int1 + var_res1;
  real<lower=0> var_ID_int2 = sigma_I[2]^2;
  real<lower=0> var_res2 = sd_R[2]^2;
  real<lower=0> var_P2 = var_ID_int2 + var_res2;
  // Correlations & Covariances
  matrix[2,2] Omega_I = L * L'; // Correlation matrix
  matrix[2,2] D_I = diag_matrix(sigma_I); // Diagonal SD matrix
  matrix[2,2] S_I = D_I*Omega_I*D_I; // Covariance matrix
  }"
  , file = "Mod3.1.stan") 

# Fit the model
fit_mod3.1 <- stan("Mod3.1.stan", data = stan_data_mod3.1, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F)

pars_mod3.1 <- c("B_0", "B_1", 
                 "var_ID_int1", "var_res1", "var_P1",
                 "var_ID_int2", "var_res2", "var_P2",
                 "Omega_I", "S_I")

# Look at the summary estimates
round(rstan::summary(fit_mod3.1, pars = pars_mod3.1)$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
# Model 3.2
#===============================================================================

# Prepare the data for the stan model
stan_data_mod3.2 <- list(N_z = nrow(df),
                         zs = as.matrix(cbind(as.vector(scale(df$Exploration)),
                                              as.vector(scale(df$Aggression)))),
                         x = as.vector(scale(df$Density)),
                         N_I = length(unique(df$Individual)),
                         ID = df$Individual)

# Write the stan model
write(
  temp <- "data{
  int<lower=1> N_z; // Total number of observations of phenotype (Z)
  int<lower=1> N_I; // Total number of unique individuals (I)
  int<lower=1> ID[N_z];  //  Individual ID repeated obs
  array[N_z] vector[2] zs;  // phenotypic observations for phenotype 1 & 2
  vector[N_z] x; // Environmental covariate
  }
  parameters {
  // Fixed effects
  vector[2] B_0; // Population intercept
  vector[2] B_1; // Population slope for covariate x
  // Random effects
  matrix[N_I,4] zI; // intercepts & slopes for each individual
  vector<lower=0>[4] sigma_I; // sd individual intercepts & slopes
  cholesky_factor_corr[4] L; // Factor to estimate correlations
  cholesky_factor_corr[2] LR; // Cholesky corr matrix for residuals
  vector<lower=0>[2] sd_R; // temporary indvidual standard deviation (i.e. residual variance)
  }
  transformed parameters {
  vector[N_z] e_z1; // Predicted values for phenotype
  vector[N_z] e_z2; // Predicted values for phenotype
  matrix[N_I,4] I  = zI * diag_pre_multiply(sigma_I, L)'; // get the unscaled value
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
  matrix[4,4] Omega_I = L * L'; // Correlation matrix
  matrix[4,4] D_I = diag_matrix(sigma_I); // Diagonal SD matrix
  matrix[4,4] S_I = D_I*Omega_I*D_I; // Covariance matrix
  }"
  , file = "Mod3.2.stan") 

# Fit the model
fit_mod3.2 <- stan("Mod3.2.stan", data = stan_data_mod3.2, 
                   chains = 4, iter = 4000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F)

pars_mod3.2 <- c("B_0", "B_1", 
                 "var_ID_int1", "var_ID_slopes1", "var_res1", "var_P1",
                 "var_ID_int2", "var_ID_slopes2", "var_res2", "var_P2",
                 "Omega_I", "S_I")

# Look at the summary estimates
round(rstan::summary(fit_mod3.2, pars = pars_mod3.2)$summary[,c(1,4,6,8,9,10)],3)

# Model 4 Error-in-variable (EIV) & selection gradient analysis 
#===============================================================================
# Model 4.1: z ~ B_0 + B1 * x
#===============================================================================

# Prepare the data for the stan model
stan_data_mod4.1 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         N_I = length(unique(df$Individual)),
                         ID = df$Individual,
                         w = as.vector(scale(unique(df$LRS))),
                         N_w = length(unique(df$LRS)),
                         ID_w = unique(df$Individual))

# Write the stan model
write(
  temp <- "data{
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
  }"
  , file = "Mod4.1.stan") 

# Fit the model
fit_mod4.1 <- stan("Mod4.1.stan", data = stan_data_mod4.1, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F)

pars_mod4.1 <- c("B_0", "Bw_0", "Bw",
                 "var_ID", "var_res", "var_P", "var_res_w")

# Look at the summary estimates
round(rstan::summary(fit_mod4.1, pars = pars_mod4.1)$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
# Model 4.2
#===============================================================================

# Prepare the data for the stan model
stan_data_mod4.2 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         N_I = length(unique(df$Individual)),
                         ID = df$Individual,
                         partnerID = df$Partner, 
                         w = as.vector(scale(df$Feeding_rate)))

# Write the stan model
write(
  temp <- "data{
  int<lower=1> N_z; // Total number of observations of phenotype (Z)
  int<lower=1> N_I; // Total number of unique individuals (I)
  int<lower=1> ID[N_z];  //  Individual ID repeated obs
  int<lower=1> partnerID[N_z];  //  Social partner ID repeated obs
  vector[N_z] z; // Phenotypic observations
  vector[N_z] w; // Fitness 
  }
  parameters {
  // Fixed effects
  real B_0; // Population intercept
  real Bw_0; // Fitness intercept
  vector[4] Bw; // Selection gradient slopes
  // Random effects
  vector[N_I] zI; // (intercepts for each individual)
  real<lower=0> sigma_I; // sd individual intercepts
  vector[N_I] zWI; // (fitness intercepts for each individual)
  real<lower=0> sigma_WI; // fitness sd individual intercepts
  real<lower=0> sigma_e;// Standard deviation of model likelihood
  real<lower=0> sigma_ew;// Standard deviation of fitness model likelihood
  }
  transformed parameters {
  vector[N_z] e_z; // Predicted values for phenotype
  vector[N_z] e_w; // Predicted values for fitness
  vector[N_I] I  =  zI * sigma_I; // Get the unscaled values for random effect
  vector[N_I] WI  =  zWI * sigma_WI; // Get the unscaled values for random effect
  // model equation 1 phenotypic model
  e_z = B_0 + I[ID];
  // model equation 2 selection gradient
  e_w = Bw_0 + Bw[1] * I[ID] +
               Bw[2] * I[partnerID] +
               Bw[3]*(I[ID].*I[ID]) + 
               Bw[3]*(I[ID].*I[ID]) +
               WI[ID];
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
  to_vector(zWI) ~ normal(0,1);
  sigma_WI ~ exponential(3);
  sigma_e ~ exponential(3);
  sigma_ew ~ exponential(3);
  }
  generated quantities{
  // Variance
  real<lower=0> var_ID = sigma_I^2;
  real<lower=0> var_res = sigma_e^2;
  real<lower=0> var_P = var_ID + var_res;
  real<lower=0> var_ID_w = sigma_WI^2;
  real<lower=0> var_res_w = sigma_ew^2;
  real<lower=0> var_P_w = var_ID_w + var_res_w;
  }"
  , file = "Mod4.2.stan") 

# Fit the model
fit_mod4.2 <- stan("Mod4.2.stan", data = stan_data_mod4.2, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F)

pars_mod4.2 <- c("B_0", "Bw_0", "Bw",
                 "var_ID", "var_res", "var_P", 
                 "var_ID_w", "var_res_w", "var_P_w")

# Look at the summary estimates
round(rstan::summary(fit_mod4.2, pars = pars_mod4.2)$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
# Model 4.3
#===============================================================================

# Prepare the data for the stan model
stan_data_mod4.3 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         x = as.vector(scale(df$Density)),
                         N_I = length(unique(df$Individual)),
                         ID = df$Individual,
                         w = as.vector(scale(unique(df$LRS))),
                         N_w = length(unique(df$LRS)),
                         ID_w = unique(df$Individual))

# Write the stan model
write(
  temp <- "data{
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
  }"
  , file = "Mod4.3.stan") 

# Fit the model
fit_mod4.3 <- stan("Mod4.3.stan", data = stan_data_mod4.3, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F)

pars_mod4.3 <- c("B_0", "Bw_0", "Bw",
                 "var_ID_int", "var_ID_slopes", "var_res", "var_P", "var_res_w",
                 "Omega_I", "S_I")

# Look at the summary estimates
round(rstan::summary(fit_mod4.3, pars = pars_mod4.3)$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
# Bonus: working with shiny stan & posterior distribution 
#===============================================================================

launch_shinystan(fit_mod4.3)

# Still have to do this


