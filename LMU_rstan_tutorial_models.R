# Rstan tutorial models

# Install rstan if necessary: 

# run the next line if you already have rstan installed
# remove.packages(c("StanHeaders", "rstan"))

install.packages("rstan", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# Load the necessary packages
library(rstan);library(shinystan)

rstan::stan_version()# We used version"2.32.2"

# Set working directory (sets it where the script is located)
#script.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(script.dir)

# load the datafile
df <- read.csv("fake_SPI_birds.csv")

# Take a quick look at the structure of the data frame
head(df)

# stan does not work with character classes, therefore transform all characters to integers 
# Use "//" to annotate text within rstan files insetad of "#"
# Denote the end of a line with ";"


# Stan models have the following chunks
# However, not all are necessary for each model
write(
  temp <- "data{} 
  transformed data{}
  parameters {}
  transformed parameters {}
  model {}
  generated quantities {}"
  , file = "Models/Mod.stan") 
# .stan is necessary to save it as a stan file!

# Model 1 Linear regression 
#===============================================================================
# Model 1.1: 
# Simple regression model
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
  int<lower=1> N_z;   // Total number of observations of phenotype (Z)
  vector[N_z] z;     // Phenotypic observations
  vector[N_z] x;    // Environmental covariate
  }
  
  parameters {
  // Fixed effects
  real B_0;               // Population intercept
  real B_1;               // Population slope for covariate x
  
  // Random effects
  real<lower=0> sigma_e;  // Standard deviation of model likelihood
  }
  
  transformed parameters {
  vector[N_z] e_z;        // Predicted values for phenotype
  
  // model equation
  e_z = B_0 + B_1 * x;
  // The model equation can also be specified in the model block, 
  // but this way you can save the predicted values to derive the residuals
  }
  
  model {
  // Likelihood function
  z ~ normal(e_z, sigma_e);
  // Alternatively: z ~ normal(B_0 + B_1 * x, sigma_e)
  
  // Priors
  B_0 ~ normal(0,1);        // Because the data is standardised this represents the mean and the sd   
  B_1 ~ normal(0,1);        // Idem
  sigma_e ~ exponential(3); // This is a half prior, because sigma_e is constraint to be >0
  }"
  , file = "Models/Mod1.1.stan") 

# Fit the model
fit_mod1.1 <- stan("Models/Mod1.1.stan", data = stan_data_mod1.1, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F,
                   seed = 26032025)

# Get the random seed of the model for reproducibility,
# Because reproducibility is important!!!
rstan::get_seed(fit_mod1.1)

# Look at the summary estimates
round(rstan::summary(fit_mod1.1)$summary[,c(1,4,6,8,9,10)],3)

# We're only interested in the intercept, slope and sd estimation for now
round(rstan::summary(fit_mod1.1, pars = c("B_0", "B_1", "sigma_e"))$summary[,c(1,4,6,8,9,10)],3)

# To save the model or summary to work with it later:
saveRDS(fit_mod1.1, "fit_mod1.1.RDS")

#===============================================================================
# Model 1.2
# Assignment 1: Add the fixed effect Sex
#===============================================================================

# Tip: copy-paste model 1.1 here and make edits

# Having trouble? See folder Solutions

#===============================================================================
# Model 1.3 
# Add an interaction effect
#===============================================================================

# Prepare the data for the stan model
stan_data_mod1.3 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         x = as.vector(scale(df$Density)),
                         sex = ifelse(df$Sex == "F", -0.5, 0.5))

# used different notation to show alternative option for multiple betas

# Write the stan model
write(
  temp <- "data{
  int<lower=1> N_z;   // Total number of observations of phenotype (Z)
  vector[N_z] z;      // Phenotypic observations
  vector[N_z] x;      // Environmental covariate
  vector[N_z] sex;    // Sex of every individual
  }
  
  parameters {
  // Fixed effects
  real B_0;                 
  vector[2] B;              // Vector of betas for fixed effects
  vector[1] B_int;          // Vector of n interactions
  // You can also use this vector notation to minimise code

  // Random effects
  real<lower=0> sigma_e;    // Standard deviation of model likelihood
  }
  
  transformed parameters {
  vector[N_z] e_z;          
  
  // model equation
  e_z = B_0 + B[1]* x + B[2] * sex + B_int[1].*x.*sex;
  // Mind the .* when working with vectors
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
  , file = "Models/Mod1.3.stan") 

# Fit the model
fit_mod1.3 <- stan("Models/Mod1.3.stan", data = stan_data_mod1.3, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F,
                   seed = 26032025)

# Look at the summary estimates
round(rstan::summary(fit_mod1.3, pars = c("B_0", "B", "B_int", "sigma_e"))$summary[,c(1,4,6,8,9,10)],3)

# Model 2 Generalised linear regression 
#===============================================================================
# Model 2.1: 
# Add a random intercept
#===============================================================================

# Prepare the data for the stan model
stan_data_mod2.1 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         x = as.vector(scale(df$Density)),
                         N_I = length(unique(df$Individual)),
                         ID = as.integer(as.factor(df$Individual)))

# Write the stan model
write(
  temp <- "data{
  int<lower=1> N_z;       // Total number of observations of phenotype (Z)
  int<lower=1> N_I;       // Total number of unique individuals (I)
  int<lower=1> ID[N_z];   //  Individual ID repeated obs
  vector[N_z] z;          // Phenotypic observations
  vector[N_z] x;          // Environmental covariate
  }
  
  parameters {
  // Fixed effects
  real B_0; 
  real B_1; 
  
  // Random effects
  vector[N_I] zI;           // Intercepts for each individual
  real<lower=0> sigma_I;    // sd individual intercepts
  real<lower=0> sigma_e;    // Standard deviation of model likelihood
  }
  
  transformed parameters {
  vector[N_z] e_z; 
  
  // This is an optimisation trick in stan for random effect computation
  vector[N_I] I  =  zI * sigma_I; // Get the unscaled values for random effect
  
  // model equation
  e_z = B_0 + B_1 * x + I[ID];
  // Note the notation for random effects i.e. + I[ID]
  }
  
  model {
  // Likelihood function
  z ~ normal(e_z, sigma_e);
  
  // Priors
  B_0 ~ normal(0,1);
  B_1 ~ normal(0,1);
  to_vector(zI) ~ normal(0,1);    // Prior for the standardised individual intercepts
  sigma_I ~ exponential(3);       // Prior for the sd of individual intercepts
  sigma_e ~ exponential(3);
  }
  
  generated quantities{
  // Variance
  real<lower=0> var_ID = sigma_I^2;         // Get the indivual intercept variance
  real<lower=0> var_res = sigma_e^2;        // Get the residual variance
  real<lower=0> var_P = var_ID + var_res;   // Get the total estimated variance
  }"
  , file = "Models/Mod2.1.stan") 

# Fit the model
fit_mod2.1 <- stan("Models/Mod2.1.stan", data = stan_data_mod2.1, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F,
                   seed = 26032025)

pars_mod2.1 <- c("B_0", "B_1", "var_ID", "var_res", "var_P")

# Look at the summary estimates
round(rstan::summary(fit_mod2.1, pars = pars_mod2.1)$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
# Model 2.2
# Add a random slope 
#===============================================================================

# Prepare the data for the stan model
stan_data_mod2.2 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         x = as.vector(scale(df$Density)),
                         N_I = length(unique(df$Individual)),
                         ID = as.integer(as.factor(df$Individual)))

# Write the stan model
write(
  temp <- "data{
  int<lower=1> N_z;       // Total number of observations of phenotype (Z)
  int<lower=1> N_I;       // Total number of unique individuals (I)
  int<lower=1> ID[N_z];   //  Individual ID repeated obs
  vector[N_z] z;          // Phenotypic observations
  vector[N_z] x;          // Environmental covariate
  }
  
  parameters {
  // Fixed effects
  real B_0; 
  real B_1; 
  
  // Random effects
  matrix[N_I,2] zI;           // Note the 2 columns for intercepts & slopes
  vector<lower=0>[2] sigma_I; // sd individual intercepts & slopes
  cholesky_factor_corr[2] L;  // Factor to estimate correlation int-slopes, equal to n-trait
  real<lower=0> sigma_e;
  }
  
  transformed parameters {
  vector[N_z] e_z; // Predicted values for phenotype
  
  // Note the use of the cholesky factor here when obtaining the unscaled values
  matrix[N_I,2] I = zI * diag_pre_multiply(sigma_I, L)'; 
  
  // model equation
  e_z = B_0 + (B_1 + I[ID,2]) .* x + I[ID, 1];
  // Note the use of the indices to assign the individual deviations to the matrix
  // 1: Individual intercepts
  // 2: Individual deviations from the population slope
  }
  
  model {
  // Likelihood function
  z ~ normal(e_z, sigma_e);
  
  // Priors
  B_0 ~ normal(0,1);
  B_1 ~ normal(0,1);
  to_vector(zI) ~ normal(0,1);
  sigma_I ~ exponential(3);
  L ~ lkj_corr_cholesky(3);   // Prior for correlations following a Cholesky distribution
  sigma_e ~ exponential(3);
  }
  
  generated quantities{
  // Variance
  real<lower=0> var_ID_int = sigma_I[1]^2;    // Note the use of matrix indexing
  real<lower=0> var_ID_slopes = sigma_I[2]^2;
  real<lower=0> var_res = sigma_e^2;
  real<lower=0> var_P = var_ID_int + var_ID_slopes + var_res;
  
  // Correlations & Covariances
  matrix[2,2] rho_I = L * L';             // Correlation matrix
  matrix[2,2] cov_I = diag_matrix(sigma_I)*rho_I*diag_matrix(sigma_I); // Covariance matrix
  }"
  , file = "Mod2.2.stan") 

# Fit the model
fit_mod2.2 <- stan("Mod2.2.stan", data = stan_data_mod2.2, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F,
                   seed = 26032025)

pars_mod2.2 <- c("B_0", "B_1", 
                 "var_ID_int", "var_ID_slopes", "var_res", "var_P",
                 "rho_I", "cov_I")

# Look at the summary estimates
round(rstan::summary(fit_mod2.2, pars = pars_mod2.2)$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
# Model 2.3
# Assignment 2: Add a random intercept & slope for population
#===============================================================================

# Tip: see 2.1 and 2.2 

# Model 3 Mulivariate linear mixed model 
#===============================================================================
# Model 3.1: 
# Extend model 2.1 to a multivariate model
#===============================================================================

# Prepare the data for the stan model
stan_data_mod3.1 <- list(N_z = nrow(df),
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
  int<lower=1> ID[N_z];       //  Individual ID repeated obs
  array[N_z] vector[2] zs;    // phenotypic observations for phenotype 1 & 2
  vector[N_z] x;              // Environmental covariate
  }
  
  parameters {
  // Fixed effects
  vector[2] B_0;    // Note the use of vectors and their dimension here
  vector[2] B_1; 
  
  // Random effects
  matrix[N_I,2] zI; 
  vector<lower=0>[2] sigma_I; 
  cholesky_factor_corr[2] L;
  cholesky_factor_corr[2] LR;   // Cholesky factor for residual correlation
  vector<lower=0>[2] sd_R;      // Temporary vector of sd's 
  }
  
  transformed parameters {
  vector[N_z] e_z1; // Predicted values for phenotype1
  vector[N_z] e_z2; // Predicted values for phenotype2
  
  matrix[N_I,2] I  = zI * diag_pre_multiply(sigma_I, L)'; 
  
  // model equations
  e_z1 = B_0[1] + B_1[1] * x + I[ID, 1];
  e_z2 = B_0[2] + B_1[2] * x + I[ID, 2];
  // Note the 2 model equations and the use of vectors for the beta's and
  // the matrix for individual intercepts for z1 and z2
  }
  
  model {
  // Transform expected values to an array
  array[N_z] vector[2] mus;
  for (o in 1:N_z) 
  mus[o] = [e_z1[o], e_z2[o]]';
  
  // Derive the sd's and residual correlation
  matrix[2,2] L_sigma = diag_pre_multiply(sd_R, LR);
  
  // Cholesky multinormal likelihood function to estimate residual correlations
  zs ~ multi_normal_cholesky(mus, L_sigma);  
  
  // Priors
  B_0 ~ normal(0,1);
  B_1 ~ normal(0,1);
  to_vector(zI) ~ normal(0,1);
  sigma_I ~ exponential(3);
  L ~ lkj_corr_cholesky(3);
  LR ~ lkj_corr_cholesky(3);    // Prior residual correlation
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
  matrix[2,2] rho_I = L * L';
  matrix[2,2] cov_I = diag_matrix(sigma_I)*rho_I*diag_matrix(sigma_I); 
  matrix[2,2] rho_R = LR * LR';   // Residual correlation matrix
  }"
  , file = "Models/Mod3.1.stan") 

# Fit the model
fit_mod3.1 <- stan("Models/Mod3.1.stan", data = stan_data_mod3.1, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F,
                   seed = 26032025)

pars_mod3.1 <- c("B_0", "B_1", 
                 "var_ID_int1", "var_res1", "var_P1",
                 "var_ID_int2", "var_res2", "var_P2",
                 "rho_I", "cov_I")

# Look at the summary estimates
round(rstan::summary(fit_mod3.1, pars = pars_mod3.1)$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
# Model 3.2
# Assignment 3: Add random slopes for individuals and correlate them across traits
#===============================================================================

# Tip: look at how the random intercepts for individuals are estimated and allowed to correlate. 
# Do the same for the individual random slopes

# Model 4 Error-in-variable (EIV) & selection gradient analysis 
#===============================================================================
# Model 4.1: 
# Add an EIV selection gradient to model 2.1 
#===============================================================================

# Prepare the data for the stan model
stan_data_mod4.1 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
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
  int<lower=1> ID[N_z];     //  Individual ID repeated obs
  int<lower=1> ID_w[N_w];   //  Individual ID repeated obs
  vector[N_z] z;            // Phenotypic observations
  vector[N_w] w;            // Fitness 
  }
  
  parameters {
  // Fixed effects
  real B_0; 
  real Bw_0;      // Fitness population intercept
  vector[2] Bw;   // Selection gradient slopes
  
  // Random effects
  vector[N_I] zI; 
  real<lower=0> sigma_I; 
  real<lower=0> sigma_e;  
  real<lower=0> sigma_ew;   // Standard deviation of fitness model likelihood
  }
  
  transformed parameters {
  vector[N_z] e_z; 
  vector[N_w] e_w;    // Predicted values for fitness
  
  vector[N_I] I  =  zI * sigma_I; 
  
  // model equation 1 phenotypic model
  e_z = B_0 + I[ID];
  // model equation 2 selection gradient
  e_w = Bw_0 + Bw[1] * I[ID_w] + Bw[2]*(I[ID_w].*I[ID_w]);
  // Note the use of the previously estimated intercepts I[ID_w]
  // Note the use of (I[ID_w].*I[ID_w]) to estimate a quadratic slope
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
  , file = "Models/Mod4.1.stan") 

# Fit the model
fit_mod4.1 <- stan("Models/Mod4.1.stan", data = stan_data_mod4.1, 
                   chains = 4, iter = 3000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F,
                   seed = 26032025)

pars_mod4.1 <- c("B_0", "Bw_0", "Bw",
                 "var_ID", "var_res", "var_P", "var_res_w")

# Look at the summary estimates
round(rstan::summary(fit_mod4.1, pars = pars_mod4.1)$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
# Model 4.2
# Change model 4.1 for repeated fitness measure and add social selection gradient (partner phenotype effect)
#===============================================================================

# Prepare the data for the stan model
stan_data_mod4.2 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         N_I = length(unique(df$Individual)),
                         ID = as.integer(as.factor(df$Individual)),
                         partnerID = as.integer(as.factor(df$Partner)), 
                         w = as.vector(scale(df$Feeding_rate)))

# Write the stan model
write(
  temp <- "data{
  int<lower=1> N_z;               // Total number of observations of phenotype (Z)
  int<lower=1> N_I;               // Total number of unique individuals (I)
  int<lower=1> ID[N_z];           //  Individual ID repeated obs
  int<lower=1> partnerID[N_z];    //  Social partner ID repeated obs
  vector[N_z] z;                  // Phenotypic observations
  vector[N_z] w;                  // Fitness 
  }
  
  parameters {
  // Fixed effects
  real B_0; 
  real Bw_0; 
  vector[4] Bw;             // Note number of slopes
  // Random effects
  vector[N_I] zI; 
  real<lower=0> sigma_I; 
  vector[N_I] zWI;          // Fitness intercepts for each individual
  real<lower=0> sigma_WI;   // Fitness sd individual intercepts
  real<lower=0> sigma_e;
  real<lower=0> sigma_ew;
  }
  
  transformed parameters {
  vector[N_z] e_z; 
  vector[N_z] e_w; 
  
  vector[N_I] I  =  zI * sigma_I; 
  vector[N_I] WI  =  zWI * sigma_WI; // Get the unscaled values for random effect
  
  // model equation 1 phenotypic model
  e_z = B_0 + I[ID];
  // model equation 2 selection gradient
  e_w = Bw_0 + Bw[1] * I[ID] +
               Bw[2] * I[partnerID] +
               Bw[3] * (I[ID].*I[ID]) + 
               Bw[4] * (I[ID].*I[partnerID]) +
               WI[ID];
  // Note the use of the partner's trait values based on the partner ID indices
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
  , file = "Models/Mod4.2.stan") 

# Fit the model
fit_mod4.2 <- stan("Models/Mod4.2.stan", data = stan_data_mod4.2, 
                   chains = 4, iter = 4000,  
                   warmup = 1500, thin = 1, 
                   cores = 4,
                   refresh = 250,
                   save_warmup = F,
                   seed = 26032025)

pars_mod4.2 <- c("B_0", "Bw_0", "Bw",
                 "var_ID", "var_res", "var_P", 
                 "var_ID_w", "var_res_w", "var_P_w")

# Look at the summary estimates
round(rstan::summary(fit_mod4.2, pars = pars_mod4.2)$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
# Model 4.3
# Assignment 4: Add individual random slopes to model 4.1 (not 4.2!)
# Setup up direct, quadratic and correlation selection gradients for individual intercepts and slopes 
#===============================================================================

# Tip: you got this! Look closely how to model quadratic and correlational selection. What does it mean?

#===============================================================================
# Bonus (assignments): working with shiny stan & posterior distribution 
#===============================================================================

# Bonus assignment 1
# Do some posterior checks with any of the models
launch_shinystan(fit_mod4.3)

# Bonus assignment 2
# Extract the posterior for further calculations
library(tidyverse)

# Extract posteriors
posterior_mod2.3 <- as.data.frame(rstan::extract(fit_mod2.3, pars = c("B_0", "B_1", 
                                                                      "var_ID_int", "var_ID_slopes",
                                                                      "var_pop_int", "var_pop_slopes",
                                                                      "var_res", "var_P")))

# Back transform variables
posterior_mod2.3$B_0 <- (posterior_mod2.3$B_0*sd(df$Exploration)) + mean(df$Exploration)
posterior_mod2.3$B_1 <- posterior_mod2.3$B_1*(sd(df$Exploration)/sd(df$Density))

posterior_mod2.3 <- posterior_mod2.3 %>%
  mutate(across(contains("var_"), ~ . * var(df$Exploration)))

# Derive repeatabilities 
posterior_mod2.3 <- posterior_mod2.3 %>%
  mutate(across(contains("var_"), 
                ~ . /var_P, 
                .names = "Rep_{.col}"))

# Derive quantiles of posterior per variable 
Sum_mod2.3 <- as.data.frame(t(round(apply(posterior_mod2.3[1:ncol(posterior_mod2.3)],2, quantile,probs = c(0.025, 0.5, 0.975)),3)))
Sum_mod2.3 <- Sum_mod2.3 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Sum_mod2.3 <- Sum_mod2.3 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Sum_mod2.3 <- Sum_mod2.3 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Sum_mod2.3 <- Sum_mod2.3 %>% relocate("Estimate")

# Make a figure with the BLUPS of a model
library(ggplot2)

BLUPS_I <- as.data.frame(round(summary(fit_mod4.2, pars = "I")$summary[,c(1,4,6, 8, 9,10)],3))
BLUPS_WI <- as.data.frame(round(summary(fit_mod4.2, pars = "WI")$summary[,c(1,4,6, 8, 9,10)],3))
W_slope <- as.data.frame(round(summary(fit_mod4.2, pars = "Bw")$summary[,c(1,4,6, 8, 9,10)],3))

ggplot() +
  geom_point(aes(x = (BLUPS_I$"50%"), 
                 y = (BLUPS_WI$"50%")), alpha = 0.5) +
  geom_abline(intercept = 0, slope = W_slope[1,"50%"], linewidth = 1.25, color = "red") +
  geom_errorbar(aes(
    x = (BLUPS_I$"50%"),
    ymin = (BLUPS_WI$"2.5%"),
    ymax = (BLUPS_WI$"97.5%")), alpha = 0.1) +
  geom_errorbarh(aes(
    y = (BLUPS_WI$"50%"),
    xmin = (BLUPS_I$"2.5%"),
    xmax = (BLUPS_I$"97.5%")), alpha = 0.1) + 
  labs(y = "Fitness (w)", x = "Average explorarion behaviour") +
  scale_x_continuous(limits = c(-1.6,1.6), breaks = seq(-1.5,1.5,0.5)) +
  scale_y_continuous(limits = c(-3.2,3.2), breaks = seq(-3,3,1)) +
  theme_classic(base_size = 17) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

#### Congrats you're a stan expert now!!! 8-)

#===============================================================================
