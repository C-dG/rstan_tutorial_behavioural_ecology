#===============================================================================
# Model 1.2
# Assignment 1: Add the fixed effect Sex
#===============================================================================

# Prepare the data for the stan model
stan_data_mod1.2 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         x = as.vector(scale(df$Density)),
                         sex = ifelse(df$Sex == "F", -0.5, 0.5))

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
# Model 1.3 
# Add an interaction effect
#===============================================================================

# Prepare the data for the stan model
stan_data_mod1.3 <- list(N_z = nrow(df),
                         z = as.vector(scale(df$Exploration)),
                         x = as.vector(scale(df$Density)),
                         sex = ifelse(df$Sex == "F", -0.5, 0.5))