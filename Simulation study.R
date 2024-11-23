# --------------------------------------------------------------------------------------
# Data Generation using the skewlmm Package
# This section generates simulated a longitudinal data based on a random slope model using Contaminated-normal distribution
# --------------------------------------------------------------------------------------

# Number of iterations, subjects (IDs), and time points
N_iter <- 500  # Number of iterations
N_id <- 100   # Number of subjects (IDs)
N_time <- 5   # Number of time points per subject

# Create a data frame with repeated measures
data <- expand.grid(iter = 1:N_iter, id = 1:N_id, time = 1:N_time)

# Parameters for simulations
Niter <- length(unique(data$iter)) # Number of iterations
N <- length(unique(data$id))       # Number of subjects
Ti <- length(unique(data$time))    # Number of time points per subject
Ntot <- nrow(data)                 # Total number of observations
iter <- data$iter
id <- data$id
time <- data$time

# Initialize fixed and interaction effects
set.seed(123)
t0 <- c(1, 2, 3, 4, 5)  # Time points
x1 <- rbinom(N, 1, 0.5) # Binary predictor variable

# Arrays to store fixed effects and their interactions
x11 <- array(0, dim = c(Niter, N, Ti))
x22 <- array(0, dim = c(Niter, N, Ti))
x33 <- array(0, dim = c(Niter, N, Ti))

# Populate the arrays
for (it in 1:Niter) {
  for (i in 1:N) {
    for (t in 1:Ti) {
      x11[it, i, t] <- x1[i]
      x22[it, i, t] <- t0[t]
      x33[it, i, t] <- x11[it, i, t] * x22[it, i, t]
    }
  }
}

# Generate outcome variable 
y1 <- array(0, dim = c(Niter, N, Ti))
X1 <- array(0, dim = c(Niter, N, Ti))
X2 <- array(0, dim = c(Niter, N, Ti))

# Covariance matrix for random effects
D0 <- matrix(c(2, 0.7, 0.7, 1), nrow = 2, ncol = 2)

# Consider an AR(1) structure with phi = 0.5 for within-subject error
set.seed(123)
for (it in 1:Niter) {
  for (i in 1:N) {
    y0 <- rsmsn.lmm(1:Ti, x1 = cbind(1, x11[it, i, ], x22[it, i, ], x33[it, i, ]),
                    z1 = cbind(1, x22[it, i, ]),
                    sigma2 = 0.25, D1 = D0,
                    beta = c(0.5, 1.5, 2, 0.7), lambda = c(0, 0),
                    depStruct = "ARp", phi = 0.5, distr = "scn", nu = c(0.3, 0.3))
    
    # Store results
    y1[it, i, ] <- y0$y
    X1[it, i, ] <- y0$x.1
    X2[it, i, ] <- y0$x.2
  }
}

# Add generated data to the data frame
data$y <- as.vector(y1)
data$x1 <- as.vector(X1)
data$x2 <- as.vector(X2)

# Sort data by iteration, ID, and time
data <- data[order(data$iter, data$id, data$time), ]

# Save the generated data to a CSV file
write.csv(data, file = "generated_data.csv", row.names = FALSE)

# --------------------------------------------------------------------------------------
# Model Fitting using Marginal Huber Model with Stan
# --------------------------------------------------------------------------------------

# Load required libraries
library(rstan)
library(dplyr)
library(tidyr)

# Load generated data
data <- read.csv("generated_data.csv")


# Recompute parameters for simulations
Niter <- length(unique(data$iter)) # Number of iterations
N <- length(unique(data$id))       # Number of subjects
Ti <- length(unique(data$time))    # Number of time points per subject
iter <- data$iter
id <- data$id
time <- data$time

# Stan model code 
model <- 	
  paste("
  data {
int  N; 
int  Ti;
matrix[N,Ti] y;
matrix[N,Ti] x1;
matrix[N,Ti] x2;
matrix[N,Ti] x1x2;
}


transformed data {
real pi=pi();
}

parameters {
  vector[4] beta; 
 vector<lower=0>[N] c;
real<lower=0> sigma;
real<lower=0> omega1;
real<lower=0> omega0;
real omega01;
}

transformed parameters {
matrix[N,Ti] mu;
matrix[N,Ti] u;
matrix[Ti,Ti] S;
matrix[Ti,Ti] R;
matrix[Ti,Ti] Rinv;
matrix[N,Ti] z;
matrix[N,Ti] q1;
matrix[N,Ti] h1;
matrix[N,Ti] h2;
matrix[N,Ti] h3;
real de;
vector[N] A;

for(i in 1:N) {
A[i] = 2*sqrt(2*pi)*(Phi(c[i]) + 1/sqrt(2*pi)*exp(-0.5*c[i]*c[i]) /c[i] - 0.5);}

mu= beta[1] + beta[2]*x1+beta[3]*x2+beta[4]*x1x2;
u=y-mu;

# covariance matrix
for (i in 1:Ti) {for (j in 1:Ti) {
    if (i == j) {
      S[i, j] = pow(omega0,2) + pow(omega1,2) * i*i + 2 *i * omega01 +sigma*sigma;
    } else {
      S[i, j] =
     pow(omega0,2) + omega01 * (i+j) + pow(omega1,2) * i*j; }}}

R=cholesky_decompose(S);
Rinv=inverse(R);
z=u*Rinv;

for(i in 1:N){for(t in 1:Ti){q1[i,t]= step(c[i]-fabs(z[i,t]));}}

for(i in 1:N){for(t in 1:Ti){h1[i,t]=q1[i,t]*pow(z[i,t],2);}}

for(i in 1:N){for(t in 1:Ti){h2[i,t]= c[i]*(1-q1[i,t])*fabs(z[i,t]);}}

for(i in 1:N){for(t in 1:Ti){h3[i,t]= pow(c[i],2)*(1-q1[i,t]);}}

de=determinant(R);
}


model {
vector[N] loglike;
for(i in 1:N) {
loglike[i]= -log(de)-Ti*log(A[i])-0.5*sum(h1[i,1:Ti])-sum(h2[i,1:Ti])+0.5*sum(h3[i,1:Ti]);}
target += loglike;


  for(i in 1:4){ beta[i] ~ normal(0, 100);} 
  sigma ~ normal(0,1)T[0,];
  omega0 ~ normal(0,1)T[0,];
 omega1 ~ normal(0,1)T[0,];


 # omega01~ normal(0,1);
 for(t in 1:Ti){c[t]~normal(0,100)T[0,];}
}
")
writeLines(model,"Huber.stan")

# Prepare data for Stan for each iteration
data_list <- list()
for (iter_val in unique(data$iter)) {
  iter_data <- subset(data, iter == iter_val)
  N <- length(unique(iter_data$id))
  Ti <- length(unique(iter_data$time))
  y_matrix <- matrix(iter_data$y, nrow = N, byrow = TRUE)
  x1_matrix <- matrix(iter_data$x1, nrow = N, byrow = TRUE)
  x2_matrix <- matrix(iter_data$x2, nrow = N, byrow = TRUE)
  x1x2_matrix <- matrix(iter_data$x1 * iter_data$x2, nrow = N, byrow = TRUE)
  
  data_list[[iter_val]] <- list(
    N = N,
    Ti = Ti,
    y = y_matrix,
    x1 = x1_matrix,
    x2 = x2_matrix,
    x1x2 = x1x2_matrix
  )
}

# Compile the Stan model
stan_model <- stan_model("Huber.stan")

# Store parameter estimates for all iterations
all_results <- list()

# Run the model for each iteration
for (iter_val in seq_along(data_list)) {
  # Fit the model
  stan_fit <- sampling(
    stan_model,
    data = data_list[[iter_val]],
    iter =5000, chains = 4, seed = 123  
  )
  
  param_names <- c("beta", "sigma", "omega0", "omega1", "omega01", "c")
  fit_summary <- summary(stan_fit, pars = param_names, probs = c(0.1, 0.5, 0.9))$summary
  
  # Store estimates for all iterations
  all_results[[iter_val]] <- fit_summary
}


param_names <- rownames(all_results[[1]])

# Initialize vectors to store results
mean <- numeric()
sd <- numeric()
lower <- numeric()
median <- numeric()
upper <- numeric()
max_Rhat <- numeric()
min_n_eff<- numeric()

# Iterate through each parameter name
for (param in param_names) {
  # Extract the mean values for the current parameter across all iterations
  param_means <- sapply(all_results, function(iter_result) {
    iter_result[param, "mean"]
  })
  
  # Extract the SD values for the current parameter across all iterations
  param_sds <- sapply(all_results, function(iter_result) {
    iter_result[param, "sd"]
  })
  
  # Extract the Rhat values for the current parameter across all iterations
  param_rhats <- sapply(all_results, function(iter_result) {
    iter_result[param, "Rhat"]
  })
  
  # Extract the n_eff values for the current parameter across all iterations
  param_ess <- sapply(all_results, function(iter_result) {
    iter_result[param, "n_eff"]
  })
  
  # Compute the mean estimate across all iterations
  mean[param] <- mean(param_means)
  
  # Compute the mean SD across all iterations
  sd[param] <- mean(param_sds, na.rm = TRUE)
  
  # Calculate the credible intervals (10%, 50%, 90%) across iterations
  lower[param] <- quantile(param_means, 0.10)
  median[param] <- quantile(param_means, 0.50)
  upper[param] <- quantile(param_means, 0.90)
  
  # Compute the maximum Rhat across iterations
  max_Rhat[param] <- max(param_rhats, na.rm = TRUE)
  
  # Compute the minimum n_eff across iterations
  min_n_eff[param] <- min(param_ess, na.rm = TRUE)
}

# Convert the results into a data frame for easy viewing
summary_df <- data.frame(
  Parameter = names(mean),
  Mean = mean,
  SD = sd,
  Lower = lower,
  Median = median,
  Upper = upper,
  Max_Rhat = max_Rhat,
  Min_n_eff = min_n_eff
)

