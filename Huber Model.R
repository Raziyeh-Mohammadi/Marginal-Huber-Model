# Load the required library
library(lme4)
library(dplyr)
library(tidyr)
library(rstan)

# Set working directory 
setwd("Path/To/Your/Project")

# Load the sleepstudy data set
sleep_df <- as.data.frame(sleepstudy)

# Reshape data into wide matrices for Stan modeling
reaction_matrix <- sleep_df %>%
  select(Subject, Days, Reaction) %>%
  pivot_wider(names_from = Days, values_from = Reaction) %>%
  select(-Subject) %>%
  as.matrix()

days_matrix <- sleep_df %>%
  select(Subject, Days) %>%
  pivot_wider(names_from = Days, values_from = Days) %>%
  select(-Subject) %>%
  as.matrix()

# Define constants
N <- nrow(reaction_matrix)  # Number of subjects
Ti <- ncol(reaction_matrix) # Number of time points

#Model
model <- 	
  paste("
  data {
int  N; 
int  Ti;
matrix[N,Ti] reaction;
matrix[N,Ti] days;
}


transformed data {
real pi=pi();
}

parameters {
  vector[2] beta; 
 real<lower=-1, upper=1> rho;
 vector<lower=0>[N] c;
real<lower=0> sigma;
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

mu= beta[1] + beta[2]*days;
u=reaction-mu;

for(i in 1:Ti){for(j in 1:Ti){S[i,j]=sigma*sigma*pow(rho,fabs(i-j));}}
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
 // Likelihood
   vector[N] loglike;
  for(i in 1:N) {
  loglike[i]= -log(de)-Ti*log(A[i])-0.5*sum(h1[i,1:Ti])-sum(h2[i,1:Ti])+0.5*sum(h3[i,1:Ti]);}
   target += loglike;

\\Priors
  for(i in 1:2){ beta[i] ~ normal(0, 100);} 
  sigma ~ normal(0,1)T[0,];
  rho ~ uniform(-1,1);
  for(t in 1:Ti){c[t]~normal(0,100)T[0,];}
  }
")
writeLines(model, "sleep.stan")

# Create the Stan data list
sleep.data<- list( reaction = reaction_matrix,days = days_matrix,N=N,Ti=Ti)

#Stan
fit1 <- stan(file = "sleep.stan",data = sleep.data,iter = 10000)
print(fit1, pars=c("beta", "sigma","rho","c"))
traceplot(fit1, par=c("beta", "sigma","rho","c"))

#MLE
sm <- stan_model(model_code = model)
mle<-optimizing(sm, data =sleep.data,algorithm = "BFGS", hessian = TRUE,init = 'random')
print(mle, digits=2)







