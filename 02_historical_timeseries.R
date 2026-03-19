# "02_historical_timeseries.R
# Script for fitting historical data using Bayesian state-space model
# Corresponds to 2/20/2026 "Pulling and Visualizing Data" Project Milestone



#### load packages and data ####
library(rjags)
library(daymetr)
library(ecoforecastR)

disease_url = 'https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/project_id=bu4cast/targets/tropical-disease-targets.csv'

disease_targets = read.csv(disease_url)


#### Define JAGS code ####

RandomWalk = "
model{
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs)
  }
  
  #### Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
}


#### Define ####

