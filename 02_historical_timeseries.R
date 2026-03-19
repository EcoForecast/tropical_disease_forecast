# "02_historical_timeseries.R
# Script for fitting historical data using Bayesian state-space model
# Corresponds to 2/20/2026 "Pulling and Visualizing Data" Project Milestone



#### load packages and data ####
library(rjags)
library(daymetr)
library(ecoforecastR)

disease_url = 'https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/project_id=bu4cast/targets/tropical-disease-targets.csv'

disease_targets = read.csv(disease_url)

sp <- subset(disease_targets, site_id >= 350000 & site_id < 360000)  # SP = 35xxxx IBGE prefix
sp_monthly <- aggregate(observation ~ datetime, data = sp, sum, na.rm = TRUE)

# Format input data
time = as.Date(sp_monthly$datetime)
y = sp_monthly$observation

plot(time, y, type='l', ylab="Cases", lwd=2)

#### Define JAGS code ####

RandomWalk = "
model{
  
 #### Data Model (Negative Binomial)
  for(t in 1:n){
    y[t] ~ dnegbin(p[t], r)
    p[t] <- r / (r + mu[t])
    log(mu[t]) <- x[t]
  }
  
  #### Process Model (same random walk)
  for(t in 2:n){
    x[t] ~ dnorm(x[t-1], tau_add)
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic, tau_ic)
  
  tau_add ~ dgamma(a_add, r_add)
  
  # Dispersion parameter (size)
  r ~ dgamma(0.001, 0.001)
}
"
# Notes
  # JAGS parameterizes NB with r (size/dispersion) and p (success probability)
  # Ensure mu[t] > 0
  # Convert mean -> probability

# Define data/priors as a list
data <- list(
  y = as.integer(y), # has to be integer for negative binary
  n = length(y),
  x_ic = log(mean(y) + 1),  # based on mean counts
  tau_ic = 1,
  a_add = 1,
  r_add = 1
)

# Define initial model state
nchain = 3
init <- list()

for(i in 1:nchain){
  y.samp = sample(y, length(y), replace=TRUE)
  
  init[[i]] <- list(
    tau_add = 1 / var(diff(log(y.samp + 1))),  
    r = 10   # starting guess for dispersion
  )
}

# Send info to JAGS
j.model <- jags.model(
  file = textConnection(RandomWalk),
  data = data,
  inits = init,
  n.chains = nchain
)

# MCMC chain samples
## burn-in
jags.out <- coda.samples(
  model = j.model,
  variable.names = c("tau_add","r"),
  n.iter = 2000
)

plot(jags.out)
dic.samples(j.model, 2000)

# Larger samples
jags.out <- coda.samples(
  model = j.model,
  variable.names = c("x","tau_add","r"),
  n.iter = 10000
)

# Visualize
time.rng = c(1,length(time))
out <- as.matrix(jags.out)

x.cols <- grep("^x", colnames(out))

ci <- apply(exp(out[,x.cols]), 2, quantile,
            c(0.025,0.5,0.975))

plot(time, ci[2,], type='n',
     ylim=range(y, na.rm=TRUE),
     ylab="Cases",
     xlim=time[time.rng])

if(diff(time.rng) < 100){ 
  axis.Date(1,
            at=seq(time[time.rng[1]],
                   time[time.rng[2]],
                   by='month'),
            format="%Y-%m")
}

ecoforecastR::ciEnvelope(
  time, ci[1,], ci[3,],
  col=ecoforecastR::col.alpha("lightBlue",0.75)
)

points(time, y, pch="+", cex=0.5)