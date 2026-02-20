# "01_daily_download_and_vis.R'
# Script for daily target data download and visualization
# Corresponds to 2/20/2026 "Pulling and Visualizing Data" Project Milestone

#### part 1: load packages ####
# install.packages("readr")
# install.packages("ggplot2")
library(readr)
library(ggplot2)

#### part 2: bring in data ####
disease_url = 'https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/targets/project_id=bu4cast/tropical-disease-targets.csv'

disease_targets = read.csv(disease_url)

#### part 3: plot time-series visualization of data ####

# split observations by time to create timeseries for each site
obs_by_time <- split(disease_targets$observation, disease_targets$datetime) 

# quantiles per time point to get a sense of the distribution of observations across sites at each time point
n.stats <- sapply(obs_by_time, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

time <- as.Date(names(obs_by_time)) # convert from character to Date for plotting
ylo  <- n.stats[1, ]
ymed <- n.stats[2, ]
yhi  <- n.stats[3, ]

# function to draw a confidence interval ribbon between ylo and yhi across time points in x
# borrowed from EF_Activities/Exercise_02_Logistic
ciEnvelope <- function(x, ylo, yhi, col = rgb(0,0,0,0.2), ...) {
  polygon(c(x, rev(x), x[1]),
          c(ylo, rev(yhi), ylo[1]),
          border = NA, col = col, ...)
}

par(mfrow = c(2, 1),mar = c(3, 4, 3, 1))

plot(time, ymed, type = "l",
     ylim = range(c(ylo, yhi), na.rm = TRUE),
     ylab = "VL Cases",
     main = "Median trajectory across sites with 95% CI")

ciEnvelope(time, ylo, yhi) # add confidence interval
lines(time, ymed, lwd = 2) # redraw median on top

plot(time, ymed, type = "l",
     ylim = range(c(ylo, yhi, disease_targets$observation), na.rm = TRUE),
     ylab = "VL Cases",
     main = "Median trajectory with 95% CI and site observations")

ciEnvelope(time, ylo, yhi) # add confidence interval
lines(time, ymed, lwd = 2) # redraw median on top

# add points to see observations outside the quantiles
points(as.Date(disease_targets$datetime), disease_targets$observation, pch = 16, cex = 0.2)
