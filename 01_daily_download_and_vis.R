# "01_daily_download_and_vis.R'
# Script for daily target data download and visualization
# Corresponds to 2/20/2026 "Pulling and Visualizing Data" Project Milestone

#### part 1: load packages ####
install.packages("readr")
library(readr)

#### part 2: bring in data ####
disease_url = 'https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/targets/project_id=bu4cast/tropical-disease-targets.csv'

disease_targets = read.csv(disease_url)

#### part 3: plot time-series visualization of data ####


