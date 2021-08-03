# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : stan_testing.r
# Programmer Name    : Luis Campos
#                     soyluiscampos@gmail.com
#
# Purpose            : An example data analysis using DLM code
#
# Input              : None
# Output             : None
#
# References         : Measuring Effects of Medication Adherence on Time-Varying
# 					           Health Outcomes using Bayesian Dynamic Linear Models
# 					           Luis F. Campos, Mark E. Glickman, Kristen B. Hunter
# 					           https://arxiv.org/abs/1811.11072
#
# Platform           : R
# Version            : v3.3.0
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

library(ggplot2)
library(rstan)
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores() - 1)

setwd('/Users/khunter/Dropbox/Medication-Adherence')
source('infer/state_space_functions.R')
source('infer/stan_functions.R')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Read in and reformat the data to a form needed for STAN
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
dat <- readRDS(file = 'data/dat.RDS')
adherence_dat = arrange.SS_stan(dat, is.sim = TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Run analysis using STAN
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
fit <- stan(file = 'infer/state_space_adherence.stan', data = adherence_dat, iter = 100, chains = 2)
ex <- extract(fit)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Save out data and STAN fit
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

saveRDS(ex, file = 'data/stan_fit.RDS')
covariate.cols = c('gender')
saveRDS(covariate.cols, file = 'data/covariate_cols.RDS')
