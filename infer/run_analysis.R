# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : stan_testing.r
# Programmer Name    : Luis Campos
#                     soyluiscampos@gmail.com
#
# Purpose            : An example data anlaysis using DLM code
#
# Input              : None
# Output             : None
# 
# References         : Measuring Effects of Medication Adherence on Time-Varying
# 					   Health Outcomes using Bayesian Dynamic Linear Models
# 					   Luis F. Campos, Mark E. Glickman, Kristen B. Hunter
# 					   https://arxiv.org/abs/1811.11072
#
# Platform           : R
# Version            : v3.3.0
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


library('ggplot2')
library('rstan')
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores() - 1)

setwd('/Users/khunter/Dropbox/Medication-Adherence')
source('state_space_functions.R')
source('stan_functions.R')



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Set parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
params = list(
  rho = c(0.8, 0.9), phi = c(0.5,0.7), sig = c(2.75, 2.2),
  sig.nu = c(1, 0.8),  sig.0 = c(1.5, 1.3), beta = c(120, -1.2, 80, -1.2), cor = 0.55
)
Sigma.BP = matrix(
  c(params[['sig']][1]^2, params[['cor']]*params[['sig']][1]*params[['sig']][2],
    params[['cor']]*params[['sig']][1]*params[['sig']][2], params[['sig']][2]^2), nrow = 2
 )
Sigma.nu = matrix(c(params[['sig.nu']][1]^2, 0, 0, params[['sig.nu']][2]^2), nrow = 2)
Sigma.0 = matrix(c(params[['sig.0']][1]^2, 0, 0, params[['sig.0']][2]^2), nrow = 2)
beta = matrix(params[['beta']], nrow = 2)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Data Parameters
# n: the number of patients
# T: the length of the follow-up period
# I.y: the days on which outcome is to be observed (defaults to 1:T)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
n = 10
T = 20
I.y = seq(1, 20, 5)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Simulate Data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

dat = replicate(n, {
		simulate.SS_2D(rho = params[[1]], phi = params[[2]], 
			Sigma.nu = Sigma.nu, Sigma.BP = Sigma.BP, 
			Sigma.0 = Sigma.0, T = T, beta = beta, I.y = I.y)
}, simplify = FALSE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# The way this is set up, all patients have the same follow-up period (T)
#  and outcome observation days I.y/
# We can randomize this as
# n.obs: number of observations per patient
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
n.obs = 4
dat = replicate(n, {
		simulate.SS_2D(rho = params[[1]], phi = params[[2]], 
			Sigma.nu = Sigma.nu, Sigma.BP = Sigma.BP, 
			Sigma.0 = Sigma.0, T = T, beta = beta, I.y = sort(sample(1:T, n.obs)))
}, simplify = FALSE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# we can similarly randomize the follow-up period T or number of 
#  observations n.obs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

dat = replicate(n, {
		simulate.SS_2D(rho = params[[1]], phi = params[[2]], 
			Sigma.nu = Sigma.nu, Sigma.BP = Sigma.BP, 
			Sigma.0 = Sigma.0, T = sample(80:120, 1), beta = beta, I.y = sort(sample(1:T, rpois(1, 2) + 1)))
}, simplify = FALSE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# This reformats the data to a form needed for STAN
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
adherence_dat = arrange.SS_stan(dat, is.sim = TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Run analysis using STAN
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
fit <- stan(file = 'state_space_adherence.stan', data = adherence_dat, iter = 200, chains = 4)
ex <- extract(fit)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Save out data and STAN fit
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
dir.create('data')
save(dat, file = 'data/dat.Rda')
save(ex, file = 'data/stan_fit.Rda')
# save out names of covariate columns
covariate.cols = c('gender')
save(covariate.cols, file = 'data/covariate_cols.Rda')
