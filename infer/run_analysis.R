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

n.train = 70
n.test = 30
n.iter = 8000
n.chains = 4

# setwd('/Users/khunter/Dropbox/Medication-Adherence')
setwd('/n/home01/khunter33/Medication-Adherence/')
source('infer/state_space_functions.R')
source('infer/stan_functions.R')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Read in and reformat the data to a form needed for STAN
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
dat <- readRDS(file = 'data/dat.RDS')

# sample down
set.seed(100)
ids <- 1:length(dat)
train.ids <- sample(ids, n.train)
test.ids <- sample(ids[!(ids %in% train.ids)], n.test)
all.ids <- c(train.ids, test.ids)
dat <- dat[all.ids]

adherence_dat <- arrange.SS_stan(dat, is.sim = TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Run analysis using STAN
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

init.params <- list(
  'chain1' = list(
    'rho' = c(0.01, 0.01),
    'phi' = c(0, 0),
    'sigma' = c(0.2, 0.2),
    'cor' = 0,
    'sigma_nu' = c(0.2, 0.2),
    'sigma_0' = c(0.2, 0.2),
    'beta' = matrix(c(120, 80, 0, 0), nrow = 2, ncol = 2, byrow = TRUE)
  ),
  'chain2' = list(
    'rho' = c(0.5, 0.5),
    'phi' = c(2, 2),
    'sigma' = c(0.5, 0.5),
    'cor' = -0.5,
    'sigma_nu' = c(0.5, 0.5),
    'sigma_0' = c(0.5, 0.5),
    'beta' = matrix(c(120, 80, 1, 1), nrow = 2, ncol = 2, byrow = TRUE)
  ),
  'chain3' = list(
    'rho' = c(0.8, 0.8),
    'phi' = c(-2, -2),
    'sigma' = c(1, 1),
    'cor' = 0.5,
    'sigma_nu' = c(1, 1),
    'sigma_0' = c(1, 1),
    'beta' = matrix(c(120, 80, -1, -1), nrow = 2, ncol = 2, byrow = TRUE)
  ),
  'chain4' = list(
    'rho' = c(0.2, 0.2),
    'phi' = c(2, -2),
    'sigma' = c(2, 2),
    'cor' = 0,
    'sigma_nu' = c(0.1, 0.1),
    'sigma_0' = c(2, 2),
    'beta' = matrix(c(120, 80, 0.2, -3), nrow = 2, ncol = 2, byrow = TRUE)
  )
)


fit <- stan(
  file = 'infer/state_space_adherence.stan',
  data = adherence_dat,
  init = init.params,
  iter = n.iter,
  chains = n.chains
)

ex <- extract(fit)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Save out data and STAN fit
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

saveRDS(ex, file = 'data/stan_fit.RDS')
covariate.cols <- c('gender')
saveRDS(covariate.cols, file = 'data/covariate_cols.RDS')
