setwd('/Users/khunter/Dropbox/Medication-Adherence')
source('simulate/sim_data_functions.R')
source('infer/state_space_functions.R')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Data Parameters
# n: the number of patients
# T: the length of the follow-up period
# n.obs: number of observations per patient
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
n = 500
T = 90
n.obs = 30

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
# Simulate Data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# dat = replicate(n, {
#   simulate.SS_2D(rho = params[[1]], phi = params[[2]],
#                  Sigma.nu = Sigma.nu, Sigma.BP = Sigma.BP,
#                  Sigma.0 = Sigma.0, T = T, beta = beta, I.y = I.y)
# }, simplify = FALSE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# The way this is set up, all patients have the same follow-up period (T)
#  and outcome observation days I.y/
# We can randomize this as

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

dat = replicate(n, {
  simulate.SS_2D(rho = params[[1]], phi = params[[2]],
                 Sigma.nu = Sigma.nu, Sigma.BP = Sigma.BP,
                 Sigma.0 = Sigma.0, T = T, beta = beta, I.y = sort(sample(1:T, n.obs)))
}, simplify = FALSE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# we can similarly randomize the follow-up period T or number of
#  observations n.obs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# dat = replicate(n, {
#   simulate.SS_2D(rho = params[[1]], phi = params[[2]],
#                  Sigma.nu = Sigma.nu, Sigma.BP = Sigma.BP,
#                  Sigma.0 = Sigma.0, T = sample(80:120, 1), beta = beta, I.y = sort(sample(1:T, rpois(1, 2) + 1)))
# }, simplify = FALSE)


if(!dir.exists('data')) { dir.create('data') }
saveRDS(dat, file = 'data/dat.RDS')
