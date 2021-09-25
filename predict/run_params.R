# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : run_params.R
# Programmer Name    : Kristen Hunter
#                      kristenbhunter@gmail.com
#
# Last Updated       : Jan 2020
#
# Purpose            : Allows you to easily change run parameters, such as
#                      number of patients in test set, number of particles, etc.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

run.params = NULL

########################################################################
# type of run
########################################################################
# for a quick code test, automatically runs a quick version, see below
run.params[['test.code']] = FALSE
# running in the cloud? if running locally, will not use all cores
run.params[['cloud']] = FALSE
########################################################################

########################################################################
# number of iterations etc
########################################################################
# number of posterior samples from training model
# corresponds to "f" in manuscript
run.params[['npostsamp']] = 1
# number of patients to predict adherence for. We have
# 503 total, and 400 are fixed in the training set
run.params[['train.size']] = 30
run.params[['test.size']] = 10
# number of particles
# 'P' in manuscript
run.params[['nparticles']] = 8
# number of particle filter iterations
# corresponds to "e" in manuscript
run.params[['npf']] = 1
# particle filter burnin
run.params[['burnin']] = round(run.params[['npf']]/5)
# MCMC length
run.params[['mcmc.chains']] = 2
run.params[['theta.h.mcmc.length']] = 500
run.params[['theta.a.mcmc.length']] = 500
########################################################################

if(run.params[['test.code']])
{
  run.params[['test.size']] = 2
  run.params[['train.size']] = 3
  run.params[['npostsamp']] = 2
  run.params[['nparticles']] = 4
  run.params[['npf']] = 3
  run.params[['burnin']] = round(run.params[['npf']]/5)
  run.params[['mcmc.chains']] = 2
  run.params[['theta.h.mcmc.length']] = 120
  run.params[['theta.a.mcmc.length']] = 120
}
