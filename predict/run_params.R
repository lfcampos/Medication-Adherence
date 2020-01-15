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
run.params[['test.code']] = TRUE
# running in the cloud? if running locally, will not use all cores
run.params[['cloud']] = FALSE
########################################################################

########################################################################
# number of iterations etc
########################################################################
# number of imputations
run.params[['imputations']] = 20
# number of posterior samples from training model
# corresponds to "h" in manuscript
run.params[['n.post.samp']] = 100
# number of patients to predict adherence for
run.params[['test.size']] = 2
# number of particles
run.params[['nparticles']] = 128
# number of particle filter iterations
# corresponds to "g" in manuscript
run.params[['npf']] = 1000
# particle filter burnin
run.params[['burnin']] = run.params[['npf']]/5
########################################################################

if(run.params[['test.code']])
{
  run.params[['imputations']] = 1
  run.params[['n.post.samp']] = 3
  run.params[['nparticles']] = 4
  run.params[['npf']] = 10
  run.params[['test.size']] = 2
  run.params[['burnin']] = run.params[['npf']]/5
}