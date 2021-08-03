# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : predict_onestep.r
# Programmer Name    : Kristen Hunter
#                      kristenbhunter@gmail.com
#
# Last Updated       : July 2021
#
# Purpose            : An example data anlaysis to predict adherence values using
#                      simulated data. Does predictions in one step instead of 2.
#
# Input              : Takes in the base directory where github repo is located
# Output             : Saves out output and plots in output directory
#
# References         :
#
# Platform           : R
# Version            : v3.6.1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

########################################################
# CHANGE THIS to github repo directory
########################################################
base.dir = '/Users/khunter/Dropbox/Medication-Adherence/'

# source necessary files
predict.dir = paste(base.dir, 'predict/', sep = '')
source(paste(predict.dir, 'setup.R', sep = ''))
run.dir = base.setup(base.dir)

predict.adherence = function()
{
  # setup the data and parameters
  setup.output = setup.all(base.dir)
  datasets = setup.output[['datasets']]
  run.params = setup.output[['run.params']]

  # predict adherence
  predict.draws = draw.c.star.onestep(datasets, run.params, base.dir)

  return(predict.draws)
}

################################################
# run the code
################################################

draws = predict.adherence()
saveRDS(draws, file = paste0(run.dir, 'draws_onestep.rds'))
summary(warnings())

