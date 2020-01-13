# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : predict.r
# Programmer Name    : Kristen Hunter
#                      kristenbhunter@gmail.com
#
# Last Updated       : Jan 2020
#
# Purpose            : An example data anlaysis to predict adherence values using
#                      simulated data
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
  ################################################
  # setup
  ################################################
  
  # setup the data and simulation
  setup.output = data.setup(base.dir)
  
  # save out important parameters
  datasets = setup.output[['datasets']]
  run.params = setup.output[['run.params']]
  # posterior draws from training set for health outcomes model
  theta.h = setup.output[['theta.h']]
  
  ################################################
  # posterior draws from training set for unconditional adherence model
  ################################################
  
  theta.a.model = get.theta.a.draws.stan(datasets, run.params, predict.dir, run.dir)
  theta.a = theta.a.model[['theta.a']]
  save(theta.a, file = paste(run.dir, 'theta_a.RData', sep = ''))
  theta = save.theta(theta.a, theta.h, run.params)

  ################################################
  # predictions of adherence
  ################################################
  
  if(run.params[['test.code']])
  {
    cl = NULL
  } else
  {
    # initiate cluster (if running locally, don't use all available cores)
    num.cores = ifelse(run.params[['cloud']], detectCores(), detectCores() - 1)
    outfile = paste(run.dir, 'out.txt', sep = '')
    cl = makeCluster(num.cores, outfile = outfile)
  }
  
  # predict mean adherence values for each patient
  predict.draws = draw.ad.star(datasets, theta, run.params, cl)
  
  if(!run.params[['test.code']])
  {
    stopCluster(cl)
  }
  
  return(predict.draws)
}

################################################
# generate predicted values
################################################
impute.draws = predict.adherence()

# extract which imputation number if parallelizing across imputations
impute.index = as.numeric(as.character(Sys.getenv("INDEX_VAR")))
if(is.na(impute.index)){ impute.index = 1 }

# save out data
saveRDS(impute.draws, file = paste(run.dir, 'draws_', impute.index, '.rds', sep = ''))
