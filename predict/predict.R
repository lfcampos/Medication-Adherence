# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : predict.r
# Programmer Name    : Kristen Hunter
#                      kristenbhunter@gmail.com
#
# Last Updated       : Jan 2021
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
# base.dir = '/Users/khunter/Dropbox/archive/Medication-Adherence/'
base.dir = '/Users/khunter/Dropbox/Medication-Adherence/'

# source necessary files
predict.dir = paste(base.dir, 'predict/', sep = '')
source(paste(predict.dir, 'setup.R', sep = ''))
run.dir = base.setup(base.dir)

predict.adherence = function()
{
  # setup the data and simulation
  setup.output = setup.all(base.dir, run.dir)

  # save out important parameters
  datasets = setup.output[['datasets']]

  if(run.params[['test.code']])
  {
    cl = NULL
  } else
  {
    # initiate cluster
    num.cores = ifelse(run.params[['cloud']], detectCores(), detectCores() - 1)
    outfile = paste(run.dir, 'out.txt', sep = '')
    cl = makeCluster(num.cores, outfile = outfile)
  }

  # calculate how many days to predict for (i.e. the maximum any patient has)
  run.params[['npredictdays']] = max(datasets$test$data.melt$mday)

  # these are posterior draws from training set for blood pressure model
  # provided by Luis
  theta.h = setup.output[['theta.h']]

  ################################################
  # posterior draws from training set for unconditional adherence model
  ################################################
  theta.a = get.theta.a.draws(
    theta.a.stan.dat = datasets[['train']][['data.melt']],
    covariate.cols = datasets[['train']][['covariate.cols']],
    run.params, base.dir
  )
  save(theta.a, file = paste(run.dir, 'theta_a.RData', sep = ''))
  theta = save.theta.draws(theta.a, theta.h, run.params)

  ################################################
  # predictions of adherence
  ################################################

  set.seed(05242021)

  predict.draws = draw.c.star(datasets = datasets[['test']], theta, run.params, base.dir, cl)

  if(!run.params[['test.code']])
  {
    stopCluster(cl)
  }

  return(predict.draws)
}

################################################
# generate predicted values
################################################
draws = predict.adherence()
# save out data
saveRDS(draws, file = paste(run.dir, 'draws.rds', sep = ''))
