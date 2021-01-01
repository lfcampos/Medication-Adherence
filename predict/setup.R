# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : setup.R
# Programmer Name    : Kristen Hunter
#                      kristenbhunter@gmail.com
# 
# Last Updated       : Jan 2021
#
# Purpose            : General setup, including installing packages, loading
#                      in data, formatting data, etc.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# install missing packages
list.of.packages = c(
  'abind', 'bayesm', 'devtools', 'ggplot2', 
  'knitr', 
  'parallel', 'plyr', 'Rcpp', 'reshape2', 
  'rstan', 'sqldf', 'wesanderson'
)
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]
if(length(new.packages)) install.packages(new.packages, repos = 'http://cran.r-project.org')

library(abind)
library(bayesm)
library(devtools)
library(ggplot2)
library(knitr)
library(parallel)
library(plyr)
library(Rcpp)
library(reshape2)
library(rstan)
library(sqldf)
library(wesanderson)

if(!('CoupledCPF' %in% installed.packages()[,'Package'])){
  install_github('pierrejacob/CoupledCPF')
}
library(CoupledCPF)

# sets up directories, cleans up data, etc.
base.setup = function(base.dir)
{
  # source all necessary code files
  predict.dir = paste(base.dir, 'predict/', sep = '')
  source(paste(predict.dir, 'run_params.R', sep = ''))
  source(paste(predict.dir, 'functions.R', sep = ''))
  source(paste(predict.dir, 'functions_pf.R', sep = ''))
  source(paste(predict.dir, 'plots.R', sep = ''))
  
  # create output directory if doesn't exist
  out.dir = paste(predict.dir, 'out/', sep = '')
  if(!dir.exists(out.dir)){ dir.create(run.dir) }
  
  # create a unique directory for each simulation run based on timestamp
  Sys.setenv(TZ = 'America/New_york')
  time.start = Sys.time()
  run = format(time.start, format = '%Y%m%d_%H')
  run.dir = paste(out.dir, 'run_', run, '/', sep = '')
  if(!dir.exists(run.dir)){ dir.create(run.dir) }
  
  # save out run parameters
  save(run.params, file = paste(run.dir, 'run_params_', run, '.Rdata', sep = ''))
  
  # also save out version of stan file for future reference
  file.copy(from = paste(predict.dir, 'ran_ef.stan', sep = ''),
            to = paste(run.dir, 'ran_ef.stan', sep = ''))
  
  return(run.dir)
}

# loads in and sets up data
data.setup = function(base.dir)
{
  data.dir = paste(base.dir, 'data/', sep = '')
  
  # load in data and posterior draws of blood pressure model
  # theta h = blood pressure model, provided by previous analysis
  load(paste(data.dir, 'dat.Rda', sep = ''))
  load(paste(data.dir, 'covariate_cols.Rda', sep = ''))
  load(paste(data.dir, 'stan_fit.Rda', sep = ''))
  theta.h = ex
  dimnames(theta.h$beta) = list('iterations' = seq(1, dim(theta.h$beta)[1]),
                                'covariates' = c('intercept', covariate.cols),
                                'blood pressure' = c('sbp', 'dbp'))
  
  # preprocess data into useful formats
  output = setup.data(run.params, theta.h, dat, covariate.cols)
  datasets = output[['datasets']]
  
  # calculate how many days to predict for (i.e. the maximum any patient has)
  run.params[['n.predict.days']] = max(datasets$test.melt$adday)
  
  return(list(datasets = datasets,
              run.params = run.params,
              theta.h = theta.h))
}

#################################
# preprocess data
#################################
preprocess.data = function(data, covariate.cols)
{
  # create melted data frame
  
  # first, melt covariates
  melted1 = unique(
    data[,c('id', 'intercept', covariate.cols), drop = FALSE]
  )[ , c('intercept', covariate.cols), drop = FALSE]
  
  # summarize adherence for non-missing values
  melted2 = ddply(data[!is.na(data$adhere),], c('id'), summarise,
                  # number of days tracked
                  adday = length(adhere),
                  # yes adhere day
                  yesadday = sum(adhere > 0),
                  # no adhere day
                  noadday = sum(adhere <= 0),
                  # overall probability of adherence
                  padhere = sum((adhere > 0)/length(adhere)),
                  # orginal probablity of adherence
                  ogpadhere = sum((adhere > 0 & missing == 0)/sum(missing == 0))
  )
  
  # summarize blood pressure
  melted3 = ddply(data, c('id'), summarise,
                  # number of bp measurements
                  bpday = sum(!is.na(dbp) & !is.na(sbp))
  )
  
  # combine melted data frames
  melted = data.frame(melted2, bpday = melted3$bpday, melted1)
  
  # add back in these new measures for full data frame
  query.adhere = "
  SELECT data.*,
  melted.adday, melted.yesadday, melted.noadday,
  melted.padhere, melted.ogpadhere, melted.bpday
  FROM data
  LEFT JOIN melted ON data.id = melted.id
  "
  full = sqldf(query.adhere)
  
  return(list('full' = full, 'melted' = melted))
}

#################################
# process data from Luis Campos (LC)
#################################

process.data.lc = function(dat, covariate.cols)
{
  data = ldply(lapply(seq_along(dat), function(i, x){
    T = x[[i]]$T
    sbp = dbp = rep(NA, T)
    sbp[x[[i]]$I.y] = x[[i]]$y1[x[[i]]$I.y]
    dbp[x[[i]]$I.y] = x[[i]]$y2[x[[i]]$I.y]     
    cov.matrix = matrix(x[[i]]$x, nrow = T, ncol = length(x[[i]]$x), byrow = TRUE)
    matrix.row = cbind(
      rep(i, T),
      seq(1, T),
      x[[i]]$c_t,
      cov.matrix,
      sbp = sbp,
      dbp = dbp,
      missing = is.na(x[[i]]$c_t)
    )
    colnames(matrix.row) = c(
      'id', 'mday', 'adhere', 'intercept', covariate.cols, 'sbp', 'dbp', 'missing'
    )
    
    return(data.frame(matrix.row))
  }, x = dat))
  
  return(data)
}

#################################
# rearrange theta.a and theta.h params for convenience and readability
# only select a subsample for final inference
#################################

save.theta = function(theta.a, theta.h, run.params)
{
  # only select a subsample for inference
  theta.h.subsample = sample(
    seq(1, nrow(theta.h$rho)),
    run.params[['n.post.samp']],
    replace = TRUE
  )
  theta.a.subsample = sample(
    seq(1, nrow(theta.a[['fixef']])),
    run.params[['n.post.samp']],
    replace = TRUE
  )
  
  params.h = data.frame(
    rho.s = theta.h$rho[theta.h.subsample,1],
    rho.d = theta.h$rho[theta.h.subsample,2],
    phi.s = theta.h$rho[theta.h.subsample,1],
    phi.d = theta.h$rho[theta.h.subsample,2],
    sigma.s = theta.h$sigma[theta.h.subsample,1],
    sigma.d = theta.h$sigma[theta.h.subsample,2],
    gamma.sd = theta.h$cor[theta.h.subsample],
    sigma.s.nu = theta.h$sigma_nu[theta.h.subsample,1],
    sigma.d.nu = theta.h$sigma_nu[theta.h.subsample,2],
    sigma.s.0 = theta.h$sigma_0[theta.h.subsample,1],
    sigma.d.0 = theta.h$sigma_0[theta.h.subsample,2],
    beta.s = theta.h$beta[theta.h.subsample,,1],
    beta.d = theta.h$beta[theta.h.subsample,,2]
  )
  params.a = data.frame(
    beta.a = theta.a[['fixef']][theta.a.subsample,],
    sd.a.ranef = theta.a[['ranef.sd']][theta.a.subsample],
    alpha.new = theta.a[['alpha.new']][theta.a.subsample,]
  )
  
  params = cbind(params.h, params.a)
  return(params)
}

#################################
# setup data
#################################

setup.data = function(run.params, theta.h, dat, covariate.cols)
{
  # process data into right format
  data = process.data.lc(dat, covariate.cols)
  
  # data processing
  datasets = preprocess.data(data, covariate.cols)
  full = datasets[['full']]
  melted = datasets[['melted']]
  
  # split into training/testing
  ids = melted$id
  n.test = run.params[['test.size']]
  n.train = length(ids) - n.test
  set.seed(100)
  train.ids = sample(ids, n.train)
  test.ids = sample(ids[!(ids %in% train.ids)], n.test)
  
  # save out data
  datasets = list()
  datasets[['train.full']] = data[data$id %in% train.ids,]
  datasets[['test.full']] = data[data$id %in% test.ids,]
  
  datasets[['train.melt']] = melted[melted$id %in% train.ids,]
  datasets[['test.melt']] = melted[melted$id %in% test.ids,]
  datasets[['covariate.cols']] = covariate.cols
  
  return(list(datasets = datasets))
}
