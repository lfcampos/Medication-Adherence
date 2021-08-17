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
  if(!dir.exists(out.dir)){ dir.create(out.dir) }

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
setup.all = function(base.dir, run.dir, onestep = FALSE)
{
  data.dir = paste(base.dir, 'data/', sep = '')

  # load in data and posterior draws of blood pressure model
  # theta h = blood pressure model, provided by previous analysis
  dat = readRDS(paste0(data.dir, 'dat.RDS'))
  covariate.cols = readRDS(paste0(data.dir, 'covariate_cols.RDS'))

  # copy over to run dir for saving
  saveRDS(dat, paste0(run.dir, 'data.RDS'))
  saveRDS(covariate.cols, paste0(run.dir, 'covariate_cols.RDS'))

  if(!onestep)
  {
    theta.h = readRDS(paste0(data.dir, 'stan_fit.RDS'))
    dimnames(theta.h$beta) = list('iterations' = seq(1, dim(theta.h$beta)[1]),
                                  'covariates' = c('intercept', covariate.cols),
                                  'blood pressure' = c('sbp', 'dbp'))
  } else
  {
    theta.h = NULL
  }

  # preprocess data into useful formats
  output = setup.data(run.params = run.params, dat = dat, covariate.cols = covariate.cols)
  datasets = output[['datasets']]

  # calculate how many days to predict for (i.e. the maximum any patient has)
  run.params[['npredictdays']] = max(datasets$full$data.melt$mday)

  # save out number of covariates
  run.params[['ncovariates']] = length(covariate.cols)

  return(list(datasets = datasets,
              run.params = run.params,
              theta.h = theta.h))
}

#################################
# preprocess data
#################################
melt.data = function(data, covariate.cols)
{
  # create melted data frame

  # first, melt covariates
  melted1 = unique(
    data[,c('id', covariate.cols), drop = FALSE]
  )[ , covariate.cols, drop = FALSE]

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
                  # orginal probablity of ahderence
                  ogpadhere = sum((adhere > 0 & missing == 0)/sum(missing == 0))
  )

  # other summaries
  melted3 = ddply(data, c('id'), summarise,
                  # number of bp measurements
                  bpday = sum(!is.na(dbp) & !is.na(sbp)),
                  # including missing days
                  mday = max(mday)
  )

  melted = data.frame(melted2, bpday = melted3$bpday, mday = melted3$mday, melted1)

  return(melted)
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
# setup data
#################################

setup.data = function(run.params, covariate.cols, dat)
{
  # process data into right format
  data = process.data.lc(dat, covariate.cols)

  # data processing
  melted = melt.data(data, covariate.cols)

  # split into training/testing
  ids = melted$id
  n.test = run.params[['test.size']]
  n.train = length(ids) - n.test
  set.seed(100)
  train.ids = sample(ids, n.train)
  test.ids = sample(ids[!(ids %in% train.ids)], n.test)
  all.ids = c(train.ids, test.ids)

  # subset down to only ids of interest
  data.full = data[data$id %in% all.ids,]
  data.melt = melted[melted$id %in% all.ids,]

  # properly set things to NA
  data.full$adhere[data.full$id %in% test.ids] = NA

  # mark out test patients
  data.melt$new.patient = data.melt$id %in% test.ids
  data.full$new.patient = data.full$id %in% test.ids

  # for fitting models
  data.stan = setup.stan.data(
    data.full, data.melt, covariate.cols, impute.missing = FALSE
  )
  data.train.stan = setup.stan.data(
    data.full[data.full$id %in% train.ids,],
    data.melt[data.melt$id %in% train.ids,],
    covariate.cols,
    impute.missing = TRUE
  )

  # save out data
  datasets = list(
    'train' = list(
      'data.full' = data.full[data.full$id %in% train.ids,],
      'data.melt' = data.melt[data.melt$id %in% train.ids,],
      'data.stan' = data.train.stan,
      'covariate.cols' = covariate.cols
    ),
    'test' = list(
      'data.full' = data.full[data.full$id %in% test.ids,],
      'data.melt' = data.melt[data.melt$id %in% test.ids,],
      'covariate.cols' = covariate.cols
    ),
    'full' = list(
      'data.full' = data.full,
      'data.melt' = data.melt,
      'data.stan' = data.stan,
      'covariate.cols' = covariate.cols
    )
  )

  return(list(datasets = datasets))
}

#################################
# setup stan data
#################################

setup.stan.data = function(data.full, data.melt, covariate.cols, impute.missing = FALSE)
{
  all.ids = unique(data.full$id)
  N = length(all.ids)

  # list form
  dat.list = list()
  for(k in 1:N)
  {
    dat.list[[k]] = data.full[data.full$id == all.ids[k],]
  }

  # calculate various quantities
  T = sapply(dat.list, function(x) { return(max(x$mday)) })
  ni = sapply(dat.list, function(x) { return(sum(!is.na(x$sbp))) })
  y = do.call('c', lapply(dat.list, function(x) return(c(x$sbp[!is.na(x$sbp)], x$dbp[!is.na(x$sbp)])) ))
  t = do.call('c', lapply(dat.list, function(x){ return(c(which(!is.na(x$sbp)), which(!is.na(x$sbp)))) }))
  yid = do.call('c', lapply(dat.list, function(x){ return(c(rep(1, sum(!is.na(x$sbp))), rep(2,sum(!is.na(x$dbp))))) }))
  X = as.matrix(cbind(intercept = 1, data.melt[,covariate.cols]))
  maxT = max(data.full$mday)

  # adherence
  c_matrix = matrix(NA, nrow = N, ncol = maxT)
  for(k in 1:N){
    c_t = dat.list[[k]]$adhere
    c_t[c_t == 0] = -1
    if(sum(is.na(c_t)) > 0 & impute.missing){
      t1 = sum(c_t > 0, na.rm = TRUE)
      t2 = sum(c_t < 0, na.rm = TRUE)
      # sample rate
      theta = rbeta(1, t1 + 1, t2 + 1)
      need.rep = is.na(c_t)
      # sample adherence
      c_t[need.rep] = 2*rbinom(sum(need.rep), 1, theta) - 1
    }
    c_matrix[k, 1:length(c_t)] = c_t
    if(length(c_t) < maxT){
      c_matrix[k, (length(c_t)+1):maxT] = 9999
    }
  }

  # compile into stan data
  data.stan = list(
    N = N,
    J = ncol(X),
    T = T,
    maxT = maxT,
    ni = ni,
    n = sum(ni),
    X = X,
    y = y,
    id = rep(1:N, 2*ni),
    t = t,
    yid = yid,
    c = c_matrix
  )
  return(data.stan)
}
