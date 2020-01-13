# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : functions.R
# Programmer Name    : Kristen Hunter
#                      kristenbhunter@gmail.com
# 
# Last Updated       : Jan 2020
#
# Purpose            : General supporting functions for prediction analysis
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


#' Copy arguments into env and re-bind any function's lexical scope to bindTargetEnv .
#' 
#' See http://winvector.github.io/Parallel/PExample.html for example use.
#' 
#' 
#' Used to send data along with a function in situations such as parallel execution 
#' (when the global environment would not be available).  Typically called within 
#' a function that constructs the worker function to pass to the parallel processes
#' (so we have a nice lexical closure to work with).
#' 
#' @param bindTargetEnv environment to bind to
#' @param objNames additional names to lookup in parent environment and bind
#' @param names of functions to NOT rebind the lexical environments of
bindToEnv <- function(bindTargetEnv=parent.frame(),objNames,doNotRebind=c()) {
  # Bind the values into environment
  # and switch any functions to this environment!
  for(var in objNames) {
    val <- get(var,envir=parent.frame())
    if(is.function(val) && (!(var %in% doNotRebind))) {
      # replace function's lexical environment with our target (DANGEROUS)
      environment(val) <- bindTargetEnv
    }
    # assign object to target environment, only after any possible alteration
    assign(var,val,envir=bindTargetEnv)
  }
}

#################################
# inverse logit
#################################
inv.logit = function(x) { return(exp(x)/(exp(x) + 1))}

#################################
# construct sigma matrix for multivariate normal
#################################

construct.sigma = function(param.row)
{
  e.matrix = matrix(c(
    param.row['sigma.s']^2,
    param.row['gamma.sd'] * param.row['sigma.s'] * param.row['sigma.d'],
    param.row['gamma.sd'] * param.row['sigma.s'] * param.row['sigma.d'],
    param.row['sigma.d']^2
  ), nrow = 2, ncol = 2)
  return(e.matrix)
}

#################################
# get draws of theta a using stan
#################################

get.theta.a.draws.stan = function(datasets, run.params, predict.dir, run.dir)
{
  # construct dataset to feed into stan
  train.data = datasets[['train.melt']]
  N_test = run.params[['test.size']]
  # feed in binomial data, number of adherence days and number of total days
  y = train.data$yesadday
  total = train.data$adday
  x = train.data[,datasets[['covariate.cols']], drop = FALSE]
  N = nrow(x)
  P = ncol(x)
  # put into correct format for stan
  stan_data = list(N = N, N_test = N_test, P = P, y = y, x = x, total = total)
  
  # if running locally, don't use all cores
  rstan_options(auto_write = TRUE)
  if(run.params[['cloud']])
  { 
    options(mc.cores = parallel::detectCores())
  } else
  {
    options(mc.cores = parallel::detectCores() - 1)
  }
  
  # if testing the code, try a mini chain
  if(run.params[['test.code']])
  {
    mcmc.length = 100
    mcmc.chains = 1
  } else
  {
    mcmc.length = 80000
    mcmc.chains = 4
  }
  
  # fit stan model
  fitted.model = stan(
    file = paste(predict.dir, 'ran_ef.stan', sep = ''),
    model_name = 'ran_ef',
    data = stan_data,
    iter = mcmc.length * 1.2,
    warmup = mcmc.length/5,
    chains = mcmc.chains,
    control = list(max_treedepth = 15)
  )
  
  # save out results
  fitted.values = extract(fitted.model)
  theta.a.fixef = fitted.values$beta
  colnames(theta.a.fixef) = c('intercept', datasets[['covariate.cols']])
  theta.a = list(
    'fixef' = theta.a.fixef,
    'ranef.sd' = fitted.values$sigma_alpha,
    'alpha.new' = fitted.values$alpha_new
  )
  
  # diagnostics
  model.summary = summary(fitted.model)$summary
  pars.keep = 
    grepl('beta', rownames(model.summary)) |
    grepl('sigma_alpha', rownames(model.summary))
  model.diag = model.summary[pars.keep,]
  write.csv(model.diag,
    file = paste(run.dir, 'model_diagnostics.csv', sep = ''))
  pars = rownames(model.summary)[pars.keep]
  
  return(list(theta.a = theta.a, fitted.model = fitted.model, pars = pars))
}

#################################
# calculate quantiles for posterior adherence intervals
#################################

get.quantiles = function(ad.means, ad.prior, test.melt)
{
  # calculate quantiles
  ad.quantiles = data.frame(
    id = test.melt$id,
    padhere = test.melt$ogpadhere,
    bpday = test.melt$bpday,
    adday = test.melt$adday
  )
  
  # quantiles
  for(i in 1:nrow(ad.means))
  {
    ad.quantiles$prior.mean[i] = mean(ad.prior[i,])
    ad.quantiles$lower.prior.95[i] = quantile(x = ad.prior[i,], probs = 0.025)
    ad.quantiles$upper.prior.95[i] = quantile(x = ad.prior[i,], probs = 0.975)
    ad.quantiles$lower.prior.80[i] = quantile(x = ad.prior[i,], probs = 0.1)
    ad.quantiles$upper.prior.80[i] = quantile(x = ad.prior[i,], probs = 0.9)
    ad.quantiles$lower.prior.50[i] = quantile(x = ad.prior[i,], probs = 0.25)
    ad.quantiles$upper.prior.50[i] = quantile(x = ad.prior[i,], probs = 0.75)
    
    ad.quantiles$lower.95[i] = quantile(x = ad.means[i,], probs = 0.025)
    ad.quantiles$upper.95[i] = quantile(x = ad.means[i,], probs = 0.975)
    ad.quantiles$lower.80[i] = quantile(x = ad.means[i,], probs = 0.1)
    ad.quantiles$upper.80[i] = quantile(x = ad.means[i,], probs = 0.9)
    ad.quantiles$lower.50[i] = quantile(x = ad.means[i,], probs = 0.25)
    ad.quantiles$upper.50[i] = quantile(x = ad.means[i,], probs = 0.75)
  }
  
  # coverage
  calc.cover = function(true, lower, upper)
  {
    cover =
      round(lower, 2) <= round(true, 2) &
      round(true, 2) <= round(upper, 2)
    return(cover)
  }
  ad.quantiles$cover.95 = calc.cover(
    ad.quantiles$padhere, ad.quantiles$lower.95, ad.quantiles$upper.95
  )
  ad.quantiles$cover.80 = calc.cover(
    ad.quantiles$padhere, ad.quantiles$lower.80, ad.quantiles$upper.80
  )
  ad.quantiles$cover.50 = calc.cover(
    ad.quantiles$padhere, ad.quantiles$lower.50, ad.quantiles$upper.50
  )
  
  # coverage for 'prior' intervals
  ad.quantiles$cover.prior.95 = calc.cover(
    ad.quantiles$padhere, ad.quantiles$lower.prior.95, ad.quantiles$upper.prior.95
  )
  ad.quantiles$cover.prior.80 = calc.cover(
    ad.quantiles$padhere, ad.quantiles$lower.prior.80, ad.quantiles$upper.prior.80
  )
  ad.quantiles$cover.prior.50 = calc.cover(
    ad.quantiles$padhere, ad.quantiles$lower.prior.50, ad.quantiles$upper.prior.50
  )
  
  # interval width
  ad.quantiles$width.95 = ad.quantiles$upper.95 - ad.quantiles$lower.95
  ad.quantiles$width.80 = ad.quantiles$upper.80 - ad.quantiles$lower.80
  ad.quantiles$width.50 = ad.quantiles$upper.50 - ad.quantiles$lower.50
  
  # interval width for 'prior' intervals
  ad.quantiles$width.prior.95 = ad.quantiles$upper.prior.95 - ad.quantiles$lower.prior.95
  ad.quantiles$width.prior.80 = ad.quantiles$upper.prior.80 - ad.quantiles$lower.prior.80
  ad.quantiles$width.prior.50 = ad.quantiles$upper.prior.50 - ad.quantiles$lower.prior.50

  return(ad.quantiles)
}