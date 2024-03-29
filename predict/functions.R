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
    param.row['sigma.s.eps']^2,
    param.row['rho.eps'] * param.row['sigma.s.eps'] * param.row['sigma.d.eps'],
    param.row['rho.eps'] * param.row['sigma.s.eps'] * param.row['sigma.d.eps'],
    param.row['sigma.d.eps']^2
  ), nrow = 2, ncol = 2)
  return(e.matrix)
}

#################################
# get draws of theta a
#################################

get.theta.a.draws = function(theta.a.stan.dat, covariate.cols, run.params, base.dir)
{
  rstan_options(auto_write = TRUE)
  if(run.params[['cloud']])
  {
    options(mc.cores = parallel::detectCores())
  } else
  {
    options(mc.cores = parallel::detectCores() - 1)
  }

  # feed in binomial data, number of adherence days and number of total days
  N_test = run.params[['test.size']]
  y = theta.a.stan.dat$yesadday
  total = theta.a.stan.dat$adday
  x = theta.a.stan.dat[,covariate.cols, drop = FALSE]
  N = nrow(x)
  P = ncol(x)
  stan_data = list(N = N, N_test = N_test, P = P, y = y, x = x, total = total)

  # fit stan model
  fit = stan(
    file = paste0(base.dir, 'predict/adherence.stan'),
    data = stan_data,
    iter = run.params[['theta.a.mcmc.length']],
    chains = run.params[['mcmc.chains']]
  )

  # save out results
  fitted.values = extract(fit)
  theta.a.fixef = fitted.values$beta
  colnames(theta.a.fixef) = c('intercept', covariate.cols)
  theta.a = list(
    'fixef' = theta.a.fixef,
    'sigma.delta' = fitted.values$sigma_delta,
    'delta.train' = fitted.values$delta,
    'delta.new' = fitted.values$delta_new,
    'rhat' = summary(fit)$summary[,'Rhat'],
    'divergent' = get_num_divergent(fit)
  )

  theta.a$post.means = get_posterior_mean(fit)

  remove(fit)
  gc()

  return(theta.a)
}

#################################
# get draws of theta h
#################################

get.theta.h.draws = function(theta.h.stan.dat, covariate.cols,
                             run.params, base.dir)
{
  rstan_options(auto_write = TRUE)

  if(run.params[['cloud']])
  {
    options(mc.cores = parallel::detectCores())
  } else
  {
    options(mc.cores = parallel::detectCores() - 1)
  }

  start.time = Sys.time()

  if(run.params[['theta.h.fixed.init']])
  {
    init.params = initialize.chains.theta.h()[1:run.params[['mcmc.chains']]]

    fit = stan(
      file = paste0(base.dir, 'infer/state_space_adherence.stan'),
      data = theta.h.stan.dat,
      init = init.params,
      iter = run.params[['theta.h.mcmc.length']],
      warmup = run.params[['theta.h.mcmc.length']]/2,
      chains = run.params[['mcmc.chains']]
    )

  } else
  {
    fit = stan(
      file = paste0(base.dir, 'infer/state_space_adherence.stan'),
      data = theta.h.stan.dat,
      iter = run.params[['theta.h.mcmc.length']],
      warmup = run.params[['theta.h.mcmc.length']]/2,
      chains = run.params[['mcmc.chains']]
    )
  }

  end.time = Sys.time()

  total.time = difftime(end.time, start.time, units = 'secs')
  print(paste('State space model:', round(total.time), 'seconds or', round(total.time/60), 'minutes'))

  theta.h = extract(fit)
  dimnames(theta.h$beta) = list('iterations' = seq(1, dim(theta.h$beta)[1]),
                                'covariates' = c('intercept', covariate.cols),
                                'blood pressure' = c('sbp', 'dbp'))

  theta.h$rhat = summary(fit)$summary[,'Rhat']
  theta.h$divergent = get_num_divergent(fit)
  theta.h$model.time = difftime(end.time, start.time, units = 'secs')

  theta.h$post.means = get_posterior_mean(fit)

  remove(fit)
  gc()

  return(theta.h)
}

#################################
# rearrange theta.a and theta.h params for convenience and readability
# only select a subsample for final inference
#################################

save.theta.draws = function(theta.a, theta.h, run.params, delta.new = TRUE)
{
  # only select a subsample for inference
  theta.h.subsample = sample(
    seq(1, nrow(theta.h$rho)),
    run.params[['npostsamp']],
    replace = TRUE
  )
  theta.a.subsample = sample(
    seq(1, nrow(theta.a[['fixef']])),
    run.params[['npostsamp']],
    replace = TRUE
  )

  params.b = data.frame(
    rho.s = theta.h$rho[theta.h.subsample,1, drop = FALSE],
    rho.d = theta.h$rho[theta.h.subsample,2, drop = FALSE],
    phi.s = theta.h$phi[theta.h.subsample,1, drop = FALSE],
    phi.d = theta.h$phi[theta.h.subsample,2, drop = FALSE],
    sigma.s.eps = theta.h$sigma[theta.h.subsample,1, drop = FALSE],
    sigma.d.eps = theta.h$sigma[theta.h.subsample,2, drop = FALSE],
    rho.eps = theta.h$cor[theta.h.subsample],
    sigma.s.nu = theta.h$sigma_nu[theta.h.subsample,1, drop = FALSE],
    sigma.d.nu = theta.h$sigma_nu[theta.h.subsample,2, drop = FALSE],
    sigma.s.0 = theta.h$sigma_0[theta.h.subsample,1, drop = FALSE],
    sigma.d.0 = theta.h$sigma_0[theta.h.subsample,2, drop = FALSE],
    beta.s = theta.h$beta[theta.h.subsample,,1, drop = FALSE],
    beta.d = theta.h$beta[theta.h.subsample,,2, drop = FALSE]
  )
  params.a = data.frame(
    beta.a = theta.a[['fixef']][theta.a.subsample, , drop = FALSE],
    sigma.delta = theta.a[['sigma.delta']][theta.a.subsample],
    delta.train = theta.a[['delta.train']][theta.a.subsample, , drop = FALSE]
  )

  if(delta.new)
  {
    params.a$delta.new = theta.a[['delta.new']][theta.a.subsample, , drop = FALSE]
  }

  params = cbind(params.b, params.a)
  return(params)
}

#################################
# calculate quantiles for posterior adherence intervals
#################################

get.quantiles = function(c.means, p.prior, test.melt)
{
  # calculate quantiles
  ad.quantiles = data.frame(
    id = test.melt$id,
    padhere = test.melt$ogpadhere,
    bpday = test.melt$bpday,
    adday = test.melt$adday
  )

  # quantiles
  for(i in 1:nrow(c.means))
  {
    ad.quantiles$prior.mean[i] = mean(p.prior[i,])
    ad.quantiles$lower.prior.95[i] = quantile(x = p.prior[i,], probs = 0.025)
    ad.quantiles$upper.prior.95[i] = quantile(x = p.prior[i,], probs = 0.975)
    ad.quantiles$lower.prior.80[i] = quantile(x = p.prior[i,], probs = 0.1)
    ad.quantiles$upper.prior.80[i] = quantile(x = p.prior[i,], probs = 0.9)
    ad.quantiles$lower.prior.50[i] = quantile(x = p.prior[i,], probs = 0.25)
    ad.quantiles$upper.prior.50[i] = quantile(x = p.prior[i,], probs = 0.75)

    ad.quantiles$lower.95[i] = quantile(x = c.means[i,], probs = 0.025)
    ad.quantiles$upper.95[i] = quantile(x = c.means[i,], probs = 0.975)
    ad.quantiles$lower.80[i] = quantile(x = c.means[i,], probs = 0.1)
    ad.quantiles$upper.80[i] = quantile(x = c.means[i,], probs = 0.9)
    ad.quantiles$lower.50[i] = quantile(x = c.means[i,], probs = 0.25)
    ad.quantiles$upper.50[i] = quantile(x = c.means[i,], probs = 0.75)
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

#################################
# inititalize chains for theta h model
#################################

initialize.chains.theta.h = function()
{
  init.params = list(
    'chain1' = list(
      'rho' = c(0.01, 0.01),
      'phi' = c(0, 0),
      'sigma' = c(0.2, 0.2),
      'cor' = 0,
      'sigma_nu' = c(0.2, 0.2),
      'sigma_0' = c(0.2, 0.2),
      'beta' = matrix(c(120, 80, 0, 0, 0, 0), nrow = 3, ncol = 2, byrow = TRUE)
    ),
    'chain2' = list(
      'rho' = c(0.5, 0.5),
      'phi' = c(2, 2),
      'sigma' = c(0.5, 0.5),
      'cor' = -0.5,
      'sigma_nu' = c(0.5, 0.5),
      'sigma_0' = c(0.5, 0.5),
      'beta' = matrix(c(120, 80, 1, 1, 1, 1), nrow = 3, ncol = 2, byrow = TRUE)
    ),
    'chain3' = list(
      'rho' = c(0.8, 0.8),
      'phi' = c(-2, -2),
      'sigma' = c(1, 1),
      'cor' = 0.5,
      'sigma_nu' = c(1, 1),
      'sigma_0' = c(1, 1),
      'beta' = matrix(c(120, 80, -1, -1, -1, -1), nrow = 3, ncol = 2, byrow = TRUE)
    ),
    'chain4' = list(
      'rho' = c(0.2, 0.2),
      'phi' = c(2, -2),
      'sigma' = c(2, 2),
      'cor' = 0,
      'sigma_nu' = c(0.1, 0.1),
      'sigma_0' = c(2, 2),
      'beta' = matrix(c(120, 80, 0.2, -3, -3, 0.2), nrow = 3, ncol = 2, byrow = TRUE)
    )
  )

  return(init.params)
}

theta.true.params = c(
  'rho.s' = 0.8,
  'rho.d' = 0.2,
  'phi.s' = 0.5,
  'phi.d' = 0.7,
  'sigma.s.eps' = 2.75,
  'sigma.d.eps' = 2.2,
  'rho.eps' = 0.55,
  'sigma.s.nu' = 1,
  'sigma.d.nu' = 0.8,
  'sigma.s.0' = 1.5,
  'sigma.d.0' = 1.3,
  'beta.s.intercept.sbp' = 120,
  'beta.s.gender.sbp' = -1.2,
  'beta.d.intercept.dbp' = 80,
  'beta.d.gender.dbp' = -1.2,
  'beta.a.intercept' = 0,
  'beta.a.gender' = -0.4,
  'sigma.delta' = 1.1
)
