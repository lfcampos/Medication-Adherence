# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : functions_pf.R
# Programmer Name    : Pierre Jacob
#                      with edits by Kristen Hunter
#                      kristenbhunter@gmail.com
#
# Source             : CoupledCPF repo by pierrejacob
#                    : https://github.com/pierrejacob/CoupledCPF
#
# Last Updated       : Jan 2020
#
# Purpose            : Supporting functions for particle filter steps
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

#################################
# apply multivariate normal density across a mean vector
#################################

apply.dmvnorm = function(mean, x, sigma)
{
  return(bayesm::lndMvn(x = x, mu = mean, rooti = backsolve(chol(sigma),diag(2))))
}

#################################
# initialize x
#################################

rinit = function(xobs.init, nparticles, xdimension, theta.row, p.particles)
{
  x1 = matrix(ncol = nparticles, nrow = xdimension)
  # blood pressure
  x1[1,] = rnorm(nparticles, mean = 0, sd = theta.row[['sigma.s.0']])
  x1[2,] = rnorm(nparticles, mean = 0, sd = theta.row[['sigma.d.0']])
  # adherence
  if(is.na(xobs.init))
  {
    x1[3,] = rbinom(nparticles, size = 1, prob = p.particles)
  } else
  {
    x1[3,] = rep(xobs.init, nparticles)
  }
  return(x1)
}

#################################
# draw transition of x
#################################

rtransition = function(xobs.time, xparticles, theta.row, p.particles)
{
  nparticles = ncol(xparticles)
  xdimension = nrow(xparticles)
  nextx = matrix(ncol = nparticles, nrow = xdimension)

  # adherence
  if(!is.na(xobs.time))
  {
    nextx[3,] = rep(xobs.time, nparticles)
  } else
  {
    nextx[3,] = rbinom(nparticles, size = 1, prob = p.particles)
  }

  # blood pressure
  nextx[1,] = theta.row['rho.s'] * xparticles[1,] +
    theta.row['phi.s'] * nextx[3,, drop = FALSE] +
    rnorm(nparticles, mean = 0, sd = theta.row[['sigma.s.nu']])

  nextx[2,] = theta.row['rho.d'] * xparticles[2,] +
    theta.row['phi.d'] * nextx[3,] +
    rnorm(nparticles, mean = 0, sd = theta.row[['sigma.d.nu']])
  return(nextx)
}

#################################
# density of transition x
#################################

dtransition = function(next_x, xparticles, theta.row, p.particles)
{
  mean.s = theta.row['rho.s'] * xparticles[1,] + theta.row['phi.s'] * next_x[3]
  mean.d = theta.row['rho.d'] * xparticles[2,] + theta.row['phi.d'] * next_x[3]
  d.s = dnorm(next_x[1], mean.s, sd = theta.row[['sigma.s.nu']], log = TRUE)
  d.d = dnorm(next_x[2], mean.d, sd = theta.row[['sigma.d.nu']], log = TRUE)
  d.c = dbinom(next_x[3], size = 1, prob = p.particles, log = TRUE)
  return(d.s + d.d + d.c)
}

#################################
# density of measurement
#################################

dmeasurement = function(observation, xparticles, theta.row,
                        xbeta, construct.sigma, apply.dmvnorm)
{
  mean.vec = rbind(xparticles[1,] + xbeta[['xsbeta']], xparticles[2,] + xbeta[['xdbeta']])
  sigma = construct.sigma(theta.row)
  return(apply(mean.vec, 2, apply.dmvnorm, x = observation, sigma = sigma))
}

#################################
# save out model
#################################

model = list(rinit = rinit, rtransition = rtransition,
             dtransition = dtransition,
             dmeasurement = dmeasurement, xdimension = 3)

#################################
# main function
#################################

PF = function(nparticles, model, theta.row, p.particles, xobs, Y, xbeta,
              ref_trajectory = NULL, with_as = FALSE)
{

  datalength = nrow(Y)
  # store ess
  ess = rep(NA, datalength)
  # create tree representation of the trajectories
  Tree = new(TreeClass, nparticles, 10*nparticles*model$xdimension, model$xdimension)

  # initialization
  xparticles = model$rinit(xobs[1], nparticles, model$xdimension, theta.row, p.particles)
  if (!is.null(ref_trajectory)){
    xparticles[ , nparticles] = ref_trajectory[,1]
  }
  Tree$init(xparticles)

  # if ancestor sampling, needs to keep the last generation of particles at hand
  if (with_as){
    xlast = xparticles
  }
  # weight with Y1
  if (is.na(Y[1,1])){
    normweights = rep(1/nparticles, nparticles)
  } else {
    logw = model$dmeasurement(Y[1,], xparticles, theta.row, xbeta,
                              construct.sigma, apply.dmvnorm)
    maxlw = max(logw)
    w = exp(logw - maxlw)
    normweights = w / sum(w)
  }
  # step t > 1
  for (time in 2:datalength){
    ancestors = multinomial_resampling_n(normweights, nparticles)
    xparticles = xparticles[,ancestors]
    xparticles = model$rtransition(xobs.time = xobs[time], xparticles, theta.row, p.particles)

    if (!is.null(ref_trajectory)){
      xparticles[,nparticles] = ref_trajectory[,time]
      if (with_as){
        # Ancestor sampling
        logm = model$dtransition(ref_trajectory[,time], xlast, theta.row, p.particles)
        logm = log(normweights) + logm
        w_as = exp(logm - max(logm))
        w_as = w_as / sum(w_as)
        ancestors[nparticles] = systematic_resampling_n(w_as, 1, runif(1))
        x_last = xparticles
      } else {
        ancestors[nparticles] = nparticles
      }
    }
    if (is.na(Y[time, 1])){
      normweights = rep(1/nparticles, nparticles)
    } else {
      logw = model$dmeasurement(Y[time, ], xparticles, theta.row, xbeta,
                                construct.sigma, apply.dmvnorm)
      maxlw = max(logw)
      w = exp(logw - maxlw)
      normweights = w / sum(w)
    }
    Tree$update(xparticles, ancestors - 1)
  }
  ess = 1/sum(normweights^2)
  new_trajectory = Tree$get_path(systematic_resampling_n(normweights, 1, runif(1)) - 1)

  return(list(new_trajectory = new_trajectory, ess = ess))
}

#################################
# run particle filter for one draw of posterior sample
#################################

post.samp.draw.pf = function(theta.row, datasets, run.params)
{
  # theta.row = unlist(theta[1,])

  # set up tree for storage
  library(CoupledCPF)
  module_tree <<- Rcpp::Module("module_tree", PACKAGE = "CoupledCPF")
  TreeClass <<- module_tree$Tree

  # preparation
  ids = unique(datasets$data.full$id)
  X = as.matrix(cbind(intercept = 1, datasets$data.melt[,datasets$covariate.cols]))

  # calculate xbetas to save time later
  beta.s = as.matrix(theta.row[grep('beta.s', names(theta.row))])
  beta.d = as.matrix(theta.row[grep('beta.d', names(theta.row))])
  beta.a = as.matrix(theta.row[grep('beta.a', names(theta.row))])
  xsbeta = t(X %*% beta.s)
  xdbeta = t(X %*% beta.d)
  xabeta = t(X %*% beta.a)

  n.patients = nrow(X)

  # calculate prior probability
  delta = matrix(NA, ncol = nrow(X))
  # if an existing patient, random intercept is known
  delta[!datasets$data.melt$new.patient] = theta.row[grep('delta.train', names(theta.row))]
  delta[datasets$data.melt$new.patient] = theta.row[grep('delta.new', names(theta.row))]
  p.prior = exp(xabeta + delta)/(1 + exp(xabeta + delta))

  # prepare to save out main paths
  c.star = array(
    NA, dim = c(n.patients, run.params[['npf']], run.params[['npredictdays']])
  )
  ess = matrix(
    NA, nrow = n.patients, ncol = run.params[['npf']]
  )
  xpaths = list()

  for(e in 1:run.params[['npf']])
  {
    print('---------------------------------------------------------------------------------------------')
    print(paste('Particle filter iteration:', e, 'out of', run.params[['npf']]))
    if(e > 2)
    {
      time.now = Sys.time()
      print(paste('Time elapsed:', timediff(time.now, time.start.pf)))
    }
    print('---------------------------------------------------------------------------------------------')

    if(e == 2)
    {
      time.start.pf = Sys.time()
    }

    for(k in 1:n.patients)
    {
      patient = datasets$data.full[datasets$data.full$id == ids[k],]
      xbeta = list(xsbeta = xsbeta[k], xdbeta = xdbeta[k])

      Y = cbind(patient$sbp, patient$dbp)
      if(patient$new.patient[1])
      {
        xobs = rep(NA, run.params[['npredictdays']])
      } else
      {
        xobs = patient$adhere
        xobs[!is.na(xobs) & xobs == -1] = 0
      }
      times = patient$mday

      p.particles = rep(p.prior[k], run.params[['nparticles']])

      if(e == 1)
      {
        pf.output = PF(run.params[['nparticles']], model, theta.row, p.particles, xobs, Y, xbeta)
      } else
      {
        pf.output = PF(run.params[['nparticles']], model, theta.row, p.particles,
                       xobs, Y, xbeta, ref_trajectory = xpaths[[k]], with_as = TRUE)
      }
      xpaths[[k]] = pf.output[['new_trajectory']]
      ess[k,e] = pf.output[['ess']]
      remove(pf.output)
      c.star[k,e,] = c(xpaths[[k]][3,], rep(NA, run.params[['npredictdays']] - length(xpaths[[k]][3,])))
    }

    if(e == 2)
    {
      time.end.pf = Sys.time()
      print('---------------------------------------------------------------------------------------------')
      print('Runtime statistics')
      print(paste(
        'One iteration took',
        round(difftime(time.end.pf, time.start.pf, units = 'secs')[[1]], 2),
        'seconds'
      ))
      print(paste(
        'Expected total run time:',
        round( (run.params[['npf']]*difftime(time.end.pf, time.start.pf, units = 'secs')[[1]])/60),
        'minutes'
      ))
      print('---------------------------------------------------------------------------------------------')
    }
  }

  # calculate means for each path
  selected.paths = c.star[,(run.params[['burnin']] + 1):run.params[['npf']], , drop = FALSE]
  c.star.means = rowMeans(selected.paths, na.rm = TRUE, dims = 2)
  # remove burn-in for ess
  ess = ess[,(run.params[['burnin']] + 1):run.params[['npf']]]
  # only keep test patients
  c.star.means = c.star.means[datasets$data.melt$new.patient, , drop = FALSE]

  return(list(c.star.means = c.star.means, p.prior = p.prior, ess = ess))
}


#################################
# run particle filter for one draw of posterior sample
#################################

post.samp.draw.pf.onestep = function(theta.row, datasets, run.params, base.dir)
{
  # set up tree for storage
  library(CoupledCPF)
  module_tree <<- Rcpp::Module("module_tree", PACKAGE = "CoupledCPF")
  TreeClass <<- module_tree$Tree

  # preparation
  ids = unique(datasets$full$data.full$id)
  X = as.matrix(cbind(intercept = 1, datasets$full$data.melt[,datasets$full$covariate.cols]))
  n.patients = nrow(X)

  # prepare to save out main paths
  c.star = array(
    NA, dim = c(n.patients, run.params[['npf']], run.params[['npredictdays']])
  )
  ess = matrix(
    NA, nrow = n.patients, ncol = run.params[['npf']]
  )
  p.prior.all = matrix(
    NA, nrow = n.patients, ncol = run.params[['npf']]
  )
  theta.all = matrix(
    NA, nrow = length(theta.row), ncol = run.params[['npf']]
  )
  xpaths = list()

  # split into old and new patients
  new.patient.rows = datasets$full$data.melt$new.patient
  old.patient.rows = !datasets$full$data.melt$new.patient

  # save out convergence
  n.theta.h.params = 12 + 2 * (1 + run.params[['ncovariates']])
  n.theta.a.params = 2 + (1 + run.params[['ncovariates']]) + n.patients
  theta.a.rhat = matrix(NA, nrow = n.theta.a.params, ncol = run.params[['npf']])
  theta.h.rhat = matrix(NA, n.theta.h.params, ncol = run.params[['npf']])

  theta.a.divergent = rep(NA, run.params[['npf']])
  theta.h.divergent = rep(NA, run.params[['npf']])

  for(e in 1:run.params[['npf']])
  {
    print('---------------------------------------------------------------------------------------------')
    print(paste('Particle filter iteration:', e, 'out of', run.params[['npf']]))
    print('---------------------------------------------------------------------------------------------')

    if(e == 2)
    {
      time.start.pf = Sys.time()
    }

    # calculate xbetas to save time later
    beta.s = as.matrix(theta.row[grep('beta.s', names(theta.row))])
    beta.d = as.matrix(theta.row[grep('beta.d', names(theta.row))])
    beta.a = as.matrix(theta.row[grep('beta.a', names(theta.row))])
    xsbeta = t(X %*% beta.s)
    xdbeta = t(X %*% beta.d)
    xabeta = t(X %*% beta.a)

    # calculate prior probability
    delta = matrix(NA, ncol = nrow(X))
    # if an existing patient, random intercept is known
    if(e == 1)
    {
      delta[old.patient.rows] = theta.row[grep('delta.train', names(theta.row))]
      delta[new.patient.rows] = theta.row[grep('delta.new', names(theta.row))]
    } else
    {
      delta = theta.row[grep('delta.train', names(theta.row))]
    }

    p.prior = exp(xabeta + delta)/(1 + exp(xabeta + delta))
    p.prior.all[,e] = p.prior

    for(k in 1:n.patients)
    {
      patient = datasets$full$data.full[datasets$full$data.full$id == ids[k],]
      xbeta = list(xsbeta = xsbeta[k], xdbeta = xdbeta[k])

      Y = cbind(patient$sbp, patient$dbp)
      if(patient$new.patient[1])
      {
        xobs = rep(NA, run.params[['npredictdays']])
      } else
      {
        xobs = patient$adhere
        xobs[!is.na(xobs) & xobs == -1] = 0
      }
      times = patient$mday

      p.particles = rep(p.prior[k], run.params[['nparticles']])

      if(e == 1)
      {
        pf.output = PF(run.params[['nparticles']], model,
                       theta.row, p.particles, xobs, Y, xbeta)
      } else
      {
        pf.output = PF(run.params[['nparticles']], model, theta.row, p.particles,
                       xobs, Y, xbeta, ref_trajectory = xpaths[[k]], with_as = TRUE)
      }
      xpaths[[k]] = pf.output[['new_trajectory']]
      ess[k,e] = pf.output[['ess']]
      remove(pf.output)
      c.star[k,e,] = c(xpaths[[k]][3,], rep(NA, run.params[['npredictdays']] - length(xpaths[[k]][3,])))
    }

    if(e == 2)
    {
      time.end.pf = Sys.time()
    }

    if(e < 10)
    {
      theta_b_stan_dat = datasets$train$data.stan
      theta_a_stan_dat = datasets$train$data.melt

      theta.h = get.theta.h.draws(theta_b_stan_dat, covariate.cols = datasets$full[['covariate.cols']], run.params, base.dir)
      theta.a = get.theta.a.draws(theta_a_stan_dat, covariate.cols = datasets$full[['covariate.cols']], run.params, base.dir)
      theta.iter = save.theta(theta.a, theta.h, run.params, delta.new = TRUE)

      theta.h.rhat[,e] = theta.h$rhat
      rownames(theta.h.rhat) = names(theta.h$rhat)

      theta.a.rhat.iter = theta.a$rhat
      param.order = names(theta.a$rhat)
      param.order = c(
        param.order[grepl('beta', param.order)],
        param.order[grepl('delta', param.order) & !(grepl('sigma', param.order))],
        'sigma_delta', 'lp__'
      )
      theta.a.rhat[,e] = theta.a.rhat.iter[param.order]
      rownames(theta.a.rhat) = param.order

      theta.a.divergent[e] = theta.a[['divergent']]
      theta.h.divergent[e] = theta.h[['divergent']]
    } else
    {
      c.star.pf = c.star[,e,]
      theta_b_stan_dat = datasets$full$data.stan
      c_matrix = matrix(NA, nrow = n.patients, ncol = run.params[['npredictdays']])
      for(k in 1:n.patients){
        c_t = c.star.pf[k,]
        c_t[c_t == 0] = -1
        c_t[is.na(c_t)] = 9999
        c_matrix[k,] = c_t
      }
      theta_b_stan_dat[['c']] = c_matrix

      theta_a_stan_dat = datasets$full$data.melt
      theta_a_stan_dat$adday = theta_a_stan_dat$mday
      theta_a_stan_dat$yesadday = apply(c.star.pf, 1, sum, na.rm = TRUE)

      theta.h = get.theta.h.draws(theta_b_stan_dat, covariate.cols = datasets$full[['covariate.cols']], run.params, base.dir)
      theta.a = get.theta.a.draws(theta_a_stan_dat, covariate.cols = datasets$full[['covariate.cols']], run.params, base.dir)
      theta.iter = save.theta(theta.a, theta.h, run.params, delta.new = FALSE)

      theta.h.rhat[,e] = theta.h$rhat
      rownames(theta.h.rhat) = names(theta.h$rhat)

      theta.a.rhat.iter = theta.a$rhat
      param.order = names(theta.a$rhat)
      param.order = c(
        param.order[grepl('beta', param.order)],
        param.order[grepl('delta', param.order) & !(grepl('sigma', param.order))],
        'sigma_delta', 'lp__'
      )
      theta.a.rhat[,e] = theta.a.rhat.iter[param.order]
      rownames(theta.a.rhat) = param.order

      theta.a.divergent[e] = theta.a[['divergent']]
      theta.h.divergent[e] = theta.h[['divergent']]
    }

    theta.row = unlist(theta.iter[1,])
    theta.all[,e] = theta.row
    remove(theta.a, theta.h)
    gc()

    # Keep track of simulation runtime
    if(e == 2)
    {
      time.end.learn = Sys.time()
      print('---------------------------------------------------------------------------------------------')
      print('Runtime statistics')
      print(paste('One iteration took', round(difftime(time.end.learn, time.start.pf, units = 'secs')[[1]], 2), 'seconds'))
      print(paste('PF time:', round(difftime(time.end.pf, time.start.pf, units = 'secs')[[1]], 2), 'seconds'))
      print(paste('Learn time:', round(difftime(time.end.learn, time.end.pf, units = 'secs')[[1]], 2), 'seconds'))
      print(paste('Expected total run time:', round( (run.params[['npf']]*difftime(time.end.learn, time.start.pf, units = 'secs')[[1]])/60), 'minutes'))
      print('---------------------------------------------------------------------------------------------')
    }

    # save out intermediate results
    c.star.means.iter = rowMeans(c.star[,1:e,, drop = FALSE], na.rm = TRUE, dims = 2)
    c.means.iter = c.star.means.iter[new.patient.rows, , drop = FALSE]
    p.prior.iter = p.prior.all[new.patient.rows, 1:e, drop = FALSE]
    ess.iter = ess[new.patient.rows, 1:e, drop = FALSE]

    draws = list(
      'c.means' = c.means.iter,
      'p.prior' = p.prior.iter,
      'ess' = ess.iter,
      'theta' = theta.all[,1:e, drop = FALSE],
      'theta.a.rhat' = theta.a.rhat[,1:e, drop = FALSE],
      'theta.h.rhat' = theta.h.rhat[,1:e, drop = FALSE],
      'theta.a.divergent' = theta.a.divergent[1:e],
      'theta.h.divergent' = theta.h.divergent[1:e]
    )
    saveRDS(draws, file = paste0(run.dir, 'draws_onestep.rds'))
  }

  # calculate means for each path
  keep.iter = (run.params[['burnin']] + 1):run.params[['npf']]
  selected.paths = c.star[,keep.iter, ,drop = FALSE]
  c.star.means = rowMeans(selected.paths, na.rm = TRUE, dims = 2)
  # remove burn-in and only keep test patients
  ess = ess[new.patient.rows, keep.iter]
  c.star.means = c.star.means[new.patient.rows, , drop = FALSE]
  p.prior.all = p.prior.all[new.patient.rows, keep.iter, drop = FALSE]
  theta.all = theta.all[, keep.iter, drop = FALSE]

  return(list(c.star.means = c.star.means, p.prior = p.prior.all, ess = ess, theta = theta.all,
              theta.a.rhat = theta.a.rhat, theta.h.rhat = theta.h.rhat,
              theta.a.divergent = theta.a.divergent, theta.h.divergent = theta.h.divergent))
}

#################################
# set up parallelization
#################################

mkWorker.pf = function(datasets, run.params) {

  bindToEnv(objNames = c('datasets', 'run.params',
                         'post.samp.draw.pf', 'PF', 'construct.sigma', 'apply.dmvnorm',
                         'rinit', 'rtransition', 'dtransition', 'dmeasurement',
                         'model'))

  worker = function(theta.row) {
    post.samp.draw.pf(
      theta.row, datasets, run.params
    )
  }
  return(worker)
}

#################################
# draw predicted adherence vector for all patients
#################################

draw.c.star = function(datasets, theta, run.params, cl = NULL)
{
  if(!run.params[['test.code']])
  {
    post.output = parApply(
      cl, theta, 1,
      mkWorker.pf(datasets, run.params)
    )
  } else
  {
    post.output = apply(
      theta, 1,
      mkWorker.pf(datasets, run.params)
    )
  }

  # predicted probability with blood pressure
  c.star.mean.list = lapply(post.output, function(x) x[['c.star.mean']])
  c.star.mean = abind(c.star.mean.list, along = 3)
  c.star.mean = aperm(c.star.mean, c(1,3,2))
  dimnames(c.star.mean) = list(
    'patient.id' = unique(datasets$test.melt$id),
    'post.sample.id' = seq(1, run.params[['n.post.samp']], 1),
    'pf.id' = seq(1, run.params[['npf']] - run.params[['burnin']])
  )

  # predicted probability before blood pressure
  ad.prior.list = lapply(post.output, function(x) x[['ad.prior']])
  ad.prior = t(abind(ad.prior.list, along = 1))
  dimnames(ad.prior) = list(
    'id' = unique(datasets$test.melt$id),
    'post.sample.id' = seq(1, run.params[['n.post.samp']], 1)
  )

  # re-arrange to make into a 2-D array instead of a 3-D array
  c.star.mean.new = apply(c.star.mean, 1, rbind)
  c.star.mean.new = aperm(c.star.mean.new, c(2,1))

  return(list('ad.means' = c.star.mean.new, 'ad.prior' = ad.prior))
}


#################################
# draw predicted adherence vector for all patients
# one-step procedure
#################################

draw.c.star.onestep = function(datasets, run.params, base.dir)
{
  # begin with a warm start: initial draw of theta parameters
  theta.h.stan.dat = datasets[['train']]$data.stan
  theta.a.stan.dat = datasets[['train']]$data.melt
  covariate.cols = datasets[['train']][['covariate.cols']]

  theta.h = get.theta.h.draws(theta.h.stan.dat, covariate.cols, run.params, base.dir)
  theta.a = get.theta.a.draws(theta.a.stan.dat, covariate.cols, run.params, base.dir)
  theta = save.theta(theta.a, theta.h, run.params, delta.new = TRUE)
  theta.row = unlist(theta[1,])

  # draw posterior samples
  post.output = post.samp.draw.pf.onestep(theta.row, datasets, run.params, base.dir)

  # predicted probability with blood pressure
  c.star.mean = post.output[['c.star.means']]
  p.prior = post.output[['p.prior']]
  ess = post.output[['ess']]
  dimnames(ess) = dimnames(p.prior) = dimnames(c.star.mean) = list(
    'id' = unique(datasets$data.melt$id),
    'pf.id' = seq(1, run.params[['npf']] - run.params[['burnin']])
  )


  return(list(
    'c.means' = c.star.mean,
    'p.prior' = p.prior,
    'ess' = ess,
    'theta' = post.output[['theta']],
    'theta.a.rhat' = post.output[['theta.a.rhat']],
    'theta.h.rhat' = post.output[['theta.h.rhat']],
    'theta.a.divergent' = post.output[['theta.a.divergent']],
    'theta.h.divergent' = post.output[['theta.h.divergent']]
  ))
}
