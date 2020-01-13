# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : functions_pf.R
# Programmer Name    : Kristen Hunter
#                      kristenbhunter@gmail.com
# 
# Last Updated       : Jan 2020
#
# Purpose            : Supporting functions for particle filter steps
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

#################################
# apply multivariate normal density across a mean vector
#################################

apply.dmvnorm = function(mean, x, sigma){
  return(bayesm::lndMvn(x = x, mu = mean, rooti = backsolve(chol(sigma),diag(2))))
}

#################################
# initialize x
#################################

rinit = function(nparticles, xdimension, theta.row, p.particles){
  
  # xdimension = model$xdimension
  x1 = matrix(ncol = nparticles, nrow = xdimension)
  # blood pressure
  x1[1,] = rnorm(nparticles, mean = 0, sd = theta.row[['sigma.s.0']])
  x1[2,] = rnorm(nparticles, mean = 0, sd = theta.row[['sigma.d.0']])
  # adherence
  x1[3,] = rbinom(nparticles, size = 1, prob = p.particles)
  return(x1)
}

#################################
# draw transition of x
#################################

rtransition = function(xparticles, theta.row, p.particles){
  nparticles = ncol(xparticles)
  xdimension = nrow(xparticles)
  nextx = matrix(ncol = nparticles, nrow = xdimension)
  
  # adherence
  nextx[3,] = rbinom(nparticles, size = 1, prob = p.particles)
  
  # blood pressure
  nextx[1,] = theta.row['rho.s'] * xparticles[1,] +
    theta.row['phi.s'] * nextx[3,] +
    rnorm(nparticles, mean = 0, sd = theta.row[['sigma.s.nu']])
  
  nextx[2,] = theta.row['rho.d'] * xparticles[2,] +
    theta.row['phi.d'] * nextx[3,] +
    rnorm(nparticles, mean = 0, sd = theta.row[['sigma.d.nu']])
  return(nextx)
}

#################################
# density of transition x
#################################

dtransition = function(next_x, xparticles, theta.row, p.particles){
  
  # next_x = ref_trajectory[,time]; xparticles = xlast
  
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

dmeasurement = function(observation, xparticles, theta.row, xbeta, construct.sigma, apply.dmvnorm){
  mean.vec = rbind(xparticles[1,] + xbeta[['xsbeta']], xparticles[2,] + xbeta[['xdbeta']])
  sigma = construct.sigma(theta.row)
  return(apply(mean.vec, 2, apply.dmvnorm, x = observation, sigma = sigma))
}

#################################
# save out model
#################################

model = list(rinit = rinit, rtransition = rtransition, dtransition = dtransition,
             dmeasurement = dmeasurement, xdimension = 3)

#################################
# main function
#################################

PF = function(nparticles, model, theta.row, p.particles, Y, xbeta,
              ref_trajectory = NULL, with_as = FALSE){

  # nparticles = sim.params[['nparticles']]; with_as = TRUE;
  # ref_trajectory = xpath

  datalength = nrow(Y)
  # create tree representation of the trajectories
  Tree = new(TreeClass, nparticles, 10*nparticles*model$xdimension, model$xdimension)
  
  # initialization
  xparticles = model$rinit(nparticles, model$xdimension, theta.row, p.particles)
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
    xparticles = model$rtransition(xparticles, theta.row, p.particles)
    
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
  new_trajectory = Tree$get_path(systematic_resampling_n(normweights, 1, runif(1)) - 1)

  return(new_trajectory)
}

#################################
# run particle filter for one draw of posterior sample
#################################

post.samp.draw.pf = function(theta.row, datasets, sim.params)
{
  # theta.row = unlist(theta[1,])

  # save stuff out
  ids = unique(datasets$test.melt$id)
  X = as.matrix(cbind(intercept = 1, datasets$test.melt[,datasets$covariate.cols, drop = FALSE]))
  
  # set up tree for storage
  library(CoupledCPF)
  module_tree <<- Rcpp::Module("module_tree", PACKAGE = "CoupledCPF")
  TreeClass <<- module_tree$Tree
  
  # calculate xbetas to save time later
  beta.s = as.matrix(theta.row[grep('beta.s', names(theta.row))])
  beta.d = as.matrix(theta.row[grep('beta.d', names(theta.row))])
  beta.a = as.matrix(theta.row[grep('beta.a', names(theta.row))])
  xsbeta = t(X %*% beta.s)
  xdbeta = t(X %*% beta.d)
  xabeta = t(X %*% beta.a)
  
  # calculate prior probability
  alpha = t(as.matrix(theta.row[grep('alpha.new', names(theta.row))]))
  p.prior = exp(xabeta + alpha)/(1 + exp(xabeta + alpha))
  
  # prepare to save
  ad.prior = matrix(NA, ncol = sim.params[['test.size']], nrow = 1)
  ad.star.mean = matrix(NA,
    nrow = sim.params[['test.size']], ncol = sim.params[['npf']] - sim.params[['burnin']]
  )
  paths = NULL
  
  for(k in 1:sim.params[['test.size']])
  {
    xbeta = list(xsbeta = xsbeta[k], xdbeta = xdbeta[k])
    patient = datasets$test.full[datasets$test.full$id == ids[k],]
    
    Y = cbind(patient$sbp, patient$dbp)
    
    datalength = nrow(Y)
    times = patient$mday

    p.particles = rep(p.prior[k], sim.params[['nparticles']])
    ad.prior[,k] = p.prior[k]
    
    alpha.sbp.paths = matrix(nrow = sim.params[['npf']], ncol = datalength)
    alpha.dbp.paths = matrix(nrow = sim.params[['npf']], ncol = datalength)
    apaths = matrix(nrow = sim.params[['npf']], ncol = datalength)
    xpath = PF(sim.params[['nparticles']], model, theta.row, p.particles, Y, xbeta)
    alpha.sbp.paths[1,] = xpath[1,]
    alpha.dbp.paths[1,] = xpath[2,]
    apaths[1,] = xpath[3,]
    for (imcmc in 2:sim.params[['npf']]){
      xpath = PF(sim.params[['nparticles']], model, theta.row, p.particles,
                  Y, xbeta, ref_trajectory = xpath, with_as = TRUE)
      alpha.sbp.paths[imcmc,] = xpath[1,]
      alpha.dbp.paths[imcmc,] = xpath[2,]
      apaths[imcmc,] = xpath[3,]
    }

    # calculate means for each path
    selected.paths = apaths[(sim.params[['burnin']] + 1):sim.params[['npf']], , drop = FALSE]
    ad.star.mean[k,] = rowMeans(selected.paths, na.rm = TRUE)
  }
  
  return(list(ad.star.mean = ad.star.mean, ad.prior = ad.prior))
}

#################################
# set up parallelization
#################################

mkWorker.pf = function(datasets, sim.params) {
  
  bindToEnv(objNames = c('datasets', 'sim.params',
                         'post.samp.draw.pf', 'PF', 'construct.sigma', 'apply.dmvnorm',
                         'rinit', 'rtransition', 'dtransition', 'dmeasurement',
                         'model'))
  
  worker = function(theta.row) {
    post.samp.draw.pf(
      theta.row, datasets, sim.params
    )
  }
  return(worker)
}

#################################
# draw predicted adherence vector for a patient
#################################

draw.ad.star = function(datasets, theta, sim.params, cl = NULL)
{
  if(!sim.params[['test.code']])
  {
    post.output = parApply(
      cl, theta, 1,
      mkWorker.pf(datasets, sim.params)
    )
  } else
  {
    post.output = apply(
      theta, 1,
      mkWorker.pf(datasets, sim.params)
    )
  }
  
  # predicted probability with blood pressure
  ad.star.mean.list = lapply(post.output, function(x) x[['ad.star.mean']])
  ad.star.mean = abind(ad.star.mean.list, along = 3)
  ad.star.mean = aperm(ad.star.mean, c(1,3,2))
  dimnames(ad.star.mean) = list(
    'patient.id' = unique(datasets$test.melt$id),
    'post.sample.id' = seq(1, sim.params[['n.post.samp']], 1),
    'pf.id' = seq(1, sim.params[['npf']] - sim.params[['burnin']])
  )
  
  # predicted probability before blood pressure
  ad.prior.list = lapply(post.output, function(x) x[['ad.prior']])
  ad.prior = t(abind(ad.prior.list, along = 1))
  dimnames(ad.prior) = list(
    'id' = unique(datasets$test.melt$id),
    'post.sample.id' = seq(1, sim.params[['n.post.samp']], 1)
  )
  
  # re-arrange to make into a 2-D array instead of a 3-D array
  ad.star.mean.new = apply(ad.star.mean, 1, rbind)
  ad.star.mean.new = aperm(ad.star.mean.new, c(2,1))
  
  return(list('ad.means' = ad.star.mean.new, 'ad.prior' = ad.prior))
}
