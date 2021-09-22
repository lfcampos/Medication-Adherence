# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Simulates from bivariate state-space model
#	y_t = x*beta + alpha_t + eps_t               eps_t ~ N([0,0], Sigma.BP)
#	alpha_t = rho*alpha_{t-1} + phi*c_t + nu_t           nu_t ~ N([0,0], Sigma.nu)
#   alpha_1 ~ N([0,0], Sigma.0)
#
# input:
#	rho = [rho.S, rho.D]    - numeric: auto-correlation
#	phi = [phi.S, phi.D]    - numeric: time-dependent input coefficient
#	Sigma.nu = [sig.nu_S^2,      0    ]   - matrix: innovation variance
#              [     0    , sig.nu_D^2]
#	Sigma.BP = [sig.BP_S^2,   cov.BP  ]   - matrix: measurement error variance
#              [  cov.BP  , sig.BP_D^2]
#	Sigma.0  = [sig0.S^2,    0    ]   - matrix: initial state variance
#              [   0    , sig0.S^2]
#	beta   - numeric: covariates [Jx2]
#	T      - integer: number of days followed (length of c_t)
# n.obs  - integer: number of days y_t to be observed (< T):
#
# output (list):
#	T     - integer: number of days followed (length of c_t)
#	c_t   - vector (-1/1): observed adherence
#	y     - vector (numeric): complete observed response
#	x     - vector (numeric): covariates
#	mu    - vector (numeric): mean vector for y
#	Sigma - matrix (numeric): covariance matrix for y
#	I.y   - vecor(integer): indexes of t = 1, 2, ..., T to keep
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

simulate.SS_2D = function(T, rho, phi, Sigma.nu, Sigma.BP, Sigma.0, beta,
                          sigma.delta, beta.a,
                          c_t = NULL, n.obs = NULL, I.y = NULL){

  # Make binary x vector the length of beta with intercept
  x = c(1, sample(c(0, 1), nrow(beta) - 1, replace = TRUE))

  # simulate adherence
  delta = rnorm(1, 0, sigma.delta)
  eta = delta + x[1] * beta.a[1] + x[2] * beta.a[2]
  p = exp(eta)/(1 + exp(eta))
  if(is.null(c_t)) c_t = matrix(2*rbinom(T, 1, p) - 1, T, 1)

  # calculate complete mean vector
  mu1 = mu_T(T, x, c_t, rho[1], phi[1], Sigma.nu[1,1], Sigma.BP[1,1], Sigma.0[1,1], beta[,1])
  mu2 = mu_T(T, x, c_t, rho[2], phi[2], Sigma.nu[2,2], Sigma.BP[2,2], Sigma.0[2,2], beta[,2])

  # calculate complete covariance matrix
  Sigma1 = Sigma_T(T, rho[1], phi[1], Sigma.nu[1,1], Sigma.BP[1,1], Sigma.0[1,1], beta[,1])
  Sigma2 = Sigma_T(T, rho[2], phi[2], Sigma.nu[2,2], Sigma.BP[2,2], Sigma.0[2,2], beta[,2])
  Sigma12 = diag(Sigma.BP[1,2], nrow = T)

  mu = c(mu1, mu2)
  Sigma = rbind(cbind(Sigma1, Sigma12), cbind(Sigma12, Sigma2))

  # simulate complete observation
  y = rmvnorm(1, mean = mu, sigma = Sigma)[1,]

  y1 = y[1:T]
  y2 = y[-(1:T)]

  # simulate which to keep
  if(is.null(n.obs)) n.obs = T
  if(is.null(I.y)) I.y = sort(c(1, sample(2:T, n.obs-1)))

  # missing values: for now, no missingngess
  missing = rep(0, T)

  return(list(T = T, c_t = c_t, y1 = y1, y2 = y2, x = x, I.y = I.y, missing = missing))
}
