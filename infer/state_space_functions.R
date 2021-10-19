# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : state_space_functions.R
# Programmer Name    : Luis Campos
#                     soyluiscampos@gmail.com
#
# Purpose            : Supporting functions for simulation
#
# Input              : None
# Output             : None
#
# References         : Measuring Effects of Medication Adherence on Time-Varying
# 					   Health Outcomes using Bayesian Dynamic Linear Models
# 					   Luis F. Campos, Mark E. Glickman, Kristen B. Hunter
# 					   https://arxiv.org/abs/1811.11072
#
# Platform           : R
# Version            : v3.3.0
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# State-space Model Notation
#	y_t = x*beta + alpha_t + eps_t               eps_t ~ N([0,0], Sigma.BP)
#	alpha_t = rho*alpha_{t-1} + phi*c_t + nu_t           nu_t ~ N([0,0], Sigma.nu)
#   alpha_1 ~ N([0,0], Sigma.0)
#
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
# theta = (rho, phi, Sigma.nu, Sigma.BP, Sigma.0, beta)
# alpha = (alpha[1], ..., alpha[T])
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


suppressMessages(library(mvtnorm))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# AR(1) correlation matrix
#  - square matrix with rho^(|i-j|) in i, jth location
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
autocorr.mat <- function(p = 100, rho = 0.9) {
    mat <- diag(p)
    return(rho^abs(row(mat)-col(mat)))
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Calculate Covariance matrix for marginal of y_1, ..., y_T | theta
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

Sigma_T = function(T, rho, phi, sig.nu, sig.BP, sig0, beta){
	alpha_cov = autocorr.mat(T, rho)
	v_t = rep(NA, T)

	v_t[1] = sig0^2
	if(T>1){
		for(t in 2:T){
			v_t[t] = sum(rho^(2*(t-t:2)))*sig.nu^2 + rho^(2*(t-1)) * sig0^2
		}
	}
	v_mat = alpha_cov*NA
	for(i in 1:T){
		for(j in 1:T){
			if(j >= i){
				v_mat[i, j] = v_t[i]
				v_mat[j, i] = v_t[i]
			}
		}
	}
	Sigma = sig.BP^2*diag(1, T) + (v_mat*alpha_cov)

	Sigma
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Calculate Covariance matrix for marginal of alpha_1, ..., alpha_T | theta
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

Sigma_alpha = function(T, rho, phi, sig.nu, sig.BP, sig0, beta){
	alpha_cov = autocorr.mat(T, rho)
	v_t = rep(NA, T)

	v_t[1] = sig0^2
	if(T>1){
		for(t in 2:T){
			v_t[t] = sum(rho^(2*(t-t:2)))*sig.nu^2 + rho^(2*(t-1)) * sig0^2
		}
	}
	v_mat = alpha_cov*NA
	for(i in 1:T){
		for(j in 1:T){
			if(j >= i){
				v_mat[i, j] = v_t[i]
				v_mat[j, i] = v_t[i]
			}
		}
	}
	Sigma =  (v_mat*alpha_cov)

	Sigma
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Calculate Mean vector for marginal of (y_1, ..., y_T | theta)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

mu_T = function(T, x, c_t, rho, phi, sig.nu, sig.BP, sig0, beta){
	alpha_mu = rep(NA, T)
	alpha_mu[1] = 0
	if(T>1){
		for(t in 2:T){
			alpha_mu[t] = phi*sum(c_t[t:2,1]*rho^(t-t:2))
		}
	}
	mu = sum(x * beta) + alpha_mu
	mu
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Sample the unobserved states given a patient and a parameter vector
# 	bivariate version
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

ralpha_post = function(dat, theta, is.data = FALSE){

	n = length(dat)

	rho <- theta[1]
	phi <- theta[2]
	sig.nu <- theta[3]
	sig.BP <- theta[4]
	sig0 <- theta[5]
	beta <-  theta[-c(1:5)]

	# calculate covariance matrix
	max.T = max(sapply(dat, function(d) d$T))
	Sigma_all = Sigma_alpha(max.T, rho, phi, sig.nu, sig.BP, sig0, beta)
	# calculate mean function by observation
	mu_all = lapply(dat, function(d){
		x = d$x
		c_t = d$c_t
		mu_T(d$T, x, c_t, rho, phi, sig.nu, sig.BP, sig0, beta)
	})

	alphas = vector('list', n)
	for(i in 1:n){
		d = dat[[i]]
		# subset data to observed outcomes: data comes subset already, so no need to subset
		if(!is.data) y = d$y[d$I.y]
		if(is.data) y = d$y
		mu = mu_all[[i]][d$I.y]
		Sigma = as.matrix(Sigma_all[d$I.y, d$I.y])

		mu1 = mu + Sigma %*% solve(sig.BP^2*diag(d$T) + Sigma) %*% t(t(y - mu))
		Sigma1 = Sigma - Sigma %*% solve(sig.BP^2*diag(d$T) + Sigma) %*% Sigma
		# calculate negative log likelihood
		alphas[[i]] = rmvnorm(1, mean = mu1, sigma = Sigma1)
	}

	alphas
}


