# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : state_space_functions.R
# Programmer Name    : Luis Campos
#                     soyluiscampos@gmail.com
#
# Purpose            : Supporting functions for analysis
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
# This function arranges a dataset (list form) into the format needed
#  to analyze data in STAN.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

arrange.SS_stan = function(dat, is.sim = TRUE){
	dat = lapply(dat, function(x) {
		x$ni = length(x$I.y);
		if(is.sim){
			x$y1 = x$y1[x$I.y]
			x$y2 = x$y2[x$I.y]
		}
		x
	})

	X = do.call('rbind', lapply(dat, function(x) x$x))
	y = do.call('c', lapply(dat, function(x) c(x$y1, x$y2)))
	yid = do.call('c', lapply(dat, function(x) c(rep(1, length(x$y1)), rep(2, length(x$y2)))))

	N = length(dat)
	J = ncol(X)

	T = do.call('c', lapply(dat, function(x) x$T))
	maxT = max(T)
	ni = do.call('c', lapply(dat, function(x) x$ni))
	n = sum(ni)

	c_matrix = matrix(NA, nrow = N, ncol = maxT)
	for(i in 1:N){
		c_t = dat[[i]]$c_t[,1]
		c_matrix[i, 1:length(c_t)] = c_t
		if(length(c_t) < maxT){
			c_matrix[i, (length(c_t)+1):maxT] = 9999
		}
	}

	id = rep(1:N, 2*ni)
	t = do.call('c', lapply(dat, function(x) c(x$I.y, x$I.y)))

	adherence_dat = list(N = N, J = J, T = T, maxT = maxT, ni = ni, n = n,
	                     X = X, c = c_matrix, y = y, id = id, t = t, yid = yid)

	return(adherence_dat)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# general plotting function takes name, parameters and a label
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

plot_densities = function(fit.extracts, which_par, which_par_id, par, lab, xlim = NULL, ylim = NULL, make.log = FALSE, lwd = 1){
	if(make.log){
		densities = lapply(fit.extracts, function(x){density(log(x[[which_par]][,which_par_id]))})
		par = log(par)
		lab = paste('log', lab)
	}else{
		densities = lapply(fit.extracts, function(x){density(x[[which_par]][,which_par_id])})
	}
	if(is.null(xlim)) xlim = range(c(sapply(densities, function(x) x$x), par), na.rm = TRUE)
	if(is.null(ylim)) ylim = range(sapply(densities, function(x) x$y))


	plot(NA, pch = 'n', xlim = xlim, ylim = ylim, ylab = 'Density', xlab = lab)
	tmp = lapply(densities, function(x) lines(x, lwd = lwd))
	abline(v = par, col ='red', lwd =2)
}


get_lims = function(fit.extracts, which_par, which_par_id, par, lab, make.log = FALSE){
	if(make.log){
		densities = lapply(fit.extracts, function(x){density(log(x[[which_par]][,which_par_id]))})
		par = log(par)
		lab = paste('log', lab)
	}else{
		densities = lapply(fit.extracts, function(x){density(x[[which_par]][,which_par_id])})
	}
	xlim = range(c(sapply(densities, function(x) x$x), par), na.rm = TRUE)
	ylim = range(sapply(densities, function(x) x$y))


	list(xlim = xlim, ylim = ylim)

}




plot_densities_data = function(fit.extracts, which_par, which_par_id, lab, make.log = FALSE){
	if(make.log){
		densities = lapply(fit.extracts, function(x){density(log(x[[which_par]][,which_par_id]))})
		par = sapply(fit.extracts, function(x){mean((x[[which_par]][,which_par_id]))})
		par = log(par)
		lab = paste('log', lab)
	}else{
		densities = lapply(fit.extracts, function(x){density(x[[which_par]][,which_par_id])})
		par = sapply(fit.extracts, function(x){mean((x[[which_par]][,which_par_id]))})
	}
	xlim = range(c(sapply(densities, function(x) x$x), par), na.rm = TRUE)
	ylim = range(sapply(densities, function(x) x$y))


	plot(NA, pch = 'n', xlim = xlim, ylim = ylim, ylab = 'Density', xlab = lab)
	tmp = lapply(densities, function(x) lines(x))
	abline(v = par, col ='red', lwd =2)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Traceplot
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
plot_trace = function(trace, par, lab, points = FALSE){
	if(is.null(par)) par = mean(trace)
	if(points){
		plot(trace, pch = 19, cex = 0.1, ylim = range(c(trace, par)), xlab = 'sample', ylab = lab, col = rgb(0.1, 0.1, 0.1, alpha = 0.5))
	}else{
		plot(trace, type = 'l', ylim = range(c(trace, par)), xlab = 'sample', ylab = lab)
	}
	abline(h = par, col ='red', lwd = 2)
}
