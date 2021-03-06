# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : functions_pf.R
# Programmer Name    : Kristen Hunter
#                      kristenbhunter@gmail.com
# 
# Last Updated       : Jan 2021
#
# Purpose            : Functions to summarize output and produce plots in paper
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

################################################
# CHANGE THIS to output directory to be summarized
################################################
run.date = '20210101_16'

# setup
base.dir = '/Users/khunter/Dropbox/archive/Medication-Adherence/'
predict.dir = paste(base.dir, 'predict/', sep = '')
run.dir = paste(predict.dir, 'out/run_', run.date, '/', sep = '')

################################################
# plots and save out summary statistics
################################################

summarize.adherence = function(ad.means, ad.prior, run.params)
{
  # get test data once
  setup.output = data.setup(base.dir)
  test.melt = setup.output[['datasets']][['test.melt']]
  
  # get adherence quantiles
  ad.quantiles = get.quantiles(ad.means, ad.prior, test.melt)
  
  # 95% intervals
  interval.plot.95 = plot.coverage(ad.quantiles, run.params, exp.coverage = 0.95)
  interval.plot.file.95 = paste(run.dir, 'coverage_95.png', sep = '')
  png(interval.plot.file.95, width = 1000, height = 1000)
  print(interval.plot.95)
  dev.off()
  
  # 80% intervals
  interval.plot.80 = plot.coverage(ad.quantiles, run.params, exp.coverage = 0.80)
  interval.plot.file.80 = paste(run.dir, 'coverage_80.png', sep = '')
  png(interval.plot.file.80, width = 1000, height = 1000)
  print(interval.plot.80)
  dev.off()
  
  # 50% intervals
  interval.plot.50 = plot.coverage(ad.quantiles, run.params, exp.coverage = 0.5)
  interval.plot.file.50 = paste(run.dir, 'coverage_50.png', sep = '')
  png(interval.plot.file.50, width = 1000, height = 1000)
  print(interval.plot.50)
  dev.off()
  
  # 95% ad intervals
  interval.plot.ad.95 = plot.coverage(ad.quantiles, run.params, exp.coverage = 0.95, type = 'ad')
  interval.plot.file.ad.95 = paste(run.dir, 'coverage_ad_95.png', sep = '')
  png(interval.plot.file.ad.95, width = 1000, height = 1000)
  print(interval.plot.ad.95)
  dev.off()
  
  # 80% ad intervals
  interval.plot.ad.80 = plot.coverage(ad.quantiles, run.params, exp.coverage = 0.8, type = 'ad')
  interval.plot.file.ad.80 = paste(run.dir, 'coverage_ad_80.png', sep = '')
  png(interval.plot.file.ad.80, width = 1000, height = 1000)
  print(interval.plot.ad.80)
  dev.off()
  
  # 50% ad intervals
  interval.plot.ad.50 = plot.coverage(ad.quantiles, run.params, exp.coverage = 0.5, type = 'ad')
  interval.plot.file.ad.50 = paste(run.dir, 'coverage_ad_50.png', sep = '')
  png(interval.plot.file.ad.50, width = 1000, height = 1000)
  print(interval.plot.ad.50)
  dev.off()
  
  # summarize performance
  interval.summary = data.frame(
    expected.coverage = c('95', '80', '50'),
    actual.coverage = c(mean(ad.quantiles$cover.95), mean(ad.quantiles$cover.80), mean(ad.quantiles$cover.50)),
    mean.width = c(mean(ad.quantiles$width.95), mean(ad.quantiles$width.80), mean(ad.quantiles$width.50)),
    max.width = c(max(ad.quantiles$width.95), max(ad.quantiles$width.80), max(ad.quantiles$width.50))
  )
  
  # save out data
  ad.quantile.file = paste(run.dir, 'ad_quantiles.csv', sep = '')
  write.csv(ad.quantiles, file = ad.quantile.file)
  
  return(interval.summary)
}

# get params
load(paste(run.dir, 'run_params_', run.date, '.Rdata', sep = ''))

# source important functions
source(paste(predict.dir, 'functions.R', sep = ''))
source(paste(predict.dir, 'plots.R', sep = ''))
source(paste(predict.dir, 'setup.R', sep = ''))

# read and concatenate in all the data
draws = readRDS(file = paste(run.dir, 'draws.rds', sep = ''))
ad.means = draws[['ad.means']]
ad.prior = draws[['ad.prior']]
interval.summary = summarize.adherence(ad.means, ad.prior, sim.params)
print(interval.summary)

# produce table summary of adherence model
summarize.theta.a = function(theta.a, param_titles)
{
  fixef.summary = cbind(
    apply(theta.a$fixef, 2, mean),
    apply(theta.a$fixef, 2, quantile, probs = 0.025),
    apply(theta.a$fixef, 2, quantile, probs = 0.975)
  )
  ranef.sd.summary = c(mean(theta.a$ranef.sd), quantile(theta.a$ranef.sd, 0.025), quantile(theta.a$ranef.sd, 0.975))
  alpha.summary = c(mean(theta.a$alpha.new), quantile(theta.a$alpha.new, 0.025), quantile(theta.a$alpha.new, 0.975))
  theta_summary = data.frame(rbind(fixef.summary, ranef.sd.summary, alpha.summary))
  colnames(theta_summary) = c('mean', 'quantile1', 'quantile2')
  theta_summary$parameter = rownames(theta_summary)
  
  if(!is.null(param_titles))
  {
    add.titles.query = "
    SELECT 
      param_titles.nice AS parameter, theta_summary.mean,
      theta_summary.quantile1, theta_summary.quantile2
    FROM
      theta_summary
    LEFT JOIN
      param_titles ON theta_summary.parameter = param_titles.orig
    "
    theta_summary = sqldf(add.titles.query)
  }
  return(theta_summary)
}

load(paste0(run.dir, 'theta_a.RData'))
# maps out original variable names to "nice" ones you may want to print in a table
param_titles = data.frame(
  orig = c(
    'intercept', 'gender', 'ranef.sd.summary', 'alpha.summary'
  ),
  nice = c(
    'Intercept', 'Male', 'Random effect standard deviation', 'Alpha'
  ))
theta.a.summary = summarize.theta.a(theta.a, param_titles)
print(theta.a.summary)
