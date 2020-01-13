# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : functions_pf.R
# Programmer Name    : Kristen Hunter
#                      kristenbhunter@gmail.com
# 
# Last Updated       : Jan 2020
#
# Purpose            : Functions to summarize output and produce plots in paper
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

################################################
# CHANGE THIS to output directory to be summarized
################################################
run.date = '20200110_16'

# setup
base.dir = '/Users/khunter/Dropbox/Medication-Adherence/'
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
  
  # save out data
  ad.quantile.file = paste(run.dir, 'ad_quantiles.csv', sep = '')
  write.csv(ad.quantiles, file = ad.quantile.file)
}

# get params
load(paste(run.dir, 'run_params_', run.date, '.Rdata', sep = ''))

# source important functions
source(paste(predict.dir, 'functions.R', sep = ''))
source(paste(predict.dir, 'plots.R', sep = ''))
source(paste(predict.dir, 'setup.R', sep = ''))

# read and concantenate in all the data
ad.means = NULL
ad.prior = NULL
for(i in 1:run.params[['imputations']])
{
  draws = readRDS(file = paste(run.dir, 'draws_', i, '.rds', sep = ''))
  ad.means = cbind(ad.means, draws[['ad.means']])
  ad.prior = cbind(ad.prior, draws[['ad.prior']])
}
summarize.adherence(ad.means, ad.prior, run.params)
