# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# File Name          : plots.R
# Programmer Name    : Kristen Hunter
#                      kristenbhunter@gmail.com
#
# Last Updated       : Jan 2020
#
# Purpose            : Code to produce plots in paper
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

library(wesanderson)

# standardized colors for all plots
plot.colors = c(wes_palette('Darjeeling2')[2], wes_palette('Darjeeling2')[3])

#################################
# standard legend for all plots
#################################

standard.legend = function(point.data)
{
  if(all(point.data$cover))
  {
    legend.color = scale_color_manual(
      name = 'Cover',
      values = c('TRUE' = plot.colors[1]),
      labels = c('True')
    )
  } else if(all(!point.data$cover))
  {
    legend.color = scale_color_manual(
      name = 'Cover',
      values = c('TRUE' = plot.colors[2]),
      labels = c('False')
    )
  } else
  {
    legend.color = scale_color_manual(
      name = 'Cover',
      values = c('FALSE' = plot.colors[2], 'TRUE' = plot.colors[1]),
      labels = c('False', 'True')
    )
  }
  legend.shape = scale_shape_manual(
    name = 'Average\nadherence',
    values = c('padhere' = 16, 'prior.mean' = 15),
    labels = c('padhere' = 'True Value', 'prior.mean' = 'Prior Mean')
  )

  return(list(legend.color, legend.shape))

}

#################################
# create interval plot
#################################

gen.interval.plot = function(coverage.data, point.data, title, square.dots = FALSE)
{
  # coverage.data = plot.data[['coverage.data']]; point.data = plot.data[['point.data']];

  legend = standard.legend(point.data)

  if(square.dots)
  {
    interval.plot = ggplot(coverage.data, aes(x = value, y = factor(id), color = cover)) +
      geom_line(size = 2.5) +
      geom_point(data = point.data, aes(x = value, shape = factor(variable)), size = 5) +
      legend[[1]] + legend[[2]] +
      ylab('Patient') +
      xlab('Predicted percent adherence') +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      ggtitle(title) +
      theme(text = element_text(size = 20)) +
      theme(plot.title = element_text(size = 20)) +
      xlim(0, 1)
  } else
  {
    interval.plot = ggplot(coverage.data, aes(x = value, y = factor(id), color = cover)) +
      geom_line(size = 2.5) +
      geom_point(data = point.data, aes(x = value), size = 5) +
      legend[[1]] +
      ylab('Patient') +
      xlab('Predicted percent adherence') +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      ggtitle(title) +
      theme(text = element_text(size = 20)) +
      theme(plot.title = element_text(size = 20)) +
      xlim(0, 1)
  }
  return(interval.plot)
}

#################################
# create plot data
#################################

# take in expected coverage
create.plot.data = function(ad.quantiles, exp.coverage, type = 'bp', square.dots = FALSE)
{
  # type = 'bp'
  # type = 'ad'

  # rearrange data
  cover.text = round(100 * exp.coverage)
  if(type == 'ad')
  {
    cover.text = paste('prior', cover.text, sep = '.')
  }
  cover.var = paste('cover', cover.text, sep = '.')
  width.var = paste('width', cover.text, sep = '.')
  vars = c(
    paste('lower', cover.text, sep = '.'),
    paste('upper', cover.text, sep = '.')
  )
  id.vars = c('id', 'bpday', 'adday', 'true', cover.var, width.var)
  ad.quantiles$true = ad.quantiles$padhere

  # standardize colnames
  colnames = c('id', 'bpday', 'adday', 'true', 'cover', 'width', 'variable', 'value')

  # coverage data
  coverage.data = melt(ad.quantiles[ ,c(id.vars, vars)], id.vars)
  colnames(coverage.data) = colnames

  if(!square.dots)
  {
    # point data
    point.data = melt(ad.quantiles[ ,c(id.vars, 'padhere')], id.vars)
    colnames(point.data) = colnames
  } else
  {
    # point data
    point.data = melt(ad.quantiles[ ,c(id.vars, 'prior.mean', 'padhere')], id.vars)
    colnames(point.data) = colnames
  }

  return(list(point.data = point.data, coverage.data = coverage.data))
}

#################################
# generate coverage text
#################################

gen.coverage.text = function(ad.quantiles, cover.var, exp.coverage){
  return(paste(
    'Coverage: ',
    round(100 * mean(ad.quantiles[ , cover.var])),
    '%, Expected: ',
    round(100 * exp.coverage),
    '%',
    sep = ''
  ))
}

#################################
# coverage plot
#################################

plot.coverage = function(ad.quantiles,
                         sim.params,
                         exp.coverage = 0.95,
                         sample.size = 25,
                         type = 'bp',
                         square.dots = FALSE)
{
  # sample.size = 25; exp.coverage = 0.95; type = 'bp';

  # first calculate coverage and generate title
  cover.text = round(100 * exp.coverage)
  if(type == 'ad')
  {
    cover.text = paste('prior', cover.text, sep = '.')
  }
  cover.var = paste('cover', cover.text, sep = '.')
  coverage.text = gen.coverage.text(ad.quantiles, cover.var, exp.coverage)
  title = paste('Predicted average adherence\n', coverage.text, sep = '')

  # sample down to only show some patients on the final plot
  sample.size.min = min(sample.size, nrow(ad.quantiles))
  set.seed(3333)
  id.sample = sample(unique(ad.quantiles$id), sample.size.min)
  ad.quantiles = ad.quantiles[ad.quantiles$id %in% id.sample,]

  plot.data = create.plot.data(ad.quantiles, exp.coverage, type = type,
                               square.dots = square.dots)

  interval.plot = gen.interval.plot(
    coverage.data = plot.data[['coverage.data']],
    point.data = plot.data[['point.data']],
    title = title,
    square.dots = square.dots
  )

  return(interval.plot)
}


#################################
# generate a table summarizing theta
#################################

summarize.theta = function(theta)
{
  theta.plot.summary = rbind(
    apply(theta.plot.data, 2, mean),
    apply(theta.plot.data, 2, quantile, probs = 0.025),
    apply(theta.plot.data, 2, quantile, probs = 0.975)
  )
  rownames(theta.plot.summary) = c('Mean', 'Lower 95% Quantile', 'Upper 95% Quantile')
  theta.plot.summary = apply(theta.plot.summary, 1, round, digits = 2)
  kable(theta.plot.summary, format = 'latex')
}

