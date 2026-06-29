set.seed(12345)
par(mfrow=c(1, 1), mar = c(4,4,2,2))

nsamples <- length(samples[['theta1']])




nonmasting <- c()
masting <- c()
for(t in 1:data$N_trees){
  idxs <- seq(1+(t-1)*data$N_max_years, t*data$N_max_years, 1)
  for(i in idxs){
    masting <-  c(masting, na.omit(as.numeric(ifelse(samples[[paste0('states_pred[', i,']')]] == 2,samples[[paste0('seed_counts_pred[', i,']')]],NA))))
    nonmasting <-  c(nonmasting, na.omit(as.numeric(ifelse(samples[[paste0('states_pred[', i,']')]] == 1,samples[[paste0('seed_counts_pred[', i,']')]],NA))))
  }
}
masting <- as.numeric(masting)
nonmasting <- as.numeric(nonmasting)

pdf(file = file.path(figpath, 'model_overview.pdf'),
    width = 10, height = 6)
layout( matrix(c(1,1,1,4,4,4,2,2,3,3,5,5), ncol=2, byrow = F))
# layout( matrix(c(2,2,3,3,5,5,1,1,1,4,4,4), ncol=2, byrow = F))
par(mar = c(3,2,2,2))


maxx <- 3000
hist(nonmasting, breaks=seq(0,maxx,l=600), col = '#48754860', border = 'white', 
     xlab = '', main = '', prob = TRUE, yaxt = "n", ylab = '', xaxt = 'n', xlim = c(0,400))

usr <- par("usr")
segments(x0 = 1,  x1 = 500, y0 = usr[3], y1 = usr[3], lwd = 1.5)
axis(1, at=seq(0, 500,100), labels=seq(0, 500,100), 
     lwd = 0, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))   
title(xlab = 'Number of seeds', line = 2)

hist(masting, breaks=seq(0,maxx,l=600), col = "#8F272760", border = 'white', 
     xlab = 'Number of seeds', main = '', add = TRUE, prob = TRUE, xlim = c(0,400))
text(x = 25, y = 0.03, labels = 'Low reproductive \nstate', col = '#487548', font = 1, adj = 0, cex = 1.2)
text(x = 100, y = 0.01, labels = 'High reproductive state', col = util$c_mid_highlight, font = 1, adj = 0, cex = 1.2)

mtext(LETTERS[1], side = 3, line = 0.5, adj = -0.05, cex = 1.2, font = 2, col = 'grey30')


# Two modalities
# nonmasting <- c()
# masting <- c()
# for(i in 1:nsamples){
#   theta <- samples[['theta1']][i]
#   nonmas <- sapply(1:100, function(x){
#     ifelse(rbinom(1, 1, theta) == 1, 0, MASS::rnegbin(1, samples[['lambda1']][i], 1/samples[['psi1']][i]))})
#   nonmasting <- c(nonmasting, nonmas)
#   mas <- MASS::rnegbin(100, samples[['lambda20']][i], 1/samples[['psi2']][i])
#   masting <- c(masting, mas)
# }
# 
# hist(nonmasting, breaks=seq(0,1070,l=80), col = '#48754860', border = 'white', 
#      xlab = '', main = '', prob = TRUE, yaxt = "n", ylab = '', xaxt = 'n', xlim = c(0,500))
# 
# usr <- par("usr")
# segments(x0 = 1,  x1 = 500, y0 = usr[3], y1 = usr[3], lwd = 1.5)
# axis(1, at=seq(0, 500,100), labels=seq(0, 500,100), 
#      lwd = 0, lwd.ticks = 1, tck = -0.03, 
#      cex.axis = 0.9, mgp = c(2, 0.7, 0))   
# title(xlab = 'Number of seeds', line = 2)
# 
# masting <- MASS::rnegbin(10000, 160, 2)
# hist(masting, breaks=seq(0,1070,l=80), col = "#8F272760", border = 'white', 
#      xlab = 'Number of seeds', main = '', add = TRUE, prob = TRUE, xlim = c(0,500))
# text(x = 25, y = 0.03, labels = 'Low-reproduction\nstate', col = '#487548', font = 1, adj = 0)
# text(x = 100, y = 0.008, labels = 'High-reproduction state', col = util$c_mid_highlight, font = 1, adj = 0)
# 
# mtext(LETTERS[1], side = 3, line = 0.5, adj = -0.05, cex = 1.2, font = 2, col = 'grey30')

# Tree-level alternation in states
par(mar = c(3,3.5,2,2))
t <- 58

idxs_tree <-(1+N_max_years*(t-1)):(N_max_years*t)
observed_idxs_tree <- tree_start_idxs[t]:tree_end_idxs[t]
observed_flags <- years[tree_start_idxs[t]:tree_end_idxs[t]]

names <- sapply(idxs_tree,
                function(n) paste0('states_pred[', n, ']'))

ylab="Latent state"
display_ylim=c(0, 1)
main=''

# Construct bins
N <- length(names)
bin_min <- 0.5
bin_max <- N + 0.5
bin_delta <- 1
breaks <- seq(bin_min, bin_max, bin_delta)

plot_config <- util$configure_bin_plotting(breaks)
plot_idxs <- plot_config[[1]]
plot_xs <- plot_config[[2]]

# Construct marginal quantiles
# probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[names[n]]]-1, probs)
}
quantiles <- sapply(1:N, calc)
plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                        function(n) quantiles[1:5, n]))

delta <- 0.05 * (display_ylim[2] - display_ylim[1])
display_ylim[1] <- display_ylim[1] - delta  
display_ylim[2] <- display_ylim[2] + delta      

plot(1, type="n", main=main,
     xlim=c(bin_min, bin_max), xlab='', xaxt="n",
     ylim=display_ylim, ylab='', yaxt="n",, frame = FALSE)

usr <- par("usr")
segments(x0 = 1,  x1 = 43, y0 = usr[3], y1 = usr[3], lwd = 1.5)
axis(1, at=seq(1,N,10), labels=seq(1980,2020,10), 
     lwd = 0, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))    

axis(2, at=seq(0,1,1), c('Low', 'High'), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))  
title(ylab = ylab, line = 2.5, cex.lab = 1.1)

for(n in 1:N){
  
  idplot <- c(n*2-1, n*2)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[1,idplot], rev(plot_quantiles[5,idplot])),
          col = ifelse(n%in%observed_flags, "#c8d8e5", 'grey90'), border = NA)
  # polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
  #         c(plot_quantiles[2,idplot], rev(plot_quantiles[8,idplot])),
  #         col = ifelse(n%in%observed_flags, util$c_light_highlight, 'grey80'), border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[2,idplot], rev(plot_quantiles[4,idplot])),
          col = ifelse(n%in%observed_flags, "#a4bed5", 'grey80'), border = NA)
  # polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
  #         c(plot_quantiles[4,idplot], rev(plot_quantiles[6,idplot])),
  #         col = ifelse(n%in%observed_flags, util$c_mid_highlight, 'grey60'), border = NA)
  
  lines(plot_xs[idplot], plot_quantiles[3, idplot],
        col= ifelse(n%in%observed_flags,"#728595", 'grey60'), lwd=2)
}


mtext(LETTERS[3], side = 3, line = 0.5, adj = -0.1, cex = 1.2, font = 2, col = 'grey30')

# Tree-level alternation in seeds
par(mar = c(3,3.5,0,2))

idxs_tree <-(1+N_max_years*(t-1)):(N_max_years*t)
observed_idxs_tree <- tree_start_idxs[t]:tree_end_idxs[t]
observed_flags <- years[tree_start_idxs[t]:tree_end_idxs[t]]

names <- sapply(idxs_tree,
                function(n) paste0('seed_counts_pred[', n, ']'))

ylab="Number of seeds"
display_ylim=c(0, 400)
main=''

# Construct bins
N <- length(names)
bin_min <- 0.5
bin_max <- N + 0.5
bin_delta <- 1
breaks <- seq(bin_min, bin_max, bin_delta)

plot_config <- util$configure_bin_plotting(breaks)
plot_idxs <- plot_config[[1]]
plot_xs <- plot_config[[2]]

# Construct marginal quantiles
# probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[names[n]]]-1, probs)
}
quantiles <- sapply(1:N, calc)
plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                        function(n) quantiles[1:5, n]))

delta <- 0.05 * (display_ylim[2] - display_ylim[1])
display_ylim[1] <- display_ylim[1] - delta  
display_ylim[2] <- display_ylim[2] + delta      

plot(1, type="n", main=main,
     xlim=c(bin_min, bin_max), xlab='', xaxt="n",
     ylim=display_ylim, ylab='', yaxt="n",, frame = FALSE)

usr <- par("usr")
segments(x0 = 1,  x1 = 43, y0 = usr[3], y1 = usr[3], lwd = 1.5)
axis(1, at=seq(1,N,10), labels=seq(1980,2020,10), 
     lwd = 0, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))    

axis(2, at=seq(0,400,100), seq(0,400,100), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))  
title(ylab = ylab, line = 2.5, cex.lab = 1.1)

for(n in 1:N){
  
  idplot <- c(n*2-1, n*2)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[1,idplot], rev(plot_quantiles[5,idplot])),
          col = ifelse(n%in%observed_flags, "#c8d8e5", 'grey90'), border = NA)
  # polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
  #         c(plot_quantiles[2,idplot], rev(plot_quantiles[8,idplot])),
  #         col = ifelse(n%in%observed_flags, util$c_light_highlight, 'grey80'), border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[2,idplot], rev(plot_quantiles[4,idplot])),
          col = ifelse(n%in%observed_flags, "#a4bed5", 'grey80'), border = NA)
  # polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
  #         c(plot_quantiles[4,idplot], rev(plot_quantiles[6,idplot])),
  #         col = ifelse(n%in%observed_flags, util$c_mid_highlight, 'grey60'), border = NA)
  
  lines(plot_xs[idplot], plot_quantiles[3, idplot],
        col= ifelse(n%in%observed_flags,"#728595", 'grey60'), lwd=2)
}

for(i in observed_idxs_tree) {
  lines(c(years[i] - 0.5, years[i] + 0.5),
        rep(seed_counts[i], 2),
        col="white", lwd=4)
  lines(c(years[i] - 0.5, years[i] + 0.5),
        rep(data$seed_counts[i], 2),
        col="black", lwd=2)
}  
# mtext(LETTERS[2], side = 3, line = 0.5, adj = -0.15, cex = 1.2, font = 2, col = 'grey30')

# Probabilities

meansummertemp <- mean(clim_data$meantmax_ja[clim_data$year %in% c(1980:2022)])
samples_transitiontohigh <- invlogit(logit(samples[['tau_nm_m0']]) + samples[['beta_nm_m']] * (meansummertemp-summertemp0))
samples_persistenceinlow <- 1 - samples_transitiontohigh
samples_persistenceinhigh <- samples[['tau_m_m0']]
samples_transitiontolow <- 1-samples_persistenceinhigh

average_transitiontohigh <- util$ensemble_mcmc_quantile_est(samples_transitiontohigh, c(0.05, 0.5, 0.95))
average_persistenceinlow <- util$ensemble_mcmc_quantile_est(samples_persistenceinlow, c(0.05, 0.5, 0.95))
prob_ratio <- util$ensemble_mcmc_quantile_est(samples_transitiontohigh/samples_persistenceinlow, c(0.05, 0.5, 0.95))
average_transitiontolow <- util$ensemble_mcmc_quantile_est(1-samples[['tau_m_m0']], c(0.05, 0.5, 0.95))
average_persistenceinhigh <- util$ensemble_mcmc_quantile_est(samples[['tau_m_m0']], c(0.05, 0.5, 0.95))


par(mar = c(1.5,2,2,2))

plot(x = NULL, y = NULL,
     xlim = c(-0.9, 2.2 + 0.8),
     ylim = c(-0.2, 2.2 + 0.5),
     yaxt = 'n', xaxt = 'n',
     xlab = '',
     ylab = '',
     main = '',
     frame.plot	= FALSE)
mtext(LETTERS[2], side = 3, line = 0.5, adj = -0.05, cex = 1.2, font = 2, col = 'grey30')
param <- hist(samples_persistenceinlow, breaks = seq(0,1, 0.025), plot = FALSE)
param_counts <- param$counts / max(param$counts)*0.9
rect(ybottom = param$breaks[1:(length(param$breaks)-1)] + 1.2,
     xleft = 0,
     ytop = param$breaks[2:(length(param$breaks))] + 1.2,
     xright =  param_counts,
     col = '#b8d2b8',
     #border = ifelse(i == 4, '#bfe3d8', ifelse(i == 10, '#bfdce3', "grey90"))
     border = '#b8d2b8')

param <- hist(samples_transitiontohigh, breaks = seq(0,1, 0.025), plot = FALSE)
param_counts <- param$counts / max(param$counts)*0.9
rect(ybottom = param$breaks[1:(length(param$breaks)-1)] + 1.2,
     xleft = 1.2,
     ytop = param$breaks[2:(length(param$breaks))] + 1.2, 
     xright =  param_counts + 1.2,
     col = '#e4abab',
     #border = ifelse(i == 4, '#bfe3d8', ifelse(i == 10, '#bfdce3', "grey90"))
     border = '#e4abab')

param <- hist(samples_transitiontolow, breaks = seq(0,1, 0.025), plot = FALSE)
param_counts <- param$counts / max(param$counts)*0.9
rect(ybottom = param$breaks[1:(length(param$breaks)-1)],
     xleft = 0,
     ytop = param$breaks[2:(length(param$breaks))],
     xright =  param_counts,
     col = '#b8d2b8',
     #border = ifelse(i == 4, '#bfe3d8', ifelse(i == 10, '#bfdce3', "grey90"))
     border = '#b8d2b8')


param <- hist(samples_persistenceinhigh, breaks = seq(0,1, 0.025), plot = FALSE)
param_counts <- param$counts / max(param$counts)*0.9
rect(ybottom = param$breaks[1:(length(param$breaks)-1)],
     xleft = 1.2,
     ytop = param$breaks[2:(length(param$breaks))],
     xright =  param_counts + 1.2,
     col = '#e4abab',
     #border = ifelse(i == 4, '#bfe3d8', ifelse(i == 10, '#bfdce3', "grey90"))
     border = '#e4abab')

segments(x0 = 0, y0 = 0, y1 = 1)
segments(x0 = 0, x1 = -0.05, y0 = 0)
segments(x0 = 0, x1 = -0.05, y0 = 0.5)
segments(x0 = 0, x1 = -0.05, y0 = 1)
text(x = -0.08, y = c(0, .5, 1), labels = c(0, .5, 1), cex = 0.9, adj = 1)

segments(x0 = 0, y0 = 1.2, y1 = 2.2)
segments(x0 = 0, x1 = -0.05, y0 = 1.2)
segments(x0 = 0, x1 = -0.05, y0 = 1.2 + 0.5)
segments(x0 = 0, x1 = -0.05, y0 = 1.2 + 1)
text(x = -0.08, y = 1.2+ c(0, .5, 1), labels = c(0, .5, 1), cex = 0.9, adj = 1)

segments(x0 = 1.2, y0 = 0, y1 = 1)
segments(x0 = 1.2, x1 = 1.2-0.05, y0 = 0)
segments(x0 = 1.2, x1 = 1.2-0.05, y0 = 0.5)
segments(x0 = 1.2, x1 = 1.2-0.05, y0 = 1)
text(x = 1.2-0.08, y = c(0, .5, 1), labels = c(0, .5, 1), cex = 0.9, adj = 1)

segments(x0 = 1.2, y0 = 1.2, y1 = 2.2)
segments(x0 = 1.2, x1 = 1.2-0.05, y0 = 1.2)
segments(x0 = 1.2, x1 = 1.2-0.05, y0 = 1.2 + 0.5)
segments(x0 = 1.2, x1 = 1.2-0.05, y0 = 1.2 + 1)
text(x = 1.2-0.08, y = 1.2+ c(0, .5, 1), labels = c(0, .5, 1), cex = 0.9, adj = 1)

segments(x0 = -0.3, y0 = -0.25, y1 = 2.45, lwd = 1.7)
segments(x0 = -0.3, x1 = -0.2, y0 = -0.25, y1 = -0.25, lwd = 1.7)
segments(x0 = -0.3, x1 = -0.2, y0 = 2.45, y1 = 2.45, lwd = 1.7)
segments(x0 = 2.4, y0 = -0.25, y1 = 2.45, lwd = 1.7)
segments(x0 = 2.4, x1 = 2.3, y0 = -0.25, y1 = -0.25, lwd = 1.7)
segments(x0 = 2.4, x1 = 2.3, y0 = 2.45, y1 = 2.45, lwd = 1.7)

text(x = -0.55, y = c(.5, 0.5+1.2), labels = c('High', 'Low'), cex = 1.4, adj = 0.5, col = c('#A25050', '#487548'), srt = 90)
text(x = c(0.5+1.2, 0.5), y = 2.6, labels = c('High', 'Low'), cex = 1.4, adj = 0.5, col = c('#A25050', '#487548'))
text(x = 2.6, y = 1.15, labels = 'Transition matrix', adj = 0.5, srt = 270, cex = 1.4)

# Add the probability values...
text(x = 0.5, y = 0.6, 
     labels = paste0(round(average_transitiontolow['50%']*100,0), '% (',round(average_transitiontolow['5%']*100,0), ' \u2013 ', round(average_transitiontolow['95%']*100,0), '%)'), 
     cex = 1, adj = 0.5, col = c('#487548'))
text(x = 0.5+1.2, y = 0.4, 
     labels = paste0(round(average_persistenceinhigh['50%']*100,0), '% (',round(average_persistenceinhigh['5%']*100,0), ' \u2013 ', round(average_persistenceinhigh['95%']*100,0), '%)'), 
     cex = 1, adj = 0.5, col = c('#A25050'))

text(x = 0.5, y = 0.5+1.2, 
     labels = paste0(round(average_persistenceinlow['50%']*100,0), '% (',round(average_persistenceinlow['5%']*100,0), ' \u2013 ', round(average_persistenceinlow['95%']*100,0), '%)'), 
     cex = 1, adj = 0.5, col = c('#487548'))
text(x = 0.5+1.2, y = 0.47+1.2, 
     labels = paste0(round(average_transitiontohigh['50%']*100,0), '% (',round(average_transitiontohigh['5%']*100,0), ' \u2013 ', round(average_transitiontohigh['95%']*100,0), '%)'), 
     cex = 1, adj = 0.5, col = c('#A25050'))



# All trees
par(mar = c(3,3.5,2,2))
trees <- trees_tokeep
sum_states <- list()
for(y in 1:N_max_years){
  
  idxs_tree <- (y+N_max_years*(trees-1))
  names <- sapply(idxs_tree,
                  function(n) paste0('states_pred[', n, ']'))
  
  sum_states_y <- Reduce("+", samples[names])
  
  sum_states[[y]] <- (sum_states_y - length(idxs_tree))/length(idxs_tree)
  
  
}
names(sum_states) <- paste0('sum_states[', 1:N_max_years, ']')
names <- names(sum_states)

ylab="Masting trees (%)"
display_ylim=c(0, 100)
main=''

# Construct bins
N <- length(names)
bin_min <- 0.5
bin_max <- N + 0.5
bin_delta <- 1
breaks <- seq(bin_min, bin_max, bin_delta)

plot_config <- util$configure_bin_plotting(breaks)
plot_idxs <- plot_config[[1]]
plot_xs <- plot_config[[2]]

# Construct marginal quantiles
# probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(sum_states[[names[n]]]*100, probs)
}
quantiles <- sapply(1:N, calc)
plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                        function(n) quantiles[1:5, n]))

delta <- 0.05 * (display_ylim[2] - display_ylim[1])
display_ylim[1] <- display_ylim[1] - delta   
display_ylim[2] <- display_ylim[2] + delta      

plot(1, type="n", main=main,
     xlim=c(bin_min, bin_max), xlab='', xaxt="n", yaxt="n",
     ylim=display_ylim, ylab='', frame = FALSE)

usr <- par("usr")
segments(x0 = 1,  x1 = 43, y0 = usr[3], y1 = usr[3], lwd = 1.5)
axis(1, at=seq(1,N,10), labels=seq(1980,2020,10), 
     lwd = 0, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))    

axis(2, at=seq(0,100,25), seq(0,100,25), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))  
title(ylab = ylab, line = 2.5, cex.lab = 1.1)



for(n in 1:N){
  
  idplot <- c(n*2-1, n*2)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[1,idplot], rev(plot_quantiles[5,idplot])),
          col = "#c8d8e5", border = NA)
  # polygon(c(plot_xs[idplot], rev(plot_xs[idplot])), 
  #         c(plot_quantiles[2,idplot], rev(plot_quantiles[8,idplot])),
  #         col = util$c_light_highlight, border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[2,idplot], rev(plot_quantiles[4,idplot])),
          col = "#a4bed5" , border = NA)
  # polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
  #         c(plot_quantiles[4,idplot], rev(plot_quantiles[6,idplot])),
  #         col = util$c_mid_highlight, border = NA)
  
  lines(plot_xs[idplot], plot_quantiles[3, idplot],
        col= "#728595" , lwd=2)
}

mtext(LETTERS[4], side = 3, line = 0.5, adj = -0.1, cex = 1.2, font = 2, col = 'grey30')

dev.off()

