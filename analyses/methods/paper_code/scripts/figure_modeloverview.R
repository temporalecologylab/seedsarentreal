set.seed(123456)
par(mfrow=c(1, 1), mar = c(4,4,2,2))

nsamples <- length(samples[['theta1']])

pdf(file = file.path(figpath, 'model_overview.pdf'),
    width = 9, height = 5)
layout( matrix(c(1,2,3,4), ncol=2, byrow = T))
par(mar = c(3,2,2,2))

# Two modalities
nonmasting <- c()
masting <- c()
for(i in 1:nsamples){
  theta <- samples[['theta1']][i]
  nonmas <- sapply(1:100, function(x){
    ifelse(rbinom(1, 1, theta) == 1, 0, MASS::rnegbin(1, samples[['lambda1']][i], 1/samples[['psi1']][i]))})
  nonmasting <- c(nonmasting, nonmas)
  mas <- MASS::rnegbin(100, samples[['lambda20']][i], 1/samples[['psi2']][i])
  masting <- c(masting, mas)
}

hist(nonmasting, breaks=seq(0,800,l=80), col = '#48757560', border = 'white', 
     xlab = '', main = '', prob = TRUE, yaxt = "n", ylab = '', xaxt = 'n', xlim = c(0,500))

usr <- par("usr")
segments(x0 = 1,  x1 = 500, y0 = usr[3], y1 = usr[3])
axis(1, at=seq(0, 500,100), labels=seq(0, 500,100), 
     lwd = 0, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))   
title(xlab = 'Number of seeds', line = 2)

masting <- MASS::rnegbin(10000, 160, 2)
hist(masting[masting < 650], breaks=seq(0,800,l=80), col = "#8F272760", border = 'white', 
     xlab = 'Number of seeds', main = '', add = TRUE, prob = TRUE, xlim = c(0,500))
text(x = 30, y = 0.03, labels = 'Low-reproduction\nstate', col = util$c_mid_teal, font = 1, adj = 0)
text(x = 100, y = 0.008, labels = 'High-reproduction state', col = util$c_mid, font = 1, adj = 0)

mtext(LETTERS[1], side = 3, line = 0.5, adj = -0.05, cex = 1.2, font = 2, col = 'grey30')


# Tree-level alternation
par(mar = c(3,3.5,2,2))
t <- 13

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
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[names[n]]], probs)
}
quantiles <- sapply(1:N, calc)
plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                        function(n) quantiles[1:9, n]))

delta <- 0.05 * (display_ylim[2] - display_ylim[1])
display_ylim[1] <- -10
display_ylim[2] <- display_ylim[2] + delta      

plot(1, type="n", main=main,
     xlim=c(bin_min, bin_max), xlab='', xaxt="n",
     ylim=display_ylim, ylab='', yaxt="n",, frame = FALSE)

usr <- par("usr")
segments(x0 = 1,  x1 = 43, y0 = usr[3], y1 = usr[3])
axis(1, at=seq(1,N,10), labels=seq(1980,2020,10), 
     lwd = 0, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))    

axis(2, at=seq(0,400,100), seq(0,400,100), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))  
title(ylab = ylab, line = 2.5)

for(n in 1:N){
  
  idplot <- c(n*2-1, n*2)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[1,idplot], rev(plot_quantiles[9,idplot])),
          col = ifelse(n%in%observed_flags, util$c_light, 'grey90'), border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[2,idplot], rev(plot_quantiles[8,idplot])),
          col = ifelse(n%in%observed_flags, util$c_light_highlight, 'grey80'), border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[3,idplot], rev(plot_quantiles[7,idplot])),
          col = ifelse(n%in%observed_flags, util$c_mid, 'grey70'), border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[4,idplot], rev(plot_quantiles[6,idplot])),
          col = ifelse(n%in%observed_flags, util$c_mid_highlight, 'grey60'), border = NA)
  
  lines(plot_xs[idplot], plot_quantiles[5, idplot],
        col= ifelse(n%in%observed_flags,util$c_dark, 'grey50'), lwd=2)
}

for(i in observed_idxs_tree) {
  lines(c(years[i] - 0.5, years[i] + 0.5),
        rep(seed_counts[i], 2),
        col="white", lwd=4)
  lines(c(years[i] - 0.5, years[i] + 0.5),
        rep(data$seed_counts[i], 2),
        col="black", lwd=2)
}  
mtext(LETTERS[2], side = 3, line = 0.5, adj = -0.15, cex = 1.2, font = 2, col = 'grey30')

# Probabilities
par(mar = c(3,2,2,2))
plot.new()
mtext(LETTERS[3], side = 3, line = 0.5, adj = -0.05, cex = 1.2, font = 2, col = 'grey30')

# All trees
par(mar = c(3,3.5,2,2))
trees <- 1:data$N_trees
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
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(sum_states[[names[n]]], probs)
}
quantiles <- sapply(1:N, calc)
plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                        function(n) quantiles[1:9, n]))

delta <- 0.05 * (display_ylim[2] - display_ylim[1])
display_ylim[1] <- display_ylim[1] - delta   
display_ylim[2] <- display_ylim[2] + delta      

plot(1, type="n", main=main,
     xlim=c(bin_min, bin_max), xlab='', xaxt="n", yaxt="n",
     ylim=display_ylim, ylab='', frame = FALSE)

usr <- par("usr")
segments(x0 = 1,  x1 = 43, y0 = usr[3], y1 = usr[3])
axis(1, at=seq(1,N,10), labels=seq(1980,2020,10), 
     lwd = 0, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))    

axis(2, at=seq(0,1,0.25), seq(0,1,0.25), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))  
title(ylab = ylab, line = 2.5)



for(n in 1:N){
  
  idplot <- c(n*2-1, n*2)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[1,idplot], rev(plot_quantiles[9,idplot])),
          col = util$c_light, border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[2,idplot], rev(plot_quantiles[8,idplot])),
          col = util$c_light_highlight, border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[3,idplot], rev(plot_quantiles[7,idplot])),
          col = util$c_mid, border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[4,idplot], rev(plot_quantiles[6,idplot])),
          col = util$c_mid_highlight, border = NA)
  
  lines(plot_xs[idplot], plot_quantiles[5, idplot],
        col= util$c_dark, lwd=2)
}

mtext(LETTERS[4], side = 3, line = 0.5, adj = -0.15, cex = 1.2, font = 2, col = 'grey30')

dev.off()