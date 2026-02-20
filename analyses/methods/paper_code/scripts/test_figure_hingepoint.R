

# Manipulate posteriors
s <- 1
for(y in 1:data$N_max_newyears_cw){
   
  newname <- paste0('withinstand_synchrony_ww[',y,']')
  samples[[newname]] <- 0
   
  name1 <- paste0('masting_new_plot_ww[', s,',',y,']')
  name2 <- paste0('nonmasting_new_plot_ww[', s,',',y,']')
    
  maxsync <- pmax(samples[[name1]], samples[[name2]])
    
  samples[[newname]] <- maxsync/data$N_newtrees_cw
  
  newname <- paste0('withinstand_synchrony_cw[',y,']')
  samples[[newname]] <- 0
  
  name1 <- paste0('masting_new_plot_cw[', s,',',y,']')
  name2 <- paste0('nonmasting_new_plot_cw[', s,',',y,']')
  
  maxsync <- pmax(samples[[name1]], samples[[name2]])
  
  samples[[newname]] <- maxsync/data$N_newtrees_cw
  
}


par(mfrow = c(1,1), mar = c(4.5,3.5,2,2))

 
names <- sapply(10:data$N_max_newyears_cw,
                function(n) paste0('withinstand_synchrony_cw[', n, ']'))

ylab="Synchrony"
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
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[names[n]]], probs)
}
quantiles <- sapply(1:N, calc)
plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                        function(n) quantiles[1:5, n]))

delta <- 0.05 * (display_ylim[2] - display_ylim[1])
display_ylim[1] <- 0
display_ylim[2] <- display_ylim[2] + delta      

plot(1, type="n", main=main,
     xlim=c(1, N), xlab='Year', xaxt="n",
     ylim=c(0.5,1), ylab='', yaxt="n",, frame = FALSE)

usr <- par("usr")
segments(x0 = 1,  x1 = 20, y0 = usr[3], y1 = usr[3], lwd = 1)
axis(1, at= seq(1,N,1),
     labels=years_to_predict[10:data$N_max_newyears_cw], 
     lwd = 0, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))    

axis(2, at=seq(0.5,1,0.25), seq(0.5,1,0.25), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))  
title(ylab = ylab, line = 2.5)

for(n in 1:N){
  
  idplot <- c(n*2-1, n*2)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[1,idplot], rev(plot_quantiles[5,idplot])),
          col = "#bcdccc90", border = NA)
  # polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
  #         c(plot_quantiles[2,idplot], rev(plot_quantiles[8,idplot])),
  #         col = "#99c7b090", border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[2,idplot], rev(plot_quantiles[4,idplot])),
          col = "#7cb99b90", border = NA)
  # polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
  #         c(plot_quantiles[4,idplot], rev(plot_quantiles[6,idplot])),
  #         col = "#50a27990", border = NA)
  
  lines(plot_xs[idplot], plot_quantiles[3, idplot],
        col= '#278f5b', lwd=2)
}


names <- sapply(10:data$N_max_newyears_cw,
                function(n) paste0('withinstand_synchrony_ww[', n, ']'))

# Construct marginal quantiles
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[names[n]]], probs)
}
quantiles <- sapply(1:N, calc)
plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                        function(n) quantiles[1:5, n]))


for(n in 1:N){
  
  idplot <- c(n*2-1, n*2)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[1,idplot], rev(plot_quantiles[5,idplot])),
          col = "#e4c19d90", border = NA)
  # polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
  #         c(plot_quantiles[2,idplot], rev(plot_quantiles[8,idplot])),
  #         col = "#deb48990", border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[2,idplot], rev(plot_quantiles[4,idplot])),
          col = "#d49b6190", border = NA)
  # polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
  #         c(plot_quantiles[4,idplot], rev(plot_quantiles[6,idplot])),
  #         col = "#cd8e4d90", border = NA)
  
  lines(plot_xs[idplot], plot_quantiles[3, idplot],
        col= '#c98139', lwd=2)
}


legend("topleft",
       legend = c("Alternation of cold and warm summers", "Only warm summers"),
       col = c("#bcdccc", "#e4c19d"),
       lwd = 8,
       cex = 0.85,
       text.col = "black",
       horiz = F,
       bty = "n",
       inset = c(0.05, 0.025))
legend("topleft",
       legend = c("", ""),
       col = c("#278f5b", "#c98139"),
       lwd = 2,
       bty = "n",
       cex = 0.85,
       text.col = "black",
       horiz = F,
       inset = c(0.05, 0.025))

par(mar = c(4,4,1,1))
plot(summertemp_cw[10:data$N_max_newyears_cw]~years_to_predict[10:data$N_max_newyears_cw], type = 'l', 
     col = '#278f5b', ylim = c(15,27),
     xlab = '', ylab = 'Summer temperature',
     lwd = 2)
lines(summertemp_ww[10:data$N_max_newyears_cw]~years_to_predict[10:data$N_max_newyears_cw], col = '#c98139', lwd = 2)

# util$plot_disc_pushforward_quantiles(samples, paste0('masting_new_plot_cw[', 1,',',1:data$N_max_newyears_cw,']'))
# 
# util$plot_disc_pushforward_quantiles(samples, paste0('masting_new_plot_ww[', 1,',',1:data$N_max_newyears_cw,']'))
# 
# 

