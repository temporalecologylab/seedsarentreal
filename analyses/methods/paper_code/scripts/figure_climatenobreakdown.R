


pdf(file = file.path(figpath, 'climate_nobreakdown.pdf'),
    width = 6.5, height = 4)
par(mfrow = c(1,1), mar = c(4.5,3.5,2,2))


idxs_tree <- 1:data$N_max_newyears
names <- sapply(idxs_tree,
                function(n) paste0('states_new_plot[1,', n, ']'))

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
  util$ensemble_mcmc_quantile_est(samples_brk[[names[n]]]/data$N_newtrees, probs)
}
quantiles <- sapply(1:N, calc)
plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                        function(n) quantiles[1:9, n]))

delta <- 0.05 * (display_ylim[2] - display_ylim[1])
display_ylim[1] <- 0
display_ylim[2] <- display_ylim[2] + delta      

plot(1, type="n", main=main,
     xlim=c(10, N), xlab='Previous summer temperature (°C)', xaxt="n",
     ylim=display_ylim, ylab='', yaxt="n",, frame = FALSE)

usr <- par("usr")
segments(x0 = 10,  x1 = 121, y0 = usr[3], y1 = usr[3], lwd = 1.5)
axis(1, at= c(19, 34, 48, 62, 76, 91, 105, 119),
     labels=seq(18,25,1), 
     lwd = 0, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))    

axis(2, at=seq(0,1,0.25), seq(0,100,25), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))  
title(ylab = ylab, line = 2.5)

for(n in 10:N){
  
  idplot <- c(n*2-1, n*2)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[1,idplot], rev(plot_quantiles[9,idplot])),
          col = "#bcdccc90", border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[2,idplot], rev(plot_quantiles[8,idplot])),
          col = "#99c7b090", border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[3,idplot], rev(plot_quantiles[7,idplot])),
          col = "#7cb99b90", border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[4,idplot], rev(plot_quantiles[6,idplot])),
          col = "#50a27990", border = NA)
  
  lines(plot_xs[idplot], plot_quantiles[5, idplot],
        col= '#278f5b', lwd=2)
}


names <- sapply(idxs_tree,
                function(n) paste0('states_new_brkd_plot[1,', n, ']'))

# Construct marginal quantiles
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples_brk[[names[n]]]/data$N_newtrees, probs)
}
quantiles <- sapply(1:N, calc)
plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                        function(n) quantiles[1:9, n]))


for(n in 10:N){
  
  idplot <- c(n*2-1, n*2)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[1,idplot], rev(plot_quantiles[9,idplot])),
          col = "#e4c19d90", border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[2,idplot], rev(plot_quantiles[8,idplot])),
          col = "#deb48990", border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[3,idplot], rev(plot_quantiles[7,idplot])),
          col = "#d49b6190", border = NA)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[4,idplot], rev(plot_quantiles[6,idplot])),
          col = "#cd8e4d90", border = NA)
  
  lines(plot_xs[idplot], plot_quantiles[5, idplot],
        col= '#c98139', lwd=2)
}


legend("topleft",
       legend = c("Constraints buffer warming", "Masting breakdown"),
       col = c("#bcdccc", "#e4c19d"),
       lwd = 8,
       cex = 0.85,
       text.col = "black",
       horiz = F,
       bty = "n",
       inset = c(0.55, 0.75))
legend("topleft",
       legend = c("", ""),
       col = c("#278f5b", "#c98139"),
       lwd = 2,
       bty = "n",
       cex = 0.85,
       text.col = "black",
       horiz = F,
       inset = c(0.55, 0.75))



dev.off()
