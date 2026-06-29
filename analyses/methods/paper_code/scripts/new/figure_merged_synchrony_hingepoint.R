

pdf(file = file.path(figpath, 'merged_synchrony_hingepoint.pdf'),
    width = 10, height = 3.3)
par(mfrow = c(1,2))
par(mar = c(3,4,0.25,1))
# Panel A: synchrony


names <- paste0('acrossstands_synchrony_', names(intervals))

# Construct bins
N <- length(names)
bin_min <- 0.5
bin_max <- N + 0.5
bin_delta <- 1
breaks <- seq(bin_min, bin_max, bin_delta)

plot_config <- util$configure_bin_plotting(breaks)
plot_idxs <- plot_config[[1]]
plot_xs <- plot_config[[2]]

# par(mar = c(0.5,4,0.25,0))
plot(1, type="n", main=main,
     xlim=c(bin_min, bin_max), xlab='', xaxt="n", yaxt = 'n',
     ylim=c(0.49,1), ylab='',
     frame = FALSE, cex.axis = 0.8)
par(xpd=TRUE)

axis(2, at= seq(0.5,1,0.1), labels=seq(50,100,10), cex.axis = 0.8, mgp = c(2,0.75,0), tck = -0.03)
title(ylab = 'Synchrony (% of trees)', line = 2.5)

mtext(LETTERS[1], side = 3, line = -1, adj = -0.17, cex = 1.2, font = 2, col = 'grey30')

# Construct marginal quantiles
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)

calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[names[n]]], probs)
}
quantiles_a <- sapply(1:N, calc)

plot_quantiles_a <- do.call(cbind, lapply(plot_idxs,
                                          function(n) quantiles_a[1:5, n]))


names <- paste0('withinstand_synchrony_', names(intervals))

calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[names[n]]], probs)
}
quantiles_w <- sapply(1:N, calc)

plot_quantiles_w <- do.call(cbind, lapply(plot_idxs,
                                          function(n) quantiles_w[1:5, n]))


lines(c(1:N)+0.2, quantiles_w[3, ],
      col='#754875', lwd=1, lty = 2)

lines(c(1:N)-0.2, quantiles_a[3, ],
      col='#466f8a', lwd=1, lty = 2)

for(i in seq(1, 2*N, 2)){
  
  rect(xleft = plot_xs[i]+0.1, xright = plot_xs[i+1]-0.5,
       ybottom = plot_quantiles_a[1,i], ytop = plot_quantiles_a[5,i+1],
       col = '#b7cddb', border = NA)
  
  rect(xleft = plot_xs[i]+0.1, xright = plot_xs[i+1]-0.5,
       ybottom = plot_quantiles_a[2,i], ytop = plot_quantiles_a[4,i+1],
       col = '#7fa4bb', border = NA)
  
  lines(c(plot_xs[i]+0.1,  plot_xs[i+1]-0.5), plot_quantiles_a[3, i:(i+1)],
        col="#1b3f5a", lwd=2)
  
}

#1b3f5a  #466f8a  #7fa4bb  #b7cddb

for(i in seq(1, 2*N, 2)){
  
  rect(xleft = plot_xs[i]+0.5, xright = plot_xs[i+1]-0.1,
       ybottom = plot_quantiles_w[1,i], ytop = plot_quantiles_w[5,i+1],
       col = '#bba6bb', border = NA)
  
  rect(xleft = plot_xs[i]+0.5, xright = plot_xs[i+1]-0.1,
       ybottom = plot_quantiles_w[2,i], ytop = plot_quantiles_w[4,i+1],
       col = '#a488a4', border = NA)
  
  lines(c(plot_xs[i]+0.5,  plot_xs[i+1]-0.1), plot_quantiles_w[3, i:(i+1)],
        col='#4f1d4b', lwd=2)
  
  lines(c(plot_xs[i]+0.1,  plot_xs[i+1]-0.1), c(0.48,0.48),
        col='black', lwd=1)
  
}

text(y = 0.515, x = 0.3, labels = 'Minimum synchrony', cex = 0.75, col = 'grey70', adj = 0)
segments(y0 = 0.5, x0 = 0.25, x1 = 9.55, lty = 2, col = 'grey70')

# text(x = 1:8, y = 0.45, labels = 
#        c('1980\nto\n1984', '1985\nto\n1989', '1990\nto\n1994', '1995\nto\n1999',
#          '2000\nto\n2004', '2005\nto\n2009', '2010\nto\n2014', '2015\nto\n2022'),
#      cex = 0.8)
par(xpd=NA)
text(x = 1:N, y = 0.44, labels = 
       c('1980-\n1984', '1985-\n1989', '1990-\n1994', '1995-\n1999',
         '2000-\n2004', '2005-\n2009', '2010-\n2014', '2015-\n2019',
         '2020-\n2024'),
     cex = 0.77)

legend("topleft",
       legend = c("Within population", "Between populations"),
       col = c("#bba6bb", "#b7cddb"),
       lwd = 8,
       cex = 0.85,
       text.col = "black",
       horiz = F,
       bty = "n",
       inset = c(0.55, 0.6))
legend("topleft",
       legend = c("", ""),
       col = c("#4f1d4b", "#1b3f5a"),
       lwd = 2,
       bty = "n",
       cex = 0.85,
       text.col = "black",
       horiz = F,
       inset = c(0.55, 0.6))

par(xpd=TRUE)



#Panel B : future synchrony?

years_to_plot <- years_to_predict[10:data$N_max_newyears_cw]
names <- sapply(10:data$N_max_newyears_cw,
                function(n) paste0('withinstand_synchrony_cw[', n, ']'))

ylab="Synchrony (% of trees)"
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
  util$ensemble_mcmc_quantile_est(newsamples[[names[n]]], probs)
}
quantiles <- sapply(1:N, calc)
plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                        function(n) quantiles[1:5, n]))

delta <- 0.05 * (display_ylim[2] - display_ylim[1])
display_ylim[1] <- 0
display_ylim[2] <- display_ylim[2] + delta  

par(mar = c(3,4,0.25,0))
plot(1, type="n", main=main,
     xlim=c(1, N), xlab='', xaxt="n",
     ylim=c(0.49,1), ylab='', yaxt="n",, frame = FALSE)


mtext(LETTERS[2], side = 3, line = -1, adj = -0.17, cex = 1.2, font = 2, col = 'grey30')


usr <- par("usr")
segments(x0 = 0.5,  x1 = N+0.5, y0 = usr[3], y1 = usr[3], lwd = 1)
axis(1, at= seq(0.5,N+0.5,1),
     # labels=years_to_plot, 
     labels = NA,
     lwd = 0, lwd.ticks = 1, tck = -0.013, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))  
axis(1, at= c(0.5,10.5,N+0.5),
     # labels=years_to_plot,
     labels = NA,
     lwd = 0, lwd.ticks = 1.2, tck = -0.04,
     cex.axis = 0.9, mgp = c(2, 0.7, 0))
#title(xlab = "Year", line = 1.5)


axis(2, at= seq(0.5,1,0.1), labels=seq(50,100,10), cex.axis = 0.8, mgp = c(2,0.75,0), tck = -0.03)
title(ylab = ylab, line = 2.5)



for(n in 1:N){
  
  col <- ifelse(years_to_plot[n] %in% c(1989:1998), "#bba6bb", '#bcdccc90')
  idplot <- c(n*2-1, n*2)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[1,idplot], rev(plot_quantiles[5,idplot])),
          col = col, border = NA)
  # polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
  #         c(plot_quantiles[2,idplot], rev(plot_quantiles[8,idplot])),
  #         col = "#99c7b090", border = NA)
  col <- ifelse(years_to_plot[n] %in%  c(1989:1998), "#a488a4", '#7cb99b90')
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[2,idplot], rev(plot_quantiles[4,idplot])),
          col = col, border = NA)
  # polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
  #         c(plot_quantiles[4,idplot], rev(plot_quantiles[6,idplot])),
  #         col = "#50a27990", border = NA)
  
  col <- ifelse(years_to_plot[n] %in%  c(1989:1998), "#4f1d4b", '#278f5b')
  lines(plot_xs[idplot], plot_quantiles[3, idplot],
        col= col, lwd=2)
}


names <- sapply(10:data$N_max_newyears_cw,
                function(n) paste0('withinstand_synchrony_ww[', n, ']'))

# Construct marginal quantiles
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(newsamples[[names[n]]], probs)
}
quantiles <- sapply(1:N, calc)
plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                        function(n) quantiles[1:5, n]))



par(xpd=TRUE)
legend("bottomleft",
       legend = c("Alternation of cold\nand warm summers"),
       col = c("#bba6bb"),
       lwd = 8,
       cex = 0.85,
       text.col = "black",
       horiz = F,
       bty = "n",
       inset = c(-0.06, -0.18))
legend("bottomleft",
       legend = c("\n"),
       col = c("#4f1d4b"),
       lwd = 2,
       bty = "n",
       cex = 0.85,
       text.col = "black",
       horiz = F,
       inset = c(-0.06, -0.18))


legend("bottomleft",
       legend = c("Only warm summers"),
       col = c("#bcdccc"),
       lwd = 8,
       cex = 0.85,
       text.col = "black",
       horiz = F,
       bty = "n",
       inset = c(0.45, -0.18))
legend("bottomleft",
       legend = c(''),
       col = c("#278f5b"),
       lwd = 2,
       bty = "n",
       cex = 0.85,
       text.col = "black",
       horiz = F,
       inset = c(0.45, -0.18))

par(xpd=FALSE)

text(y = 0.515, x = 0.5, labels = 'Minimum synchrony', cex = 0.75, col = 'grey70', adj = 0)
segments(y0 = 0.5, x0 = -0.5, x1 = length(years_to_predict)+0.5, lty = 2, col = 'grey70')

dev.off()
