
# Make new predictions with summer temp. alternating
# Create hypothetical trees c/w. summer alternations
years_to_predict <- seq(1980,  2010, 1)
summertemp_ww <- rnorm(length(years_to_predict), 24, 2/2.57)
summertemp_cw <- c(rnorm(9, 24, 2/2.57), rnorm(1, 15, 1/2.57), rnorm(1, 24, 2/2.57),  rnorm(1, 15, 1/2.57),
                   rnorm(1, 24, 2/2.57),  rnorm(1, 15, 1/2.57), rnorm(1, 24, 2/2.57),  rnorm(1, 15, 1/2.57),
                   rnorm(1, 24, 2/2.57),  rnorm(1, 15, 1/2.57),  rnorm(13, 24, 2/2.57))


trees_per_stand <- 50
unique_stands <- "Benwell"
newtree_cw_stand_idxs <- rep(which(unique_stands==unique_stands), each = trees_per_stand)
N_newtrees_cw <- length(newtree_cw_stand_idxs)
N_max_newyears_cw <- length(years_to_predict)
first_newyear_cw <- min(years_to_predict)
newyears_cw <- c()
newprevsummer_temps_cw <-c()
newprevsummer_temps_ww <-c()
N_newyears_cw <- c()
idx <- 1
newtree_cw_start_idxs <- c()
newtree_cw_end_idxs <- c()
for (i in 1:N_newtrees_cw) {
  
  newyears_tree <- years_to_predict-first_newyear+1
  newyears_cw <- c(newyears_cw, newyears_tree)
  
  N_newyears_tree <- length(newyears_tree)
  N_newyears_cw <- c(N_newyears_cw, N_newyears_tree)
  
  newprevsummer_temps_tree <- summertemp_cw
  newprevsummer_temps_cw <- c(newprevsummer_temps_cw, newprevsummer_temps_tree)
  
  newprevsummer_temps_tree <- summertemp_ww
  newprevsummer_temps_ww <- c(newprevsummer_temps_ww, newprevsummer_temps_tree)
  
  newtree_cw_start_idxs <- c(newtree_cw_start_idxs, idx)
  idx <- idx + N_newyears_tree
  newtree_cw_end_idxs <- c(newtree_cw_end_idxs, idx - 1)
}

# Collect data
N <- length(years)
Nnew <- length(newyears)
Nnew_cw <- length(newyears_cw)
data <- mget(c('N', 'N_trees', 'N_max_years',
               'N_years', 'tree_start_idxs', 'tree_end_idxs',
               'seed_counts', 'years', 
               'prevsummer_temps', 'spring_temps', 'gdd_lastfrost',
               # for predictions
               'Nnew', 'N_newtrees', 'N_max_newyears', 'newtree_stand_idxs',
               'N_newyears', 'newtree_start_idxs', 'newtree_end_idxs','newyears', 
               'newprevsummer_temps', 'newgdd_lastfrost', 'newspring_temps',
               
               # for predictions
               'Nnew_cw', 'N_newtrees_cw', 'N_max_newyears_cw', 'newtree_cw_stand_idxs',
               'N_newyears_cw', 'newtree_cw_start_idxs', 'newtree_cw_end_idxs','newyears_cw', 
               'newprevsummer_temps_cw', 'newprevsummer_temps_ww'
))


sm_gq <- stan_model(file = file.path(wd, "stan", "paper_models", "model2_treelevel_zinb_missing_rescaled_summertemp_springfrost_springtemp_wpred_standlevel.stan"))
newgq <- gqs(sm_gq, data = data, as.matrix(fit_biol))
newsamples <- util$extract_expectand_vals(newgq)


# Manipulate posteriors
s <- 1
for(y in 1:data$N_max_newyears_cw){
  
  newname <- paste0('withinstand_synchrony_ww[',y,']')
  newsamples[[newname]] <- 0
  
  name1 <- paste0('masting_new_plot_ww[', s,',',y,']')
  name2 <- paste0('nonmasting_new_plot_ww[', s,',',y,']')
  
  maxsync <- pmax(newsamples[[name1]], newsamples[[name2]])
  
  newsamples[[newname]] <- maxsync/data$N_newtrees_cw
  
  newname <- paste0('withinstand_synchrony_cw[',y,']')
  newsamples[[newname]] <- 0
  
  name1 <- paste0('masting_new_plot_cw[', s,',',y,']')
  name2 <- paste0('nonmasting_new_plot_cw[', s,',',y,']')
  
  maxsync <- pmax(newsamples[[name1]], newsamples[[name2]])
  
  newsamples[[newname]] <- maxsync/data$N_newtrees_cw
  
}


pdf(file = file.path(figpath, 'coldsummer_hingepoint.pdf'),
    width = 6, height = 3.3)
par(mfrow = c(1,1), mar = c(3,4,0.25,0))

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

plot(1, type="n", main=main,
     xlim=c(1, N), xlab='', xaxt="n",
     ylim=c(0.49,1), ylab='', yaxt="n",, frame = FALSE)

usr <- par("usr")
segments(x0 = 0.5,  x1 = N+0.5, y0 = usr[3], y1 = usr[3], lwd = 1.5)
axis(1, at= seq(0.5,N+0.5,1),
     # labels=years_to_plot, 
     labels = NA,
     lwd = 0, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))  
axis(1, at= c(0.5,10.5,N+0.5),
     # labels=years_to_plot, 
     labels = NA,
     lwd = 0, lwd.ticks = 1.2, tck = -0.07, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))    
#title(xlab = "Year", line = 1.5)

axis(2, at=seq(0.5,1,0.1), seq(50,100,10), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))  
title(ylab = ylab, line = 2.5)

for(n in 1:N){
  
  col <- ifelse(years_to_plot[n] %in% c(1989:1998), "#ccbcdc90", '#bcdccc90')
  idplot <- c(n*2-1, n*2)
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[1,idplot], rev(plot_quantiles[5,idplot])),
          col = col, border = NA)
  # polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
  #         c(plot_quantiles[2,idplot], rev(plot_quantiles[8,idplot])),
  #         col = "#99c7b090", border = NA)
  col <- ifelse(years_to_plot[n] %in%  c(1989:1998), "#9b7cb990", '#7cb99b90')
  polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
          c(plot_quantiles[2,idplot], rev(plot_quantiles[4,idplot])),
          col = col, border = NA)
  # polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
  #         c(plot_quantiles[4,idplot], rev(plot_quantiles[6,idplot])),
  #         col = "#50a27990", border = NA)
  
  col <- ifelse(years_to_plot[n] %in%  c(1989:1998), "#5b278f", '#278f5b')
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
       col = c("#ccbcdc"),
       lwd = 8,
       cex = 0.85,
       text.col = "black",
       horiz = F,
       bty = "n",
       inset = c(0.05, -0.18))
legend("bottomleft",
       legend = c("\n"),
       col = c("#5b278f"),
       lwd = 2,
       bty = "n",
       cex = 0.85,
       text.col = "black",
       horiz = F,
       inset = c(0.05, -0.18))


legend("bottomleft",
       legend = c("Only warm summers"),
       col = c("#bcdccc"),
       lwd = 8,
       cex = 0.85,
       text.col = "black",
       horiz = F,
       bty = "n",
       inset = c(0.55, -0.18))
legend("bottomleft",
       legend = c(''),
       col = c("#278f5b"),
       lwd = 2,
       bty = "n",
       cex = 0.85,
       text.col = "black",
       horiz = F,
       inset = c(0.55, -0.18))

par(xpd=FALSE)

text(y = 0.515, x = 0.5, labels = 'Minimum synchrony', cex = 0.75, col = 'grey70', adj = 0)
segments(y0 = 0.5, x0 = -0.5, x1 = 22.5, lty = 2, col = 'grey70')

dev.off()
