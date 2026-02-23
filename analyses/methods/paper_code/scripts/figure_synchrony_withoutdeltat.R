

# Manipulate posteriors
sites <- stringr::str_split_i(uniq_tree_ids, '_', 1)
uniq_sites <- unique(sites)
synchrony_within <- synchrony_across <- c()
trees <- 1:data$N_trees
for(y in 1980:2022){
  
  ## Look at all trees within one year
  year <- y-1980+1
  
  idxs_tree <- (year+N_max_years*(trees-1))
  names <- sapply(idxs_tree,
                  function(n) paste0('states_pred[', n, ']'))
  
  for(s in 1:length(uniq_sites)){
    name <- paste0('masting_site[', s, ']')
    samples[[name]] <- 0
    name <- paste0('nonmasting_site[', s, ']')
    samples[[name]] <- 0
  }
  
  trees_persite <- rep(0, length(uniq_sites))
  for(t in 1:data$N_trees){
    name <- paste0('masting_site[', which(uniq_sites == sites[t]), ']')
    samples[[name]] <- samples[[name]] + (samples[[names[t]]]-1)
    
    name <- paste0('nonmasting_site[', which(uniq_sites == sites[t]), ']')
    samples[[name]] <- samples[[name]] + abs(samples[[names[t]]]-2)
    
    trees_persite[which(uniq_sites == sites[t])] <-
      trees_persite[which(uniq_sites == sites[t])] + 1
  }
  
  
  newname <- paste0('withinstand_synchrony[',year,']')
  samples[[newname]] <- 0
  for(s in 1:length(uniq_sites)){
    
    name1 <- paste0('masting_site[', s, ']')
    name2 <- paste0('nonmasting_site[', s, ']')
    
    maxsync <- pmax(samples[[name1]], samples[[name2]])
    
    samples[[newname]] <- samples[[newname]] + maxsync/trees_persite[s]
    
    newname_here <- paste0('withinstand_masting_synchrony[',s, ",", year,']')
    samples[[newname_here]] <- samples[[name1]]/trees_persite[s]
    
    newname_here <- paste0('withinstand_nonmasting_synchrony[',s, ",", year,']')
    samples[[newname_here]] <- samples[[name2]]/trees_persite[s]
    
  }
  samples[[newname]] <- samples[[newname]]/length(uniq_sites)
  
  newname <- paste0('acrossstands_masting_synchrony[',year,']')
  samples[[newname]] <- 0
  for(s in 1:length(uniq_sites)){
    name<- paste0('masting_site[', s, ']')
    samples[[newname]] <- samples[[newname]] + samples[[name]]/trees_persite[s]
  }
  samples[[newname]] <- samples[[newname]]/length(uniq_sites)
  
  newname <- paste0('acrossstands_nonmasting_synchrony[',year,']')
  samples[[newname]] <- 0
  for(s in 1:length(uniq_sites)){
    name<- paste0('nonmasting_site[', s, ']')
    samples[[newname]] <- samples[[newname]] + samples[[name]]/trees_persite[s]
  }
  samples[[newname]] <- samples[[newname]]/length(uniq_sites)
  
  name1 <- paste0('acrossstands_nonmasting_synchrony[',year,']')
  name2 <- paste0('acrossstands_masting_synchrony[',year,']')
  newname <- paste0('acrossstands_synchrony[',year,']')
  samples[[newname]] <- pmax(samples[[name1]], samples[[name2]])
  
}

# Summarize by intervals
intervals <- list('1980s' = 1980:1984, '1985s' = 1985:1989,  
                  '1990s' = 1990:1994, '1995s' = 1995:1999, 
                  '2000s' = 2000:2004, '2005s' = 2005:2009, 
                  '2010s' = 2010:2014, '2015s' = 2015:2019,
                  '2020s' = 2020:2022)

for(i in 1:length(intervals)){
  newname <- paste0('withinstand_synchrony_', names(intervals)[i])
  samples[[newname]] <- 0
  for(y in intervals[[i]]){
    name <- paste0('withinstand_synchrony[', y-1980+1, ']')
    samples[[newname]] <- samples[[newname]] + samples[[name]]
  }
  samples[[newname]] <- samples[[newname]]/length(intervals[[i]])
}

# par(mfrow=c(1, 1), mar=c(2,5,0.5,0.5), cex.lab = 1)
names <- paste0('withinstand_synchrony_', names(intervals))
# util$plot_disc_pushforward_quantiles(samples, names, xticklabs = names(intervals),
#                                      ylab = 'Within-stand synchrony (% stands)',
#                                      display_ylim = c(0,1))


for(i in 1:length(intervals)){
  newname <- paste0('acrossstands_synchrony_', names(intervals)[i])
  samples[[newname]] <- 0
  for(y in intervals[[i]]){
    name <- paste0('acrossstands_synchrony[', y-1980+1, ']')
    samples[[newname]] <- samples[[newname]] + samples[[name]]
  }
  samples[[newname]] <- samples[[newname]]/length(intervals[[i]])
}

names <- paste0('acrossstands_synchrony_', names(intervals))
# util$plot_disc_pushforward_quantiles(samples, names, xticklabs = names(intervals),
#                                      ylab = 'Across-stand synchrony (% stands)',
#                                      display_ylim = c(0,1))



# Plot
pdf(file = file.path(figpath, 'synchrony_withoutdeltat.pdf'),
    width = 6, height = 3.3)

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

par(mar = c(0.5,4,0.25,0))
plot(1, type="n", main=main,
     xlim=c(bin_min, bin_max), xlab='', xaxt="n", yaxt = 'n',
     ylim=c(0.43,1.05), ylab='Synchrony (% of trees)',
     frame = FALSE, cex.axis = 0.9)

axis(2, at= seq(0.5,1,0.1), labels=seq(50,100,10), cex.axis = 0.9)

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



# text(x = 1:8, y = 0.45, labels = 
#        c('1980\nto\n1984', '1985\nto\n1989', '1990\nto\n1994', '1995\nto\n1999',
#          '2000\nto\n2004', '2005\nto\n2009', '2010\nto\n2014', '2015\nto\n2022'),
#      cex = 0.8)

text(x = 1:N, y = 0.44, labels = 
       c('1980-\n1984', '1985-\n1989', '1990-\n1994', '1995-\n1999',
         '2000-\n2004', '2005-\n2009', '2010-\n2014', '2015-\n2019',
         '2020-\n2022'),
     cex = 0.77)

legend("topleft",
       legend = c("Within population", "Between populations"),
       col = c("#bba6bb", "#b7cddb"),
       lwd = 8,
       cex = 0.85,
       text.col = "black",
       horiz = F,
       bty = "n",
       inset = c(0.05, 0.07))
legend("topleft",
       legend = c("", ""),
       col = c("#754875", "#466f8a"),
       lwd = 2,
       bty = "n",
       cex = 0.85,
       text.col = "black",
       horiz = F,
       inset = c(0.05, 0.07))

text(y = 0.515, x = 0.3, labels = 'Minimum synchrony', cex = 0.75, col = 'grey70', adj = 0)
segments(y0 = 0.5, x0 = 0, x1 = 9.5, lty = 2, col = 'grey70')


dev.off()
