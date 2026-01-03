
# Retrodictive check
pdf(file = file.path('/home/victor/projects/seedsarentreal/docs/reports/external_figures', 'retrodictive_check.pdf'),
    width = 6, height = 4)
observed_idxs <- c()
for (t in 1:data$N_trees){
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  observed_idxs_tree <- N_max_years*(t-1)+years[idxs]
  observed_idxs <- c(observed_idxs, observed_idxs_tree)
}
par(mfrow=c(1, 1), mar = c(4,4,2,2))
names <- sapply(observed_idxs, function(n) paste0('seed_counts_pred[',n,']'))
util$plot_hist_quantiles(samples[names], 'seed_counts_pred', 0, 340, 20,
                         baseline_values=data$seed_counts, 
                         xlab="Number of seeds")

dev.off()

# Tree-level retrodictive checks
pdf(file = file.path('/home/victor/projects/seedsarentreal/docs/reports/external_figures', 'tree_checks.pdf'),
    width = 10, height = 6)
par(mfrow = c(2,2), mar = c(3,3.5,2,2))
for(t in c(1,13,34,50)){
  
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
  segments(x0 = 1,  x1 = 43, y0 = usr[3], y1 = usr[3], lwd = 1.5)
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
  
  
}
dev.off()

# Diagram quantiles
pdf(file = file.path('/home/victor/projects/seedsarentreal/docs/reports/external_figures', 'quantiles.pdf'),
    width = 6, height = 5)
par(mfrow = c(1,1), mar = c(0.5,0.5,0.5,0.5))
plot(1, type="n", main=main,
     xlim=c(0, 13), xlab='', xaxt="n",
     ylim=c(100,300), ylab='', yaxt="n",, frame = FALSE)

usr <- par("usr")
# segments(x0 = 1,  x1 = 43, y0 = usr[3], y1 = usr[3], lwd = 1.5)
# axis(1, at=seq(1,N,10), labels=seq(1980,2020,10), 
#      lwd = 0, lwd.ticks = 1, tck = -0.03, 
#      cex.axis = 0.9, mgp = c(2, 0.7, 0))    

# axis(2, at=seq(0,400,100), seq(0,400,100), 
#      lwd = 1, lwd.ticks = 1, tck = -0.03, 
#      cex.axis = 0.9, mgp = c(2, 0.7, 0))  
# title(ylab = ylab, line = 2.5)

for(n in 8){
  
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
        col= ifelse(n%in%observed_flags,util$c_dark, 'grey50'), lwd=4)
}

for(i in observed_idxs_tree[8]) {
  lines(c(years[i] - 0.5, years[i] + 0.5),
        rep(seed_counts[i], 2),
        col="white", lwd=8)
  lines(c(years[i] - 0.5, years[i] + 0.5),
        rep(data$seed_counts[i], 2),
        col="black", lwd=4)
}  

lines(c(years[i] - 1, years[i] -1),
      c(plot_quantiles[4,idplot[1]], plot_quantiles[6,idplot[1]]),
      col="grey30", lwd=2)
text(y =  mean(c(plot_quantiles[4,idplot[1]], plot_quantiles[6,idplot[1]])), x = years[i] - 1.5,
     labels = c('40th-60th quantile'), srt = 90, col = util$c_mid_highlight)

lines(c(years[i] - 3, years[i] -3),
      c(plot_quantiles[3,idplot[1]], plot_quantiles[7,idplot[1]]),
      col="grey30", lwd=2)
text(y =  mean(c(plot_quantiles[3,idplot[1]], plot_quantiles[7,idplot[1]])), x = years[i] - 3.5,
     labels = c('30th-70th quantile'), srt = 90, col = util$c_mid)

lines(c(years[i] - 5, years[i] -5),
      c(plot_quantiles[2,idplot[1]], plot_quantiles[8,idplot[1]]),
      col="grey30", lwd=2)
text(y =  mean(c(plot_quantiles[2,idplot[1]], plot_quantiles[8,idplot[1]])), x = years[i] - 5.5,
     labels = c('20th-80th quantile'), srt = 90, col = util$c_light_highlight)

lines(c(years[i] - 7, years[i] -7),
      c(plot_quantiles[1,idplot[1]], plot_quantiles[9,idplot[1]]),
      col="grey30", lwd=2)
text(y =  mean(c(plot_quantiles[1,idplot[1]], plot_quantiles[9,idplot[1]])), x = years[i] - 7.5,
     labels = c('10th-90th quantile'), srt = 90, col = util$c_light)

text(y =  c(plot_quantiles[5,idplot[1]]), x = years[i] + 3.1,
     labels = c('50th quantile (median)'), srt = 0, col = util$c_dark)

text(y =  c(seed_counts[i]), x = years[i] + 2.3,
     labels = c('Observations'), srt = 0, col = 'black')

dev.off()

# Tree-level retrodictive checks
pdf(file = file.path('/home/victor/projects/seedsarentreal/docs/reports/external_figures', 'tree_states.pdf'),
    width = 10, height = 6)
par(mfrow = c(2,2), mar = c(3,3.5,2,2))
for(t in c(1,13,34,50)){
  
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
  probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  calc <- function(n) {
    util$ensemble_mcmc_quantile_est(samples[[names[n]]]-1, probs)
  }
  quantiles <- sapply(1:N, calc)
  plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                          function(n) quantiles[1:9, n]))
  
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
  
  axis(2, at=seq(0,1,1), c('off-year', 'on-year'), 
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
  
  # for(i in observed_idxs_tree) {
  #   lines(c(years[i] - 0.5, years[i] + 0.5),
  #         rep(seed_counts[i], 2),
  #         col="white", lwd=4)
  #   lines(c(years[i] - 0.5, years[i] + 0.5),
  #         rep(data$seed_counts[i], 2),
  #         col="black", lwd=2)
  # }  
  
  
}
dev.off()


# All trees!
pdf(file = file.path('/home/victor/projects/seedsarentreal/docs/reports/external_figures', 'masting_trees.pdf'),
    width = 6, height = 4)
par(mar = c(3,3.5,1,2))
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

ylab="Trees in an on-year (%)"
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
segments(x0 = 1,  x1 = 43, y0 = usr[3], y1 = usr[3], lwd = 1.5)
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
dev.off()




# Synchrony
pdf(file = file.path('/home/victor/projects/seedsarentreal/docs/reports/external_figures', 'synchrony.pdf'),
    width = 10, height = 6)
# Manipulate posteriors
sites <- stringr::str_split_i(uniq_tree_ids, '_', 1)
uniq_sites <- unique(sites)
synchrony_within <- synchrony_across <- c()
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
  
  obstrees_persite <- rep(0, length(uniq_sites))
  trees_persite <- rep(0, length(uniq_sites))
  for(t in 1:data$N_trees){
    name <- paste0('masting_site[', which(uniq_sites == sites[t]), ']')
    samples[[name]] <- samples[[name]] + (samples[[names[t]]]-1)
    
    name <- paste0('nonmasting_site[', which(uniq_sites == sites[t]), ']')
    samples[[name]] <- samples[[name]] + abs(samples[[names[t]]]-2)
    
    obstrees_persite[which(uniq_sites == sites[t])] <- 
      obstrees_persite[which(uniq_sites == sites[t])] +
      ifelse(idxs_tree[t]%in%observed_flags,1, 0)
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
                  '2010s' = 2010:2014, '2015s' = 2015:2022)

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
     ylim=c(0.45,1.05), ylab='Synchrony',
     frame = FALSE)

axis(2, at= seq(0.5,1,0.1), labels=seq(0.5,1,0.1), cex.axis = 0.9)

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


lines(c(1:8)+0.2, quantiles_w[3, ],
      col=util$c_mid, lwd=1, lty = 2)

lines(c(1:8)-0.2, quantiles_a[3, ],
      col=util$c_mid_teal, lwd=1, lty = 2)

for(i in seq(1, 2*N, 2)){
  
  rect(xleft = plot_xs[i]+0.1, xright = plot_xs[i+1]-0.5,
       ybottom = plot_quantiles_a[1,i], ytop = plot_quantiles_a[5,i+1],
       col = '#a6bbbb', border = NA)
  
  rect(xleft = plot_xs[i]+0.1, xright = plot_xs[i+1]-0.5,
       ybottom = plot_quantiles_a[2,i], ytop = plot_quantiles_a[4,i+1],
       col = '#88a4a4', border = NA)
  
  lines(c(plot_xs[i]+0.1,  plot_xs[i+1]-0.5), plot_quantiles_a[3, i:(i+1)],
        col=util$c_dark_teal, lwd=2)
  
}


for(i in seq(1, 2*N, 2)){
  
  rect(xleft = plot_xs[i]+0.5, xright = plot_xs[i+1]-0.1,
       ybottom = plot_quantiles_w[1,i], ytop = plot_quantiles_w[5,i+1],
       col = util$c_light, border = NA)
  
  rect(xleft = plot_xs[i]+0.5, xright = plot_xs[i+1]-0.1,
       ybottom = plot_quantiles_w[2,i], ytop = plot_quantiles_w[4,i+1],
       col = util$c_light_highlight, border = NA)
  
  lines(c(plot_xs[i]+0.5,  plot_xs[i+1]-0.1), plot_quantiles_w[3, i:(i+1)],
        col=util$c_dark, lwd=2)
  
  lines(c(plot_xs[i]+0.1,  plot_xs[i+1]-0.1), c(0.5,0.5),
        col='black', lwd=1)
  
}

# text(x = 1:8, y = 0.45, labels = 
#        c('1980\nto\n1984', '1985\nto\n1989', '1990\nto\n1994', '1995\nto\n1999',
#          '2000\nto\n2004', '2005\nto\n2009', '2010\nto\n2014', '2015\nto\n2022'),
#      cex = 0.8)

text(x = 1:8, y = 0.47, labels = 
       c('1980-\n1984', '1985-\n1989', '1990-\n1994', '1995-\n1999',
         '2000-\n2004', '2005-\n2009', '2010-\n2014', '2015-\n2022'),
     cex = 0.77)

legend("topleft", 
       legend = c("Within stand", "Between stands"), 
       col = c(util$c_mid_highlight, util$c_mid_teal), 
       lwd = 2, 
       bty = "n", 
       cex = 0.9, 
       text.col = "black", 
       horiz = T, 
       inset = c(0.05, 0.07))
dev.off()
