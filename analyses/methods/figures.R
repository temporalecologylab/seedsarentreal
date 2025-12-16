
samples <- util$extract_expectand_vals(fit_full_const)

# select... all trees
trees_s <- which(grepl( paste(unique(raw_data$site.ID),collapse="|"), uniq_tree_ids))
idxs_tree <-(1+N_max_years*(min(trees_s)-1)):(N_max_years*max(trees_s))

# which predictions correspond to observations (we will not show missing observations here)
observed_flags <- unlist(sapply(1:N_trees, function(t){
  idxs_tree <-(1+N_max_years*(t-1)):(N_max_years*t)
  observed_flags_tree <- years[tree_start_idxs[t]:tree_end_idxs[t]]
  return(idxs_tree[observed_flags_tree])
}))

# prepare data for Mike's function
sum_states <- list()
for(y in 1:N_max_years){
  
  idxs_tree <- (y+N_max_years*(trees_s-1))
  idxs_tree_observed <- idxs_tree[idxs_tree %in% observed_flags]
  names <- sapply(idxs_tree_observed,
                  function(n) paste0('states_pred[', n, ']'))
  
  sum_states_y <- Reduce("+", samples[names])
  
  sum_states[[y]] <- (sum_states_y - length(idxs_tree_observed))/length(idxs_tree_observed)
  
  
}
names(sum_states) <- paste0('sum_states[', 1:N_max_years, ']')

layout(matrix(c(1,2), ncol = 1), height = c(1,1))
par(mfrow=c(1, 1), mar=c(2,5,0.5,0.5), cex.lab = 1)
util$plot_conditional_median_quantiles(sum_states, names(sum_states), obs_xs = first_year:last_year,
                                       bin_min = first_year-0.5, bin_max = last_year+0.5, bin_delta = 1,
                                       ylab = 'Proportion of masting trees (observed only)', main = '')
nb_trees <- data.frame()
for(y in 1:N_max_years){
  idxs_tree <- (y+N_max_years*(trees_s-1))
  idxs_tree_observed <- idxs_tree[idxs_tree %in% observed_flags]
  # text(x = 1979+y, y = util$ensemble_mcmc_quantile_est(sum_states[[y]], 0.9) + 0.05, 
  #      labels = length(idxs_tree_observed), cex = 0.75)
  nb_trees <- rbind(
    nb_trees,
    data.frame(year = 1979+y, count = length(idxs_tree_observed))
  )
}

par(mfrow=c(1, 1), mar=c(0.5,5, 1,0.5), cex.lab = 1.5)
barplot(height = nb_trees$count, names = '', ylim = c(0,60),
        ylab = '# trees', col = 'grey80', border = NA)

## Same figure than above, but will all trees all years
# select... all trees
trees_s <- which(grepl( paste(unique(raw_data$site.ID),collapse="|"), uniq_tree_ids))
idxs_tree <-(1+N_max_years*(min(trees_s)-1)):(N_max_years*max(trees_s))


# prepare data for Mike's function
sum_states <- list()
for(y in 1:N_max_years){
  
  idxs_tree <- (y+N_max_years*(trees_s-1))
  names <- sapply(idxs_tree,
                  function(n) paste0('states_pred[', n, ']'))
  
  sum_states_y <- Reduce("+", samples[names])
  
  sum_states[[y]] <- (sum_states_y - length(idxs_tree))/length(idxs_tree)
  
  
}
names(sum_states) <- paste0('sum_states[', 1:N_max_years, ']')

layout(matrix(c(1,2), ncol = 1), height = c(1,1))
par(mfrow=c(1, 1), mar=c(2,5,0.5,0.5), cex.lab = 1)
util$plot_conditional_median_quantiles(sum_states, names(sum_states), obs_xs = first_year:last_year,
                                       bin_min = first_year-0.5, bin_max = last_year+0.5, bin_delta = 1,
                                       ylab = 'Proportion of masting trees (all)', main = '')


## Look at all trees within one year
sites <- stringr::str_split_i(uniq_tree_ids, '_', 1)
par(mfrow=c(1, 1), mar=c(2,4, 1,1), cex.lab = 1)
year <- 2018-1980+1
{
  xlab="Tree"
  xticklabs=1:data$N_trees
  ylab="Tree-level state"
  display_ylim=c(0, 1.2)
  main=paste("Tree", uniq_tree_ids[t])  
  
  idxs_tree <- (year+N_max_years*(trees_s-1))
  names <- sapply(idxs_tree,
                  function(n) paste0('states_pred[', n, ']'))
  
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
  
  plot(1, type="n", main='',
       xlim=c(bin_min, bin_max), xlab=xlab, xaxt="n",
       ylim=display_ylim, ylab=ylab, xaxs = 'i', yaxt="n",)
  # axis(1, at=1:N, labels=xticklabs[1:N])    
  axis(2, at=c(0,1), labels=c('Non-mast', 'Mast'))   

  
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
          col= ifelse(n%in%observed_flags,util$c_dark, 'grey40'), lwd=2)
    
  }
  
  abline(v=as.numeric(cumsum(table(sites)))+0.5, col="black")
  text(y = 1.15, x = cumsum(table(sites))-table(sites)/2, label = unique(sites), cex = 0.75, adj = 0.5)
}


## FIGURE, one tree only ## 
t <- which(uniq_tree_ids == 'Buckholt_6')
idxs_tree <-(1+N_max_years*(t-1)):(N_max_years*t)
observed_idxs_tree <- tree_start_idxs[t]:tree_end_idxs[t]
observed_flags <- years[tree_start_idxs[t]:tree_end_idxs[t]]

names <- sapply(idxs_tree,
                function(n) paste0('seed_counts_pred[', n, ']'))
xlab="Year"
xticklabs=first_year:last_year
ylab="Seed counts"
display_ylim=c(0, 300)
main=paste("Tree", uniq_tree_ids[t])

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
display_ylim[1] <- display_ylim[1] - delta
display_ylim[2] <- display_ylim[2] + delta      

plot(1, type="n", main='',
     xlim=c(bin_min+2, bin_max), xlab=xlab, xaxt="n",
     ylim=display_ylim, ylab=ylab)
axis(1, at=3:N, labels=xticklabs[3:N])    

for(n in 3:N){
  
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

## FIGURE WITH CLIMATE##
logit <- function(p){
  return(log(p/(1-p)))
}

invlogit <- function(r){
  return(1/(1+exp(-r)))
}

temp0 <- 15
tau_nm_m_df <- data.frame()
for(temp_prevsummer in seq(10,30,0.1)){
  tau_nm_m <- matrix(nrow = nrow(samples[['tau_nm_m0']]), ncol = ncol(samples[['tau_nm_m0']]))
  for(c in 1:nrow(tau_nm_m)){
    tau_nm_m[c,] <- sapply(1:ncol(tau_nm_m), function(i) invlogit(logit(samples[['tau_nm_m0']][c,i])+samples[['beta_nm_m']][c,i]*(temp_prevsummer-temp0)))
  }
  tau_nm_m_df <- rbind(tau_nm_m_df,
                       data.frame(temp_prevsummer, t(util$ensemble_mcmc_quantile_est(tau_nm_m, probs = c(0.1,0.5,0.9)))))
}

temp0 <- 15
tau_m_m_df <- data.frame()
for(temp_prevsummer in seq(10,30,0.1)){
  tau_m_m <- matrix(nrow = nrow(samples[['tau_m_m0']]), ncol = ncol(samples[['tau_m_m0']]))
  for(c in 1:nrow(tau_nm_m)){
    tau_m_m[c,] <- sapply(1:ncol(tau_nm_m), function(i) invlogit(logit(samples[['tau_m_m0']][c,i])+samples[['beta_m_m']][c,i]*(temp_prevsummer-temp0)))
  }
  tau_m_m_df <- rbind(tau_m_m_df,
                       data.frame(temp_prevsummer, t(util$ensemble_mcmc_quantile_est(tau_m_m, probs = c(0.1,0.5,0.9)))))
}


gdd0 <- 15
lambda2_df <- data.frame()
for(gdd_lastfrost in seq(0,40,0.1)){
  lambda2 <- matrix(nrow = nrow(samples[['lambda20']]), ncol = ncol(samples[['lambda20']]))
  for(c in 1:nrow(lambda2)){
    lambda2[c,] <- sapply(1:ncol(lambda2), function(i) exp(log(samples[['lambda20']][c,i]) +samples[['beta_lambda2_frost']][c,i]*(gdd_lastfrost-gdd0)))
  }
  lambda2_df <- rbind(lambda2_df,
                      data.frame(gdd_lastfrost, t(util$ensemble_mcmc_quantile_est(lambda2, probs = c(0.1,0.5,0.9)))))
}

layout(matrix(c(1,2))) # layout(matrix(c(1,1,2, 3, 4, 4), nrow = 2, byrow = T)) if I want to add parameter plot?
par(mar=c(4,6,0.5,0.5), cex.lab = 1.2, cex.axis = 1)
plot(1, type="n", main='',
     xlim=c(10, 30), xlab='Average previous summer max. temperature (degC)',
     ylim=c(0,1), ylab='Probability of transtion\nfrom non-masting to masting')
lines(tau_nm_m_df$X50. ~ tau_nm_m_df$temp_prevsummer, col = util$c_dark, lwd = 2)
lines(tau_nm_m_df$X10. ~ tau_nm_m_df$temp_prevsummer, col = util$c_mid, lty = 2, lwd = 1.2)
lines(tau_nm_m_df$X90. ~ tau_nm_m_df$temp_prevsummer, col = util$c_mid, lty = 2, lwd = 1.2)
par(new = TRUE, mar=c(4,6,10.5,0.5))
boxplot( data$prevsummer_temps, horizontal = TRUE, axes = FALSE,
        col = util$c_light, border = util$c_mid, ylim = c(10, 30), cex = 0.5)
par(mar=c(4,6,0.5,0.5), cex.lab = 1.2, cex.axis = 1)
plot(1, type="n", main='',
     xlim=c(0, 40), xlab='GDD until the last frost day (x10degC)',
     ylim=c(100,220), ylab='Masting intensity\n(seed count)')
lines(lambda2_df$X50. ~ lambda2_df$gdd_lastfrost, col = util$c_dark, lwd = 2)
lines(lambda2_df$X10. ~ lambda2_df$gdd_lastfrost, col = util$c_mid, lty = 2, lwd = 1.2)
lines(lambda2_df$X90. ~ lambda2_df$gdd_lastfrost, col = util$c_mid, lty = 2, lwd = 1.2)
par(new = TRUE, mar=c(4,6,10.5,0.5))
boxplot( data$gdd_lastfrost, horizontal = TRUE, axes = FALSE,
         col = util$c_light, border = util$c_mid, ylim = c(0, 40), cex = 0.5)



layout(matrix(c(1,2))) # layout(matrix(c(1,1,2, 3, 4, 4), nrow = 2, byrow = T)) if I want to add parameter plot?
par(mar=c(4,6,0.5,0.5), cex.lab = 1.2, cex.axis = 1)
plot(1, type="n", main='',
     xlim=c(10, 30), xlab='Average previous summer max. temperature (degC)',
     ylim=c(0 ,1), ylab='Probability of staying in masting')
lines(tau_m_m_df$X50. ~ tau_m_m_df$temp_prevsummer, col = util$c_dark, lwd = 2)
lines(tau_m_m_df$X10. ~ tau_m_m_df$temp_prevsummer, col = util$c_mid, lty = 2, lwd = 1.2)
lines(tau_m_m_df$X90. ~ tau_m_m_df$temp_prevsummer, col = util$c_mid, lty = 2, lwd = 1.2)
par(new = TRUE, mar=c(3.5,6,11.5,0.5))
boxplot( data$prevsummer_temps, horizontal = TRUE, axes = FALSE,
         col = util$c_light, border = util$c_mid, ylim = c(10, 30), cex = 0.5)

