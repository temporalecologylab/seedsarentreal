idxs_tree <-(1+N_max_years*(t-1)):(N_max_years*t)
observed_idxs_tree <- tree_start_idxs[t]:tree_end_idxs[t]
observed_flags <- years[tree_start_idxs[t]:tree_end_idxs[t]]


par(mfrow=c(1, 1))
names <- sapply(idxs_tree,
                function(n) paste0('states_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, names, obs_xs = first_year:last_year,
                                       bin_min = first_year, bin_max = last_year+1, bin_delta = 1,
                                       ylab = 'Probability of masting')



par(mfrow=c(4, 2))
for(s in unique(raw_data$site.ID)){
  
  trees_s <- which(grepl(s, uniq_tree_ids))
  
  sum_states <- list()
  for(y in 1:N_max_years){
    
    idxs_tree <- (y+N_max_years*(trees_s-1))
    names <- sapply(idxs_tree,
                    function(n) paste0('states_pred[', n, ']'))
    
    sum_states_y <- Reduce("+", samples[names])
    
    sum_states[[y]] <- sum_states_y - length(trees_s)
    
    
  }
  names(sum_states) <- paste0('sum_states[', 1:N_max_years, ']')
  
  util$plot_conditional_median_quantiles(sum_states, names(sum_states), obs_xs = first_year:last_year,
                                         bin_min = first_year, bin_max = last_year+1, bin_delta = 1,
                                         ylab = 'Number of masting trees', main = s)
  points(2013.5, -0.5, pch = 16, col = '#008080', cex = 2)
  
  
}

par(mfrow=c(1, 1))
trees_s <- which(grepl( paste(unique(raw_data$site.ID),collapse="|"), uniq_tree_ids))
idxs_tree <-(1+N_max_years*(min(trees_s)-1)):(N_max_years*max(trees_s))


sum_states <- list()
for(y in 1:N_max_years){
  
  idxs_tree <- (y+N_max_years*(trees_s-1))
  names <- sapply(idxs_tree,
                  function(n) paste0('states_pred[', n, ']'))
  
  sum_states_y <- Reduce("+", samples[names])
  
  sum_states[[y]] <- sum_states_y - length(trees_s)
  
  
}
names(sum_states) <- paste0('sum_states[', 1:N_max_years, ']')
util$plot_conditional_median_quantiles(sum_states, names(sum_states), obs_xs = first_year:last_year,
                                       bin_min = first_year, bin_max = last_year+1, bin_delta = 1,
                                       ylab = 'Number of masting trees', main = '')


observed_flags <- unlist(sapply(1:N_trees, function(t){
  idxs_tree <-(1+N_max_years*(t-1)):(N_max_years*t)
  observed_flags_tree <- years[tree_start_idxs[t]:tree_end_idxs[t]]
  return(idxs_tree[observed_flags_tree])
}))

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
par(mfrow=c(1, 1), mar=c(2,4,0.5,0.5))
util$plot_conditional_median_quantiles(sum_states, names(sum_states), obs_xs = first_year:last_year,
                                       bin_min = first_year-0.5, bin_max = last_year+0.5, bin_delta = 1,
                                       ylab = 'Proportion of masting trees', main = '')
for(y in 1:N_max_years){
  idxs_tree <- (y+N_max_years*(trees_s-1))
  idxs_tree_observed <- idxs_tree[idxs_tree %in% observed_flags]
  text(x = 1979+y, y = util$ensemble_mcmc_quantile_est(sum_states[[y]], 0.9) + 0.05, 
       labels = length(idxs_tree_observed), cex = 0.75)
}





sapply(1:N_max_years, function(y){
  
  return(length(idxs_tree_observed))
})


par(mfrow=c(4, 2))
for(s in unique(raw_data$site.ID)){http://127.0.0.1:13865/graphics/0f8b3f05-bc7d-41c4-a0b5-426485b13168.png
  
  trees_s <- which(grepl(s, uniq_tree_ids))
  display_ylim=c(0, 1)
  
  # Construct bins
  N <- N_max_years
  bin_min <- 0.5
  bin_max <- N + 0.5
  bin_delta <- 1
  breaks <- seq(bin_min, bin_max, bin_delta)
  
  plot_config <- util$configure_bin_plotting(breaks)
  plot_idxs <- plot_config[[1]]
  plot_xs <- plot_config[[2]]
  
  plot(1, type="n", main= s,
       xlim=c(bin_min, bin_max), xlab= '', xaxt="n",
       ylim=display_ylim, ylab='Proportion of masting trees')
  axis(1, at=1:N, labels=first_year:last_year)    
  
  sum_states <- list()
  for(y in 1:N_max_years){
    
    idxs_tree <- (y+N_max_years*(trees_s-1))
    idxs_tree_observed <- idxs_tree[idxs_tree %in% observed_flags]
    idplot <- c(y*2-1, y*2)
    
    if(length(idxs_tree_observed) > 0){
      names <- sapply(idxs_tree_observed,
                      function(n) paste0('states_pred[', n, ']'))
      sum_states_y <- Reduce("+", samples[names])
      sum_states[[y]] <- (sum_states_y - length(idxs_tree_observed))/length(idxs_tree_observed)
      
      
      probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
      plot_quantiles <- util$ensemble_mcmc_quantile_est(sum_states[[y]], probs)
      
      polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
              c(rep(plot_quantiles[1],2), rep(plot_quantiles[9],2)),
              col = util$c_light, border = NA)
      polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
              c(rep(plot_quantiles[2],2), rep(plot_quantiles[8],2)),
              col = util$c_light_highlight, border = NA)
      polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
              c(rep(plot_quantiles[3],2), rep(plot_quantiles[7],2)),
              col = util$c_mid, border = NA)
      polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
              c(rep(plot_quantiles[4],2), rep(plot_quantiles[6],2)),
              col = util$c_mid_highlight, border = NA)
      
      lines(plot_xs[idplot], rep(plot_quantiles[5],2),
            col= util$c_dark, lwd=2)
      
    }else{
      
      points(mean(plot_xs[idplot]), 0.5, pch = 4)
      sum_states[[y]] <- NA

    }
    
  }

}

par(mfrow=c(1, 1))
trees_s <- which(grepl( paste(unique(raw_data$site.ID),collapse="|"), uniq_tree_ids))
idxs_tree <-(1+N_max_years*(min(trees_s)-1)):(N_max_years*max(trees_s))
