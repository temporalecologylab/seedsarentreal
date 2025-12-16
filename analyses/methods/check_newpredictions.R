



par(mfrow = c(length(unique_stands),2), mar = c(4,4,1,1))
for(s in 1:length(unique_stands)){
  trees_s <- which(data$newtree_stand_idxs==s)
  names <- sapply(1:N_max_newyears, function(n) paste0('states_new_plot[',s, ',', n,']'))
  util$plot_conditional_median_quantiles(samples, names, obs_xs = years_to_predict,
                                         bin_min = years_to_predict[1], bin_max = years_to_predict[length(years_to_predict)]+1, bin_delta = 1,
                                         ylab = 'Number of masting trees')
  util$plot_conditional_median_quantiles(samples, names, obs_xs = data$newprevsummer_temps[data$newtree_start_idxs[min(trees_s)]:data$newtree_end_idxs[min(trees_s)]],
                                         bin_min = 16, bin_max = 29, bin_delta = 0.5,
                                         ylab = 'Number of masting trees', xlab = 'Summer temperatures')
}


par(mfrow = c(length(unique_stands),2), mar = c(4,4,1,1))
for(s in 1:length(unique_stands)){
  trees_s <- which(data$newtree_stand_idxs==s)
  names <- sapply(1:N_max_newyears, function(n) paste0('seed_counts_new_plot[',s, ',', n,']'))
  util$plot_conditional_median_quantiles(samples, names, obs_xs = years_to_predict,
                                         bin_min = years_to_predict[1], bin_max = years_to_predict[length(years_to_predict)]+1, bin_delta = 1,
                                         ylab = 'Seed production')
  util$plot_conditional_median_quantiles(samples, names, obs_xs = data$newprevsummer_temps[data$newtree_start_idxs[min(trees_s)]:data$newtree_end_idxs[min(trees_s)]],
                                         bin_min = 16, bin_max = 29, bin_delta = 0.25,
                                         ylab = 'Seed production', xlab = 'Summer temperatures')
}


par(mfrow = c(1,2), mar = c(4,4,1,1))
names <- sapply(1:N_max_years, function(n) paste0('seed_counts_new_plot[',n,']'))
util$plot_conditional_median_quantiles(samples, names, obs_xs = first_year:last_year,
                                       bin_min = first_year, bin_max = last_year+1, bin_delta = 1,
                                       ylab = 'Seed production')
# util$plot_conditional_median_quantiles(samples, names, obs_xs = data$spring_temps[1:N_max_years],
#                                        bin_min = 5, bin_max = 10, bin_delta = 0.1,
#                                        ylab = 'Seed production')
# util$plot_conditional_median_quantiles(samples, names, obs_xs = data$prevsummer_temps[1:N_max_years],
#                                        bin_min = 15, bin_max = 21, bin_delta = 0.25,
#                                        ylab = 'Seed production')
util$plot_conditional_median_quantiles(samples, names, obs_xs = data$gdd_lastfrost[1:N_max_years],
                                       bin_min = 0, bin_max = 25, bin_delta = 1,
                                       ylab = 'Seed production', xlab = 'Frost risk')
