
N_newtrees <- 100
newpreds <- data.frame()
for(y in N_max_years){
  newstates_y <- c()
  for(t in 1:N_newtrees){
    i <- y + (N_max_years)*(t-1)
    states_y_t <- samples[['states_new[1]']]-1
    newstates_y <- rbind(newstates_y, states_y_t)
  }
  
  colSums(newstates_y)
  
  newpreds <- rbind(
    newpreds,
    data.frame(y = y, ntrees_masting = )
  )
  
  
  
}


probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[['states_new_plot[1]']], probs)
}
quantiles <- sapply(1:N, calc)

par(mfrow = c(1,2), mar = c(4,4,1,1))
names <- sapply(1:N_max_years, function(n) paste0('states_new_plot[',n,']'))
util$plot_conditional_median_quantiles(samples, names, obs_xs = first_year:last_year,
                                       bin_min = first_year, bin_max = last_year+1, bin_delta = 1,
                                       ylab = 'Number of masting trees')
util$plot_conditional_median_quantiles(samples, names, obs_xs = data$prevsummer_temps[1:N_max_years],
                                       bin_min = 16, bin_max = 21, bin_delta = 0.25,
                                       ylab = 'Number of masting trees', xlab = 'Summer temperatures')

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
