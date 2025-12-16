


wd <- "~/projects/seedsarentreal/analyses/methods"
setwd(wd)
util <- new.env()
source(file.path(wd, "mcmc_analysis_tools_rstan.R"), local=util) 
source(file.path(wd, "mcmc_visualization_tools.R"), local=util)

fit_full_const <- readRDS(file.path(wd, 'output', 'fit_25sept2025_fullmodel_constrainedslopes.rds'))
samples_full_const <- util$extract_expectand_vals(fit_full_const)
fit_biol <- readRDS(file = file.path(wd, 'output', 'fit_25sept2025_biologymodel.rds'))
samples_biol <- util$extract_expectand_vals(fit_biol)

yearstoplot <- 1960:2100-1949

par(mfrow = c(3,1), mar = c(4,4,1,1))
s <- 1
names <- sapply(yearstoplot, function(n) paste0('seed_counts_new_brkd_plot[',s, ',', n,']'))
util$plot_conditional_median_quantiles(samples_full_const, names, obs_xs = data$newprevsummer_temps[yearstoplot],
                                       bin_min = 18, bin_max = 27, bin_delta = 0.5,
                                       ylab = 'Number of masting trees', xlab = '',
                                       main = 'Breakdown predictions', display_ylim = c(700, 1800))

s <- 1
trees_s <- which(data$newtree_stand_idxs==s)
names <- sapply(yearstoplot, function(n) paste0('seed_counts_new_plot[',s, ',', n,']'))
util$plot_conditional_median_quantiles(samples_full_const, names, obs_xs = data$newprevsummer_temps[yearstoplot],
                                       bin_min = 18, bin_max = 27, bin_delta = 0.5,
                                       ylab = 'Number of masting trees', xlab = '',
                                       main = 'Effect of summer temperature on both transitions', display_ylim = c(700, 1800))

s <- 1
trees_s <- which(data$newtree_stand_idxs==s)
names <- sapply(yearstoplot, function(n) paste0('seed_counts_new_plot[',s, ',', n,']'))
util$plot_conditional_median_quantiles(samples_biol, names, obs_xs = data$newprevsummer_temps[yearstoplot],
                                       bin_min = 18, bin_max = 27, bin_delta = 0.5,
                                       ylab = 'Number of masting trees', xlab = 'Summer temperatures',
                                       main = 'Biological model', display_ylim = c(700, 1800))


names <- sapply(yearstoplot, function(n) paste0('seed_counts_new_plot[',s, ',', n,']'))
predic_biol <- data.frame(
  q5 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_biol[[n]], c(0.1))}),
  q50 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_biol[[n]], c(0.5))}),
  q95 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_biol[[n]], c(0.9))}),
  prevsummertemp = data$newprevsummer_temps[yearstoplot]
)

names <- sapply(yearstoplot, function(n) paste0('seed_counts_new_plot[',s, ',', n,']'))
predic_full <- data.frame(
  q5 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.1))}),
  q50 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.5))}),
  q95 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.9))}),
  prevsummertemp = data$newprevsummer_temps[yearstoplot]
)

names <- sapply(yearstoplot, function(n) paste0('seed_counts_new_brkd_plot[',s, ',', n,']'))
predic_brkd <- data.frame(
  q5 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.1))}),
  q50 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.5))}),
  q95 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.9))}),
  prevsummertemp = data$newprevsummer_temps[yearstoplot]
)

ggplot() +
  geom_line(data = predic_biol, aes(x = prevsummertemp, y = q5), linetype = 'dashed', col = 'darkblue') +
  geom_line(data = predic_biol, aes(x = prevsummertemp, y = q50), col = 'darkblue') +
  geom_line(data = predic_biol, aes(x = prevsummertemp, y = q95), linetype = 'dashed', col = 'darkblue') +
  geom_line(data = predic_full, aes(x = prevsummertemp, y = q5), linetype = 'dashed', col = 'darkgreen') +
  geom_line(data = predic_full, aes(x = prevsummertemp, y = q50), col = 'darkgreen') +
  geom_line(data = predic_full, aes(x = prevsummertemp, y = q95), linetype = 'dashed', col = 'darkgreen') +
  geom_line(data = predic_brkd, aes(x = prevsummertemp, y = q5), linetype = 'dashed', col = 'darkred') +
  geom_line(data = predic_brkd, aes(x = prevsummertemp, y = q50), col = 'darkred') +
  geom_line(data = predic_brkd, aes(x = prevsummertemp, y = q95), linetype = 'dashed', col = 'darkred') +
  theme_classic() +
  labs(y = 'Seed production at the stand-level (100 trees)')


names <- sapply(yearstoplot, function(n) paste0('states_new_plot[',s, ',', n,']'))
predic_biol <- data.frame(
  q5 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_biol[[n]], c(0.1))}),
  q50 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_biol[[n]], c(0.5))}),
  q95 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_biol[[n]], c(0.9))}),
  prevsummertemp = data$newprevsummer_temps[yearstoplot]
)

names <- sapply(yearstoplot, function(n) paste0('states_new_plot[',s, ',', n,']'))
predic_full <- data.frame(
  q5 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.1))}),
  q50 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.5))}),
  q95 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.9))}),
  prevsummertemp = data$newprevsummer_temps[yearstoplot]
)

names <- sapply(yearstoplot, function(n) paste0('states_new_brkd_plot[',s, ',', n,']'))
predic_brkd <- data.frame(
  q5 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.1))}),
  q50 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.5))}),
  q95 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.9))}),
  prevsummertemp = data$newprevsummer_temps[yearstoplot]
)


ggplot() +
  # geom_line(data = predic_biol, aes(x = prevsummertemp, y = q5), linetype = 'dashed', col = 'darkblue') +
  # geom_line(data = predic_biol, aes(x = prevsummertemp, y = q50), col = 'darkblue') +
  # geom_line(data = predic_biol, aes(x = prevsummertemp, y = q95), linetype = 'dashed', col = 'darkblue') +
  geom_line(data = predic_full, aes(x = prevsummertemp, y = q5), linetype = 'dashed', col = 'darkgreen') +
  geom_line(data = predic_full, aes(x = prevsummertemp, y = q50), col = 'darkgreen', linewidth = 1.2) +
  geom_line(data = predic_full, aes(x = prevsummertemp, y = q95), linetype = 'dashed', col = 'darkgreen') +
  geom_line(data = predic_brkd, aes(x = prevsummertemp, y = q5), linetype = 'dashed', col = 'darkred') +
  geom_line(data = predic_brkd, aes(x = prevsummertemp, y = q50), col = 'darkred', linewidth = 1.2) +
  geom_line(data = predic_brkd, aes(x = prevsummertemp, y = q95), linetype = 'dashed', col = 'darkred') +
  theme_classic() +
  labs(y = 'Number of masting trees (out of 100 trees)', x = 'Average previous summer max. temperature (degC)')



