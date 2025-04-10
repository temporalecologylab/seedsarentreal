


# Setup
rm(list = ls())
library(rstan)
library(ggplot2)
rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains

wd <- "~/projects/seedsarentreal/analyses/methods"
util <- new.env()
source(file.path(wd, "mcmc_analysis_tools_rstan.R"), local=util) 
source(file.path(wd, "mcmc_visualization_tools.R"), local=util) 

# Data exploration
raw_data <- read.csv(file.path(wd, 'methods', 'data',  'comast', 'individual_seed_production.csv'))

raw_data_filt <- raw_data[raw_data$general_method == 'TRAP' & raw_data$species_name == 'Acer_saccharum',]
raw_data_filt <- raw_data_filt[raw_data_filt$site_name == 'AEC',] # HBR site contains non integer seed counts....
raw_data_filt <- raw_data_filt[raw_data_filt$trap %in% c("958", "959", "960"),] # HBR site contains non integer seed counts....

ggplot(data = raw_data_filt) + 
  geom_line(aes(x = year, y = count, group = paste0(plot, trap)),
            linewidth = 0.1, alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(data = raw_data_filt) + 
  facet_wrap(~ trap) +
  geom_histogram(aes(x = count), bins = 50,
            linewidth = 0.1, alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Sizes
raw_data_filt$trapID <- paste0(raw_data_filt$plot, raw_data_filt$trap)
uniq_trap_ids <- unique(raw_data_filt$trapID)
N_traps <- length(uniq_trap_ids)

# Format data into ragged arrays
seed_counts <- c()
sizes <- c()
years <- c()
N_years <- c()

idx <- 1
trap_start_idxs <- c()
trap_end_idxs <- c()

for (tid in uniq_trap_ids) {
  raw_data_trap <- raw_data_filt[raw_data_filt$trapID == tid,]
  
  years_trap <- raw_data_trap$year
  N_years_trap <- length(years_trap)
  
  size_trap <- raw_data_trap$trap_area_m2[1]
  seed_count_trap <- raw_data_trap$count
  
  seed_counts <- c(seed_counts, seed_count_trap)
  sizes <- c(sizes, size_trap)
  years <- c(years, years_trap)
  N_years <- c(N_years, N_years_trap)
  
  trap_start_idxs <- c(trap_start_idxs, idx)
  idx <- idx + N_years_trap
  trap_end_idxs <- c(trap_end_idxs, idx - 1)
}

# Cross check sizes
length(seed_counts)
length(years)
length(sizes)
length(N_years)
length(trap_start_idxs)
length(trap_end_idxs)


# Plot all trees using new data format
par(mfrow=c(2, 3))

for (t in 1:N_traps) {
  idxs <- trap_start_idxs[t]:trap_end_idxs[t]
  
  years_trap <- years[idxs]
  seed_counts_trap <- seed_counts[idxs]
  
  plot(years_trap, seed_counts_trap, pch=16, cex=1.0, 
       xlab="Year", ylab="Seed Counts", ylim=c(0, 100),
       main=paste("Trap", uniq_trap_ids[t]))
}

# Format data
N <- length(years)

data <- mget(c('N', 'N_traps', 'sizes',
               'seed_counts', 'years', 'N_years',
               'trap_start_idxs', 'trap_end_idxs'))

# Posterior quantification
fit <- stan(file= file.path(wd, "stan", "model2_inv.stan"),
            data=data, seed=5838299,
            warmup=1000, iter=2024, refresh=0)

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('lambda1', 'lambda2', 
                                         'psi1', 'rho0',
                                         'tau_nm_m', 'tau_m_nm'))
util$check_all_expectand_diagnostics(base_samples)

# Retrodicive Check
par(mfrow=c(1, 1))

util$plot_hist_quantiles(samples, 'seed_count_pred', 0, 100, 1,
                         baseline_values=data$seed_counts, 
                         xlab="Seed Counts")

par(mfrow=c(2, 3))

for (t in 1:data$N_traps) {
  idxs <- trap_start_idxs[t]:trap_end_idxs[t]
  names <- sapply(idxs,
                  function(n) paste0('seed_count_pred[', n, ']'))
  util$plot_disc_pushforward_quantiles(samples, names, 
                                       baseline_values=data$seed_counts[idxs],
                                       xlab="Year",
                                       xticklab=data$years[idxs],
                                       ylab="Seed Counts",
                                       display_ylim=c(0, 100),
                                       main=paste("Trap", uniq_trap_ids[t]))
}

par(mfrow=c(2, 3))

for (n in c(10, 20, 30, 40, 50, 60)) {
  tid <- tail(which(data$trap_start_idxs <= n), 1)
  name <- paste0('seed_count_pred[', n, ']')
  util$plot_expectand_pushforward(samples[[name]],
                                  50, flim=c(0, 300),
                                  baseline=data$seed_count[n],
                                  display_name="Predicted Seed Counts",
                                  main=paste0("Trap ", uniq_trap_ids[tid],
                                              ", Year ", data$years[n]))
}

# Posterior inferences
par(mfrow=c(2, 3))

util$plot_expectand_pushforward(samples[['lambda1']], 25,
                                display_name="lambda1")
xs <- seq(0, 100, 0.01)
ys <- dnorm(xs, 0, 500 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['lambda2']], 25,
                                display_name="lambda2")
xs <- seq(000, 500, 0.5)
ys <- dnorm(xs, 0, 500 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['psi1']], 25,
                                display_name="psi1")
xs <- seq(0, 5, 0.1)
ys <- dnorm(xs, 0, 5 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['rho0']], 
                                22, flim=c(0, 1.1),
                                display_name="rho0",)
xs <- seq(0, 1.1, 0.01)
ys <- rep(1, length(xs))
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['tau_nm_m']], 
                                25, flim=c(0, 1),
                                display_name="tau_nm_m",)
xs <- seq(0, 1, 0.01)
ys <- rep(1, length(xs))
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['tau_m_nm']], 
                                25, flim=c(0, 1),
                                display_name="tau_m_nm",)
xs <- seq(0, 1, 0.01)
ys <- rep(1, length(xs))
lines(xs, ys, lwd=2, col=util$c_light)




par(mfrow=c(2, 3))

for (t in 1:data$N_traps) {
  idxs <- trap_start_idxs[t]:trap_end_idxs[t]
  names <- sapply(idxs,
                  function(n) paste0('p_masting[', n, ']'))
  util$plot_disc_pushforward_quantiles(samples, names,
                                       xlab="Year",
                                       xticklab=data$years[idxs],
                                       ylab="Posterior Masting Probability",
                                       display_ylim=c(-0.01, 1),
                                       main=paste("Trap", uniq_trap_ids[t]))
}
