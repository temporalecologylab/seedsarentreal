######################################################################
#
# Setup
#
######################################################################

par(family="serif", las=1, bty="l",
    cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 1))

library(rstan)
rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

######################################################################
#
# Data Exploration
#
######################################################################

raw_data <- read.csv('data/TSHE.csv')

# Visualize data
par(mfrow=c(2, 1))

util$plot_line_hist(raw_data$allseeds, 0, 270, 10,
                    xlab="Seed Count", main="Trap 1")

plot(raw_data$year, raw_data$allseeds, 
     pch=16, cex=1.0, xlab="Year", ylab="Seed Count")

# Sizes
uniq_trap_ids <- unique(raw_data$trap)
N_traps <- length(uniq_trap_ids)

# Format data into ragged arrayss
seed_counts <- c()
sizes <- c()
years <- c()
N_years <- c()

idx <- 1
trap_start_idxs <- c()
trap_end_idxs <- c()

for (tid in uniq_trap_ids) {
  raw_data_trap <- raw_data[raw_data$trap == tid,]
  
  years_trap <- raw_data_trap$year
  N_years_trap <- length(years_trap)
  
  size_trap <- raw_data_trap$size[1]
  seed_count_trap <- raw_data_trap$allseeds
 
  seed_counts <- c(seed_counts, seed_count_trap)
  sizes <- c(sizes, size_trap)
  years <- c(years, years_trap)
  N_years <- c(N_years, N_years_trap)
  
  trap_start_idxs <- c(trap_start_idxs, idx)
  idx <- idx + N_years_trap
  trap_end_idxs <- c(trap_end_idxs, idx - 1)
}

# Cross check sizes
N_traps

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
       xlab="Year", ylab="Seed Counts", ylim=c(0, 200),
       main=paste("Trap", uniq_trap_ids[t]))
}

# Format data
N <- length(years)

data <- mget(c('N', 'N_traps', 'sizes',
               'seed_counts', 'years', 'N_years',
               'trap_start_idxs', 'trap_end_idxs'))

######################################################################
#
# Model 1
#
######################################################################

l <- log(50)
u <- log(200)

0.5 * (u + l)
0.5 * (u - l) / 2.32

# Posterior Quantification
fit <- stan(file='stan_programs/model1.stan',
            data=data, seed=5838299,
            warmup=1000, iter=2024, refresh=0)

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('lambda1', 'lambda2', 
                                         'psi1', 'rho'))
util$check_all_expectand_diagnostics(base_samples)

# Retrodicive Check
par(mfrow=c(1, 1))

util$plot_hist_quantiles(samples, 'seed_count_pred', 0, 300, 25,
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
                                       display_ylim=c(0, 270),
                                       main=paste("Trap", uniq_trap_ids[t]))
}

par(mfrow=c(2, 3))

for (n in c(10, 20, 30, 40, 50, 60)) {
  tid <- tail(which(data$trap_start_idxs <= n), 1)
  name <- paste0('seed_count_pred[', n, ']')
  util$plot_expectand_pushforward(samples[[name]],
                                  50, flim=c(0, 250),
                                  baseline=data$seed_count[n],
                                  display_name="Predicted Seed Counts",
                                  main=paste0("Trap ", uniq_trap_ids[tid],
                                             ", Year ", data$years[n]))
}

# Posterior inferences
par(mfrow=c(2, 2))

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

util$plot_expectand_pushforward(samples[['rho']], 
                                22, flim=c(0, 1.1),
                                display_name="rho",)
xs <- seq(0, 1.1, 0.01)
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
                                       display_ylim=c(0, 1),
                                       main=paste("Trap", uniq_trap_ids[t]))
}

dnegbin <- function(x, mu, psi) {
  choose(x + 1 / psi - 1, x) *
  (mu * psi / (mu * psi + 1))**x * 
  (1 / (mu * psi + 1))**(1 / psi)
}

xs <- 0:300
ys <- dnegbin(xs, 60, 0.75)
plot(xs, ys)

######################################################################
#
# Model 2
#
######################################################################

# Posterior Quantification
fit <- stan(file='stan_programs/model2.stan',
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

util$plot_hist_quantiles(samples, 'seed_count_pred', 0, 300, 25,
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
                                       display_ylim=c(0, 270),
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