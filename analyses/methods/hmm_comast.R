
# Setup
rm(list = ls())
options(stringsAsFactors = FALSE)
library(rstan)
library(ggplot2)
rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains

wd <- "~/projects/seedsarentreal/analyses/methods"
setwd(wd)
util <- new.env()
source(file.path(wd, "mcmc_analysis_tools_rstan.R"), local=util) 
source(file.path(wd, "mcmc_visualization_tools.R"), local=util) 

kippenberger <- c("#8B174DFF", "#AE2565FF", "#C1447EFF", "#D06C9BFF", "#DA9FB8FF", 
                  "#ADBE7CFF", "#8BA749FF", "#6E8537FF", "#4F5F28FF", "#343D1FFF")


# Data exploration
raw_data <- read.csv(file.path(wd, 'data',  'comast', 'individual_seed_production.csv'))

raw_data_filt <- raw_data[raw_data$general_method == 'TRAP',]
# View(data.frame(table(raw_data_filt$species_name)))
raw_data_filt <- raw_data[raw_data$general_method == 'TRAP' & raw_data$species_name == 'Fagus_grandifolia',]

raw_data_filt$uniqueID <- paste0(raw_data_filt$plot, raw_data_filt$trap)

raw_data_filt <- raw_data_filt[raw_data_filt$year %in% c(1995:2020),]
obscount <- data.frame(table(raw_data_filt$uniqueID))
trapstokeep <- obscount[obscount$Freq == 26, "Var1"]
raw_data_filt <- raw_data_filt[raw_data_filt$uniqueID %in% trapstokeep,]

s <- sample(unique(raw_data_filt$uniqueID), 1)
ggplot(data = raw_data_filt[raw_data_filt$uniqueID == "TF_Upper24",]) + 
  geom_line(aes(x = year, y = count, group = paste0(plot, trap)),
            linewidth = 0.1, alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(data = raw_data_filt) + 
  # facet_wrap(~ trap) +
  geom_histogram(aes(x = count), bins = 100,
                 linewidth = 0.1, alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank())

# raw_data_filt <- raw_data_filt[raw_data_filt$uniqueID == "TF_Upper24",]

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
  seed_count_trap <- round(raw_data_trap$count,0) # aarrrg
  
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

N_years <- as.array(N_years)
sizes <- as.array(sizes)
trap_start_idxs <- as.array(trap_start_idxs)
trap_end_idxs <- as.array(trap_end_idxs)

# Format data
N <- length(years)
data <- mget(c('N', 'N_traps', 'sizes',
               'seed_counts', 'years', 'N_years',
               'trap_start_idxs', 'trap_end_idxs'))

fit <- stan(file= file.path(wd, "stan", "model2_alt.stan"),
            data=data, seed=5838299,
            warmup=1000, iter=4000, refresh=100)


diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('lambda1', 'lambda2', 
                                         'psi1', 'rho0',
                                         'tau_nm_m', 'tau_m_nm'))
util$check_all_expectand_diagnostics(base_samples)



# Posterior inferences
par(mfrow=c(1, 4))

util$plot_expectand_pushforward(samples[['lambda1']], 25,
                                display_name="lambda1")
xs <- seq(0, 100, 0.01)
ys <- dnorm(xs, 0, 1 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['psi1']], 25,
                                display_name="psi1")
xs <- seq(0, 100, 0.1)
ys <- dnorm(xs, 0, 5 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['lambda2']], 25,
                                display_name="lambda2")
xs <- seq(0, 800, 0.5)
ys <- dnorm(xs, 0, 500 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)


util$plot_expectand_pushforward(samples[['psi2']], 25,
                                display_name="psi2")
xs <- seq(0, 50, 0.5)
ys <- dnorm(xs, 0, 5 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)



par(mfrow=c(1, 2))
hist(MASS::rnegbin(10000, mean(samples$lambda1), 1/mean(samples$psi1)), 
     main = "Non-masting", xlab = "Seed count", , breaks = 50)
 hist(rpois(10000, mean(samples$lambda2)), main = "Masting", xlab = "Seed count")
hist(MASS::rnegbin(10000, mean(samples$lambda2), 1/mean(samples$psi2)), 
     main = "Non-masting", xlab = "Seed count", breaks = 50)

par(mfrow=c(1, 3))
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

par(mfrow=c(3, 3))
for (n in seq(1,10,1)) {
  tid <- tail(which(data$tree_start_idxs <= n), 1)
  name <- paste0('seed_count_pred[', n, ']')
  util$plot_expectand_pushforward(samples[[name]],
                                  50, flim=c(0, 100),
                                  baseline=data$seed_count[n],
                                  display_name="Predicted Seed Counts",
                                  main=paste0("Tree ", uniq_trap_ids[tid],
                                              ", Year ", data$years[n]))
}


par(mfrow=c(1, 1))

for (t in 1:data$N_traps) {
  idxs <- trap_start_idxs[t]:trap_end_idxs[t]
  names <- sapply(idxs,
                  function(n) paste0('seed_count_pred[', n, ']'))
  util$plot_disc_pushforward_quantiles(samples, names, 
                                       baseline_values=data$seed_counts[idxs],
                                       xlab="Year",
                                       xticklab=data$years[idxs],
                                       ylab="Seed Counts",
                                       display_ylim=c(0, 20),
                                       main=paste("Tree", uniq_trap_ids[t]))
}




