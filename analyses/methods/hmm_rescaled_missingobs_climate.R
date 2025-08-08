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
raw_data <- read.csv(file.path(wd, 'data',  'ebms', 'Beech_tree-ring_masting_data.csv'))
clim_data <- readRDS(file.path(wd, 'data',  'ebms', 'era5land_sitesextract.rds'))

length(unique(paste0(raw_data$site.ID,raw_data$tree.ID))) # 57 unique trees
length(unique(raw_data$site.ID)) # 7 unique sites
raw_data$uniqueID <- paste0(raw_data$site.ID, "_", raw_data$tree.ID)
raw_data <- na.omit(raw_data[,c('site.ID', 'uniqueID', 'year', 'seeds')])

# Sizes
uniq_tree_ids <- unique(raw_data$uniqueID)
N_trees <- length(uniq_tree_ids)
first_year <- min(raw_data$year)
last_year <- max(raw_data$year)
max_years <- first_year:max(raw_data$year)
N_max_years <- length(max_years)

# Format data into ragged arrays
seed_counts <- c()
years <- c()
prevsummer_temps <-c()
spring_temps <-c()
gdd_lastfrost <-c()
N_years <- c()

idx <- 1
tree_start_idxs <- c()
tree_end_idxs <- c()

for (tid in uniq_tree_ids) {
  
  raw_data_tree <- raw_data[raw_data$uniqueID == tid,]
  
  years_tree <- raw_data_tree$year-first_year+1
  years <- c(years, years_tree)
  
  N_years_tree <- length(years_tree)
  N_years <- c(N_years, N_years_tree)
  
  seed_count_tree <- raw_data_tree$seeds
  seed_counts <- c(seed_counts, seed_count_tree)
  
  prevsummer_temps_tree <- clim_data[clim_data$site.ID == raw_data_tree$site.ID[1] & clim_data$year %in% (max_years-1), 'meantmax_ja'] # all years, even those unobserved
  prevsummer_temps <- c(prevsummer_temps, prevsummer_temps_tree)
  
  spring_temps_tree <- clim_data[clim_data$site.ID == raw_data_tree$site.ID[1] & clim_data$year %in% (max_years-1), 'meantmean_am'] # all years, even those unobserved
  spring_temps <- c(spring_temps, spring_temps_tree)
  
  gdd_lastfrost_tree <- clim_data[clim_data$site.ID == raw_data_tree$site.ID[1] & clim_data$year %in% (max_years-1), 'gdd_b5_tolastfrost']/10 # all years, even those unobserved, in *10degC
  gdd_lastfrost <- c(gdd_lastfrost, gdd_lastfrost_tree)
  
  tree_start_idxs <- c(tree_start_idxs, idx)
  idx <- idx + N_years_tree
  tree_end_idxs <- c(tree_end_idxs, idx - 1)
  
}

# Format data
N <- length(years)
data <- mget(c('N', 'N_trees', 'N_max_years',
               'N_years', 'tree_start_idxs', 'tree_end_idxs',
               'seed_counts', 'years', 
               'prevsummer_temps', 'spring_temps', 'gdd_lastfrost'))

# Posterior quantification
fit <- stan(file= file.path(wd, "stan", "model2_treelevel_zinb_missing_rescaled_summertempfull_springfrost_springtemp.stan"),
            data=data, seed=5838299, chain = 4, cores = 4,
            warmup=1000, iter=2024, refresh=100)

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('lambda1', 'theta1', 'psi1',
                                         'lambda20', 'psi2', 
                                         'beta_lambda2_frost', 'beta_lambda2_spring',
                                         'rho0', 'beta_nm_m',
                                         'tau_nm_m0', 'tau_m_nm0'))
util$check_all_expectand_diagnostics(base_samples)

# Retrodictive check
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
                         xlab="Seed counts")

par(mfrow=c(3, 5), mar = c(4,4,2,2))
for (t in 1:data$N_trees) {
  
  idxs_tree <-(1+N_max_years*(t-1)):(N_max_years*t)
  observed_idxs_tree <- tree_start_idxs[t]:tree_end_idxs[t]
  observed_flags <- years[tree_start_idxs[t]:tree_end_idxs[t]]
  
  names <- sapply(idxs_tree,
                  function(n) paste0('seed_counts_pred[', n, ']'))
  xlab="Year"
  xticklabs=first_year:last_year
  ylab="Seed Counts"
  display_ylim=c(0, 400)
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
  
  plot(1, type="n", main=main,
       xlim=c(bin_min, bin_max), xlab=xlab, xaxt="n",
       ylim=display_ylim, ylab=ylab)
  axis(1, at=1:N, labels=xticklabs)    
  
  for(n in 1:N){
    
    idplot <- c(n*2-1, n*2)
    polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
            c(plot_quantiles[1,idplot], rev(plot_quantiles[9,idplot])),
            col = ifelse(n%in%observed_flags, util$c_light, '#b9d9b9'), border = NA)
    polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
            c(plot_quantiles[2,idplot], rev(plot_quantiles[8,idplot])),
            col = ifelse(n%in%observed_flags, util$c_light_highlight, '#96c796'), border = NA)
    polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
            c(plot_quantiles[3,idplot], rev(plot_quantiles[7,idplot])),
            col = ifelse(n%in%observed_flags, util$c_mid, '#72b472'), border = NA)
    polygon(c(plot_xs[idplot], rev(plot_xs[idplot])),
            c(plot_quantiles[4,idplot], rev(plot_quantiles[6,idplot])),
            col = ifelse(n%in%observed_flags, util$c_mid_highlight, '#50a250'), border = NA)
    
    lines(plot_xs[idplot], plot_quantiles[5, idplot],
          col= ifelse(n%in%observed_flags,util$c_dark, '#278f27'), lwd=2)
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


names <- sapply(1:(N_trees*N_max_years), function(n) paste0('tau_nm_m_pred[',n,']'))
util$plot_conditional_median_quantiles(samples, names, obs_xs = rep(first_year:last_year, N_trees),
                                       bin_min = first_year, bin_max = last_year+1, bin_delta = 1,
                                       ylab = 'Probability of masting')


# Posterior inferences
par(mfrow=c(2, 5))

util$plot_expectand_pushforward(samples[['lambda1']], 25,
                                display_name="lambda1")

util$plot_expectand_pushforward(samples[['psi1']], 25,
                                display_name="psi1")

util$plot_expectand_pushforward(samples[['theta1']], 25,
                                display_name="theta1")

util$plot_expectand_pushforward(samples[['lambda20']], 25,
                                display_name="lambda20")

util$plot_expectand_pushforward(samples[['beta_lambda2_frost']], 25,
                                display_name="beta_lambda2_frost")

util$plot_expectand_pushforward(samples[['beta_lambda2_spring']], 25,
                                display_name="beta_lambda2_spring")

util$plot_expectand_pushforward(samples[['psi2']], 25,
                                display_name="psi2")


util$plot_expectand_pushforward(samples[['rho0']], 
                                24, flim=c(0, 1.1),
                                display_name="rho0",)


util$plot_expectand_pushforward(samples[['tau_nm_m0']], 
                                25, flim=c(0, 1),
                                display_name="tau_nm_m0",)

util$plot_expectand_pushforward(samples[['beta_nm_m']], 25,
                                display_name="beta_nm_m")

util$plot_expectand_pushforward(samples[['tau_m_nm0']], 
                                25, flim=c(0, 1),
                                display_name="tau_m_nm0")


logit <- function(p){
  return(log(p/(1-p)))
}

invlogit <- function(r){
  return(1/(1+exp(-r)))
}

temp_summer <- seq(10,25,0.1)
tau_nm_m <- sapply(temp_summer, function(x) invlogit(logit(0.3)+0.4*(x-15)))
par(mfrow=c(1, 1))
plot(tau_nm_m ~ temp_summer, type = 'l')


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

library(ggplot2)
ggplot() +
  geom_line(data = tau_nm_m_df, aes(x = temp_prevsummer, y = X10.), color = 'black', linetype = 'dashed') +
  geom_line(data = tau_nm_m_df,aes(x = temp_prevsummer, y = X90.), color = 'black', linetype = 'dashed') +
  geom_line(data = tau_nm_m_df,aes(x = temp_prevsummer, y = X50.), color = 'black', linewidth = 1.5) +
  geom_line(data = tau_nm_m_df,aes(x = temp_prevsummer, y = X50.), color = 'white', linewidth = 0.6) +
  geom_boxplot(aes(x = data$prevsummer_temps, y = -0), inherit.aes = FALSE, width = 0.05) +
  theme_classic() +
  labs(y = 'Probability of transtion from non-masting to masting', x = 'Average previous summer max. temperature (degC)')


gdd0 <- 15
lambda2_df <- data.frame()
for(gdd_lastfrost in seq(0,35,0.1)){
  lambda2 <- matrix(nrow = nrow(samples[['lambda20']]), ncol = ncol(samples[['lambda20']]))
  for(c in 1:nrow(lambda2)){
    lambda2[c,] <- sapply(1:ncol(lambda2), function(i) exp(log(samples[['lambda20']][c,i]) +samples[['beta_lambda2_frost']][c,i]*(gdd_lastfrost-gdd0)))
  }
  lambda2_df <- rbind(lambda2_df,
                      data.frame(gdd_lastfrost, t(util$ensemble_mcmc_quantile_est(lambda2, probs = c(0.1,0.5,0.9)))))
}
ggplot() +
  geom_line(data = lambda2_df, aes(x = gdd_lastfrost, y = X10.), color = 'black', linetype = 'dashed') +
  geom_line(data = lambda2_df,aes(x = gdd_lastfrost, y = X90.), color = 'black', linetype = 'dashed') +
  geom_line(data = lambda2_df,aes(x = gdd_lastfrost, y = X50.), color = 'black', linewidth = 1.5) +
  geom_line(data = lambda2_df,aes(x = gdd_lastfrost, y = X50.), color = 'white', linewidth = 0.6) +
  geom_boxplot(aes(x = data$gdd_lastfrost, y = 120), inherit.aes = FALSE, width = 5) +
  theme_classic() +
  labs(y = 'Masting intensity', x = 'GDD until the last frost day (x10degC)')




gdd0 <- 15
lambda2_df <- data.frame()
for(gdd_lastfrost in seq(0,35,0.1)){
  lambda2 <- matrix(nrow = nrow(samples[['lambda20']]), ncol = ncol(samples[['lambda20']]))
  for(c in 1:nrow(lambda2)){
    lambda2[c,] <- sapply(1:ncol(lambda2), function(i) exp(log(samples[['lambda20']][c,i]) + (-0.005) *(gdd_lastfrost-gdd0)))
  }
  lambda2_df <- rbind(lambda2_df,
                      data.frame(gdd_lastfrost, t(util$ensemble_mcmc_quantile_est(lambda2, probs = c(0.1,0.5,0.9)))))
}
ggplot() +
  geom_line(data = lambda2_df, aes(x = gdd_lastfrost, y = X10.), color = 'black', linetype = 'dashed') +
  geom_line(data = lambda2_df,aes(x = gdd_lastfrost, y = X90.), color = 'black', linetype = 'dashed') +
  geom_line(data = lambda2_df,aes(x = gdd_lastfrost, y = X50.), color = 'black', linewidth = 1.5) +
  geom_line(data = lambda2_df,aes(x = gdd_lastfrost, y = X50.), color = 'white', linewidth = 0.6) +
  geom_boxplot(aes(x = data$gdd_lastfrost, y = 120), inherit.aes = FALSE, width = 5) +
  theme_classic() +
  labs(y = 'Masting intensity', x = 'GDD until the last frost day (x10degC)')

