
logseed0 <- log(100)
beta.prev <- -0.035
beta.summer <- 0.06
beta.frost <- -0.05

heuristic_simulations <- data.frame()
for(t in unique(raw_data$uniqueID)){
  s <- unique(raw_data[raw_data$uniqueID ==t, 'site.ID'])
  seeds0 <- runif(1,1,10) # seed production in 1979
  seed.n1 <- seeds0
  seeds <- c()
  for(y in 1980:2022){
    seedhat <- exp(logseed0 + beta.prev*seed.n1 + 
                     beta.summer * (clim_data[clim_data$site.ID == s & clim_data$year == y, 'meantmax_ja']-15) +
                     beta.frost * (clim_data[clim_data$site.ID == s & clim_data$year == y, 'gdd_b5_tolastfrost']/10-15))
    seed <- rpois(1, lambda = seedhat)
    seed.n1 <- seed
    seeds <- c(seeds,seed)
  }
  
  heuristic_simulations <- rbind(
    heuristic_simulations,
    data.frame(site.ID =s, uniqueID = t, year = 1980:2022, seeds = seeds)
  )
  
}

ggplot(data = heuristic_simulations) +
  geom_line(aes(x = year, y = seeds, group = uniqueID)) +
  theme_bw()


# Sizes
uniq_tree_ids <- unique(heuristic_simulations$uniqueID)
N_trees <- length(uniq_tree_ids)
first_year <- min(heuristic_simulations$year)
last_year <- max(heuristic_simulations$year)
max_years <- first_year:max(heuristic_simulations$year)
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
  
  raw_data_tree <- heuristic_simulations[heuristic_simulations$uniqueID == tid,]
  
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

# Format data for Stan
N <- length(years)
Nnew <- length(newyears)
data <- mget(c('N', 'N_trees', 'N_max_years',
               'N_years', 'tree_start_idxs', 'tree_end_idxs',
               'seed_counts', 'years', 
               'prevsummer_temps', 'spring_temps', 'gdd_lastfrost'
))


fit <- stan(file= file.path(wd, "stan", "paper_models", "model2_treelevel_zinb_missing_rescaled_summertemp_springfrost_springtemp_constrainedslopes.stan"),
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
                                         'tau_nm_m0', 'tau_m_m0'))
util$check_all_expectand_diagnostics(base_samples)


# Posterior inferences
par(mfrow=c(2, 6))

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

util$plot_expectand_pushforward(samples[['tau_m_m0']], 
                                25, flim=c(0, 1),
                                display_name="tau_m_m0")


par(mfrow=c(1, 1), mar = c(4,4,2,2))
names <- sapply(1:data$N, function(n) paste0('seed_counts_pred[',n,']'))
util$plot_hist_quantiles(samples[names], 'seed_counts_pred', 0, 340, 20,
                         baseline_values=data$seed_counts, 
                         xlab="Seed counts")

for (t in 1:data$N_trees) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  names <- sapply(idxs,
                  function(n) paste0('seed_counts_pred[', n, ']'))
  util$plot_disc_pushforward_quantiles(samples, names, 
                                       baseline_values=data$seed_counts[idxs],
                                       xlab="Year",
                                       xticklab=data$years[idxs],
                                       ylab="Seed Counts",
                                       display_ylim=c(0, 270))
}

