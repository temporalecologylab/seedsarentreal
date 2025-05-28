
raw_data_filt <- raw_data[raw_data$general_method == 'tree',]

raw_data_filt <- raw_data[raw_data$general_method == 'PARTIALCONECOUNT' & raw_data$species_name == 'Abies_amabilis',]
raw_data_filt$uniqueID <- paste0(raw_data_filt$plot, raw_data_filt$plant_ID)
raw_data_filt <- raw_data_filt[raw_data_filt$year %in% c(1962:2016) & !is.na(raw_data_filt$count),]
obscount <- data.frame(table(raw_data_filt$uniqueID))
tokeep <- obscount[obscount$Freq == length(c(1962:2016)), "Var1"]
length(tokeep)
raw_data_filt <- raw_data_filt[raw_data_filt$uniqueID %in% tokeep,]

ggplot(data = raw_data_filt) + 
  geom_line(aes(x = year, y = count, group = uniqueID),
            linewidth = 0.1, alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank())


# raw_data_filt <- raw_data_filt[raw_data_filt$plant_ID %in% c('CNCT_51ABAM19', 'CNCT_01ABAM11', 'CNCT_51ABAM12'),]


# Sizes
uniq_tree_ids <- unique(raw_data_filt$uniqueID)
N_trees<- length(uniq_tree_ids)

# Format data into ragged arrays
seed_counts <- c()
sizes <- c()
years <- c()
N_years <- c()

idx <- 1
tree_start_idxs <- c()
tree_end_idxs <- c()

for (tid in uniq_tree_ids) {
  raw_data_tree <- raw_data_filt[raw_data_filt$uniqueID == tid,]
  
  years_tree <- raw_data_tree$year
  N_years_tree <- length(years_tree)
  
  size_tree <- 1
  seed_count_tree <- raw_data_tree$count
  
  seed_counts <- c(seed_counts, seed_count_tree)
  sizes <- c(sizes, size_tree)
  years <- c(years, years_tree)
  N_years <- c(N_years, N_years_tree)
  
  tree_start_idxs <- c(tree_start_idxs, idx)
  idx <- idx + N_years_tree
  tree_end_idxs <- c(tree_end_idxs, idx - 1)
}

# Cross check sizes
length(seed_counts)
length(years)
length(sizes)
length(N_years)
length(tree_start_idxs)
length(tree_end_idxs)

N_years <- as.array(N_years)
sizes <- as.array(sizes)
tree_start_idxs <- as.array(tree_start_idxs)
tree_end_idxs <- as.array(tree_end_idxs)

# Format data
N <- length(years)
data <- mget(c('N', 'N_trees', 'sizes',
               'seed_counts', 'years', 'N_years',
               'tree_start_idxs', 'tree_end_idxs'))

fit <- stan(file= file.path(wd, "stan", "model2_treelevel_diffdist.stan"),
            data=data, seed=5838299, cores =4,
            warmup=1000, iter=4000, refresh=100)

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('rho0',
                                         'tau_nm_m', 'tau_m_nm',
                                         'lambda1', 'psi1',
                                         'lambda2'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

# Retrodicive Check
par(mfrow=c(1, 1))

util$plot_hist_quantiles(samples, 'seed_count_pred', 0, 500, 5,
                         baseline_values=data$seed_counts, 
                         xlab="Seed Counts")

par(mfrow=c(3, 3))

for (t in sample(1:data$N_trees, 9)) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  names <- sapply(idxs,
                  function(n) paste0('seed_count_pred[', n, ']'))
  util$plot_disc_pushforward_quantiles(samples, names, 
                                       baseline_values=data$seed_counts[idxs],
                                       xlab="Year",
                                       xticklab=data$years[idxs],
                                       ylab="Seed Counts",
                                       display_ylim=c(0, 300),
                                       main=paste("tree", uniq_tree_ids[t]))
}

par(mfrow=c(2, 3))

for (n in c(10, 20, 30, 40, 50, 60)) {
  tid <- tail(which(data$tree_start_idxs <= n), 1)
  name <- paste0('seed_count_pred[', n, ']')
  util$plot_expectand_pushforward(samples[[name]],
                                  50, flim=c(0, 300),
                                  baseline=data$seed_count[n],
                                  display_name="Predicted Seed Counts",
                                  main=paste0("tree ", uniq_tree_ids[tid],
                                              ", Year ", data$years[n]))
}


# Posterior inferences
par(mfrow=c(1, 4))

util$plot_expectand_pushforward(samples[['lambda1']], 25,
                                display_name="lambda1")
xs <- seq(0, 100, 0.01)
ys <- dnorm(xs, 0, 500 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['psi1']], 25,
                                display_name="psi1", flim = c(0,8))
xs <- seq(0, 100, 0.1)
ys <- dnorm(xs, 0, 5 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['lambda2']], 25,
                                display_name="lambda2")
xs <- seq(0, 800, 0.5)
ys <- dnorm(xs, 0, 500 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['psi2']], 25,
                                display_name="psi2", flim = c(0,8))
xs <- seq(0, 100, 0.1)
ys <- dnorm(xs, 0, 5 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

par(mfrow=c(1, 3))
util$plot_expectand_pushforward(samples[['rho0[1]']], 25,
                                display_name="rho", flim = c(0,1))
util$plot_expectand_pushforward(samples[['rho0[2]']], 25,
                                display_name="rho", flim = c(0,1))
util$plot_expectand_pushforward(samples[['rho0[3]']], 25,
                                display_name="rho", flim = c(0,1))

par(mfrow=c(1, 3))
util$plot_expectand_pushforward(samples[['lambda1[1]']], 25,
                                display_name="rho", flim = c(0,200))
util$plot_expectand_pushforward(samples[['lambda1[2]']], 25,
                                display_name="rho", flim = c(0,200))
util$plot_expectand_pushforward(samples[['lambda1[3]']], 25,
                                display_name="rho", flim = c(0,200))

par(mfrow=c(3, 2))
util$plot_expectand_pushforward(samples[['tau_nm_m[1]']], 25,
                                display_name="tau_nm_m", flim = c(0,1))
util$plot_expectand_pushforward(samples[['tau_m_nm[1]']], 25,
                                display_name="tau_m_nm", flim = c(0,1))
util$plot_expectand_pushforward(samples[['tau_nm_m[2]']], 25,
                                display_name="tau_nm_m", flim = c(0,1))
util$plot_expectand_pushforward(samples[['tau_m_nm[2]']], 25,
                                display_name="tau_m_nm", flim = c(0,1))
util$plot_expectand_pushforward(samples[['tau_nm_m[3]']], 25,
                                display_name="tau_nm_m", flim = c(0,1))
util$plot_expectand_pushforward(samples[['tau_m_nm[3]']], 25,
                                display_name="tau_m_nm", flim = c(0,1))

