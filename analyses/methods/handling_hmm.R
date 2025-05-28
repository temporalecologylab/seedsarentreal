
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

length(unique(paste0(raw_data$site.ID,raw_data$tree.ID))) # 57 unique trees
length(unique(raw_data$site.ID)) # 7 unique sites

raw_data$uniqueID <- paste0(raw_data$site.ID, "_", raw_data$tree.ID)

ggplot(data = raw_data) + 
  facet_wrap(~ uniqueID) +
  geom_line(aes(x = year, y = seeds, group = uniqueID, color = uniqueID),
            linewidth = 0.2, alpha = 1) +
  geom_vline(data = raw_data[is.na(raw_data$seeds),], aes(xintercept = year), color = "grey95",
             linewidth = 1.5, alpha = 1) +
  scale_color_manual(values = colorRampPalette(kippenberger)(length(unique(raw_data$uniqueID)))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        strip.text.x.top = element_blank(),
        axis.text = element_text(size = 7))

ggplot(data = raw_data) + 
  facet_wrap(~ uniqueID) +
  geom_histogram(aes(x = seeds, fill = uniqueID), bins = 20, alpha = 1, color = "white") +
  scale_fill_manual(values = colorRampPalette(kippenberger)(length(unique(raw_data$uniqueID)))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        strip.text.x.top = element_blank())

ggplot(data = raw_data) + 
  geom_histogram(aes(x = seeds), bins = 40, alpha = 1, color = "white") +
  scale_fill_manual(values = colorRampPalette(kippenberger)(length(unique(raw_data$uniqueID)))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        strip.text.x.top = element_blank())

# Missing data? Not mentioned in HP et al. 2025... nor in https://onlinelibrary.wiley.com/doi/10.1111/gcb.16730...

## as first brute-force approach, we keep only trees with all obs. between 1985 and 2020
length(1985:2020)
data.frame(table(raw_data[,c("uniqueID")]))
raw_data_filter <- na.omit(raw_data[raw_data$year %in% c(1985:2020),])
trees_nobs <- data.frame(table(raw_data_filter[,c("uniqueID")]))
names(trees_nobs) <- c("uniqueID", "nobs")
raw_data_filter <- raw_data_filter[raw_data_filter$uniqueID %in% trees_nobs[trees_nobs$nobs == length(1985:2020), "uniqueID"],]

ggplot(data = raw_data_filter) + 
  geom_line(aes(x = year, y = seeds, group = uniqueID, color = uniqueID),
            linewidth = 0.2, alpha = 1) +
  geom_vline(data = raw_data_filter[is.na(raw_data_filter$seeds),], aes(xintercept = year), color = "grey95",
             linewidth = 1.5, alpha = 1) +
  scale_color_manual(values = colorRampPalette(kippenberger)(length(unique(raw_data_filter$uniqueID)))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        strip.text.x.top = element_blank(),
        axis.text = element_text(size = 7))

ggplot(data = raw_data_filter) + 
  geom_histogram(aes(x = seeds), bins = 60, alpha = 1, color = "white", fill = "grey80") +
  scale_fill_manual(values = colorRampPalette(kippenberger)(length(unique(raw_data$uniqueID)))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        strip.text.x.top = element_blank())



# Sizes
uniq_tree_ids <- unique(raw_data_filter$uniqueID)
N_trees <- length(uniq_tree_ids)

# Format data into ragged arrays
seed_counts <- c()
sizes <- c()
years <- c()
N_years <- c()

idx <- 1
tree_start_idxs <- c()
tree_end_idxs <- c()

for (tid in uniq_tree_ids) {
  raw_data_tree <- raw_data_filter[raw_data_filter$uniqueID == tid,]
  
  years_tree <- raw_data_tree$year
  N_years_tree <- length(years_tree)
  
  size_tree <- 1 # the ground below each tree was searched for seeds, no info on area
  seed_count_tree <- raw_data_tree$seeds
  
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


# Plot all trees using new data format
par(mfrow=c(3, 3))

for (t in 1:N_trees) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  
  years_tree <- years[idxs]
  seed_counts_tree <- seed_counts[idxs]
  
  plot(years_tree, seed_counts_tree, pch=16, cex=1.0, 
       xlab="Year", ylab="Seed Counts", ylim=c(0, 100),
       main=paste("Trap", uniq_tree_ids[t]))
}

# Format data
N <- length(years)

data <- mget(c('N', 'N_trees', 'sizes',
               'seed_counts', 'years', 'N_years',
               'tree_start_idxs', 'tree_end_idxs'))

# Posterior quantification
fit <- stan(file= file.path(wd, "stan", "model2_treelevel_alt.stan"),
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
par(mfrow=c(1, 1), mar = c(4,4,2,2))

util$plot_hist_quantiles(samples, 'seed_count_pred', 0, 40, 1,
                         baseline_values=data$seed_counts, 
                         xlab="Seed counts")

par(mfrow=c(3, 3))

for (t in 1:data$N_traps) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  names <- sapply(idxs,
                  function(n) paste0('seed_count_pred[', n, ']'))
  util$plot_disc_pushforward_quantiles(samples, names, 
                                       baseline_values=data$seed_counts[idxs],
                                       xlab="Year",
                                       xticklab=data$years[idxs],
                                       ylab="Seed Counts",
                                       display_ylim=c(0, 400),
                                       main=paste("Tree", uniq_tree_ids[t]))
}

par(mfrow=c(3, 3))

for (n in c(10, 20, 30, 40, 50, 60, 70, 80, 90)) {
  tid <- tail(which(data$tree_start_idxs <= n), 1)
  name <- paste0('seed_count_pred[', n, ']')
  util$plot_expectand_pushforward(samples[[name]],
                                  50, flim=c(0, 300),
                                  baseline=data$seed_count[n],
                                  display_name="Predicted Seed Counts",
                                  main=paste0("Tree ", uniq_tree_ids[tid],
                                              ", Year ", data$years[n]))
}

# Posterior inferences
par(mfrow=c(2, 4))

util$plot_expectand_pushforward(samples[['lambda1']], 25,
                                display_name="lambda1")
xs <- seq(0, 100, 0.01)
ys <- dnorm(xs, 0, 10 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['lambda2']], 25,
                                display_name="lambda2")
xs <- seq(000, 500, 0.5)
ys <- dnorm(xs, 100, 500 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['psi1']], 25,
                                display_name="psi1")
xs <- seq(0, 5, 0.1)
ys <- dnorm(xs, 0, 2 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['psi2']], 25,
                                display_name="psi2")
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


par(mfrow=c(3, 3))

for (t in 1:data$N_trees) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  names <- sapply(idxs,
                  function(n) paste0('p_masting[', n, ']'))
  util$plot_disc_pushforward_quantiles(samples, names,
                                       xlab="Year",
                                       xticklab=data$years[idxs],
                                       ylab="Posterior Masting Probability",
                                       display_ylim=c(-0.01, 1),
                                       main=paste("Tree", uniq_tree_ids[t]))
}
