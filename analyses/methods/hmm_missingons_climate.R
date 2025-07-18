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



# Data processing
raw_data <- read.csv(file.path(wd, 'data',  'ebms', 'Beech_tree-ring_masting_data.csv'))
clim_data <- readRDS(file.path(wd, 'data',  'ebms', 'era5land_sitesextract.rds'))

length(unique(paste0(raw_data$site.ID,raw_data$tree.ID))) # 57 unique trees
length(unique(raw_data$site.ID)) # 7 unique sites
raw_data$uniqueID <- paste0(raw_data$site.ID, "_", raw_data$tree.ID)
raw_data <- na.omit(raw_data[,c('site.ID', 'uniqueID', 'year', 'seeds')])

# Sizes
uniq_tree_ids <- unique(raw_data$uniqueID)
N_trees <- length(uniq_tree_ids)
max_years <- min(raw_data$year):max(raw_data$year)
N_max_years <- length(max_years)

# Format data into ragged arrays
seed_counts <- c()
springtemp <- c()
prevsummertemp <-c()
observed_years <- c()
years <- c()
N_years <- c()

idx <- 1
tree_start_idxs <- c()
tree_end_idxs <- c()

for (tid in uniq_tree_ids) {
  
  raw_data_tree <- raw_data[raw_data$uniqueID == tid,]
  
  years_tree <- raw_data_tree$year
  years <- c(years, years_tree)
  
  N_years_tree <- length(years_tree)
  N_years <- c(N_years, N_years_tree)
  
  observed_years_tree <- match(years_tree, max_years)
  observed_years <- c(observed_years, observed_years_tree)
  
  seed_count_tree <- raw_data_tree$seeds
  seed_counts <- c(seed_counts, seed_count_tree)
  
  springtemp_tree <- clim_data[clim_data$site.ID == raw_data_tree$site.ID[1] & clim_data$year %in% (max_years), 'meantmean_mam'] # all years, even those unobserved
  springtemp <- c(springtemp, springtemp_tree)
  
  prevsummertemp_tree <- clim_data[clim_data$site.ID == raw_data_tree$site.ID[1] & clim_data$year %in% (max_years-1), 'maxtmax_ja'] # all years, even those unobserved
  prevsummertemp <- c(prevsummertemp, prevsummertemp_tree)
  
  tree_start_idxs <- c(tree_start_idxs, idx)
  idx <- idx + N_years_tree
  tree_end_idxs <- c(tree_end_idxs, idx - 1)
  
}

# Format data
N <- length(years)
data <- mget(c('N', 'N_trees', 'N_max_years',
               'seed_counts', 'observed_years', 
               'springtemp', 'prevsummertemp',
               'years', 'N_years',
               'tree_start_idxs', 'tree_end_idxs'))

# Posterior quantification
fit <- stan(file= file.path(wd, "stan", "model2_treelevel_zinb_missing_summertemp.stan"),
            data=data, seed=5838299, chain = 4, cores = 4,
            warmup=1000, iter=2024, refresh=100)

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('lambda1', 'theta1', 'psi1',
                                         'lambda2', 'psi2', 
                                         'beta_spring', 'beta_prsumm',
                                         'rho0',
                                         'tau_nm_m', 'tau_m_nm'))
util$check_all_expectand_diagnostics(base_samples)



util$plot_pairs_by_chain(samples[['lambda1']], 'lambda1', samples[['psi1']], 'psi1')
util$plot_pairs_by_chain(samples[['lambda2']], 'lambda2', samples[['lambda1']], 'lambda1')
util$plot_pairs_by_chain(samples[['lambda2']], 'lambda2', samples[['psi2']], 'psi2')
util$plot_pairs_by_chain(samples[['lambda2']], 'lambda2', samples[['beta_prsumm']], 'beta_prsumm')
util$plot_pairs_by_chain(samples[['lambda2']], 'lambda2', samples[['beta_spring']], 'beta_spring')
util$plot_pairs_by_chain(samples[['beta_prsumm']], 'beta_prsumm', samples[['beta_spring']], 'beta_spring')

util$plot_pairs_by_chain(samples[['tau_nm_m']], 'tau_nm_m', samples[['tau_m_nm']], 'tau_m_nm')



# Posterior inferences
par(mfrow=c(2, 5))

util$plot_expectand_pushforward(samples[['lambda1']], 25,
                                display_name="lambda1")

util$plot_expectand_pushforward(samples[['psi1']], 25,
                                display_name="psi1")

util$plot_expectand_pushforward(samples[['theta1']], 25,
                                display_name="theta1")

util$plot_expectand_pushforward(samples[['lambda2']], 25,
                                display_name="lambda2")

util$plot_expectand_pushforward(samples[['beta_spring']], 25,
                                display_name="beta_spring")

util$plot_expectand_pushforward(samples[['beta_prsumm']], 25,
                                display_name="beta_prsumm")

util$plot_expectand_pushforward(samples[['psi2']], 25,
                                display_name="psi2")

util$plot_expectand_pushforward(samples[['rho0']], 
                                22, flim=c(0, 1.1),
                                display_name="rho0",)


util$plot_expectand_pushforward(samples[['tau_nm_m']], 
                                25, flim=c(0, 1),
                                display_name="tau_nm_m",)

util$plot_expectand_pushforward(samples[['tau_m_nm']], 
                                25, flim=c(0, 1),
                                display_name="tau_m_nm",)



