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

# Load data
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

# Format observed data into ragged arrays
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

# Create some data to make some predictions!
years_to_predict <- (1950:2100)
## summer temperature
clim_df <- clim_data[clim_data$site.ID == 'Benwell',]
baseline_2000 <- mean(clim_df[clim_df$year %in% c(1980:2000), 'meantmax_ja']) # average summer temp. from 1980-2000 in Benwell
clim_df$year <- clim_df$year - 1978
trend <- 0.05 # 5degC warming in 100 years
summertemp <- baseline_2000 + trend * (years_to_predict-2000)
ggplot() +
  geom_point(data = clim_df, aes(x = year +1978, y = meantmax_ja)) +
  geom_line(data = data.frame(years_to_predict, summertemp), aes(x = years_to_predict, y = summertemp), color = 'red') +
  theme_minimal()
## frost risk
simplelm <- lm(gdd_b5_tolastfrost ~ year, data = clim_df)
intercept <- summary(simplelm)$coefficients['(Intercept)','Estimate']
trend <- summary(simplelm)$coefficients['year','Estimate']
baseline_2000 <- intercept + (2000-1978)*trend
frostgdd <- baseline_2000 + trend * (years_to_predict-2000)
frostgdd <- rnorm(length(years_to_predict), mean = predict(simplelm, newdata = data.frame(year = years_to_predict-1978)), sd = 0)
ggplot() +
  geom_point(data = clim_df, aes(x = year +1978, y = gdd_b5_tolastfrost)) +
  geom_line(data = data.frame(years_to_predict, frostgdd), aes(x = years_to_predict, y = frostgdd), color = 'red') +
  theme_minimal()
## spring temp
simplelm <- lm(meantmean_am ~ year, data = clim_df)
intercept <- summary(simplelm)$coefficients['(Intercept)','Estimate']
trend <- summary(simplelm)$coefficients['year','Estimate']
baseline_2000 <- intercept + (2000-1978)*trend
springtemp <- baseline_2000 + trend * (years_to_predict-2000)
springtemp <- rnorm(length(years_to_predict), mean = predict(simplelm, newdata = data.frame(year = years_to_predict-1978)), sd = 0)
ggplot() +
  geom_point(data = clim_df, aes(x = year +1978, y = meantmean_am)) +
  geom_line(data = data.frame(years_to_predict, springtemp), aes(x = years_to_predict, y = springtemp), color = 'red') +
  theme_minimal()

## format new data
trees_per_stand <- 100
unique_stands <- "Benwell"
newtree_stand_idxs <- rep(which(unique_stands==unique_stands), each = trees_per_stand)
N_newtrees <- length(newtree_stand_idxs)
N_max_newyears <- length(years_to_predict)
first_newyear <- min(years_to_predict)
newyears <- c()
newprevsummer_temps <-c()
newspring_temps <-c()
newgdd_lastfrost <-c()
N_newyears <- c()
idx <- 1
newtree_start_idxs <- c()
newtree_end_idxs <- c()
for (i in 1:N_newtrees) {
  
  newyears_tree <- years_to_predict-first_newyear+1
  newyears <- c(newyears, newyears_tree)
  
  N_newyears_tree <- length(newyears_tree)
  N_newyears <- c(N_newyears, N_newyears_tree)
  
  newprevsummer_temps_tree <- summertemp
  newprevsummer_temps <- c(newprevsummer_temps, newprevsummer_temps_tree)
  
  newspring_temps_tree <- springtemp
  newspring_temps <- c(newspring_temps, newspring_temps_tree)

  newgdd_lastfrost_tree <- frostgdd/10
  newgdd_lastfrost <- c(newgdd_lastfrost, newgdd_lastfrost_tree)
  
  newtree_start_idxs <- c(newtree_start_idxs, idx)
  idx <- idx + N_newyears_tree
  newtree_end_idxs <- c(newtree_end_idxs, idx - 1)
}

N <- length(years)
Nnew <- length(newyears)
data <- mget(c('N', 'N_trees', 'N_max_years',
               'N_years', 'tree_start_idxs', 'tree_end_idxs',
               'seed_counts', 'years', 
               'prevsummer_temps', 'spring_temps', 'gdd_lastfrost',
               # for predictions
               'Nnew', 'N_newtrees', 'N_max_newyears', 'newtree_stand_idxs',
               'N_newyears', 'newtree_start_idxs', 'newtree_end_idxs','newyears', 
               'newprevsummer_temps', 'newgdd_lastfrost', 'newspring_temps'
))

# Posterior quantification
if(FALSE){
  fit <- stan(file= file.path(wd, "stan", "paper_models", "model2_treelevel_zinb_missing_rescaled_summertempfull_springfrost_springtemp_constrainedslopes_wpred_breakdownpred_standlevel.stan"),
              data=data, seed=5838299, chain = 4, cores = 4,
              warmup=1000, iter=2024, refresh=100)
  saveRDS(fit, file = file.path(wd, 'output', 'fit_25sept2025_fullmodel_constrainedslopes.rds'))
  diagnostics <- util$extract_hmc_diagnostics(fit)
  util$check_all_hmc_diagnostics(diagnostics)
  samples <- util$extract_expectand_vals(fit)
  base_samples <- util$filter_expectands(samples,
                                         c('lambda1', 'theta1', 'psi1',
                                           'lambda20', 'psi2', 
                                           'beta_lambda2_frost', 'beta_lambda2_spring',  
                                           'rho0', 'beta_nm_m', 'beta_m_m',
                                           'tau_nm_m0', 'tau_m_m0'))
  util$check_all_expectand_diagnostics(base_samples)
  
  fit <- stan(file= file.path(wd, "stan", "paper_models", "model2_treelevel_zinb_missing_rescaled_summertemp_springfrost_springtemp_constrainedslopes_wpred.stan"),
              data=data, seed=5838299, chain = 4, cores = 4,
              warmup=1000, iter=2024, refresh=100)
  saveRDS(fit, file = file.path(wd, 'output', 'fit_25sept2025_biologymodel_constrainedslopes.rds'))
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
  
  fit <- stan(file= file.path(wd, "stan", "paper_models", "model2_treelevel_zinb_missing_rescaled_summertemp_springfrost_springtemp_wpred_standlevel.stan"),
              data=data, seed=5838299, chain = 4, cores = 4,
              warmup=1000, iter=2024, refresh=100)
  saveRDS(fit, file = file.path(wd, 'output', 'fit_25sept2025_biologymodel.rds'))
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
}else{
  fit_full_const <- readRDS(file.path(wd, 'output', 'fit_25sept2025_fullmodel_constrainedslopes.rds'))
  fit_biol_const <- readRDS(file = file.path(wd, 'output', 'fit_25sept2025_biologymodel_constrainedslopes.rds'))
  fit_biol <- readRDS(file = file.path(wd, 'output', 'fit_25sept2025_biologymodel.rds'))
}






