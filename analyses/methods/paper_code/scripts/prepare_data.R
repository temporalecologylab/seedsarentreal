
# Load data
raw_data <- read.csv(file.path(wd, 'data',  'ebms', 'Beech_tree-ring_masting_data.csv'))
clim_data <- readRDS(file.path(wd, 'data',  'ebms', 'era5land_sitesextract.rds'))

raw_data$uniqueID <- paste0(raw_data$site.ID, "_", raw_data$tree.ID)
raw_data <- na.omit(raw_data[,c('site.ID', 'uniqueID', 'year', 'seeds')])

uniq_tree_ids <- unique(raw_data$uniqueID)
N_trees <- length(uniq_tree_ids)
N_sites <- length(unique(raw_data$site.ID))
first_year <- min(raw_data$year)
last_year <- max(raw_data$year)
max_years <- first_year:max(raw_data$year)
N_max_years <- length(max_years)

# Format data
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

# Create hypothetical future climates
years_to_predict <- (1980:2100)
## summer temperature
clim_df <- clim_data[clim_data$site.ID == 'Benwell',]
baseline_2000 <- mean(clim_df[clim_df$year %in% c(1980:2000), 'meantmax_ja']) # average summer temp. from 1980-2000 in Benwell
clim_df$year <- clim_df$year - 1978
trend <- 0.05 # 5degC warming in 100 years
summertemp <- baseline_2000 + trend * (years_to_predict-2000)
## frost GDD
simplelm <- lm(gdd_b5_tolastfrost ~ year, data = clim_df)
intercept <- summary(simplelm)$coefficients['(Intercept)','Estimate']
trend <- summary(simplelm)$coefficients['year','Estimate']
baseline_2000 <- intercept + (2000-1978)*trend
frostgdd <- baseline_2000 + trend * (years_to_predict-2000)
frostgdd <- rnorm(length(years_to_predict), mean = predict(simplelm, newdata = data.frame(year = years_to_predict-1978)), sd = 0)
## spring temperature
simplelm <- lm(meantmean_am ~ year, data = clim_df)
intercept <- summary(simplelm)$coefficients['(Intercept)','Estimate']
trend <- summary(simplelm)$coefficients['year','Estimate']
baseline_2000 <- intercept + (2000-1978)*trend
springtemp <- baseline_2000 + trend * (years_to_predict-2000)
springtemp <- rnorm(length(years_to_predict), mean = predict(simplelm, newdata = data.frame(year = years_to_predict-1978)), sd = 0)

# Create hypotehtical trees to predict
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

# Collect data
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