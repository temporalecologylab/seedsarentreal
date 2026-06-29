
# Pipeline to fit HMM and make figures
## V. Van der Meersch, 2025-2026

rm(list = ls())
wd <- '/home/victor/projects/seedsarentreal/analyses/methods'
figpath <- '/home/victor/projects/seedsarentreal/docs/beechmasting/figures/new'
library(rstan)
library(terra)

# Load functions created by M. Betancourt
setwd(wd)
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)


logit <- function(p){
  return(log(p/(1-p)))
}

invlogit <- function(r){
  return(1/(1+exp(-r)))
}

# 1. Prepare data
source(file.path(wd, 'paper_code/scripts/new', 'prepare_data.R'))

# 2. Posterior quantification
source(file.path(wd, 'paper_code/scripts/new', 'posterior_quantification.R'))

# Baselines in the model
summertemp0 <- 17
frostgdd0 <- 5
springtemp0 <- 9

# For the figures, should we keep only trees with at least 20 obs.?
# sum(data$N_years < 20) # we would lost 97 trees
# trees_tokeep <- which(data$N_years >= 20)
trees_tokeep <- 1:data$N_trees # we keep all trees (makes no difference anyway)

# 3. Figures 
source(file.path(wd, 'paper_code/scripts/new', 'figure_modeloverview_extended.R'))
source(file.path(wd, 'paper_code/scripts', 'figure_synchrony_withoutdeltat.R'))
# source(file.path(wd, 'paper_code/scripts', 'figure_climatepredictors.R'))
source(file.path(wd, 'paper_code/scripts', 'figure_climatepredictors_withdeltat.R'))
source(file.path(wd, 'paper_code/scripts', 'figure_climatenobreakdown.R'))
# source(file.path(wd, 'paper_code/scripts', 'figure_synchrony.R'))
source(file.path(wd, 'paper_code/scripts/new', 'figure_coldsummer_hingepoint.R'))
source(file.path(wd, 'paper_code/scripts/new','figure_merged_synchrony_hingepoint.R'))


# 4. Associated quantities (to include in the text)
meansummertemp <- mean(clim_data$meantmax_ja[clim_data$year %in% c(1980:2024)])


# Average seed production predicted in each state
samples[['sumseeds_inmasting']] <- 0
samples[['numtrees_inmasting']] <- 0
samples[['sumseeds_innonmasting']] <- 0
samples[['numtrees_innonmasting']] <- 0
for(t in 1:data$N_trees){
  
  idxs <- seq(1+(t-1)*data$N_max_years, t*data$N_max_years, 1)
  
  for(i in idxs){
    
    samples[['sumseeds_inmasting']] <- samples[['sumseeds_inmasting']] +
      ifelse(samples[[paste0('states_pred[', i,']')]] == 2,samples[[paste0('seed_counts_pred[', i,']')]],0)
    samples[['numtrees_inmasting']] <- samples[['numtrees_inmasting']] +
      ifelse(samples[[paste0('states_pred[', i,']')]] == 2,1,0)
    
    samples[['sumseeds_innonmasting']] <- samples[['sumseeds_innonmasting']] +
      ifelse(samples[[paste0('states_pred[', i,']')]] == 1,samples[[paste0('seed_counts_pred[', i,']')]],0)
    samples[['numtrees_innonmasting']] <- samples[['numtrees_innonmasting']] +
      ifelse(samples[[paste0('states_pred[', i,']')]] == 1,1,0)
    
  }
}
samples[['meanseeds_inmasting']] <- samples[['sumseeds_inmasting']]/samples[['numtrees_inmasting']] 
samples[['meanseeds_innonmasting']] <- samples[['sumseeds_innonmasting']]/samples[['numtrees_innonmasting']] 
meanseeds_inmasting <- util$ensemble_mcmc_quantile_est(samples[['meanseeds_inmasting']], c(0.05, 0.5, 0.95))
meanseeds_innonmasting <- util$ensemble_mcmc_quantile_est(samples[['meanseeds_innonmasting']], c(0.05, 0.5, 0.95))


warmsummertemp <- 22
samples_transition_warm <- invlogit(logit(samples[['tau_nm_m0']]) + samples[['beta_nm_m']] * (warmsummertemp-summertemp0))
warmsummer_prob <- util$ensemble_mcmc_quantile_est(samples_transition_warm, c(0.05, 0.5, 0.95))
coldsummertemp <- 18
samples_transition_cold <- invlogit(logit(samples[['tau_nm_m0']]) + samples[['beta_nm_m']] * (coldsummertemp-summertemp0))
coldsummer_prob <- util$ensemble_mcmc_quantile_est(samples_transition_cold, c(0.05, 0.5, 0.95))
warmcold_prob_ratio <- util$ensemble_mcmc_quantile_est(samples_transition_warm/samples_transition_cold, c(0.05, 0.5, 0.95))

lowfrostrisk <- 100
samples_lowfrostrisk <- exp(log(samples[['lambda20']]) + samples[['beta_lambda2_frost']] * (lowfrostrisk/10-frostgdd0))
highfrostrisk <- 3*lowfrostrisk
samples_highfrostrisk <- exp(log(samples[['lambda20']]) + samples[['beta_lambda2_frost']] * (highfrostrisk/10-frostgdd0))
frost_decrease <- util$ensemble_mcmc_quantile_est((samples_highfrostrisk-samples_lowfrostrisk)/samples_lowfrostrisk*100, c(0.05, 0.5, 0.95))


samples[['withinstand_synchrony_mean']] <- 0
for(n in paste0('withinstand_synchrony[', 1:data$N_max_years, ']')){
  samples[['withinstand_synchrony_mean']] <- samples[['withinstand_synchrony_mean']] +
    samples[[n]]
}
samples[['withinstand_synchrony_mean']] <- samples[['withinstand_synchrony_mean']]/data$N_max_years
withinstand_synchrony <- util$ensemble_mcmc_quantile_est(samples[['withinstand_synchrony_mean']]*100, c(0.05, 0.5, 0.95))


samples[['acrossstands_synchrony_mean']] <- 0
for(n in paste0('acrossstands_synchrony[', 1:data$N_max_years, ']')){
  samples[['acrossstands_synchrony_mean']] <- samples[['acrossstands_synchrony_mean']] +
    samples[[n]]
}
samples[['acrossstands_synchrony_mean']] <- samples[['acrossstands_synchrony_mean']]/data$N_max_years
acrossstands_synchrony <- util$ensemble_mcmc_quantile_est(samples[['acrossstands_synchrony_mean']]*100, c(0.05, 0.5, 0.95))

samples[['acrossstands_synchrony_mean_before2015']] <- 0
for(n in paste0('acrossstands_synchrony[', which(max_years %in% 1980:2014), ']')){
  samples[['acrossstands_synchrony_mean_before2015']] <- samples[['acrossstands_synchrony_mean_before2015']] +
    samples[[n]]
}
samples[['acrossstands_synchrony_mean_before2015']] <- samples[['acrossstands_synchrony_mean_before2015']]/length(1980:2014)
acrossstands_synchrony_before2015 <- util$ensemble_mcmc_quantile_est(samples[['acrossstands_synchrony_mean_before2015']]*100, c(0.05, 0.5, 0.95))

samples[['acrossstands_synchrony_mean_2015onwards']] <- 0
for(n in paste0('acrossstands_synchrony[', which(max_years %in% 2015:2022), ']')){
  samples[['acrossstands_synchrony_mean_2015onwards']] <- samples[['acrossstands_synchrony_mean_2015onwards']] +
    samples[[n]]
}
samples[['acrossstands_synchrony_mean_2015onwards']] <- samples[['acrossstands_synchrony_mean_2015onwards']]/length(2015:2022)
acrossstands_synchrony_2015onwards <- util$ensemble_mcmc_quantile_est(samples[['acrossstands_synchrony_mean_2015onwards']]*100, c(0.05, 0.5, 0.95))


samples[['withinstand_synchrony_mean_before2015']] <- 0
for(n in paste0('withinstand_synchrony[', which(max_years %in% 1980:2014), ']')){
  samples[['withinstand_synchrony_mean_before2015']] <- samples[['withinstand_synchrony_mean_before2015']] +
    samples[[n]]
}
samples[['withinstand_synchrony_mean_before2015']] <- samples[['withinstand_synchrony_mean_before2015']]/length(1980:2014)
withinstand_synchrony_before2015 <- util$ensemble_mcmc_quantile_est(samples[['withinstand_synchrony_mean_before2015']]*100, c(0.05, 0.5, 0.95))

samples[['withinstand_synchrony_mean_2015onwards']] <- 0
for(n in paste0('withinstand_synchrony[', which(max_years %in% 2015:2022), ']')){
  samples[['withinstand_synchrony_mean_2015onwards']] <- samples[['withinstand_synchrony_mean_2015onwards']] +
    samples[[n]]
}
samples[['withinstand_synchrony_mean_2015onwards']] <- samples[['withinstand_synchrony_mean_2015onwards']]/length(2015:2022)
withinstand_synchrony_2015onwards <- util$ensemble_mcmc_quantile_est(samples[['withinstand_synchrony_mean_2015onwards']]*100, c(0.05, 0.5, 0.95))

n <- data$N_max_newyears
mastingtrees_warmfuture <- util$ensemble_mcmc_quantile_est(samples_brk[[paste0('states_new_plot[1,', n, ']')]]/data$N_newtrees*100, c(0.05, 0.5, 0.95))

# Distance between sites
sites_coords <- read.csv(file.path(wd, 'data',  'ebms', 'sites_ebms.csv'))
sites_vect <-  vect(sites_coords, geom=c("Longitude", "Latitude"))
D <- distance(sites_vect, sites_vect, unit="km") 
D[upper.tri(D)] <- NA
D[D==0] <- NA
meandistance <- mean(D, na.rm = T)
mindistance <- min(D, na.rm = T)

# Years observed by sites:
maxyears_persite <- rep(0, length(uniq_sites))
for(s in unique(sites)){
  for(t in which(sites == s)){
    maxyears_persite[which( unique(sites) == s)] <- ifelse(data$N_years[t] > maxyears_persite[which( unique(sites) == s)], 
                                                            data$N_years[t], maxyears_persite[which( unique(sites) == s)])
  }
}

key_quantities <- list(
  "meansummertemp" = meansummertemp,
  "average_transitiontohigh" = average_transitiontohigh,
  "average_persistenceinlow" = average_persistenceinlow,
  "prob_ratio" = prob_ratio,
  "average_transitiontolow" = average_transitiontolow,
  "average_persistenceinhigh" = average_persistenceinhigh,
  "meanseeds_inmasting" = meanseeds_inmasting,
  "meanseeds_innonmasting" = meanseeds_innonmasting,
  "warmsummertemp" = warmsummertemp,
  "warmsummer_prob" = warmsummer_prob,
  "coldsummertemp" = coldsummertemp,
  "coldsummer_prob" = coldsummer_prob,
  "warmcold_prob_ratio" = warmcold_prob_ratio,
  "lowfrostrisk" = lowfrostrisk,
  "highfrostrisk" = highfrostrisk,
  "frost_decrease" = frost_decrease,
  "withinstand_synchrony" = withinstand_synchrony,
  "acrossstands_synchrony" = acrossstands_synchrony,
  "acrossstands_synchrony_before2015" = acrossstands_synchrony_before2015,
  "acrossstands_synchrony_2015onwards" = acrossstands_synchrony_2015onwards,
  "mastingtrees_warmfuture" = mastingtrees_warmfuture,
  "meandistance" = meandistance,
  "mindistance" = mindistance,
  "ntrees" = data$N_trees,
  "nsites" = length(uniq_sites),
  "ntrees_persite" = range(trees_persite),
  "nyears_persite" = range(maxyears_persite),
  "avg_nyears_pertree" = round(mean(data$N_years),1)
  
)
saveRDS(key_quantities, file.path(figpath, 'key_quantities.rds'))

# Supplementary figures







