
# Pipeline to fit HMM and make figures
## V. Van der Meersch, 2025-2026

rm(list = ls())
wd <- '/home/victor/projects/seedsarentreal/analyses/methods'
figpath <- '/home/victor/projects/seedsarentreal/docs/beechmasting/figures'
library(rstan)

# Load functions created by M. Betancourt
setwd(wd)
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

# 1. Prepare data
source(file.path(wd, 'paper_code/scripts', 'prepare_data.R'))

# 2. Posterior quantification
source(file.path(wd, 'paper_code/scripts', 'posterior_quantification.R'))

# 3. Figures 
source(file.path(wd, 'paper_code/scripts', 'figure_modeloverview_extended.R'))
source(file.path(wd, 'paper_code/scripts', 'figure_synchrony_withoutdeltat.R'))
# source(file.path(wd, 'paper_code/scripts', 'figure_climatepredictors.R'))
source(file.path(wd, 'paper_code/scripts', 'figure_climatepredictors_withdeltat.R'))
source(file.path(wd, 'paper_code/scripts', 'figure_climatenobreakdown.R'))
# source(file.path(wd, 'paper_code/scripts', 'figure_synchrony.R'))
source(file.path(wd, 'paper_code/scripts', 'figure_coldsummer_hingepoint.R'))


# 4. Associated quantities (to include in the text)
summertemp0 <- 15
meansummertemp <- mean(clim_data$meantmax_ja[clim_data$year %in% c(1980:2022)])
samples_transitiontohigh <- invlogit(logit(samples[['tau_nm_m0']]) + samples[['beta_nm_m']] * (meansummertemp-summertemp0))
average_transitiontohigh <- util$ensemble_mcmc_quantile_est(samples_transitiontohigh, c(0.05, 0.5, 0.95))
samples_persistenceinlow <- 1 - samples_transitiontohigh
average_persistenceinlow <- util$ensemble_mcmc_quantile_est(samples_persistenceinlow, c(0.05, 0.5, 0.95))
prob_ratio <- util$ensemble_mcmc_quantile_est(samples_transitiontohigh/samples_persistenceinlow, c(0.05, 0.5, 0.95))
average_transitiontolow <- util$ensemble_mcmc_quantile_est(1-samples[['tau_m_m0']], c(0.05, 0.5, 0.95))
average_persistenceinhigh <- util$ensemble_mcmc_quantile_est(samples[['tau_m_m0']], c(0.05, 0.5, 0.95))


warmsummertemp <- 23
samples_transition_warm <- invlogit(logit(samples[['tau_nm_m0']]) + samples[['beta_nm_m']] * (warmsummertemp-summertemp0))
warmsummer_prob <- util$ensemble_mcmc_quantile_est(samples_transition_warm, c(0.05, 0.5, 0.95))
coldsummertemp <- 17
samples_transition_cold <- invlogit(logit(samples[['tau_nm_m0']]) + samples[['beta_nm_m']] * (coldsummertemp-summertemp0))
coldsummer_prob <- util$ensemble_mcmc_quantile_est(samples_transition_cold, c(0.05, 0.5, 0.95))
warmcold_prob_ratio <- util$ensemble_mcmc_quantile_est(samples_transition_warm/samples_transition_cold, c(0.05, 0.5, 0.95))

frostgdd0 <- 15
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

n <- data$N_max_newyears
mastingtrees_warmfuture <- util$ensemble_mcmc_quantile_est(samples_brk[[paste0('states_new_plot[1,', n, ']')]]/data$N_newtrees*100, c(0.05, 0.5, 0.95))


key_quantities <- list(
  "meansummertemp" = meansummertemp,
  "average_transitiontohigh" = average_transitiontohigh,
  "average_persistenceinlow" = average_persistenceinlow,
  "prob_ratio" = prob_ratio,
  "average_transitiontolow" = average_transitiontolow,
  "average_persistenceinhigh" = average_persistenceinhigh,
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
  "mastingtrees_warmfuture" = mastingtrees_warmfuture
)
saveRDS(key_quantities, file.path(figpath, 'key_quantities.rds'))

