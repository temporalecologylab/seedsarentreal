
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
source(file.path(wd, 'paper_code/scripts', 'figure_modeloverview.R'))
source(file.path(wd, 'paper_code/scripts', 'figure_climatepredictors.R'))
source(file.path(wd, 'paper_code/scripts', 'figure_climatenobreakdown.R'))
source(file.path(wd, 'paper_code/scripts', 'figure_synchrony.R'))

# 4. Associated quantities (to include in the text)
summertemp0 <- 15
meansummertemp <- mean(clim_data$meantmax_ja[clim_data$year %in% c(1980:2022)])
average_transition <- util$ensemble_mcmc_quantile_est(invlogit(logit(samples[['tau_nm_m0']]) + 
                                                                 samples[['beta_nm_m']] * (meansummertemp-summertemp0)), c(0.025, 0.5, 0.975))
average_persistence <- util$ensemble_mcmc_quantile_est(samples[['tau_m_m0']], c(0.025, 0.5, 0.975))

warmsummertemp <- 23
samples_transition_warm <- invlogit(logit(samples[['tau_nm_m0']]) + samples[['beta_nm_m']] * (warmsummertemp-summertemp0))
coldsummertemp <- 17
samples_transition_cold <- invlogit(logit(samples[['tau_nm_m0']]) + samples[['beta_nm_m']] * (coldsummertemp-summertemp0))
probability_ratio <- util$ensemble_mcmc_quantile_est(samples_transition_warm/samples_transition_cold, c(0.025, 0.5, 0.975))

frostgdd0 <- 15
lowfrostrisk <- 100
samples_lowfrostrisk <- exp(log(samples[['lambda20']]) + samples[['beta_lambda2_frost']] * (frostrisk/10-frostgdd0))
highfrostrisk <- 3*lowfrostrisk
samples_highfrostrisk <- exp(log(samples[['lambda20']]) + samples[['beta_lambda2_frost']] * (highfrostrisk/10-frostgdd0))
frost_decrease <- util$ensemble_mcmc_quantile_est((samples_highfrostrisk-samples_lowfrostrisk)/samples_lowfrostrisk*100, c(0.025, 0.5, 0.975))

n <- data$N_max_newyears
mastingtrees_warmfuture <- util$ensemble_mcmc_quantile_est(samples_brk[[paste0('states_new_plot[1,', n, ']')]]/data$N_newtrees*100, c(0.025, 0.5, 0.975))


key_quantities <- list(
  "meansummertemp" = meansummertemp,
  "average_transition" = average_transition,
  "average_persistence" = average_persistence,
  "warmsummertemp" = warmsummertemp,
  "coldsummertemp" = coldsummertemp,
  "probability_ratio" = probability_ratio,
  "lowfrostrisk" = lowfrostrisk,
  "highfrostrisk" = highfrostrisk,
  "frost_decrease" = frost_decrease,
  "mastingtrees_warmfuture" = mastingtrees_warmfuture
)
saveRDS(key_quantities, file.path(figpath, 'key_quantities.rds'))
