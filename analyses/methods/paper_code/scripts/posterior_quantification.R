
# Posterior quantification
if(FALSE){
  fit <- stan(file= file.path(wd, "stan", "paper_models", "model2_treelevel_zinb_missing_rescaled_summertempfull_springfrost_springtemp_constrainedslopes_wpred_breakdownpred_standlevel.stan"),
              data=data, seed=5838299, chain = 4, cores = 4,
              warmup=1000, iter=2024, refresh=100)
  saveRDS(fit, file = file.path(wd, 'output', 'fit_03jan2026_fullmodel_constrainedslopes.rds'))
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
  
  # fit <- stan(file= file.path(wd, "stan", "paper_models", "model2_treelevel_zinb_missing_rescaled_summertemp_springfrost_springtemp_constrainedslopes_wpred.stan"),
  #             data=data, seed=5838299, chain = 4, cores = 4,
  #             warmup=1000, iter=2024, refresh=100)
  # saveRDS(fit, file = file.path(wd, 'output', 'fit_25sept2025_biologymodel_constrainedslopes.rds'))
  # diagnostics <- util$extract_hmc_diagnostics(fit)
  # util$check_all_hmc_diagnostics(diagnostics)
  # samples <- util$extract_expectand_vals(fit)
  # base_samples <- util$filter_expectands(samples,
  #                                        c('lambda1', 'theta1', 'psi1',
  #                                          'lambda20', 'psi2', 
  #                                          'beta_lambda2_frost', 'beta_lambda2_spring', 
  #                                          'rho0', 'beta_nm_m', 
  #                                          'tau_nm_m0', 'tau_m_m0'))
  # util$check_all_expectand_diagnostics(base_samples)
  
  fit <- stan(file= file.path(wd, "stan", "paper_models", "model2_treelevel_zinb_missing_rescaled_summertemp_springfrost_springtemp_wpred_standlevel.stan"),
              data=data, seed=5838299, chain = 4, cores = 4,
              warmup=1000, iter=2024, refresh=100)
  saveRDS(fit, file = file.path(wd, 'output', 'fit_03jan2026_biologymodel.rds'))
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
  fit_full_const <- readRDS(file.path(wd, 'output', 'fit_03jan2026_fullmodel_constrainedslopes.rds'))
  # fit_biol_const <- readRDS(file = file.path(wd, 'output', 'fit_25sept2025_biologymodel_constrainedslopes.rds'))
  fit_biol <- readRDS(file = file.path(wd, 'output', 'fit_03jan2026_biologymodel.rds'))
}

samples <- util$extract_expectand_vals(fit_biol)
