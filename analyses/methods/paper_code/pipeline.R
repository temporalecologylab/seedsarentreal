
# Pipeline to fit HMM and make figures
## V. Van der Meersch, 2025-2026

wd <- '/home/victor/projects/seedsarentreal/analyses/methods'
library(rstan)

# Load functions created by M. Betancourt
setwd(wd)
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

# 1. Prepare data
source(file.path(wd, 'paper_code/scripts', 'prepare_data.R'))

# 2. Prepare data
source(file.path(wd, 'paper_code/scripts', 'posterior_quantification.R'))