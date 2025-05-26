set.seed(04111954)

# Simulate data
seeds <- data.frame()
params <- data.frame()

ids <- 20

for(id in 1:ids){
  
  # years <- round(runif(1, 40,80),0) # no of obs.
  years <- 70
  
  rho0 <- runif(1) # initial masting probability
  state0 <- ifelse(runif(1) < rho0, 1, 0) # initial state
  tau_nm_m <- runif(1) # non-masting to masting probability
  tau_m_nm <- runif(1) # masting to non-masting probability
  
  # distribution parameters
  lambda1 <- abs(rnorm(1, 0, 10 / 2.57))
  psi1 <- abs(rnorm(1, 0, 5 / 2.57))
  lambda2 <- abs(rnorm(1, 50, 200 / 2.57))
  while(lambda2 < lambda1){
    lambda2 <- abs(rnorm(1, 50, 200 / 2.57))
  }
  psi2 <- abs(rnorm(1, 0, 5 / 2.57))
  
  state <- state0
  seeds_tree <- data.frame()
  for(y in 1:years){
    
    if(state == 1){
      # masting year 
      nseeds <- MASS::rnegbin(1, lambda2, 1/psi2)
      statenext <- ifelse(runif(1) < tau_m_nm, 0, 1)
    }else if (state == 0){
      # non-masting year 
      nseeds <- MASS::rnegbin(1, lambda1, 1/psi1)
      statenext <- ifelse(runif(1) < tau_nm_m, 1, 0)
    }
    
    seeds_tree <- rbind(seeds_tree, data.frame(id = id, year = y, n = nseeds, state = state))
    state <- statenext
  }
  
  seeds <- rbind(seeds, seeds_tree)
  params <- rbind(params, data.frame(rho0, lambda1, psi1, lambda2, psi2, tau_m_nm, tau_nm_m))
  
}

ggplot(data = seeds) + 
  geom_line(aes(x = year, y = n, group = id),
            linewidth = 0.1, alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank())

# sizes
uniq_tree_ids <- unique(seeds$id)
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
  seeds_tree <- seeds[seeds$id == tid,]
  
  years_tree <- seeds_tree$year
  N_years_tree <- length(years_tree)
  
  size_tree <- 1
  seed_count_tree <- seeds_tree$n
  
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

# Posterior quantification
fit <- stan(file= file.path(wd, "stan", "model2_treelevel_alt_diffparam.stan"),
            data=data, seed=04111954, cores =4,
            warmup=1000, iter=4000, refresh=100)

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('rho0',
                                         'tau_nm_m', 'tau_m_nm',
                                         'lambda1', 'psi1',
                                         'lambda2', 'psi2'),
                                       check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

# Compare simulated and predicted parameters
lambda1 <- data.frame()
for(t in 1:N_trees){
  tsamp <- c()
  for(c in 1:4){
    csamp <- samples[[paste0('lambda1[',t,']')]][c,]
    tsamp <- c(tsamp, csamp)
  }
  tsamp <- data.frame(t = t, sim = params[t, 'lambda1'], samp = tsamp)
  lambda1 <- rbind(lambda1, tsamp)
}
psi1 <- data.frame()
for(t in 1:N_trees){
  tsamp <- c()
  for(c in 1:4){
    csamp <- samples[[paste0('psi1[',t,']')]][c,]
    tsamp <- c(tsamp, csamp)
  }
  tsamp <- data.frame(t = t, sim = params[t, 'psi1'], samp = tsamp)
  psi1 <- rbind(psi1, tsamp)
}
lambda2 <- data.frame()
for(t in 1:N_trees){
  tsamp <- c()
  for(c in 1:4){
    csamp <- samples[[paste0('lambda2[',t,']')]][c,]
    tsamp <- c(tsamp, csamp)
  }
  tsamp <- data.frame(t = t, sim = params[t, 'lambda2'], samp = tsamp)
  lambda2 <- rbind(lambda2, tsamp)
}
psi2 <- data.frame()
for(t in 1:N_trees){
  tsamp <- c()
  for(c in 1:4){
    csamp <- samples[[paste0('psi1[',t,']')]][c,]
    tsamp <- c(tsamp, csamp)
  }
  tsamp <- data.frame(t = t, sim = params[t, 'psi2'], samp = tsamp)
  psi2 <- rbind(psi2, tsamp)
}
rho0 <- data.frame()
for(t in 1:N_trees){
  tsamp <- c()
  for(c in 1:4){
    csamp <- samples[[paste0('rho0[',t,']')]][c,]
    tsamp <- c(tsamp, csamp)
  }
  tsamp <- data.frame(t = t, sim = params[t, 'rho0'], samp = tsamp)
  rho0 <- rbind(rho0, tsamp)
}
tau_m_nm <- data.frame()
for(t in 1:N_trees){
  tsamp <- c()
  for(c in 1:4){
    csamp <- samples[[paste0('tau_m_nm[',t,']')]][c,]
    tsamp <- c(tsamp, csamp)
  }
  tsamp <- data.frame(t = t, sim = params[t, 'tau_m_nm'], samp = tsamp)
  tau_m_nm <- rbind(tau_m_nm, tsamp)
}

ggplot() +
  geom_boxplot(data = tau_m_nm, aes(y = samp, x = t, group = t), outliers = FALSE, width = 0.5,
               fill = 'grey90', color = 'grey80') +
  geom_point(data = unique(tau_m_nm[,c('sim', 't')]), aes(y = sim, x = t, group = t), color = 'white', shape = '-',
             size = 10) +
  geom_point(data = unique(tau_m_nm[,c('sim', 't')]), aes(y = sim, x = t, group = t), color = 'red', shape = '-',
             size = 6) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        plot.margin = margin(r = 5, l = 5),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = TeX(r"( $\lambda_{non-masting}$ )"))




ggplot(data = seeds[seeds$id == 2,]) + 
  geom_line(aes(x = year, y = n, group = id),
            linewidth = 0.1, alpha = 0.7) +
  geom_point(aes(x = year, y = n, group = id, color = as.character(state))) +
  scale_color_manual(name = "", values = c("#6E8537FF",'#C1447EFF'),
                     labels = c('Non-masting', 'Masting')) + 
  theme_bw() +
  theme(panel.grid = element_blank())

