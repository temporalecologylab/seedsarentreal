

library(latex2exp)
library(ggplot2)

set.seed(06111954)
wd <- "~/projects/seedsarentreal/analyses/methods"

# Simulate data


ids <- 20

tau_nm_m <- runif(1) # non-masting to masting probability
tau_m_nm <- runif(1, min = 0.8, max = 1) # masting to non-masting probability

tau_nm_m <- 0.5 # non-masting to masting probability
tau_m_nm <- 0.9 # masting to non-masting probability

# rho0 <- runif(1) # initial masting probability
rho0 <- 0
state0 <- ifelse(runif(1) < rho0, 1, 0) # initial state

params <- data.frame()
for(id in 1:ids){
  
  # distribution parameters
  lambda1 <- abs(rnorm(1, 0, 10 / 2.57))
  psi1 <- abs(rnorm(1, 0, 5 / 2.57))
  lambda2 <- abs(rnorm(1, 50, 200 / 2.57))
  while(lambda2 < lambda1){
    lambda2 <- abs(rnorm(1, 50, 200 / 2.57))
  }
  params <- rbind(params, data.frame(lambda1, psi1, lambda2))
}

seeds <- data.frame()
years <- 70
state <- state0
for(y in 1:years){
  
  for(id in 1:ids){
    
    if(state == 1){
      # masting year 
      nseeds <- rpois(1, params[id, 'lambda2'])
    }else if (state == 0){
      # non-masting year 
      nseeds <- MASS::rnegbin(1, params[id, 'lambda1'], 1/params[id, 'psi1'])
    }
    
    seeds <- rbind(seeds, data.frame(id = id, year = y, n = nseeds, state = state))
    
  }
  
  if(state == 1){
    state <- ifelse(runif(1) < tau_m_nm, 0, 1)
  }else if (state == 0){
    state <- ifelse(runif(1) < tau_nm_m, 1, 0)
  }
  
  
}

ggplot(data = seeds) + 
  geom_line(aes(x = year, y = n, group = id),
            linewidth = 0.1, alpha = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank())
