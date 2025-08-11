

par(mfrow=c(2, 4), mar = c(4,4,1,1))

lambda1 <- rnorm(1e6, 0,20/2.57)
hist(abs(lambda1), prob = TRUE,
     main='', xlab = 'lambda non-masting')

lambda20 <- rnorm(1e6, 0,500/2.57)
hist(abs(lambda20), prob = TRUE,
     main='', xlab = 'lambda masting')

psi <- rnorm(1e6, 0,5/2.57)
hist(abs(psi), prob = TRUE,
     main='', xlab = 'psi1 or psi2')

beta <- rnorm(1e6, 0, log(1.3)/2.57)
hist(-abs(beta), prob = TRUE,
     main='', xlab = 'beta_frost')
hist(abs(beta), prob = TRUE,
     main='', xlab = 'beta_spring')

theta <- rbeta(1e6, 2,2)
hist(theta, prob = TRUE,
     main='', xlab = 'theta')

tau0 <- rbeta(1e6, 1,4)
hist(tau0, prob = TRUE,
     main='', xlab = 'tau_nm_m0 or tau_m_m0')

beta_summer <- rnorm(1e6, 0,0.6/2.57)
hist(abs(beta_summer), prob = TRUE,
     main='', xlab = 'beta_summer (nm_m or m_m)')

# Implications of our priors for tau
logit <- function(p){
  return(log(p/(1-p)))
}
invlogit <- function(r){
  return(1/(1+exp(-r)))
}

par(mfrow=c(2, 5), mar = c(4,4,1,1))
temp_summer <- seq(10,25,1)
tau0s <- rbeta(10, 1,1.5)
beta_summers <- abs(rnorm(100, 0, 0.6/2.57))
probs <- data.frame()
for(tau0 in tau0s){
  plot.new()
  plot.window(xlim = c(10, 25), ylim = c(0, 1))
  axis(side = 1, at = seq(10,25,5))
  axis(side = 2, at = c(0, 0.5, 1))
  tauchange <- c()
  for(beta_summer in beta_summers){
    tau_nm_m <- sapply(temp_summer, function(x) invlogit(logit(tau0)+beta_summer*(x-15)))
    probs <- rbind(probs, data.frame(tau0, beta_summer, temp_summer, tau_nm_m))
    lines(tau_nm_m ~ temp_summer, col = 'grey80')
    tauchange <- c(tauchange, invlogit(logit(tau0) + beta_summer)/invlogit(logit(tau0)))
  }
  text(x = 10, y = 0.95, labels = paste0('Max. chg. for 1deg:\n',round(max(tauchange-1)*100,1),'%'), adj = 0)
  text(x = 23, y = 0.05, labels = paste0('Tau0: ',round(tau0,1), adj = 1))
}


