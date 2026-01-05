
logit <- function(p){
  return(log(p/(1-p)))
}

invlogit <- function(r){
  return(1/(1+exp(-r)))
}

pdf(file = file.path(figpath, 'climate_predictors.pdf'),
    width = 10, height = 3.3)
layout( matrix(c(1,2), ncol=2, byrow = T))

temp0 <- 15
tau_nm_m_df <- data.frame()
for(temp_prevsummer in seq(12,26,0.1)){
  tau_nm_m <- matrix(nrow = nrow(samples[['tau_nm_m0']]), ncol = ncol(samples[['tau_nm_m0']]))
  for(c in 1:nrow(tau_nm_m)){
    tau_nm_m[c,] <- sapply(1:ncol(tau_nm_m), function(i) invlogit(logit(samples[['tau_nm_m0']][c,i])+samples[['beta_nm_m']][c,i]*(temp_prevsummer-temp0)))
  }
  tau_nm_m_df <- rbind(tau_nm_m_df,
                       data.frame(temp_prevsummer, t(util$ensemble_mcmc_quantile_est(tau_nm_m, probs = c(0.1,0.5,0.9)))))
}

par(mar = c(4.5,4.5,2,2))
plot(1, type="n", main=main,
     xlim=c(12, 26), xlab='', xaxt="n", yaxt="n",
     ylim=c(0,1), ylab='', frame = FALSE)

usr <- par("usr")
axis(1, at=seq(12,26,2), labels=seq(12,26,2), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.85, mgp = c(2, 0.7, 0))
title(xlab = 'Previous summer temperature (°C)', line = 2.4, cex.lab = 0.85)

axis(2, at=seq(0,1,0.25), seq(0,1,0.25), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.85, mgp = c(2, 0.7, 0))  
title(ylab = 'Probability of transition from\nlow- to high-reproduction state', line = 2.5, cex.lab = 0.85)

lines(tau_nm_m_df$X50. ~ tau_nm_m_df$temp_prevsummer, col = util$c_dark, lwd = 2)
lines(tau_nm_m_df$X10. ~ tau_nm_m_df$temp_prevsummer, col = util$c_mid, lty = 2)
lines(tau_nm_m_df$X90. ~ tau_nm_m_df$temp_prevsummer, col = util$c_mid, lty = 2)

mtext(LETTERS[1], side = 3, line = 0.5, adj = -0.2, cex = 1.2, font = 2, col = 'grey30')

gdd0 <- 15
lambda2_df <- data.frame()
for(gdd_lastfrost in seq(0,35,0.1)){
  lambda2 <- matrix(nrow = nrow(samples[['lambda20']]), ncol = ncol(samples[['lambda20']]))
  for(c in 1:nrow(lambda2)){
    lambda2[c,] <- sapply(1:ncol(lambda2), function(i) exp(log(samples[['lambda20']][c,i]) +samples[['beta_lambda2_frost']][c,i]*(gdd_lastfrost-gdd0)))
  }
  lambda2_df <- rbind(lambda2_df,
                      data.frame(gdd_lastfrost, t(util$ensemble_mcmc_quantile_est(lambda2, probs = c(0.1,0.5,0.9)))))
}


par(mar = c(4.5,4.5,2,2))
plot(1, type="n", main=main,
     xlim=c(0, 350), xlab='', xaxt="n", yaxt="n",
     ylim=c(120,220), ylab='', frame = FALSE)

usr <- par("usr")
axis(1, at=seq(0,350,50), labels=seq(0,350,50), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.85, mgp = c(2, 0.7, 0))
title(xlab = 'GDD until the last frost day (°C)', line = 2.4, cex.lab = 0.85)

axis(2, at=seq(120,220,20), seq(120,220,20), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.85, mgp = c(2, 0.7, 0))  
title(ylab = 'Average seed production in\na high-reproduction state', line = 2.5, cex.lab = 0.85)

lines(lambda2_df$X50. ~ c(lambda2_df$gdd_lastfrost*10), col = util$c_dark, lwd = 2)
lines(lambda2_df$X10. ~ c(lambda2_df$gdd_lastfrost*10), col = util$c_mid, lty = 2)
lines(lambda2_df$X90. ~ c(lambda2_df$gdd_lastfrost*10), col = util$c_mid, lty = 2)

mtext(LETTERS[2], side = 3, line = 0.5, adj = -0.2, cex = 1.2, font = 2, col = 'grey30')

par(fig = c(0,0.5, 0, 0.6), new = T)  
plot(1, type="n", main=main,
     xlim=c(12, 26), xlab='', xaxt="n", yaxt="n",
     ylim=c(0.4,2.4), ylab='', frame = FALSE)
boxplot(data$prevsummer_temps[1:data$N_max_years], horizontal = TRUE, add = T, frame = F, axes = F)

par(fig = c(0.5,1, 0, 0.6), new = T)  
plot(1, type="n", main=main,
     xlim=c(0, 350), xlab='', xaxt="n", yaxt="n",
     ylim=c(0.4,2.4), ylab='', frame = FALSE)
boxplot(data$gdd_lastfrost[1:data$N_max_years]*10, horizontal = TRUE, add = T, frame = F, axes = F)

dev.off()
