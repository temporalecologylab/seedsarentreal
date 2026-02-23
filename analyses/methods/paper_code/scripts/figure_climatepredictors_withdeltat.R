
logit <- function(p){
  return(log(p/(1-p)))
}

invlogit <- function(r){
  return(1/(1+exp(-r)))
}

pdf(file = file.path(figpath, 'climate_predictors_withdeltat.pdf'),
    width = 10, height = 5.3)
layout( matrix(c(1,2, 3, 4), ncol=2, byrow = F))

temp0 <- 15
tau_nm_m_df <- data.frame()
for(temp_prevsummer in seq(12,26,0.1)){
  tau_nm_m <- matrix(nrow = nrow(samples[['tau_nm_m0']]), ncol = ncol(samples[['tau_nm_m0']]))
  for(c in 1:nrow(tau_nm_m)){
    tau_nm_m[c,] <- sapply(1:ncol(tau_nm_m), function(i) invlogit(logit(samples[['tau_nm_m0']][c,i])+samples[['beta_nm_m']][c,i]*(temp_prevsummer-temp0)))
  }
  tau_nm_m_df <- rbind(tau_nm_m_df,
                       data.frame(temp_prevsummer, t(util$ensemble_mcmc_quantile_est(tau_nm_m, probs = c(0.05,0.5,0.95)))))
}

### Part 1: transition

par(mar = c(2.5,4.5,2,2))
plot(1, type="n", main=main,
     xlim=c(12, 26), xlab='', xaxt="n", yaxt="n",
     ylim=c(0,1), ylab='', frame = FALSE)

usr <- par("usr")
axis(1, at=seq(12,26,2), labels=seq(12,26,2), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.85, mgp = c(0, 0.3, 0))
title(xlab = 'Previous summer temperature (°C)', line = 1.5, cex.lab = 0.9)

axis(2, at=seq(0,1,0.25), seq(0,1,0.25), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.85, mgp = c(2, 0.7, 0))  
title(ylab = 'Probability of transition from\nlow to high reproductive state', line = 2.5, cex.lab = 0.9)

lines(tau_nm_m_df$X50. ~ tau_nm_m_df$temp_prevsummer, col = util$c_dark, lwd = 2)
lines(tau_nm_m_df$X5. ~ tau_nm_m_df$temp_prevsummer, col = util$c_mid, lty = 2)
lines(tau_nm_m_df$X95. ~ tau_nm_m_df$temp_prevsummer, col = util$c_mid, lty = 2)

mtext(LETTERS[1], side = 3, line = 0.5, adj = -0.18, cex = 1.2, font = 2, col = 'grey30')

### Part 2: frost

gdd0 <- 15
lambda2_df <- data.frame()
for(gdd in seq(0,35,0.1)){
  lambda2 <- matrix(nrow = nrow(samples[['lambda20']]), ncol = ncol(samples[['lambda20']]))
  for(c in 1:nrow(lambda2)){
    lambda2[c,] <- sapply(1:ncol(lambda2), function(i) exp(log(samples[['lambda20']][c,i]) +samples[['beta_lambda2_frost']][c,i]*(gdd-gdd0)))
  }
  lambda2_df <- rbind(lambda2_df,
                      data.frame(gdd, t(util$ensemble_mcmc_quantile_est(lambda2, probs = c(0.05,0.5,0.95)))))
}

par(mar = c(3.5,4.5,1,2))
plot(1, type="n", main=main,
     xlim=c(0, 350), xlab='', xaxt="n", yaxt="n",
     ylim=c(120,220), ylab='', frame = FALSE)

usr <- par("usr")
axis(1, at=seq(0,350,50), labels=seq(0,350,50), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.85, mgp = c(0, 0.3, 0))
title(xlab = 'GDD until the last frost day (°C)', line = 1.5, cex.lab = 0.9)

axis(2, at=seq(120,220,20), seq(120,220,20), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.85, mgp = c(2, 0.7, 0))  
title(ylab = 'Average seed production in\na high reproductive state', line = 2.5, cex.lab = 0.9)

lines(lambda2_df$X50. ~ c(lambda2_df$gdd*10), col = util$c_dark, lwd = 2)
lines(lambda2_df$X5. ~ c(lambda2_df$gdd*10), col = util$c_mid, lty = 2)
lines(lambda2_df$X95. ~ c(lambda2_df$gdd*10), col = util$c_mid, lty = 2)

mtext(LETTERS[2], side = 3, line = 0.5, adj = -0.18, cex = 1.2, font = 2, col = 'grey30')

# par(fig = c(0,0.5, 0, 0.6), new = T)  
# plot(1, type="n", main=main,
#      xlim=c(12, 26), xlab='', xaxt="n", yaxt="n",
#      ylim=c(0.4,2.4), ylab='', frame = FALSE)
# boxplot(clim_data$meantmax_ja[clim_data$year %in% c(1980:2022)], horizontal = TRUE, add = T, frame = F, axes = F,
#         outline=FALSE)
# 
# par(fig = c(0.5,1, 0, 0.6), new = T)  
# plot(1, type="n", main=main,
#      xlim=c(0, 350), xlab='', xaxt="n", yaxt="n",
#      ylim=c(0.4,2.4), ylab='', frame = FALSE)
# boxplot(clim_data$gdd_b5_tolastfrost[clim_data$year %in% c(1980:2022)], horizontal = TRUE, add = T, frame = F, axes = F,
#         outline=FALSE)


### Part 3: deltaT

names <- paste0('acrossstands_masting_synchrony[', 1:N_max_years, ']')
probs <- c(0.1, 0.5, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[names[n]]], probs)
}
quantiles_m <- sapply(1:N_max_years, calc)

names <- paste0('acrossstands_nonmasting_synchrony[', 1:N_max_years, ']')
probs <- c(0.1, 0.5, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[names[n]]], probs)
}
quantiles_nm <- sapply(1:N_max_years, calc)


prevsummer_n2 <- aggregate(meantmax_ja ~ year, data = clim_data[clim_data$year %in% c(1978:2020),], FUN = mean)$meantmax_ja
prevsummer_n1 <- aggregate(meantmax_ja ~ year, data = clim_data[clim_data$year %in% c(1979:2021),], FUN = mean)$meantmax_ja

delta <- prevsummer_n1-prevsummer_n2

par(mar=c(1.5,5,2,0.5), cex.lab = 0.85)
plot(1, type="n", main=main,
     xlim=c(-4, 4), xlab='', 
     yaxt = 'n', xaxt = 'n',
     ylim=c(-0.05,1.05), ylab='',
     frame = FALSE)
title(ylab = 'Masting synchrony\nbetween populations (% of trees)', line = 2, cex.lab = 0.9)

# axis(2, at= seq(0,1,0.1), labels=seq(0,1,0.1))
axis(2, at=seq(0,1,0.25), seq(0,100,25), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.85, mgp = c(2, 0.7, 0)) 

segments(x0 = delta, y0 = quantiles_m['10%',], y1 = quantiles_m['90%',], col=util$c_light, lwd = 0.75)
points(quantiles_m['50%',] ~ delta, col = util$c_light, pch = 20)
fit <- lm(quantiles_m['50%',] ~ delta)
xlim <- c(-4, 4)
x <- seq(xlim[1], xlim[2], length.out = 100)
lines(x, predict(fit, newdata = data.frame(delta = x)), col = util$c_mid, lty = 2, lwd = 1.2)



mtext(LETTERS[3], side = 3, line = 0.5, adj = -0.18, cex = 1.2, font = 2, col = 'grey30')


par(mar=c(3.5,5,0,0.5), cex.lab = 0.85)
plot(1, type="n", main=main,
     xlim=c(-4, 4), xlab='', 
     yaxt = 'n', xaxt = 'n',
     ylim=c(-0.05,1.05), ylab='',
     frame = FALSE)
title(ylab = 'Non-masting synchrony\nbetween populations (% of trees)', line = 2, cex.lab = 0.9)

axis(2, at=seq(0,1,0.25), seq(0,100,25), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.85, mgp = c(2, 0.7, 0)) 

axis(1, at=seq(-4,4,1), seq(-4,4,1), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.85, mgp = c(0, 0.3, 0)) 
title(xlab = 'Temperature diff. between two previous summers (°C)', line = 1.5, cex.lab = 0.9)

segments(x0 = delta, y0 = quantiles_nm['10%',], y1 = quantiles_nm['90%',], col='#96b096', lwd = 0.75)
points(quantiles_nm['50%',] ~ delta, col = '#96b096', pch = 20)
fit <- lm(quantiles_nm['50%',] ~ delta)
xlim <- c(-4, 4)
x <- seq(xlim[1], xlim[2], length.out = 100)
lines(x, predict(fit, newdata = data.frame(delta = x)), col = "#487548", lty = 2, lwd = 1.2)

dev.off()
