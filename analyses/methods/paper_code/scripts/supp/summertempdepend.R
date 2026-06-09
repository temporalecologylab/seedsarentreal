

pdf(file = file.path(figpath, 'supp', 'summertempdepend.pdf'),
    width = 6.5, height = 4)

tau_nm_m_df <- data.frame()
for(temp_prevsummer in seq(12,26,0.1)){
  tau_nm_m <- matrix(nrow = nrow(samples_brk[['tau_nm_m0']]), ncol = ncol(samples_brk[['tau_nm_m0']]))
  for(c in 1:nrow(tau_nm_m)){
    tau_nm_m[c,] <- sapply(1:ncol(tau_nm_m), function(i) invlogit(logit(samples_brk[['tau_nm_m0']][c,i])+
                                                                    samples_brk[['beta_nm_m']][c,i]*(temp_prevsummer-summertemp0)))
  }
  tau_nm_m_df <- rbind(tau_nm_m_df,
                       data.frame(temp_prevsummer, t(util$ensemble_mcmc_quantile_est(tau_nm_m, probs = c(0.05,0.5,0.95)))))
}

tau_m_m_df <- data.frame()
for(temp_prevsummer in seq(12,26,0.1)){
  tau_m_m <- matrix(nrow = nrow(samples_brk[['tau_m_m0']]), ncol = ncol(samples_brk[['tau_m_m0']]))
  for(c in 1:nrow(tau_m_m)){
    tau_m_m[c,] <- sapply(1:ncol(tau_m_m), function(i) invlogit(logit(samples_brk[['tau_m_m0']][c,i])+
                                                                    samples_brk[['beta_m_m']][c,i]*(temp_prevsummer-summertemp0)))
  }
  tau_m_m_df <- rbind(tau_m_m_df,
                       data.frame(temp_prevsummer, t(util$ensemble_mcmc_quantile_est(tau_m_m, probs = c(0.05,0.5,0.95)))))
}

par(mar = c(2.5,4.5,2,2))
plot(1, type="n", main=main,
     xlim=c(12, 26), xlab='', xaxt="n", yaxt="n",
     ylim=c(0,1), ylab='', frame = FALSE)

lines(tau_nm_m_df$X5. ~ tau_nm_m_df$temp_prevsummer, col = "#a4bed5", lty = 2)
lines(tau_nm_m_df$X95. ~ tau_nm_m_df$temp_prevsummer, col = "#a4bed5", lty = 2)
lines(tau_nm_m_df$X50. ~ tau_nm_m_df$temp_prevsummer, col = '#728595', lwd = 2)

lines(tau_m_m_df$X5. ~ tau_m_m_df$temp_prevsummer, col = "#F0A24A", lty = 2)
lines(tau_m_m_df$X95. ~ tau_m_m_df$temp_prevsummer, col = "#F0A24A", lty = 2)
lines(tau_m_m_df$X50. ~ tau_m_m_df$temp_prevsummer, col = '#B85C2E', lwd = 2)

usr <- par("usr")
axis(1, at=seq(12,26,2), labels=seq(12,26,2), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.85, mgp = c(0, 0.3, 0))
title(xlab = 'Previous summer temperature (°C)', line = 1.5, cex.lab = 0.9)

axis(2, at=seq(0,1,0.25), seq(0,1,0.25), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.85, mgp = c(2, 0.7, 0))  
title(ylab = 'Probability of transition', line = 2.5, cex.lab = 0.9)

legend(
  x = 12, y = 1,
  lty = 1, lwd = 1.5, 
  col = c('#728595', '#B85C2E'),
  legend = c('from low to high state',
            'from high to high state\n(persistence)'),
  box.lwd = NA)

dev.off()


