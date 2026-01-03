

pdf(file = file.path(figpath, 'climate_nobreakdown.pdf'),
    width = 11, height = 4)
layout( matrix(c(1,2), ncol=2, byrow = T))

plot.new()

yearstoplot <- 1960:2100-1949
s <- 1

samples_full_const <- util$extract_expectand_vals(fit_full_const)
names <- sapply(yearstoplot, function(n) paste0('states_new_plot[',s, ',', n,']'))
predic_full <- data.frame(
  q10 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.1))}),
  q50 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.5))}),
  q90 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.9))}),
  prevsummertemp = data$newprevsummer_temps[yearstoplot]
)

names <- sapply(yearstoplot, function(n) paste0('states_new_brkd_plot[',s, ',', n,']'))
predic_brkd <- data.frame(
  q10 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.1))}),
  q50 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.5))}),
  q90 = sapply(names, function(n){ util$ensemble_mcmc_quantile_est(samples_full_const[[n]], c(0.9))}),
  prevsummertemp = data$newprevsummer_temps[yearstoplot]
)

ggplot() +
  geom_line(data = predic_full, aes(x = prevsummertemp, y = q10), linetype = 'dashed', col = 'darkgreen') +
  geom_line(data = predic_full, aes(x = prevsummertemp, y = q50), col = 'darkgreen') +
  geom_line(data = predic_full, aes(x = prevsummertemp, y = q90), linetype = 'dashed', col = 'darkgreen') +
  geom_line(data = predic_brkd, aes(x = prevsummertemp, y = q10), linetype = 'dashed', col = 'darkred') +
  geom_line(data = predic_brkd, aes(x = prevsummertemp, y = q50), col = 'darkred') +
  geom_line(data = predic_brkd, aes(x = prevsummertemp, y = q90), linetype = 'dashed', col = 'darkred') +
  theme_classic() +
  labs(y = 'Seed production at the stand-level (100 trees)')


plot(1, type="n", main=main,
     xlim=c(10, 35), xlab='', xaxt="n", yaxt="n",
     ylim=c(0,100), ylab='', frame = FALSE)


usr <- par("usr")
axis(1, at=seq(0,35,5), labels=seq(0,35,5), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))
title(xlab = 'GDD until the last frost day (°C)', line = 2.6)

axis(2, at=seq(0,100,20), seq(0,100,20), 
     lwd = 1, lwd.ticks = 1, tck = -0.03, 
     cex.axis = 0.9, mgp = c(2, 0.7, 0))  
title(ylab = 'Seed production in\nhigh-reproduction state', line = 2.5)

lines(predic_full$q50 ~ c(predic_full$prevsummertemp), col = util$c_mid, lwd = 2)
lines(predic_full$q10 ~ c(predic_full$prevsummertemp), col = util$c_mid, lty = 2)
lines(predic_full$q90 ~ c(predic_full$prevsummertemp), col = util$c_mid, lty = 2)
