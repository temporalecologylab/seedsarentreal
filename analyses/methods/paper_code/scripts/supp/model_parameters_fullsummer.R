


pdf(file = file.path(figpath, 'supp', 'negbin_parameters_fullsummer.pdf'),
    width = 9, height = 4)

par(mfrow = c(2,4), mar = c(4,4,1,1))
util$plot_expectand_pushforward(samples_brk[['theta1']], 50, 
                                bquote(theta~'(probability of zeros)'), 
                                flim = c(0,1))
# prior <- rbeta(1e6, 2, 4)
# lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(samples_brk[['lambda1']], 50, 
                                bquote(lambda[low]~'(non-masting mean seed count)'), 
                                flim = c(0,100))
# prior <- rnorm(1e6, 0, 100/2.57)
# lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
util$plot_expectand_pushforward(samples_brk[['psi1']], 50, 
                                bquote(psi[low]~'(non-masting dispersion)'), 
                                flim = c(0,3))
# prior <- rnorm(1e6, 3, 2.5 / 2.57)
# lines(density(prior), col = util$c_light_teal, lwd = 2, lty = 2)
plot.new()

util$plot_expectand_pushforward(samples_brk[['lambda20']], 50, 
                                bquote(lambda[high]~'(masting mean seed count)'),
                                flim = c(100,300))
util$plot_expectand_pushforward(samples_brk[['psi2']], 50, 
                                bquote(psi[high]~'(masting dispersion)'), 
                                flim = c(0,1))
util$plot_expectand_pushforward(samples_brk[['beta_lambda2_spring']], 40, 
                                bquote(beta[high]^{spring}~'(effect of spring temp.)'),
                                flim = c(-0.05, 0.07))
util$plot_expectand_pushforward(samples_brk[['beta_lambda2_frost']], 50, 
                                bquote(beta[high]^{frost}~'(effect of spring frost)'),
                                flim = c(-0.05, 0.07))

dev.off()


pdf(file = file.path(figpath, 'supp', 'transition_parameters_fullsummer.pdf'),
    width = 9, height = 2.5)

par(mfrow = c(1,4), mar = c(5,4,1,1))
util$plot_expectand_pushforward(samples_brk[['tau_m_m0']], 50,
                                expression(logit^{-1}*(alpha[2 %->% 2])~'(persistence in high)'),
                                flim = c(0,1))
util$plot_expectand_pushforward(samples_brk[['beta_m_m']], 40,
                                "",
                                flim = c(0,1.5))

util$plot_expectand_pushforward(samples_brk[['tau_nm_m0']], 50,
                                expression(logit^{-1}*(alpha[1 %->% 2])~'\n(transition low to high)'),
                                flim = c(0,1))
util$plot_expectand_pushforward(samples_brk[['beta_nm_m']], 40,
                                "",
                                flim = c(0,1.5))
par(fig = c(0, 1, 0, 1), new = TRUE, xpd = NA)
text(1.45, -1.6,
     expression(beta[1 %->% 2]^{summer}~'(effect of prev. summer temp.,'),
     xpd = NA)
text(1.45, -2.1,
     "transition low to high)",
     xpd = NA)

par(fig = c(0, 1, 0, 1), new = TRUE, xpd = NA)
text(.3, -1.6,
     expression(beta[2 %->% 2]^{summer}~'(effect of prev. summer temp.,'),
     xpd = NA)
text(.3, -2.1,
     "persistence in high)",
     xpd = NA)

dev.off()