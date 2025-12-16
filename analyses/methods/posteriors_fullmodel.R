
# Posterior inferences
par(mfrow=c(2, 7))

util$plot_expectand_pushforward(samples[['lambda1']], 25,
                                flim = c(0,80), 
                                display_name="lambda1")
xs <- seq(0,80,1)
ys <- dnorm(xs, 0, 20 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['psi1']], 25,
                                flim = c(0,2), 
                                display_name="psi1")
xs <- seq(0,2,0.05)
ys <- dnorm(xs, 0, 5 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['theta1']], 25,
                                flim = c(0,1),
                                display_name="theta1")
xs <- seq(0,1,0.05)
ys <- dbeta(xs, 2, 2)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['lambda20']], 25,
                                flim = c(0,300), 
                                display_name="lambda20")
xs <- seq(0,300,1)
ys <- dnorm(xs, 0, 500 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)


util$plot_expectand_pushforward(samples[['psi2']], 25,
                                flim = c(0,2), 
                                display_name="psi2")
xs <- seq(0,2,0.05)
ys <- dnorm(xs, 0, 5 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_lambda2_frost']], 25,
                                flim = c(-0.10,0), 
                                display_name="beta_lambda2_frost")
xs <- seq(-0.1,0.1,0.001)
ys <- dnorm(xs, 0, log(1.3)/2.57)
lines(xs, ys, lwd=2, col=util$c_light)


util$plot_expectand_pushforward(samples[['beta_lambda2_spring']], 25,
                                flim = c(0,0.1), 
                                display_name="beta_lambda2_spring")
xs <- seq(-0.1,0.1,0.001)
ys <- dnorm(xs, 0, log(1.3)/2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['rho0']], 25,
                                flim = c(0,1),
                                display_name="rho0",)
xs <- seq(0,1,0.05)
ys <- dunif(xs, 0, 1)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['tau_nm_m0']], 25, 
                                flim=c(0, 1),
                                display_name="tau_nm_m0",)
xs <- seq(0,1,0.05)
ys <- dbeta(xs, 1, 1.5)
lines(xs, ys, lwd=2, col=util$c_light)


util$plot_expectand_pushforward(samples[['beta_nm_m']], 25,
                                flim = c(0,0.8),
                                display_name="beta_nm_m")
xs <- seq(0,0.6, 0.001)
ys <- dnorm(xs, 0, 0.6/2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['tau_m_m0']], 25, 
                                flim=c(0, 1),
                                display_name="tau_m_m0")
xs <- seq(0,1,0.05)
ys <- dbeta(xs, 1, 1.5)
lines(xs, ys, lwd=2, col=util$c_light)


util$plot_expectand_pushforward(samples[['beta_m_m']], 25,
                                flim=c(0, 0.8),
                                display_name="beta_m_m")
xs <- seq(0,0.6, 0.001)
ys <- dnorm(xs, 0, 0.6/2.57)
lines(xs, ys, lwd=2, col=util$c_light)

