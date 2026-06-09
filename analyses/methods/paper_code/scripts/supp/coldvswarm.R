

pdf(file = file.path(figpath, 'supp', 'coldvswarmsummer.pdf'),
    width = 5, height = 3.5)

par(mfrow = c(1,1), mar = c(4,4,1,1))
util$plot_expectand_pushforward(samples_transition_warm, 50, 
                                "Probability of transition from low to high state", 
                                flim = c(0,1), col = util$c_mid)
util$plot_expectand_pushforward(samples_transition_cold, 50, 
                                bquote(lambda[low]~'(non-masting mean seed count)'), 
                                flim = c(0,1), add = T, col = util$c_light_teal)

legend(
  x = 0, y = 25,
  lty = 1, lwd = 1.5, 
  col = c(util$c_light_teal, util$c_mid),
  legend = c(paste0('Cold summer (', coldsummertemp, '°C)'),
             paste0('Warm summer (', warmsummertemp, '°C)')),
  box.lwd = NA)

dev.off()