

pdf(file = file.path(figpath, 'supp', 'coldvswarmsummer.pdf'),
    width = 9, height = 2.5)

layout(matrix(c(1,2,3), nrow = 1),
       widths = c(0.5,1, 0.5))
par(mar = c(4,4,1,1))
plot.new()
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