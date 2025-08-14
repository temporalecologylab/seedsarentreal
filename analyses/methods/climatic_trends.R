
library(ggplot2)
library(ggpubr)

ggplot(data = clim_data, ) +
  facet_wrap(~site.ID, ncol = 2) +
  geom_point(aes(x = year-1978, y = meantmax_ja), col = '#ADBE7CFF', size = 0.8) +
  geom_smooth(aes(x = year-1978, y = meantmax_ja), method = 'lm', se = FALSE, col = '#D06C9BFF') +
  stat_regline_equation(aes(x = year-1978, y = meantmax_ja),
                        label.x = -1, label.y = 23.5, hjust = 0, col = '#D06C9BFF') +  
  scale_x_continuous(breaks = seq(2,42,10), labels = 1978+seq(2,42,10)) +
  labs(x = '') + 
  theme_classic()

ggplot(data = clim_data, ) +
  facet_wrap(~site.ID, ncol = 2) +
  geom_point(aes(x = year-1978, y = gdd_b5_tolastfrost), col = '#ADBE7CFF', size = 0.8) +
  geom_smooth(aes(x = year-1978, y = gdd_b5_tolastfrost), method = 'lm', se = FALSE, col = '#D06C9BFF') +
  stat_regline_equation(aes(x = year-1978, y = gdd_b5_tolastfrost),
                        label.x = -1, label.y = 320, hjust = 0, col = '#D06C9BFF') +  
  scale_x_continuous(breaks = seq(2,42,10), labels = 1978+seq(2,42,10)) +
  geom_smooth(aes(x = year-1978, y = gdd_b5_tolastfrost), method = 'loess', se = FALSE, 
              col = '#DA9FB8FF', linewidth = 1, linetype = 'dashed') +
  labs(x = '') + 
  theme_classic()

newclim_data <- data.frame()
for(s in unique(clim_data$site.ID)){
  clim_s <- clim_data[clim_data$site.ID == s,]
  clim_s$year <- clim_s$year - 1978
  # summer temperature
  intercept <- summary(lm(meantmax_ja ~ year, data = clim_s))$coefficients['(Intercept)','Estimate']
  trend <- summary(lm(meantmax_ja ~ year, data = clim_s))$coefficients['year','Estimate']
  
  delta2025 <- trend * (2025-1978)
  newclim_s <- clim_s
  newclim_s$year <- newclim_s$year + (2025-1978)
  newclim_s$meantmax_ja <- newclim_s$meantmax_ja + delta2025
  
  delta2025 <- trend * (2025-1978)
  newclim_s_2 <- newclim_s
  newclim_s_2$year <- newclim_s_2$year + (2025-1978)
  newclim_s_2$meantmax_ja <- newclim_s_2$meantmax_ja + delta2025
  
  
  newclim_data <- rbind(newclim_data, rbind(clim_s, newclim_s, newclim_s_2))
}
newclim_data$year <- newclim_data$year + 1978
ggplot(data = newclim_data, ) +
  facet_wrap(~site.ID, ncol = 2) +
  geom_point(aes(x = year-1978, y = meantmax_ja), col = '#ADBE7CFF', size = 0.4) +
  geom_smooth(aes(x = year-1978, y = meantmax_ja, linetype = year > 2024), method = 'lm', se = FALSE, col = '#D06C9BFF') +
  stat_regline_equation(aes(x = year-1978, y = meantmax_ja),
                        label.x = -1, label.y = 25.5, hjust = 0, col = '#D06C9BFF') +  
  scale_x_continuous(breaks = seq(2,122,30), labels = 1978+seq(2,122,30)) +
  labs(x = '') + 
  theme_classic() +
  theme(legend.position = 'none')

newclim_data <- data.frame()
for(s in unique(clim_data$site.ID)){
  clim_s <- clim_data[clim_data$site.ID == s,]
  clim_s$year <- clim_s$year - 1978
  # summer temperature
  intercept <- summary(lm(meantmax_ja ~ year, data = clim_s))$coefficients['(Intercept)','Estimate']
  trend <- summary(lm(meantmax_ja ~ year, data = clim_s))$coefficients['year','Estimate']
  newtrend <- 0.02
  delta2025 <- trend * (2025-1978)
  newclim_s <- clim_s
  newclim_s$year <- newclim_s$year + (2025-1978)
  newclim_s$meantmax_ja <- newclim_s$meantmax_ja + delta2025 + newtrend * (newclim_s$year-(2025-1978))
  
  delta2025 <- trend * (2025-1978)
  newclim_s_2 <- newclim_s
  newclim_s_2$year <- newclim_s_2$year + (2025-1978)
  newclim_s_2$meantmax_ja <- newclim_s_2$meantmax_ja + delta2025 + newtrend * (newclim_s$year-(2025-1978))
  
  
  newclim_data <- rbind(newclim_data, rbind(clim_s, newclim_s, newclim_s_2))
}
newclim_data$year <- newclim_data$year + 1978
ggplot(data = newclim_data, ) +
  facet_wrap(~site.ID, ncol = 2) +
  geom_point(aes(x = year-1978, y = meantmax_ja), col = '#ADBE7CFF', size = 0.4) +
  geom_smooth(aes(x = year-1978, y = meantmax_ja, linetype = year > 2024), method = 'lm', se = FALSE, col = '#D06C9BFF') +
  stat_regline_equation(aes(x = year-1978, y = meantmax_ja),
                        label.x = -1, label.y = 26.5, hjust = 0, col = '#D06C9BFF') +  
  scale_x_continuous(breaks = seq(2,122,30), labels = 1978+seq(2,122,30)) +
  labs(x = '') + 
  theme_classic() +
  theme(legend.position = 'none')

