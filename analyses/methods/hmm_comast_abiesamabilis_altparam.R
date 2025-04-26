
kippenberger <- c("#8B174DFF", "#AE2565FF", "#C1447EFF", "#D06C9BFF", "#DA9FB8FF", 
                  "#ADBE7CFF", "#8BA749FF", "#6E8537FF", "#4F5F28FF", "#343D1FFF")

tau_nm_m <- data.frame()
for(t in 1:N_trees){
  tsamp <- c()
  for(c in 1:4){
    csamp <- samples[[paste0('tau_nm_m[',t,']')]][c,]
    tsamp <- c(tsamp, csamp)
  }
  tsamp <- data.frame(t = t, s = tsamp)
  tau_nm_m <- rbind(tau_nm_m, tsamp)
}

ggplot(data = tau_nm_m) +
  geom_density(aes(x = s), linewidth = 0.3) +
  scale_color_gradientn(colors = kippenberger) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        plot.margin = margin(r = 10, l = 5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
  labs(y = "Density", x = TeX(r"( P(non-masting$\rightarrow$masting) )"))

ggplot(data = tau_nm_m) +
  geom_boxplot(aes(x = s, y = t, group = t, color = t), outliers = FALSE, width = 0.75) +
  scale_color_gradientn(colors = kippenberger) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))

pnmm <- tau_nm_m %>%
  group_by(t) %>%
  summarise(smed = median(s))

tau_m_nm <- data.frame()
for(t in 1:N_trees){
  tsamp <- c()
  for(c in 1:4){
    csamp <- samples[[paste0('tau_m_nm[',t,']')]][c,]
    tsamp <- c(tsamp, csamp)
  }
  tsamp <- data.frame(t = t, s = tsamp)
  tau_m_nm <- rbind(tau_m_nm, tsamp)
}

ggplot(data = tau_m_nm) +
  geom_density(aes(x = s, group = t, color = t), linewidth = 0.3) +
  scale_color_gradientn(colors = kippenberger) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        plot.margin = margin(r = 10, l = 5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
  labs(y = "Density", x = TeX(r"( P(masting$\rightarrow$non-masting) )"))

ggplot(data = tau_m_nm) +
  geom_boxplot(aes(x = s, y = t, group = t, color = t), outliers = FALSE, width = 0.75) +
  scale_color_gradientn(colors = kippenberger) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))



lambda2 <- data.frame()
for(t in 1:N_trees){
  tsamp <- c()
  for(c in 1:4){
    csamp <- samples[[paste0('lambda2[',t,']')]][c,]
    tsamp <- c(tsamp, csamp)
  }
  tsamp <- data.frame(t = t, s = tsamp)
  lambda2 <- rbind(lambda2, tsamp)
}

ggplot(data = lambda2) +
  geom_density(aes(x = s, group = t, color = t), linewidth = 0.3) +
  scale_color_gradientn(colors = kippenberger) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        plot.margin = margin(r = 10, l = 5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(expand = c(0,0), limits = c(0,200)) +
  labs(y = "Density", x = TeX(r"( $\lambda_{masting}$ )"))

lambda1 <- data.frame()
for(t in 1:N_trees){
  tsamp <- c()
  for(c in 1:4){
    csamp <- samples[[paste0('lambda1[',t,']')]][c,]
    tsamp <- c(tsamp, csamp)
  }
  tsamp <- data.frame(t = t, s = tsamp)
  lambda1 <- rbind(lambda1, tsamp)
}

ggplot(data = lambda1) +
  geom_density(aes(x = s, group = t, color = t), linewidth = 0.3) +
  scale_color_gradientn(colors = kippenberger) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        plot.margin = margin(r = 10, l = 5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(expand = c(0,0), limits = c(0,2)) +
  labs(y = "Density", x = TeX(r"( $\lambda_{non-masting}$ )"))



rho0 <- data.frame()
for(t in 1:N_trees){
  tsamp <- c()
  for(c in 1:4){
    csamp <- samples[[paste0('rho0[',t,']')]][c,]
    tsamp <- c(tsamp, csamp)
  }
  tsamp <- data.frame(t = t, s = tsamp)
  rho0 <- rbind(rho0, tsamp)
}

ggplot(data = rho0) +
  geom_density(aes(x = s, group = t, color = t), linewidth = 0.3) +
  scale_color_gradientn(colors = kippenberger) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        plot.margin = margin(r = 10, l = 5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
  labs(y = "Density", x = TeX(r"( $\rho_0 )"))


year <- 2000
pmasting <- data.frame()
for(t in 1:N_trees){
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  n <- idxs[year-1961]
  tsamp <- c()
  for(c in 1:4){
    csamp <- samples[[paste0('p_masting[',n,']')]][c,]
    tsamp <- c(tsamp, csamp)
  }
  tsamp <- data.frame(t = t, s = tsamp)
  pmasting <- rbind(pmasting, tsamp)
}

ggplot(data = left_join(pmasting, pnmm)) +
  geom_vline(xintercept = 0.5, color = 'grey', linetype = 'dashed') +
  geom_boxplot(aes(x = s, y = reorder(t, s, FUN = median), group = t, color = smed), outliers = FALSE, width = 0.65) +
  scale_colour_gradient2(midpoint = 0.5) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.margin = margin(r = 10, l = 5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1.01)) +
  labs(y = "", x = TeX(r"( Posterior $P_{masting}$ in 2016 )"), color = TeX(r"( $\lambda_{masting}$ )"))





ggplot(data = pmasting) +
  facet_wrap(~ t, scales = 'free_y') +
  geom_density(aes(x = s, group = t, color = t)) +
  scale_color_gradientn(colors = kippenberger) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        plot.margin = margin(r = 10, l = 5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1.01)) +
  labs(y = "Density", x = TeX(r"( Posterior $P_{masting}$ in 2000 )"))
