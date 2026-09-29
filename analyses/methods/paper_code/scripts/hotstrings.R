
prob_vals <- readRDS(file = file.path(wd, 'output', 'meantmaxja_probproj_ukcp18.rds'))

hotsummer_strings <- data.frame()
rcps <- c('rcp45', 'rcp85')
base <- c(2020:2039)
decades <- list(
  base, base+10, base+20, base+30, base+40, base+50, base+60
)

for(rcp in rcps){
  for(i in 1:length(decades)){
    dec <- decades[[i]]
    print(mean(dec))
    for(site in unique(prob_vals$site)){
      dat <- prob_vals[prob_vals$site == site & prob_vals$year %in% dec & prob_vals$rcp == rcp,]
      samples <- unique(dat$sample)
      for(s in samples){
        
        dat_s <- dat[dat$sample == s,]
        means <- dat_s[order(dat_s$year), 'mean']
        # above3 <- rle(means > 3)
        # above5 <- rle(means > 5)
        # above7 <- rle(means > 7)
        
        # n_above3 <- sum(above3$values & above3$lengths >= 5)
        # n_above5 <- sum(above5$values & above5$lengths >= 5)
        # n_above7 <- sum(above7$values & above7$lengths >= 5)
        
        above3 <- means > 3
        n_above3 <- sum(rowSums(embed(above3, 5)) == 5)
        
        above5 <- means > 5
        n_above5 <- sum(rowSums(embed(above5, 5)) == 5)
        
        above7 <- means > 7
        n_above7 <- sum(rowSums(embed(above7, 5)) == 5)
        
        if(n_above7 > n_above3 | n_above7 > n_above5 | n_above5 > n_above3){stop()}
        
        hss <- data.frame(rcp = rcp, pos = i, decade = mean(dec), site = site, n_above3, n_above5, n_above7)
        hotsummer_strings <- rbind(hotsummer_strings, hss)
      }
    }
  }
}

saveRDS(hotsummer_strings, file = file.path(wd, 'output', 'hotsummer_strings.rds'))
