


# Create hypothetical trees c/w. summer alternations
years_to_predict <- seq(1980,  2010, 1)
summertemp_ww <- rnorm(length(years_to_predict), 24, 2/2.57)
summertemp_cw <- c(rnorm(9, 24, 2/2.57), rnorm(1, 15, 1/2.57), rnorm(1, 24, 2/2.57),  rnorm(1, 15, 1/2.57),
                   rnorm(1, 24, 2/2.57),  rnorm(1, 15, 1/2.57), rnorm(1, 24, 2/2.57),  rnorm(1, 15, 1/2.57),
                   rnorm(1, 24, 2/2.57),  rnorm(1, 15, 1/2.57),  rnorm(13, 24, 2/2.57))
                   

trees_per_stand <- 50
unique_stands <- "Benwell"
newtree_cw_stand_idxs <- rep(which(unique_stands==unique_stands), each = trees_per_stand)
N_newtrees_cw <- length(newtree_cw_stand_idxs)
N_max_newyears_cw <- length(years_to_predict)
first_newyear_cw <- min(years_to_predict)
newyears_cw <- c()
newprevsummer_temps_cw <-c()
newprevsummer_temps_ww <-c()
N_newyears_cw <- c()
idx <- 1
newtree_cw_start_idxs <- c()
newtree_cw_end_idxs <- c()
for (i in 1:N_newtrees_cw) {
  
  newyears_tree <- years_to_predict-first_newyear+1
  newyears_cw <- c(newyears_cw, newyears_tree)
  
  N_newyears_tree <- length(newyears_tree)
  N_newyears_cw <- c(N_newyears_cw, N_newyears_tree)
  
  newprevsummer_temps_tree <- summertemp_cw
  newprevsummer_temps_cw <- c(newprevsummer_temps_cw, newprevsummer_temps_tree)
  
  newprevsummer_temps_tree <- summertemp_ww
  newprevsummer_temps_ww <- c(newprevsummer_temps_ww, newprevsummer_temps_tree)
  
  newtree_cw_start_idxs <- c(newtree_cw_start_idxs, idx)
  idx <- idx + N_newyears_tree
  newtree_cw_end_idxs <- c(newtree_cw_end_idxs, idx - 1)
}

# Collect data
N <- length(years)
Nnew <- length(newyears)
Nnew_cw <- length(newyears_cw)
data <- mget(c('N', 'N_trees', 'N_max_years',
               'N_years', 'tree_start_idxs', 'tree_end_idxs',
               'seed_counts', 'years', 
               'prevsummer_temps', 'spring_temps', 'gdd_lastfrost',
               # for predictions
               'Nnew', 'N_newtrees', 'N_max_newyears', 'newtree_stand_idxs',
               'N_newyears', 'newtree_start_idxs', 'newtree_end_idxs','newyears', 
               'newprevsummer_temps', 'newgdd_lastfrost', 'newspring_temps',
               
               # for predictions
               'Nnew_cw', 'N_newtrees_cw', 'N_max_newyears_cw', 'newtree_cw_stand_idxs',
               'N_newyears_cw', 'newtree_cw_start_idxs', 'newtree_cw_end_idxs','newyears_cw', 
               'newprevsummer_temps_cw', 'newprevsummer_temps_ww'
))


sm_gq <- stan_model(file = file.path(wd, "stan", "paper_models", "model2_treelevel_zinb_missing_rescaled_summertemp_springfrost_springtemp_wpred_standlevel.stan"))
test <- gqs(sm_gq, data = data, as.matrix(fit_biol))
samples <- util$extract_expectand_vals(test)
