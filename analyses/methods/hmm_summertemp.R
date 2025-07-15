


# Setup
rm(list = ls())
options(stringsAsFactors = FALSE)
library(rstan)
library(ggplot2)
rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains

wd <- "~/projects/seedsarentreal/analyses/methods"
setwd(wd)
util <- new.env()
source(file.path(wd, "mcmc_analysis_tools_rstan.R"), local=util) 
source(file.path(wd, "mcmc_visualization_tools.R"), local=util) 

kippenberger <- c("#8B174DFF", "#AE2565FF", "#C1447EFF", "#D06C9BFF", "#DA9FB8FF", 
                  "#ADBE7CFF", "#8BA749FF", "#6E8537FF", "#4F5F28FF", "#343D1FFF")



# Data exploration
raw_data <- read.csv(file.path(wd, 'data',  'ebms', 'Beech_tree-ring_masting_data.csv'))

length(unique(paste0(raw_data$site.ID,raw_data$tree.ID))) # 57 unique trees
length(unique(raw_data$site.ID)) # 7 unique sites

raw_data$uniqueID <- paste0(raw_data$site.ID, "_", raw_data$tree.ID)

ggplot(data = raw_data) + 
  facet_wrap(~ uniqueID) +
  geom_line(aes(x = year, y = seeds, group = uniqueID, color = uniqueID),
            linewidth = 0.2, alpha = 1) +
  geom_vline(data = raw_data[is.na(raw_data$seeds),], aes(xintercept = year), color = "grey95",
             linewidth = 1.5, alpha = 1) +
  scale_color_manual(values = colorRampPalette(kippenberger)(length(unique(raw_data$uniqueID)))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        strip.text.x.top = element_blank(),
        axis.text = element_text(size = 7))

ggplot(data = raw_data) + 
  facet_wrap(~ uniqueID) +
  geom_histogram(aes(x = seeds, fill = uniqueID), bins = 20, alpha = 1, color = "white") +
  scale_fill_manual(values = colorRampPalette(kippenberger)(length(unique(raw_data$uniqueID)))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        strip.text.x.top = element_blank())

ggplot(data = raw_data) + 
  geom_histogram(aes(x = seeds), bins = 40, alpha = 1, color = "white") +
  scale_fill_manual(values = colorRampPalette(kippenberger)(length(unique(raw_data$uniqueID)))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        strip.text.x.top = element_blank())

# Missing data? Not mentioned in HP et al. 2025... nor in https://onlinelibrary.wiley.com/doi/10.1111/gcb.16730...

## as first brute-force approach, we keep only trees with all obs. between 1985 and 2020
length(1986:2020)
data.frame(table(raw_data[,c("uniqueID")]))
raw_data_filter <- na.omit(raw_data[raw_data$year %in% c(1986:2020),])
trees_nobs <- data.frame(table(raw_data_filter[,c("uniqueID")]))
names(trees_nobs) <- c("uniqueID", "nobs")
raw_data_filter <- raw_data_filter[raw_data_filter$uniqueID %in% trees_nobs[trees_nobs$nobs == length(1986:2020), "uniqueID"],]

ggplot(data = raw_data_filter) + 
  geom_line(aes(x = year, y = seeds, group = uniqueID, color = uniqueID),
            linewidth = 0.2, alpha = 1) +
  geom_vline(data = raw_data_filter[is.na(raw_data_filter$seeds),], aes(xintercept = year), color = "grey95",
             linewidth = 1.5, alpha = 1) +
  scale_color_manual(values = colorRampPalette(kippenberger)(length(unique(raw_data_filter$uniqueID)))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        strip.text.x.top = element_blank(),
        axis.text = element_text(size = 7))

ggplot(data = raw_data_filter) + 
  geom_histogram(aes(x = seeds), bins = 60, alpha = 1, color = "white", fill = "grey80") +
  scale_fill_manual(values = colorRampPalette(kippenberger)(length(unique(raw_data$uniqueID)))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        strip.text.x.top = element_blank())



# Sizes
uniq_tree_ids <- unique(raw_data_filter$uniqueID)
N_trees <- length(uniq_tree_ids)

# Format data into ragged arrays
seed_counts <- c()
prev_summer_temp <- c()
sizes <- c()
years <- c()
N_years <- c()

idx <- 1
tree_start_idxs <- c()
tree_end_idxs <- c()

for (tid in uniq_tree_ids) {
  raw_data_tree <- raw_data_filter[raw_data_filter$uniqueID == tid,]
  
  years_tree <- raw_data_tree$year
  N_years_tree <- length(years_tree)
  
  size_tree <- 1 # the ground below each tree was searched for seeds, no info on area
  seed_count_tree <- raw_data_tree$seeds
  
  seed_counts <- c(seed_counts, seed_count_tree)
  sizes <- c(sizes, size_tree)
  years <- c(years, years_tree)
  N_years <- c(N_years, N_years_tree)
  
  tree_start_idxs <- c(tree_start_idxs, idx)
  idx <- idx + N_years_tree
  tree_end_idxs <- c(tree_end_idxs, idx - 1)
  
  prevsummertemptree <- raw_data[raw_data$uniqueID == tid & raw_data$year %in% c(1985:2019), 'summer.temp']
  if(length(prevsummertemptree) != N_years_tree){
    stop()
  }
  prev_summer_temp <- c(prev_summer_temp, prevsummertemptree)
  
}



