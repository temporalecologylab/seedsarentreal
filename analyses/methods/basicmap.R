# Setup
rm(list = ls())
options(stringsAsFactors = FALSE)
library(rstan)
library(ggplot2)
rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains

library(rnaturalearth)
library(terra)
library(tidyterra)
library(ggplot2)

wd <- "~/projects/seedsarentreal/analyses/methods"

sites <- read.csv(file.path(wd, 'data', 'ebms', 'sites_ebms.csv'))
sites <- vect(sites, geom = c('longitude', 'latitude'))
world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
crs(sites) <- crs(world)


ggplot() +
  geom_spatvector(data = crop(world, ext(c(-7, 3, 49,57))), fill = 'grey95', color = 'grey80') +
  geom_spatvector(data = sites, shape = 21, color = 'white', fill = '#be2d2d', size = 5, stroke = 2) + 
  scale_size_manual(values = c(2, 2.3), guide = 'none') +
  coord_sf(expand = FALSE) +
  theme_bw() +
  theme(legend.position = 'inside', legend.position.inside = c(0.1, 0.2),
        legend.background = element_rect(color = 'black', linewidth = 0.2),
        legend.text = element_text(size = 9, margin = margin(t = 5, r = 5, b = 5, l = 0, unit = "pt")), 
        legend.title = element_text(size = 9), panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank()) 
