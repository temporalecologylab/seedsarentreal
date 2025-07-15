# Prepare climate data

rm(list=ls())
library(terra)

wd <- "/home/vdmeersch/projects/masting"

# Read sites (from a file I found in the repo)
sites <- read.csv(file.path(wd, 'data', 'sites_ebms.csv'))
sites$id <- 1:nrow(sites)
# sites_vect <-  vect(sites, geom=c("longitude", "latitude"))
# sites_vect <- shift(rotate(sites_vect, longitude = 0),360)
sites_vect <- readRDS(file = file.path(wd, 'data', "files_vect.rds")) # just because Norm does not want to isntall the last version of terra.......
nsites <- nrow(sites)

years <- 1979:2024
nyears <- length(years)

# climate.data <- data.frame(Site = NA, Latitude = NA, Year = NA, DOY = NA,
#                            TAV = NA, TMIN = NA, TMAX = NA)

climdir<- "/home/vdmeersch/data/climate"

climate.data <- data.frame()
for(y in years){
  cat(paste0(y, '\n'))
  
  era5_tmax <- rast(file.path(climdir, "era5land", paste0("tmx_",y,".nc")))
  time(era5_tmax) <- seq(as.Date(paste0(y, "-01-01")), as.Date(paste0(y, "-12-31")), by = "day")
  
  e <- terra::extract(era5_tmax[[1]], sites_vect, search_radius = 20000, xy = TRUE) # 20km > 0.1deg, to get the closest cell if required
  sites_vect_e <- vect(e, geom = c('x','y'))
  
  start_ind <- match(as.Date(paste0(y, "-06-01")), time(era5_tmax))
  end_ind <- match(as.Date(paste0(y, "-07-31")), time(era5_tmax))
  
  tmax <- terra::extract(era5_tmax[[start_ind:end_ind]], sites_vect_e, ID = FALSE)
  tmax <- rowMeans(tmax) - 273.15
  
  clim <-
    data.frame(site.ID = sites$site.ID, lat = sites$latitude, year = y, meantmax_jj = tmax)

  
  climate.data <- rbind(climate.data, clim)
}

saveRDS(climate.data, file = "/home/vdmeersch/projects/masting/data/era5land_sitesextract.rds")