
library(terra)
uk <- geodata::gadm(country = "United Kingdom", level = 0, resolution = 2, path = file.path(wd, 'data'))
e <- terra::ext(-8, 2, 49.5, 56)

pdf(file = file.path(figpath, 'supp', 'map.pdf'),
    width = 6.5, height = 4)

layout(matrix(c(1,2), nrow = 1),
       widths = c(1.1,1))

plot(uk, ext = e, axes = TRUE, box = TRUE, xaxs = "i", yaxs = "i",
     col = 'grey85', border = 'white', mar = c(1.5,1.5,1,1))
points(sites_vect, col = 'white', cex = 1.5)
points(sites_vect, col = util$c_mid, cex = 1)

mtext(LETTERS[1], side = 3, line = 1, adj = 0.05, cex = 1.2, font = 2, col = 'grey30')

par(mar = c(3.5,3,2,1), mgp = c(2,0.5,0.25), cex.axis = 0.8, cex.lab = 0.9, tcl = -0.2)
D <- distance(sites_vect, sites_vect, unit="km") 
D[upper.tri(D)] <- NA
D[D==0] <- NA
hist(D, breaks = 10, xlab = 'Pairwise distance (km)', main = '',
     col = util$c_light_highlight, border = 'white')

mtext(LETTERS[2], side = 3, line = -1, adj = -0.25, cex = 1.2, font = 2, col = 'grey30')

dev.off()