library(RStoolbox)
library(raster)
library(sf)
library(ggplot2)

#k means Raster

shapeArea <- read_sf("D:/Thesis/DATA/Vector/VezBasinCOS12_WGS84_29N.shp")

shapeAreaBuff <- st_buffer(shapeArea, 2500)


fl_10m <- list.files("D:/Thesis/DATA/Raster/SentinelS2A/SENTINEL2A_20161116-112335-459_L2A_T29TNG_D",
                     pattern="_SRE_", full.names = TRUE)[c(3,4,5,9)]

rst20m_resamp <- raster("D:/Thesis/OUT/autumn16/resamp20_20161116.tif")

ndvi <- raster("D:/Thesis/OUT/autumn16/NDVI_autumn16_20160801.tif")

rst10m <- stack(fl_10m)
rst10m <- crop(rst10m, shapeAreaBuff)
rst10m <- mask(rst10m, shapeAreaBuff)

rst<- stack(rst10m, rst20m_resamp, ndvi)

# Calibrate kmeans algorithm
# nSamples - Integer. Number of random samples to draw to fit cluster map
#
uc <- unsuperClass(rst, nClasses = 10, nSamples = 1E5)

map <- predict(uc, rst)



