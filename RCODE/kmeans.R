library(RStoolbox)
library(raster)
library(sf)
library(ggplot2)

#k means Raster

shapeArea <- read_sf("D:/Thesis/DATA/Vector/VezBasinCOS12_WGS84_29N.shp")

shapeAreaBuff <- st_buffer(shapeArea, 2500)


fl_10m <- list.files("D:/Thesis/DATA/Raster/SentinelS2A/SENTINEL2A_20160124-114205-331_L2A_T29TNG_D",
                     pattern="_SRE_", full.names = TRUE)[c(3,4,5,9)]

fl_20m <- list.files("D:/Thesis/DATA/Raster/SentinelS2A/SENTINEL2A_20160124-114205-331_L2A_T29TNG_D",
                     pattern="_SRE_", full.names = TRUE)[-c(3,4,5,9)]

rst10m <- stack(fl_10m)
rst20m <- stack(fl_20m)

rst10m <- crop(rst10m, shapeAreaBuff)
rst10m <- mask(rst10m, shapeAreaBuff)

rst20m <- crop(rst20m, shapeAreaBuff)
rst20m <- mask(rst20m, shapeAreaBuff)

rst20m_resamp <- resample(rst20m, rst10m, method="ngb")

rst<- stack(rst10m, rst20m_resamp)


## Files for functional variables

fl <- c("./OUT/AnnualMetrics/Non-smoothed/NDVI_MAD_2016.tif",
        #"./OUT/AnnualMetrics/Non-smoothed/NDVI_Max_2016.tif",
        "./OUT/AnnualMetrics/Non-smoothed/NDVI_Median_2016.tif")


rst<- stack(fl)

#rst <- scale(rst)

# Calibrate kmeans algorithm
# nSamples - Integer. Number of random samples to draw to fit cluster map
#
uc1 <- unsuperClass(rst, nClasses = 10, nSamples = 1E6)

map1 <- predict(uc1, rst, "./DataToShare/kmeans_NDVI_10c_v2.tif")

uc2 <- unsuperClass(rst, nClasses = 20, nSamples = 1E6)

map2 <- predict(uc2, rst, "./DataToShare/kmeans_NDVI_20c_v2.tif")



