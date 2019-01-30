

library(raster)
library(sf)

shapeArea <- read_sf("D:/Thesis/DATA/Vector/VezBasinCOS12_WGS84_29N.shp")

shapeAreaBuff <- st_buffer(shapeArea, 2500)

fl_10m <- list.files("D:/Thesis/DATA/Raster/SentinelS2A/SENTINEL2A_20161229-113456-457_L2A_T29TNG_D_V1-4",
                 pattern="_SRE_", full.names = TRUE)[c(3,4,5,9)]

fl_20m <- list.files("D:/Thesis/DATA/Raster/SentinelS2A/SENTINEL2A_20161229-113456-457_L2A_T29TNG_D_V1-4",
                     pattern="_SRE_", full.names = TRUE)[-c(3,4,5,9)]


rst10m <- stack(fl_10m)
rst20m <- stack(fl_20m)

rst10m <- crop(rst10m, shapeAreaBuff)
rst10m <- mask(rst10m, shapeAreaBuff)

rst20m <- crop(rst20m, shapeAreaBuff)
rst20m <- mask(rst20m, shapeAreaBuff)

rst20m_resamp <- resample(rst20m, rst10m, method="ngb")

rst<- stack(rst10m, rst20m_resamp)


ndvi <- (rst10m[[4]] - rst10m[[3]]) / (rst10m[[4]] + rst10m[[3]])

writeRaster(ndvi, "./OUT/NDVI_test.tif")

