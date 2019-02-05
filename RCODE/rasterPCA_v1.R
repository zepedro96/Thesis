library(RStoolbox)
library(raster)
library(sf)

#PCA Raster

shapeArea <- read_sf("D:/Thesis/DATA/Vector/VezBasinCOS12_WGS84_29N.shp")

shapeAreaBuff <- st_buffer(shapeArea, 2500)


fl_10m <- list.files("D:/Thesis/DATA/Raster/SentinelS2A/SENTINEL2A_20160124-114205-331_L2A_T29TNG_D",
                     pattern="_SRE_", full.names = TRUE)[c(3,4,5,9)]

rst20m_resamp <- raster("D:/Thesis/OUT/Winter15/resamp20_20160124.tif")

ndvi <- raster("D:/Thesis/OUT/Winter15/NDVI_winter15_20160124.tif")

rst10m <- stack(fl_10m)
rst10m <- crop(rst10m, shapeAreaBuff)
rst10m <- mask(rst10m, shapeAreaBuff)

rst<- stack(rst10m, rst20m_resamp, ndvi)


ggRGB(rst, 1,2,3)
set.seed(112)
RstPCA <- rasterPCA(rst,  nSamples = NULL, nComp = 10, spca = FALSE, maskCheck = TRUE)

summary(RstPCA$model)

test_PCA <- ggRGB(RstPCA$map,1,2,3, stretch="lin", q=0)

