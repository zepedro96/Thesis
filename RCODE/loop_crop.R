

library(raster)
library(sf)

shapeArea <- read_sf("D:/Thesis/DATA/Vector/VezBasinCOS12_WGS84_29N.shp")

shapeAreaBuff <- st_buffer(shapeArea, 2500)


files <- list.files(path="DATA/Raster/SentinelS2A",  full.names=TRUE, recursive=FALSE)

lapply(files, fucntion( {
  fl_10m <- list.files(pattern="_SRE_", full.names = TRUE)[c(3,4,5,9)]
  fl_20m <- list.files(pattern="_SRE_", full.names = TRUE)[-c(3,4,5,9)]
  
  rst10m <- stack(fl_10m)
  rst20m <- stack(fl_20m)
  
  rst10m <- crop(rst10m, shapeAreaBuff)
  rst10m <- mask(rst10m, shapeAreaBuff)
  
  rst20m <- crop(rst20m, shapeAreaBuff)
  rst20m <- mask(rst20m, shapeAreaBuff)
  
  rst20m_resamp <- resample(rst20m, rst10m, method="ngb")
  
  rst<- stack(rst10m, rst20m_resamp)
  writeRaster(rst, "./DATA/Raster/SentinelS2A_crop_rspl/rst.tif", overwrite=FALSE)
  
  ndvi <- (rst10m[[4]] - rst10m[[3]]) / (rst10m[[4]] + rst10m[[3]])
  writeRaster(ndvi, "./DATA/Raster/SentinelS2A_crop_rspl/_NDVI.tif", overwrite=FALSE)
})
