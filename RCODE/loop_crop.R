

library(raster)
library(sf)

shapeArea <- read_sf("./DATA/Vector/VezBasinCOS12_WGS84_29N.shp")

shapeAreaBuff <- st_buffer(shapeArea, 2500)


#files <- list.files(path="DATA/Raster/SentinelS2A",  full.names=TRUE, recursive=FALSE)

fileList <- list.files(path="DATA/Raster/SentinelS2A",  full.names=TRUE, recursive=FALSE)

for(fname in fileList){
  
  print(fname)
  
  fl_10m <- list.files(fname, pattern="_SRE_", full.names = TRUE)[c(3,4,5,9)] # Input here is each Sentinel2 scene/folder
  fl_20m <- list.files(fname, pattern="_SRE_", full.names = TRUE)[-c(3,4,5,9)]
  
  
  print("Stacking and croping 20m bands!")
  rst10m <- stack(fl_10m)
  rst10m <- crop(rst10m, shapeAreaBuff)
  rst10m <- mask(rst10m, shapeAreaBuff)
  
  print("Stacking, croping and resampling 20m bands!")
  rst20m <- stack(fl_20m)
  rst20m <- crop(rst20m, shapeAreaBuff)
  rst20m <- mask(rst20m, shapeAreaBuff)
  # Resample 20 m bands
  rst20m_resamp <- resample(rst20m, rst10m, method="ngb")
  
  print("Stacking and writing all bands!")
  # Stack 10 m and resampled 20 m variables
  rst<- stack(rst10m, rst20m_resamp)
  
  ##
  ##
  # The output folder needs to be changed??! |||| <----------------------
  # Change this: "./DATA/RASTER/S2_L2A_CropResamp/"
  fnameOut <- paste("./DATA/Raster/S2A_crop/",basename(fname),"_crop_r10m.tif",sep="")
  fnameNDVI <- paste("./DATA/Raster/S2A_crop/",basename(fname),"_NDVI_crop_r10m.tif",sep="")
  
  # Write full sentinel multiband image
  writeRaster(rst, fnameOut, overwrite=FALSE)
  
  # Write NDVI image
  ndvi <- (rst10m[[4]] - rst10m[[3]]) / (rst10m[[4]] + rst10m[[3]])
  writeRaster(ndvi, fnameNDVI, overwrite=FALSE)
  print("Finished!!")
}

  
