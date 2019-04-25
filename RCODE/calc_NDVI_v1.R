#median, MaxVal & mad

library(raster)
library(sf)
library(dplyr)
library(rgdal)
fileList <- list.files(path="DATA/Raster/S2A_crop",  full.names=TRUE, recursive=FALSE, pattern = "_NDVI")

stackNDVI <- stack(fileList)

rstMedian <- calc(stackNDVI, fun = median, na.rm=TRUE)

rstMax <- calc(stackNDVI, fun = max, na.rm=TRUE)

rstMAD <- calc(stackNDVI, fun = mad, na.rm=TRUE)


writeRaster(rstMedian,"./OUT/AnnualMetrics/Non-smoothed/NDVI_Median_2016.tif")
writeRaster(rstMax,"./OUT/AnnualMetrics/Non-smoothed/NDVI_Max_2016.tif")
writeRaster(rstMAD,"./OUT/AnnualMetrics/Non-smoothed/NDVI_MAD_2016.tif")


#extract 

birdsgrid <- readOGR("D:/Dropbox/MasterThesisJoséSilva/DATA/FieldData/Birds/Grids/Grid200mEffectSampSSU_WGS84UTM29N_VEZ_v1.shp")

extNDVI_med_mn <- extract(rstMedian, birdsgrid, fun=mean)
extNDVI_med_sd <- extract(rstMedian, birdsgrid, fun=sd)

extNDVI_max_mn <- extract(rstMax, birdsgrid, fun=mean)
extNDVI_max_sd <- extract(rstMax, birdsgrid, fun=sd)

extNDVI_MAD_mn <- extract(rstMAD, birdsgrid, fun=mean)
extNDVI_MAD_sd <- extract(rstMAD, birdsgrid, fun=sd)

birdsDF <- birdsgrid %>% 
  st_drop_geometry() %>% 
  select(ID_PSU,ID_SSU,PSU_SSU_ID)

birdsDFall <- data.frame(birdsDF,
           NDVI_med_mn = extNDVI_med_mn,
           NDVI_med_sd = extNDVI_med_sd,
           NDVI_max_mn = extNDVI_max_mn,
           NDVI_max_sd = extNDVI_max_sd,
           NDVI_MAD_mn = extNDVI_MAD_mn,
           NDVI_MAD_sd = extNDVI_MAD_sd)

birdsDFall <- write.csv(birdsDFall, "OUT/AnnualMetrics/Non-smoothed/NS_birds_all.csv")

