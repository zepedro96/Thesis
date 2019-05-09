

library(raster)
library(lubridate)

fl <- list.files("./DATA/RASTER/S2_Vez/S2A_crop/NDVI", pattern=".tif$", full.names = TRUE)

dt <- substr(basename(fl), 12, 19)
yr <- substr(dt,1,4)
mo <- substr(dt,5,6)
dy <- substr(dt,7,8)

dts <- date(paste(yr,mo,dy,sep="-"))

monthNum <- month(dts)
monthNum[1:2] <- 0
monthNum <- monthNum + 1

rstNDVI <- stack(fl)

rstNDVImax <- stackApply(rstNDVI, monthNum, fun=max)

#plot(rstNDVImax[[10]])

for(i in 1:nlayers(rstNDVImax)){
  
  writeRaster(rstNDVImax[[i]], filename = paste("./DATA/RASTER/S2_Vez/S2A_crop/NDVImax/NDVI_Max_m",i-1,".tif",sep=""))
  print(i)
  
}


