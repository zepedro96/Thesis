

library(raster)
library(rgdal)
library(fasterize)
library(sf)
library(dplyr)

#birdsgrid <- readOGR("D:/Dropbox/MasterThesisJoséSilva/DATA/FieldData/Birds/Grids/Grid200mEffectSampSSU_WGS84UTM29N_VEZ_v1.shp")
birdsgrid <- read_sf("./DATA/FieldData/Birds/Grids/Grid200mEffectSampSSU_WGS84UTM29N_VEZ_v1.shp")


r <- raster("./DATAtoShare/RASTER/kmeans_NDVI_20c_v2.tif")


birdGridRst <- fasterize(birdsgrid,rst[[1]],field = "ID_SSU")

rstDF <- na.omit(values(stack(birdGridRst, rst)))

head(rstDF)


