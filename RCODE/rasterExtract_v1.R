
library(raster)
library(rgdal)
library(readxl)

# read_excel("",x)

#imgPath <-"C:/Users/Utilizador/AppData/Local/Temp/processing_18c0361eafc14768be4578fa0990aa0d/24d9c78479a24524979610bc97101ea6/OUTPUT.tif"

fl <- list.files("D:/Thesis/DATA/Raster/SentinelS2A/SENTINEL2A_20161116-112335-459_L2A_T29TNG_D",
                     pattern="_FRE_", full.names = TRUE)[c(3,4,5,9)]

ndvi <- raster("D:/Thesis/OUT/NDVI_autumn16_20160801.tif")

rst <- stack(fl)

names(rst) <- paste("autumn16",c("b02","b03","b04","b08"),sep="_")

#bands

birdsgrid <- readOGR("D:/Dropbox/MasterThesisJoséSilva/DATA/FieldData/Birds/Grids/Grid200mEffectSampSSU_WGS84UTM29N_VEZ_v1.shp")

extRst_AVG <- extract(rst, birdsgrid, fun=mean)

colnames(extRst_AVG) <- paste(colnames(extRst_AVG),"_avg",sep="")

extRst_STD <- extract(rst, birdsgrid, fun=sd)

colnames(extRst_STD) <- paste(colnames(extRst_STD),"_std",sep="")

#ndvi

extRst_AVG_NDVI <- extract(ndvi, birdsgrid, fun=mean)

colnames(extRst_AVG_NDVI) <- paste(colnames(extRst_AVG_NDVI),"autumn16_NDVI_avg",sep="")

extRst_STD_NDVI <- extract(ndvi, birdsgrid, fun=sd)

colnames(extRst_STD_NDVI) <- paste(colnames(extRst_STD_NDVI),"autumn16_NDVI_std",sep="")

#save data

extRst <- cbind(birdsgrid@data[,c("ID_PSU","ID_SSU","PSU_SSU_ID")], 
                    extRst_AVG, extRst_STD, extRst_AVG_NDVI, extRst_STD_NDVI)

write.csv(extRst, "./OUT/autumn16_extrct.csv",row.names = FALSE)



