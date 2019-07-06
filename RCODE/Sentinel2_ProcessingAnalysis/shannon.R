library(raster)
library(rgdal)
library(fasterize)
library(sf)
library(dplyr)

#shannon

shannon <- function (x, na.rm=TRUE){
  x <- na.omit(x)
  tb <- table(x)
  tbrel <- tb / sum(tb)
  z <- tbrel*log(tbrel)
  out <- -1*sum(z)
  return(out)
}




#birdsgrid <- readOGR("D:/Dropbox/MasterThesisJoséSilva/DATA/FieldData/Birds/Grids/Grid200mEffectSampSSU_WGS84UTM29N_VEZ_v1.shp")
birdsgrid <- read_sf("D:/Dropbox/MasterThesisJoséSilva/DATA/FieldData/Birds/Grids/Grid200mEffectSampSSU_WGS84UTM29N_VEZ_v1.shp")

length(unique(birdsgrid$ID_SSU))


rst <- stack(c("./DATAtoShare/kmeans_winter15_10c.grd", 
               "./DATAtoShare/kmeans_winter15_20c.grd"))

names(rst) <- c("wint15_km10c","wint15_km20c")

birdGridRst <- fasterize(birdsgrid,rst[[1]],field = "ID_SSU")

rstDF <- na.omit(values(stack(birdGridRst, rst)))

head(rstDF)

rstDFsummary <- rstDF %>% 
  as.data.frame() %>% 
  group_by(layer) %>% 
  summarise_all(.funs = list(shn=shannon), na.rm=TRUE) 
  #mutate(ID_SSU = layer)

colnames(rstDFsummary)[1] <- "ID_SSU"

write.csv(rstDFsummary, "winter15_shn.csv",row.names = FALSE)



gridBirds <- birdsgrid %>% left_join(rstDFsummary, by = "ID_SSU")

plot(gridBirds)

write_sf(gridBirds, "./OUT/shn_wint15.shp")

hist(rstDFsummary$wint15_km10c_shn)

#kmeans 5c and 15c

birdsgrid <- read_sf("D:/Dropbox/MasterThesisJoséSilva/DATA/FieldData/Birds/Grids/Grid200mEffectSampSSU_WGS84UTM29N_VEZ_v1.shp")

length(unique(birdsgrid$ID_SSU))


rst <- stack(c("./DataToShare/kmeans_NDVI_5c_v2.tif", 
               "./DataToShare/kmeans_NDVI_15c_v2.tif"))

names(rst) <- c("NDVI_km5c","NDVI_km15c")

birdGridRst <- fasterize(birdsgrid,rst[[1]],field = "ID_SSU")

rstDF <- na.omit(values(stack(birdGridRst, rst)))

head(rstDF)

rstDFsummary <- rstDF %>% 
  as.data.frame() %>% 
  group_by(layer) %>% 
  summarise_all(.funs = list(shn=shannon), na.rm=TRUE) 
#mutate(ID_SSU = layer)

colnames(rstDFsummary)[1] <- "ID_SSU"

write.csv(rstDFsummary, "NDVI_shn_v2.csv",row.names = FALSE)



gridBirds <- birdsgrid %>% left_join(rstDFsummary, by = "ID_SSU")

plot(gridBirds)

write_sf(gridBirds, "./OUT/shn_NDVI_v2.shp")


# 
# extRst_AVG <- extract(rst[[1:3]], birdsgrid, fun = mean, na.rm=TRUE, method="simple")
# extRst_STD <- extract(rst[[1:3]], birdsgrid, fun = sd, na.rm=TRUE, method="simple")
# 
# colnames(extRst_AVG) <- paste("Wint15_",colnames(extRst_AVG),"_avg",sep="")
# colnames(extRst_STD) <- paste("Wint15_",colnames(extRst_STD),"_std",sep="")
# 
# extRst <- cbind(birdsgrid@data[,c("ID_PSU","ID_SSU","PSU_SSU_ID")], 
#                 extRst_AVG, extRst_STD)
# 
# write.csv(extRst, "winter15_extrctPCA.csv",row.names = FALSE)
# 
