
library(raster)
library(rgdal)
library(fasterize)
library(sf)
library(dplyr)

#birdsgrid <- readOGR("D:/Dropbox/MasterThesisJoséSilva/DATA/FieldData/Birds/Grids/Grid200mEffectSampSSU_WGS84UTM29N_VEZ_v1.shp")
birdsgrid <- read_sf("D:/Dropbox/MasterThesisJoséSilva/DATA/FieldData/Birds/Grids/Grid200mEffectSampSSU_WGS84UTM29N_VEZ_v1.shp")

length(unique(birdsgrid$ID_SSU))


rst <- stack("./OUT/PCA_S2a_autumn2016.tif")

birdGridRst <- fasterize(birdsgrid,rst[[1]],field = "ID_SSU")

rstDF <- na.omit(values(stack(birdGridRst, rst[[1:3]])))

head(rstDF)

rstDFsummary <- rstDF %>% 
  as.data.frame() %>% 
  group_by(layer) %>% 
  summarise_all(.funs = list(avg=mean, std=sd), na.rm=TRUE) 
  #mutate(ID_SSU = layer)

colnames(rstDFsummary)[1] <- "ID_SSU"
 
write.csv(rstDFsummary, "autumn16_extrctPCA.csv",row.names = FALSE)


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
