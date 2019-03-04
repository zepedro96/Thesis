
library(raster)
library(rgdal)


birdsgrid <- readOGR("C:/Users/Utilizador/Desktop/TESE/MasterThesisJoseSilva/MasterThesisJoseSilva/DATA/fielddata/Birds/Grids","Grid200mEffectSampSSU_WGS84UTM29N_VEZ_v1")

rst <- stack("./OUT/PCA_S2a_Winter2015.tif")

extRst_AVG <- extract(rst[[1:3]], birdsgrid, fun = mean, na.rm=TRUE, method="simple")
extRst_STD <- extract(rst[[1:3]], birdsgrid, fun = sd, na.rm=TRUE, method="simple")

colnames(extRst_AVG) <- paste("Wint15_",colnames(extRst_AVG),"_avg",sep="")
colnames(extRst_STD) <- paste("Wint15_",colnames(extRst_STD),"_std",sep="")

extRst <- cbind(birdsgrid@data[,c("ID_PSU","ID_SSU","PSU_SSU_ID")], 
                extRst_AVG, extRst_STD)

write.csv(extRst, "Spring16_extrctBands.csv",row.names = FALSE)

