
library(raster)
library(rgdal)
library(readxl)

# read_excel("",x)

#imgPath <-"C:/Users/Utilizador/AppData/Local/Temp/processing_18c0361eafc14768be4578fa0990aa0d/24d9c78479a24524979610bc97101ea6/OUTPUT.tif"

fl <- list.files("C:/Users/Utilizador/Downloads/S2A_MSIL1C_20160430T112122_N0201_R037_T29TNG_20160430T112639.SAFE/GRANULE/L1C_T29TNG_A004467_20160430T112639/IMG_DATA"
                 , pattern=".jp2$", full.names = TRUE)[c(2,3,4,8)]

rst <- stack(fl)

names(rst) <- paste("Spring16",c("b02","b03","b04","b08"),sep="_")


birdsgrid <- readOGR("C:/Users/Utilizador/Desktop/TESE/MasterThesisJoseSilva/MasterThesisJoseSilva/DATA/fielddata/Birds/Grids","Grid200mEffectSampSSU_WGS84UTM29N_VEZ_v1")

extRst_AVG <- extract(rst, birdsgrid, fun=mean)

colnames(extRst_AVG) <- paste(colnames(extRst_AVG),"_avg",sep="")

extRst_STD <- extract(rst, birdsgrid, fun=sd)

colnames(extRst_STD) <- paste(colnames(extRst_STD),"_std",sep="")


extRst <- cbind(birdsgrid@data[,c("ID_PSU","ID_SSU","PSU_SSU_ID")], 
                    extRst_AVG, extRst_STD)

write.csv(extRst, "Spring16_extrctBands.csv",row.names = FALSE)



