
library(raster)
library(sf)

vezDF1km <- read.csv("D:/MyDocs/GeoData/ByProject/IND_CHANGE/FieldWorkData/CommonBirds_FieldData_MainTable/birdDiv_multiSARavg_preds_Vez_v2.csv")


fpaths <- 
  paste("./DATA/RASTER/Landsat/NDVI-32day/",
        c("LT8_NDVI_2013-19_EFAamp.tif",
          "LT8_NDVI_2013-19_EFAmin.tif",
          "LT8_NDVI_2013-19_EFAmax.tif",
          "LT8_NDVI_2013-19_EFAavg.tif",
          "LT8_NDVI_2013-19_EFAstd.tif",
          "LT8_NDVI_2013-19_EFAmed.tif",
          "LT8_NDVI_2013-19_EFAcv.tif",
          "LT8_NDVI_2013-19_EFAsprg.tif",
          "LT8_NDVI_2013-19_EFAwint.tif",
          "LT8_NDVI_2013-19_ESPI.tif"
        ),sep=""
  )


for(i in 1:length(fpaths)){
  
  print(i)
  tmp <- stack(fpaths[[i]])[[2]]
  if(i==1)
    rst <- stack(tmp)
  else
    rst <- stack(rst, tmp)
}

names(rst) <- c("EFAamp","EFAmin","EFAmax","EFAavg","EFAstd",
                "EFAmed","EFAcv","EFAsprg","EFAwint","ESPI")


#grid1kmall <- read_sf("D:/MyDocs/GeoData/ByProject/IND_CHANGE/SampDesign/SampDesignVEZ_0_INDCHANGE_v1/SampGridVez_1km_Vez_All_v1.shp") 
grid1kmall <- read_sf("C:/Users/JG/Desktop/birdData_GBIF/Grid1km_Vez_GBIFdata_v1.shp")

gridSub1km <- grid1kmall %>% 
  filter(ID_PSU %in% DF$ID_PSU)

birdGridRst1km <- fasterize::fasterize(gridSub1km, rst[[1]], field = "ID_PSU")

plot(birdGridRst1km)

rstDF1km <- na.omit(values(stack(birdGridRst1km, rst))) %>% 
  as.data.frame %>% 
  rename("ID_PSU" = "layer") %>% 
  group_by(ID_PSU) %>% 
  summarize_all(.funs=list(avg = mean, std = sd))


vezDF_vars1km <- DF %>% 
  left_join(rstDF1km, by=c("ID_PSU"))

cm <- cor(vezDF_vars1km[,-1], method="pearson")

data.frame(cn=colnames(cm)[-c(1:2)],
           corVal=cm[-c(1:2),1,drop=FALSE]) %>% arrange(desc(abs(SpRich)))


plot(vezDF_vars1km$EFAcv_avg, vezDF_vars1km$SpRich)

cor.test(vezDF_vars1km$EFAavg_std, vezDF_vars1km$SpRich, method="pearson")



